# SnaQ submission files 
# This script call snaq_onedata.jl to run SnaQ. 
# It allows to run specific replicate (red_start to rep_end). 

using ArgParse 
using TimerOutputs 
using Dates
using TimeZones 
using RCall 
using CSV
using DataFrames
@everywhere using Printf 
@everywhere using Distributed  
@everywhere include("utilities.jl")

const to = TimerOutput()  
tz = TimeZone("America/Chicago") 
current_time_tz = ZonedDateTime(now(), tz) 
time = Dates.format(current_time_tz, "yyyy-mm-dd HH:MM:SS zzz") 

function parse_commandline()
    s = ArgParseSettings() 
    @add_arg_table s begin 
      
      # Specify which rep and which parameter sets to run SnaQ: 
      "--dup_rate"
        help = "Paramater setting (duplication rate) to run snaq"
        arg_type = Float64
        required = true
      "--loss_rate"
        help = "Paramater setting (gene loss rate) to run snaq"
        arg_type = Float64
        required = true
      "--ratevar"
        help ="Parameter setting (variation rate) to run snaq"
        arg_type = String
        required = true
      "--n_reps"
        help = "How many reps this parameter set was run in total"
        arg_type = Int
        default = 100
      "--rep_start"
        help = "Parameter setting (start of the range of replicates) to run snaq"
        arg_type = Int
        default = 1 # default is 1 if not specified
      "--rep_end"
        help = "Parameter setting (end of the range of replicates) to run snaq"
        arg_type = Int
        default = -1 # temporary default. 
      "--n_inds"
        help = "Paramater setting (gene loss rate) to run snaq"
        arg_type = Int
        default = 1
      "--SF"
        help = "Scaling factor to scale effective population Ne
          (Default = 1.0, no scaling)"
        arg_type = Float64
        default = 1.0
      "--debug_mode" 
        help = "Debugging mode: Keep all temporary folders"
        default = false
        arg_type = Bool
      "--gene_len" 
        help = "The length of simulated gene sequences (Default = 1000 bp)" 
        arg_type = Int 
        default = 1000  

      # Specify arguments for SnaQ: 
      "--runs"
        help = "Number of runs in SnaQ"
        arg_type = Int
        required = true
      "--n_snaqboot_rep" # legancy parameters for snaq_bootstrap 
        help = "Number of replicate for SnaQ boostrap for Hmax = 1" 
        arg_type = Int
        default = 30 
      "--method" # snaq_bootstrap is a lagancy method, kept for compatibility 
        help = """Method in selecting the best networks using SnaQ. 
                  Options: QuartetNetworkGoodnessFit / snaq_bootstrap / Both"""
        arg_type = String
        default = "QuartetNetworkGoodnessFit" 
    end 
    
    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1 # If rep_end is not specified, set it to n_reps
        parsed_args["rep_end"] = parsed_args["n_reps"]
    end

    return parsed_args
  end
  
parsed_args = parse_commandline()

# Parse arguments:  
dup_rate = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"]
n_reps = parsed_args["n_reps"]
n_inds = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
runs = parsed_args["runs"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]
SF = parsed_args["SF"] 
debug_mode = parsed_args["debug_mode"]
gene_len = parsed_args["gene_len"] 

# method to compute networks between different Hmax: 
method = lowercase(parsed_args["method"])
if !(method in ["snaq_bootstrap", "quartetnetworkgoodnessfit", "both"])
    error("method should be either 'snaq_bootstrap' or 'quartetnetworkgoodnessfit' or 'both'")
end

# set up folders: 
paramname_root = set_up_paramname_root(dup_rate, loss_rate, ratevar, 
                                      n_inds, SF, gene_len)  
# Above: specify to find the folder

outfolder = "output/$paramname_root"

#-----------------------------------------------#       
#    Check pre-existing files and remove them
#-----------------------------------------------# 
folder_path_list = []
for simulation_rep in rep_start:rep_end
  rep_number_string = pad_number(simulation_rep, n_reps)
  rep_folder_path = joinpath(outfolder, "rep$rep_number_string") 
  push!(folder_path_list, rep_folder_path) 
end 

snaqfolder_list = []
index_length = rep_end - rep_start + 1 # match how we set up rep_start:rep_end, 
for ind in 1:index_length 
  snaqfolder = setup_rep_output_folders(folder_path_list, ind, "snaqfolder")
  push!(snaqfolder_list, snaqfolder)
end 
check_existing_dir(snaqfolder_list) # see utilies.jl. 

#-----------------------------------------------#       
#       Estimate Networks using SnaQ
#-----------------------------------------------#
#--------------- Set up seeds -------------------#
#= Seed generation: 
# 1. Generate a master seed for each parameter setting using paramname_root
# 2. Generate a seed for each software using the master seed
# 3. Generate a n_rep x 4 seed array for SnaQ using the seed for SnaQ
# 4. For each replicate, select the corresponding row in the seed array to run SnaQ (see below)
=#
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
# unique master seed for each parameter setting: 
master_seed = generate_master_seed(params_dict_for_seed_setting) 

software_names = ["snaq"] 
seed_dic = generate_software_seeds(master_seed, software_names) # see utility.jl 
# This seed used to generate m (n_reps) x 2 seed array: 
seed_snaq = seed_dic["snaq"] 

# Generate a n_rep x 4 vector with random seeds generated from the master seed snaq_seed. 
# 1st seed -> infer net 0 
# 2nd seed -> infer net1
# 3rd seed -> QuartetNetworkGoodnessFit for H = 0 
# 4th seed -> QuartetNetworkGoodnessFit for H = 1 
# If using snaq_bootstrap: 
# snaq_boostrapping might not be used in the end, 
# kept it for legacy reason. 
# 5th seed -> infer net0 boostrapping (might not use it)
# 6th seed -> infer net1 boostrapping (might not use it)
seed_array = seed_generator(seed_snaq, n_reps, 6, outfolder, "random_seed_snaq.txt") 
# Here, seed_array is generated using n_reps x 6. 
# For each rep, seed is selected as [repID, i for i in 1:6] 

#-----------------------------------------------# 
#       Run SNaQ for each replicate
#-----------------------------------------------# 
# broadcast variables to all processors: 
@everywhere seed_array = $seed_array 
@everywhere snaqfolder_list = $snaqfolder_list 
@everywhere folder_path_list = $folder_path_list 
@everywhere outfolder = $outfolder
@everywhere paramname_root = $paramname_root
@everywhere runs = $runs
@everywhere n_snaqboot_rep = $n_snaqboot_rep
@everywhere n_inds = $n_inds
@everywhere method = $method
@everywhere rep_start = $rep_start
@everywhere rep_end = $rep_end
@everywhere method = $method 

@everywhere begin 
"""
  run_snaq_for_replicate(ind::Int)
Run SNaQ analysis for a specific replicate index.
# Arguments
- `ind::Int`: The index of the replicate to process.
# Description
This function constructs and executes a command to run the SNaQ analysis
  for the specified replicate. It utilizes various parameters and settings
  defined in the global scope, including output directories, seeds, and
  analysis options. The function is designed to be called in a parallel
  processing context, allowing multiple replicates to be processed concurrently.
# Note
- Ensure that all required global variables are defined and accessible
  in the scope where this function is called.
- The function assumes that the SNaQ script (`snaq_1rep.jl`) is located
  in the `./scripts/` directory relative to the current working directory.
"""
  function run_snaq_for_replicate(ind)

    iqtreefolder  = setup_rep_output_folders(folder_path_list, ind, "iqtreefolder")
    astralfolder = setup_rep_output_folders(folder_path_list, ind, "astralfolder")
    snaqfolder = snaqfolder_list[ind] 

    # Identify the correct rep to start with 
    ind_in_seed_array = ind + rep_start - 1 

    seed_net0 = seed_array[ind_in_seed_array, 1]
    seed_net1 = seed_array[ind_in_seed_array, 2]
    seed_qgof0 = seed_array[ind_in_seed_array, 3]
    seed_qgof1 = seed_array[ind_in_seed_array, 4]
    seed_net0_bs = seed_array[ind_in_seed_array, 5]
    seed_net1_bs = seed_array[ind_in_seed_array, 6]

    #= For snaq, each replicate should have a different seed.  
    Explanations for selecting seed_array[ind_in_seed_array: i for i in i:4] 
    For example 1: 
    For rep2, and rep_start = 2 and rep_end = 5 (index_length = 5 - 2 + 1 = 4)
    ind for rep2 in this loop is 1, so the actual repID = 1 + 2 - 1 = 2

    example 2: 
    For rep 5 and rep_start = 3 and rep_end = 10 (index_length = 10 - 3 + 1 = 8) 
    ind for rep5 in this loop is 3, so actual repID = 3 + 3 - 1 = 6 
    
    This method makes sure we selects the same seed from seed_array no matter of rep_start and rep_end 
    =# 
    run(`julia ./scripts/snaq_1rep.jl \
        --outfolder $outfolder \
        --iqtreefolder $iqtreefolder \
        --astralfolder $astralfolder \
        --snaqfolder $snaqfolder \
        --paramname_root $paramname_root \
        --seed_net0 $seed_net0 \
        --seed_net1 $seed_net1 \
        --seed_qgof0 $seed_qgof0 \
        --seed_qgof1 $seed_qgof1 \
        --seed_net0_bs $seed_net0_bs \
        --seed_net1_bs $seed_net1_bs \
        --runs $runs \
        --n_snaqboot_rep $n_snaqboot_rep \
        --n_inds $n_inds \
        --method $method`) 
  end 
end 

data_frames = []  # A list of data frames to collect results from each replicate 

@timeit to "Running SNaQ from rep$rep_start to rep$rep_end" begin  
  
  pmap(ind -> begin
      println("Worker $(myid()): Starting task $ind")
      run_snaq_for_replicate(ind)
  end, 1:index_length)
  # index_length = rep_end - rep_start + 1 
  # pmap will automatically distribute 
  
end

#-----------------------------------------------#       
#       Ouput the running time into .log
#-----------------------------------------------#
host_name = gethostname() 

# Define the number of processors based on the Julia parallel environment
processors = nprocs()

# Update the SNAQ_arguments to include the processors variable
SNAQ_arguments = """
  #=====================================================#
  #------------------------SNaQ-------------------------#
  #=====================================================#
  ---Arguments used to run SNAQ and SNaQ bootstrap---
  runs = $runs, number of runs for snaq on the original data & each bootstrap reps;
  seed_snaq = $seed_snaq, Master seed to generate seed array for SNaQ;
  n_snaqboot_rep = $n_snaqboot_rep, Number of bootstrap replicates for SNaQ (Hmax=0,1);
  --- Other Information ---
  processors = $processors, Number of processors to run SNaQ (and bootsnaq);
  Server for running the script = $host_name.stat.wisc.edu
  Time of running the script = $time;
  --- Running time ---
  """

argument_file = joinpath(outfolder, "arguments-$paramname_root.log")

open(argument_file, "a") do io
    println(io, SNAQ_arguments)  # Append SNAQ arguments
    show(io, to)  # Append timer output
end

println("=============================================")
println("SNAQ Analysis Completed!") 
println("SNAQ have been saved to $argument_file ")