# SnaQ submission files 
# This script call snaq_onedata.jl to run SnaQ. It allows to run specific replicate (red_start to rep_end). 

using ArgParse 
using TimerOutputs 
include("utilities.jl")

const to = TimerOutput()  

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
        required = true
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
      "--processors"
        help= "Number of processors to multi-process snaq"
        arg_type = Int
        default = nprocs() - 1 # -p will naturally add 1 to this maximum processors 
      
      # Specify arguments for SnaQ: 
      "--runs"
        help = "Number of runs in SnaQ"
        arg_type = Int
        required = true
      "--seed_snaq"
        help = "Seed to run SnaQ"
        arg_type = Int
        required = true
      "--n_snaqboot_rep" 
        help = "Number of replicate for SnaQ boostrap for Hmax = 1" 
        arg_type = Int
        required = true
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
seed_snaq = parsed_args["seed_snaq"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]
processors = parsed_args["processors"]

# set up folders: 
paramname_root = "DUP$dup_rate-LOS$loss_rate-RV$ratevar-N_ind$n_inds" # Specify to find the folder
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
index_length = rep_end - rep_start + 1 # see utilities.jl to match how we set up rep_start:rep_end, 
for ind in 1:index_length 
  snaqfolder = setup_rep_output_folders(folder_path_list, ind, "snaqfolder")
  push!(snaqfolder_list, snaqfolder)
end 
check_existing_dir(snaqfolder_list) # see utilies.jl. 

#-----------------------------------------------#       
#       Estimate Networks using SnaQ
#-----------------------------------------------#
# Generate a n_rep x 4 vector with random seeds generated from the master seed snaq_seed. 
# 1st seed -> infer net 0 
# 2nd seed -> infer net1
# 3rd seed -> infer net0 boostrapping 
# 4th seed -> infer net1 boostrapping 
seed_array = seed_generator(seed_snaq, n_reps, 4, outfolder, "random_seed_snaq.txt") 
# Here, seed_array is generated using n_reps x 4. 
# For each rep, seed is selected as [repID, i for i in 1:4] 

@timeit to "Running SNaQ from rep$rep_start to rep$rep_end" begin  
  for ind in 1:index_length 
    
    iqtreefolder  = setup_rep_output_folders(folder_path_list, ind, "iqtreefolder")
    astralfolder = setup_rep_output_folders(folder_path_list, ind, "astralfolder")
    snaqfolder = snaqfolder_list[ind] 

    ind_in_seed_array = ind + n_start - 1 # For snaq, each replicate should have a different seed.  
    #= Explanations for selecting seed_array[ind_in_seed_array: i for i in i:4] 
    For example 1: 
    For rep2, and rep_start = 2 and rep_end = 5 (index_length = 5 - 2 + 1 = 4)
    ind for rep2 in this loop is 1, so the actual repID = 1 + 2 - 1 = 2

    example 2: 
    For rep 5 and rep_start = 3 and rep_end = 10 (index_length = 10 - 3 + 1 = 8) 
    ind for rep5 in this loop is 3, so actual repID = 3 + 3 - 1 = 6 
    
    This method makes sure we selects the same seed from seed_array no matter of rep_start and rep_end 
    =# 

    seed_net0 = seed_array[ind_in_seed_array,1] # net0 
    seed_net1 = seed_array[ind_in_seed_array,2] # net1
    seed_net0_bs = seed_array[ind_in_seed_array, 3] # net0 boostrapping 
    seed_net1_bs = seed_array[ind_in_seed_array, 4] # net1 boostrapping 
        
    # run snaq_onedata.jl: 
    run(`julia -p $processors ./scripts/snaq_onedata.jl --outfolder $outfolder --iqtreefolder $iqtreefolder --astralfolder $astralfolder --snaqfolder $snaqfolder --paramname_root $paramname_root --seed_net0 $seed_net0 --seed_net1 $seed_net1 --seed_net0_bs $seed_net0_bs --seed_net1_bs $seed_net1_bs --runs $runs --n_snaqboot_rep $n_snaqboot_rep --processors $processors`) 
    # Are both processors are important? Need to test 
  end
end 

#-----------------------------------------------#       
#       Ouput the running time into .log
#-----------------------------------------------#
host_name = gethostname() 

SNAQ_arguments = """
  #--------------------------SNaQ--------------------------------#
  Arguments used to run SNAQ and SNaQ bootstrap:
  runs = $runs, number of independent runs for snaq on the original data & each bootstrap replicate;
  seed_snaq = $seed_snaq, Master seed to generate seed array for SNaQ;
  n_snaqboot_rep = $n_snaqboot_rep, Number of bootstrap replicates for SNaQ (Hmax=0,1);
  processors = $processors, Number of processors to run SNaQ (and bootsnaq);
  Server for running the script = $host_name.stat.wisc.edu

  Running time shown below: 
  """

argument_file = joinpath(outfolder, "arguments-$paramname_root.log")

open(argument_file, "a") do io
    println(io, SNAQ_arguments)  # Append SNAQ arguments
    show(io, to)  # Append timer output
end
