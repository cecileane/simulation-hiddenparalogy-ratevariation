#= 
Run PhyNEST on all the replicates of one parameter setting.
For each replicate, this calls phynest_1rep.jl, which runs PhyNEST with
hmax=0 (tree) and hmax=1 (one reticulation) on the concatenated alignment,
starting from the ASTRAL tree. Model selection between H=0 and H=1 uses
T = score(H0) - score(H1) (scores are negative log composite likelihoods),
compared to a threshold calibrated on the null setting (no dup/loss,
no rate variation), like the recalibrated WR threshold of find_graphs.
Input, in output/<paramname>/rep<id>/: seqgenfolder/concate_alignment_*.fasta
  and astralfolder/astral.tre
Output, in output/<paramname>/rep<id>/phynestfolder/:
  H0_output/H0_hc.{log,out}, H1_output/H1_hc.{log,out}, phynest_results.csv
  and for the whole setting: output/<paramname>/phynest_summary_results.csv
Usage: julia -p 100 scripts/phynest.jl --dup_rate 0.0003 --loss_rate 0.0003 \
       --ratevar G --n_reps 100 --n_inds 1 --runs 100
--rep_start / --rep_end to run or resume part of the replicates,
the seeds do not depend on that.
PhyNEST v0.1.x needs PhyloNetworks < 1.0, so phynest_1rep.jl runs in its
own environment, envs/phynest (in the repo). The worker activates it and
runs Pkg.instantiate() itself. To build it by hand:
  julia --project=envs/phynest -e 'using Pkg; Pkg.instantiate()'
PhyNEST, SNaQ and find_graphs are independent, all three only need the
output of simulation.jl 
=# 

using ArgParse
using TimerOutputs
using Dates
using TimeZones
using CSV
using DataFrames
using Statistics
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

      # Specify which rep and which parameter sets to run PhyNEST:
      "--dup_rate"
        help = "Paramater setting (duplication rate) to run phynest"
        arg_type = Float64
        required = true
      "--loss_rate"
        help = "Paramater setting (gene loss rate) to run phynest"
        arg_type = Float64
        required = true
      "--ratevar"
        help = "Parameter setting (variation rate) to run phynest"
        arg_type = String
        required = true
      "--n_reps"
        help = "How many reps this parameter set was run in total"
        arg_type = Int
        default = 100
      "--rep_start"
        help = "Start of the replicate range to run phynest"
        arg_type = Int
        default = 1 # default is 1 if not specified
      "--rep_end"
        help = "Parameter setting (end of the range of replicates) to run"
        arg_type = Int
        default = -1 # temporary default.
      "--n_inds"
        help = "Number of individuals per species"
        arg_type = Int
        default = 1
      "--SF"
        help = "Scaling factor to scale effective population Ne
          (Default = 1.0, no scaling)"
        arg_type = Float64
        default = 1.0
      "--gene_len"
        help = "The length of simulated gene sequences (Default = 1000 bp)"
        arg_type = Int
        default = 1000

      # Specify arguments for PhyNEST:
      "--runs"
        help = "Number of independent runs per PhyNEST search
          (default = 100; PhyNEST's own default is 10). The same
          number is used for Hmax=0 and Hmax=1 so the composite
          likelihood comparison stays fair; 100 runs saturate both
          searches (P(all starts NNI-perturbed) = 0.75^100)."
        arg_type = Int
        default = 100
      "--max_steps"
        help = "Max steps per PhyNEST search (default = 250000)"
        arg_type = Int
        default = 250000
      "--tree"
        help = "Species tree used in the simulation: reptile or plant"
        arg_type = String
        default = "reptile"
      "--outgroup"
        help = "Outgroup taxon name used in PhyNEST (default: A)"
        arg_type = String
        default = "A" # Homo, see speciestree_reptile.jl
    end

    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1 # default: run all reps
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
max_steps = parsed_args["max_steps"]
outgroup = parsed_args["outgroup"]
SF = parsed_args["SF"]
tree = parsed_args["tree"]
gene_len = parsed_args["gene_len"]

# Define true species tree in Newick format (reptile or plant, as simulated)
true_tree_newick = (startswith(tree, "plant") ?
    "(A,(B,((C,D),(E,(H,(F,G))))));" :
    "(A,((((B,C),(D,E)),F),(G,H)));")

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

phynestfolder_list = []
index_length = rep_end - rep_start + 1 # match rep_start:rep_end
for ind in 1:index_length
  phynestfolder = setup_rep_output_folders(
      folder_path_list, ind, "phynestfolder")
  push!(phynestfolder_list, phynestfolder)
end
check_existing_dir(phynestfolder_list) # see utilities.jl.

#-----------------------------------------------#
#       Estimate Networks using PhyNEST
#-----------------------------------------------#
#--------------- Set up seeds -------------------#
#= Seed generation:
# 1. Generate a master seed for each parameter setting using paramname_root
# 2. Generate a seed for each software using the master seed
# 3. Generate a n_rep x 2 seed array for PhyNEST using the PhyNEST seed
# 4. For each replicate, select the row in seed_array matching the rep index
=#
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
# unique master seed for each parameter setting:
master_seed = generate_master_seed(params_dict_for_seed_setting)

software_names = ["phynest"]
seed_dic = generate_software_seeds(master_seed, software_names)
# This seed used to generate m (n_reps) x 2 seed array:
seed_phynest = seed_dic["phynest"]

# n_reps x 2 seed array:
# 1st seed -> infer tree Hmax = 0
# 2nd seed -> infer network Hmax = 1
# phyne! has no seed argument, so phynest_1rep.jl seeds the global RNG
# of its own Julia process before each search.
seed_array = seed_generator(
    seed_phynest, n_reps, 2, outfolder, "random_seed_phynest.txt")
# For each rep, seed is selected as [repID, i for i in 1:2]

#-----------------------------------------------#
#       Run PhyNEST for each replicate
#-----------------------------------------------#
# broadcast variables to all processors:
@everywhere seed_array = $seed_array
@everywhere phynestfolder_list = $phynestfolder_list
@everywhere folder_path_list = $folder_path_list
@everywhere outfolder = $outfolder
@everywhere paramname_root = $paramname_root
@everywhere runs = $runs
@everywhere max_steps = $max_steps
@everywhere outgroup = $outgroup
@everywhere n_inds = $n_inds
@everywhere n_reps = $n_reps
@everywhere rep_start = $rep_start
@everywhere rep_end = $rep_end
@everywhere true_tree_newick = $true_tree_newick

@everywhere begin
"""
Run phynest_1rep.jl for one replicate, in the envs/phynest environment
Input: ind (Int): index of the replicate in phynestfolder_list (1 = rep_start)
Notes: uses the globals sent to the workers above (seed_array,
    phynestfolder_list, folder_path_list, runs, ...)
"""
  function run_phynest_for_replicate(ind)

    seqgenfolder = setup_rep_output_folders(
        folder_path_list, ind, "seqgenfolder")
    astralfolder = setup_rep_output_folders(
        folder_path_list, ind, "astralfolder")
    phynestfolder = phynestfolder_list[ind]

    # Identify the correct rep to start with
    ind_in_seed_array = ind + rep_start - 1
    rep_id = pad_number(ind_in_seed_array, n_reps)

    seed_net0 = seed_array[ind_in_seed_array, 1]
    seed_net1 = seed_array[ind_in_seed_array, 2]
    # seed selection does not depend on rep_start / rep_end,
    # see snaq.jl for an example of the indexing

    run(`julia --project=envs/phynest ./scripts/phynest_1rep.jl \
        --seqgenfolder $seqgenfolder \
        --astralfolder $astralfolder \
        --phynestfolder $phynestfolder \
        --rep_id $rep_id \
        --seed_net0 $seed_net0 \
        --seed_net1 $seed_net1 \
        --runs $runs \
        --n_inds $n_inds \
        --max_steps $max_steps \
        --outgroup $outgroup \
        --true_tree_newick $true_tree_newick`)
  end
end

@timeit to "Running PhyNEST from rep$rep_start to rep$rep_end" begin

  pmap(ind -> begin
      println("Worker $(myid()): Starting task $ind")
      run_phynest_for_replicate(ind)
  end, 1:index_length)
  # index_length = rep_end - rep_start + 1
  # pmap will automatically distribute

end

#----------------------------------------------------#
#       Organize all results
#----------------------------------------------------#
"""
Concatenate the phynest_results.csv of all replicates into one table
Input: phynestfolder_list (Vector), outfolder (String)
Output: writes <outfolder>/phynest_summary_results.csv and returns the table
Notes: also prints the 95th percentile of T = score(H0) - score(H1) in this
    setting. On the null setting this is the calibrated threshold T*, on the
    other settings compare with the null T* to get the false positive rate
"""
function summarize_phynest_results(phynestfolder_list::Vector,
        outfolder::String)
    results = DataFrame()
    for phynestfolder in phynestfolder_list
        result_path = joinpath(phynestfolder, "phynest_results.csv")
        if isfile(result_path)
            df = CSV.read(result_path, DataFrame)
            results = vcat(results, df)
        else
            @warn "Missing PhyNEST results file: $result_path"
        end
    end

    output_file = joinpath(outfolder, "phynest_summary_results.csv")
    CSV.write(output_file, results)
    println("PhyNEST summary saved to $output_file")

    if nrow(results) > 0
        T95 = quantile(results.T, 0.95)
        println("Within-setting 95th percentile of T: $T95")
    end
    return results
end

@timeit to "Organize all results" begin
    summarize_phynest_results(phynestfolder_list, outfolder)
end

#-----------------------------------------------#
#       Ouput the running time into .log
#-----------------------------------------------#
host_name = gethostname()

# Define the number of processors based on the Julia parallel environment
processors = nprocs()

PHYNEST_arguments = """
  #=====================================================#
  #----------------------PhyNEST------------------------#
  #=====================================================#
  ---Arguments used to run PhyNEST---
  runs = $runs, independent runs per search (Hmax=0 and Hmax=1);
  max_steps = $max_steps, max steps per hill-climbing search;
  outgroup = $outgroup, outgroup taxon used to root networks;
  seed_phynest = $seed_phynest, master seed for the PhyNEST seed array;
  --- Other Information ---
  processors = $processors, Number of processors to run PhyNEST;
  Server for running the script = $host_name.stat.wisc.edu
  Time of running the script = $time;
  --- Running time ---
  """

argument_file = joinpath(outfolder, "arguments-$paramname_root.log")

open(argument_file, "a") do io
    println(io, PHYNEST_arguments)  # Append PhyNEST arguments
    show(io, to)  # Append timer output
end

println("=============================================")
println("PhyNEST Analysis Completed!")
println("PhyNEST arguments have been saved to $argument_file")
