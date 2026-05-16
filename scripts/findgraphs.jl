# ============================================================================
# scripts/findgraphs.jl
#
# Purpose : Entry point for find_graphs (admixture-graph) inference across
#           replicates of one parameter setting. For each replicate it
#           launches findgraphs_1rep.R, which (1) extracts SNPs via snp-sites,
#           (2) converts to eigenstrat, and (3) runs find_graphs() at h=0
#           and h=1 from the admixtools R package.
# Inputs  : output/<paramname>/rep<id>/iqtreefolder/concatenated.fasta
#           (and per-rep ancillary files written by simulation.jl)
# Outputs : output/<paramname>/rep<id>/findgraph/rep<id>_admix{0,1}_*.{rds,txt}
#           output/<paramname>/rep<id>/findgraph/rep<id>.vcf, eigenstrat_*.*
# Usage   : julia -p 100 --project=. scripts/findgraphs.jl \
#               --dup_rate 0.0003 --loss_rate 0.0003 \
#               --ratevar G --n_reps 100 --n_inds 1 --runs 100 --block 1000
# Note    : Steps 2 (SNaQ) and 3 (find_graphs) are independent and can run in
#           parallel after Step 1 (simulation.jl).
# ============================================================================

using ArgParse
using TimerOutputs 
using Distributed
using Dates, TimeZones 
using PhyloNetworks
using RCall

# Set R environment for all workers to avoid R_HOME conflicts
# This is important when using RCall in a distributed setting 
@everywhere begin
    # Get R_HOME if not already set
    if !haskey(ENV, "R_HOME") || ENV["R_HOME"] == ""
        try
            r_home = strip(read(`R RHOME`, String))
            ENV["R_HOME"] = r_home
        catch
            @warn "Could not detect R_HOME. Please set it manually."
        end
    end
end

@everywhere using Printf 
@everywhere using PhyloNetworks
@everywhere using RCall
@everywhere using DataFrames, CSV 
@everywhere include("utilities.jl")
@everywhere include("../third_party_scripts/interop_admixtools.jl")

const to = TimerOutput()  
tz = TimeZone("America/Chicago")
current_time_tz = ZonedDateTime(now(), tz) 
time = Dates.format(current_time_tz, "yyyy-mm-dd HH:MM:SS zzz") 

function parse_commandline()
    s = ArgParseSettings() 
    @add_arg_table s begin 
      
      # Specify which rep and which parameter sets to run SnaQ: 
      "--dup_rate"
        help = "Paramater setting (duplication rate) to run findgraph"
        arg_type = Float64
        required = true
      "--loss_rate"
        help = "Paramater setting (gene loss rate) to run findgraph"
        arg_type = Float64
        required = true
      "--ratevar"
        help ="Parameter setting (variation rate) to run findgraph"
        arg_type = String
        required = true
      "--n_reps" 
        help = "How many reps this parameter set was run in total"
        arg_type = Int
        default = 100 
      "--rep_start"
        help = "Parameter setting (start of replicates) to run findgraph"
        arg_type = Int
        default = 1 # default is 1 if not specified
      "--rep_end"
        help = "Parameter setting (end of replicates) to run findgraph"
        arg_type = Int
        default = -1 # temporary default. 
      "--n_inds"
        help = "Paramater setting (gene loss rate) to run findgraph"
        arg_type = Int
        default = 1
      "--processors"
        help= "Number of processors to multi-process findgraph"
        arg_type = Int
        default = nprocs() - 1 #-p will naturally add 1 to this max processors 
      "--SF"
        help = "Scaling factor to scale effective population Ne
        (Default = 1.0, no scaling)"
        arg_type = Float64
        default = 1.0
      "--debug_mode" 
        help = "Debugging mode: Keep all temporary folders"
        default = false
        arg_type = Bool
      
      # Specify arguments for findgraph: 
      "--runs"
        help = "Number of runs in findgraph"
        arg_type = Int
        default = 100
      "--stop_gen" 
        help = "Number of generation to stop, default in findgraphs = 100"
        arg_type = Int 
        default = 100
      "--outgroup"
        help = "Outgroup population name used in findgraphs (default: A)"
        arg_type = String
        default = "A" # Homo, see speciestree.jl 
      "--block"
        help = "Number of blocks (approximate) used in find_graphs" 
        arg_type = Int
        default = 1000 # we typically simulate genes with 1000 bp 
      "--selection_method"
        help = "Method to select graphs"
        arg_type = String
        default = "select_threshold_within_ll" # always use this 
        # the other option is a legacy option 
      "--threshold"
        help = "Threshold for graph selection, default: NULL"
        arg_type = Int
        default = 2 # corresponds to 4 AIC

      # Other arguments: 
      "--log"
        help = "Enable logging: Save all screen output to worker-specific logs"
        default = true
        arg_type = Bool
      "--vcf_transfer_scripts" 
        # Multiallelic is a legacy option. 
        help = "Script to transfer vcf to eigenstrat.\n
                Biallelic: standard, no missing alleles.\n
                Multiallelic: legacy, treats missing alleles."
        arg_type = String 
        default = "biallelic"
      "--gene_len" 
        help = "The length of simulated gene sequences (Default = 1000 bp)" 
        arg_type = Int 
        default = 1000  

    end 
    
    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1 # If rep_end not specified, then n_reps
        parsed_args["rep_end"] = parsed_args["n_reps"]
    end

    return parsed_args
  end

#-----------------------------------------------#       
#               Initialization
#-----------------------------------------------# 
parsed_args = parse_commandline()

#-------------- Parse arguments #--------------#
dup_rate = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"]
n_reps = parsed_args["n_reps"]
n_inds = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
processors = parsed_args["processors"]
# Find graph specific arguments: 
runs = parsed_args["runs"]
stop_gen = parsed_args["stop_gen"]
outgroup = parsed_args["outgroup"]
block = parsed_args["block"] 
selection_method = parsed_args["selection_method"]  
SF = parsed_args["SF"]
gene_len = parsed_args["gene_len"] 
debug_mode = parsed_args["debug_mode"]
log_output = parsed_args["log"]
vcf_transfer_scripts = parsed_args["vcf_transfer_scripts"] 


#= This threshold is used to select graphs within Δ LL of the best graph
Under the same complexity (same num_admix), 
LRT = -2log(L0/L1) ~ χ² with df = difference in number of parameters 
Because complexity difference = 0 when under the same num_admix 
LRT = 2 * Δ LL which is the same as AIC or BIC difference 
Because differences of Δ ≤ 2 in AIC/BIC are widely regarded as 
indicating that models have essentially the same level of empirical support, 
we use default threshold = 2 here (although the threshold can be configured).
=# 
threshold = parsed_args["threshold"]

# if selection_method not in ["both", ""]
#     error("Invalid selection method. Must be 'both' or empty string.")
# end

# output suffix -> hard-code to avoid too many parameters 
output_graph_suffix = "_unique_graphs.rds" 
output_f2_suffix = "_f2.rds" 
output_summary_table_suffix = "_summary_table.txt"

# Define true species tree in Newick format
true_tree_newick = "(A,((((B,C),(D,E)),F),(G,H)));" 

#--------------- Set up folders ----------------# 
paramname_root = set_up_paramname_root(dup_rate, loss_rate, ratevar, 
                                      n_inds, SF, gene_len)  
outfolder = "output/$paramname_root"

# Ensure outfolder exists for log files
mkpath(outfolder) 

#--------------- Set up seeds -------------------#
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
# unique master seed for each parameter setting: 
master_seed = generate_master_seed(params_dict_for_seed_setting) 

software_names = ["findgraphs", "qpgraph"] 
seed_dic = generate_software_seeds(master_seed, software_names)
# This seed used to generate m (n_reps) x 2 seed array: 
seed_findgraphs = seed_dic["findgraphs"] 
# This generates seeds for qpgraph: 
seed_qpgraph = seed_dic["qpgraph"] 

# Set up seed arrays for findgraphs: n_reps x runs matrix  
# Each replicate gets 'runs' number of seeds for independent randomization
seed_array_findgraphs = seed_generator(
                        seed_findgraphs, 
                        n_reps, # number of replicates (rows)
                        runs,   # number of runs per replicate (columns) 
                        outfolder, 
                        "findgraph_seed.txt") 

# Keep the old 2-column format for qpgraph compatibility
# seed_array_qpgraph = seed_generator(
#                       seed_qpgraph, 
#                       n_reps, 
#                       2, 
#                       outfolder, 
#                       "random_seed_qpgraph.txt")

seed_array_findgraphs_path = joinpath(outfolder, "findgraph_seed.txt") 
# seed_array_qpgraph_path = joinpath(outfolder, "random_seed_qpgraph.txt")

#-----------------------------------------------#       
#    Check pre-existing files and remove them
#-----------------------------------------------# 
folder_path_list = []
for simulation_rep in rep_start:rep_end
  rep_number_string = pad_number(simulation_rep, n_reps)
  rep_folder_path = joinpath(outfolder, "rep$rep_number_string") 
  push!(folder_path_list, rep_folder_path) 
end 

findgraph_folderlist = []
index_length = rep_end - rep_start + 1 # To match rep_start:rep_end 
for ind in 1:index_length 
  fingraph_folder = setup_rep_output_folders(folder_path_list, ind, "findgraph")
  push!(findgraph_folderlist, fingraph_folder)
end 

check_existing_dir(findgraph_folderlist) # see utilies.jl. 

#-----------------------------------------------#       
#  Push global params to all processors
#-----------------------------------------------#
@everywhere global index_length = $index_length
@everywhere global n_reps = $n_reps
@everywhere global rep_start = $rep_start
@everywhere global rep_end = $rep_end
@everywhere global findgraph_folderlist = $findgraph_folderlist 
@everywhere global folder_path_list = $folder_path_list
@everywhere global runs = $runs
@everywhere global stop_gen = $stop_gen
@everywhere global block = $block
@everywhere global outgroup = $outgroup
@everywhere global seed_array_findgraphs_path = $seed_array_findgraphs_path
# @everywhere global seed_array_qpgraph_path = $seed_array_qpgraph_path
@everywhere global output_graph_suffix = $output_graph_suffix 
@everywhere global output_f2_suffix = $output_f2_suffix
@everywhere global output_summary_table_suffix = $output_summary_table_suffix
@everywhere global true_tree_newick = $true_tree_newick
@everywhere global selection_method = $selection_method
@everywhere global threshold = $threshold
@everywhere global log_output = $log_output
@everywhere global outfolder = $outfolder
@everywhere global paramname_root = $paramname_root
@everywhere global vcf_transfer_scripts = $vcf_transfer_scripts 
@everywhere global debug_mode = $debug_mode

# Initialize worker-specific logging variable
@everywhere global worker_logfile = nothing

#-----------------------------------------------#       
#         Variant Calling + File Conversion
#    Goal 1: snp-sites calls SNP from fasta files 
#   Goal 2:  Convert VCF to eigenstrat using codes 
#           from https://github.com/mathii/gdc 
# Goal 3: Modify the .ind files from vcf2eigenstrat_modified_py3.py 
#               so each sample gets unique ID
#-----------------------------------------------# 

@everywhere function callVariant_convertVCF_1rep(
            simulation_rep::Int, 
            findgraph_folderlist::Vector, 
            folder_path_list::Vector, 
            rep_start::Int, 
            n_reps::Int, 
            block::Int,
            vcf_transfer_scripts::String = "biallelic"
            ) 
    
    global debug_mode
    # concatenated fasta files are stored in Rep$id/seqgenfolder
    seqgenfolder = setup_rep_output_folders(folder_path_list, 
                                          simulation_rep, 
                                          "seqgenfolder")

    findgraph_folder = findgraph_folderlist[simulation_rep]
    mkpath(findgraph_folder)
   
    # Goal 1: snp-sites calls SNP from fasta files 
    match_result = match(r"rep(\d+)/", seqgenfolder) 
    # Below: Find the file name in the concatenated fasta file  
    rep_id = match_result.captures[1] 
    
    # This could also be listed as simulation_rep + rep_start - 1:
    # Having a debug-part to check if the rep_id is correct: 
    rep_id_check = pad_number(simulation_rep + rep_start - 1, n_reps)
    if rep_id != rep_id_check # Exit if rep_id doesn't match
      println("rep_id {$rep_id} and {$rep_id_check}")
      error("Check rep ID!")
    end # Only for debugging and double checking 

    fasta_file_name = "concate_alignment_rep$(rep_id).fasta" 
    fasta_file = joinpath(seqgenfolder, fasta_file_name) 
    vcf_file = joinpath(findgraph_folder, "rep$(rep_id).vcf")

    # Analyze alignment length using utility function
    alignment_length = get_alignment_length(fasta_file)
    worker_println("Rep$rep_id: Alignment length is $alignment_length bp")
    
    # estimate the blgsize based on alignment length 
    # blgsize should be the size of each individual gene locus 
    # block is the number of blocks (approximate) used in find_graphs 
    # blgsize = alignment_length / block = average block size
    blgsize = Int(floor(alignment_length / block)) 

    if debug_mode 
      println("="^60)
      println("DEBUG INFO for Rep$rep_id:")
      println("Alignment length: $alignment_length bp")
      println("Block parameter: $block") 
      println("$alignment_length / $block = $(alignment_length / block)")
      println("Calculated blgsize: $blgsize")
      println("="^60)
    end

    # run snp-sites: 
    run(`./executables/snp-sites $fasta_file -v -o $vcf_file`)
    
    # Convert VCF to eigenstrat (github.com/mathii/gdc)
    eigenstrat_file = joinpath(findgraph_folder, "eigenstrat_rep$(rep_id)") 

    if vcf_transfer_scripts == "biallelic"
      # Use the biallelic version which doesn't treat missing alleles 
      run(`python third_party_scripts/vcf2eigenstrat_biallelic.py \
        -v $(vcf_file) \
        -o $(eigenstrat_file)`) 
    elseif vcf_transfer_scripts == "multiallelic" 
      #= multiallelic: legacy, treats missing alleles (unused by default) =#
      run(`python third_party_scripts/vcf2eigenstrat_multiallelic.py \
        -v $(vcf_file) \
        -o $(eigenstrat_file)`)
    else
      error("Invalid vcf_transfer_scripts: " *
          "must be 'biallelic' or 'multiallelic'.")
    end

    # Goals 3: output .ind files is re-assigned a unique pop to each taxon
    eigenstrat_ind_file = "$(eigenstrat_file).ind"
    modified_ind_file = "$(eigenstrat_file)_modified.ind"
    replace_pop_inInd_file(eigenstrat_ind_file, modified_ind_file) 
    # f2_from_geno requires a single prefix; rename in place
    # The below process is not done in a function to avoid over-writing 
    # Choose not to use indAsPop from vcf2eigenstrat_modified_py3.py, 
    # since we might have multiple individual from one pop (eg. A_0, A_1) 
    run(`sh -c "rm $eigenstrat_ind_file && \
        mv $modified_ind_file $eigenstrat_ind_file"`)
    worker_println("Rep$rep_id: vcf → eigenstrat conversion done")

    return blgsize # return the estimated blgsize 
end 

#-----------------------------------------------#       
#         Run findgraphs
#    Goal 1: Run findgraphs with runs 
#   Goal 2: Merge all graphs from a single runs
#   Goal 3: Calculate wr and make a summary table for each run 
#-----------------------------------------------# 

@everywhere function run_findgraphs_1rep(
                      simulation_rep::Int, 
                      findgraph_folderlist::Vector, 
                      rep_start::Int, n_reps::Int, 
                      num_admix::Int, stop_gen::Int, 
                      outgroup::String, runs::Int, 
                      blgsize::Int, seed_array_findgraphs_path::String, 
                      output_graph_suffix::String, 
                      output_f2_suffix::String, 
                      output_summary_table_suffix::String,
                      true_tree_newick::String
                      ) 
  
  findgraph_folder = findgraph_folderlist[simulation_rep]
  rep_id = pad_number(simulation_rep + rep_start - 1, n_reps)

  prefix = "eigenstrat_rep$rep_id"

  # Construct command arguments as a vector to avoid shell escaping issues
  cmd_args = [
    "Rscript", "./scripts/findgraphs_1rep.R",
    "-i", findgraph_folder,
    "-o", findgraph_folder,
    "-p", prefix,
    "-k", string(num_admix),
    "--stop_gen", string(stop_gen),
    "--outgroup", outgroup,
    "--runs", string(runs),
    "-b", string(blgsize),
    "-r", rep_id,
    "-s", seed_array_findgraphs_path,
    "--output_graph_suffix", output_graph_suffix,
    "--output_f2_suffix", output_f2_suffix,
    "--output_summary_table_suffix", output_summary_table_suffix,
    "--true_tree_newick", true_tree_newick,
    "--selection_method", selection_method,
    "--threshold", string(threshold)
  ]
  # All things passed to Rscript should be in a vector 
  run(Cmd(cmd_args)) 
  
end 

#-----------------------------------------------#       
#         Summarize results  
#-----------------------------------------------# 
function summarize_findgraphs_results(
          simulation_rep::Int, 
          findgraph_folderlist::Vector, 
          rep_start::Int, n_reps::Int, 
          output_summary_table_suffix::String
          )
  
  findgraph_folder = findgraph_folderlist[simulation_rep]
  rep_id = pad_number(simulation_rep + rep_start - 1, n_reps)

  tables = Dict(
      0 => "rep$(rep_id)_admix0$output_summary_table_suffix",
      1 => "rep$(rep_id)_admix1$output_summary_table_suffix"
  )

  wr_found = Dict(0 => false, 1 => false)

  for k in sort(collect(keys(tables)))
      table_path = joinpath(findgraph_folder, tables[k])
      df = CSV.read(table_path, DataFrame; delim=' ', 
                    ignorerepeated=true, quotechar='"')
      valid_graphs = df[df.WR_smaller_than_3 .== true, :]

      if nrow(valid_graphs) > 0
          wr_found[k] = true
      end
  end

  best_k = if wr_found[0]
      0
  elseif wr_found[1]
      1
  else
      ">1"
  end

  return Dict(
      "rep_id" => rep_id,
      "best_k" => best_k
  )
end

#-----------------------------------------------#       
#           Run and time the process
#-----------------------------------------------# 
function main()
  
  # Initialize worker logs if logging is enabled
  if log_output
    @everywhere open_worker_log()
  end

  # Store blgsize for each replicate
  blgsize_list = Int[]
  
  @timeit to "Variant Calling + Convert VCF to eigenstrat" begin 
    blgsize_list = pmap(simulation_rep -> callVariant_convertVCF_1rep(
                            simulation_rep, 
                            findgraph_folderlist, 
                            folder_path_list, 
                            rep_start, 
                            n_reps, 
                            block 
                            ), 1:index_length)
  end 

  @timeit to (
    "Run findgraphs with num_admix=0 for $runs runs "*
    "from rep$rep_start to rep$rep_end"
  ) begin 
    num_admix = 0
    pmap(simulation_rep -> run_findgraphs_1rep(
                          simulation_rep, 
                          findgraph_folderlist, 
                          rep_start, n_reps, 
                          num_admix, stop_gen, outgroup, 
                          runs, 
                          blgsize_list[simulation_rep], 
                          seed_array_findgraphs_path, 
                          output_graph_suffix, 
                          output_f2_suffix, 
                          output_summary_table_suffix,
                          true_tree_newick
                          ), 1:index_length)
  end 

  @timeit to (
    "Run findgraphs with num_admix=1 for $runs runs " *
    "from rep$rep_start to rep$rep_end"
  ) begin
    num_admix = 1
    pmap(simulation_rep -> run_findgraphs_1rep(
                          simulation_rep, 
                          findgraph_folderlist, 
                          rep_start, n_reps, 
                          num_admix, stop_gen, 
                          outgroup, runs, blgsize_list[simulation_rep], 
                          seed_array_findgraphs_path, 
                          output_graph_suffix, 
                          output_f2_suffix, 
                          output_summary_table_suffix,
                          true_tree_newick
                          ), 1:index_length)
  end

  @timeit to "Organize all results" begin

    results = DataFrame(rep_id = String[], best_k = Any[]) 

    for simulation_rep in 1:index_length 
        result = summarize_findgraphs_results(
                  simulation_rep, 
                  findgraph_folderlist, 
                  rep_start, 
                  n_reps, 
                  output_summary_table_suffix) 
        push!(results, (rep_id = result["rep_id"], best_k = result["best_k"]))
    end 
    output_file = joinpath(outfolder, "findgraphs_summary_results.csv")
    CSV.write(output_file, results)
  end
  
  # Close worker logs if logging is enabled
  if log_output
    @everywhere close_worker_log()
  end
end


if abspath(PROGRAM_FILE) == @__FILE__
  main()

  n_processors = nprocs() 
  host_name = gethostname()  

  # save all parameters and write into the output dir
  arguments = """
  #=====================================================#
  #-----------------AdmixtureTools2---------------------#
  #=====================================================#
  ---Arguments used for simulating this output dataset---
  master seed = $master_seed, unique master_seed for this parameter setting 
  Seed for findgraph= $seed_findgraphs (rep x 2 seed array)
  This experiment runs from rep$rep_start to rep$rep_end
  ---'findgraphs' specific patameters---
  Number of runs to run findgraphs = $runs;
  Number of blocks used in findgraphs = $block;
  Number of generations to stop running findgraphs = $stop_gen;
  ---Other information---
  Number of processors used = $n_processors;
  Server for running the script = $host_name.stat.wisc.edu;
  Time of running the script = $time;
  ---Running time---
  """
  argument_file = joinpath(outfolder, "arguments-$paramname_root.log")

  if isfile(argument_file)
      open(argument_file, "a") do io
          show(io, to)
      end
  else
      write(argument_file, arguments)
      open(argument_file, "a") do io
          show(io, to)
      end
  end
end

println("===============================================") 
println("Findgraphs analysis completed for parameter setting: $paramname_root") 
println("Output folder: $outfolder")

