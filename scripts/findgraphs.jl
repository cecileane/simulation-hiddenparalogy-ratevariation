# This script documents using find_graphs to estimate graph 

using ArgParse 
using TimerOutputs 
using Distributed
using CSV, DataFrames
@everywhere include("utilities.jl")

const to = TimerOutput()  

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
        required = true
      "--rep_start"
        help = "Parameter setting (start of the range of replicates) to run findgraph"
        arg_type = Int
        default = 1 # default is 1 if not specified
      "--rep_end"
        help = "Parameter setting (end of the range of replicates) to run findgraph"
        arg_type = Int
        default = -1 # temporary default. 
      "--n_inds"
        help = "Paramater setting (gene loss rate) to run findgraph"
        arg_type = Int
        default = 1
      "--processors"
        help= "Number of processors to multi-process findgraph"
        arg_type = Int
        default = nprocs() - 1 # -p will naturally add 1 to this maximum processors 
      
      # Specify arguments for findgraph: 
      "--runs"
        help = "Number of runs in findgraph"
        arg_type = Int
        default = 100
      "--stop_gen" 
        help = "Number of generation to stop, default in findgraphs = 100"
        arg_type = Int 
        default = 100
      "--blgsize"
        help = "Block size used in findgraphs"
        arg_type = Int
        default = 200
      "--outgroup"
        help = "Outgroup population name used in findgraphs, default is homo sapiens"
        arg_type = String
        default = "A"
    end 
    
    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1 # If rep_end is not specified, set it to n_reps
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
blgsize = parsed_args["blgsize"]
outgroup = parsed_args["outgroup"]

# output suffix -> hard-code to avoid too many parameters 
output_graph_suffix = "_unique_graphs.rds"
output_f2_suffix = "_f2.rds" 
output_summary_table_suffix = "_summary_table.txt" 

#--------------- Set up folders ----------------# 
paramname_root = "DUP$dup_rate-LOS$loss_rate-RV$ratevar-N_ind$n_inds" # Specify to find the folder
outfolder = "output/$paramname_root" 

#--------------- Set up seeds -------------------#
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
# unique master seed for each parameter setting: 
master_seed = generate_master_seed(params_dict_for_seed_setting) 

software_names = ["findgraphs", "qpgraph"] 
seed_dic = generate_software_seeds(master_seed, software_names) # see utility.jl 
seed_findgraphs = seed_dic["findgraphs"] # This seed used to generate m (n_reps) x 2 seed array
seed_qpgraph = seed_dic["qpgraph"] 

# Set up seed arrays for num_admix = 0 and 1 
seed_array_findgraphs = seed_generator(seed_findgraphs, n_reps, 2, outfolder, "random_seed_findgraphs.txt") # here, use n_reps since in findgraph_1rep seed is selected based on the actual simulation_rep
seed_arrary_qpgraph = seed_generator(seed_qpgraph, n_reps, 2, outfolder, "random_seed_qpgraph.txt")

seed_array_findgraphs_path = joinpath(outfolder, "random_seed_findgraphs.txt") 
seed_array_qpgraph_path = joinpath(outfolder, "random_seed_qpgraph.txt")

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
index_length = rep_end - rep_start + 1 # see utilities.jl to match how we set up rep_start:rep_end, 
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
@everywhere global blgsize = $blgsize
@everywhere global outgroup = $outgroup
@everywhere global seed_array_findgraphs_path = $seed_array_findgraphs_path
@everywhere global seed_array_qpgraph_path = $seed_array_qpgraph_path
@everywhere global output_graph_suffix = $output_graph_suffix 
@everywhere global output_f2_suffix = $output_f2_suffix
@everywhere global output_summary_table_suffix = $output_summary_table_suffix

#-----------------------------------------------#       
#         Variant Calling + File Conversion
#    Goal 1: snp-sites calls SNP from fasta files 
#   Goal 2:  Convert VCF to eigenstrat using codes from https://github.com/mathii/gdc 
# Goal 3: Modify the .ind files from vcf2eigenstrat_modified_py3.py so each sample gets unique ID
#-----------------------------------------------# 

@everywhere function callVariant_convertVCF_1rep(simulation_rep::Int, findgraph_folderlist::Vector, folder_path_list::Vector, rep_start::Int, n_reps::Int) 
    # concatenated fasta files are stored in Rep$id/seqgenfolder
    seqgenfolder = setup_rep_output_folders(folder_path_list, simulation_rep, "seqgenfolder")
    findgraph_folder = findgraph_folderlist[simulation_rep]
    mkpath(findgraph_folder)
   
    # Goal 1: snp-sites calls SNP from fasta files 
    match_result = match(r"rep(\d+)/", seqgenfolder) 
    rep_id = match_result.captures[1]  # Find the file name in the concatenated fasta file: 
    
    # This could also be listed as simulation_rep + rep_start - 1:
    # Having a debug-part to check if the rep_id is correct: 
    rep_id_check = pad_number(simulation_rep + rep_start - 1, n_reps)
    if rep_id != rep_id_check # Exit if rep_id doesn't match
      println("rep_id {$rep_id} and {$rep_id_check}")
      error("Check rep ID!")
    end # This part is useless but I want to keep it as a debugging function and double check. 

    fasta_file_name = "concate_alignment_rep$(rep_id).fasta" 
    fasta_file = joinpath(seqgenfolder, fasta_file_name) 
    vcf_file = joinpath(findgraph_folder, "rep$(rep_id).vcf")
    run(`./executables/snp-sites $fasta_file -v -o $vcf_file`) # run snp-sites
    
    # Goal 2:  Convert VCF to eigenstrat using codes from https://github.com/mathii/gdc  
    eigenstrat_file = joinpath(findgraph_folder, "eigenstrat_rep$(rep_id)") 
    run(`python scripts/vcf2eigenstrat_modified_py3.py -v $(vcf_file) -o $(eigenstrat_file)`)

    # Goals 3: output .ind files is re-assigned a unique pop to each taxon
    eigenstrat_ind_file = "$(eigenstrat_file).ind"
    modified_ind_file = "$(eigenstrat_file)_modified.ind"
    replace_pop_inInd_file(eigenstrat_ind_file, modified_ind_file) 
    # rm and rename modified_ind_file, because f2_from_geno can only recognize one prefix 
    # The below process is not done in a function to avoid over-writing 
    # Choose not to use indAsPop from vcf2eigenstrat_modified_py3.py, 
    # since we might have multiple individual from one pop (eg. A_0, A_1) 
    run(`sh -c "rm $eigenstrat_ind_file && mv $modified_ind_file $eigenstrat_ind_file"`)
    println("Finish converting .vcf to eigenstrate and adding pop names")
end 

#-----------------------------------------------#       
#         Run findgraphs
#    Goal 1: Run findgraphs with runs 
#   Goal 2: Merge all graphs from a single runs
#   Goal 3: Calculate wr and make a summary table for each run 
#-----------------------------------------------# 

@everywhere function run_findgraphs_1rep(simulation_rep::Int, findgraph_folderlist::Vector, rep_start::Int, n_reps::Int, num_admix::Int, stop_gen::Int, outgroup::String, runs::Int, blgsize::Int, seed_array_findgraphs_path::String, output_graph_suffix::String, output_f2_suffix::String, output_summary_table_suffix::String) 
  
  findgraph_folder = findgraph_folderlist[simulation_rep]
  rep_id = pad_number(simulation_rep + rep_start - 1, n_reps)

  prefix = "eigenstrat_rep$rep_id"

  run(`Rscript ./scripts/findgraphs_1rep.R -i $findgraph_folder -o $findgraph_folder -p $prefix -k $num_admix --stop_gen $stop_gen --outgroup $outgroup --runs $runs -b $blgsize -r $simulation_rep -s $seed_array_findgraphs_path --output_graph_suffix $output_graph_suffix --output_f2_suffix $output_f2_suffix --output_summary_table_suffix $output_summary_table_suffix`) 
  
end 

#-----------------------------------------------#       
#         Summarize results  
#-----------------------------------------------# 
function summarize_findgraphs_results(simulation_rep::Int, findgraph_folderlist::Vector, rep_start::Int, n_reps::Int, output_summary_table_suffix::String)
  findgraph_folder = findgraph_folderlist[simulation_rep]
  rep_id = pad_number(simulation_rep + rep_start - 1, n_reps)

  tables = Dict(0 => "rep$(simulation_rep)_admix0$output_summary_table_suffix",
                  1 => "rep$(simulation_rep)_admix1$output_summary_table_suffix")
  
  result = Dict("rep_id" => rep_id, "best_K" => "None", "reject_H0" => false)

  for k in sort(collect(keys(tables)))
    table_path = joinpath(findgraph_folder, tables[k])
    df = CSV.read(table_path, DataFrame, delim='\t')
    valid_graphs = df[df.WR_smaller_than_3 .== true .&& df.is_outgroup_correct .== true, :]
    if nrow(valid_graphs) > 0
        result["best_K"] = k
        result["reject_H0"] = true
        break
    end
  end
  return result
end 

#-----------------------------------------------#       
#           Run and time the process
#-----------------------------------------------# 
function main() 

  @timeit to "Variant Calling + Convert VCF to eigenstrat" begin 
    pmap(simulation_rep -> callVariant_convertVCF_1rep(simulation_rep, findgraph_folderlist, folder_path_list, rep_start, n_reps), 1:index_length)
  end 

  @timeit to "Run findgraphs with num_admix=0 for $runs runs from rep$rep_start to rep$rep_end" begin 
    num_admix = 0
    pmap(simulation_rep -> run_findgraphs_1rep(simulation_rep, findgraph_folderlist, rep_start, n_reps, num_admix, stop_gen, outgroup, runs, blgsize, seed_array_findgraphs_path, output_graph_suffix, output_f2_suffix, output_summary_table_suffix), 1:index_length)
  end 

  @timeit to "Run findgraphs with num_admix=1 for $runs runs from rep$rep_start to rep$rep_end" begin 
    num_admix = 1
    pmap(simulation_rep -> run_findgraphs_1rep(simulation_rep, findgraph_folderlist, rep_start, n_reps, num_admix, stop_gen, outgroup, runs, blgsize, seed_array_findgraphs_path, output_graph_suffix, output_f2_suffix, output_summary_table_suffix), 1:index_length)
  end

  @timeit to "Organize all results" begin 
    for simulation_rep in 1:index_length 
      result = summarize_findgraphs_results(simulation_rep, findgraph_folderlist, rep_start, n_reps, output_summary_table_suffix)
      output_file = joinpath(outfolder, "findgraphs_summary_results.csv") 
      CSV.write(output_file, result)
    end 
  end 

end 

if abspath(PROGRAM_FILE) == @__FILE__
  main()

  n_processors = nprocs() 
  host_name = gethostname()  

  # save all parameters and write into the output dir
  arguments = """
  #-----------------AdmixtureTools2---------------------#
  Arguments used for simulating this output dataset: 
  master seed = $master_seed, unique master_seed for this parameter setting 
  Seed for findgraph= $seed_findgraphs, used to generate seed array (replicate x 2) used for findgraphs
  Find graph specific patameters: 
  Number of runs to run findgraphs = $runs
  Block size (blgsize) used in findgraphs = $blgsize 
  Number of generations to stop running findgraphs = $stop_gen
  Number of processors used = $n_processors;
  Server for running the script = $host_name.stat.wisc.edu

  Running time shown below: 
  """
  argument_file = joinpath(outfolder, "arguments-$paramname_root.log")
  write(argument_file, arguments)
  open(argument_file, "a") do io
    show(io, to)  
  end

end







