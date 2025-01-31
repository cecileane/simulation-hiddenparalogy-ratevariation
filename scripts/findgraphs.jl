# This script documents using find_graphs to estimate graph 

# SnaQ submission files 
# This script call snaq_onedata.jl to run SnaQ. It allows to run specific replicate (red_start to rep_end). 

using ArgParse 
include("utilities.jl")

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

findgraph_folder_list = []
index_length = rep_end - rep_start + 1 # see utilities.jl to match how we set up rep_start:rep_end, 
for ind in 1:index_length 
  findgraph_folder = setup_rep_output_folders(folder_path_list, ind, "findgraphfolder")
  push!(findgraph_folder_list, findgraph_folder)
end 
check_existing_dir(findgraph_folder_list) # see utilies.jl. 

#-----------------------------------------------#       
#         Variant Calling using snp-sites 
#       snp-sites calls SNP from fasta files 
#-----------------------------------------------# 
for ind in 1:index_length 
    # concatenated fasta files are stored in Rep$id/seqgenfolder
    seqgenfolder = setup_rep_output_folders(folder_path_list, ind, "seqgenfolder")
    findgraphfolder = findgraph_folder_list[ind]
    mkpath(findgraphfolder)
    # Find the file name in the concatenated fasta file: 
    match_result = match(r"rep(\d+)/", seqgenfolder) 
    rep_id = match_result.captures[1] 
    fasta_file_name = "concate_alignment_rep$(rep_id).fasta" 
    fasta_file = joinpath(seqgenfolder, fasta_file_name) 
    println("for $seqgenfolder, here is the file name: $fasta_file_name")
    run(`./executables/snp-sites $fasta_file -v -o $findgraphfolder/snps.vcf`) # run snp-sites 
end 


