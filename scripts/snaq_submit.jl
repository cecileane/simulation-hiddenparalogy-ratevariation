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
for ind in 1:index_length 
  iqtreefolder  = setup_rep_output_folders(folder_path_list, ind, "iqtreefolder")
  astralfolder = setup_rep_output_folders(folder_path_list, ind, "astralfolder")
  snaqfolder = snaqfolder_list[ind] 
  # run snaq_onedata.jl: 
  run(`julia ./scripts/snaq_onedata.jl --outfolder $outfolder --iqtreefolder $iqtreefolder --astralfolder $astralfolder --snaqfolder $snaqfolder --paramname_root $paramname_root --seed_snaq $seed_snaq --runs $runs --n_snaqboot_rep $n_snaqboot_rep`) 
end
