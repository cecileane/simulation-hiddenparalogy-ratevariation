# SnaQ submission files 
# This script call snaq_onedata.jl to run SnaQ. It allows to run specific replicate (red_start to rep_end). 

using ArgParse 

function parse_commandline()
    s = ArgParseSettings() 
    @add_arg_table s begin 
      
      # Specify which rep and which parameter sets to run SnaQ: 
      "--duploss"
        help = "Paramater setting (duplication rate) to run snaq"
        arg_type = Float64
        required = true
      "--ratevar"
        help ="Parameter setting (variation rate) to run snaq"
        arg_type = String
        required = true
      "--rep_start"
        help = "Parameter setting (start of the range of replicates) to run snaq"
        arg_type = Int
        required = true
      "--rep_end"
        help = "Parameter setting (end of the range of replicates) to run snaq"
        arg_type = Int
        required = true
      
      # Specify arguments for SnaQ: 
      "--runs"
        help = "Number of runs in SnaQ"
        arg_type = Int
        required = true
      "--threads"
        help = "Number of threads to run SnaQ. If unspecified, the defult is set to the maximum number of available threads"
        arg_type = Int
        default = Sys.CPU_THREADS # default sets to max number of CPUs
      "--seed_snaq"
        help = "Seed to run SnaQ"
        arg_type = Int
        required = true
      "--n_snaqboot_rep" 
        help = "Number of replicate for SnaQ boostrap for Hmax = 1" 
        arg_type = Int
        required = true
    end 
    return parse_args(s)
  end
  
parsed_args = parse_commandline()

# Parse arguments:  
duploss = parsed_args["duploss"]
ratevar = parsed_args["ratevar"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
runs = parsed_args["runs"]
threads = parsed_args["threads"]
seed_snaq = parsed_args["seed_snaq"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]

# set up folders: 
paramname_root = "DL$duploss-RV$ratevar" # Specify to find the folder
outfolder = "output/$paramname_root"
raxmlfolder  = joinpath(outfolder, "raxml-outfiles")
astralfolder = joinpath(outfolder, "astral-outfiles")
snaqfolder = joinpath(outfolder, "snaq-outfiles")
mkpath(snaqfolder)

#-----------------------------------------------#       
#       Estimate Networks using SnaQ
#-----------------------------------------------#
for simulation_rep in rep_start:rep_end
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(rep_end+1)), '0')
  # Some args used were passed down to snaq.jl script to create folder
  run(`julia ./scripts/snaq_onedata.jl --outfolder $outfolder --raxmlfolder $raxmlfolder --astralfolder $astralfolder --snaqfolder $snaqfolder --paramname_root $paramname_root --simulation_rep $simulation_rep --seed_snaq $seed_snaq --threads $threads --runs $runs --n_snaqboot_rep $n_snaqboot_rep`) 
end
