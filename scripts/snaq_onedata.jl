# The script to run snaq 

using PhyloNetworks 
using Distributed
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        # outfolder, raxmlfolder, astralfolder, sanqfolder are simply passed down from simulation.jl to create meaingful filenames for snaq output 

        "--outfolder" 
        help = "root folder to be specified" 
        arg_type = String
        required = true

        "--raxmlfolder"
        help = "raxml output folder"
        arg_type = String
        required = true

        "--astralfolder"
        help = "astral output folder"
        arg_type = String
        required = true

        "--snaqfolder"
        help = "snaq output folder"
        arg_type = String
        required = true

        "--paramname_root" # This is important to create meaningful file names for snaq outputs
        help = "Parameter set used to simulate sim-phy"
        arg_type = String
        required = true

        "--simulation_rep" # This is important to create meaingful file names for snaq outputs
        help = "Rep ID for simulation"
        arg_type = Int
        required = true

        "--runs"
        help = "Number of runs to run SnaQ"
        arg_type = Int
        required = true
        
        "--threads" 
        help = "Number of threads to be used"
        arg_type = Int
        required = true

        "--seed_snaq"
        help = "Seed for SnaQ"
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

# Get arguments: 
outfolder = parsed_args["outfolder"]
raxmlfolder = parsed_args["raxmlfolder"]
astralfolder = parsed_args["astralfolder"]
snaqfolder = parsed_args["snaqfolder"]
paramname_root = parsed_args["paramname_root"]
simulation_rep = parsed_args["simulation_rep"]

# Get args used for SnaQ: 
runs = parsed_args["runs"]
threads = parsed_args["threads"]
seed_snaq = parsed_args["seed_snaq"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]

# Specify folder path: 
species_tree = joinpath(astralfolder,"astral-outfiles-$paramname_root-$simulation_rep", "astral.tre") # Astral tree 
gene_trees = joinpath(raxmlfolder,"raxml-outfiles-$paramname_root-$simulation_rep", "besttrees.tre") # RAxML trees
outdir = joinpath(snaqfolder, "rep$simulation_rep") # output directory for snaq. Each rep gets its own folder
bslist = joinpath(astralfolder, "astral-outfiles-$paramname_root-$simulation_rep", "BSlistfiles") # BS files to run boostrap 

# Make output directories for H=0 and H=1
H0folder = joinpath(outdir, "H0_output")
H1folder = joinpath(outdir, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

# Multi-thread the pipeline
addprocs(threads)
@everywhere using PhyloNetworks

# Load starting topology, gene trees, and CFs
species_tree = readMultiTopology(species_tree)[102] # starting topology -- the last (102th) tree
gene_trees = readMultiTopology(gene_trees)
q,t = countquartetsintrees(gene_trees)
df = writeTableCF(q,t)
raxmlCF = readTableCF(df)

# hmax = 0 or 1
# Question: do I need two different seeds for H=1 and H=0 and bootnet? 
# Below estimates network 
println("Parameter setting for Snaq Hmax = 0:\n runs = $runs, seed = $seed_snaq, threads = $threads")
net0 = snaq!(species_tree, raxmlCF, hmax=0, filename=joinpath(H0folder, "H0"), seed=seed_snaq, runs=runs)
net1 = snaq!(net0, raxmlCF, hmax=1, filename=joinpath(H1folder, "H1"), seed=seed_snaq, runs=runs)

# Incorporate bootrapping for hmax = 1
# In simulation.jl the output of raxml and astral was moved into the folder so the filenames listed in BSListfiles starts from $root instead of the correct file 
println("running snaq boostrap for $n_snaqboot_rep reps and Hmax = 1")
# Below setting bslist = path doesn't work. 
write("bslist_revised.txt", join((joinpath(outfolder, "raxml-outfiles", line) for line in readlines(bslist)), "\n")) # revise the bslist 
bootTrees = readBootstrapTrees("bslist_revised.txt") # The bslist_revised store the path from the root to raxml boostrap file
bootnet = bootsnaq(net0, bootTrees, hmax=1, nrep= n_snaqboot_rep, runs= runs, filename=joinpath(H1folder, "bootnet_H1"), seed=seed_snaq)

rm("bslist_revised.txt") # removed intermediate file

