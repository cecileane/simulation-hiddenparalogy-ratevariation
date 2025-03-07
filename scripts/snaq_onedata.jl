# The script to run snaq 
using CSV
using DataFrames
using Distributed
using ArgParse
using PhyloNetworks
include("utilities.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        # outfolder, raxmlfolder, astralfolder, sanqfolder are simply passed down from simulation.jl to create meaingful filenames for snaq output 

        "--outfolder" 
        help = "root folder to be specified" 
        arg_type = String
        required = true

        "--iqtreefolder"
        help = "iqtree output folder"
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

        "--paramname_root" # create meaningful file names for snaq outputs
        help = "Parameter set used to simulate sim-phy"
        arg_type = String
        required = true

        "--seed_net0"
        help = "Seed for estimating tree Hmax = 0"
        arg_type = Int
        required = true

        "--seed_net1"
        help = "Seed for estimating network Hmax = 1"
        arg_type = Int
        required = true 

        "--seed_net0_bs"
        help = "bootnet seed Hmax = 0"
        arg_type = Int 
        required = true 

        "--seed_net1_bs"
        help = "bootnet seed Hmax = 1" 
        arg_type = Int
        required = true

        "--runs"
        help = "Number of runs to run SnaQ"
        arg_type = Int
        required = true

        "--n_snaqboot_rep"
        help = "Number of replicate for SnaQ boostrap for Hmax = 1" 
        arg_type = Int
        required = true

        "--processors"
        help = "Number of processors"
        arg_type = Int
        required = true

    end 
    return parse_args(s)
end

parsed_args = parse_commandline()

# Get arguments: 
outfolder = parsed_args["outfolder"]
iqtreefolder = parsed_args["iqtreefolder"]
astralfolder = parsed_args["astralfolder"]
snaqfolder = parsed_args["snaqfolder"]
paramname_root = parsed_args["paramname_root"]
processors = parsed_args["processors"]

# Get args used for SnaQ: 
runs = parsed_args["runs"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]
seed_net0 = parsed_args["seed_net0"]
seed_net0_bs = parsed_args["seed_net0_bs"]
seed_net1 = parsed_args["seed_net1"]
seed_net1_bs = parsed_args["seed_net1_bs"]

# Specify folder path: 
species_tree = joinpath(astralfolder,"astral.tre") # Astral tree 
gene_trees = joinpath(iqtreefolder,"besttrees.tre") # IQtree trees
outdir = snaqfolder # output directory for snaq. Each rep gets its own folder

# Make output directories for H=0 and H=1
H0folder = joinpath(outdir, "H0_output")
H1folder = joinpath(outdir, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

# Load starting topology, gene trees, and CFs
# species_tree = readMultiTopology(species_tree)[102] # starting topology -- the last (102th) tree
species_tree = readMultiTopology(species_tree)[1]
gene_trees = readMultiTopology(gene_trees)
q,t = countquartetsintrees(gene_trees)
df = writeTableCF(q,t)
CSV.write("$outdir/CF_results.csv", df) # save CFs into a file 
iqtreeCF = readTableCF(df)

@everywhere using PhyloNetworks 
addprocs(processors)

# hmax = 0 --> The starting tree for Hmax = 1 
# Below estimates network 
println("Parameter setting for Snaq Hmax = 0:\n runs = $runs, seed = $seed_net0. Running...")
net0out = joinpath(H0folder, "H0")
net0 = snaq!(species_tree, iqtreeCF, hmax=0, filename=net0out, seed=seed_net0, runs=runs)

println("Parameter setting for Snaq Hmax = 1:\n runs = $runs, seed = $seed_net1. Running...")
net1out = joinpath(H1folder, "H1")
net1 = snaq!(net0, iqtreeCF, hmax=1, filename=net1out, seed=seed_net1, runs=runs)

# Boostrapping
println("bootsnaq Hmax = 0:\n runs = $runs, seed = $seed_net0_bs. Running...")
bslist = joinpath(iqtreefolder, "bslist.txt") 
gene_tree_list = [ [tree] for tree in gene_trees ] # Convert to Vector{Vector{HybridNetwork}}

# Boostrapping H = 0:  
bootnet0 = bootsnaq(species_tree, gene_tree_list, hmax=0, nrep=n_snaqboot_rep, filename=net1out, seed=seed_net0_bs, runs=runs)
net0 = readTopology("$net0out.out")
BSe_tree0, tree0 = treeEdgesBootstrap(bootnet0,net0)
bootnet0_output = joinpath(H0folder,"bootnet0_onTree0.csv")
CSV.write(bootnet0_output, BSe_tree0)

# Boostrapping H = 1: 
println("bootsnaq Hmax = 0:\n runs = $runs, seed = $seed_net1_bs. Running...")
bootnet1 = bootsnaq(net0, gene_tree_list, hmax=1, nrep=n_snaqboot_rep, filename=net1out, seed=seed_net1_bs, runs=runs)
net1 = readTopology("$net1out.out")
BSe_tree1, tree1 = treeEdgesBootstrap(bootnet1,net1)
bootnet1_output = joinpath(H1folder,"bootnet1_onNet1.csv")
CSV.write(bootnet1_output, BSe_tree1)


