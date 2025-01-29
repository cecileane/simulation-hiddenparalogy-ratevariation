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

        "--runs"
        help = "Number of runs to run SnaQ"
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
iqtreefolder = parsed_args["iqtreefolder"]
astralfolder = parsed_args["astralfolder"]
snaqfolder = parsed_args["snaqfolder"]
paramname_root = parsed_args["paramname_root"]
# simulation_rep = parsed_args["simulation_rep"]

# Get args used for SnaQ: 
runs = parsed_args["runs"]
seed_snaq = parsed_args["seed_snaq"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]

# Specify folder path: 
species_tree = joinpath(astralfolder,"astral.tre") # Astral tree 
gene_trees = joinpath(iqtreefolder,"besttrees.tre") # IQtree trees
outdir = snaqfolder # output directory for snaq. Each rep gets its own folder

# Make output directories for H=0 and H=1
H0folder = joinpath(outdir, "H0_output")
H1folder = joinpath(outdir, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

# Multi-thread the pipeline
addprocs(Sys.CPU_THREADS)
@everywhere using PhyloNetworks

# Load starting topology, gene trees, and CFs
# species_tree = readMultiTopology(species_tree)[102] # starting topology -- the last (102th) tree
species_tree = readMultiTopology(species_tree)[1]
gene_trees = readMultiTopology(gene_trees)
q,t = countquartetsintrees(gene_trees)
df = writeTableCF(q,t)
iqtreeCF = readTableCF(df)

# hmax = 0 or 1
# Question: do I need two different seeds for H=1 and H=0 and bootnet? 
# Below estimates network 
println("Parameter setting for Snaq Hmax = 0:\n runs = $runs, seed = $seed_snaq. Running...")
net0 = snaq!(species_tree, iqtreeCF, hmax=0, filename=joinpath(H0folder, "H0"), seed=seed_snaq, runs=runs)
net1out = joinpath(H1folder, "H1")
net1 = snaq!(net0, iqtreeCF, hmax=1, filename=net1out, seed=seed_snaq, runs=runs)

#=  Boostrapping: 
bootnetout_dir = joinpath(H1folder, "bootnet_H1")
bootnet = bootsnaq(net0, iqtreeCF, hmax=1, nrep=n_snaqboot_rep, filename=bootnetout_dir, seed=seed_snaq, runs=runs)
bootnet_networks = readMultiTopology(bootnetout);
net1 = readTopology("$net1out.out")
BSe_tree, tree1 = treeEdgesBootstrap(bootnet,net1)

bslist = joinpath(astralfolder, "astral-outfiles-$paramname_root-$simulation_rep", "BSlistfiles") # BS files to run boostrap 
# Incorporate bootrapping for hmax = 1
# In simulation.jl the output of raxml and astral was moved into the folder so the filenames listed in BSListfiles starts from $root instead of the correct file 
println("running snaq boostrap for $n_snaqboot_rep reps and Hmax = 1")
# Below setting bslist = path doesn't work. 
write("bslist_revised.txt", join((joinpath(outfolder, "raxml-outfiles", line) for line in readlines(bslist)), "\n")) # revise the bslist 
bootTrees = readBootstrapTrees("bslist_revised.txt") # The bslist_revised store the path from the root to raxml boostrap file
bootnet = bootsnaq(net0, bootTrees, hmax=1, nrep= n_snaqboot_rep, runs= runs, filename=joinpath(H1folder, "bootnet_H1"), seed=seed_snaq)

rm("bslist_revised.txt") # removed intermediate file
=# 
