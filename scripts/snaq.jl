# The script to run snaq 

using PhyloNetworks 
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--species_tree"
        help = "Path to species tree in Newich file"
        arg_type = String
        required = true

        "--gene_trees"
        help = "path to the file with all gene trees"
        arg_type = String
        required = true 

        "--hmax"
        help = "Max number of reticulation for SnaQ: 0 or 1 only"
        arg_type = Int # 0 or 1 
        required = true

        "--outdir"
        help = "Path to the output directory"
        arg_type = String
        required = true
    end 
    return parse_args(s)
end

parsed_ages = parde_commandline()

# Load species tree
species_tree = parsed_args["species_tree"]
gene_trees = parsed_agrs["--gene_trees"]
hmax = parsed_args["hmax"] 
outdir = parsed_args["outdir"]

# Load tree topology 
species_tree = readTopology(secies_tree)
gene_trees = readMultiTopology(gene_trees)

# If hmax = 1 then start to infer from hmax = 0 


