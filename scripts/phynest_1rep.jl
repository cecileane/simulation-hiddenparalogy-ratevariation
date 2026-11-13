#= 
PhyNEST for one replicate, called by phynest.jl (not run directly)
Runs phyne! with hmax=0 (tree) and then hmax=1 (one reticulation) 
on the concatenated alignment,
then saves the negative log composite likelihood of both, 
of the true species tree,
and T = score(H0) - score(H1) used for model selection
Input: seqgenfolder/concate_alignment_rep<id>.fasta
     astralfolder/astral.tre, the starting topology (tips A..H)
Output, all in phynestfolder: rep<id>.phy (converted alignment),
    H0_output/H0_hc.{log,out}, H1_output/H1_hc.{log,out}, 
    phynest_results.csv (one row)

Environment setting see phynest.jl, and it cannot be the same env as utilities, so that 
the few functions needed here are copied below.
PhyNEST scores are negative log composite likelihoods (smaller = better fit),
so T = score_H0 - score_H1 >= 0 and a larger T means the reticulation helps more
=# 

# activate envs/phynest so --project is not needed on the command line 

using Pkg
phynest_env = normpath(joinpath(@__DIR__, "..", "envs", "phynest"))
Pkg.activate(phynest_env; io=devnull)
Pkg.instantiate(; io=devnull)

using ArgParse
using CSV
using DataFrames
using PhyNEST
using Random

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        # folders, passed from phynest.jl

        "--seqgenfolder"
        help = "seq-gen output folder with the concatenated fasta"
        arg_type = String
        required = true

        "--astralfolder"
        help = "astral output folder (starting topology)"
        arg_type = String
        required = true

        "--phynestfolder"
        help = "phynest output folder"
        arg_type = String
        required = true

        "--rep_id"
        help = "Zero-padded replicate ID string (e.g. 01)"
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

        "--runs"
        help = "Number of independent runs per PhyNEST search"
        arg_type = Int
        required = true

        "--n_inds"
        help = "Number of individuals per species"
        arg_type = Int
        required = true

        "--max_steps"
        help = "Max steps per search (PhyNEST default = 250000)"
        arg_type = Int
        default = 250000

        "--outgroup"
        help = "Outgroup taxon used to root the network (default: A)"
        arg_type = String
        default = "A"

        "--true_tree_newick"
        help = "True species tree in Newick format"
        arg_type = String
        default = "(A,((((B,C),(D,E)),F),(G,H)));"
    end
    return parse_args(s)
end

#-------------------------------------#
# Local helper functions, see phynest.jl for understanding why helper functions need to
# be copied here.
 
#-------------------------------------#
# same as clean_newick_nan in utilities.jl, which we can't include here (see above)
"""
Input: newick_str (String): newick from ASTRAL, can have :nan branch lengths
Output: same string, with :nan replaced by a small branch length
"""
function clean_newick_nan(newick_str::String)
    return replace(newick_str, ":nan" => ":0.000001")
end

"""
Input: newick_str (String)
Output: newick without internal node labels, 
    (ASTRAL writes support values right after a closing parenthesis)
Notes: PhyNEST's runs returns Inf on a topology with such labels,
and the labels of the starting tree get into every network of the search,
 so remove them first
"""
function strip_support_labels(newick_str::String)
    return replace(newick_str, r"\)[0-9.eE+-]+" => ")")
end

"""
Input:
-net: topology (HybridNetwork of the PhyloNetworks version inside PhyNEST)
-outgroup (String): taxon to root at
Output: same topology rooted at the outgroup, branch lengths dropped
Notes: same unroot / re-root steps as 
    PhyNEST does on its own proposals (preNNI in Searches.jl)
ASTRAL trees have an arbitrary root, 
    and PhyNEST's clock model depends on the root of the starting tree
"""
function root_at_outgroup(net, outgroup::String)
    PN = PhyNEST.PhyloNetworks
    unrooted = PN.readTopologyUpdate(PN.writeTopologyLevel1(net))
    for e in unrooted.edge
        e.length = -1.0
    end
    return readTopology(PN.writeTopologyLevel1(unrooted, outgroup))
end

"""
Input: logpath (String): .log file written by phyne!
Output: smallest negative log composite likelihood over all runs,
that's the score of the network that phyne! returns. 
"""
function best_score_from_log(logpath::String)
    scores = Float64[]
    pattern = r"in this run: ([-+0-9.eE]+)"
    for line in eachline(logpath)
        m = match(pattern, line)
        m === nothing && continue
        s = tryparse(Float64, m.captures[1])
        s === nothing || isnan(s) || push!(scores, s)
    end
    return isempty(scores) ? NaN : minimum(scores)
end

"""
Input: net (topology), p (Phylip data read by readPhylip)
Output: negative log composite likelihood of that topology,
after optimizing branch lengths, gammas and theta with do_optimization
"""
function score_network(net, p)
    res, _ = do_optimization(net, p)
    return res.minimum
end

#-------------------------------------#
# Parse commandline arguments
#-------------------------------------#
parsed_args = parse_commandline()

seqgenfolder = parsed_args["seqgenfolder"]
astralfolder = parsed_args["astralfolder"]
phynestfolder = parsed_args["phynestfolder"]
rep_id = parsed_args["rep_id"]
seed_net0 = parsed_args["seed_net0"]
seed_net1 = parsed_args["seed_net1"]
runs = parsed_args["runs"]
n_inds = parsed_args["n_inds"]
max_steps = parsed_args["max_steps"]
outgroup = parsed_args["outgroup"]
true_tree_newick = parsed_args["true_tree_newick"]

# output folders for H=0 and H=1, same layout as snaqfolder
H0folder = joinpath(phynestfolder, "H0_output")
H1folder = joinpath(phynestfolder, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

#----------------------------------#
# Convert alignment and parse inputs
#-----------------------------------#
fasta_file = joinpath(seqgenfolder, "concate_alignment_rep$(rep_id).fasta")
phylip_file = joinpath(phynestfolder, "rep$(rep_id).phy")

# fasta to sequential phylip with Biopython (scripts/fasta2phylip.py), 
# same library as concatenate_seq.py
run(`python3 scripts/fasta2phylip.py \
    -i $fasta_file \
    -o $phylip_file \
    --n_inds $n_inds`)

p = readPhylip(phylip_file, showProgress=false)
# readPhylip returns nothing on a parse error instead of failing, so check here
p isa PhyNEST.Phylip ||
    error("readPhylip failed on $phylip_file for rep$rep_id")

# starting tree = ASTRAL tree (tips already A..H), support labels removed, rooted at the outgroup
species_tree_content = read(joinpath(astralfolder, "astral.tre"), String)
species_tree_content = clean_newick_nan(species_tree_content)
species_tree_content = strip_support_labels(species_tree_content)
start_tree = root_at_outgroup(readTopology(species_tree_content), outgroup)

#-------------------------------------#
# Run PhyNEST
#-------------------------------------#
# phyne! has no seed argument, so seed the global RNG before each search
# each replicate is its own julia process, so this is stable
# can't use StableRNG here: phyne! calls rand() on the default RNG
# results are reproducible for a fixed julia version 
# (same as snaq! with its seed argument)
# only the seed generation in utilities.jl needs StableRNG
# same number of runs for H=0 and H=1 so that the T comparison is fair

# hmax = 0 --> tree search; also the starting point for Hmax = 1
println("PhyNEST Hmax=0: runs=$runs, seed=$seed_net0. Running...") 

Random.seed!(seed_net0) # see above comments for the reason here 

net0 = phyne!(start_tree, p, outgroup; hmax=0,
    number_of_runs=runs, maximum_number_of_steps=max_steps,
    filename=joinpath(H0folder, "H0"))
net0 isa PhyNEST.HybridNetwork ||
    error("PhyNEST Hmax=0 failed for rep$rep_id: $net0")

# hmax = 1 --> the network with one hybrid edge
println("PhyNEST Hmax=1: runs=$runs, seed=$seed_net1. Running...")

Random.seed!(seed_net1) # see above comments for the reason here   

net1 = phyne!(net0, p, outgroup; hmax=1,
    number_of_runs=runs, maximum_number_of_steps=max_steps,
    filename=joinpath(H1folder, "H1"))
net1 isa PhyNEST.HybridNetwork ||
    error("PhyNEST Hmax=1 failed for rep$rep_id: $net1")

#-------------------------------------#
# Post-PhyNEST analyses
# Extract scores for H = 0 and H = 1
# T = score_H0 - score_H1 (both are negative log composite
#     likelihoods, so T >= 0 up to search noise; large T favors H = 1)
# Score the true species tree (fixed topology) for reference
# Extract gamma values from the H = 1 network

#= scores from the per-run lines of the phyne! .log 
(what the search actually reached)
re-optimizing the returned network is only the backup: 
    do_optimization can give Inf on some networks
=# 
score_net0 = best_score_from_log(joinpath(H0folder, "H0_hc.log"))
isnan(score_net0) && (score_net0 = score_network(net0, p))

score_net1 = best_score_from_log(joinpath(H1folder, "H1_hc.log"))
isnan(score_net1) && (score_net1 = score_network(net1, p))

T_stat = score_net0 - score_net1 # get the test statistic 

true_tree = readTopology(true_tree_newick)
score_true = score_network(true_tree, p)

println("score_H0 = $score_net0")
println("score_H1 = $score_net1")
println("T = $T_stat")
println("score_truetree = $score_true")

gammas = [e.gamma for e in net1.edge if e.hybrid]
println("Number of hybrid edges found: ", length(gammas))
println("Gamma values: ", gammas)

if length(gammas) >= 2
    gamma_1, gamma_2 = sort(gammas, rev=true)[1:2]
else
    gamma_1, gamma_2 = NaN, NaN
end

#-----------------------------------#
# Save one-row results table
#------------------------------------#
results = DataFrame(
        rep_id = [rep_id],
        score_H0 = [score_net0],
        score_H1 = [score_net1],
        T = [T_stat],
        score_truetree = [score_true],
        num_hybrid_H1 = [net1.numHybrids],
        gamma_1 = [gamma_1],
        gamma_2 = [gamma_2],
        net_H0 = [PhyNEST.writeTopologyLevel1(net0)],
        net_H1 = [PhyNEST.writeTopologyLevel1(net1)],
    )

results_path = joinpath(phynestfolder, "phynest_results.csv")
CSV.write(results_path, results)
println("Rep$rep_id: PhyNEST results saved to $results_path")
