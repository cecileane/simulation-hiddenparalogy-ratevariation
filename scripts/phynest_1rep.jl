# ============================================================================
# scripts/phynest_1rep.jl
#
# Purpose : Per-replicate PhyNEST worker. Runs PhyNEST at hmax=0 (tree) and
#           hmax=1 (one reticulation) on the concatenated alignment of a
#           single replicate, then records the negative log composite
#           likelihood of both models, of the true species tree, and the
#           statistic T = score(H0) - score(H1) used for model selection
#           (analogous to SNaQ GoF p-values and find_graphs WR thresholds).
# Inputs  : <seqgenfolder>/concate_alignment_rep<id>.fasta
#           <astralfolder>/astral.tre  (starting topology, tips A..H)
# Outputs : <phynestfolder>/rep<id>.phy            converted alignment
#           <phynestfolder>/H0_output/H0_hc.{log,out}
#           <phynestfolder>/H1_output/H1_hc.{log,out}
#           <phynestfolder>/phynest_results.csv    one-row score summary
# Usage   : Not run directly; called by phynest.jl across replicates as
#               julia scripts/phynest_1rep.jl ...
#           (the script activates envs/phynest itself, so the
#           --project flag is optional; passing it is harmless)
# Note    : PhyNEST v0.1.x pins PhyloNetworks < 1.0, which conflicts with
#           the SNaQ/PhyloSummaries stack, so this worker runs in its own
#           Julia environment (envs/phynest) and does not include
#           utilities.jl (written against PhyloNetworks >= 1.0). The few
#           helpers needed here are defined locally instead.
#           PhyNEST scores are negative log composite likelihoods
#           (smaller = better fit), so T = score_H0 - score_H1 >= 0 and
#           larger T means the reticulation improves the fit more.
# ============================================================================

# Activate the dedicated PhyNEST environment (envs/phynest) so this
# script also works without `--project=envs/phynest` on the command
# line. The path is resolved relative to this file, so it works from
# any working directory. Pkg.instantiate() installs the exact
# versions pinned in envs/phynest/Manifest.toml on first use and is
# a fast no-op afterwards.
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
        # Folders passed down from phynest.jl for meaningful filenames

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
# Local helper functions
#-------------------------------------#

# Same fix as clean_newick_nan in utilities.jl (not includable here,
# see header): ASTRAL may write ":nan" branch lengths.
"""
    clean_newick_nan(newick_str)

Replace `:nan` in ASTRAL newick output with a small branch length.
"""
function clean_newick_nan(newick_str::String)
    return replace(newick_str, ":nan" => ":0.000001")
end

"""
    strip_support_labels(newick_str)

Remove internal node labels (e.g. ASTRAL support values written
directly after a closing parenthesis). PhyNEST's optimization returns
Inf on topologies carrying such labels, and labels propagate from the
starting tree into every searched network, so they must be removed.
"""
function strip_support_labels(newick_str::String)
    return replace(newick_str, r"\)[0-9.eE+-]+" => ")")
end

"""
    root_at_outgroup(net, outgroup)

Root a topology at `outgroup` and drop branch lengths, using the same
unroot/re-root path PhyNEST applies to its own search proposals
(see preNNI in PhyNEST's Searches.jl). ASTRAL trees come arbitrarily
rooted, and PhyNEST's clock-based model is sensitive to the rooting
of the starting topology.
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
    best_score_from_log(logpath)

Smallest per-run negative log composite likelihood reported in a
phyne! .log file, i.e. the score of the network phyne! returns.
Returns NaN if no per-run score line is found.
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
    score_network(net, p)

Negative log composite likelihood of a fixed topology, with branch
lengths, gammas and theta re-optimized by `do_optimization`.
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

# Make output directories for H=0 and H=1 (mirrors snaqfolder layout)
H0folder = joinpath(phynestfolder, "H0_output")
H1folder = joinpath(phynestfolder, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

#-------------------------------------#
# Convert alignment and parse inputs
#-------------------------------------#
fasta_file = joinpath(seqgenfolder, "concate_alignment_rep$(rep_id).fasta")
phylip_file = joinpath(phynestfolder, "rep$(rep_id).phy")

# FASTA -> sequential PHYLIP via Biopython (scripts/fasta2phylip.py),
# the same library concatenate_seq.py uses; see that script's docstring
run(`python3 scripts/fasta2phylip.py \
    -i $fasta_file \
    -o $phylip_file \
    --n_inds $n_inds`)

p = readPhylip(phylip_file, showProgress=false)
# readPhylip swallows parse errors and returns nothing; fail loudly here
p isa PhyNEST.Phylip ||
    error("readPhylip failed on $phylip_file for rep$rep_id")

# Starting topology: ASTRAL species tree (tips are already A..H),
# support labels removed and rooted at the outgroup (see helpers above)
species_tree_content = read(joinpath(astralfolder, "astral.tre"), String)
species_tree_content = clean_newick_nan(species_tree_content)
species_tree_content = strip_support_labels(species_tree_content)
start_tree = root_at_outgroup(readTopology(species_tree_content), outgroup)

#-------------------------------------#
# Run PhyNEST
#-------------------------------------#
# phyne! has no seed argument, so seed the global RNG before each
# search; each replicate runs in its own process, so this is stable.
# StableRNG cannot be used here: phyne!'s internal rand() calls draw
# from Julia's task-default RNG, which a StableRNG object cannot
# replace. Results are exactly reproducible for a fixed Julia
# version (same property as snaq!'s seed argument); only the seed
# *generation* in utilities.jl needs StableRNG, and it has it.
# The same number of runs is used for H=0 and H=1 to keep the
# composite likelihood comparison (T statistic) fair.

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
# --> Extract scores for H = 0 and H = 1
# --> T = score_H0 - score_H1 (both are negative log composite
#     likelihoods, so T >= 0 up to search noise; large T favors H = 1)
# --> Score the true species tree (fixed topology) for reference
# --> Extract gamma values from the H = 1 network
#-------------------------------------#
# Scores come from the per-run lines of the phyne! .log (the exact
# values the search achieved); re-optimizing the returned network is
# only a fallback because do_optimization can fail (Inf) on some
# network shapes.
score_net0 = best_score_from_log(joinpath(H0folder, "H0_hc.log"))
isnan(score_net0) && (score_net0 = score_network(net0, p))
score_net1 = best_score_from_log(joinpath(H1folder, "H1_hc.log"))
isnan(score_net1) && (score_net1 = score_network(net1, p))
T_stat = score_net0 - score_net1

true_tree = readTopology(true_tree_newick)
score_true = score_network(true_tree, p)

println("score_H0 = $score_net0")
println("score_H1 = $score_net1")
println("T = $T_stat")
println("score_truetree = $score_true")

# With hmax=1, we expect 1 hybrid node, i.e. 2 hybrid edges with
# gamma and 1-gamma. The search may still return fewer (e.g. a tree);
# record NaN gammas in that case instead of failing the replicate.
gammas = [e.gamma for e in net1.edge if e.hybrid]
println("Number of hybrid edges found: ", length(gammas))
println("Gamma values: ", gammas)

if length(gammas) >= 2
    gamma_1, gamma_2 = sort(gammas, rev=true)[1:2]
else
    gamma_1, gamma_2 = NaN, NaN
end

#-------------------------------------#
# Save one-row results table
#-------------------------------------#
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
