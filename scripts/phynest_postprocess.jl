#!/usr/bin/env julia
# ============================================================================
# scripts/phynest_postprocess.jl
#
# Purpose : Aggregate per-replicate PhyNEST results for one parameter set
#           into a single per-setting table, and add what the inference step
#           does not record: RF distances of the H=0 tree and of the two
#           trees displayed in the H=1 network to the true species tree
#           (full, and with F, G or H pruned) and to the three alternative
#           placements of F; the hybrid clade and its donors in the H=1
#           network; and, from the per-run search logs, where the true tree
#           sits in each search: how many runs reached the best score, how
#           many runs found (H=0) or displayed (H=1) the true tree, and the
#           score gap and rank of the best such run; and the number of
#           sites M of the concatenated alignment, read from the PhyNEST
#           log, and the per-site score difference T_per_site = T / M: the
#           raw difference T grows linearly with M, and M is smaller in
#           replicates where fewer loci were retained.
# Inputs  : <output_dir>/<paramname>/rep*/phynestfolder/phynest_results.csv
#           <output_dir>/<paramname>/rep*/phynestfolder/H?_output/H?_hc.log
# Outputs : <output_dir>/<paramname>/PhyNEST-<paramname>-summary.csv
#           (copied into phynest_summary/ by run_postprocessing.jl)
# Usage   : julia scripts/phynest_postprocess.jl --dup_rate 0.0 \
#               --loss_rate 0.0 --ratevar N --n_inds 1 --SF 0.5 \
#               --n_reps 100 [--output_dir output]
#           or, for every parameter set at once:
#           julia scripts/run_postprocessing.jl --mode phynest \
#               [--output_dir DIR]
# Note    : No threshold on T is applied here. Model selection (naive T > 0
#           and the calibrated T*) is done later in summary_phynest.jl, so the
#           threshold can change without re-running this script.
#           RF distances use mudistance_semidirected, as for SNaQ/find_graphs.
# ============================================================================

# Self-activate the repo root environment so this script never depends on
# the caller's default Julia environment. Pkg.instantiate() installs the
# exact versions pinned in Manifest.toml on first use and is a fast no-op
# afterwards.
using Pkg
repo_env = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(repo_env; io=devnull)
Pkg.instantiate(; io=devnull)

using ArgParse
using CSV
using DataFrames
using PhyloNetworks
include("utilities.jl") # pad_number, set_up_paramname_root, get_hybrid_info

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        # Parameter settings to identify the correct output folder
        "--dup_rate"
            help = "Parameter setting (duplication rate) used in phynest"
            arg_type = Float64
            required = true
        "--loss_rate"
            help = "Parameter setting (gene loss rate) used in phynest"
            arg_type = Float64
            required = true
        "--ratevar"
            help = "Parameter setting (variation rate) used in phynest"
            arg_type = String
            required = true
        "--n_reps"
            help = "Total number of reps for this parameter set"
            arg_type = Int
            default = 100
        "--rep_start"
            help = "Start of the range of replicates to process"
            arg_type = Int
            default = 1
        "--rep_end"
            help = "End of the range of replicates to process"
            arg_type = Int
            default = -1
        "--n_inds"
            help = "Number of individuals per species"
            arg_type = Int
            default = 1
        "--SF"
            help = "Scaling factor to scale effective population Ne
                (Default = 1.0, no scaling)"
            arg_type = Float64
            default = 1.0
        "--gene_len"
            help = "Gene length used in simulation (default = 1000)"
            arg_type = Int
            default = 1000
        "--output_dir"
            help = "Directory holding the parameter folders (default = output)"
            arg_type = String
            default = "output"
    end

    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1
        parsed_args["rep_end"] = parsed_args["n_reps"]
    end

    return parsed_args
end

parsed_args = parse_commandline()

# Parse arguments
dup_rate = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"]
n_reps = parsed_args["n_reps"]
n_inds = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
SF = parsed_args["SF"]
gene_len = parsed_args["gene_len"]
output_dir = parsed_args["output_dir"]

# Set up folders
paramname_root = set_up_paramname_root(dup_rate, loss_rate,
                                    ratevar, n_inds, SF,
                                    gene_len)
outfolder = joinpath(output_dir, paramname_root)

if !isdir(outfolder)
    error("Output folder does not exist: $outfolder")
end

println("PhyNEST Postprocessing Script")
println("Parameter set: $paramname_root")
println("Processing replicates $rep_start to $rep_end")
println("Output folder: $outfolder")

#-----------------------------------------------#
#   Reference topologies
#-----------------------------------------------#

true_species_tree = readnewick("(A,((((B,C),(D,E)),F),(G,H)));")

#= Same alternatives as in snaq_postprocess.jl: when an estimated tree is
   not the species tree, check whether F was simply moved.
   alter1: F with (B,C); alter2: F with (D,E); alter3: F outside all =#
alter1 = readnewick("(A,((G,H),(((B,C),F),(D,E))));")
alter2 = readnewick("(A,((G,H),((B,C),(F,(D,E)))));")
alter3 = readnewick("(A,(((G,H),((D,E),(B,C))),F));")

# True species tree with one taxon pruned
species_tree_pruned = Dict(
    "F" => readnewick("(A,((G,H),((B,C),(D,E))));"),
    "G" => readnewick("(A,((((B,C),(D,E)),F),H));"),
    "H" => readnewick("(A,((((B,C),(D,E)),F),G));"))

# Column suffixes, in the order returned by compare_to_truth
rf_names = ["true", "true_noF", "true_noG", "true_noH",
            "alter1", "alter2", "alter3"]
rf_columns(prefix) = Tuple(Symbol.(prefix * "_" .* rf_names))

rf(tree, ref) = mudistance_semidirected(tree, ref, preorder=true)

"""
RF distance to the true tree after pruning `taxon` from a copy of `tree`.
Returns NaN if the pruning fails.
"""
function rf_pruned(tree::HybridNetwork, taxon::String)
    pruned = deepcopy(tree)
    try
        deleteleaf!(pruned, taxon)
        return rf(pruned, species_tree_pruned[taxon])
    catch e
        println("Warning: could not prune $taxon: $e")
        return NaN
    end
end

"""
All topology comparisons for one tree: RF distance to the true tree,
to the true tree without F/G/H, and to the 3 alternative placements of F.
"""
function compare_to_truth(tree::HybridNetwork)
    return (rf(tree, true_species_tree),
            rf_pruned(tree, "F"), rf_pruned(tree, "G"), rf_pruned(tree, "H"),
            rf(tree, alter1), rf(tree, alter2), rf(tree, alter3))
end

no_comparison = Tuple(fill(NaN, length(rf_names)))

"""
Read a PhyNEST hill-climbing log and return the estimated topology (Newick)
and composite likelihood of every run, in run order. A run whose
optimization failed prints NaN and is dropped.
"""
function read_run_results(logfile::String)
    topologies = Dict{Int,String}()
    scores = Dict{Int,Float64}()
    isfile(logfile) || return (String[], Float64[])
    for line in eachline(logfile)
        m = match(r"^\((\d+)/\d+\) Estimated topology in this run: (.*)$", line)
        m === nothing || (topologies[parse(Int, m[1])] = String(m[2]))
        m = match(r"^\((\d+)/\d+\) Composite likelihood of the estimated " *
                  r"topology in this run: (\S+)$", line)
        m === nothing || (scores[parse(Int, m[1])] =
            something(tryparse(Float64, m[2]), NaN))
    end
    runs = sort([k for k in keys(scores)
                 if haskey(topologies, k) && !isnan(scores[k])])
    return ([topologies[k] for k in runs], [scores[k] for k in runs])
end

"""
Number of sites in the alignment analyzed by PhyNEST, from the header of
its log ("Sequence length: 951000"). Returns 0 if the line is missing.
"""
function read_sequence_length(logfile::String)
    isfile(logfile) || return 0
    for line in eachline(logfile)
        m = match(r"^Sequence length: (\d+)", line)
        m === nothing || return parse(Int, m[1])
    end
    return 0
end

# Scores within 1 unit of each other are treated as ties: differences below
# that come from the optimizer, not from the topology.
score_tol = 1.0

"""
Summarize one search (H=0 or H=1) from its per-run scores, the flag saying
which runs found/displayed the true tree, and the score of the true tree:
number of runs, runs at the best score, runs with the true tree, score gap
between the best such run and the best run (0 when the best run has it),
rank of the best such run (1 = best), and the percentage of runs that
scored better than the true tree itself.
"""
function summarize_runs(scores::Vector{Float64}, has_true::Vector{Bool},
                        score_truetree::Float64)
    n = length(scores)
    n == 0 && return (0, 0, 0, NaN, NaN, NaN)
    best = minimum(scores)
    if any(has_true)
        best_true = minimum(scores[has_true])
        gap = best_true - best
        rank = 1 + count(s -> s < best_true - score_tol, scores)
    else
        gap, rank = NaN, NaN
    end
    return (n, count(s -> s - best < score_tol, scores), count(has_true),
            gap, rank,
            100 * count(s -> s < score_truetree - score_tol, scores) / n)
end

displays_true_tree(net) =
    any(t -> rf(t, true_species_tree) == 0, displayedtrees(net, 0.0))

#-----------------------------------------------#
#   Process each replicate
#-----------------------------------------------#

rows = []
missing_reps = String[]

for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    phynestfolder = joinpath(outfolder, "rep$rep_number_string",
                             "phynestfolder")
    results_file = joinpath(phynestfolder, "phynest_results.csv")

    if !isfile(results_file)
        println("Warning: results file not found for rep $rep_number_string")
        push!(missing_reps, rep_number_string)
        continue
    end

    # one row per replicate, written by phynest_1rep.jl
    res = CSV.read(results_file, DataFrame)[1, :]
    net0 = readnewick(res.net_H0)
    net1 = readnewick(res.net_H1)

    # H=0: estimated tree vs. true species tree
    rf_net0 = compare_to_truth(net0)

    # H=1: the two displayed trees; the first one is the major tree
    displayed_trees = displayedtrees(net1, 0.0)
    if length(displayed_trees) != 2
        println("Warning: $(length(displayed_trees)) displayed tree(s) " *
            "for rep $rep_number_string, expected 2")
    end
    rf_net1_1 = compare_to_truth(displayed_trees[1])
    rf_net1_2 = length(displayed_trees) >= 2 ?
        compare_to_truth(displayed_trees[2]) : no_comparison

    # hybrid clade and donors; deepcopy because get_hybrid_info reroots
    hybrid_info = net1.numhybrids == 1 ?
        get_hybrid_info(deepcopy(net1), "A") :
        (hybrid_taxon="NA", major_donor="NA", minor_donor="NA")

    # make sure gamma_1 is the major inheritance probability
    gamma_1 = max(res.gamma_1, res.gamma_2)
    gamma_2 = min(res.gamma_1, res.gamma_2)

    # all runs of both searches: is the true tree found (H=0) / displayed
    # (H=1), and where does it rank
    tops0, scores0 = read_run_results(
        joinpath(phynestfolder, "H0_output", "H0_hc.log"))
    tops1, scores1 = read_run_results(
        joinpath(phynestfolder, "H1_output", "H1_hc.log"))
    runs_H0 = summarize_runs(scores0,
        [rf(readnewick(t), true_species_tree) == 0 for t in tops0],
        res.score_truetree)
    runs_H1 = summarize_runs(scores1,
        [displays_true_tree(readnewick(t)) for t in tops1],
        res.score_truetree)

    # alignment length M, the same for both searches of a replicate
    n_sites = read_sequence_length(
        joinpath(phynestfolder, "H0_output", "H0_hc.log"))
    n_sites == read_sequence_length(
        joinpath(phynestfolder, "H1_output", "H1_hc.log")) ||
        println("Warning: H0 and H1 alignment lengths differ " *
            "for rep $rep_number_string")
    n_sites > 0 ||
        println("Warning: sequence length not found for rep $rep_number_string")

    push!(rows, merge(
        (repID = rep_number_string,
         score_H0 = res.score_H0,
         score_H1 = res.score_H1,
         T = res.T,
         n_sites = n_sites,
         T_per_site = n_sites > 0 ? res.T / n_sites : NaN,
         score_truetree = res.score_truetree,
         num_hybrid_H1 = res.num_hybrid_H1,
         gamma_1 = gamma_1,
         gamma_2 = gamma_2),
        NamedTuple{rf_columns("RF_net0")}(rf_net0),
        NamedTuple{rf_columns("RF_net1_1")}(rf_net1_1),
        NamedTuple{rf_columns("RF_net1_2")}(rf_net1_2),
        hybrid_info,
        NamedTuple{(:n_runs_H0, :n_runs_at_best_H0, :n_runs_true_H0,
                    :gap_best_true_H0, :rank_best_true_H0,
                    :pct_runs_better_than_truetree_H0)}(runs_H0),
        NamedTuple{(:n_runs_H1, :n_runs_at_best_H1, :n_runs_true_H1,
                    :gap_best_true_H1, :rank_best_true_H1,
                    :pct_runs_better_than_truetree_H1)}(runs_H1)))

    println("Processed rep$rep_number_string: T=$(round(res.T, digits=1)) " *
        "RF_H0=$(rf_net0[1]) RF_H1=($(rf_net1_1[1]),$(rf_net1_2[1])) " *
        "hybrid=$(hybrid_info.hybrid_taxon)")
end

#-----------------------------------------------#
#   Write the per-setting summary
#-----------------------------------------------#

summary_df = DataFrame(rows)
summary_file = joinpath(outfolder, "PhyNEST-$paramname_root-summary.csv")
CSV.write(summary_file, summary_df)

println("\nWrote $summary_file")
println("Summary statistics:")
println("  Replicates processed: $(nrow(summary_df))")
println("  Replicates missing: $(length(missing_reps)) " *
    (isempty(missing_reps) ? "" : join(missing_reps, ",")))
println("  H0 tree = true tree: $(count(==(0), summary_df.RF_net0_true))")
println("  H1 displays true tree: " *
    "$(count((summary_df.RF_net1_1_true .== 0) .|
             (summary_df.RF_net1_2_true .== 0)))")
println("  T > 0: $(count(>(0), summary_df.T))")
println("  Sites per alignment: " *
    "$(minimum(summary_df.n_sites))-$(maximum(summary_df.n_sites))")
println("  Runs reaching the best score (median): " *
    "H0 $(median(summary_df.n_runs_at_best_H0)) / " *
    "$(median(summary_df.n_runs_H0)), " *
    "H1 $(median(summary_df.n_runs_at_best_H1)) / " *
    "$(median(summary_df.n_runs_H1))")
println("  Runs with the true tree (median): " *
    "H0 $(median(summary_df.n_runs_true_H0)), " *
    "H1 $(median(summary_df.n_runs_true_H1))")
println("  Rank of the best run with the true tree (median): " *
    "H0 $(median(summary_df.rank_best_true_H0)), " *
    "H1 $(median(summary_df.rank_best_true_H1))")
