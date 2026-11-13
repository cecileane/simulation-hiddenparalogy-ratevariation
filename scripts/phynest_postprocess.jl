
#= Postprocessing for PhyNEST, one parameter set at a time.
Reads the per-replicate results written by phynest_1rep.jl 
then writes one summary csv per parameter set 
run_postprocessing.jl --mode phynest calls this for all parameter sets 
and copies the csv to phynest_summary/

What gets added here that phynest_1rep.jl does not record:
 - RF distances of the H=0 tree and of the two trees displayed in the H=1 network
    to the true species tree (full, and with one taxon pruned) 
and to the three alternative trees,
    same alternatives as snaq_postprocess.jl
 - hybrid clade and its donors in the H=1 network
 - from the search logs: how many runs reached the best score, 
how many found (H=0) or displayed (H=1)
    the true tree, score gap and rank of the best such run
 - number of sites M (from the PhyNEST log) and T_per_site = T / M
    T grows linearly with M, and M is smaller when fewer loci were retained

Usage: julia scripts/phynest_postprocess.jl
       --dup_rate 0.0 --loss_rate 0.0 --ratevar N
       --n_inds 1 --SF 0.5 --n_reps 100
       --tree reptile (default) or plant: picks the true tree, the alternative trees
       and the default output folder (output/ or output_plant/)

Note: no threshold on T here, because 
model selection (T > 0 or T > T*) is done later in the summary scripts,
so T* can change without re-running this.
RF distances: mudistance_semidirected, same as for SNaQ and find_graphs
=# 
using Pkg
repo_env = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(repo_env; io=devnull)
Pkg.instantiate(; io=devnull)

using ArgParse
using CSV
using DataFrames
using PhyloNetworks
include("utilities.jl") # for pad_number, set_up_paramname_root, get_hybrid_info

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        # parameter settings, to find the right output folder
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
        "--tree"
            help = "Species tree used in the simulation: reptile or plant " *
                   "(selects the reference trees and the default output folder)"
            arg_type = String
            default = "reptile"
        "--output_dir"
            help = "Directory holding the parameter folders (default = " *
                   "output, or output_plant for the plant tree)"
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

# get the arguments
dup_rate = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"]
n_reps = parsed_args["n_reps"]
n_inds = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
SF = parsed_args["SF"]
gene_len = parsed_args["gene_len"]

tree = parsed_args["tree"]
tree in ("reptile", "plant") ||
  error("--tree must be reptile or plant")
output_dir = parsed_args["output_dir"]
# plant runs are in output_plant/ unless --output_dir says otherwise

if tree != "reptile" && output_dir == "output"
    output_dir = "output_$tree"
end

# folder of this parameter set
paramname_root = set_up_paramname_root(dup_rate, loss_rate,
                                    ratevar, n_inds, SF,
                                    gene_len)
outfolder = joinpath(output_dir, paramname_root)

if !isdir(outfolder)
    error("Output folder does not exist: $outfolder")
end

println("PhyNEST Postprocessing Script")
println("Species tree: $tree")
println("Parameter set: $paramname_root")
println("Processing replicates $rep_start to $rep_end")
println("Output folder: $outfolder")

#-----------------------------------------------#
#   Reference topologies
#-----------------------------------------------#

#= true tree, alternative trees and pruned taxa: 
    same as in snaq_postprocess.jl
   Idea: if the estimated tree is not the species tree, 
    was the taxon on the short branches (or some suspecious taxa) just moved?
   reptile tree, taxon F:
       alter1 = F with (B,C), alter2 = F with (D,E), 
       alter3 = F outside of (B,C),(D,E),(G,H)
       pruned taxa: F, G, H
   plant tree, taxon B:
       alter1 = B with (C,D); alter2 = B with the E clade; alter3 = H moved next to G
       pruned taxa: B, F, G
   B is the hard one in the plant tree. With lineage rate variation, 
   the wrong findgraphs trees are alter2 (B and (C,D) exchanged) 
   but the wrong PhyNEST H=0 trees are alter2 + an E/H swap,
   so pruning B fixes the findgraphs trees but not the PhyNEST trees 
   (Note: might want to double check this if reviewers want to know the inferred topology of phynest).
   Column names keep the reptile suffixes noF, noG, noH. 
   For the plant tree they mean without B, F, G 
=#

if startswith(tree, "plant")
    true_species_tree = readnewick("(A,(B,((C,D),(E,(H,(F,G))))));")
    alter1 = readnewick("(A,((B,(C,D)),(E,(H,(F,G)))));")
    alter2 = readnewick("(A,((C,D),(B,(E,(H,(F,G))))));")
    alter3 = readnewick("(A,(B,((C,D),(E,(F,(G,H))))));")
    pruned_taxa = ["B", "F", "G"]
else
    true_species_tree = readnewick("(A,((((B,C),(D,E)),F),(G,H)));")
    alter1 = readnewick("(A,((G,H),(((B,C),F),(D,E))));")
    alter2 = readnewick("(A,((G,H),((B,C),(F,(D,E)))));")
    alter3 = readnewick("(A,(((G,H),((D,E),(B,C))),F));")
    pruned_taxa = ["F", "G", "H"]
end

# true tree with one taxon pruned, one per taxon in pruned_taxa
species_tree_pruned = Dict{String,HybridNetwork}()
for taxon in pruned_taxa
    species_tree_pruned[taxon] = deepcopy(true_species_tree)
    deleteleaf!(species_tree_pruned[taxon], taxon)
end

# column suffixes, same order as what compare_to_truth returns
rf_names = ["true", "true_noF", "true_noG", "true_noH",
            "alter1", "alter2", "alter3"]
rf_columns(prefix) = Tuple(Symbol.(prefix * "_" .* rf_names))

rf(tree, ref) = mudistance_semidirected(tree, ref, preorder=true)

"""
Input:
-tree (HybridNetwork): estimated tree, or a tree displayed in the estimated network
-taxon (String): tip to prune, one of pruned_taxa
Output: RF distance between tree without taxon and the true species tree without taxon
Notes: deleteleaf! can fail on some trees, then we return NaN and move on
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
All the topology comparisons for one tree
Input: tree (HybridNetwork)
Output: Tuple of 7 RF distances, in the order of rf_names:
to the true tree, 
to the true tree without each of the 3 pruned taxa, 
to alter1, alter2, alter3
"""
function compare_to_truth(tree::HybridNetwork)
    return (rf(tree, true_species_tree),
            rf_pruned(tree, pruned_taxa[1]), rf_pruned(tree, pruned_taxa[2]),
            rf_pruned(tree, pruned_taxa[3]),
            rf(tree, alter1), rf(tree, alter2), rf(tree, alter3))
end

no_comparison = Tuple(fill(NaN, length(rf_names)))

"""
Read the hill climbing log of PhyNEST (H0_hc.log or H1_hc.log)
Input: logfile (String): path to the log files 
Output: Tuple (topologies, scores)
    -topologies (Vector{String}): estimated topology (newick) of each run
    -scores (Vector{Float64}): composite likelihood of each run, same order
Notes: runs where there is a failure with NaN in the log, those are dropped
    if the log is missing, returns two empty vectors. This is okay and tested. 
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
Number of sites M in the alignment that PhyNEST analyzed
Input: logfile (String): path to H0_hc.log or H1_hc.log
Output: Int, from the line "Sequence length: 951000" in the log header. 
    0 if the line is not there
"""
function read_sequence_length(logfile::String)
    isfile(logfile) || return 0
    for line in eachline(logfile)
        m = match(r"^Sequence length: (\d+)", line)
        m === nothing || return parse(Int, m[1])
    end
    return 0
end

score_tol = 1.0 # scores within 1 unit = tie, that's optimizer noise and not the topology

"""
Summarize the runs of one search (H=0 or H=1)
Input:
    -scores (Vector{Float64}): 
    composite likelihood of each run, from read_run_results
    -has_true (Vector{Bool}): did this run find (H=0) or display (H=1) the true tree or not 
    -score_truetree (Float64): score of the true tree, from phynest_results.csv
Output: Tuple (n_runs, n_at_best, n_true, gap, rank, pct_better)
    -n_at_best: 
        number of runs within score_tol of the best score
    -gap: 
        best score among the runs with the true tree
    -rank: 
        rank of the best run with the true tree, 1 = best
    -pct_better: 
        % of runs that scored better than the true tree itself
    gap and rank are NaN if no run has the true tree
"""
function summarize_runs(scores::Vector{Float64}, 
        has_true::Vector{Bool},
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



#---------------------------------------------#
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

    # phynest_results.csv has one row, which is written by phynest_1rep.jl
    res = CSV.read(results_file, DataFrame)[1, :]
    net0 = readnewick(res.net_H0)
    net1 = readnewick(res.net_H1)

    # H=0: estimated tree vs true species tree
    rf_net0 = compare_to_truth(net0)

    # H=1: the two displayed trees, the first one is the major tree
    displayed_trees = displayedtrees(net1, 0.0)
    if length(displayed_trees) != 2
        println("Warning: $(length(displayed_trees)) displayed tree(s) " *
            "for rep $rep_number_string, expected 2")
    end
    rf_net1_1 = compare_to_truth(displayed_trees[1])
    rf_net1_2 = length(displayed_trees) >= 2 ?
        compare_to_truth(displayed_trees[2]) : no_comparison

    # hybrid clade and donors. deepcopy because get_hybrid_info reroots the network
    hybrid_info = net1.numhybrids == 1 ?
        get_hybrid_info(deepcopy(net1), "A") :
        (hybrid_taxon="NA", major_donor="NA", minor_donor="NA")

    # gamma_1 should be the major one
    gamma_1 = max(res.gamma_1, res.gamma_2)
    gamma_2 = min(res.gamma_1, res.gamma_2)

    # go through all runs of both searches: 
    # which runs found (H=0) or displayed (H=1) the true tree, and their rank
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

    # alignment length M, should be the same in H0 and H1 logs
    n_sites = read_sequence_length(
        joinpath(phynestfolder, "H0_output", "H0_hc.log"))
    n_sites == read_sequence_length(
        joinpath(phynestfolder, "H1_output", "H1_hc.log")) ||
        println("Warning: H0 and H1 alignment lengths differ " *
            "for rep $rep_number_string")

    n_sites > 0 ||
        println("Warning: sequence length not found for rep $rep_number_string")

    # Get all results together: 
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

end

#-----------------------------------------------#
#   Write the per-setting summary
#-----------------------------------------------#
summary_df = DataFrame(rows)
summary_file = joinpath(outfolder, "PhyNEST-$paramname_root-summary.csv")
CSV.write(summary_file, summary_df)

# Print a lot of stuff to potentially see if there could be an issue: 
println("Replicates processed: $(nrow(summary_df))")
println("H0 tree = true tree: $(count(==(0), summary_df.RF_net0_true))")
println("H1 displays true tree: " *
"$(count((summary_df.RF_net1_1_true .== 0) .| (summary_df.RF_net1_2_true .== 0)))")
println("T > 0: $(count(>(0), summary_df.T))")
println("Sites per alignment: " * "$(minimum(summary_df.n_sites))-$(maximum(summary_df.n_sites))")

# When running, could be good to save in a log file and might check later... 
println("Runs reaching the best score (median): " *
    "H0 $(median(summary_df.n_runs_at_best_H0)) / " *
     "$(median(summary_df.n_runs_H0)), " *
    "H1 $(median(summary_df.n_runs_at_best_H1)) / " * "$(median(summary_df.n_runs_H1))")
println("Runs with the true tree (median): " *
        "H0 $(median(summary_df.n_runs_true_H0)), " *
            "H1 $(median(summary_df.n_runs_true_H1))")
println("Rank of the best run with the true tree (median): " *
    "H0 $(median(summary_df.rank_best_true_H0)), " *
        "H1 $(median(summary_df.rank_best_true_H1))")
