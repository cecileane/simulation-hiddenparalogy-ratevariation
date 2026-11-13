
#= For one setting, count how many branches of each gene tree disagree 
with the species tree (true SimPhy gene tree) 
and how many are misestimated by
IQ-TREE (true vs estimated gene tree), by locus category. 
This gives the supplementary tables on gene tree discordance and estimation error, 
with the % of gene trees with 0, 1, 2, ... such branches. 

This script is used in one of the appendix table for better understanding 
    of what a discordant branch is and what normalized RF at certain level represent. 

The per-gene normalized RF from simulation_postprocess.jl is not enough with two
individuals per taxon: 
there the stored distance is the mean over the 8-taxon subtrees of the gene tree 
(one per choice of individual in each non-monophyletic species). 

Here the subtrees are enumerated again and each one counts 1/(number of subtrees) 
in its branch count.

Input: genetrees_stats_<setting>.csv from simulation_postprocess.jl and,
  with two individuals per taxon, the SimPhy gene trees it points to
Output: <output dir>/<setting>/genetree_tiers_<setting>.csv, one row per
  locus category (all, ortholog, falseHP, weakHP, strongHP, lossonly):
  n_genes, mean_RF_species, d0..d5 (fraction of gene trees with k
  discordant branches), mean_RF_estimated, e0..e13 (fraction with k
  misestimated branches: 5 internal branches with one individual per
  taxon, 13 with two)
Usage: julia scripts/genetree_tiers.jl --tree reptile --dup_rate 0.0003 \
  --loss_rate 0.0003 --ratevar N --n_inds 2 --SF 0.5 [--gene_len 1000]
  [--output_dir output]

the per-gene mean recomputed from the subtrees is checked against the
stored value, mismatches are printed
=# 

using ArgParse
using PhyloNetworks
using CSV, DataFrames
using Statistics

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--dup_rate"
            help = "Duplication rate of the setting"
            arg_type = Float64
            required = true
        "--loss_rate"
            help = "Loss rate of the setting"
            arg_type = Float64
            required = true
        "--ratevar"
            help = "Rate variation: N, G or L"
            arg_type = String
            required = true
        "--n_inds"
            help = "Number of individuals per species"
            arg_type = Int
            default = 1
        "--SF"
            help = "Scaling factor of the species tree branch lengths"
            arg_type = Float64
            default = 1.0
        "--gene_len"
            help = "Length of the simulated genes"
            arg_type = Int
            default = 1000
        "--tree"
            help = "Species tree used in the simulation: reptile or plant"
            arg_type = String
            required = true
        "--output_dir"
            help = "Folder holding the setting folders (default: output, " *
                   "or output_plant for the plant tree)"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

# the species tree, topology only, and the species of a tip (A_0 -> A)
species_newick(tree) = startswith(tree, "plant") ?
    "(A,(B,((C,D),(E,(H,(F,G))))));" : "(A,((((B,C),(D,E)),F),(G,H)));"
species_of(tip) = String(first(split(tip, "_")))

function category(r)
    r.nothing_or_not == "Y" && return "ortholog"
    r.false_HP == "Y" && return "falseHP"
    r.weak_HP == "Y" && return "weakHP"
    r.strong_HP == "Y" && return "strongHP"
    return "lossonly"
end

function issister(tree, a, b)
    na = tree.node[findfirst(n -> n.leaf && n.name == a, tree.node)]
    children = getchildren(getparent(na))
    return length(children) == 2 && any(n -> n.leaf && n.name == b, children)
end

"""
Input: genetree (8 taxa, tips get renamed to their species), species_tree
Output: number of discordant branches between the two, the species tree is
    pruned to the species present in the gene tree
"""
function discordant_branches(genetree, species_tree)
    for n in genetree.leaf; n.name = species_of(n.name); end
    sp = species_tree
    missing_taxa = setdiff(tipLabels(sp), tipLabels(genetree))
    if !isempty(missing_taxa)
        sp = deepcopy(sp)
        for t in missing_taxa; deleteleaf!(sp, t); end
    end
    return round(Int, mudistance_semidirected(genetree, sp, preorder=true) / 2)
end

"""
For two individuals per species
Input: file (gene tree newick), species_tree
Output: fraction of the 8-taxon subtrees of the gene tree with k = 0..5
    discordant branches (vector of length 6), and their mean normalized RF,
    as in simulation_postprocess.jl
"""
function subtree_counts(file, species_tree)
    tree = readnewick(file)
    tips = tipLabels(tree)
    species = unique(species_of.(tips))
    for s in species # merge conspecific sister tips
        ts = filter(t -> species_of(t) == s, tips)
        length(ts) == 2 && issister(tree, ts...) && deleteleaf!(tree, ts[2])
    end
    tips = tipLabels(tree)
    choices = [filter(t -> species_of(t) == s, tips) for s in species]
    w = zeros(6)
    ncombo = 0
    for combo in Iterators.product(choices...)
        sub = deepcopy(tree)
        for t in setdiff(tips, collect(combo)); deleteleaf!(sub, t); end
        w[discordant_branches(sub, species_tree) + 1] += 1
        ncombo += 1
    end
    w ./= ncombo
    return w, sum(w .* (0:5)) / 5
end

function main()
    args = parse_commandline()
    tree = args["tree"]
    tree in ("reptile", "plant") ||
        error("--tree must be reptile or plant")
    n_inds = args["n_inds"]
    paramname = "DUP$(args["dup_rate"])-LOS$(args["loss_rate"])" *
        "-RV$(args["ratevar"])-N_ind$(n_inds)-SF$(args["SF"])" *
        "-genelen$(args["gene_len"])"
    output_dir = args["output_dir"]
    if tree != "reptile" && output_dir == "output"
        output_dir = "output_$tree"
    end
    outfolder = joinpath(output_dir, paramname)
    statsfile = joinpath(outfolder, "genetrees_stats_$paramname.csv")
    isfile(statsfile) || error("missing $statsfile")
    stats = CSV.read(statsfile, DataFrame)
    species_tree = readnewick(species_newick(tree))

    # per category: 
    # n genes, sum of RF to the species tree, 
    # counts of discordant branches 0..5, 
    # sum of RF to the estimated tree, 
    # counts of misestimated branches 0..13
    acc = Dict{String,Vector{Float64}}()
    add!(c, v) = (acc[c] = get(acc, c, zeros(23)) .+ v)
    nmismatch = 0
    for r in eachrow(stats)
        v = zeros(23)
        v[1] = 1
        v[2] = r.RF_true_genetrees_vs_species_tree
        nbranch = r.num_taxa_gene_tree - 3 # internal branches of the gene tree
        if n_inds == 1
            v[3 + round(Int, v[2] * nbranch)] += 1
        else
            f = r.gene_tree_file_path
            # the stats file may point to the folder of another server
            isfile(f) ||
                (f = joinpath(outfolder, split(f, paramname * "/")[2]))
            isfile(f) || (println("missing ", f); continue)
            w, m = subtree_counts(f, species_tree)
            v[3:8] .= w
            if abs(m - v[2]) > 1e-9
                nmismatch += 1
                println("mismatch ", f, " ", m, " ", v[2])
            end
        end
        v[9] = r.RF_true_vs_estimated_genetrees
        v[10 + round(Int, v[9] * nbranch)] += 1
        add!("all", v)
        add!(category(r), v)
    end
    nmismatch > 0 && println(nmismatch, " genes with a recomputed mean RF ",
                             "different from the stored one")

    cols = ["d$k" for k in 0:5]
    append!(cols, ["e$k" for k in 0:13])
    out = DataFrame(setting = String[], category = String[], n_genes = Int[],
                    mean_RF_species = Float64[])
    for c in cols[1:6]; out[!, c] = Float64[]; end
    out[!, "mean_RF_estimated"] = Float64[]
    for c in cols[7:end]; out[!, c] = Float64[]; end
    for c in ("all", "ortholog", "falseHP", "weakHP", "strongHP", "lossonly")
        haskey(acc, c) || continue
        v = acc[c]
        n = v[1]
        push!(out, vcat(paramname, c, Int(n), v[2] / n, v[3:8] ./ n,
                        v[9] / n, v[10:23] ./ n))
    end
    outfile = joinpath(outfolder, "genetree_tiers_$paramname.csv")
    CSV.write(outfile, out)
    println("written ", outfile)
end

main()
