#=
code used to build a second "true" species tree, from angiosperms,
used by SimPhy to simulate gene trees.
The idea here is that angiosperms that diverge earlier
are known to have lots of hidden paralogs,
so we want to see how that affects the gene tree distributions.

Same procedure as speciestree_reptile.jl
(reptile tree, Crawford et al. 2012 UCE data), applied to the Angiosperms353
genes of the Kew Tree of Life Explorer, release 4.0 (Baker et al. 2022;
Zuntini et al. 2024), data under CC BY 4.0 at
https://sftp.kew.org/pub/paftol/current_release
Previsou trials tried 1KP but it retrieved too few genes for the 8 taxa.

Here, branch lengths of the species tree are in coalescent units (CU);
per-branch multipliers give the lineage-specific substitution rates.

Run from the repo root (each step is skipped if its output exists):
  julia --project=. scripts/speciestree/speciestree_angiosperm.jl
needs: executables/iqtree2 (v2.4.0), executables/Astral/astral.5.7.8.jar,
Rscript with ape (for scripts/speciestree/calibrate.r)
These set of softwares are compatible with the rest of the repo.

taxa: 1 gymnosperm outgroup + 7 basal angiosperms (Zuntini et al. 2024 without
monocots and eudicots). Kew "Sequence_ID" of each sample in parentheses.
  A = Pinus ponderosa (ERR2040837) -> outgroup
  B = Amborella trichopoda (GCF_000471905.2)
  C = Nymphaea nouchali var. zanzibariensis (GCF_008831285.2)
      = the N. colorata genome
  D = Brasenia schreberi (GCA_030142015.1)
  E = Schisandra chinensis (SRR26404306)
  F = Chloranthus sessilifolius (GCA_021018995.1)
  G = Liriodendron chinense (GCA_003013855.2)
  H = Ceratophyllum demersum var. quadrispinum (SRR6374676)

data: Kew's published per-gene alignments
(fasta/alignments/<gene>.dna.aln.fasta).
Note:
1. Our samples are got as aligned;
2. columns that are gaps in all of them are removed;
3. only genes with all 8 taxa are used (260 of 353).
We can only use 260 genes because the other 93 are missing one or more of
our samples. If the gene trees miss one or more taxa, we cannot accurately
estimate the median tree height (see below), so we decided to remove those
genes from the analysis.

gene trees: IQ-TREE with model chosen by ModelFinder for each gene,
1000 ultrafast bootstraps, then ASTRAL with the bootstrap trees,
same as the reptile tree.

ASTRAL estimates branch lengths in CU for internal branches only. With one
sample per species, the external branches and the ingroup stem cannot be
estimated. They are calibrated against substitutions per site, which we
estimate on every branch from the gene trees (step 4): CU = k * substitutions
with one ratio k = 22, from the stem of the water lilies (Nymphaea + Brasenia).
This is done in scripts/speciestree/calibrate.r (step 5), which gives the
water-lily tips, the ingroup stem from the path Pinus -> ingroup crown, and
the other tips from ultrametricity. The same tree is rebuilt here with
ultrametrize! from QuartetNetworkGoodnessFit, as a check.

substitution rates as in ratevariation.jl (reptiles repo):
1. locus height d_l = median distance from the outgroup to the ingroup taxa
    in gene tree l; This is why we need to remove genes with missing taxa.
2. log-normal fit of d_l/mean(d_l) gives -hl;
3. gene trees rescaled by mean(d_l)/d_l, pairwise distances averaged
across genes and fitted to the species tree by least squares give d_i;
4. bar_r = sum(d_i)/sum(tau_i),
5. -su = bar_r/2Ne,
6. m_i = d_i/(tau_i bar_r).

2Ne needs to be decided by the overall CU. See below.

=#

using PhyloNetworks
using QuartetNetworkGoodnessFit # for ultrametrize!
using PhyloPlots, RCall
using CSV, DataFrames
using Distributions, Statistics
using Downloads

datadir = "data/kew_a353" # run from the repo root
mkpath(datadir)

#= Simphy takes the number of generations as input,
so we need to decide on a value for 2Ne.
If the tree have the same or approxiamately the same number of generations
as the reptile tree,
then we don't need to re-calibrate the duplication and loss rates.
However, since the angiosperm tree is much deeper than the reptile tree,
we need to use a smaller 2Ne to keep the same number of generations.
=#
# The number of generations of the reptile tree is 3440 generations in total
# so here the eff_pop for the angiosperm tree should be
# 3440 / CU length of the tree

# height of the reptile tree in generations: 3.44 CU (Homo branch) x 2Ne = 1000
reptile_ngen = 3440

# 2Ne is set in step 6: 2Ne = 3440 / (height of the tree in CU), rounded to the
# nearest 50, so that the height in generations is about the same as for the
# reptile tree and the duplication/loss rates per generation can be reused.
# This 2Ne is passed to simulation.jl with --Ne (SimPhy -sp), together with
# -su.
kew = "https://sftp.kew.org/pub/paftol/current_release"
iqtree = abspath("executables/iqtree2")
astral = abspath("executables/Astral/astral.5.7.8.jar")

taxa = DataFrame(
    letter = ["A","B","C","D","E","F","G","H"],
    code = ["Pinus","Amborella","Nymphaea","Brasenia",
            "Schisandra","Chloranthus","Liriodendron","Ceratophyllum"],
    sample = ["ERR2040837","GCF_000471905.2","GCF_008831285.2",
                "GCA_030142015.1","SRR26404306",
              "GCA_021018995.1","GCA_003013855.2","SRR6374676"])
codes = taxa.code
code2letter = Dict(taxa.code .=> taxa.letter)
sample2code = Dict(taxa.sample .=> taxa.code)
CSV.write(joinpath(datadir, "taxa.csv"), taxa)

#--------------------------------------------------------------------#
# 1. data: our 8 samples from Kew's alignments
# use genes with all 8 taxa
#--------------------------------------------------------------------#
alndir = joinpath(datadir, "alignments")
if !isdir(alndir)
    mkpath(alndir)
    listing = Downloads.download(kew * "/fasta/alignments/", IOBuffer())
    listing = String(take!(listing))
    genes = unique(m.captures[1] for m in
                   eachmatch(r"href=\"(\d+)\.dna\.aln\.fasta\"", listing))
    println(length(genes), " gene alignments at Kew") # 353
    for g in genes
        f = Downloads.download("$kew/fasta/alignments/$g.dna.aln.fasta")
        seqs = Dict{String,String}(); name = ""
        # keep the sequences of our samples: "Sequence_ID:xxx" in the header
        for line in eachline(f)
            if startswith(line, ">")
                m = match(r"Sequence_ID:(\S+)", line)
                name = (m !== nothing && haskey(sample2code, m.captures[1])) ?
                       sample2code[m.captures[1]] : ""
                name == "" || (seqs[name] = "")
            elseif name != ""
                seqs[name] *= uppercase(strip(line))
            end
        end
        rm(f)
        # genes with all 8 taxa only
        all(haskey(seqs, c) for c in codes) || continue
        # remove the columns that are gaps in all 8 samples
        L = length(seqs[codes[1]])
        keep = [j for j in 1:L if
                any(seqs[c][j] ∉ ('-', '?', 'N') for c in codes)]
        open(joinpath(alndir, "$g.fasta"), "w") do io
            for c in codes
                println(io, ">", c); println(io, seqs[c][keep])
            end
        end
    end
end
println("genes with all 8 taxa: ", length(readdir(alndir))) # 260

#--------------------------------------------------------------------#
# 2. gene trees (IQ-TREE), bootstrap trees split by gene, ASTRAL
#--------------------------------------------------------------------#
iqdir = joinpath(datadir, "iqtree")
if !isfile(joinpath(iqdir, "loci.treefile"))
    mkpath(iqdir) # iqtree does not create it
    cd(datadir) do
        run(`$iqtree -S alignments --prefix iqtree/loci -T 8
             -B 1000 -wbtl -seed 1`)
    end
end
# loci in the order used by IQ-TREE, from the charset lines
locusnames = String[]
for line in readlines(joinpath(iqdir, "loci.best_scheme.nex"))
    m = match(r"charset (\S+) =", line)
    m === nothing || push!(locusnames, m.captures[1])
end
# 1 bootstrap file per gene (as separate-boot-bygene.jl in the reptiles repo)
# and their list for ASTRAL
if !isfile(joinpath(datadir, "BSlistfiles"))
    # 1000 trees per gene, in the same order
    boot = readlines(joinpath(iqdir, "loci.ufboot"))
    mkpath(joinpath(iqdir, "bootstrap"))
    open(joinpath(datadir, "BSlistfiles"), "w") do list
        for (i, name) in enumerate(locusnames)
            open(joinpath(iqdir, "bootstrap", "$name.ufboot"), "w") do io
                for t in boot[(i-1)*1000+1 : i*1000]; println(io, t); end
            end
            println(list, "iqtree/bootstrap/$name.ufboot")
        end
    end
end
astralfile = joinpath(datadir, "astral", "species.tre")
if !isfile(astralfile)
    mkpath(joinpath(datadir, "astral"))
    cd(datadir) do
        run(`java -jar $astral -i iqtree/loci.treefile -b BSlistfiles
             -r 1000 -o astral/species.tre`)
    end
end
# species.tre: 1000 bootstrap trees, the ASTRAL tree with support,
# then the same tree with lengths in CU
speciestree_string = readlines(astralfile)[end]
println("ASTRAL tree:\n", speciestree_string)
#= we get:
(Pinus,(((Nymphaea,Brasenia)100.0:3.7727609380946374,(Schisandra,(Ceratophyllum,(Chloranthus,Liriodendron)100.0:0.44234373849353004)100.0:0.9469277013360451)100.0:0.7874999644079429)90.0:0.15847034097177473,Amborella):0.0);
same topology as the Kew 4.0 species tree pruned to these 8 samples.
=#

#--------------------------------------------------------------------#
# 3. tree in coalescent units: internal branches from ASTRAL
#--------------------------------------------------------------------#
tree = readnewick(speciestree_string)
# remove bootstrap values from node names (SimPhy doesn't accept them)
for n in tree.node
    isleaf(n) || (n.name = "")
end
rootatnode!(tree, "Pinus")
for e in tree.edge # round to 2 digits: not significant for simulations
    e.length == -1.0 || (e.length = round(e.length, digits=2))
end
# the tree lacks external edge lengths, and the ingroup stem is 0.0
# (not estimable), as for the reptile tree
stem = tree.edge[findfirst(e -> getparent(e) === getroot(tree) &&
                                !getchild(e).leaf, tree.edge)]
stem.length = -1.0 # unknown, like the external edges
println(writenewick(tree))
# (Pinus,(((Nymphaea,Brasenia):3.77,(Schisandra,(Ceratophyllum,(Chloranthus,Liriodendron):0.44):0.95):0.79):0.16,Amborella));
# distance from the root to a node, and the node of a given tip: used below
function depth(node, tree)
    d = 0.0
    while node !== getroot(tree)
        d += getparentedge(node).length
        node = getparent(node)
    end
    return d
end
tipnode(name, tree) = tree.node[findfirst(x -> x.name == name, tree.node)]

#--------------------------------------------------------------------#
# 4. substitution rates: across genes (-hl) and per branch (d_i)
#--------------------------------------------------------------------#
genetrees = readmultinewick(joinpath(iqdir, "loci.treefile"))
# pairwise distances of a gene tree, taxa in the order of `codes` (Pinus first)
function distances(gt)
    D = pairwisetaxondistancematrix(gt)
    order = [findfirst(==(c), tiplabels(gt)) for c in codes]
    return D[order, order]
end
# locus height d_l: median distance from the outgroup to the ingroup taxa
d_l = [median(distances(gt)[1, 2:end]) for gt in genetrees]
mean_d = mean(d_l)
println("locus heights d_l: mean ", round(mean_d, digits=4),
        ", range ", round(minimum(d_l), digits=3),
        " to ", round(maximum(d_l), digits=3))
# rate variation across genes: log-normal fitted to d_l / mean(d_l),
# then mean 1 imposed
rel = d_l ./ mean_d
for family in [Normal, LogNormal, Gamma]
    d = fit(family, rel)
    println("loglik=$(round(sum(logpdf.(d, rel)), digits=1)), $d")
end
sigma_gene = fit(LogNormal, rel).σ
# mu so that the multiplier has mean exp(mu + sigma^2/2) = 1
mu_gene = -sigma_gene^2 / 2
hl = "ln:$(round(mu_gene, digits=6)),$(round(sigma_gene, digits=6))"
println("-hl ", hl) # -hl ln:-0.03634,0.269593
# in reptile: -hl ln:-0.19,0.6164414002968976
# we get: d_l mean 0.7173 (range 0.376 to 3.234),
# log-normal best (loglik -17.1 vs -33.6 gamma, -92.8 normal)
# -hl ln:-0.03634,0.269593 : much less variation across genes than for the
# reptile UCEs
# substitutions per site on each branch: rescale each gene tree by
# mean(d_l)/d_l, average the pairwise distances, fit them to the species tree
# by least squares
D_mean = zeros(length(codes), length(codes))
for (gt, h) in zip(genetrees, d_l)
    global D_mean .+= distances(gt) .* (mean_d / h)
end
D_mean ./= length(genetrees)
sub_tree = deepcopy(tree) # same topology and edge order as `tree`
for e in sub_tree.edge; e.length = 0.1; end # starting values
# no clock: lineage rates differ
calibratefrompairwisedistances!(sub_tree, D_mean, codes; ultrametric=false)
println("substitutions per site:\n", writenewick(sub_tree))
# (Pinus:0.2272,(((Nymphaea:0.1314,Brasenia:0.1359):0.1721,(Schisandra:0.1671,(Ceratophyllum:0.3200,(Chloranthus:0.1473,Liriodendron:0.1345):0.0224):0.0330):0.0311):0.0065,Amborella:0.2536):0.2272);
# note: the outgroup branch and the stem get the same length: only their sum
# is identifiable from pairwise distances (same for the reptile tree:
# Homo 0.0347, stem 0.0347)
# These values (3 digits) are the ones typed into calibrate.r:
for (ecu, esub) in zip(tree.edge, sub_tree.edge)
    println(getchild(ecu).leaf ? getchild(ecu).name : "internal", "\t",
            ecu.length, "\t", round(esub.length, digits=4))
end

#--------------------------------------------------------------------#
# 5. external branches and stem in CU: calibrate.r,
#    then check with ultrametrize!
#--------------------------------------------------------------------#
# CU / substitutions on the 5 internal branches: 19.6, 28.8, 25.4, 24.6, 21.9.
# calibrate.r uses k = 22 (the water-lily stem: 3.77 / 0.172), sets the
# water-lily tips to k * their substitutions, the stem from the path
# Pinus -> ingroup crown (k * 0.454 = 2 * stem + ingroup height), the other
# tips from ultrametricity. It writes the tree to
# calibration/speciestree_cu_from_subs.tre
calfile = joinpath(datadir, "calibration", "speciestree_cu_from_subs.tre")
run(`Rscript scripts/speciestree/calibrate.r`)
tree_r = readnewick(readline(calfile))
println("tree in CU from calibrate.r:\n", writenewick(tree_r))

# same thing with PhyloNetworks, as in speciestree_reptile.jl: assign the
# water-lily tips and the stem, then ultrametrize! assigns the other
# external edges
k = 22
# water-lily stem, 3.77 CU, and its substitutions
e_nb = getparentedge(getparent(tipnode("Nymphaea", tree)))
sub_nb = sub_tree.edge[findfirst(e -> e === e_nb, tree.edge)].length
println("k = CU / substitutions on the water-lily stem: ",
        round(e_nb.length / sub_nb, digits=2)) # 21.9
subs(name) = sub_tree.edge[findfirst(e -> getchild(e).name == name,
                                     sub_tree.edge)].length
# 2.94, same for both tips
cu_nbtip = k * (subs("Nymphaea") + subs("Brasenia")) / 2
for name in ("Nymphaea", "Brasenia")
    getparentedge(tipnode(name, tree)).length = cu_nbtip
end
# angiosperm crown (child of the stem) to the tips: 0.16 + 3.77 + 2.94 = 6.87
H_ingroup = depth(tipnode("Nymphaea", tree), tree) -
            depth(getchild(stem), tree)
# path Pinus tip -> ingroup crown: Pinus + stem = 0.454 subs = k * 0.454 CU,
# and the tree being ultrametric, Pinus = stem + H_ingroup.
# So 2 * stem + H_ingroup = k * 0.454
sub_stem = sub_tree.edge[findfirst(e -> e === stem, tree.edge)].length
stem.length = (k * (subs("Pinus") + sub_stem) - H_ingroup) / 2 # 1.56
# verbose: lists the edges assigned
QuartetNetworkGoodnessFit.ultrametrize!(tree, true)
for e in tree.edge; e.length = round(e.length, digits=2); end
println(writenewick(tree))
# (Pinus:8.43,(((Nymphaea:2.94,Brasenia:2.94):3.77,(Schisandra:5.92,(Ceratophyllum:4.97,(Chloranthus:4.53,Liriodendron:4.53):0.44):0.95):0.79):0.16,Amborella:6.87):1.56);

hardwiredclusterdistance(tree, tree_r, true) == 0 # are they the same? Yes!
for e in tree.edge # same edge lengths?
    er = tree_r.edge[findfirst(x -> hardwiredcluster(x, codes) ==
                                    hardwiredcluster(e, codes), tree_r.edge)]
    abs(e.length - er.length) < 0.005 ||
        error("edge lengths differ from calibrate.r: ",
              "$(e.length) vs $(er.length)")
end

#= we get, from both:
(Pinus:8.43,(((Nymphaea:2.94,Brasenia:2.94):3.77,(Schisandra:5.92,(Ceratophyllum:4.97,(Chloranthus:4.53,Liriodendron:4.53):0.44):0.95):0.79):0.16,Amborella:6.87):1.56);
height 8.43 CU, stem 1.56 CU
=#

#--------------------------------------------------------------------#
# 6. 2Ne, genome-wide rate, multipliers, SimPhy tree
#--------------------------------------------------------------------#
# substitutions on the stem and the Pinus branch: only their sum is known
# (0.454). Same rate on both (same k), so split in proportion to CU:
e_pinus = getparentedge(tipnode("Pinus", tree))
i_stem = findfirst(e -> e === stem, tree.edge)
i_pinus = findfirst(e -> e === e_pinus, tree.edge)
sub_path = sub_tree.edge[i_stem].length + sub_tree.edge[i_pinus].length
cu_path = stem.length + e_pinus.length
sub_tree.edge[i_stem].length  = sub_path * stem.length / cu_path    # 0.071
sub_tree.edge[i_pinus].length = sub_path * e_pinus.length / cu_path # 0.383

for n in tree.node; n.leaf && (n.name = code2letter[n.name]); end
for n in sub_tree.node; n.leaf && (n.name = code2letter[n.name]); end
cu_total = sum(e.length for e in tree.edge)
sub_total = sum(e.length for e in sub_tree.edge)
height_cu = depth(tipnode("A", tree), tree) # root to tip (Pinus), in CU: 8.43
# 2Ne: 3440 / 8.43 = 408 -> 400
eff_pop = 50 * round(Int, reptile_ngen / height_cu / 50)
bar_r = sub_total / cu_total # substitutions per site per coalescent unit
su = bar_r / eff_pop         # per site per generation, SimPhy -su
println("height ", height_cu, " CU, total ", round(cu_total, digits=2),
        " CU; 2Ne = ", eff_pop, ": height ", round(Int, height_cu * eff_pop),
        " generations (reptile: 3440), sum of all branches ",
        round(Int, cu_total * eff_pop), " generations (reptile: 17480)")
println("bar_r = ", round(bar_r, digits=6),
        " substitutions/site/CU, -su f:", su)

# tree with edge lengths in number of generations: gen = cu * 2Ne
tree_ngen = deepcopy(tree)
for e in tree_ngen.edge; e.length = round(e.length * eff_pop); end
# per-branch multiplier m_i = d_i / (tau_i bar_r), unit-less,
# as in speciestree_reptile.jl
tree_mult = deepcopy(tree)
for (e, ecu, esub) in zip(tree_mult.edge, tree.edge, sub_tree.edge)
    e.length = round(esub.length / (ecu.length * bar_r), digits=6)
end
println("branch\tCU\tsubs/site\tm_i")
for (ecu, esub, em) in zip(tree.edge, sub_tree.edge, tree_mult.edge)
    println(getchild(ecu).leaf ? getchild(ecu).name : "internal", "\t",
            ecu.length, "\t", round(esub.length, digits=4), "\t", em.length)
end
# SimPhy tree with "generations*multiplier" on each edge: the two newick
# strings (generations, multipliers) have the same structure, so their
# numbers are paired
function simphynewick(tree_gen, tree_mult)
    g = split(writenewick(tree_gen), ":")
    m = split(writenewick(tree_mult), ":")
    s = g[1]
    for i in 2:length(g)
        # same regex for both: lengths may be written as 1.0e8
        ng = match(r"^[0-9.eE+-]+", g[i]).match
        nm = match(r"^[0-9.eE+-]+", m[i]).match
        s *= ":" * string(round(Int, parse(Float64, ng))) * "*" * nm *
             g[i][length(ng)+1:end]
    end
    return s
end
simphy_flat = replace(writenewick(tree_ngen), ".0" => "")
simphy_mult = simphynewick(tree_ngen, tree_mult)
println("SimPhy tree, no lineage variation:\n", simphy_flat)
println("SimPhy tree with multipliers:\n", simphy_mult)
open(joinpath(datadir, "simphy_tree_angiosperm.txt"), "w") do io
    println(io, "# taxa: ",
            join(["$(r.letter)=$(r.code)" for r in eachrow(taxa)], ", "))
    println(io, "# 2Ne = ", eff_pop, " (--Ne in simulation.jl), so that the ",
            "height is ", round(Int, height_cu * eff_pop),
            " generations, about as for the reptile tree (3440)")
    println(io, "# -su f:", su)
    println(io, "# -hl ", hl)
    println(io, "# species tree, coalescent units:\n", writenewick(tree))
    println(io, "# species tree, substitutions/site:\n", writenewick(sub_tree))
    println(io, "# SimPhy species tree, no lineage rate variation ",
            "(generations, 2Ne=$eff_pop):\n", simphy_flat)
    println(io, "# SimPhy species tree with lineage-specific multipliers:\n",
            simphy_mult)
end

#= we get:
height 8.43 CU, total 48.8 CU; 2Ne = 400: height 3372 generations
(reptile: 3440), sum of all branches 19520 (reptile: 17480)
bar_r = 0.041172 substitutions/site/CU, -su f:0.00010293039257500044
(A:3372,(((C:1176,D:1176):1508,(E:2368,(H:1988,(F:1812,G:1812):176):380):316):64,B:2748):624);
(A:3372*1.104565,(((C:1176*1.085755,D:1176*1.12242):1508*1.108897,(E:2368*0.685743,(H:1988*1.564038,(F:1812*0.789595,G:1812*0.721226):176*1.233907):380*0.842807):316*0.954714):64*0.993533,B:2748*0.896529):624*1.104565);
-su is larger than for the reptile tree (5.6e-5 in the paper) because the same
number of substitutions accumulates over about the same number of generations
but on a tree with a smaller total length (2Ne is smaller).
The multipliers are all between 0.69 and 1.56: the stem and Pinus get the same
one (1.10) because their substitutions were split in proportion to CU.
=#

#--------------------------------------------------------------------#
# 7. parameters for seq-gen: from the loci for which ModelFinder chose HKY
#    (same as notes/choice-seqgen-parameters.md for the reptile loci)
#--------------------------------------------------------------------#
models = readlines(joinpath(iqdir, "loci.best_model.nex"))
hky = filter(l -> occursin("HKY", l), models)
kappa = [parse(Float64, m.captures[1]) for l in hky
         for m in eachmatch(r"HKY\{([0-9.]+)\}", l)]
alpha = [parse(Float64, m.captures[1]) for l in hky
         for m in eachmatch(r"G4\{([0-9.]+)\}", l)]
freqs = [parse.(Float64, split(m.captures[1], ",")) for l in hky
         for m in eachmatch(r"F\{([0-9.,]+)\}", l)]
println("\nloci with HKY: ", length(kappa), " of ", length(genetrees))
println("kappa: ", fit(LogNormal, kappa),
        " (mean ", round(mean(kappa), digits=3), ")")
println("gamma shape: ", fit(LogNormal, alpha),
        " (mean ", round(mean(alpha), digits=3), ")")
println("base frequencies: ", fit(Dirichlet, reduce(hcat, freqs)))
#= we get, from 27 HKY loci (TPM2u+F+G4 chosen for 52 loci, TIM2+F+G4 for 26):
kappa: LogNormal(1.3330, 0.2079), mean 3.876
gamma shape: LogNormal(-0.4045, 0.6540), mean 0.833
base frequencies: Dirichlet([109.99, 75.31, 88.96, 106.00])
=#

#--------------------------------------------------------------------#
# 8. plot: left: the tree in coalescent units, "CU (x multiplier)" on each
#    edge; right: branch lengths in substitutions per site as simulated,
#    i.e. generations x su x lineage multiplier (= d_i)
#--------------------------------------------------------------------#
tree_subs = deepcopy(tree) # as simulated
for (e, eg, em) in zip(tree_subs.edge, tree_ngen.edge, tree_mult.edge)
    e.length = eg.length * su * em.length
end
letter2code = Dict(taxa.letter .=> taxa.code)
for t in (tree, tree_subs)
    for n in t.node; n.leaf && (n.name = letter2code[n.name]); end
end
R"png"(joinpath(datadir, "speciestree_angiosperm.png"),
       width=2400, height=1000, res=110)
R"layout"([1 2])
R"par(mar=c(1,1,3,1), cex.main=1.6)"
lab = DataFrame(number = [e.number for e in tree.edge],
                label = ["$(e.length) (x$(round(em.length, digits=2)))"
                         for (e, em) in zip(tree.edge, tree_mult.edge)])
plot(tree, useedgelength=true, edgelabel=lab, tipcex=1.6, edgecex=1.4,
     xlim=[0, 13])
R"title"("coalescent units (x lineage multiplier), height $height_cu, " *
         "2Ne = $eff_pop")
lab = DataFrame(number = [e.number for e in tree_subs.edge],
                label = [string(round(e.length, digits=3))
                         for e in tree_subs.edge])
# PhyloPlots puts the root at x=1
plot(tree_subs, useedgelength=true, edgelabel=lab, tipcex=1.6, edgecex=1.4,
     xlim=[1, 1.9])
R"title"("substitutions/site after the multipliers")
R"dev.off"()
println("figure: ", joinpath(datadir, "speciestree_angiosperm.png"))
