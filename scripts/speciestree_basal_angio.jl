#=
code used to build a second "true" species tree, from basal angiosperms,
used by SimPhy to simulate gene trees. Same procedure as speciestree.jl
(reptile tree, Crawford et al. 2012 UCE data), applied to nuclear genes from
the One Thousand Plant Transcriptomes Initiative (1KP, 2019, Nature 574:679,
doi:10.1038/s41586-019-1693-2).
Branch lengths of the species tree are in coalescent units;
per-branch multipliers give the lineage-specific substitution rates.

Two taxon sets are built, differing only in the gymnosperm outgroup:
  "ginkgo": A = Ginkgo biloba (SGTW)      -> 106 loci with all 8 taxa
  "pinus" : A = Pinus ponderosa (JBND)    -> 132 loci with all 8 taxa
  B = Amborella trichopoda (URDJ)         Amborellales
  C = Nymphaea sp. (PZRT)                 Nymphaeales
  D = Nuphar advena (WTKZ)                Nymphaeales
  E = Illicium floridanum (VZCI)          Austrobaileyales
  F = Ascarina rubricaulis (WZFE)         Chloranthales
  G = Michelia maudiae (XQWC)             magnoliids (Magnoliaceae)
  H = Ceratophyllum demersum (NPND)       Ceratophyllales
Codes are 1KP 4-letter sample codes (data/onekp/annotations.csv).

- data: 1KP capstone single-copy gene alignments, nucleotide version with all
  three codon positions ("FNA2AA", 374 genes, ~850 taxa), downloaded from
  the CyVerse Data Commons (anonymous):
  base=https://data.cyverse.org/dav-anon/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/alignments_and_trees
  curl -o data/onekp/alignments-FNA2AA.tar.bz $base/alignments/alignments-FNA2AA.tar.bz   # 77 MB
  curl -o data/onekp/annotations.csv          $base/misc/annotations.csv                  # code -> species
  mkdir -p data/onekp/FNA2AA && tar xjf data/onekp/alignments-FNA2AA.tar.bz -C data/onekp/FNA2AA
  (this is the version *before* 1KP's fragmentary-sequence filter; see step 4)
- subset to the 8 taxa, keeping genes present in all 8 and dropping codons
  that are gaps in all 8: (from the repo root)
  python3 scripts/angiosperm/subset_onekp.py --name ginkgo --codes SGTW URDJ PZRT WTKZ VZCI WZFE XQWC NPND
  python3 scripts/angiosperm/subset_onekp.py --name pinus  --codes JBND URDJ PZRT WTKZ VZCI WZFE XQWC NPND
  which creates data/onekp/subset_<name>/alignments/<gene>.fasta and taxa.csv
- gene trees were estimated with IQ-TREE (v3.1.3 here; the reptile trees used
  v2.4.0 with the same options), one tree per gene, 1000 ultra-fast bootstraps,
  HKY+F+G per gene so that kappa, base frequencies and the gamma shape can be
  re-fitted as in notes/choice-seqgen-parameters.md:
  cd data/onekp/subset_<name>
  iqtree -S alignments --prefix iqtree/loci -T 4 -B 1000 -wbtl -m MFP -mset HKY -mrate G -mfreq F -seed 1
  which creates in `iqtree/`
  * a single file `loci.treefile` with the ML tree, one per locus,
  * 1 file `loci.ufboot` with all bootstrap trees (1000 per locus),
  * `loci.best_model.nex` with the per-gene HKY+F+G parameters.
  followed by (same as separate-boot-bygene.jl in the reptiles repo):
  python3 ../../../scripts/angiosperm/split_bootstraps.py iqtree 1000
  to get 1 bootstrap file per gene and the list `BSlistfiles`.
- loci with any branch >= 5 substitutions/site were removed (step 4 below
  writes `iqtree/loci.filtered.treefile`, `BSlistfiles.filtered`,
  `dropped_loci.txt`): 2 loci in the ginkgo set, 1 in the pinus set,
  all caused by a fragmentary sequence (e.g. Ginkgo with 642 non-gap sites
  in locus 6405, giving IQ-TREE's maximum branch length of 10).
- species tree estimated with ASTRAL (v5.7.8) with IQ-TREE input, using
  astral -i iqtree/loci.filtered.treefile -b BSlistfiles.filtered -r 1000 -o astral/species.tre
  The output file has 1002 lines: 1000 bootstrap trees, then the ASTRAL tree
  with bootstrap support, then the same tree + edge lengths in coalescent units.
- edge lengths from ASTRAL (in coalescent units) were rounded to 2 digits
- missing edge lengths were assigned in such a way that the tree is ultrametric.
- substitution rates: as in the paper's Methods ("Substitution Rate Variation"),
  i.e. as ratevariation.jl in the reptiles repo:
  * locus height d_l = median genetic distance between ingroup taxa and the
    outgroup, in the IQ-TREE gene tree of locus l
  * each gene tree rescaled by mean(d_l)/d_l, pairwise distances averaged
    across loci, then fitted to the species tree topology by ordinary least
    squares with calibratefrompairwisedistances! (non-ultrametric) -> d_i
  * rate variation across genes: log-normal fitted to d_l / mean(d_l)
  * bar_r = sum(d_i) / sum(tau_i), SimPhy -su = bar_r / 2Ne, m_i = d_i / (tau_i bar_r)

Run from the repo root, for one taxon set at a time:
  julia --project=. scripts/speciestree_basal_angio.jl ginkgo
  julia --project=. scripts/speciestree_basal_angio.jl pinus
Expected outputs are copy-pasted in comments below, after each step.
Requires PhyloNetworks >= 1.0 (tested with 1.3.1) and QuartetNetworkGoodnessFit 1.0.
=#

using PhyloNetworks
using QuartetNetworkGoodnessFit # has a function to ultrametrize a network
using CSV, DataFrames
using Distributions, Statistics, Printf

setname = length(ARGS) >= 1 ? ARGS[1] : "ginkgo"   # "ginkgo" or "pinus"
datadir = joinpath("data", "onekp", "subset_" * setname)  # run from repo root
eff_pop = 1000 # = 2Ne, diploid effective population size, as for the reptile tree

taxa = CSV.read(joinpath(datadir, "taxa.csv"), DataFrame) # letter,code,species
codes = String.(taxa.code)
code2letter = Dict(String(r.code) => String(r.letter) for r in eachrow(taxa))
outgroup = codes[1] # SGTW (ginkgo) or JBND (pinus): letter A

# set of tip names below an edge: used to match edges across trees
function leafset(e::PhyloNetworks.Edge)
    ch = getchild(e)
    ch.leaf && return Set([ch.name])
    s = Set{String}()
    for e2 in ch.edge
        getparent(e2) === ch && union!(s, leafset(e2))
    end
    return s
end
totallength(t) = sum(e.length for e in t.edge)

#--------------------------------------------------------------------#
# 1. ASTRAL tree in coalescent units
#--------------------------------------------------------------------#
astral_lines = filter(!isempty, readlines(joinpath(datadir, "astral", "species.tre")))
speciestree_string = astral_lines[end] # last line: ASTRAL tree with CU lengths
#= ginkgo set:
"(SGTW,(((PZRT,WTKZ)100.0:2.4140805972193884,(VZCI,(NPND,(WZFE,XQWC)100.0:0.3658915990436563)100.0:0.9900459081166504)100.0:0.7329185039868412)91.60000000000001:0.18530766214824268,URDJ):0.0);"
pinus set:
"(JBND,(((PZRT,WTKZ)100.0:2.4774716713004823,(VZCI,(NPND,(WZFE,XQWC)100.0:0.24961012450237297)100.0:1.0543633370578753)100.0:0.711980128893038)99.5:0.2510479392720005,URDJ):0.0);"
Topology, for both outgroups:
 (A,(B,((C,D),(E,(H,(F,G)))))) = (outgroup,(Amborella,((Nymphaea,Nuphar),
   (Illicium,(Ceratophyllum,(Ascarina,Michelia))))))
which matches Zuntini et al. 2024 (Nature) once monocots and eudicots are
dropped. All bootstrap supports are 100 except Amborella-alone as sister to
the other angiosperms: 91.6 (ginkgo) and 99.5 (pinus).
=#
tree = readnewick(speciestree_string)
# remove bootstrap values from node names (SimPhy doesn't accept them)
for n in tree.node
    isleaf(n) && continue
    n.name = ""
end
rootatnode!(tree, outgroup)
# round edge lengths to avoid 15 digits: not significant for simulations
for e in tree.edge
    e.length == -1.0 && continue # missing
    e.length = round(e.length, digits=2)
end
println("ASTRAL tree, coalescent units, rounded:\n", writenewick(tree))
#= the tree lacks external edge lengths, and a length for the ingroup stem
   edge (ASTRAL writes 0.0 for it), as for the reptile tree.
ginkgo: (SGTW,(((PZRT,WTKZ):2.41,(VZCI,(NPND,(WZFE,XQWC):0.37):0.99):0.73):0.19,URDJ):0.0);
pinus:  (JBND,(((PZRT,WTKZ):2.48,(VZCI,(NPND,(WZFE,XQWC):0.25):1.05):0.71):0.25,URDJ):0.0);
=#

#--------------------------------------------------------------------#
# 2. assign missing lengths, ultrametrize
#--------------------------------------------------------------------#
# assign 0.5 coalescent units to the ingroup stem edge (as for reptiles)
rootnode = getroot(tree)
for e in rootnode.edge
    getchild(e).leaf && continue # outgroup edge: assigned by ultrametrize! below
    e.length = 0.5
end
# assign 0.5 coalescent unit to one external edge going to the deepest cherry
# (as was done for Anolis, in the deepest cherry of the reptile tree):
# this sets the height of the tree, and all other external edges are then
# assigned by ultrametrize! to be >= 0.5.
function depth(n) # from the root, in coalescent units
    d = 0.0
    while n !== rootnode
        e = getparentedge(n); d += e.length; n = getparent(n)
    end
    return d
end
internalnodes = filter(n -> !n.leaf && n !== rootnode, tree.node)
deepest = internalnodes[argmax(depth.(internalnodes))]
for e in deepest.edge
    if getparent(e) === deepest && getchild(e).leaf
        e.length = 0.5 # Nymphaea (PZRT) in both sets
        println("assigned 0.5 coalescent unit to the external edge to ", getchild(e).name)
        break
    end
end
QuartetNetworkGoodnessFit.ultrametrize!(tree, true) # verbose=true
for e in tree.edge e.length = round(e.length, digits=2); end
println("ultrametric tree, coalescent units:\n", writenewick(tree))
#= we get this below:
ginkgo: (SGTW:3.6,(((PZRT:0.5,WTKZ:0.5):2.41,(VZCI:2.18,(NPND:1.19,(WZFE:0.82,XQWC:0.82):0.37):0.99):0.73):0.19,URDJ:3.1):0.5);
pinus:  (JBND:3.73,(((PZRT:0.5,WTKZ:0.5):2.48,(VZCI:2.27,(NPND:1.22,(WZFE:0.97,XQWC:0.97):0.25):1.05):0.71):0.25,URDJ:3.23):0.5);
total length: 17.90 (ginkgo), 18.63 (pinus) CU; reptile tree: 17.48 CU.
height: 3.60 (ginkgo), 3.73 (pinus); reptile tree: 3.44.
=#

# tips renamed A-H in the order of taxa.csv (A = outgroup)
for n in tree.node
    n.leaf && (n.name = code2letter[n.name])
end
tree_cu = deepcopy(tree)
# tree with edge lengths in number of generations: gen = cu * 2Ne
tree_ngen = deepcopy(tree_cu)
for e in tree_ngen.edge
    e.length *= eff_pop
end
speciestree_ngen = replace(writenewick(tree_ngen), ".0" => "")
println("tree in generations:\n", speciestree_ngen)
#= this below is the SimPhy species tree without lineage rate variation:
ginkgo: (A:3600,(((C:500,D:500):2410,(E:2180,(H:1190,(F:820,G:820):370):990):730):190,B:3100):500);
pinus:  (A:3730,(((C:500,D:500):2480,(E:2270,(H:1220,(F:970,G:970):250):1050):710):250,B:3230):500);
=#

#--------------------------------------------------------------------#
# 3. gene trees: locus heights d_l, filtering of saturated loci
#--------------------------------------------------------------------#
# names of loci, in the order used by IQ-TREE (= order in loci.treefile)
locusnames = String[]
for l in readlines(joinpath(datadir, "iqtree", "loci.best_scheme.nex"))
    m = match(r"^\s*charset\s+(\S+)\s*=", l)
    isnothing(m) || push!(locusnames, m.captures[1])
end
genetrees_all = readmultinewick(joinpath(datadir, "iqtree", "loci.treefile"))
length(genetrees_all) == length(locusnames) || error("number of gene trees != number of loci")
# The 1KP nucleotide alignments are the version before 1KP's filter of
# fragmentary sequences. A sequence covering a small fraction of a locus can
# get a huge branch length (up to IQ-TREE's maximum of 10 subs/site).
# We remove any locus with a branch >= 5 substitutions/site.
maxbranch(gt) = maximum(e.length for e in gt.edge)
keep = [maxbranch(gt) < 5.0 for gt in genetrees_all]
println("loci dropped for a branch >= 5 subs/site: ", locusnames[.!keep])
#= ginkgo: ["6405.fasta", "6544.fasta"]  (Ginkgo at 10.0; Ceratophyllum at 5.7)
   pinus:  ["6544.fasta"]                (Ceratophyllum at 6.1)
=#
genetrees = genetrees_all[keep]
nloci = length(genetrees)
# files used by ASTRAL (see header): filtered gene trees & bootstrap file list
open(joinpath(datadir, "iqtree", "loci.filtered.treefile"), "w") do io
    for gt in genetrees; println(io, writenewick(gt)); end
end
bslist = filter(!isempty, readlines(joinpath(datadir, "BSlistfiles")))
open(joinpath(datadir, "BSlistfiles.filtered"), "w") do io
    for f in bslist[keep]; println(io, f); end
end
open(joinpath(datadir, "dropped_loci.txt"), "w") do io
    for l in locusnames[.!keep]; println(io, l); end
end

# locus height d_l: median genetic distance between the ingroup taxa and the
# outgroup, in the gene tree of locus l. The median is robust to lineages
# with extreme rates (e.g. Ceratophyllum).
function distancematrix(gt) # pairwise distances, rows/columns in the order of `codes`
    D = pairwisetaxondistancematrix(gt)
    tl = tiplabels(gt)
    idx = [findfirst(==(c), tl) for c in codes]
    return D[idx, idx]
end
d_l = [median(distancematrix(gt)[1, 2:end]) for gt in genetrees]
mean_d = mean(d_l)
@printf("%d loci; locus heights d_l: mean = %.6f, median = %.4f, min = %.4f, max = %.4f\n",
        nloci, mean_d, median(d_l), minimum(d_l), maximum(d_l))
#= ginkgo: 106 loci; mean = 0.821788, median 0.7345, min 0.4276, max 2.3401
   pinus:  132 loci; mean = 0.855225, median 0.7948, min 0.4755, max 2.0070
   reptile UCEs (1,145 loci): mean = 0.096825
   These are ~8.5x larger than for the reptile UCEs, although both trees span
   ~310 My from the root: coding genes with all 3 codon positions (3rd
   positions saturated: Ginkgo-ingroup p-distance = 0.54 at position 3,
   vs 0.24 and 0.16 at positions 1 and 2) vs ultra-conserved elements.
=#

#--------------------------------------------------------------------#
# 4. rate variation across genes: distribution of d_l
#--------------------------------------------------------------------#
# relative rate of locus l: d_l / mean(d_l), with mean 1 across loci
rel = d_l ./ mean_d
for family in [Normal, LogNormal, Gamma]
    d = fit(family, rel)
    loglik = sum(logpdf.(d, rel))
    println("loglik=$(round(loglik,digits=1)), $d")
end
#= ginkgo:
loglik=-49.5, Normal{Float64}(μ=1.0, σ=0.38607582927241535)
loglik=-26.7, LogNormal{Float64}(μ=-0.059443224366056056, σ=0.3319484557776711)
loglik=-32.3, Gamma{Float64}(α=8.574659775942234, θ=0.11662270295617816)
pinus:
loglik=-40.5, Normal{Float64}(μ=1.0, σ=0.32897020481715566)
loglik=-21.2, LogNormal{Float64}(μ=-0.04679504467967938, σ=0.29883314963093305)
loglik=-25.8, Gamma{Float64}(α=10.848900211058975, θ=0.092175241779866)
The log-normal fits best, as for the reptile UCEs (sigma = 0.6164 there).
=#
sigma_gene = fit(LogNormal, rel).σ
# log-normal with mean 1: mu = -sigma^2/2, so that E(r_g mu) = mu
mu_gene = -sigma_gene^2 / 2
@printf("SimPhy option for rate variation across genes: -hl ln:%.6f,%.6f\n", mu_gene, sigma_gene)
#= ginkgo: -hl ln:-0.055095,0.331948
   pinus:  -hl ln:-0.044651,0.298833
   reptile: -hl ln:-0.19,0.6164414002968976
=#

#--------------------------------------------------------------------#
# 5. tree with substitutions/site per branch, d_i
#--------------------------------------------------------------------#
# rescale each gene tree by mean(d_l)/d_l to eliminate rate variation across
# loci, then average the pairwise genetic distances across gene trees
D_mean = zeros(length(codes), length(codes))
for (gt, h) in zip(genetrees, d_l)
    global D_mean .+= distancematrix(gt) .* (mean_d / h)
end
D_mean ./= nloci
# fit these distances to the species tree topology with ordinary least squares.
# ultrametric=false (the default is true): a clock would erase the rate
# variation across lineages that we want to capture.
sub_tree = readnewick(speciestree_string)
for n in sub_tree.node
    isleaf(n) || (n.name = "")
end
rootatnode!(sub_tree, outgroup)
for e in sub_tree.edge e.length = 0.1; end # starting values
calibratefrompairwisedistances!(sub_tree, D_mean, codes; verbose=false, ultrametric=false)
for n in sub_tree.node
    n.leaf && (n.name = code2letter[n.name])
end
sub_tree_string = writenewick(sub_tree)
println("tree with substitutions/site per branch:\n", sub_tree_string)
#= ginkgo:
"(A:0.2541068363072812,(((C:0.16254254609231938,D:0.09178843588136217):0.20987812653084184,(E:0.20948517154173973,(H:0.38581531682726133,(F:0.15138465362462225,G:0.16188365438146046):0.02729876021912248):0.044884418198549395):0.03804069287399796):0.014503657212708115,B:0.3052436921289787):0.2541068363072812);"
pinus:
"(A:0.27557928636789697,(((C:0.1598895856743684,D:0.09098697456226051):0.19860449693000193,(E:0.2024878381104711,(H:0.38226689435509986,(F:0.14702765692027303,G:0.15767462145320987):0.020844552284902448):0.048377110621895425):0.03267716342537563):0.01780154505692354,B:0.29281420563980926):0.27557928636789697);"
Note: the outgroup edge and the ingroup stem edge get the same length:
their sum is identifiable from pairwise distances, not the split. The
reptile sub_tree has the same property (Homo: 0.03473, root edge: 0.03473).
=#
# goodness of fit of the tree distances to the mean pairwise distances
D_fit = pairwisetaxondistancematrix(sub_tree)
tl = tiplabels(sub_tree); idx = [findfirst(==(code2letter[c]), tl) for c in codes]
@printf("max |D_mean - D_fit| = %.5f\n", maximum(abs.(D_mean .- D_fit[idx, idx])))
# ginkgo: 0.00959 ; pinus: 0.00626

#--------------------------------------------------------------------#
# 6. genome-wide rate, and per-branch multipliers m_i
#--------------------------------------------------------------------#
cu_total = totallength(tree_cu)     # 17.90 (ginkgo), 18.63 (pinus)
sub_total = totallength(sub_tree)   # 2.3110 (ginkgo), 2.3026 (pinus)
sub_per_cu = sub_total / cu_total   # bar_r: 0.129104 (ginkgo), 0.123597 (pinus); reptile: 0.019526
@printf("bar_r = %.6f substitutions per site per coalescent unit\n", sub_per_cu)
# SimPhy baseline substitution rate, per site per generation: -su bar_r / 2Ne
@printf("SimPhy -su f:%.10e\n", sub_per_cu / eff_pop)
# ginkgo: 1.2910413e-4 ; pinus: 1.2359702e-4 ; reptile: 1.9526e-5

#= per-branch multiplier: m_i = d_i / (bar_r * tau_i), dimensionless, so that
   SimPhy's expected substitutions along branch i are
   t_gen,i * mu * m_i = (tau_i 2Ne) (bar_r / 2Ne) (d_i / (bar_r tau_i)) = d_i
   (the "updated" multiplier of speciestree.jl, makeSimPhytree_dividied_by_cu)
=#
sub_by_clade = Dict(leafset(e) => e.length for e in sub_tree.edge)
mult_by_clade = Dict{Set{String},Float64}()
println("branch\tCU\td_i\tm_i")
for e in tree_ngen.edge
    ls = leafset(e)
    tau = e.length / eff_pop
    d = sub_by_clade[ls]
    mult_by_clade[ls] = d / (sub_per_cu * tau)
    @printf("%s\t%.2f\t%.4f\t%.3f\n", join(sort(collect(ls))), tau, d, mult_by_clade[ls])
end
# SimPhy tree edges annotated as `generations*multiplier`, written directly
# from the tree (edges matched by the set of tips below them, not by their
# length: two internal edges can have the same length, e.g. pinus set).
function simphynewick(node, parentedge)
    s = node.leaf ? node.name :
        "(" * join([simphynewick(getchild(e), e) for e in node.edge if getparent(e) === node], ",") * ")"
    isnothing(parentedge) && return s * ";"
    return s * ":" * string(Int(round(parentedge.length))) * "*" *
           string(round(mult_by_clade[leafset(parentedge)], digits=6))
end
simphytree = simphynewick(getroot(tree_ngen), nothing)
println("SimPhy species tree with lineage-specific multipliers:\n", simphytree)
#= ginkgo:
branch  CU     d_i     m_i
A       3.60   0.2541  0.547   Ginkgo
C       0.50   0.1625  2.518   Nymphaea
D       0.50   0.0918  1.422   Nuphar
CD      2.41   0.2099  0.675
E       2.18   0.2095  0.744   Illicium
H       1.19   0.3858  2.511   Ceratophyllum
F       0.82   0.1514  1.430   Ascarina
G       0.82   0.1619  1.529   Michelia
FG      0.37   0.0273  0.571
FGH     0.99   0.0449  0.351
EFGH    0.73   0.0380  0.404
CDEFGH  0.19   0.0145  0.591
B       3.10   0.3052  0.763   Amborella
BCDEFGH 0.50   0.2541  3.936   ingroup stem
"(A:3600*0.546731,(((C:500*2.518008,D:500*1.421929):2410*0.674544,(E:2180*0.744315,(H:1190*2.511265,(F:820*1.429974,G:820*1.529147):370*0.57148):990*0.351172):730*0.403632):190*0.591267,B:3100*0.762685):500*3.936465);"
pinus:
"(A:3730*0.597764,(((C:500*2.587274,D:500*1.472317):2480*0.647932,(E:2270*0.721714,(H:1220*2.535123,(F:970*1.226364,G:970*1.315171):250*0.674598):1050*0.372772):710*0.372373):250*0.576116,B:3230*0.733469):500*4.459322);"
Fastest lineages: Ceratophyllum (2.5x) and Nymphaea (2.5x), plus the ingroup
stem (3.9-4.5x). Multiplier range 0.35-4.5, vs 0.19-9.2 for the reptile tree
(snake 9.2x, lizard 4.1x).
=#

# save everything needed by simulation.jl for this taxon set
open(joinpath(datadir, "simphy_tree_basal_angio.txt"), "w") do io
    println(io, "# taxon set: ", setname, "; ",
            join(["$(r.letter)=$(r.species) ($(r.code))" for r in eachrow(taxa)], ", "))
    println(io, "# loci used: ", nloci)
    @printf(io, "# -su f:%.10e\n", sub_per_cu / eff_pop)
    @printf(io, "# -hl ln:%.6f,%.6f\n", mu_gene, sigma_gene)
    println(io, "# species tree, coalescent units:\n", writenewick(tree_cu))
    println(io, "# species tree, substitutions/site:\n", sub_tree_string)
    println(io, "# SimPhy species tree, no lineage rate variation (generations, 2Ne=$eff_pop):\n", speciestree_ngen)
    println(io, "# SimPhy species tree with lineage-specific multipliers:\n", simphytree)
end
