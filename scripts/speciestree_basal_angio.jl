#=
code used to build a second "true" species tree, from basal angiosperms,
used by SimPhy to simulate gene trees. Same procedure as speciestree.jl
(reptile tree, Crawford et al. 2012 UCE data), applied to the Angiosperms353
nuclear genes of the Kew Tree of Life Explorer, release 4.0 (April 2026),
the data behind Zuntini et al. 2024, Nature 629:843,
doi:10.1038/s41586-024-07324-0.
Branch lengths of the species tree are in coalescent units;
per-branch multipliers give the lineage-specific substitution rates.

Two taxon sets are built, differing only in the gymnosperm outgroup.
Letters, species, and Kew "Sequence_ID" of the sample used:
  A = outgroup, one of
      "pinus" : Pinus ponderosa        ERR2040837       (1KP transcriptome)
      "ginkgo": Ginkgo biloba          SGTW             (1KP transcriptome)
  B = Amborella trichopoda             GCF_000471905.2  (genome) Amborellales
  C = Nymphaea nouchali                GCF_008831285.2  (genome) Nymphaeales
  D = Brasenia schreberi               GCA_030142015.1  (genome) Nymphaeales
  E = Schisandra chinensis             SRR26404306      (reads) Austrobaileyales
  F = Chloranthus sessilifolius        GCA_021018995.1  (genome) Chloranthales
  G = Liriodendron chinense            GCA_003013855.2  (genome) magnoliids
  H = Ceratophyllum demersum           SRR6374676       (reads) Ceratophyllales
Design mirrors the reptile tree: an outgroup (Homo), two cherries (C,D and
F,G; crocodilians and birds), one lone deep taxon (E; turtle) and one fast,
historically unstable taxon (H; lizard + snake).

- data: official per-gene DNA alignments of the Kew release (353 genes,
  ~20,000 samples per gene, 20-27 MB per file), anonymous access at
  https://sftp.kew.org/pub/paftol/current_release/fasta/alignments/
  Each file is downloaded, the 8 samples are pulled out (by Sequence_ID),
  columns that are gaps in all 8 are removed, and the download is deleted:
  python3 scripts/angiosperm/subset_a353_official.py --name pinus \
    --samples PINPON=ERR2040837 AMBTRI=GCF_000471905.2 \
    NYMNOU=GCF_008831285.2 BRASCH=GCA_030142015.1 SCHCHI=SRR26404306 \
    CHLSES=GCA_021018995.1 LIRCHI=GCA_003013855.2 CERDEM=SRR6374676
  (same with --name ginkgo and GINBIL=SGTW as first sample)
  which creates data/a353/<name>/alignments/<gene>.fasta, taxa.csv and
  gene_occupancy.csv. All 353 genes are kept (genes need >= 4 taxa);
  genes recovered per taxon in the pinus set: Pinus 322, Amborella 347,
  Nymphaea 349, Brasenia 346, Schisandra 353, Chloranthus 352,
  Liriodendron 351, Ceratophyllum 289; 260 genes have all 8 taxa.
- gene trees were estimated with IQ-TREE (v3.1.3 here; v2.4.0 for the
  reptile trees, same options), one run per gene, 1000 ultra-fast
  bootstraps, HKY+F+G so that kappa, base frequencies and the gamma
  shape can be re-fitted as in notes/choice-seqgen-parameters.md:
  cd data/a353/<name>
  for f in alignments/*.fasta; do g=$(basename $f .fasta)
    iqtree -s $f -st DNA -B 1000 -wbtl -m MFP -mset HKY -mrate G \
           -mfreq F -seed 1 -T 1 --prefix iqtree_pergene/$g
  done
  (-st DNA: Kew pads low-coverage regions with N, which can defeat
  IQ-TREE's automatic detection of the sequence type.
  One run per gene rather than one partitioned run `iqtree -S`: with -S,
  the bootstrap trees of genes lacking some taxa were written with wrong
  tip names by IQ-TREE 3.1.3, which corrupted ASTRAL's bootstrap support;
  the reptile UCE loci all had all taxa, so -S was fine there.)
  followed by (replaces separate-boot-bygene.jl in the reptiles repo):
  python3 ../../../scripts/angiosperm/collect_pergene.py
  which writes `loci.txt` (gene names), `iqtree/loci.treefile` (ML trees,
  same order), `iqtree/bootstrap/<gene>.ufboot`, the list `BSlistfiles`,
  and `seqgen_params.csv` (per-gene kappa, gamma shape, base frequencies).
- species tree estimated with ASTRAL (v5.7.8) with IQ-TREE input, using
  astral -i iqtree/loci.treefile -b BSlistfiles -r 1000 -o astral/species.tre
  The output file has 1002 lines: 1000 bootstrap trees, then the ASTRAL
  tree with bootstrap support, then the same tree + edge lengths in
  coalescent units.
- step 3 below checks for loci with a branch >= 5 substitutions/site
  (a symptom of a fragmentary sequence: IQ-TREE's maximum is 10) and
  writes `iqtree/loci.filtered.treefile` and `BSlistfiles.filtered`.
  No locus was removed from the Kew data. If some were, ASTRAL would need
  to be re-run on the filtered files (from data/a353/<name>) and this
  script run again.
- edge lengths from ASTRAL (in coalescent units) were rounded to 2 digits
- missing edge lengths were assigned in such a way that the tree is
  ultrametric.
- substitution rates: as in the paper's Methods ("Substitution Rate
  Variation"), i.e. as ratevariation.jl in the reptiles repo:
  * locus height d_l = median genetic distance between ingroup taxa and
    the outgroup, in the IQ-TREE gene tree of locus l
  * each gene tree rescaled by mean(d_l)/d_l, pairwise distances averaged
    across loci, then fitted to the species tree topology by ordinary
    least squares with calibratefrompairwisedistances! (non-ultrametric)
    -> d_i
  * rate variation across genes: log-normal fitted to d_l / mean(d_l)
  * bar_r = sum(d_i) / sum(tau_i); SimPhy -su = bar_r / 2Ne;
    m_i = d_i / (tau_i bar_r)

Run from the repo root, for one taxon set at a time:
  julia --project=. scripts/speciestree_basal_angio.jl pinus
  julia --project=. scripts/speciestree_basal_angio.jl ginkgo
Expected outputs (pinus set) are copy-pasted in comments after each step.
Requires PhyloNetworks >= 1.0 (tested with 1.3.1) and
QuartetNetworkGoodnessFit 1.0.
=#

using PhyloNetworks
using QuartetNetworkGoodnessFit # has a function to ultrametrize a network
using CSV, DataFrames
using Distributions, Statistics, Printf

setname = length(ARGS) >= 1 ? ARGS[1] : "pinus"   # "pinus" or "ginkgo"
datadir = joinpath("data", "a353", setname)       # run from repo root
eff_pop = 1000 # = 2Ne, diploid effective population size, as for reptiles

taxa = CSV.read(joinpath(datadir, "taxa.csv"), DataFrame) # letter,code,species
codes = String.(taxa.code)
code2letter = Dict(String(r.code) => String(r.letter) for r in eachrow(taxa))
outgroup = codes[1] # PINPON or GINBIL: letter A
ntax = length(codes)

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
astral_lines = filter(!isempty,
                      readlines(joinpath(datadir, "astral", "species.tre")))
speciestree_string = astral_lines[end] # last: ASTRAL tree with CU lengths
#= pinus set:
"(PINPON,(((NYMNOU,BRASCH)100.0:3.9281879422815016,(SCHCHI,(CERDEM,
(CHLSES,LIRCHI)100.0:0.32843173162533745)100.0:0.8933869344080152)100.0:
0.732254263805713)99.3:0.20062688289231798,AMBTRI):0.0);"
Topology, for both outgroups:
 (A,(B,((C,D),(E,(H,(F,G)))))) = (outgroup,(Amborella,((Nymphaea,Brasenia),
   (Schisandra,(Ceratophyllum,(Chloranthus,Liriodendron))))))
which matches Zuntini et al. 2024 once monocots and eudicots are dropped.
All bootstrap supports are 100 except Amborella alone as sister to the
other angiosperms: 95.7 (pinus); 100 (ginkgo).
ginkgo set: (AMBTRI,(GINBIL,((SCHCHI,(CERDEM,(CHLSES,LIRCHI)100.0:0.34)
100.0:0.93)100.0:0.74,(NYMNOU,BRASCH)100.0:4.61)100.0:0.25):0.0);
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
pinus: (PINPON,(((NYMNOU,BRASCH):3.93,(SCHCHI,(CERDEM,(CHLSES,LIRCHI):0.33)
       :0.89):0.73):0.2,AMBTRI):0.0);
=#

#--------------------------------------------------------------------#
# 2. assign missing lengths, ultrametrize
#--------------------------------------------------------------------#
# assign 0.5 coalescent units to the ingroup stem edge (as for reptiles)
rootnode = getroot(tree)
for e in rootnode.edge
    getchild(e).leaf && continue # outgroup edge: set by ultrametrize! below
    e.length = 0.5
end
# assign 0.5 coalescent unit to one external edge going to the deepest
# cherry (as was done for Anolis, in the deepest cherry of the reptile
# tree): this sets the height of the tree, and all other external edges
# are then assigned by ultrametrize! to be >= 0.5.
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
        e.length = 0.5 # Nymphaea (NYMNOU)
        println("assigned 0.5 coalescent unit to the external edge to ",
                getchild(e).name)
        break
    end
end
QuartetNetworkGoodnessFit.ultrametrize!(tree, true) # verbose=true
for e in tree.edge e.length = round(e.length, digits=2); end
println("ultrametric tree, coalescent units:\n", writenewick(tree))
#= we get this below:
pinus: (PINPON:5.13,(((NYMNOU:0.5,BRASCH:0.5):3.93,(SCHCHI:3.7,(CERDEM:2.81,
       (CHLSES:2.48,LIRCHI:2.48):0.33):0.89):0.73):0.2,AMBTRI:4.63):0.5);
ginkgo: (GINBIL:5.86,(((SCHCHI:4.37,(CERDEM:3.44,(CHLSES:3.1,LIRCHI:3.1)
       :0.34):0.93):0.74,(BRASCH:0.5,NYMNOU:0.5):4.61):0.25,AMBTRI:5.36):0.5);
total length: 28.86 CU (pinus), 33.02 (ginkgo); reptile tree: 17.48 CU.
height: 5.13 (pinus), 5.86 (ginkgo); reptile tree: 3.44.
The angiosperm tree is taller: with genome-quality sequences, gene trees
are less discordant, which ASTRAL translates into longer branches in CU.
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
pinus: (A:5130,(((C:500,D:500):3930,(E:3700,(H:2810,(F:2480,G:2480):330)
       :890):730):200,B:4630):500);
=#

#--------------------------------------------------------------------#
# 3. gene trees: check for saturated loci, locus heights d_l
#--------------------------------------------------------------------#
# names of loci, in the order of loci.treefile (see collect_pergene.py)
locusnames = filter(!isempty, readlines(joinpath(datadir, "loci.txt")))
genetrees_all = readmultinewick(joinpath(datadir, "iqtree", "loci.treefile"))
length(genetrees_all) == length(locusnames) ||
    error("number of gene trees != number of loci")
# A sequence covering a small fraction of a locus can get a huge branch
# length (up to IQ-TREE's maximum of 10 subs/site). Legitimate root-to-tip
# distances are < 1 here, so any locus with a branch >= 5 substitutions/site
# is flagged and removed. None in the Kew data (which was filtered for
# coverage by Kew); this was needed for the 1KP alignments tried earlier.
maxbranch(gt) = maximum(e.length for e in gt.edge)
keep = [maxbranch(gt) < 5.0 for gt in genetrees_all]
println("loci dropped for a branch >= 5 subs/site: ", locusnames[.!keep])
# pinus: String[]  (nothing dropped)
genetrees = genetrees_all[keep]
nloci = length(genetrees)
# files for ASTRAL if loci were dropped (see header)
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

# locus height d_l: median genetic distance between the ingroup taxa and
# the outgroup, in the gene tree of locus l. The median is robust to
# lineages with extreme rates (e.g. Ceratophyllum).
# Gene trees may lack some taxa (loci are kept if >= 4 of the 8 taxa are
# present): missing entries are NaN in the distance matrix below.
function distancematrix(gt) # rows/columns in the order of `codes`
    D = pairwisetaxondistancematrix(gt)
    tl = tiplabels(gt)
    idx = [findfirst(==(c), tl) for c in codes]
    M = fill(NaN, ntax, ntax)
    for (a, ia) in enumerate(idx), (b, ib) in enumerate(idx)
        (isnothing(ia) || isnothing(ib)) && continue
        M[a, b] = D[ia, ib]
    end
    return M
end
# the locus height needs the outgroup: loci without it are used for the
# species tree (ASTRAL) but not for the substitution-rate calibration
Dmats = [distancematrix(gt) for gt in genetrees]
hasoutgroup = [!all(isnan, D[1, 2:end]) for D in Dmats]
println("loci without the outgroup, excluded from rate calibration: ",
        count(!, hasoutgroup))
d_l = [median(filter(!isnan, D[1, 2:end])) for D in Dmats[hasoutgroup]]
mean_d = mean(d_l)
@printf("%d loci with outgroup; locus heights d_l: mean = %.6f, ",
        length(d_l), mean_d)
@printf("median = %.4f, min = %.4f, max = %.4f\n",
        median(d_l), minimum(d_l), maximum(d_l))
#= pinus:  31 loci without the outgroup; 322 loci with outgroup;
   d_l: mean = 0.713752, median 0.6629, min 0.3663, max 3.1012
   ginkgo: 160 loci without the outgroup (Ginkgo is a fragmentary 1KP
   transcriptome); 193 loci with outgroup; d_l: mean = 0.664268
   reptile UCEs (1,145 loci): mean = 0.096825
   d_l is ~7x larger than for the reptile UCEs, although both trees span
   ~310 My from the root: protein-coding genes with all 3 codon positions
   (3rd positions saturated: outgroup-ingroup p-distance ~0.5 at
   position 3, vs ~0.2 and ~0.15 at positions 1 and 2) vs ultra-conserved
   elements.
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
#= pinus:
loglik=-74.6, Normal{Float64}(μ=1.0, σ=0.3050)
loglik=0.9, LogNormal{Float64}(μ=-0.0339, σ=0.2500)
loglik=-14.7, Gamma{Float64}(α=14.89, θ=0.0671)
The log-normal fits best, as for the reptile UCEs (sigma = 0.6164 there).
=#
sigma_gene = fit(LogNormal, rel).σ
# log-normal with mean 1: mu = -sigma^2/2, so that E(r_g mu) = mu
mu_gene = -sigma_gene^2 / 2
@printf("SimPhy option for rate variation across genes: -hl ln:%.6f,%.6f\n",
        mu_gene, sigma_gene)
#= pinus:   -hl ln:-0.031126,0.249505
   ginkgo:  -hl ln:-0.034203,0.261546
   reptile: -hl ln:-0.19,0.6164414002968976
=#

#--------------------------------------------------------------------#
# 5. tree with substitutions/site per branch, d_i
#--------------------------------------------------------------------#
# rescale each gene tree by mean(d_l)/d_l to eliminate rate variation
# across loci, then average the pairwise genetic distances across gene
# trees (each pair of taxa: over the loci in which both are present)
D_sum = zeros(ntax, ntax)
D_n = zeros(Int, ntax, ntax)
for (D, h) in zip(Dmats[hasoutgroup], d_l)
    for a in 1:ntax, b in 1:ntax
        isnan(D[a, b]) && continue
        D_sum[a, b] += D[a, b] * (mean_d / h)
        D_n[a, b] += 1
    end
end
D_mean = D_sum ./ D_n
println("loci per pair of taxa used in the mean distances: min = ",
        minimum(D_n[i, j] for i in 1:ntax for j in 1:ntax if i != j),
        ", max = ", maximum(D_n))
# pinus: min = 267, max = 322 ; ginkgo: min = 164, max = 193
# fit these distances to the species tree topology with ordinary least
# squares. ultrametric=false (the default is true): a clock would erase
# the rate variation across lineages that we want to capture.
sub_tree = readnewick(speciestree_string)
for n in sub_tree.node
    isleaf(n) || (n.name = "")
end
rootatnode!(sub_tree, outgroup)
for e in sub_tree.edge e.length = 0.1; end # starting values
calibratefrompairwisedistances!(sub_tree, D_mean, codes;
                                verbose=false, ultrametric=false)
for n in sub_tree.node
    n.leaf && (n.name = code2letter[n.name])
end
sub_tree_string = writenewick(sub_tree)
println("tree with substitutions/site per branch:\n", sub_tree_string)
#= pinus:
"(A:0.22623565785742358,(((C:0.129977272463007,D:0.13407652421089206):
0.17468672998792792,(E:0.16498639269290707,(H:0.32060456087239814,
(F:0.14588649567995127,G:0.1308336895326273):0.021204504118894928):
0.03342822447644728):0.030205051866430785):0.01031117003260008,
B:0.25142906496098044):0.22623565785742358);"
Note: the outgroup edge and the ingroup stem edge get the same length:
their sum is identifiable from pairwise distances, not the split. The
reptile sub_tree has the same property (Homo: 0.03473, root: 0.03473).
=#
# goodness of fit of the tree distances to the mean pairwise distances
D_fit = pairwisetaxondistancematrix(sub_tree)
tl = tiplabels(sub_tree)
idx = [findfirst(==(code2letter[c]), tl) for c in codes]
@printf("max |D_mean - D_fit| = %.5f\n",
        maximum(abs.(D_mean .- D_fit[idx, idx])))
# pinus: 0.00730 ; ginkgo: 0.01638

#--------------------------------------------------------------------#
# 6. genome-wide rate, and per-branch multipliers m_i
#--------------------------------------------------------------------#
cu_total = totallength(tree_cu)     # 28.86 (pinus)
sub_total = totallength(sub_tree)   #  2.0036 (pinus)
sub_per_cu = sub_total / cu_total   # bar_r: 0.069393 (pinus),
                                    # 0.058034 (ginkgo);
                                    # reptile: 0.019526
@printf("bar_r = %.6f substitutions per site per coalescent unit\n",
        sub_per_cu)
# SimPhy baseline substitution rate, per site per generation: bar_r / 2Ne
@printf("SimPhy -su f:%.10e\n", sub_per_cu / eff_pop)
# pinus: 6.939262e-5 ; ginkgo: 5.803402e-5 ; reptile: 1.9526e-5

#= per-branch multiplier: m_i = d_i / (bar_r * tau_i), dimensionless,
   so that SimPhy's expected substitutions along branch i are
   t_gen,i * mu * m_i = (tau_i 2Ne) (bar_r / 2Ne) (d_i / (bar_r tau_i))
                      = d_i
   (the "updated" multiplier of speciestree.jl,
    makeSimPhytree_dividied_by_cu)
=#
sub_by_clade = Dict(leafset(e) => e.length for e in sub_tree.edge)
mult_by_clade = Dict{Set{String},Float64}()
println("branch\tCU\td_i\tm_i")
for e in tree_ngen.edge
    ls = leafset(e)
    tau = e.length / eff_pop
    d = sub_by_clade[ls]
    mult_by_clade[ls] = d / (sub_per_cu * tau)
    @printf("%s\t%.2f\t%.4f\t%.3f\n", join(sort(collect(ls))), tau, d,
            mult_by_clade[ls])
end
# SimPhy tree edges annotated as `generations*multiplier`, written
# directly from the tree (edges matched by the set of tips below them,
# not by their length: two internal edges can have the same length).
function simphynewick(node, parentedge)
    if node.leaf
        s = node.name
    else
        kids = [simphynewick(getchild(e), e)
                for e in node.edge if getparent(e) === node]
        s = "(" * join(kids, ",") * ")"
    end
    isnothing(parentedge) && return s * ";"
    return s * ":" * string(Int(round(parentedge.length))) * "*" *
           string(round(mult_by_clade[leafset(parentedge)], digits=6))
end
simphytree = simphynewick(getroot(tree_ngen), nothing)
println("SimPhy species tree with lineage-specific multipliers:\n",
        simphytree)
#= pinus:
branch  CU     d_i     m_i
A       5.13   0.2262  0.635   Pinus
C       0.50   0.1300  3.744   Nymphaea
D       0.50   0.1341  3.863   Brasenia
CD      3.93   0.1747  0.640
E       3.70   0.1650  0.642   Schisandra
H       2.81   0.3206  1.643   Ceratophyllum
F       2.48   0.1459  0.847   Chloranthus
G       2.48   0.1308  0.760   Liriodendron
FG      0.33   0.0212  0.926
FGH     0.89   0.0334  0.541
EFGH    0.73   0.0302  0.596
CDEFGH  0.20   0.0103  0.743
B       4.63   0.2514  0.782   Amborella
BCDEFGH 0.50   0.2262  6.518   ingroup stem
"(A:5130*0.635236,(((C:500*3.744456,D:500*3.86255):3930*0.640263,
(E:3700*0.6423,(H:2810*1.643443,(F:2480*0.847334,G:2480*0.759905):
330*0.925562):890*0.541022):730*0.596002):200*0.742625,B:4630*0.782214):
500*6.51752);"
Fastest lineages: the two water lilies (3.7-3.9x) and Ceratophyllum
(1.6x); slowest: Pinus and Schisandra (0.64x). The ingroup stem (6.5x)
is inflated by its 0.5 CU convention (3.6x in the reptile tree).
Multiplier range 0.54-3.9 at the tips, vs 0.19-9.2 for the reptile tree
(snake 9.2x, lizard 4.1x, alligator 0.19x).
ginkgo:
"(A:5860*0.590137,(((E:4370*0.647977,(H:3440*1.593824,(F:3100*0.805864,
G:3100*0.73324):340*1.356683):930*0.562365):740*0.463047,(D:500*4.654533,
C:500*4.524783):4610*0.677257):250*1.027253,B:5360*0.802687):500*6.916408);"
(water lilies 4.5-4.7x, Ceratophyllum 1.6x, Ginkgo 0.59x)
=#

#--------------------------------------------------------------------#
# 7. Seq-Gen parameters: distribution across genes of the HKY+G model
#--------------------------------------------------------------------#
# same recipe as notes/choice-seqgen-parameters.md for the reptile UCEs:
# per-gene estimates from IQ-TREE (seqgen_params.csv, written by
# collect_pergene.py), then a parametric family fitted across genes.
sg = CSV.read(joinpath(datadir, "seqgen_params.csv"), DataFrame)
println("Seq-Gen parameters, from ", nrow(sg), " genes:")
# shape alpha of the Gamma distribution of rates across sites
for family in [Normal, LogNormal, Gamma]
    d = fit(family, sg.alpha)
    println("  alpha: loglik=$(round(sum(logpdf.(d, sg.alpha)),digits=1)), $d")
end
alpha_fit = fit(Gamma, sg.alpha)
# transition/transversion ratio kappa
for family in [Normal, LogNormal, Gamma]
    d = fit(family, sg.kappa)
    println("  kappa: loglik=$(round(sum(logpdf.(d, sg.kappa)),digits=1)), $d")
end
kappa_fit = fit(LogNormal, sg.kappa)
# base frequencies: Dirichlet (one column per gene for `fit`)
freqmat = permutedims(Matrix{Float64}(sg[:, [:fA, :fC, :fG, :fT]]))
freqmat ./= sum(freqmat, dims=1) # renormalize rounded frequencies
freq_fit = fit(Dirichlet, freqmat)
@printf("  alpha ~ Gamma(shape=%.4f, scale=%.5f), mean %.4f\n",
        alpha_fit.α, alpha_fit.θ, mean(alpha_fit))
@printf("  kappa ~ LogNormal(mu=%.4f, sigma=%.4f), mean %.4f\n",
        kappa_fit.μ, kappa_fit.σ, mean(kappa_fit))
@printf("  base frequencies A,C,G,T ~ Dirichlet(%s), mean %s\n",
        join(round.(freq_fit.alpha, digits=2), ", "),
        join(round.(mean(freq_fit), digits=3), ", "))
#= pinus (353 genes):
   alpha ~ Gamma(shape=12.20, scale=0.0412), mean 0.503
   kappa ~ LogNormal(mu=1.314, sigma=0.158), mean 3.77
   base frequencies ~ Dirichlet(86.9, 61.7, 77.6, 84.9),
   mean 0.279,0.198,0.249,0.273  (ginkgo set: nearly identical)
   reptile UCEs: alpha ~ Gamma(3.267, 0.109), mean 0.356;
   kappa ~ LogNormal(1.4216, 0.2798), mean 4.31;
   frequencies ~ Dirichlet(66.6, 38.4, 38.6, 67.1), mean 0.316,0.182,0.183,0.319
=#

# save everything needed by simulation.jl for this taxon set
open(joinpath(datadir, "simphy_tree_basal_angio.txt"), "w") do io
    println(io, "# taxon set: ", setname, "; ",
            join(["$(r.letter)=$(r.species) ($(r.code))"
                  for r in eachrow(taxa)], ", "))
    println(io, "# loci used: ", nloci, " for the species tree, ",
            length(d_l), " for the substitution rates")
    @printf(io, "# -su f:%.10e\n", sub_per_cu / eff_pop)
    @printf(io, "# -hl ln:%.6f,%.6f\n", mu_gene, sigma_gene)
    println(io, "# species tree, coalescent units:\n", writenewick(tree_cu))
    println(io, "# species tree, substitutions/site:\n", sub_tree_string)
    println(io, "# SimPhy species tree, no lineage rate variation ",
            "(generations, 2Ne=$eff_pop):\n", speciestree_ngen)
    println(io, "# SimPhy species tree with lineage-specific multipliers:\n",
            simphytree)
    @printf(io, "# Seq-Gen: alpha ~ Gamma(shape=%.4f, scale=%.5f)\n",
            alpha_fit.α, alpha_fit.θ)
    @printf(io, "# Seq-Gen: kappa ~ LogNormal(mu=%.4f, sigma=%.4f)\n",
            kappa_fit.μ, kappa_fit.σ)
    println(io, "# Seq-Gen: base frequencies A,C,G,T ~ Dirichlet(",
            join(round.(freq_fit.alpha, digits=3), ","), ")")
end
