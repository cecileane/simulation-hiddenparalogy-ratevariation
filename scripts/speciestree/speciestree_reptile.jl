#=
code used to build a "true" species tree used by SimPhy to simulate gene trees.
Its branch lengths are in coalescent units.

- taken from the Crawford et al. (2012) data.
  nexus alignments downloaded from
  [dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8gm85)
  which corresponds to Brown & Thomson (2017).
- gene trees were estimated with IQ-TREE,
  using this command and option for 1000 ultra-fast bootstraps:
  `iqtree2 -S folder --prefix loci -T auto -B 1000 -wbtl`
  which creates in `folder`
  * a single file `loci.treefile` for the ML tree, one per locus,
  * and 1 file `loci.ufboot` containing all bootstrap trees (1000 per locus).
  followed by julia script separate-boot-bygene.jl (reptiles repo)
  to get 1 bootstrap tree file per gene, run from `crawford`:
  `julia ../scripts/separate-boot-bygene.jl -B 1000 -o iqtree/bootstrap iqtree/IQ_0??/loci.ufboot`
  Then, a file `BSlistfiles` was created with 1 row per locus,
  giving the name of each bootstrap file.
- species tree estimated with ASTRAL with IQTree input, using
  `java -jar astral.5.7.7.jar -i merged_iqtree.treefile -b BSlistfiles -r 1000`
  The output file has 1002 lines:
  * 1000 bootstrap trees first, then
  * ASTRAL tree with bootstrap support as node names
  * same ASTRAL tree + edge lengths in coalescent units
- edge lengths from ASTRAL (in coalescent units) were rounded to 2 digits
- missing edge lengths were assigned in such a way that the tree is ultrametric.

A network is ultrametric if it has a "height", such that the length
of any path from the root to any tip is equal to this height.

Note: On Jun 11 2025, PhyloPlots and QuartetNetworkGoodnessFit
does not match the newest PhyloNetwork version 1.1.0. Therefore, 
this script should be run under PhyloNetworks < 1.1.0. 
=#

using PhyloNetworks
using PhyloPlots
using QuartetNetworkGoodnessFit # has a function to ultrametrized a network
using PhyloPlots, RCall
R"layout"([1 2])
R"par"(mar=[0,0,0,0])
include("scripts/utilities.jl") # run from repo root

# below: tree from crawford data, ASTRAL on IQTree, copy-pasted from
# reptiles repo: estimatednets_collapsed.csv#L44
speciestree_string = "(Homo,(((Chrysemys,Pelomedusa)100.0:1.0288670824548658,((Crocodylus,Alligator)100.0:1.7082993293412243,(Taeniopygia,Gallus)93.2:1.6647117465735648)0.0:0.16825214855685644)0.0:0.17931560798411753,(Sphenodon,(Anolis,Pantherophis):2.2370887285658227)100.0:0.20028294826845958)0.0);"
tree = readnewick(speciestree_string)

# round edge lengths to avoid 15 digits: not significant for simulations
for e in tree.edge
    e.length = round(e.length, digits=2)
end
# the plot below shows that the tree lacks external edge lengths,
# also lacks a length for the edge leading to the ingroup:
plot(tree, showedgelength=true, useedgelength=true);
R"mtext"("edge lengths in coalescent units, estimated",
    line=-2, side=1, cex=0.5)
# assign 0.5 coalescent units to the ingroup stem edge:
# find this edge first
findfirst(x -> x.length == -1.0 && !getchild(x).leaf, tree.edge) # 18
tree.edge[18].length = 0.5 # assign 0.5 coalescent unit as length
# assign 0.5 coalescent unit to the external edge going to Anolis
# find this edge first
findfirst(x -> getchild(x).name == "Anolis", tree.edge) # 14
tree.edge[14].length = 0.5 # assign 0.5 coalescent unit as length
# ultrametrize this tree by assigning lengths to its extrenal branches,
# which don't have any length so far.
QuartetNetworkGoodnessFit.ultrametrize!(tree, true) # verbose=true
plot(tree, showedgelength=true, useedgelength=true); # looks good
R"mtext"("edge lengths in coalescent units, estimated " *
    "else assigned to make the tree ultrametric",
    line=-3, side=1, cex=0.5)

# prune Sphenodon and 1 turtle, Pelomedusa
deleteleaf!(tree, "Pelomedusa")
deleteleaf!(tree, "Sphenodon")
for e in tree.edge e.length = round(e.length, digits=2); end

writenewick(tree)

# we get this below, which was copy-pasted into the main readme file
speciestree_string = "(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88)100.0:1.71,(Taeniopygia:0.93,Gallus:0.93)93.2:1.66)0.0:0.17,Chrysemys:2.76)0.0:0.18,(Anolis:0.5,Pantherophis:0.5):2.44)0.0:0.5);"
tree = readnewick(speciestree_string)
# remove bootstrap values from node names (SimPhy doesn't accept them)
for n in tree.node
    isleaf(n) && continue
    n.name = ""
end
writenewick(tree)
"(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);"

# ratevariation.jl (reptiles repo): tree with subs/site averaged across genes
sub_tree_string = "(Homo:0.03472578056145014,((((Crocodylus:0.003700978635732742,Alligator:0.0032363141626990384):0.013424996872103507,(Taeniopygia:0.0219417335205793,Gallus:0.01858077060211097):0.013265606603973844):0.001170217330315521,Chrysemys:0.018550404584983873):0.0017655989402352748,(Anolis:0.03989844563755188,Pantherophis:0.08984619023323578):0.046478889220648016):0.0347294194947231);"
sub_tree = readnewick(sub_tree_string)

# Let's calculate substitutions / tree length in coalescent: 
cu_total = totallength(tree) # total length in coalescence = 17.48
# next: total length in substitition
sub_total = totallength(sub_tree) # 0.341315346400343 
sub_per_cu = sub_total / cu_total # 0.019526049565237014 

#= 
 Lineage-specific substitution rate per coalescent unit: 
 sub / cu = substitution_length / CU_length, per specfic branch 
 For example, for external edge to homo:
 0.03472578056145014 (subs/site) / 3.44 (CU) = 0.010094703651584344

 Calculate substitution rate per site per generation for each lineage
=# 
eff_pop = 1000 # = 2Ne, diploid effective population size 

# tree with edge lengths in number of generations: gen = cu * 2Ne
tree_ngen = deepcopy(tree)
for e in tree_ngen.edge
    e.length *= eff_pop
end
replace(writenewick(tree_ngen), ".0" => "")
"(Homo:3440,((((Crocodylus:880,Alligator:880):1710,(Taeniopygia:930,Gallus:930):1660):170,Chrysemys:2760):180,(Anolis:500,Pantherophis:500):2440):500);"

#--------------------------------------------------------------------# 
# Estimate per-branc-multiplier m_i for lineage-specific rate variation 
#--------------------------------------------------------------------#  

#= Simphy Substitution Rate Multiplier (*) 
The substitution rate multiplier in SimPhy is dimensionless. 
It scales the baseline substitution rate. 
=#

# After talk on April 17th 2026 
#= 
Below we present one out-dated and one updated approach
 we used in the simulation pipeline: 
 1) The out-dated approach used incorrect multiplier m_i_old 
 2) The updated approach used unit-less multilier m_i 

 Notation for one branch i:
   tau_i = branch length in coalescent units (from `tree`)
   d_i = branch length in substitutions/site (from `sub_tree`)
   t_i = tau_i * 2Ne = branch length in generations
   r_i = d_i / tau_i = per-branch subs/site per CU
   bar_r = sub_per_cu = sum(d_i) / sum(tau_i) = tree-wide avg
   mu = -su baseline = bar_r / 2Ne (option 3 from above) 

SimPhy computes expected subs/site per branch as T_i * mu * m_i. 

Out-dated approach described as below:  
m_i_old = makeSimPhytree(1/sub_per_cu) used m_b = s_b / bar_r.
   -> T_i * mu * m_i_old
    = tau_i * 2Ne * (bar_r / 2Ne) * (d_i / bar_r) 
    = tau_i * d_i, NOT d_i.
   Each branch's simulated substitution length is scaled by its
   own CU length. Thus, it does not recover d_i 
=# 

# SimPhy tree edges annotated as `generations*substitutions/cu`
# sub_factor: factor to multiple substitution per site / coalescent units
function makeSimPhytree(sub_factor)
    nwk = writenewick(tree_ngen) # (Homo:3440.0,((((Crocodylus:880.0,...
    for (e_g,e_s) in zip(tree_ngen.edge, sub_tree.edge)
        cn = getchild(e_g).name
        cn == getchild(e_s).name ||
            error("edges don't come in the order, an assumption made here")
        if cn == "" # internal node have no name
            cn == "\\)" # its parent edge length would be after ')'
        end
        tofind = Regex("($cn:" * string(Int(e_g.length)) * ")\\.0")
        toreplace = SubstitutionString("\\1*" * string(e_s.length * sub_factor))
        nwk = replace(nwk, tofind => toreplace)
    end
    return nwk
end

#=
The baseline substitution rate is set using the `-su` parameter (or `-SU`) 
in the SimPhy configuration or command-line. The value represents 
**substitutions per site per generation**. 
SimPhy uses this to compute branch lengths on gene trees in expected 
substitutions per site, before passing to the sequence simulator.
=# 
# sub_per_cu / eff_pop as the baseline: -su 0.000019526049565237014
sub_per_cu / eff_pop # 1.9526049565237015e-5
# Baseline rate used for root population; no root-specific multiplier needed.
# SimPhy tree for this:
simpytree = makeSimPhytree(1 / sub_per_cu)
"(Homo:3440*1.7784334944675035,((((Crocodylus:880*0.18954057365099278,Alligator:880*0.165743416346785):1710*0.6875429065797596,(Taeniopygia:930*1.1237159594044575,Gallus:930*0.951588827019626):1660*0.679379951364618):170*0.05993108469820257,Chrysemys:2760*0.9500336728638588):180*0.09042274190364852,(Anolis:500*2.043344482162159,Pantherophis:500*4.601350105819277):2440*2.380352926246597):500*1.778619857472514);"

# use simphy tree below, with -su 0.000019526049565237014:
replace_tips_with_letters(simpytree)
"(A:3440*1.7784334944675035,((((B:880*0.18954057365099278,C:880*0.165743416346785):1710*0.6875429065797596,(D:930*1.1237159594044575,E:930*0.951588827019626):1660*0.679379951364618):170*0.05993108469820257,F:2760*0.9500336728638588):180*0.09042274190364852,(G:500*2.043344482162159,H:500*4.601350105819277):2440*2.380352926246597):500*1.778619857472514);"

#= Below shows a updated multiplier, 
    which is the FINAL multiplier used in the simulation pipeline 
The updated per-branch multiplier below: 
New makeSimPhytree_dividied_by_cu(1/sub_per_cu) uses
   m_i = d_i / (bar_r * tau_i) 
   Thus, T_i * mu * m_i = 
    = (tau_i * 2Ne) * (bar_r / 2Ne) * (d_i / bar_r * tau_i) 
    = d_i  
   m_i is also the dimensionally-correct object: a unit-free
   per-branch rate divided by the tree-wide average rate, which
   matches SimPhy's interpretation of the `*` multiplier in -hs
   (a rate multiplier with mean one).

 Note: the line `cu = e_g.length / 1000` hard-codes 2Ne = 1000.
 If `eff_pop` (line 88) is changed, use `cu = e_g.length/eff_pop`
 instead so the two stay in sync. 
 This is the reference effective population size.  
=# 
# global substitution rate: 
sub_per_cu / eff_pop # 1.9526049565237015e-5: this is bar_r
function makeSimPhytree_dividied_by_cu(sub_factor)
    nwk = writenewick(tree_ngen) # (Homo:3440.0,((((Crocodylus:880.0,...
    for (e_g,e_s) in zip(tree_ngen.edge, sub_tree.edge)
        cn = getchild(e_g).name
        cu = e_g.length / 1000 
        cn == getchild(e_s).name ||
            error("edges don't come in the order, an assumption made here")
        if cn == "" # internal node have no name
            cn == "\\)" # its parent edge length would be after ')'
        end
        tofind = Regex("($cn:" * string(Int(e_g.length)) * ")\\.0")
        toreplace = SubstitutionString(
            "\\1*" * string((e_s.length * sub_factor) / cu))
        # sub * sub_factor / (cu)
        nwk = replace(nwk, tofind => toreplace)
    end
    return nwk
end

simphytree = makeSimPhytree_dividied_by_cu(1 / sub_per_cu)
"(Homo:3440*0.5169864809498557,((((Crocodylus:880*0.2153870155124918,Alligator:880*0.18834479130316475):1710*0.4020718751928419,(Taeniopygia:930*1.2082967305424273,Gallus:930*1.0232137924942215):1660*0.409265030942541):170*0.35253579234236804,Chrysemys:2760*0.344215098863717):180*0.5023485661313807,(Anolis:500*4.086688964324318,Pantherophis:500*9.202700211638554):2440*0.9755544779699168):500*3.557239714945028);"

# This below tree is used for Lineage substitution rate variation 
# in simulation.jl: 
# m_i is followed after "*" for each branch 
replace_tips_with_letters(simphytree) 
"(A:3440*0.5169864809498557,((((B:880*0.2153870155124918,C:880*0.18834479130316475):1710*0.4020718751928419,(D:930*1.2082967305424273,E:930*1.0232137924942215):1660*0.409265030942541):170*0.35253579234236804,F:2760*0.344215098863717):180*0.5023485661313807,(G:500*4.086688964324318,H:500*9.202700211638554):2440*0.9755544779699168):500*3.557239714945028);"



