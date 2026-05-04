#=
code used to build a "true" species tree used by SimPhy to simulate gene trees.
Its branch lengths are in coalescent units.

- taken from the crawford data, tree estimated with ASTRAL with IQTree input
- edge lengths were rounded to 2 digits
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
include("utilities.jl") # defines replace_tips_with_letters, used below

# below: tree from crawford data, ASTRAL on IQTree, copy-pasted from
# https://github.com/cecileane/reptiles/blob/main/estimatednets_collapsed.csv#L44
speciestree_string = "(Homo,(((Chrysemys,Pelomedusa)100.0:1.0288670824548658,((Crocodylus,Alligator)100.0:1.7082993293412243,(Taeniopygia,Gallus)93.2:1.6647117465735648)0.0:0.16825214855685644)0.0:0.17931560798411753,(Sphenodon,(Anolis,Pantherophis):2.2370887285658227)100.0:0.20028294826845958)0.0);"
tree = readnewick(speciestree_string)

# round edge lengths to avoid 15 digits: not significant for simulations
for e in tree.edge
    e.length = round(e.length, digits=2)
end
# the plot below shows that the tree lacks external edge lengths,
# also lacks a length for the edge leading to the ingroup:
plot(tree, showedgelength=true, useedgelength=true);
R"mtext"("edge lengths in coalescent units, estimated", line=-2, side=1, cex=0.5)
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
R"mtext"("edge lengths in coalescent units, estimated else assigned to make the tree ultrametric",
        line=-3, side=1, cex=0.5)

# prune Sphenodon and 1 turtle, Pelomedusa
deleteleaf!(tree, "Pelomedusa")
deleteleaf!(tree, "Sphenodon")
for e in tree.edge e.length = round(e.length, digits=2); end

writenewick(tree)
# we get this below, which was copy-pasted into the main readme file
speciestree_string = "(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88)100.0:1.71,(Taeniopygia:0.93,Gallus:0.93)93.2:1.66)0.0:0.17,Chrysemys:2.76)0.0:0.18,(Anolis:0.5,Pantherophis:0.5):2.44)0.0:0.5);"
tree = readnewick(speciestree_string)
# has bootstrap values as node names. Let's remove them: SimPhy doesn't like them
for n in tree.node
    isleaf(n) && continue
    n.name = ""
end
writenewick(tree)
"(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);"

# see reptiles repo, analysis in ratevariation.jl: tree with substitutions/site averaged across genes:
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
 0.03472578056145014 (substitutions/site) / 3.44 (coalescent unit) = 0.010094703651584344 

 Now, calculate substitution rate per site per generation specific to each lineage
=# 
eff_pop = 1000 # = 2Ne, diploid effective population size 

# tree with edge lengths in number of generations: gen = cu * 2Ne
tree_ngen = deepcopy(tree)
for e in tree_ngen.edge
    e.length *= eff_pop
end
replace(writenewick(tree_ngen), ".0" => "")
"(Homo:3440,((((Crocodylus:880,Alligator:880):1710,(Taeniopygia:930,Gallus:930):1660):170,Chrysemys:2760):180,(Anolis:500,Pantherophis:500):2440):500);"

# tree for SimPhy, with edges annotated like this: `generations*substitutions/cu`
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
makeSimPhytree(1.0)
"(Homo:3440*0.03472578056145014,((((Crocodylus:880*0.003700978635732742,Alligator:880*0.0032363141626990384):1710*0.013424996872103507,(Taeniopygia:930*0.0219417335205793,Gallus:930*0.01858077060211097):1660*0.013265606603973844):170*0.001170217330315521,Chrysemys:2760*0.018550404584983873):180*0.0017655989402352748,(Anolis:500*0.03989844563755188,Pantherophis:500*0.08984619023323578):2440*0.046478889220648016):500*0.0347294194947231);"
# Using this species tree above, the SimPhy configuration will need to specify a
# baseline rate of 1/(2Ne) substitutions/generation: -su 0.001

# alternatively, we could multiply lineage-specific rate multipliers by 1/(2Ne)
#   these numbers will be tiny. In that case, SimPhy should use: -su 1
makeSimPhytree(1/eff_pop)
"(Homo:3440*3.472578056145014e-5,((((Crocodylus:880*3.700978635732742e-6,Alligator:880*3.2363141626990382e-6):1710*1.3424996872103508e-5,(Taeniopygia:930*2.19417335205793e-5,Gallus:930*1.858077060211097e-5):1660*1.3265606603973845e-5):170*1.1702173303155208e-6,Chrysemys:2760*1.8550404584983873e-5):180*1.7655989402352749e-6,(Anolis:500*3.989844563755188e-5,Pantherophis:500*8.984619023323579e-5):2440*4.647888922064802e-5):500*3.47294194947231e-5);"

# Both species tree still have a problem: we also needs a substitution rate
# multiplier for the (infinite) root population.
# Let's use the weighted average already calculated above:
# substitions per site / coalescent units
sub_per_cu # 0.019526049565237014
sub_per_cu / eff_pop # 1.9526049565237015e-5

# we get these 2 trees, with a root-edge multiplier added manually:
# 1. assuming a baseline -su 0.001:
simphytree_op1 = "(Homo:3440*0.03472578056145014,((((Crocodylus:880*0.003700978635732742,Alligator:880*0.0032363141626990384):1710*0.013424996872103507,(Taeniopygia:930*0.0219417335205793,Gallus:930*0.01858077060211097):1660*0.013265606603973844):170*0.001170217330315521,Chrysemys:2760*0.018550404584983873):180*0.0017655989402352748,(Anolis:500*0.03989844563755188,Pantherophis:500*0.08984619023323578):2440*0.046478889220648016):500*0.0347294194947231)*0.019526049565237014;"

replace_tips_with_letters(simphytree_op1)
"(A:3440*0.03472578056145014,((((B:880*0.003700978635732742,C:880*0.0032363141626990384):1710*0.013424996872103507,(D:930*0.0219417335205793,E:930*0.01858077060211097):1660*0.013265606603973844):170*0.001170217330315521,F:2760*0.018550404584983873):180*0.0017655989402352748,(G:500*0.03989844563755188,H:500*0.08984619023323578):2440*0.046478889220648016):500*0.0347294194947231)*0.019526049565237014;"

# 2. assuming a baseline -su 1:
simphytree_op2 = "(Homo:3440*3.472578056145014e-5,((((Crocodylus:880*3.700978635732742e-6,Alligator:880*3.2363141626990382e-6):1710*1.3424996872103508e-5,(Taeniopygia:930*2.19417335205793e-5,Gallus:930*1.858077060211097e-5):1660*1.3265606603973845e-5):170*1.1702173303155208e-6,Chrysemys:2760*1.8550404584983873e-5):180*1.7655989402352749e-6,(Anolis:500*3.989844563755188e-5,Pantherophis:500*8.984619023323579e-5):2440*4.647888922064802e-5):500*3.47294194947231e-5)*1.9526049565237015e-5;"

replace_tips_with_letters(simphytree_op2) 
"(A:3440*3.472578056145014e-5,((((B:880*3.700978635732742e-6,C:880*3.2363141626990382e-6):1710*1.3424996872103508e-5,(D:930*2.19417335205793e-5,E:930*1.858077060211097e-5):1660*1.3265606603973845e-5):170*1.1702173303155208e-6,F:2760*1.8550404584983873e-5):180*1.7655989402352749e-6,(G:500*3.989844563755188e-5,H:500*8.984619023323579e-5):2440*4.647888922064802e-5):500*3.47294194947231e-5)*1.9526049565237015e-5;"

# 3. yet another alternative would be to use the average rate / generation
# sub_per_cu / eff_pop as the baseline: -su 0.000019526049565237014
sub_per_cu / eff_pop # 1.9526049565237015e-5
# This baseline rate should then be used for the root population, so there should
# be *no* need to specify a rate multiplier specific to the root population.
# SimPhy tree for this:
simpytree = makeSimPhytree(1 / sub_per_cu)
"(Homo:3440*1.7784334944675035,((((Crocodylus:880*0.18954057365099278,Alligator:880*0.165743416346785):1710*0.6875429065797596,(Taeniopygia:930*1.1237159594044575,Gallus:930*0.951588827019626):1660*0.679379951364618):170*0.05993108469820257,Chrysemys:2760*0.9500336728638588):180*0.09042274190364852,(Anolis:500*2.043344482162159,Pantherophis:500*4.601350105819277):2440*2.380352926246597):500*1.778619857472514);"

# use simphy tree below, with -su 0.000019526049565237014:
replace_tips_with_letters(simpytree)
"(A:3440*1.7784334944675035,((((B:880*0.18954057365099278,C:880*0.165743416346785):1710*0.6875429065797596,(D:930*1.1237159594044575,E:930*0.951588827019626):1660*0.679379951364618):170*0.05993108469820257,F:2760*0.9500336728638588):180*0.09042274190364852,(G:500*2.043344482162159,H:500*4.601350105819277):2440*2.380352926246597):500*1.778619857472514);"
