#=
code used to build a "true" species tree used by SimPhy to simulate gene trees.
Its branch lengths are in coalescent units.

- taken from the crawford data, tree estimated with ASTRAL with IQTree input
- edge lengths were rounded to 2 digits
- missing edge lengths were assigned in such a way that the tree is ultrametric.

A network is ultrametric if it has a "height", such that the length
of any path from the root to any tip is equal to this height.
=#

using PhyloNetworks
using PhyloPlots
using QuartetNetworkGoodnessFit # has a function to ultrametrized a network
include("utilities.jl") # Use the function to replace all tips

# below: tree from crawford data, ASTRAL on IQTree, copy-pasted from
# https://github.com/cecileane/reptiles/blob/main/estimatednets_collapsed.csv#L44
speciestree_string = "(Homo,(((Chrysemys,Pelomedusa)100.0:1.0288670824548658,((Crocodylus,Alligator)100.0:1.7082993293412243,(Taeniopygia,Gallus)93.2:1.6647117465735648)0.0:0.16825214855685644)0.0:0.17931560798411753,(Sphenodon,(Anolis,Pantherophis):2.2370887285658227)100.0:0.20028294826845958)0.0);"
tree = readTopology(speciestree_string)
# round edge lengths to avoid 15 digits: not significant for simulations
for e in tree.edge
    e.length = round(e.length, digits=2)
end
# the plot below shows that the tree lacks external edge lengths,
# also lacks a length for the edge leading to the ingroup:
plot(tree, :R, showEdgeLength=true, useEdgeLength=true);
# assign 0.5 coalescent units to the ingroup stem edge:
# find this edge first
findfirst(x -> x.length == -1.0 && !PhyloNetworks.getChild(x).leaf, tree.edge) # 18
tree.edge[18].length = 0.5 # assign 0.5 coalescent unit as length
# assign 0.5 coalescent unit to the external edge going to Anolis
# find this edge first
findfirst(x -> PhyloNetworks.getChild(x).name == "Anolis", tree.edge) # 14
tree.edge[14].length = 0.5 # assign 0.5 coalescent unit as length
plot(tree, :R, showEdgeLength=true, useEdgeLength=true);
# looks good

# ultrametrize this tree by assigning lengths to its extrenal branches,
# which don't have any length so far.
QuartetNetworkGoodnessFit.ultrametrize!(tree, true) # verbose=true
plot(tree, :R, showEdgeLength=true, useEdgeLength=true); # looks good

# prune Sphenodon and 1 turtle, Pelomedusa
deleteleaf!(tree, "Pelomedusa")
deleteleaf!(tree, "Sphenodon")

writeTopology(tree, round=true)

# we get this below, which was copy-pasted into the main readme file
"(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88)100.0:1.71,(Taeniopygia:0.93,Gallus:0.93)93.2:1.66)0.0:0.17,Chrysemys:2.76)0.0:0.18,(Anolis:0.5,Pantherophis:0.5):2.44)0.0:0.5);"
# has bootstrap values as node names. Let's remove them: SimPhy doesn't like them
"(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);"

# see reptiles repo, analysis in ratevariation.jl: tree with substitutions/site averaged across genes:
"(Homo:0.03472578056145014,((((Crocodylus:0.003700978635732742,Alligator:0.0032363141626990384):0.013424996872103507,(Taeniopygia:0.0219417335205793,Gallus:0.01858077060211097):0.013265606603973844):0.001170217330315521,Chrysemys:0.018550404584983873):0.0017655989402352748,(Anolis:0.03989844563755188,Pantherophis:0.08984619023323578):0.046478889220648016):0.0347294194947231);"

# Let's calculate substitutions / tree length in coalescent: 
cu_total = sum(e.length for e in tree.edge if e.length > 0) # total length in coalescence = 17.48 
sub_tree = readTopology(sub_tree) # total length in substitition 
sub_total = sum(e.length for e in sub_tree.edge if e.length > 0) # 0.341315346400343 
sub_per_cu = sub_total / cu_total # 0.019526049565237014 

#= 
 Overall substitution rate per coalescent unit:
 tree length in substitutions / tree length in coalescent = 0.019526049565237014
 This is calculated using the overall branch length based on substition rate / overall branch length in coalescent units, across all branches. 

 Lineage-specific substitution rate per coalescent unit: 
 sub / cu = substitution_length / CU_length, per specfic branch 
 For example, for homo, 0.03472578056145014 (substitution rate) / 3.44 (coalescece unit) for homo branch length = 0.010094703651584344 
=# 
"(Homo:0.0100947,((((Crocodylus:0.0042057,Alligator:0.0036776):0.0078509,(Taeniopygia:0.0235933,Gallus:0.0199793):0.0079913):0.0068836,Chrysemys:0.0067212):0.0098089,(Anolis:0.0797969,Pantherophis:0.1796924):0.0190487):0.0694588);"

# Replace tip labels on the final tree with letters: 
tre_string1 = "(Homo:3.44,((((Crocodylus:0.88,Alligator:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);"
replace_tips_with_letters(tre_string1) # see utilities.jl 
"(A:3.44,((((B:0.88,C:0.88):1.71,(D:0.93,E:0.93):1.66):0.17,F:2.76):0.18,(G:0.5,H:0.5):2.44):0.5);"

#= the branch length in coalescent unit by substition rate / CU
 All the information below is based on: https://github.com/adamallo/SimPhy/wiki/Manual 
 In Simphy, * = substition rate multiplier, which is in 
=# 
tre_string2 = "(Homo:3.44*0.0100947,((((Crocodylus:0.88*0.0042057,Alligator:0.88*0.0036776):1.71*0.0078509,(Taeniopygia:0.93*0.0235933,Gallus:0.93*0.0199793):1.66*0.0079913):0.17*0.0068836,Chrysemys:2.76*0.0067212):0.18*0.0098089,(Anolis:0.5*0.0797969,Pantherophis:0.5*0.1796924):2.44*0.0190487):0.5*0.0694588);"
replace_tips_with_letters(tre_string2)
"(A:3.44*0.0100947,((((B:0.88*0.0042057,C:0.88*0.0036776):1.71*0.0078509,(D:0.93*0.0235933,E:0.93*0.0199793):1.66*0.0079913):0.17*0.0068836,F:2.76*0.0067212):0.18*0.0098089,(G:0.5*0.0797969,H:0.5*0.1796924):2.44*0.0190487):0.5*0.0694588);"

# Now manually combining coalescent units and substitution multiplier in SimPhy's format: 
# The tree without lineage variations: 
"(A:3.44,((((B:0.88,C:0.88):1.71,(D:0.93,E:0.93):1.66):0.17,F:2.76):0.18,(G:0.5,H:0.5):2.44):0.5);"
"(A:3.44*0.0195260,((((B:0.88*0.0195260,C:0.88*0.0195260):1.71,(D:0.93*0.0195260,E:0.93*0.0195260):1.66*0.0195260):0.17*0.0195260,F:2.76*0.0195260):0.18*0.0195260,(G:0.5*0.0195260,H:0.5*0.0195260):2.44*0.0195260):0.5*0.0195260);"

"(A:3.44*0.0195260,((((B:0.88*0.0195260,C:0.88*0.0195260):1.71,(D:0.93*0.0195260,E:0.93*0.0195260):1.66*0.0195260):0.17*0.0195260,F:2.76*0.0195260):0.18*0.0195260,(G:0.5*0.0195260,H:0.5*0.0195260):2.44*0.0195260):0.5*0.0195260);"

# The tree with lineage variations: 
"(A:3.44*0.0100947,((((B:0.88*0.0042057,C:0.88*0.0036776):1.71*0.0078509,(D:0.93*0.0235933,E:0.93*0.0199793):1.66*0.0079913):0.17*0.0068836,F:2.76*0.0067212):0.18*0.0098089,(G:0.5*0.0797969,H:0.5*0.1796924):2.44*0.0190487):0.5*0.0694588);"
