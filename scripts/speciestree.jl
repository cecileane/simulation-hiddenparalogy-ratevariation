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
# below: tree from crawford data, ASTRAL on IQTree, copy-pasted from
# https://github.com/cecileane/reptiles/blob/main/estimatednets_collapsed.csv#L44
speciestree_string = "(Homo,(((Chrysemys,Pelomedusa)100.0:1.0288670824548658,((Crocodylus,alligator_mississippiensis)100.0:1.7082993293412243,(Taeniopygia,Gallus)93.2:1.6647117465735648)0.0:0.16825214855685644)0.0:0.17931560798411753,(Sphenodon,(Anolis,Pantherophis):2.2370887285658227)100.0:0.20028294826845958)0.0);"
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
"(Homo:3.44,((((Crocodylus:0.88,alligator_mississippiensis:0.88)100.0:1.71,(Taeniopygia:0.93,Gallus:0.93)93.2:1.66)0.0:0.17,Chrysemys:2.76)0.0:0.18,(Anolis:0.5,Pantherophis:0.5):2.44)0.0:0.5);"
# has bootstrap values as node names. Let's remove them: SimPhy doesn't like them
"(Homo:3.44,((((Crocodylus:0.88,alligator_mississippiensis:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);"
