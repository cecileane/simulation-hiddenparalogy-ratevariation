using PhyloNetworks
using PhyloPlots

#Don't know where these came from which may be an issue when caculating substitution lengths of branches will still be used as the primary tree for the time being
tree = readTopology("(outgroup:5.0,((((crocodilia:1.84,testudines:1.844):0.182,bird:2.778):0.405,squamata:3.878):1.442));")
#2.17 average coalescent unit length
#other lenghts determined by tree
#coalescent units
plot(tree, :R, showEdgeLength=true, useEdgeLength=true);
#Todo: add substitution units here


#Added extra trees from the fast tree (246) trees chiari file incase the location
#of the tree becomes an issue at a later time.

#branch lengths look appropriate
ftree1 = readTopology("(Homo:0.05637,(Monodelphis:0.08214,(Xenopus:0.29262,Ornithorhynchus:0.14720)0.816:0.04271)0.933:0.04555,((Anolis:0.08332,podarcis:0.05018)0.964:0.04358,(caiman:0.07705,((Gallus:0.11749,Taeniopygia:0.07573)0.996:0.08505,caretta:0.02965)0.898:0.02392)0.680:0.02368)0.992:0.07142);")
