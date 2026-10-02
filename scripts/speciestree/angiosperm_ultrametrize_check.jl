# check of scripts/speciestree/calibrate.r: rebuild the calibrated angiosperm tree with
# PhyloNetworks + ultrametrize!, as in scripts/speciestree_reptile.jl for the reptile tree.
# start from the ASTRAL tree (internal branches in coalescent units), set the two
# branches the calibration gives (water-lily tips 2.94, ingroup stem 1.64), and let
# ultrametrize! fill in the other external branches. The result should equal the
# tree computed by hand in scripts/speciestree/calibrate.r.
# usage, from the repo root: julia --project=. scripts/speciestree/angiosperm_ultrametrize_check.jl
using PhyloNetworks
using QuartetNetworkGoodnessFit # for ultrametrize!

# ASTRAL tree, last line of data/kew_a353/astral/species.tre, rounded to 2 digits
speciestree_string = "(Pinus,(((Nymphaea,Brasenia):3.77,(Schisandra,(Ceratophyllum,(Chloranthus,Liriodendron):0.44):0.95):0.79):0.16,Amborella):0.0);"
tree = readnewick(speciestree_string)
rootatnode!(tree, "Pinus")

# from scripts/speciestree/calibrate.r: k = 22 CU per substitution/site
k = 22
cu_NBtip = k * (0.131 + 0.136) / 2    # 2.94: Nymphaea and Brasenia tips
H_ingroup = 0.16 + 3.77 + cu_NBtip    # 6.87: angiosperm crown to tips
cu_stem = (k * 0.454 - H_ingroup) / 2   # 1.56: 2*stem + H_ingroup = k*0.454

# assign the water-lily tips
for name in ("Nymphaea", "Brasenia")
    i = findfirst(e -> getchild(e).name == name, tree.edge)
    tree.edge[i].length = cu_NBtip
end
# assign the ingroup stem: the edge from the root to the ingroup
i = findfirst(e -> getparent(e) === getroot(tree) && !getchild(e).leaf, tree.edge)
tree.edge[i].length = cu_stem

# the other external edges have no length yet: ultrametrize! assigns them
QuartetNetworkGoodnessFit.ultrametrize!(tree, true) # verbose = true
for e in tree.edge e.length = round(e.length, digits=2); end
println("julia:  ", writenewick(tree))

