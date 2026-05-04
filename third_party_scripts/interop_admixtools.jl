# The below script is adapted from: 
# https://github.com/JuliaPhylo/PhyloUtilities?tab=MIT-1-ov-file 

"""
Utilities for interoperability between JuliaPhylo packages and the R package
[admixtools](https://github.com/uqrmaie1/admixtools), which can fit
admixture graphs using f-statistics.

by: Cécile Ané, 2025, using PhyloNetworks v1.0.0.

Functions to convert a HybridNetwork (Julia) to/from an edge list (R):

    net_to_edgelist
    edgelist_to_net


See their documentation below for details and examples.
The examples assume that R and admixtools have already been installed,
like this for example:

```{r}
install.packages("devtools")
devtools::install_github("uqrmaie1/admixtools")
library("admixtools")
```
"""

"""
    net_to_edgelist(net::HybridNetwork)
    net_to_edgelist(net_newickfilename::AbstractString)

Given a HybridNetwork object `net` or the path to a file name containing a
network in extended, Newick format, `net_to_edgelist` returns a DataFrame
describing this network in the edge list format used by `admixtools`.
This data frame has 1 row per edge.
Its columns are:
1. name of the source node
2. name of the target node
3. edge type: "admix" (for hybrid edges) or "edge" (for tree edges)
4. weight: inheritance γ for hybrid edges; length for tree edges.

```@repl
julia> using DataFrames, PhyloNetworks

julia> include("scripts/interop_admixtools.jl")

julia> net = readnewick("((t4:1.504,((t3:1.268,t1:1.268):1.057,#H8:0.0::0.4):1.179):1.16,((t2:1.325)#H8:0::0.6):1.34);")

julia> adjlist = net_to_edgelist(net)
10×4 DataFrame
 Row │ source   target   type    weight   
     │ String?  String?  String  Float64? 
─────┼────────────────────────────────────
   1 │ i3       t4       edge       1.504
   2 │ i1       t3       edge       1.268
   3 │ i1       t1       edge       1.268
   4 │ i2       i1       edge       1.057
   5 │ i2       H8       admix      0.4
   6 │ i3       i2       edge       1.179
   7 │ i5       i3       edge       1.16
   8 │ H8       t2       edge       1.325
   9 │ i4       H8       admix      0.6
  10 │ i5       i4       edge       1.34

julia> using RCall

julia> R"library"("admixtools")

julia> R"plot_graph"(adjlist)

julia> using PhyloPlots

julia> plot(net, showedgelength=true, showgamma=true); # same graph
```
"""
function net_to_edgelist(net::HybridNetwork)
    ne = length(net.edge)
    adj_list = DataFrame(
        source = Vector{Union{Missing, String}}(missing, ne),
        target = Vector{Union{Missing, String}}(missing, ne),
        type   = ["edge" for _ in 1:ne],
        weight = Vector{Union{Missing, Float64}}(missing, ne),
    )
    PhyloNetworks.resetedgenumbers!(net)
    PhyloNetworks.resetnodenumbers!(net)
    nameinternalnodes!(net)
    for (i,e) in enumerate(net.edge)
        adj_list.source[i] = getparent(e).name
        adj_list.target[i] = getchild(e).name
        if e.hybrid
            adj_list.type[i] = "admix"
            if e.gamma != -1.0
                adj_list.weight[i] = e.gamma
            end
            e.length == -1 || e.length == 0 ||
                @warn "break hybrid edge number $(e.number): its length is not going to be recorded."
        else
            if e.length != -1.0
                adj_list.weight[i] = e.length
            end
        end
    end
    return adj_list
end

net_to_edgelist(netpath::AbstractString) = net_to_edgelist(readnewick(netpath))

"""
    edgelist_to_net(dataframe; kwargs...)
    edgelist_to_net(
        source::AbstractVector,
        target::AbstractVector,
        type::AbstractVector,
        weight;
        weighttol = -1e-8
    )

HybridNetwork object built from a list of edges, where each edge corresponds to a row.
If the input vectors are given in a data frame, they should be in columns named
`source`, `target`, `type` and `weight`.
- `source`: name of the parent node
- `target`: name of the child node
- `type`: "edge" for tree edges, "admix" for hybrid edges
- `weight`: edge length for tree edges and the inheritance probability for
  hybrid edges. If not given, then all of them are given the default value of 0.5.

optional argument: tolerance `weighttol`. Any weight `w` that is negative
is changed to 0,
- without any message if `w ≥ weighttol`
- with an error message if `w < weighttol`.

# example

```@repl
julia> net0 = readnewick("((t4:1.504,((t3:1.268,t1:1.268):1.057,#H8:0.0::0.4):1.179):1.16,((t2:1.325)#H8:0::0.6):1.34);")

julia> al = net_to_edgelist(net0);

julia> etype = map(x -> (startswith(x, "H") ? "admix" : "edge"), al[!,:target])

julia> all(etype .== al[!,:type]) # sanity check

julia> net1 = edgelist_to_net(al) # original γ's and edge lengths are used
HybridNetwork, Semidirected Network
10 edges
10 nodes: 4 tips, 1 hybrid nodes, 5 internal tree nodes.
tip labels: t4, t3, t1, t2
((t4:1.504,((t3:1.268,t1:1.268)i1:1.057,#H1:0.0::0.4)i2:1.179)i3:1.16,((t2:1.325)#H1:0.0::0.6)i4:1.34)i5;

julia> hardwiredclusterdistance(net0, net1, true) # true: as rooted networks
0

julia> net2 = edgelist_to_net(al.source, al.target, etype); # weights not given

julia> writenewick(net2) # edge weights default to 0.5 -- γs and lengths
"((t4:0.5,((t3:0.5,t1:0.5)i1:0.5,(t2:0.5)#H1:0.0::0.5)i2:0.5)i3:0.5,(#H1:0.0::0.5)i4:0.5)i5;"

julia> hardwiredclusterdistance(net0, net2, true)
0
```
"""
function edgelist_to_net(df; kwargs...)
    wgt = ("weight" in names(df) ? df.weight : missing)
    return edgelist_to_net(df.source, df.target, df.type, wgt; kwargs...)
end
function edgelist_to_net(
    esource::AbstractVector,
    etarget::AbstractVector,
    etype::AbstractVector,
    eweight::Union{AbstractVector,Missing}=missing,
    weighttol = -1e-8
)
    ne = length(esource)
    length(etarget) == ne || error("source & target vectors should have same length")
    length(etype) == ne || error("edge type vector should have length $ne")
    if ismissing(eweight)
        eweight = [0.5 for _ in 1:ne] # default edge lengths and γ = 0.5
    else
        length(eweight) == ne || error("edge weight vector should have length $ne")
        # check that edge weights are non-negative, with some tolerance
        for i in eachindex(eweight)
            eweight[i] < 0.0 || continue
            eweight[i] >= weighttol ||
                @error("negative edge weights $(eweight[i]) < tolerated $weighttol, for edge # i=$i")
            eweight[i] = 0.0
        end
    end
    nodenames = union(esource, etarget)
    # create all nodes
    nodevec = PhyloNetworks.Node[]
    nodename2index = Dict{eltype(nodenames),Int}()
    isrootnode = trues(length(nodenames)) # to find root later
    for (ni,nname) in enumerate(nodenames)
        nn = PhyloNetworks.Node(ni, true) # leaf by default
        nn.name = nname
        push!(nodevec, nn)
        push!(nodename2index, nname => ni)
    end
    # create all edges, and connects nodes & edges
    edgevec = PhyloNetworks.Edge[]
    ei = 0
    for (snn, tnn, et, ew) in zip(esource, etarget, etype, eweight)
        ei += 1
        sn = nodevec[nodename2index[snn]]
        tn = nodevec[nodename2index[tnn]]
        hyb = et == "admix"
        if hyb
            ee = PhyloNetworks.Edge(ei, 0.0)
            ee.hybrid = true
            ee.gamma = ew
            tn.hybrid = true
        else
            ee = PhyloNetworks.Edge(ei, ew)
        end
        sn.leaf = false # the source is not a leaf
        isrootnode[nodename2index[tnn]] = false # the target is not the root
        push!(ee.node, tn) # child first, because ischild1 is true by default
        push!(ee.node, sn)
        push!(edgevec, ee)
        push!(sn.edge, ee)
        push!(tn.edge, ee)
    end
    # rename hybrid nodes to Hi with i = i0,...,i0+h
    ihyb = 0
    for nn in nodevec # handle rare case when a tree node is already named "Hi"
        nn.hybrid && continue
        m = match(r"^H(\d+)$", nn.name)
        isnothing(m) && continue
        ihyb = max(ihyb, parse(Int, m.captures[1]))
    end
    for nn in nodevec
        nn.hybrid || continue
        ihyb += 1
        nn.name = "H$ihyb"
    end
    # for each hybrid node: normalize its parents' γs, find major parent
    for nn in nodevec
        nn.hybrid || continue
        hybparents = [ee for ee in edgevec if getchild(ee).name == nn.name]
        np = length(hybparents)
        np > 1 || error("hybrid node $nn with $np parents: $hybparents")
        gammasum = sum(ee.gamma for ee in hybparents)
        for ee in hybparents
            ee.gamma /= gammasum
            ee.ismajor = false
        end
        majorparentindex = argmax(ee.gamma for ee in hybparents)
        hybparents[majorparentindex].ismajor = true
    end
    # find the root
    sum(isrootnode) == 1 || error("found $(sum(isrootnode)) root(s)")
    rootindex = findfirst(isrootnode)
    # create the network
    net = HybridNetwork(nodevec, edgevec)
    net.rooti = rootindex
    return net
end
