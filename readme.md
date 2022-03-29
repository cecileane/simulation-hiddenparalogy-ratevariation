# goal

Evaluate if some model violations, or lack of information in gene trees
may result of overestimation of reticulation by network reconstruction methods.

Start with a simulation without any model violation: under a tree
estimated from the real data, and parameters estimated from real data.
See if the network methods display the behavior as in the real data
(with 1 reticulation inferred with hybridization proportions close to .5/.5)

parameters from real data: very low information content within a single gene,
that is, very low rate of evolution.

Later, if we have time, add a violation of assumptions
(hidden paralogy: gene duplications & losses? long branch attraction?)

## simulation plan

### get a species tree from real data

topology, and
branch lengths in coalescent units, also in substitution (substitution/unit of time) units
for this we need to process the gene trees estimated from the real data.

use the tree from the **crawford** data (10 taxa) obtained with ASTRAL + IQTree
[here](https://github.com/cecileane/reptiles/blob/main/estimatednets_collapsed.csv#L44)
because it has the most accepted topology, and
shorter internal edge lengths (than the chiari tree)
for more incomplete lineage sorting. But:
subsample taxa to match the taxon sampling from shaffer data (8 taxa),
to get faster analysis time.
**Remove**:
- sphenodon
- one turtle (Pelomedusa) if need be
This was done in file `scripts/speciestree.jl`. The final tree is
copy-pasted below as a string, with code to read it in julia.
Taxon names could be simplified for the simulation!

```julia
using PhyloNetworks
using PhyloPlots
treestring = "(Homo:3.44,((((Crocodylus:0.88,alligator_mississippiensis:0.88)100.0:1.71,(Taeniopygia:0.93,Gallus:0.93)93.2:1.66)0.0:0.17,Chrysemys:2.76)0.0:0.18,(Anolis:0.5,Pantherophis:0.5):2.44)0.0:0.5);"
tree = readTopology(treestring)
plot(tree, :R, showEdgeLength=true, useEdgeLength=true);
```

### simulate gene trees along this species tree

1000 gene trees in each data set to match the size of our real data
(chiari: 248 genes, other data: between 1113 and 1955 genes).
using the coalescent, using
[SimPhy](https://github.com/adamallo/SimPhy) and
[paper](https://dx.doi.org/10.1093%2Fsysbio%2Fsyv082)

To do this, we need the rate of evolution of each gene tree:
average rate, also variation of rates across genes.
We should get this from the real data.

Parameters to get gene tree branch lengths in substitutions per site:
1. an overall genome-wide substitution rate: in substitutions/site per coalescent unit
2. a species-specific rate to model rate variation between species:
   for example if birds evolve faster, or if turtles evolve slower,
   or if some ancestral lineage evolved faster or slower
3. a distribution of rate variation across genes: if some genes evolve faster or slower than others
4. a distribution of gene x lineage rate variation:
   if lineages evolve at different rates in genes trees,
   in a way that's independent across genes and lineages.

Parameters for rate variation across sites:
look at log files from gene trees, to get the substitution parameters
in real gene trees.

to do:
- look at the branch lengths of all IQ-Tree gene trees
- average or median tree length in substitutions / site (#1)
- distribution of tree length across gene trees (#3)
- among genes that match the species tree:
  average or median of their length for each branch (#2)

### simulate an alignment along each gene tree

using [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/)
we'll need to size of each gene: in # of base pairs

### analyze the simulated sets of alignments

using the same methods we used for real data (or a subset?)
* to estimate gene trees: IQ-Tree, RAxML, FastTree, MrBayes
* to estimate a network from a set of gene trees: SNaQ and Phylonet_MPL.

### summarize

what is the support for a reticulation? gammas?

