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

- get a species tree from real data
  * topology
  * branch lengths in coalescent units, also in substitution units
    for this we need to process the gene trees estimated from the real data. 

example to play with:
```julia
using PhyloNetworks
using PhyloPlots
tree = readTopology("(outgroup:5.0,((((crocodilia:1.84,testudines:1.844):0.182,bird:2.778):0.405,squamata:3.878):1.442));")
plot(tree, :R, showEdgeLength=true, useEdgeLength=true);
```

- simulate gene trees along this species tree, using the coalescent, using
  [SimPhy](https://github.com/adamallo/SimPhy) and
  [paper](https://dx.doi.org/10.1093%2Fsysbio%2Fsyv082)
  * need the rate of evolution of each gene tree:
    both average rate, also variation of rates across genes (get this from
    the real data)

- simulate an alignment along each gene tree using
  [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/)
  we'll need to size of each gene: in # of base pairs

- analyze the set of alignments using the same methods we used for real data
  (or a subset?)
  * to estimate gene trees: IQ-Tree, RAxML, FastTree, MrBayes
  * to estimate a network from a set of gene trees: SNaQ and Phylonet_MPL.

- summarize: what is the support for a reticulation? gammas?

