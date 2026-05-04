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

# How to run this simulation

### Installing software

- Create a new directory where you want to run the simulation,
  then run git clone with this repo

- This simulation needs four executable programs to run:
  [SimPhy](https://github.com/adamallo/SimPhy),
  [Seq-Gen](https://github.com/rambaut/Seq-Gen),
  [RAxML](https://github.com/stamatak/standard-RAxML) and
  [ASTRAL](https://github.com/smirarab/ASTRAL).  
  Download these 4 executable programs and move their respective executable file to a folder `executables` to be created within `simulation-reptiles`,
  with the following titles:

            * SimPhy
            * seq-gen
            * raxmlHPC-PTHREADS
            * astral.5.7.8.jar

- Now the simulation is ready to be run. Navigate to the
  `simulation-reptiles/scripts` directory.
  There are two julia scripts for this simulation.
  The first is `julia clean.jl`: this will remove all old output files and leave the simulation ready to be ran.
  The second is `julia simulation.jl`: this will run the simulation based
  on the parameters passed to the simphy configuration file in the
  `simulation-reptiles/simphy-configs` folder.

### Installing executables

for SimPhy, on `franklin00` server:

```shell
cd ~/private/apps # or other private folder for private installs
wget https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz
tar -xvf SimPhy_1.0.2.tar.gz
cp SimPhy_1.0.2/bin/simphy_lnx64 ~/bin
chmod ug+x ~/bin/simphy_lnx64
cd - # go back to simulation-reptiles folder
ln -s ~/bin/simphy_lnx64 executables/SimPhy
executables/SimPhy -h # just checking
```

for RAxML, still on `franklin00`:

```shell
git clone https://github.com/stamatak/standard-RAxML.git ~/private/apps/standard-RAxML
cd ~/private/apps/standard-RAxML # or other private folder for private installs
make -f Makefile.PTHREADS.gcc; rm *.o
mv raxmlHPC-PTHREADS ~/bin/
cd - # go back to simulation-reptiles folder
ln -s ~/bin/raxmlHPC-PTHREADS executables/raxmlHPC-PTHREADS
executables/raxmlHPC-PTHREADS -h # just checking
```

for ASTRAL:

```shell
cd ~/private/apps
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip
cd - # go back to simulation-reptiles folder
ln -s ~/private/apps/Astral/astral.5.7.8.jar executables/astral.5.7.8.jar
```

# Simulation plan

### get a species tree from real data

species tree topology, with branch lengths in coalescent units,
and also in substitution (substitution/unit of time) units.
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
tree = readnewick(treestring)
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

For #1 and #2: use the notation from the newick format used by [SimPhy](https://github.com/adamallo/SimPhy/wiki/Manual#521-input-files-newick-tree-format):
`: branchlength_num_generations * substitution_rate_multiplier ~ generation_time_multiplier # Ne`
although the spaces are just for readability --there shouldn't be any spaces.
We get this (same as at the end of `scripts/speciestree.jl`):

    (Homo:3.44*0.0100947,((((Crocodylus:0.88*0.0042057,Alligator:0.88*0.0036776):1.71*0.0078509,(Taeniopygia:0.93*0.0235933,Gallus:0.93*0.0199793):1.66*0.0079913):0.17*0.0068836,Chrysemys:2.76*0.0067212):0.18*0.0098089,(Anolis:0.5*0.0797969,Pantherophis:0.5*0.1796924):2.44*0.0190487):0.5*0.0694588);

For #3: the best fit the to crawford's genes rate was lognormal, so use
the `HL` option for a log-normal distribution of gene-family-specific rates
(h for heterogeneity, l for locus): `-hl l:-0.19,0.6164414002968976`
(or `ln`?) which has mean 1 because `0.6164414002968976 = sqrt(2*0.19)`.

For #4: look at gene trees some more (to do).
For now, use a Gamma distribution with mean 1 and shape 10
(SimPhy option `-hg f:10`), for little variation.

Parameters for rate variation across sites:
we looked at log files from IQTree (to estimate gene trees),
to get the substitution parameters in real gene trees:
see [choice-seqgen-parameters.md](choice-seqgen-parameters.md).

Conclusions:
- to simulate all genes with the same substitution model, use:
  * HKY (-m option) with transition/transversion ratio kappa = 4.143 (option -t)
  * base frequencies 0.316,0.182,0.183,0.319 (-f option)
  * shape alpha = 0.356 (-a option) for the Gamma distribution of rates across sites

- to simulate each gene with its own substitution model, use HKY with:
  * kappa from LogNormal(μ=1.4215, σ=0.2798)
  * frequencies from Dirichlet(66.59, 38.41, 38.61, 67.12)
  * alpha from Gamma(α=3.267, θ=0.109).


### simulate an alignment along each gene tree

using [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/)
we'll need to size of each gene: in # of base pairs

### analyze the simulated sets of alignments

using the same methods we used for real data (or a subset?)
* to estimate gene trees: IQ-Tree, RAxML, FastTree, MrBayes
  - focus on RAxML only
* to estimate a species tree from a set of gene trees:
  ASTRAL, and network methods under the constraint of h=0 hybrid nodes
* to estimate a network from a set of gene trees: SNaQ and Phylonet_MPL,
  under the constraint h≤1.

### summarize

- results of ASTRAL: how often is the true species tree recovered?
- results of network methods (e.g. SNaQ) with h=0: idem.
- results of network methods with h=1:
  what is the support for a reticulation? gammas?

## License

This repository is primarily licensed under the MIT License (see `LICENSE`).

**Third-Party Code:**  
The script `third_party/vcf2eigenstrat.py` is derived from a project by Iain Mathieson and is licensed under the Apache License 2.0.  
See `third_party/LICENSE-APACHE-2.0.txt` for details.


