##Notes

Todo This Week:
  Find substitution lengths of the species tree
  Begin running SimPhy simulations based on the data gathered from the original tree

#SimPhy Notes

  *Prioritizing species tree?
  *Go through SimPhy notation and verify that I know what I am doing haha.
  *Assuming I want to use -S to run with data based on our species tree.
  *ST for the height of the tree

Here is the example given by SimPhy so I assume that we will start with this and build off of it

-rs 100 //Number of replicates
-rl f:100 //100 locus per replicate
-s (((A:10000,B:10000):30000,C:40000):1000,D:41000); //Fixed species tree
-sg f:0.5 //Generation time
-sp f:10000 //Population size
-gt e:10000 //Genome-wide horizontal gene transfer (sampled once for each species tree, and applied for all locus trees)
-hg f:100 //Heterogeneity sampled for every gene tree branch using the same Gamma distribution with shape = rate =100
-cs 22 //Seed for the random number generator, in order to make the experiment repetible.
-om 1 //Tree mapping output
-od 1 //Database
-op 1 //Output with the general sampled options (describes the simulation run)
-oc 1 //Activates the backup of the original command line and configuration file (we recommend to always activate this option)

With this I also believe we can just use the Loci Tree as out gene data can this be confirmed?


Definitions for quick reference:

  Speciation: separation of one ancestral population into two new populations that do not interbreed.

  Extinction: disappearance of a species.

  Gene duplication: copy of one gene into a new locus in one individual of the population, which gets fixed in the sample (i.e., we assume no duplication polymorphism).

  Gene loss: deletion of one gene in one individual of the population, which gets fixed in the sample (i.e., we assume no loss polymorphism).

  Horizontal gene transfer: copy of one locus from one species to another contemporary species via replacement. The transfer initially affects one individual in the receptor species, and gets fixed in the sample (i.e., there is no transfer polymorphism).

  Gene conversion: replacement of one homolog by another within a single species. This conversion initially affects one individual and then gets fixed in the sample (i.e., there is no gene conversion polymorphism).

  Lineage sorting: consideration of the coalescent process of the sampled gene copies, allowing their history to be incompatible with the species tree history.