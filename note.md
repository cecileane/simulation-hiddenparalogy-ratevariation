##Notes

Todo This Week:

  *Use tree number two to run simulations. (crawford astral iq line 44 of git spreadsheet)
  *remove top group (turtle if need be) (Pelomedusa)
  *play with external branch lengths to make it ultrametric
  *run simphy with the "fat" tree
  - write script to run simulation of 1000, 100 times
  *pass simphy to seqgen to get gene tree data
  - write script to run this section 100 times
  *pass seq gen data to our methods to obtain understanding and complete analysis

  - try for one script that can get it done in one go using picture as reference. (100 times)
  

Seq-Gen params

-n 1000 -l 1000 -m GENERAL <tree from loop of SimPhy data>
 do we need weights? or is this all that we need

-n 1000 - l 1000 -m HKY -a 0.356 -t 4.143 -f waiting for -z 22 -on

will need to scale branch lengths 

median of distributions
transition/transversion ratio kappa=4.143 (-t in seq-gen)
alpha=0.356 for rate variation across sites (-a option in seq-gen)



SimPhy Params

-rs 2 //Number of replicates
-rl f:10 //1000 locus per replicate - Start here - work up to 1000 when we know that our settings are correct
-s (Homo:3.44,((((Crocodylus:0.88,alligator_mississippiensis:0.88):1.71,(Taeniopygia:0.93,Gallus:0.93):1.66):0.17,Chrysemys:2.76):0.18,(Anolis:0.5,Pantherophis:0.5):2.44):0.5);
-sg f:0.00001 //Generation time - look back at this
-sp f:10000 //Population size?
-gt f:0 //Genome-wide horizontal gene transfer (sampled once for each species tree, and applied for all locus trees) we want no horizontal gene transfer look into what value should go here
-hg f:2 //Heterogeneity sampled for every gene tree branch using the same Gamma distribution with shape = rate =100 - can adjust later
-cs 22 //Seed for the random number generator, in order to make the experiment repetible.
-om 1 //Tree mapping output
-od 1 //Database
-op 1 //Output with the general sampled options (describes the simulation run)
-oc 1 //Activates the backup of the original command line and configuration file (we recommend to always activate this option)

Rate lineage in SimPhy
see if can find way to fix the rate
add SU l:mean,sd






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


My start

  -rs 2 //Number of replicates
  -rl f:5 //1000 locus per replicate - Start here - work up to 1000 when we know that our settings are correct
  -lb f:0 - duplication rate - 0 duplications
  -s (outgroup:5.0,((((crocodilia:1.84,testudines:1.844):0.182,bird:2.778):0.405,squamata:3.878):1.442)); //Fixed species tree NEED TO DOUBLE CHECK THIS LINE
  -sg f:0.00001 //Generation time - look back at this
  -sp f:10000 //Population size?
  -gt f:0 //Genome-wide horizontal gene transfer (sampled once for each species tree, and applied for all locus trees) we want no horizontal gene transfer look into what value should go here
  -hg f:2 //Heterogeneity sampled for every gene tree branch using the same Gamma distribution with shape = rate =100 - can adjust later
  //Standard Below Here
  -cs 22 //Seed for the random number generator, in order to make the experiment repetible.
  -om 1 //Tree mapping output
  -od 1 //Database
  -op 1 //Output with the general sampled options (describes the simulation run)
  -oc 1 //Activates the backup of the original command line and configuration file (we recommend to always activate this option)


  -SU  average rate of substitution per coalescent unit - defaults to f:0.00001 - check with group

5 coalescent units = # generations / population size
t = years; -sg = generation time in years / generations
number generations = t / sg
5 coalescent units = (t / sg) / sp

from this we can use coalescent units as the tree branch lengths

With this I also believe we can just use the Loci Tree as out gene data can this be confirmed?


Ultrametric tree incase there is an error down the road

(outgroup:5.0,(((crocodilia:1.84,testudines:1.84):0.18,bird:2.02):0.40,squamata:2.42):2.58);




Definitions for quick reference:

  Speciation: separation of one ancestral population into two new populations that do not interbreed.

  Extinction: disappearance of a species.

  Gene duplication: copy of one gene into a new locus in one individual of the population, which gets fixed in the sample (i.e., we assume no duplication polymorphism).

  Gene loss: deletion of one gene in one individual of the population, which gets fixed in the sample (i.e., we assume no loss polymorphism).

  Horizontal gene transfer: copy of one locus from one species to another contemporary species via replacement. The transfer initially affects one individual in the receptor species, and gets fixed in the sample (i.e., there is no transfer polymorphism).

  Gene conversion: replacement of one homolog by another within a single species. This conversion initially affects one individual and then gets fixed in the sample (i.e., there is no gene conversion polymorphism).

  Lineage sorting: consideration of the coalescent process of the sampled gene copies, allowing their history to be incompatible with the species tree history.