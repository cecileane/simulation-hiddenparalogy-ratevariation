
# Jan 23 to Jan 30 
1) In simulation_iqtree.jl, after max_iter, if lower_threshold (change to a more informative name, min_gene_porportion) not met, ignore the replicate intead of stopping the pipeline. 
2) Re-organize the output folder structure so that the top layer is each rep. Inside each rep, there are outputs from genetrees_simphy, iqtree, astral, seq-gen, genetrees_singlecopy folders.  
3) Revise documentation: After running max_iter, if still not hitting min, then ignore the rep

# Jan 16 to Jan 22
To do list: 
1) Get ready for server -- Not heard from stat lab. 
  Note here: Botany server doesn't have enough memory to run iqtree parallelly. In iqtree.pl, num core needs to be set at 2 for now. 
2) git rm all deleted files -- Done (Agenda for this week's meeting: merge all commits and delete useless branches)
3) Codes for re-run simphy untill it reaches to a lower_threshold -- done. The codes have been tested for different ratevar, duploss, n_genes, n_reps and it all runs pretty fast. Running n_genes = 100 and n_reps = 100 takes about 1 min. (Agenda for this week's meeting: Having many different functions for this part. Potentially having another scripts or move them to utilities.jl?) n_genes = 100 and n_reps = 1000 took about 10 to 15 mins. 

4) Shorten dog strings -- Done. 

Some trials for testing the codes for point 3 and time to run it: 

1. julia scripts/simulation_iqtree.jl --dup_rate 0.01 --loss_rate 0.01 --ratevar G --n_genes 100 --n_reps 10 --seed_simphy 12345 -- Took a few seconds, much less than one min 
2. julia scripts/simulation_iqtree.jl --dup_rate 0.01 --loss_rate 0.01 --ratevar G --n_genes 25 --n_reps 25 --seed_simphy 54321 -- Took a few second, much less than 1 min 
3. julia scripts/simulation_iqtree.jl --dup_rate 0.01 --loss_rate 0.01 --ratevar G --n_genes 25 --n_reps 25 --seed_simphy 123 -- took a few second, much less than 1 min 
4. Try a lower number of max_iteration: julia scripts/simulation_iqtree.jl --dup_rate 0.01 --loss_rate 0.01 --ratevar G --n_genes 25 --n_reps 25 --seed_simphy 1000 --max_iteration 2 -- still much less than a min 
5. julia scripts/simulation_iqtree.jl --dup_rate 0.02 --loss_rate 0.01 --ratevar GL --seed_simphy 4321 --n_genes 20 --n_reps 100 -> This runs for about one min with n_reps = 100 
6. julia scripts/simulation_iqtree.jl --dup_rate 0.02 --loss_rate 0.01 --ratevar L --seed_simphy 4321 --n_genes 100 --n_reps 100  -> with n_reps = 100 and n_genes = 100, this runs for about 1 min. 
7. julia scripts/simulation_iqtree.jl --dup_rate 0.02 --loss_rate 0.01 --ratevar N --seed_simphy 4321 --n_genes 100 --n_reps 10 -> took a few seconds 
8. julia scripts/simulation_iqtree.jl --dup_rate 0.03 --loss_rate 0.03 --ratevar G --seed_simphy 12345 --n_genes 30 --n_reps 30 -> took a few seconds 
9. julia scripts/simulation_iqtree.jl --dup_rate 0.03 --loss_rate 0.03 --ratevar GL --seed_simphy 12345 --n_genes 30 --n_reps 30 -> took a few seconds  


Notes from meeting: 
1) SLURM (used in the industry): There is a tutorial on UW statistics (high performance). There are other tutorials. Website (States 679: Job scheduling)
2) Darwin -- need to check job scheduling (check hostname command). Private and Public. Put stuff in private (not public directory). All works go to the private directory. AFS (shared spot -- instruction on the CS website) -- don't store large output files in afs. 
3) run and store the large files in the harddrive of each machine. (nobackup). Clone the github repo each time in nobackup and run the codes to get the data 
In nobackup for each machine, create a folder with my login name and only do stuff inside 

ls -l (The first line shows the permission for owener, group, everyone else (r - read, x - exe, w - write)) 

4) Not for lunchbox, but for Darwin and Franklin, no permission after log-off. If I want to continue running the script after log-off, need to do stuff (check email sent by Cecile on Jan 16, 2025) 

5) Learn to have useful commit message. Don't create branches. 
Squash all smaller commits to a bigger, cleaner commit -- /delete the older branch after finishing minor jobs, while just leaving the bigger commit for a cleaner history 

git log --abbrev-commit --pretty=oneline --all --decorate 

6) git rm all the deleted files (next time do the squashing thing) --done. 

7) Shorten the dog-string 


# Dec 13 to Dec 18 
1) Think of changing iqtree.pl to something that can use the directory as the input in iqtree, which could better utilizes the multi-thread process and be faster. -- Finished but questions left 

Question: If using more than 2 cores, the codes will end with error "IQ-TREE CRASHES WITH SIGNAL SEGMENTATION FAULT" indicating that insufficient memory problem. If -T is 2, then the issue is solved. For now, the code is still running with using individual file as the inputs. 

2) First, need something to track the number of missing simulated trees. Then, get an estimated % of how many tree got missing with N_genes. This determines how many gene trees more need to be simulated. (before seq-gen.sh). 

Have a target N genes to hit and run a loop 
No infinit loop (max iteration = M / if the success percentage is too low, then exit out and give a warning -- <1% )
Send a info.log file to check 
Calculate the number of genes to be used 
Have a lower threshold (check the written notes) and then set a max iteration. After max iteration if it is higher than lower threshold while still lower than the target, then exit.  
Test with target 10 genes for test run (> 2)  
-- Have a test file  



3) In utilities.jl -- revise documentation: Change to """ """ just before the function (julia documentation. No line in between the dog strin and the actual function) 
-- Can also have a test file. -- Done but need to double check.  

4) Remove intermediate files created by iqtree. Need to double check but done. 

5) Alison (based on iqtree) simulator to simulate molecular sequences. (No need but interesting to check). -- Not necessary.  


# Dec 6 to Dec 12, 2024 
1) Potentially change to iqtree. Create iqtree.pl file based on raxml.pl -- Done. 

Notes from the meeting on Dec 6 
ramxl.pl line 195: remove one identical sequence and add the taxa to the same place of the tree. This is why ".reduced" got removed. 
Learn how to use perl. 
In iqtree, no boostrap but run pp on astral. -- Done 
Some notes during meeting on Dec 5: The reason we want to change to iqtree is because raxml.pl calculated BS, which is very slow. For iqtree, we want to change to no BS and later in the astral, use PP as the support. 
For record, I am running raxml.pl fir n_genes = 2, n_reps = 2 and n_inds = 2 for four hours until now on the Botany server, while it is still not finished! 

Notes when changing from raxml.pl to iqtree.pl: 
Some notes when I worked on changing to iqtree.pl -- This indeed works sooo fast! 
a. Two simulate julia files: simulation_raxml.jl uses raxml.jl (the original version) and simulation_iqtree.jl uses iqtree.jl (new version with iqtree.jl). 
b. I tested the code with dup_rate = 0 and 0.01, loss_rate = 0 and 0.01, n_reps = 2, n_genes = 2, n_inds = 2, ratevar = G, N, GXL, and seed_simphy = 12345. It seems that 
simulation_iqtree.jl and iqtree.pl are robust in creating iqtree gene trees and astral species tree. The script can handle missing gene trees from modify_newick from utilities.jl 
c. Running the test kit using simulation_iqtree.jl takes a few seconds (certainly less than 1 min) so this is much faster than running BS in RAxML. 

Something I need to do: Modify the iqtree.pl script so it can remove the unnecessary intermediate files. 
What are some intermediate files that could be deleted? 
iqtree.pl.log: This documents the command line run for each gene. Better to save this. 
geneX.treefile: The major treefile. This is the individual tree file for each locus in each rep. This will be concatenated into one BestTrees.tre file to be input into astral. However, for now, we can save this first. If space becomes a major issue, we can delete this. 
geneX.mldist: pairwise distance file. Delete this? 
geneX.log: Log for running iqtree on this gene. Delete this? 
geneX.bionj: Neighbor-joinning starting tree. Delete this. 
geneX.ckp.gz: compressed checkpoint files.  -- code finishd for this but for now got commented out. 


2) Change "taxa" to "gene copy" in the description of utilities.jl. -- Done. 

3) Remove sim_hidden argument from simulation.jl. When sim_hidden = F, this implicitly means that Dup_rate = 0 and loss_rate = 0, so no need to have this argument. -- Done. 

Fix this: When --ratevar N or --ratevar GxL, the new utilities/modify_newick function doesn't work. Check for details. When using N, the error comes from line 191, when using GxL, the error comes from line 119, indicating that the simulation.julia code cannot handle missing tree file. -- Fixed! 

4) Make detailed notes of NetRAX and drop it (try to document it and then drop it). (Need to do)  

# Nov 15 to Dec 5, 2024 
1) Have examples for each function in utilities or scripts. For example, have some strings listed in utilities.jl to run each function. Just use the test set to write the function as the test set. Done.  
2) License issue -- change to the concatenate_seq.py script -- check the file path. -- done. This current script also handle the missing taxa. For example, if one nexus file is missing one taxon, the script will use "-" to fill in the gap. Question: How to put the test set? 
3) Add different accessions to taxon -- simphy. Done. I also changed the utilities (revise newick strings function set) and concatenation code so that the current simulation code is able to handle multiple taxa. This process needs to be re-thought. BTW, after adding additional taxa, even if I set N_inds = 2, the codes are running very slow. 

Questions to double check: 
If sim_hidden = true, then simulation.jl will start to modify the newick strings so that: 1. Any tree will repeated gene copy will not be processed (more than one gene copy per indicidual). 2. Any tree with more than 4 missing taxa will not be processed (not enough left taxa). 3. The tree tips will be changed from A_(locus-id)-(accession-id) to A_(accession-id). 

If sim_hidden = false, then no modification will be made to the newick trees. Here, this should be changed that even if sim_hidden = false, tree tips will still be modified so that A_(locus-id)-(accession-id) to A_(accession-id). However, the function should not exclude trees with repeated taxa or trees with more than 4 missing taxa.  

4) I got the responses from NetRax. It seems to be a compiler issues from gcc. My old system has gcc 4+ and then I used gcc 13+. NetRax developer suggested me to downgrade my gcc version 11 as suggested. I re-installed the dependencies isl version 0.24 (wget https://libisl.sourceforge.io/isl-0.24.tar.xz). Still working on it. 
Installation is done in conda environment NetRax. 

5) Change the output file from simulation.jl. Add sim-hidden = true or false? -- done.

Not finished: 
6) Check NANUQ+ paper and potentially incorporate it: https://johnarhodes.r-universe.dev/articles/MSCquartets/NANUQplus.html Not done. 
and https://www.biorxiv.org/content/10.1101/2024.10.30.621146v1  check TINNIK and each step.  
7) For find_graph, ask Lauren -- possibly talk with her about contribution or authorship. Scheduled a meeting with Lauren and discussed how to process the file. The first molecular sequence could be used as the reference for SNP calling. 


# Nov 7 to Nov 14, 2024:
1) Next Thursday meeting online at 3pm. Send Cecile the zoom link before the meeting.
Changed to Friday at 4 online. done. 

2) Ask developer how to install NetRax in Github issues. --- has no response yet. Based on the average response speed on NetRax github, they might take about 10 to 15 days to respond, so that's fine. Let's wait. 

3) Use find_graph. <--- Reminder: the question I want to talk for today! 
Find_graph needs to calculate f-stats. I have a question here regarding calling SNPs. Now we have concatenated dna sequences for each replicate. To calculate f-stats from this file, the basic workflow is to extract snp first and then convert the data to a genotype matrix (Plink format or *EIGENSTRAT*). Using admixtool, we can simply get the f-stats from the Plink. 

Proposed workflow and script to write: 
 -- use vcftools to compare multiple alignments and find the sites that are different across taxon (extract snps without reference fasta). In bcftools, "--no-reference" can allow for no reference. Then, similar to reference-based method, use "bcftools mpileup" and "bcftools call -mv" to call snp. 
 -- convert snp vcf to plink by using: vcftools --vcf mydata.vcf --plink --out genotype_data 
 -- This genotype_data should be used for admixtools in R. 
I want to double check if this workflow is correct. 
Or, I should consider how to include multiple accessions to each taxa first since admixtools are designed for it.

4) Change the concatenated_seq.py file into the previously written file. 

Used this script and incorporate into the workflow in simulation.jl https://github.com/tkchafin/scripts/blob/master/concatenateNexus.py -- How to get the license? Next step unsure. This code needs to be modified in the end to change the resulting concate.nex to fasta files to be further used in find_graph or NetRax. -- Change to .fasta

Other notes: https://gist.github.com/cecileane/d8a7e9e4d4b4fadeb5c62e89506cf3e9 This code concatenated .phy files. It needs to have a current folder with all alignment. However, this won't work for our case. In out seq-gen output, there are .phy files only for the last rep. Based on raxml.pl: it seems that RAxML is executed in an iterative process for each gene, and the script may overwrite the input/output files such that only the last replicate's PHYLIP file gets processed. Thus, we don't have the .phy for all rep at once. This could be feasible if I change the ramxl.pl script. However, we have all .nex files so it would be better to start working on the .nex files. 

5) Write the function/script to process intermediate files. --- Done. The actual code is implemented in scripts/utilities.jl. The implementation of the function is in simulation but the actual folder structure needs to be discussed.  

Some notes from simphy tutorial: "Locus tree leaves are automatically labeled using modifications to the species name in order to indicate both the species and locus, following the scheme "species_locusid", while gene tree leaves also indicate the individual “species_locusid_individualid”. Locus trees may also contain lost leaves. In this case the scheme indicates both the evolutionary event that generated them and the corresponding species tree node id “event-streenodeid”. More than one lost leave can pertain to the same species tree node, and consequently some leaves can have the same name. When using the outgroup addition option (see 5.1.1) the outgroup species is called “0”." Therefore, it seems that the second number in the labeling is the locus id which we want to remove. In case if we have multiple individuals, we want to keep the individual ID. 

6) Add gene loss rates into the simulation.jl. --- Done. The output file is named as dup_$dup_rate_loss_$loss_rate_ratevar_$ratevar. $dup_rate and $loss_rate are two parameters and can be specified differently. -- No the loss rate and dup rate are set to have different values. However, loss rate could be linked to dup rate based on simphy tutorial by "-ld f:lb". Question: Should we have loss rate linked to dup rate (loss rate is forced to be the same as dup rate)? The current script now can have separate dup rate and loss rate. 

7) There is a bug in using ratevar = "N" in simulation.jl. When running this code, it will throw an error message. Fix it. -- Done. It turns out to be a sting manipulation issues but all resolved now. 



# Oct 31 to Nov 7, 2024: 
1) Change snaq filenames to be simpler. No need to have all the rate and rep since it is already specified in the folder name -- Done
2) SnaQ files should be listed as two: snaq_onedata or snaq_submission. Remove snaq loop from the end. The submit file is named as snaq_submit.jl.  Question: Think of how to simply run this script for a particular rep? Like rep 1 to 10 or rep 12. The last part in simulation.jl should be changed as well. -- Done. The current version of snaq_submit.jl could specify which replicate to be run by using "--rep_start" and "--rep_end" arguments. 
3) Write NetRax script -- Installing NetRax is hard ... Should I even use another server? :( 
  Write the code to concatenate the resulting nexus files from Seq-Gen to a concatenated fasta file, which is useful for NetRax. However, first the code is in Python since it is easier to write. I tried to write in bash or julia but all failed. I am not sure it is fast. It might be good to change it to a julia script, but I didn't find a julia package to read in nexus file -- should I change to .pl, .jl, or .sh? 

  This concatenation script can handle duplication but I am not sure if this is the best approach to do it. See the example in /home/bli283/simulation-reptiles/output/DL0.01-RVG/seq-gen-outfiles/rep1 and F_0_1  

  Thinking of this further, how to simulate hidden paralogy. In this case, there are two copies of F. One is a copy of F gene. If there is no F_0_0, then F_0_1 and F_0_0 will be hidden paralogy. Questions here!

Problem found when testing simulation.jl: 
when specify ratevar = N, the code doesn't work because it lacks l parameters in simphy. 

Similar to the simulation.jl. The should should be applied to NetRax 
Optional: Look at AliSim as substitute to Seq-Gen -- this is not necessary tho. 

# Oct 7 to Oct 30, 2024: 
To do list for weeek Oct 10: 
1) Read the mannual of astral and document the values on the trees (line 101 or 102), where is uncertainty? Done. 
Note: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#multi-locus-bootstrapping 
Here is the information. The first 100 trees are boostrapped replicate trees. 101st tree is a greedy consensus of the 100 bootstrapped replicate trees; this tree has support values drawn on branches based on the bootstrap replicate trees. Support values show the percentage of bootstrap replicates that contain a branch. The 102nd tree is the “main” ASTRAL tree; this is the results of running ASTRAL on the best_ml input gene trees. This main tree also includes support values, which are again drawn based on the 100 bootstrap replicate trees. The 102nd tree is the tree we should use. In snaq.jl, I hard-coded this tree 102 in the code for now. Potentially I could hard-code to be the last line of the tree.

2) Check astral v.5.7.8 (BS) and astral-pro (PP): 1) compare the outputs from both, 2) check mannual (check uncertainty), 3) which one is faster?  4) From the astral v.5.7.8, can we get both BS and PP? 
Using astral v.5.7.8, it is not possible to get both BS and PP, and the new ASTER can do PP. Based on this paper https://ar5iv.labs.arxiv.org/html/1904.03826 , PP on astral tree is the likelihood that the branch is correct, based on the frequency of quartets in gene trees around that branch. If a branch in the species tree is accurate, each associated quartet tree should occur in the gene trees with a probability of at least 1/3, which is a threshold under which no single tree would dominate randomly. This is relatively new while BS is more tranditional. However, PP is faster than BS since it uses quartets but many traditional papers still use BS. 

Conclusion: PP is faster but BS is more commonly used. Some recent papers have adpated PP though. 

3) Read the SimPhy paper -- Done but haven't yet read the paper about the ILNDEible paper. Question: Should we have HKY model, base frequencies, shape based on the README? See simphy tutorial about INDELible_wrapper.pl I think we should incorporate INDELible_wrapper.pl as a later step after running SimPhy 

4) Use ArgParse. --- Done

5) Write script for SNAQ. --- Done. script/snaq.jl is incorporated in simulation.jl in the end. It could be useful to separate simulation.jl with another script called estimated_networks.jl but that would require to sepecify the path again. Now, simulation.jl has too many argument setting. One of the old versions I pushed to github doesn't contain bootsnap. The newer version committed on oct 30 contains bootsnaq.  

My question: should I have two script with simulation.jl to do simphy, seq-gen, raxml and astral, and estimation.jl to run SnaQ and other network approaches, or the current simulation.jl is enough. I think having one simulation.jl is sufficient since we can keep using a lot of variables specified, but it also means the simulation.jl has many arguments now. 

Another question: For H=1, I included bootnet but I didn't do this for H=0. Should I include bootsrap? 

The current structure of the output files. This needs to be revised after all coding is done: 
Output:
  DL$duploss_RV$ratevar: The folder for each parameter set
    sim-phy-outfiles:
      rep$ID: 
      simphysim-conf-$parameters: configuration file
    seq-gen-outfiles: 
      rep$ID
      phylip: an intermediate folder created by raxml.pl
    raxml-outfiles:
      raxml-outfiles-DL$duloss_RV$ratevar_$repID: folder containing raxml results for each replicate
        raxml.pl.log: Log file from running raxml.pl
        contree.tgz: 
        bootstrap: Folder containing boostrapped replicates for each gene
        besttrees.tre: The gene trees used for SnaQ 
        besttrees.tgz: 
    astral-outfiles:
      astral-outfiles-DL$duploss-RV$ratevar-$repID: Folder containing astral results from each replicate 
       ...
    snaq-outfiles: 
      rep$ID: folder containing snaq results from each replicate
        H0_output: Folder containing Hmax = 0 outputs
        H1_output: Folder containing Hmax = 1 outputs
It might be better to think of potential ways to have more consistent labeling of each folder. 

6) Read the paper and start to think of the algorithum for distance-based comparison between semi-directed trees. --- It is indeed very interesting although I cannot understand all the details. I need a lot of help to get started. 

7) There are many methods to estimate networks. Besides SnaQ, what else should we use? 
My preliminary proposal mentioned two methods mainly, sanq and PhyloNet_MPL, which are both MPL methods. Maybe it would be better to compare between different algorithums. Here I list several methods that are different from SNAQ and could be used. To summarize, NetRAX doesn't account for ILS but it is fast. NANUQ doesn't require pre-specified number of reticulation and has a hypothesis testing process to test if the tree is a non-tree structure. Poolfast or find-graph is f-statistics based method, and based on previous studies, lineage variation could cause false detection of hybridization. PhyNEST (or BEAST2 but PhyNEST is a julia package already) can directly infer networks from sequences, which is really different from all the other methods. 
  a. NetrRAX: ML inference of phylogenetic networks that doesn't account for ILS.
  b. NANUQ: Using multispecies coalescent model but doesn't require a pre-specified number of reticulation events unlike snaq. It is still a quartet-based approach, but NANUQ uses hypothesis testing on the quartets to confirm the presence of network cycles. These statistical tests check for specific patterns that suggest non-tree-like evolution. The algorithm incorporates NeighborNet and Circular Network methods to create splits graphs, visual representations of the species network, which are then interpreted to infer network topologies.
  c. poolfast or find_graph (admixture tools): R package which uses f-statistics and admixture graph. Similar to admixture tools. I think it would be important to include a method like this which uses f-statistics, because based on previous research lineage variation could impact the detection of hybiridization using test statistics. Here, using such tools could help us verify this. I think picking one from poolfast or find-graph could potentially help us find better results. 
  d. Potentially try to estimate networks from sequencing data without the estimated gene trees? Would it be better? Some tools like PhyNEST. PhyNest is also a julia package so might be easier to install and run? 

Some notes for another interesting paper : Cherry-picking method: https://www.researchgate.net/publication/373979203_Constructing_phylogenetic_networks_via_cherry_picking_and_machine_learning I am not sure if I undersrand all of these but it is interesting. There is no software for this: 

They introduced Cherry Picking Heuristics (CPH), heuristics to combine a set of binary phylogenetic trees into a single binary phylogenetic network based on cherry picking. A cherry is a pair of taxa that share the same parent. This algorithum remove one taxa from a pair of cherries. If the cherry is reticulated, one approach is to replace the edges connected to the hybrid node with simpler, direct edges to preserve the hybridization information. For example, if A and B are reticulated cherries, this method reduce this reticulation into a single new edge. Sometimes, only one edge may remain between the parent and the reduced node, creating a simpler representation of the hybridization (? How? Don't understand what it doesn't lose reticulate info). The cherry-picking process distinguishes between regular cherries (which share a single parent) and reticulated cherries (where one member is descended from a reticulate or hybrid node).
Cherry-Picking Sequence (CPS): A CPS is an ordered sequence of cherries (pairs of leaves) that, when reduced step-by-step, will eventually collapse the network into a single root node with one descendant leaf. Once all cherries in the CPS are reduced, the resulting network ideally captures the main evolutionary relationships with minimal reticulation nodes, providing an efficient structure that balances accuracy and complexity.






############## Below are to-do list for Max ##############  
Todo This Week:

  *use new tree with coalecent unit and substitution unit
  *work on finalizing script with an explanation of how to run
  *work on tutorial for installing executables
  *run simulation
  *Document Document Document
  *add raxml to julia script and push to github
  *make a figure and put in the readme file to help show visual representation with sequence file for adding to resume (show pipeline)
  *phyloplots coalecent units only
  *document what you are doing for employers
  *store astral results


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