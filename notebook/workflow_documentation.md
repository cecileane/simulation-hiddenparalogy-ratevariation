# Documentation for seed generators 
Step 1: A master_seed is generated based on each unique combination of parameter settings (dup_rate, loss_rate, ratevar, n_inds) 
Step 2: Each software gets its own seed generated from this master_seed using stableRNG to make sure it is reproducible. 
Step 3: Then, each seed_"software" gets a seed_array (m x n) generated from stableRNG, which will be used for different runs (or iterations). 
    In simulation_iqtree.jl: 
        1. seed_simphy: 
            For seed_simphy, a seed_array (m x n) is used for different replicate (n_reps = m) and different iteration (n_iterations = n) (details see simulation_iqtree.jl and utilities.jl) 
            Seeds used are stored in "random_seed_simphy.txt" 
        2. seed_seqgen: 
            For seed_seqgen, a seed_array (m x 1) is used each replicate (n_reps = m) 
            Seeds used are stored in "random_seed_seqgen.txt" 
        3. seed_iqtree: 
            For iqtree, a seed_array (m x 1) is used each replicate (n_reps = m) 
            Seeds used are stored in "random_seed_iqtree.txt" 
        PS: Astral, which is also used in simulation_iqtree.jl, is deterministic, so no seed used. 
    In findgraphs.jl: 
        1. The master_seed (step one) based on parameter setting is used to generate a seed_findgraphs. Then, a seed_array is generated with m (=num of replicate) x 2 (num_admix = 0 and 1) dimension. When running findgraphs in findgraphs_1rep.R, seed_array[m, 1] is used for m replicate and num_admix = 0, seed_array[m, 2] for num_admix = 1. 
    
     
# Documentation for reasoning behind findgraph 
Workflow 1: 
1. In each replicate, run find_graphs for 100 independent runs with K = 0 and K = 1.
    Reasoning: 
    a. In Maier et al. for their simulation dataset, they knew the real K and they tested with n-1, n and n+1. This is also used in another study with simulated data (Flegontov et al. 2023). Here, we know our tree with k = 0, and we only want to reject the null hypothesis k = 0, so running with k = 0 and k = 1 is probably okay. 
    b. 100 independent runs is used in Maier et al. 

2. In each replicate and each run, find 5 graphs with the best (lowest) LL score from both Ks. Thus, there will be 10 graphs for each run (k = 0 and k = 1) and each replicate has 100 runs, so (1,000 graphs per replicate). 
    Maier et al. used the same method of running for 100 runs and chose the top 5 graphs with best LL within each run. However, this indeed generates a lot of graphs for our data. Assume we have 100 replicate for each parameter setting, this will generate 1,000 * 100 = 100,000 graphs in total for each parameter setting. 
    Thus, I think an alternative is to keep the best graphs for each run in find_graphs, which can bring the total number of graphs in each replicate (for k = 0 and k = 1) to 2 * 100 * 100 = 20,000. One worry I have is that in Maier et al., they specifically crtiques the drawbacks of only having the best fitted graphs. 

3. Pool all graphs from the same replicate together, then discard topologically redundant graphs within each replicate. To reduce searching time, we first merge graphs from all runs in the replicate first, and start to discard topologically redundant graphs. Then, we merge graphs from all replicates, and delete topologically redundant graphs iteratively. Hopefully this will reduce the number of total graphs. 

Notes from Feb 28: When pooling all the graphs together, ideally most of the runs could have many redundant graphs compared to another runs, so it is very likely to throw a lot of graphs. Eespecially k = 0. Thus, it is less likely to have too many graphs. The final list would have 5 to 500 graphs. 

4. Use qpgraph to find |WR| among the the remaining graphs, in set with both k = 0 and k = 1, if there is any graph with |WR| < 3se, then we set the lowst K as the number of K we found. Here, if k = 0 set has any graphs with |WR| <3, then no migration is found. If K = 1 has graphs with |WR| < 3 while k = 0 doesn't, then we find one migration. 

This process is described in the discussion of Maier et al. about how to chose migration number. They mentioned a threshold of WR or LL should be decided to find the lowest K number which meet this threshold. Since we used LL to select the best graphs, here we could use |WR < 3| to select the K. This threshold is also widely used in many emprical studies (eg. Gutaker et al. 2020 using |f4-statistic z-scores| < 3 with qdgraph in AdmixtureTools1 and Flegontov et al. 2023 with findgraph in AdmixtureTools2, ) and Maire et al iteself. 

If any graph places the outgroup not as the outgroup, then this graph is not possible. 

If no graph has wr < 3, then migration is higher than one. 

5. Across all replications within a particular parameter setting, we can calculate how many replicates have falsely inferred K = 1 or higher (of no graph with wr < 3), which is the type I error rate. We calculate the number of replicates which are correctly inferred as K = 0, which is the probability of if we can fit the graph correctly. This is not power, since power is the accept the alternative. 

If we only want to know if we want to reject the null, we only need to do k = 0. 
If we want to know what kind of graphs we could infer, we then need to do k = 0 and k = 1. 

Let's see how long it takes to run K = 0. Let me code k = 0 first, and then test the speed. Have the code to run k = 0 or k = 1, separately. 

--Double check the function of booststrapping. 

Workflow 2: 
1. In each replicate, run find_graphs for 100 times with K = 0 and K = 1 (same as workflow 1).  

2. In each run in each replicate, find 5 graphs with the best LL score from each K (0, 1). 10 graphs for each replicate (same as workflow 1). 

3. Pool all graphs from all runs together and discard the topologically redundant graphs. 

4. Maire et al. introduced a cross-validating boostrapping approach to compare two models with different complexity (different Ks) and a boostrap methods to compare models with the same complexity. Here, we have 5 models from K = 0 (m0 set: m0_1, m0_2, ... m0_5) and 5 models from k = 1 (m1 set: m1_1, m1_2, ..., m1_5). Our grand truth model is a tree structure without any migration. 

5. Calculate |WR| for all remaining trees. 

6. Using the boostrap methods described in Maier et al. to test between the grand truth model and all models with k = 0. This tells us if migration = 0, 

This tells us how many times we could reject our null hypothesis that we should pick K = 0. 

However, the cross-validating boostrapping approach to compare models with different complexity is documented in Mainer et al (Appendix 1.B.3 second paragraph, and Appendix 2.E last two paragraphs). I didn't think there is any documentation describing how to use it on AdmixtureTools2 Github or any webpage. I suspect this should be mannually implemented in codes perhaps but it sounds very fun. Below I listed how they conducted the boostrapping approach for models with the same complexity and models with different complexity. 

Same K model comparison: 
    a. Divide genome into n block indexed by i 
    b. Select b blocks with replacement, indexed by j 
    c. Fit both graphs b times once for each boostrap set of b SNP blocks 
    d. A set of b score difference deltaj. 
    e. Boostrap confidence interval for difference in scores is given by the distribution of deltaj. 
    d. Compyte emprical boostrap p-values (see appdenix 2.e)

Different K model comparison: 
First, more complex models have more degree of freedom and could be over-fitted. 
    a. Implement out-of-sample likelihood score, which have observed f3 and expected f3 statistics defined on mutually exclusive sets of SNP blocks. --> I need some help understanding this (appendix 2). 

The function for boostrap: https://uqrmaie1.github.io/admixtools/reference/qpgraph_resample_multi.html 


# Paper: https://www.science.org/doi/10.1126/science.ade2833 
**Seehausen, O., Meier, J. I., Marques, D. A., Wagner, C. E., Excoffier, L., & Malinsky, M. (2023). Cycles of fusion and fission enabled rapid parallel adaptive radiations in cichlid fishes. Science, 381(6656), eade2833. https://doi.org/10.1126/science.ade2833** 

First, the phylogenetic inferrance of this paper is interesting. They used SNP data only and used IQ-tree to fast the speed. This part could be checked out later. 

In their methods, the authors used admixture graph reconstruction to estimate the Congo-Nilotic admixture proportions in the Lake Victoria Region Superflock (LVRS). They applied the “qpgraph” function from the ADMIXTOOLS 2.0.0 R package, which models demographic history by incorporating both population drift and gene flow (admixture edges) into phylogenetic trees. To validate model robustness, the authors applied block bootstrapping. 

This paper indeed uses the bootstrapping method to compare models with different complexity. In their supplementary materials Figure S13 C and D, they compared models with two and one admixture edges. They found that models with two edges have higher likelihood than models with one edge, but a model comparison based on block-bootstrapping reveals that this difference is not significant (p-value: 0.157).  

This paper could be cited for using the block-bootstrapping method. 


# Paper: 
**Gutaker, R.M., Groen, S.C., Bellis, E.S. et al. Genomic history and ecology of the geographic spread of rice. Nat. Plants 6, 492–502 (2020). https://doi.org/10.1038/s41477-020-0659-6** 

*Note*: This paper used qdgraph in admixturetool but not findgraphs in admixture. When they chose models, they used |f4-statistic z-scores| < 3. This is from qd_graph in admixturetools1. 

*Paper summary*: The study reconstructs relationships among japonica and indica rice subpopulations using admixture graph analysis.They used admixturetool (version one but ADMIXTURETOOL2) and the QPgraph function to reconstruct the graph. The study reconstructs the geographic spread of rice using whole-genome resequencing data from more than 1,400 landraces.

*Methods to infer admixture graphs*: 
First, kd is defined as the number of subpopulations (clusters) without outgroup and k is number of clusters with outgroups. The authors tested different admixture graph models for japonica and indica rice by varying the number of subpopulations (kd) and migration events to identify the best-fitting models.
1. For kd = 2 to 5 (3 to 6 subpopulations including outgroup)
    * They explored all possible models with 0, 1, or 2 migration events.
    * Only models with f4-statistic z-scores <3.0 (indicating a good fit) were retained.
2. For kd = 7 to 9
    * They started with 6 subpopulations and tested models with 0, 1, or 2 migrations, again keeping only those with f4-statistic z-scores <3.0.
    * They then progressively added subpopulations to different positions in the graph and re-tested the models, keeping only those with f4-statistic z-scores <3.0.
    * This process continued until either no more subpopulations could be added or no more valid models ( f4-statistic z-scores <3.0) were found.
    * If no models with z-scores <3.0 remained, they considered those with z-scores <10.0.
    * They introduced an additional admixture event in all possible locations in the graph and tested the new models.
    * Only models with z-scores <3.0 were retained. 

Since many models met this criterion, they grouped them into three major "topology groups" and selected the best representative models for each. I don’t really understand how do they summarized their tree into topology group. Based on s.figure 16 and 22, they might just classified the topologies based on sturctural similarities, but I couldn’t tell with my bare eyes in their supplementary figures. 

Supplementary 16 and 22 are important to understand the process: 
In supplementary figures, from kd = 3 to kd = 8, they showed the summary of admixture graphs. Two things to notice: 
	1.  Models with |f4-statistic z-scores |< 3  are chosen 
	2. They summarized chosen models into three topology groups starting from kd = 5, which I don’t understand how do they summarized into three topology groups. 

# Paper: 
**Flegontov, P., Işıldak, U., Maier, R., Yüncü, E., Changmai, P., & Reich, D. (2023). Modeling of African population history using f-statistics is biased when applying all previously proposed SNP ascertainment schemes. PLOS Genetics, 19(9), e1010931. https://doi.org/10.1371/journal.pgen.1010931** 

*Notes*: They used the ADMIXTURETOOLS2 (using f4 statistics). In the paper, they argued that methods like find_graphs or TreeMix, which estimate AGs from SNPs, could be problematic due to Maier et al, saying that "f-statistics do not constrain even moderately complex topology spaces (e.g., graphs including 8–9 groups and 4–5 admixture events) well enough, and topologically diverse graphs often fit the data significantly better than true simulated histories." Instead, they used f4-statistics and qdAdm to fit AGs to calculated f4 statistics data. However, in one of their simulated datasets, they used findgraphs to find poorly fitted graphs and used them as poorly fitted AGs. 

This part is something I don't undersstand, because based on Maier et al, find_graph is also a greedy approach which explores all possible AGs and find how good those AGs could fit to the f2-statistics. 

They also use WR < 3 as the cut-off to choose well-fitted models. 

It is noted that they didn't explore graphs using find_graphs, instead they fit all possible topologies to their datasets and check if the model got rejected or accepted. I kind of understand how this could be achieved with their simulated dataset since they know the true AGs but I don't understand how this could be useful to their real datset.  

Something that could be useful to my study is that when they were using find_graphs. They run findgraphs for 300 runs in total, seeded by random graphs containing either the simulated number of admixture events (n, 100 runs), or n-1 events (100 runs), or n+1 events (100 runs). Those number of runs could be an reference. 

*Summary of the paper*: SNP ascertainment schemes refer to the methods used to select which single nucleotide polymorphisms (SNPs) are included in a genetic analysis. Because SNP ascertainment typically doesn't come from whole genome data and came from target enrichment, these schemes determine how and from which populations SNPs are identified, which can influence the accuracy of downstream population genetic analyses. The paper found that commonly used ascertainment panels in human genetics introduce biases when analyzing human groups using f4 statistics, leading to false rejection of correct demographic models and failure to reject incorrect ones. 

*AG inference for their simulated dataset*:  
Simulation 1: Simualting SNP sets based on real human data (section Method 1.1 and 1.2): 
They simulated different populations of unascertained and ascertained SNP sets. 
Simulation 2: Simulating random admixture graphs and simple trees (Method 1.3): 
They simulated random AGs of four complexity classes: including 8 or 9 sampled non-outgroup populations, one outgroup, and 4 or 5 pulse-like admixture events. They selected only simulations where pairwise FST for groups were in the range characteristic for anatomically modern and archaic humans (at least one FST value below 0.15). 

How did they fit they models? WR > 3 -> rejected models and WR < 3 -> accepted models 
They fit all possible AGs with two migrations to the three datasets they have with different ascertainment schemes. Then, they examined the fits of collections of AGs. The authors analyzed 5,000 best-fitting admixture graph topologies for each population quintuplet, using unascertained site sets to compare model performance. They evaluated two admixture graph fit metrics:
Log-likelihood (LL): More accurate but difficult to compare across datasets.
Worst residual (WR): Measures how poorly an admixture graph fits f-statistics.
To assess ascertainment bias, they compared graph fits on ascertained vs. unascertained data using Pearson correlation. Additionally, they measured:
Bias: The fraction of models rejected under ascertainment (WR >3 SE) but accepted on full data (WR <3 SE).
Power to reject incorrect models: The fraction of models accepted under ascertainment (WR <3 SE) but rejected on full data (WR >3 SE).
This approach helped them determine how SNP ascertainment influences the accuracy of population history reconstructions.

# Paper 3: 
**Robert MaierPavel FlegontovOlga FlegontovaUlaş IşıldakPiya ChangmaiDavid Reich (2023) On the limits of fitting complex models of population history to f-statistics eLife 12:e85492.** 




