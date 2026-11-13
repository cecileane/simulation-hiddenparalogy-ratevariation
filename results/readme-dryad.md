# Data from: Substitution rate variation, not hidden paralogy, drives false hybridization signal in phylogenetic network inference

Summary tables of the simulations, with the column dictionaries. The code
that produced them is on
[GitHub](https://github.com/cecileane/simulation-hiddenparalogy-ratevariation),
see section 3.

## 1. The simulations

We simulated phylogenomic data sets to test whether gene duplication and
loss (which can leave hidden paralogs in the data) and substitution rate
variation can make network inference methods (find_graphs, SNaQ, PhyNEST)
detect a reticulation when the true history is a tree. 
Two 8-taxon species trees were used: 
a reptile tree calibrated on the UCE data of Crawford et al. 2012, 
and a basal angiosperm tree calibrated on the Angiosperms353
genes of the Kew Tree of Life Explorer (Zuntini et al. 2024). 

Each tree simulation has 
36 settings per tree, 100 replicates each, 1000 genes of 1000 bp:

- `RV`: substitution rate variation, N (none), G (across genes), L (across
  lineages)

- `DUP` = `LOS`: duplication and loss rate per gene per generation. 0 (no
  hidden paralogy possible), medium and high: 3e-4 and 4e-4 for the
  reptile tree, 1e-4 and 1.4e-4 for the angiosperm tree

- `SF`: scaling factor of the population size, 1.0 (high ILS) or 0.5
  (low ILS)

- `N_ind`: 1 or 2 individuals per species

- `genelen`: 1000 in all settings

A setting is named `DUP<d>-LOS<l>-RV<r>-N_ind<n>-SF<s>-genelen<g>`, for
example `DUP0.0003-LOS0.0003-RVG-N_ind1-SF0.5-genelen1000`. PhyNEST was run
on the 18 settings per tree with one individual per species.

For each replicate: SimPhy simulates the gene trees, 
the trees with more than one copy per individual or with missing taxa are dropped, 
Seq-Gen simulates the sequences,
IQ-TREE estimates the gene trees, ASTRAL the species tree. 
Then SNaQ (on the gene trees), 
find_graphs (on the SNPs of the concatenated alignment)
and PhyNEST (on the concatenated alignment) each fit a tree (h = 0) 
and a network with one reticulation (h = 1) (see paper and git repo for details). 

## 2. Files

```
readme-dryad.md                 this file
summary_reptile/<method>/*.csv  per-replicate summaries, reptile tree, 198 files
summary_plant/<method>/*.csv    same for the angiosperm tree, 198 files
tables/*.csv                    the supplementary tables (10 files)
tables/per_setting.csv          the numbers behind Fig. 1
tables/T_star.csv               the PhyNEST thresholds
```

The data of the paper are the per-replicate summaries: 
one csv per setting and method, 
one row per replicate (36 files per method and tree,
18 for PhyNEST). 
Every figure and table of the paper is computed from
these files by the R scripts in `scripts/visuals/` of the repository.

Taxon names: 
A to H, A is the outgroup (Homo in the reptile tree, Pinus in
the angiosperm tree). 
The "alternative trees" are the three most likely
wrong trees: the true tree with the taxon on the shortest branches moved
(F for the reptile tree, B for the angiosperm tree), 
see scripts/snaq_postprocess.jl. 
The `_noF`, `_noG`, `_noH` columns repeat a
comparison after pruning one taxon from both trees: 
F, G, H for the reptile tree and B, F, G for the angiosperm tree 
(the column names were kept). 
RF distances are the Robinson-Foulds distances computed by
`mudistance_semidirected` of PhyloNetworks, 0 means the same topology.

NaN means that the quantity is not defined for that replicate or setting,
instead of a failed computation. 
For example, the `_loss_only` columns of the simulation
summary are all NaN: 
a tree with a lost taxon is dropped, so no kept tree has loss without duplication. 
Another example is that the hidden paralogy columns are NaN when DUP = LOS = 0.

### 2.1 simulation/summary_<setting>.csv

One row per setting (not per replicate), 
from scripts/simulation_postprocess.jl. 
Hidden paralogy (HP) categories of a kept gene tree: `nothing` 
(no duplication or loss), `false_HP`
(duplication and loss, but the locus tree has the topology and internal
branch lengths of the species tree), 
`weak_HP` (same topology, different
branch lengths), 
`strong_HP` (different topology). 
`dup_and_loss` is false_HP + weak_HP + strong_HP.

| column | meaning |
|---|---|
| `parameter_setting` | setting name |
| `n_genes_mean`, `n_genes_min`, `n_genes_max` | gene trees kept per replicate (1000 in all replicates except a few at the highest rate) |
| `n_iterations` | mean number of SimPhy runs needed to get 1000 trees |
| `n_repeated_taxa_removed` | mean number of trees dropped for having 2 copies of a gene in one individual |
| `n_insufficient_taxa_removed` | mean number of trees dropped for missing taxa |
| `percentage_genes_meet_min` | % of replicates that reached 1000 trees |
| `pert_trees_experiencing_gene_loss_only`, `..._gene_duplication_and_loss`, `..._nothing` | % of kept trees in each category (`pert_` stands for percent) |
| `pert_false_HP`, `pert_weak_HP`, `pert_strong_HP` | % of kept trees in each hidden paralogy category |
| `mean_RF_true_and_estimated_species_trees` | mean RF between the true tree and the ASTRAL tree |
| `num_estimated_species_tree_diff_from_truth` | replicates whose ASTRAL tree is not the true tree |
| `mean_num_taxa_<category>` | mean number of tips of the kept gene trees, by category |
| `mean_internal_bl_locus_<category>` | mean internal branch length of the true gene trees (SimPhy locus trees), by category |
| `mean_internal_bl_gene_<category>` | same for the gene trees estimated by IQ-TREE |
| `mean_RF_genetree_vs_sptree_<category>` | mean normalized RF between the true gene trees and the species tree, by category: gene tree discordance |
| `mean_RF_genetree_vs_sptree_noF_<category>` | same without taxon F (B for the angiosperm tree) |
| `mean_RF_true_vs_est_genetree_<category>` | mean normalized RF between the true and the estimated gene trees, by category: gene tree estimation error |

### 2.2 genetrees/genetree_tiers_<setting>.csv

From scripts/genetree_tiers.jl, one row per locus category (`all`,
`ortholog`, `falseHP`, `weakHP`, `strongHP`, `lossonly`). 
Behind the supplementary tables on gene tree discordance and estimation error.

| column | meaning |
|---|---|
| `setting`, `category` | setting name and locus category |
| `n_genes` | gene trees in the category, over the 100 replicates |
| `mean_RF_species` | mean normalized RF between the true gene trees and the species tree |
| `d0` to `d5` | fraction of gene trees with 0, 1, ... 5 branches discordant with the species tree |
| `mean_RF_estimated` | mean normalized RF between the true and the estimated gene trees |
| `e0` to `e13` | fraction of gene trees with 0, 1, ... 13 misestimated branches (5 internal branches with one individual per taxon, 13 with two) |

With two individuals per taxon a gene tree has up to 2^8 8-taxon subtrees
(one individual kept per species); each subtree counts 1 / (number of
subtrees) in the fractions.

### 2.3 snaq/SNaQ-<setting>-summary.csv

From scripts/snaq_postprocess.jl, one row per replicate.

| column | meaning |
|---|---|
| `repID` | replicate, 001 to 100 |
| `p_H0`, `p_H1` | p-value of the goodness of fit test (QuartetNetworkGoodnessFit) of the h = 0 tree and the h = 1 network. h = 0 is kept when p_H0 > 0.05, h = 1 when p_H0 <= 0.05 and p_H1 > 0.05, else h > 1 |
| `z_uncorrected_H0`, `z_uncorrected_H1`, `sigma_H0`, `sigma_H1` | z statistic and its correction factor behind each p-value |
| `score_H0`, `score_H1` | SNaQ pseudolikelihood score of each network (lower is better) |
| `gamma_1`, `gamma_2` | major and minor inheritance probabilities of the h = 1 network, gamma_1 + gamma_2 = 1 |
| `RF_net0_true`, `RF_net0_true_noF` | RF between the h = 0 tree and the true tree, full and with one taxon pruned |
| `RF_net1_1_true`, `RF_net1_2_true` | RF between the major (1) or minor (2) tree displayed by the h = 1 network and the true tree |
| `RF_net1_1_true_noF`, `_noG`, `_noH`, same for `RF_net1_2` | same with one taxon pruned |
| `RF_net0_alter1`, `_alter2`, `_alter3` | RF between the h = 0 tree and the three alternative trees |
| `RF_net1_1_alter1` to `RF_net1_2_alter3` | same for the two displayed trees |
| `hybrid_taxon`, `major_donor`, `minor_donor` | taxa below the hybrid node of the h = 1 network and below its major and minor parent (network rooted at A), NA if the network has no hybrid |

### 2.4 findgraphs/findgraph-<setting>.csv

From scripts/findgraphs_postprocess.jl, one row per replicate. 
find_graphs returns several graphs per search (100 searches per replicate, 
with 0 and with 1 admixture event), 
kept when within 2 log likelihood units of the best one, duplicates removed. 
"best" = lowest score. k = number of admixture events = h.

| column | meaning |
|---|---|
| `repID` | replicate |
| `best_k` | model choice with the worst residual (WR) threshold 3.0: "0" if a k = 0 tree has WR < 3, else "1" if a k = 1 graph has, else ">1" |
| `best_k_new_WR_3.7` | same with the threshold 3.7 |
| `H0_trees_found`, `H1_graphs_found` | number of distinct graphs kept with k = 0 and k = 1 |
| `true_tree_wr` | worst residual of the true species tree, fitted to the same data |
| `H0_best_tree_is_true_tree`, `_noF`, `_noG`, `_noH` | is the best k = 0 tree the true tree, full and with one taxon pruned |
| `H0_best_tree_WR` | worst residual of the best k = 0 tree |
| `H1_best_graph_displayed_true_tree`, `_noF`, `_noG`, `_noH` | does the best k = 1 graph display the true tree (as its major or minor tree) |
| `H1_best_graph_WR` | worst residual of the best k = 1 graph |
| `true_tree_found_in_H0`, `_noF`, `_noG`, `_noH` | is the true tree among all the k = 0 trees kept |
| `avg_gamma1_H1`, `avg_gamma2_H1` | mean major and minor admixture weight over the k = 1 graphs kept |
| `best_graph_gamma1`, `best_graph_gamma2` | admixture weights of the best k = 1 graph |
| `True_tree_displayed_H1_major`, `True_tree_displayed_H1_minor`, and `_noF`, `_noG`, `_noH` versions | is the true tree the major (or minor) tree of any k = 1 graph kept |
| `best_graph_hybrid_taxon`, `best_graph_major_donor`, `best_graph_minor_donor` | taxa below the admixed node of the best k = 1 graph and below its major and minor parent |

### 2.5 findgraphs_bestgraph/findgraph-bestgraph-<setting>.csv

From scripts/findgraphs_bestgraph_trees.R, one row per replicate: the two
trees displayed by the best k = 1 graph.

| column | meaning |
|---|---|
| `rep_id`, `score`, `worst_residual` | replicate, score and worst residual of the best graph |
| `gamma_major`, `gamma_minor` | its admixture weights |
| `hybrid_taxon`, `major_donor`, `minor_donor` | as in 2.4 |
| `major_tree`, `minor_tree` | newick of the two displayed trees |
| `rf_major`, `rf_minor` | their RF distance to the true tree |
| `major_wrong_clades`, `minor_wrong_clades` | the clades of each displayed tree that the true tree does not have |

### 2.6 phynest/PhyNEST-<setting>-summary.csv

From scripts/phynest_postprocess.jl, one row per replicate, 18 settings
per tree (one individual per species). 
PhyNEST scores are negative log composite likelihoods over all sites of 
the concatenated alignment, 
so the score difference grows with the number of sites M; 
the statistic used in the paper is the per-site difference `T_per_site`. 
h = 1 is selected when T_per_site > T*, 
where T* is the 95th percentile of T_per_site over the
settings without lineage rate variation, one T* per tree and ILS level
(computed by the R scripts, written to `tables/T_star.csv`).

| column | meaning |
|---|---|
| `repID` | replicate |
| `score_H0`, `score_H1` | score of the best h = 0 tree and the best h = 1 network over the 100 search runs |
| `T` | score_H0 - score_H1 |
| `n_sites` | number of sites M of the alignment, from the PhyNEST log (10^6 in most replicates, fewer when fewer loci were kept) |
| `T_per_site` | T / n_sites |
| `score_truetree` | score of the true species tree with its parameters optimized |
| `num_hybrid_H1` | hybrid nodes in the h = 1 network (1) |
| `gamma_1`, `gamma_2` | major and minor inheritance probabilities of the h = 1 network |
| `RF_net0_true`, `_noF`, `_noG`, `_noH`, `RF_net0_alter1` to `_alter3` | RF of the h = 0 tree to the true tree (full, pruned) and to the alternative trees |
| `RF_net1_1_*`, `RF_net1_2_*` | the same seven comparisons for the major (1) and minor (2) tree displayed by the h = 1 network |
| `hybrid_taxon`, `major_donor`, `minor_donor` | as in 2.3 |
| `n_runs_H0`, `n_runs_H1` | search runs with a finite score (a failed optimization gives NaN) |
| `n_runs_at_best_H0`, `n_runs_at_best_H1` | runs within 1 score unit of the best run |
| `n_runs_true_H0`, `n_runs_true_H1` | runs whose tree is the true tree (h = 0) or whose network displays it (h = 1) |
| `gap_best_true_H0`, `gap_best_true_H1` | score of the best such run minus the best score (0 when the best run has the true tree, NaN when no run has it) |
| `rank_best_true_H0`, `rank_best_true_H1` | rank of the best such run, 1 = best, NaN when no run has it |
| `pct_runs_better_than_truetree_H0`, `_H1` | % of runs scoring better than the true tree by more than 1 unit |

### 2.7 Tables of the paper

`visual_tables.R` pools the replicates of all settings sharing a factor
level (rate variation, duplication/loss, ILS, individuals per taxon) 
and writes one csv (here, in `tables/`) and one tex per table, 
with a reptile and a plant column under every statistic: 
`wr_summary` (worst residual of the true tree), 
`T_summary` and `T_star` (PhyNEST), `model_choice` 
(% of replicates selecting h = 1 and h > 1 with each
method), `genetree_discordance`, `genetree_estimation_error`,
`gamma_<method>` (minor gamma) and `bestgraph_summary`. The csv are in
long format: tree, factor, level, n, statistic, value.
`visual_main_figures.R` writes `per_setting.csv`, the numbers behind
Fig. 1, one row per setting and method.

## 3. Code

The GitHub repository has the whole pipeline: `scripts/readme.md` says how
to run each step. 
Software: Julia 1.11.7 with the packages of Project.toml
(PhyloNetworks 1.3.1, SNaQ 1.1.1, QuartetNetworkGoodnessFit 1.0.0,
PhyloPlots, RCall, CSV, DataFrames, Distributions, StableRNGs and others,
pinned in Manifest.toml), 
PhyNEST 0.1 in its own environment
(envs/phynest), SimPhy 1.0.2, Seq-Gen 1.3.5, IQ-TREE 2.4.0, ASTER
(ASTRAL-IV) 1.24.4.8, snp-sites 2.5.1, 
Python 3.12 with biopython, R with
admixtools 2.0.8 and igraph 1.6.0. All on Linux x86-64.

## 4. License

CC0 (Dryad policy).
