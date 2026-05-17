# Data from: Substitution rate variation, not hidden paralogy, drives false hybridization signal in phylogenetic network inference

This folder contains the cross-setting summary tables that underlie
figures and tables in the paper.
The code used to produce these tables from raw simulations is on
[GitHub](https://github.com/cecileane/simulation-hiddenparalogy-ratevariation)
and [Dryad](https://doi.org/10.5061/dryad.dfn2z35h7).

---

## 1. Data description

We simulated phylogenomic datasets to test whether two common model violations
— gene duplication/loss (which can produce *hidden paralogy*) and substitution
rate variation — can cause network inference methods (`find_graphs` and `SNaQ`)
to detect spurious hybridization when the true history is a tree. Branch
lengths and substitution-rate parameters were calibrated from empirical reptile
UCE data ([Crawford et al. 2012](doi:10.1098/rsbl.2012.0331)).

**Experimental design — 36 parameter settings, 100 replicates each:**

| Factor | Values | Meaning |
|---|---|---|
| `DUP` / `LOS` | 0.0, 3 × 10⁻⁴, 4 × 10⁻⁴ | Gene duplication / loss rate (per gene, per generation). Duplication and loss rates are equal in every setting. `0.0` disables both, so no hidden paralogy can arise. |
| `RV` (rate variation) | N, G, L | None / Gene-specific / Lineage-specific substitution rates |
| `N_ind` | 1, 2 | Individuals sampled per species / taxon |
| `SF` (scale factor) | 0.5, 1.0 | Effective-population-size scaling factor; higher = more incomplete lineage sorting |
| `genelen` | 1000 | Gene length, in base pairs (legacy parameter, 1000 across all settings) |

Each setting is encoded in a string of the form
`DUP<d>-LOS<l>-RV<r>-N_ind<n>-SF<s>-genelen<g>` (e.g., `DUP0.0003-LOS0.0003-RVG-N_ind1-SF0.5-genelen1000`).
This string is the row identifier.
`genlen` is a legacy parameter which was disabled during our simulation.
In our simulation, all `genlen = 1000` (number of base pairs per gene).

Per setting we ran 100 replicates of the full pipeline:
SimPhy → paralogy / missing taxa filter → Seq-Gen → IQ-TREE → ASTRAL-IV → SNaQ + find_graphs.

---

## 2. Files and variables

This folder contains 6 top-level CSV files, plus two subfolders of
per-setting per-replicate CSVs (`findgraph_summary/` and `snaq_summary/`,
36 files each — see §2.6). The six top-level files break down as: three
primary cross-setting summary tables, and three pre-computed cross-tabulations
used to make the paper's combined figure. The two subfolder contain primary 
statistics summarized from raw graphs from SNaQ and findgraphs. All together, 
those result files can be used to reproduce figures and tables in our manuscript. 

All six top-level CSV files have **one row per parameter setting (36 rows)**,
*except* `combined_hypothesis_acceptance_marginal.csv` which has 10 rows (one
per factor level). Each file uses comma separators and a single header row.

### 2.1 Missing-data convention

`NaN` (Not a Number) is used throughout to indicate **"this metric is not
defined under this parameter combination"**. It is *never* used to indicate
a measurement failure. There are two reasons it appears:

- **`*_dup_and_loss`, `*_false_HP`, `*_weak_HP`, `*_strong_HP` columns are NaN
  in the 12 rows where `DUP/LOS = 0.0`.** With no duplication or loss events,
  no gene tree can fall into the dup-and-loss bucket or any of the hidden-
  paralogy buckets. The category is structurally empty.
- **The six `*_loss_only` columns are NaN in *all 36 rows*.** These columns
  are quality assurance placeholders from an earlier design that allowed missing taxa
  in simulated gene trees. The final published design requires every simulated
  gene tree to retain all taxa (gene trees with any loss are dropped during
  paralogy filtering, and SimPhy is re-run until enough complete trees are
  obtained — see the paper, Methods). Trees with loss but no duplication
  therefore never enter the analysis, so the "loss only" bucket is always
  empty. The columns are kept for schema stability with intermediate per-rep
  outputs that still record the per-tree classification.

### 2.2 `findgraph_summary.csv` — find_graphs cross-setting summary

One row per parameter setting (36 rows). Columns:

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `paramname_root` | string | Parameter-setting identifier (`DUP*-LOS*-RV*-N_ind*-SF*-genelen*`). |
| 2 | `total_replicates` | int | Number of replicates that contributed to this row (100 in every setting). |
| 3 | `H0Accepted` | int | Replicates where the best find_graphs model chosen by the *standard* worst-residual (WR) threshold (WR ≤ 3.0) had h = 0 admixture events. |
| 4 | `H1Accepted` | int | As above, but model chosen had h = 1. |
| 5 | `BT1Accepted` | int | As above, but model chosen had h > 1 ("BT1" = "best h beyond 1"). `H0Accepted + H1Accepted + BT1Accepted = total_replicates`. |
| 6 | `pct_true_tree_H0` | float (%) | % of replicates whose best h = 0 tree equals the true species tree (Robinson-Foulds distance = 0). |
| 7 | `avg_gamma1_H1` | float ∈ [0, 1] | Mean major edge weight (γ₁) of the best h = 1 graph across replicates. |
| 8 | `avg_gamma2_H1` | float ∈ [0, 1] | Mean minor edge weight (γ₂ = 1 − γ₁) of the best H = 1 graph. |
| 9 | `avg_best_gamma1` | float | Mean γ₁ of the single best-likelihood h = 1 graph per replicate. |
| 10 | `avg_best_gamma2` | float | Mean γ₂ of the single best-likelihood h = 1 graph per replicate. |
| 11 | `pct_true_tree_H1` | float (%) | % of replicates where any h = 1 graph displays the true species tree as its major or minor backbone. |
| 12 | `pct_true_tree_H1_noF` | float (%) | Same as column 11 but counted *without taxon F* in both true species tree and the networks generated under H=1. |
| 13 | `mean_H0_trees_found` | float | Mean number of distinct H = 0 trees found per replicate after deduplication. |
| 14 | `median_H0_trees_found` | float | Median (per replicate). |
| 15 | `sd_H0_trees_found` | float | Standard deviation (per replicate). |
| 16 | `min_H0_trees_found` | int | Minimum across replicates. |
| 17 | `max_H0_trees_found` | int | Maximum across replicates. |
| 18–22 | `mean/median/sd/min/max_H1_graphs_found` | numeric | Same statistics for H = 1 graphs. |
| 23 | `count_H0_best_is_true_tree` | int | Replicates where the highest-likelihood H = 0 tree equals the true species tree (RF = 0). |
| 24 | `count_H0_best_is_true_tree_noF` | int | Same, without F taxon. |
| 25 | `count_H1_best_displays_true_tree` | int | Replicates where the highest-likelihood H = 1 graph displays the true species tree. |
| 26 | `count_H1_best_displays_true_tree_noF` | int | Same, without F taxon. |
| 27 | `H0Accepted_WR_3.7` | int | Like column 3, but using the *relaxed* WR threshold (WR ≤ 3.7). |
| 28 | `H1Accepted_WR_3.7` | int | As column 4, relaxed threshold. |
| 29 | `BT1Accepted_WR_3.7` | int | As column 5, relaxed threshold. |

### 2.3 `SNaQ_summary.csv` — SNaQ cross-setting summary

One row per parameter setting (36 rows). Columns:

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `parameter_setting` | string | Parameter-setting identifier (in this file prefixed with `SNaQ-…-summary.csv`, i.e. the source filename for the per-setting CSV). |
| 2 | `H=0Accepted` | int | Replicates accepting the H = 0 (tree-only) model. SNaQ uses a goodness-of-fit p-value; H = 0 is accepted when p<sub>H0</sub> > 0.05. |
| 3 | `H=1Accepted` | int | Replicates where H = 0 is rejected (p<sub>H0</sub> ≤ 0.05) but H = 1 is accepted (p<sub>H1</sub> > 0.05). |
| 4 | `H>1Accepted` | int | Replicates where both H = 0 and H = 1 are rejected. Sum across columns 2–4 = total replicates. |
| 5 | `mean_score_H0` | float | Mean SNaQ pseudo-likelihood score for the best H = 0 network across replicates (lower is better). |
| 6 | `mean_score_H1` | float | Mean SNaQ pseudo-likelihood score for the best H = 1 network. |
| 7 | `mean_gamma_1` | float ∈ [0, 1] | Mean major-edge weight γ₁ of the best H = 1 network. |
| 8 | `mean_gamma_2` | float ∈ [0, 1] | Mean minor-edge weight γ₂ = 1 − γ₁. |
| 9 | `find_true_net0` | int | Replicates where the best H = 0 network exactly matches the true species tree (RF = 0). |
| 10 | `find_true_net0_noF` | int | Same from column 9, with taxon F removed. |
| 11 | `find_true_net1` | int | Replicates where the best H = 1 network displays the true tree as either its major or minor backbone. |
| 12 | `find_true_net1_noF` | int | Same, with taxon F removed. |
| 13 | `mean_p_H0` | float ∈ [0, 1] | Mean goodness-of-fit p-value under H = 0 across replicates. |
| 14 | `mean_p_H1` | float ∈ [0, 1] | Mean goodness-of-fit p-value under H = 1. |
| 15 | `find_alter_net0` | int | Replicates whose best H = 0 network matches **any** alternative pre-specified non-true tree (sum across the three alternatives the paper considers). |
| 16 | `find_alter1_net1_major` | int | Replicates whose best H = 1 network's **major** backbone matches alternative tree #1. |
| 17 | `find_alter2_net1_major` | int | Same for alternative tree #2. |
| 18 | `find_alter3_net1_major` | int | Same for alternative tree #3. |
| 19 | `find_alter1_net1_minor` | int | Replicates whose best H = 1 network's **minor** backbone matches alternative tree #1. |
| 20 | `find_alter2_net1_minor` | int | Same for #2. |
| 21 | `find_alter3_net1_minor` | int | Same for #3. |

The three "alternative" topologies are pre-specified competing species trees
defined in `scripts/utilities.jl`; see the paper Methods for their definitions.

### 2.4 `summary_concatenated.csv` — simulation-side cross-setting summary

One row per parameter setting (36 rows). Records what the simulation itself
produced (gene-tree counts, hidden-paralogy categories, branch-length
statistics, gene-tree-vs-species-tree distances) *before* SNaQ/find_graphs ran.

**Hidden-paralogy (HP) categories** used in column names:

- `loss_only`     — gene trees that experienced loss but **no** duplication. Always empty in this design (see §2.1). This column is kept for quality assurance purpose for now. All values should be NaN since our simulation does not allow missing taxa so there should be no loss_only category. 
- `dup_and_loss`  — gene trees that experienced both duplication and loss. Always empty when DUP/LOS = 0.
- `nothing`       — gene trees with neither duplication nor loss.
- `false_HP`      — duplication and loss happened but neither internal branch length and topology was changed between the locus tree and the species tree .
- `weak_HP`       — duplication and loss happened and internal branch length was changed without changing the topology. 
- `strong_HP`     — the locus tree and species tree have different internal branch length and topology. 

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `parameter_setting` | string | Parameter-setting identifier. |
| 2 | `n_genes_mean` | float | Mean number of gene trees retained per replicate (target = 1000). |
| 3 | `n_genes_min` | int | Minimum per replicate. |
| 4 | `n_genes_max` | int | Maximum per replicate. |
| 5 | `n_iterations` | float | Mean number of SimPhy re-draws required to reach `n_genes_min` or target number of genes (averaged across replicates, see simulation.jl and our manuscript for details). |
| 6 | `n_repeated_taxa_removed` | float | Mean number of gene trees dropped per replicate due to within-species paralogs (the species had ≥ 2 copies of the gene). |
| 7 | `n_insufficient_taxa_removed` | float | Mean number of gene trees dropped per replicate for having too few remaining taxa to be phylogenetically informative. |
| 8 | `percentage_genes_meet_min` | float (%) | % of replicates that met the minimum-genes threshold without exhausting the SimPhy re-draw quota. |
| 9 | `pert_trees_experiencing_gene_loss_only` | float (%) | % of gene trees in the `loss_only` bucket (see HP categories above). Note: `pert_` is artifact in the column schema; read as "percent". |
| 10 | `pert_trees_experiencing_gene_duplication_and_loss` | float (%) | % in the `dup_and_loss` bucket. |
| 11 | `pert_trees_experiencing_nothing` | float (%) | % in the `nothing` bucket. |
| 12 | `pert_false_HP` | float (%) | % in the `false_HP` bucket. |
| 13 | `pert_weak_HP` | float (%) | % in the `weak_HP` bucket. |
| 14 | `pert_strong_HP` | float (%) | % in the `strong_HP` bucket. |
| 15 | `mean_RF_true_and_estimated_species_trees` | float | Mean Robinson-Foulds distance between the true species tree and the ASTRAL-IV estimate, across replicates. |
| 16 | `num_estimated_species_tree_diff_from_truth` | int | Replicates whose ASTRAL estimate differs from the true species tree (RF > 0). |
| 17 | `mean_num_taxa_all_genes` | float | Mean number of taxa per gene tree, across all gene trees in all replicates. |
| 18 | `mean_num_taxa_loss_only` | float | Same, restricted to `loss_only` trees. **All NaN by design — all NaN, see §2.1.** |
| 19 | `mean_num_taxa_dup_and_loss` | float | Same, restricted to `dup_and_loss` trees. NaN when DUP/LOS = 0. |
| 20 | `mean_num_taxa_nothing` | float | Same, `nothing` trees. |
| 21 | `mean_num_taxa_false_HP` | float | Same, `false_HP`. NaN when DUP/LOS = 0. |
| 22 | `mean_num_taxa_weak_HP` | float | Same, `weak_HP`. NaN when DUP/LOS = 0. |
| 23 | `mean_num_taxa_strong_HP` | float | Same, `strong_HP`. NaN when DUP/LOS = 0. |
| 24 | `mean_internal_bl_locus_all` | float | Mean of average internal branch length on the *locus* tree (the true SimPhy gene tree before sequence simulation), across all genes. Units: SimPhy's default branch-length units (substitutions per site). |
| 25 | `mean_internal_bl_locus_loss_only` | float | Same, restricted to `loss_only`. **All NaN by design. see 2.1** |
| 26–30 | `mean_internal_bl_locus_{dup_and_loss, nothing, false_HP, weak_HP, strong_HP}` | float | Same locus-tree branch-length statistic restricted to each HP bucket. NaN when DUP/LOS = 0 (except `_nothing` which is always populated). |
| 31 | `mean_internal_bl_gene_all` | float | Same metric but on the *estimated* gene tree (post IQ-TREE). |
| 32 | `mean_internal_bl_gene_loss_only` | float | **All NaN by design.** |
| 33–37 | `mean_internal_bl_gene_{dup_and_loss, nothing, false_HP, weak_HP, strong_HP}` | float | Estimated-gene-tree branch lengths per HP bucket. NaN when DUP/LOS = 0 (except `_nothing`). |
| 38 | `mean_RF_genetree_vs_sptree_all` | float | Mean RF between **true** locus tree and the true species tree, across all genes. Measures gene-tree discordance due to ILS + paralogy filtering. |
| 39 | `mean_RF_genetree_vs_sptree_loss_only` | float | **All NaN by design.** |
| 40–44 | `mean_RF_genetree_vs_sptree_{dup_and_loss, nothing, false_HP, weak_HP, strong_HP}` | float | Same RF restricted to each HP bucket. NaN when DUP/LOS = 0 (except `_nothing`). |
| 45 | `mean_RF_genetree_vs_sptree_noF_all` | float | Like column 38, but computed *without filtering* gene trees by the paralogy criterion (upper-bound discordance). |
| 46 | `mean_RF_genetree_vs_sptree_noF_loss_only` | float | **All NaN by design.** |
| 47–51 | `mean_RF_genetree_vs_sptree_noF_{dup_and_loss, nothing, false_HP, weak_HP, strong_HP}` | float | As 40–44 but no filtering. |
| 52 | `mean_RF_true_vs_est_genetree_all` | float | Mean RF between **true** locus tree and **estimated** gene tree (IQ-TREE). Measures gene-tree estimation error. |
| 53 | `mean_RF_true_vs_est_genetree_loss_only` | float | **All NaN by design.** |
| 54–58 | `mean_RF_true_vs_est_genetree_{dup_and_loss, nothing, false_HP, weak_HP, strong_HP}` | float | Same estimation-error metric per HP bucket. NaN when DUP/LOS = 0 (except `_nothing`). |

### 2.5 Derived tables used to build the paper's combined figure

These three CSVs are pre-computed pivots of the three tables above. They are
the direct inputs to `visualization_scripts/visual_combined.qmd`, which
produces `plots/combined_three_panel_figure.{pdf,png}` (Figure 1 of the paper).

#### `combined_graph_recovery_summary.csv` (36 rows)

Panel B of Figure 1: % of replicates that recovered the true species tree.

| # | Column | Description |
|--:|---|---|
| 1 | `dup_loss_rate` | Numeric duplication/loss rate (0, 3e-4, or 4e-4). |
| 2 | `RV` | Rate-variation level (N, G, L). |
| 3 | `N_ind` | Individuals per species (1 or 2). |
| 4 | `ILS` | Categorical ILS level: `"high (SF=1.0)"` or `"low (SF=0.5)"`. |
| 5 | `snaq_rec_all` | % of replicates where SNaQ's best network recovers the true species tree as its backbone. |
| 6 | `snaq_rec_noF` | Same, without taxon F. |
| 7 | `fg_rec_alltrees` | % of replicates where find_graphs recovers the true tree among **all** returned trees (not just the best). |
| 8 | `fg_rec_alltrees_noF` | Same, without taxon F. |
| 9 | `fg_rec_best` | % of replicates where find_graphs' **best** tree equals the true species tree. |
| 10 | `fg_rec_best_noF` | Same, without taxon F. |

#### `combined_hypothesis_acceptance_summary.csv` (36 rows)

Panel C of Figure 1 (per-setting view): hypothesis acceptance rates under
both methods. `wr30` = standard WR ≤ 3.0; `wr37` = relaxed WR ≤ 3.7.

| # | Column | Description |
|--:|---|---|
| 1–4 | `dup_loss_rate`, `RV`, `N_ind`, `ILS` | Setting identifiers (same as above). |
| 5 | `snaq_H0` | % of replicates accepting H = 0 under SNaQ. |
| 6 | `snaq_H1` | % accepting H = 1 under SNaQ. |
| 7 | `snaq_Hgt1` | % accepting H > 1 under SNaQ. |
| 8 | `fg_H0_wr30` | % accepting k = 0 admixture events under find_graphs (WR ≤ 3.0). |
| 9 | `fg_H1_wr30` | % accepting k = 1 (WR ≤ 3.0). |
| 10 | `fg_Hgt1_wr30` | % accepting k > 1 (WR ≤ 3.0). |
| 11–13 | `fg_*_wr37` | Same find_graphs columns at the relaxed WR ≤ 3.7 threshold. |

#### `combined_hypothesis_acceptance_marginal.csv` (10 rows)

Marginal view of the table above: each row is the **mean** over the subset
of settings where the named factor takes the named level. Used for Panel A
and the text statistics in the paper.

| # | Column | Description |
|--:|---|---|
| 1 | `group` | Which factor is marginalized over (`dup_loss_rate`, `RV`, `N_ind`, `ILS`). |
| 2 | `Level` | Level of that factor (e.g., `"Dup/Loss = 0"`, `"RVL"`, `"high (SF=1.0)"`). |
| 3 | `snaq_type1` | Mean **type I error** of SNaQ at this level: % of replicates falsely accepting H ≥ 1 when no reticulation exists. Identical to `snaq_H1 + snaq_Hgt1` because the true history is a tree. |
| 4 | `snaq_H1` | Mean % accepting exactly H = 1 (SNaQ). |
| 5 | `snaq_Hgt1` | Mean % accepting H > 1 (SNaQ). |
| 6 | `fg_type1_wr30` | Mean type I error of find_graphs at WR ≤ 3.0 (`fg_H1_wr30 + fg_Hgt1_wr30`). |
| 7 | `fg_H1_wr30` | Mean % accepting k = 1 (find_graphs, WR ≤ 3.0). |
| 8 | `fg_Hgt1_wr30` | Mean % accepting k > 1 (find_graphs, WR ≤ 3.0). |
| 9–11 | `fg_*_wr37` | Same find_graphs columns at WR ≤ 3.7. |

### 2.6 Per-setting per-replicate inputs — `findgraph_summary/` and `snaq_summary/`

These two subfolders hold the **per-replicate** CSVs that the cross-setting
summary tables in §2.2–§2.3 (and several paper figures) are computed from.
Each subfolder contains **36 CSV files — one per parameter setting**, with
**one row per replicate (typically 100 rows per file)**.

**Note on the `_noF` / `_noG` / `_noH` column suffixes.** Several columns in
both files come in variants suffixed with `_noF`, `_noG`, or `_noH`. These
re-evaluate the same comparison (true-tree recovery, RF distance, etc.)
after **pruning that single taxon — F, G, or H respectively — from both the
true species tree and the inferred tree / network before the comparison**.
They serve as a robustness check: a "true tree recovered" hit under the
unfiltered comparison can be lost just because one taxon is misplaced, and
the `_noF` / `_noG` / `_noH` variants reveal how much of the apparent
non-recovery is driven by a single hard-to-place taxon. Columns *without*
any of these suffixes use the full taxon set.

In this schema `k` denotes the number of admixture / reticulation edges and
is interchangeable with `h` in the cross-setting tables — `k` follows the
`find_graphs` paper's notation; `h` follows the SNaQ / network-inference
literature.

#### `findgraph_summary/findgraph-<paramname>.csv` — per-replicate find_graphs

One row per replicate (typically 100 per file).

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `repID` | string | Three-digit replicate identifier (`001` … `100`). |
| 2 | `best_k` | int / `">1"` | Best `k` (num of reticulatio, = h) accepted under the standard WR ≤ 3.0 threshold: `0`, `1`, or `">1"`. |
| 3 | `H0_trees_found` | int | Number of distinct `k = 0` (tree) candidates returned by `find_graphs` for this replicate (after deduplication). |
| 4 | `H1_graphs_found` | int | Number of distinct `k = 1` (one-reticulation) candidates returned. |
| 5 | `true_tree_wr` | float | Worst-residual value when the **true species tree** itself is scored — used to draw the true-tree WR distribution / percentile plots. |
| 6–9 | `H0_best_tree_is_true_tree`, `..._noF`, `..._noG`, `..._noH` | bool | Whether the highest-likelihood `k = 0` tree equals the true species tree (full taxon set, then pruning F/G/H — see note above). |
| 10 | `H0_best_tree_WR` | float | Worst-residual of the highest-likelihood `k = 0` tree. |
| 11–14 | `H1_best_graph_displayed_true_tree`, `..._noF`, `..._noG`, `..._noH` | bool | Whether the highest-likelihood `k = 1` graph displays the true species tree as its major or minor backbone (full / pruned). |
| 15 | `H1_best_graph_WR` | float | Worst-residual of the highest-likelihood `k = 1` graph. |
| 16–19 | `true_tree_found_in_H0`, `..._noF`, `..._noG`, `..._noH` | bool | Whether **any** `k = 0` tree returned for this replicate equals the true species tree (full / pruned). |
| 20 | `avg_gamma1_H1` | float ∈ [0, 1] | Per-replicate mean major-edge weight γ₁ across all returned `k = 1` graphs. |
| 21 | `avg_gamma2_H1` | float ∈ [0, 1] | Per-replicate mean minor-edge weight γ₂ across all returned `k = 1` graphs. |
| 22 | `best_graph_gamma1` | float ∈ [0, 1] | γ₁ of the single highest-likelihood `k = 1` graph (drives the γ-distribution plots). |
| 23 | `best_graph_gamma2` | float ∈ [0, 1] | γ₂ of the single highest-likelihood `k = 1` graph. |
| 24–31 | `True_tree_displayed_H1_major/minor`, `True_tree_noF_..._major/minor`, `..._noG_...`, `..._noH_...` | bool | Pair of indicators per taxon-pruning variant: whether the true tree is displayed as the **major** vs. **minor** backbone of the best `k = 1` graph. |
| 32 | `best_k_new_WR_3.7` | int / `">1"` | Same as `best_k`, but re-classified under the **relaxed** WR ≤ 3.7 threshold. |
| 33 | `best_graph_hybrid_taxon` | string / `"NA"` | Name of the hybrid (admixed) taxon in the best `k = 1` graph. Legacy columns, not used |
| 34 | `best_graph_major_donor` | string / `"NA"` | Major-donor parent taxon in the best `k = 1` graph.  Legacy columns, not used|
| 35 | `best_graph_minor_donor` | string / `"NA"` | Minor-donor parent taxon in the best `k = 1` graph.  Legacy columns, not used|

Columns 33–35 are used only by the legacy taxon-recovery analysis, which
is **disabled in the published pipeline** (see the LEGACY block in
`scripts/summary_findgraph.jl`). They are retained for schema stability.

#### `snaq_summary/SNaQ-<paramname>-summary.csv` — per-replicate SNaQ

One row per replicate (typically 100 per file).

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `repID` | string | Three-digit replicate identifier. |
| 2 | `p_H0` | float ∈ [0, 1] | Goodness-of-fit p-value under H = 0. H = 0 is accepted when `p_H0 > 0.05`. |
| 3 | `p_H1` | float ∈ [0, 1] | Goodness-of-fit p-value under H = 1. H = 1 is accepted when `p_H0 ≤ 0.05` *and* `p_H1 > 0.05`. |
| 4 | `z_uncorrected_H0` | float | Uncorrected goodness-of-fit z-statistic under H = 0 (used to compute `p_H0`). |
| 5 | `z_uncorrected_H1` | float | Same under H = 1. |
| 6 | `sigma_H0` | float | Standard error used to scale the H = 0 z-statistic. |
| 7 | `sigma_H1` | float | Same under H = 1. |
| 8 | `score_H0` | float | SNaQ pseudo-deviance score for the best H = 0 network (lower is better). |
| 9 | `score_H1` | float | Same for the best H = 1 network. |
| 10 | `gamma_1` | float ∈ [0, 1] | Major-edge weight γ₁ of the best H = 1 network. |
| 11 | `gamma_2` | float ∈ [0, 1] | Minor-edge weight γ₂ = 1 − γ₁. |
| 12 | `RF_net0_true` | float | Robinson-Foulds distance between the best H = 0 network and the true species tree (full taxon set). `0.0` means exact match. |
| 13 | `RF_net0_true_noF` | float | Same after pruning taxon F from both trees. |
| 14 | `RF_net1_1_true` | float | RF distance between the **major** backbone of the best H = 1 network and the true species tree (full taxon set). |
| 15 | `RF_net1_2_true` | float | RF distance between the **minor** backbone of the best H = 1 network and the true species tree. |
| 16–17 | `RF_net1_1_true_noF`, `RF_net1_2_true_noF` | float | Same after pruning taxon F. |
| 18–19 | `RF_net1_1_true_noG`, `RF_net1_2_true_noG` | float | Same after pruning taxon G. |
| 20–21 | `RF_net1_1_true_noH`, `RF_net1_2_true_noH` | float | Same after pruning taxon H. |
| 22–24 | `RF_net0_alter1`, `..._alter2`, `..._alter3` | float | RF distances from the best H = 0 network to the three pre-specified non-true alternative topologies (defined in `scripts/utilities.jl`). |
| 25–27 | `RF_net1_1_alter1/2/3` | float | RF distances from the H = 1 **major** backbone to each alternative topology. |
| 28–30 | `RF_net1_2_alter1/2/3` | float | RF distances from the H = 1 **minor** backbone to each alternative topology. |
| 31 | `hybrid_taxon` | string / `"NA"` | Hybrid taxon in the best H = 1 network (used only by the legacy taxon-recovery analysis, **disabled in publication** — see comments in `scripts/summary_snaq.jl`). |
| 32 | `major_donor` | string / `"NA"` | Major-donor parent in the best H = 1 network (legacy). |
| 33 | `minor_donor` | string / `"NA"` | Minor-donor parent in the best H = 1 network (legacy). |

**Why these are included in the folder.** The cross-setting tables in §2.2–§2.5
only carry *means*, *medians*, and *counts* across replicates — they do not
preserve the per-replicate values needed to redraw the distributions
of worst-residual (WR) of of the minor γ (admixture weights / inheritance prob).
Those plots are produced by `scripts/summary_findgraph.jl` and
`scripts/summary_snaq.jl`, which read these subfolders directly:

- `scripts/summary_findgraph.jl` defaults to `--input_dir findgraph_summary`
  and uses its files to compute the cross-setting `findgraph_summary.csv`
  *and* to drive the R routines in `scripts/visual_utilities.R` that produce
  the gamma-overlap plots (`plot_overlapping_ratevar_by_n_inds_sf`, reading
  `best_graph_gamma1` / `best_graph_gamma2`), the true-tree WR distribution
  and percentile plots (`plot_WR_distributions`,
  `plot_WR_percentiles_by_rate_variation`,
  `plot_WR_percentiles_jitter_combined`, all reading `true_tree_wr`), and
  the per-replicate statistic-by-tree-display plots (reading
  `best_graph_gamma2`, `H1_best_graph_WR`, and
  `H1_best_graph_displayed_true_tree`).
- `scripts/summary_snaq.jl` defaults to `--input_dir snaq_summary` and uses
  its files to compute `SNaQ_summary.csv`, the filtered
  `SNaQ_minor_gamma_filtered.csv` (replicates where neither H1 backbone
  matches the true tree), and the corresponding γ-distribution plots
  (`plot_overlapping_ratevar_by_n_inds_sf` reading `gamma_1` / `gamma_2`,
  and `plot_snaq_minor_gamma_by_tree_display`).

**Where these files come from.** The two folders are produced
upstream by `scripts/findgraphs_postprocess.jl` and `scripts/snaq_postprocess.jl`
(which collect each replicate's raw find_graphs / SNaQ output into one CSV
per setting) and then assembled into the central folders by
`scripts/run_postprocessing.jl`, whose `--saved_path` defaults to
`<project_root>/findgraph_summary` and `<project_root>/snaq_summary`
respectively.
**In the repository, these two folders have been relocated from the project root
into `results`**, as `results/findgraph_summary` and `results/snaq_summary`,
to share key results files, those of reasonable size.
If we re-run the upstream postprocessing scripts,
they will write their output to the project root by default
(and will not overwrite the files in this `results/` folder).
To keep paths consistent,
we may either move the folders to `results/` afterwards, or pass
`--saved_path results/<mode>_summary` and `--input_dir results/<mode>_summary`
to the summary scripts.

---

## 3. Code and software

The [GitHub](https://github.com/cecileane/simulation-hiddenparalogy-ratevariation)
repository (also available on [Dryad](https://doi.org/10.5061/dryad.dfn2z35h7))
contains the full simulation pipeline, the per-setting and
cross-setting summary scripts that produced the CSVs in this folder, and the
Quarto notebooks that produce figures.

**Pipeline summary** (from `readme.md` in the repository):

1. `scripts/simulation.jl`         — SimPhy → paralogy filter → Seq-Gen → IQ-TREE → ASTRAL-IV (per replicate, parallel)
2. `scripts/snaq.jl`               — SNaQ at H = 0 and H = 1 (per replicate)
3. `scripts/findgraphs.jl`         — find_graphs / qpgraph at k = 0 and k = 1 (per replicate)
4. `scripts/run_postprocessing.jl` — drives `scripts/snaq_postprocess.jl` and
   `scripts/findgraphs_postprocess.jl` over every parameter set and collects
   their per-setting per-replicate CSVs into `snaq_summary/` and
   `findgraph_summary/` (the two subfolders described in §2.6;
   in this folder they live under `results/`)
5. `scripts/summary_simulation.jl`, `scripts/summary_snaq.jl`,
   `scripts/summary_findgraph.jl` — read `snaq_summary/` and `findgraph_summary/`
   (along with the simulation outputs) to produce both the cross-setting tables
   `summary_concatenated.csv`, `SNaQ_summary.csv`, `findgraph_summary.csv`
   and the WR / γ distribution plots used in the paper (see §2.6)
6. `visualization_scripts/visual_combined.qmd` — produces the paper's combined
   figure from the `combined_*.csv` files

**Software versions** used to produce these data (Linux x86-64). External
binaries are listed first; Julia packages — pinned via `Project.toml` /
`Manifest.toml` in the repository — and the R package follow.

| Software | Version | Role in the pipeline |
|---|---|---|
| Julia | 1.11.7 | Driver language for the simulation pipeline and the postprocessing / summary scripts. |
| SimPhy | 1.0.2 | Simulates species and locus trees under the multispecies coalescent with duplication / loss. |
| Seq-Gen | 1.3.5 | Simulates sequence evolution along the locus trees. |
| IQ-TREE | 2.4.0 | Estimates per-locus gene trees from the simulated sequences. |
| ASTER (ASTRAL-IV) | 1.24.4.8 | Estimates the species tree from the per-locus gene trees. |
| snp-sites | 2.5.1 | Extracts variant sites from each locus alignment for `find_graphs` / qpgraph input. |
| Python | 3.12 | Used by helper scripts in `scripts/` (e.g. `concatenate_seq.py`). |
| R | system R | Drives `admixtools` for `find_graphs` / qpgraph, and the visualization routines in `scripts/visual_utilities.R`. |
| R `admixtools` | 2.0.8 (igraph pinned to 1.6.0) | Provides `find_graphs` / qpgraph, the network-inference method tested in this study (h = 0 and h = 1). |
| `PhyloNetworks.jl` | 1.3.1 | Julia network/tree data structures and I/O used throughout. SNaQ depends on it. |
| `SNaQ.jl` | 1.1.1 | The pseudolikelihood-based network-inference method tested in this study (H = 0 and H = 1). |
| `QuartetNetworkGoodnessFit.jl` | 1.0.0 | Goodness-of-fit test for SNaQ networks; produces the `p_H0`, `p_H1`, and `z_uncorrected_*` values in `snaq_summary/`. |
| `PhyloPlots.jl` | 2.1.0 | Plotting of inferred / consensus networks. |
| `RCall.jl` | 0.14.12 | Bridges Julia to R for the visualization layer and for `admixtools`. |
| `CSV.jl` / `DataFrames.jl` | 0.10.16 / 1.8.1 | Tabular I/O and manipulation for every summary CSV in this folder. |
| `Distributions.jl` | 0.25.123 | Distributions used during simulation. |
| `StableRNGs.jl` | 1.0.4 | Reproducible random-number generation (the deterministic seed scheme). |
| `StatsBase.jl` / `StatsPlots.jl` / `Plots.jl` / `CairoMakie.jl` | 0.34.10 / 0.15.8 / 1.41.6 / 0.15.9 | Summary statistics and the Julia-side plotting. |
| `KernelDensity.jl` | 0.6.11 | Kernel-density estimates for the γ and WR distribution diagnostics. |
| `ArgParse.jl` / `Glob.jl` / `IterTools.jl` / `Primes.jl` / `PrettyTables.jl` / `TimerOutputs.jl` / `TimeZones.jl` / `InlineStrings.jl` | 1.2.0 / 1.4.0 / 1.10.0 / 0.5.7 / 3.3.1 / 0.5.29 / 1.22.2 / 1.4.5 | Auxiliary utilities used by the scripts. |

The Julia versions above are the resolved versions recorded in
`Manifest.toml` in the repository; `Project.toml` lists
the direct dependencies. To reproduce the exact environment, instantiate
the project with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`,
which honors `Manifest.toml`.

Per-pipeline-stage runtime, dependency versions, and the deterministic seed
scheme are documented in `scripts/readme.md` and `plots/seed_control.png`.

---

## 4. Access information

**License**: CC0 (per Dryad policy).
