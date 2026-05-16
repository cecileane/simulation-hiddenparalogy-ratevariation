# Data from: Substitution rate variation, not hidden paralogy, drives false hybridization signal in phylogenetic network inference

This Dryad deposit accompanies the manuscript submitted to *Systematic Biology*
(USYB-2026-124). It contains the cross-setting summary tables that underlie
every figure and statistic in the paper. The code used to produce these tables
from raw simulations is on GitHub (link below).

---

## 1. Data description

We simulated phylogenomic datasets to test whether two common model violations
— gene duplication/loss (which can produce *hidden paralogy*) and substitution
rate variation — can cause network inference methods (`find_graphs` and SNaQ)
to detect spurious hybridization when the true history is a tree. Branch
lengths and substitution-rate parameters were calibrated from empirical reptile
UCE data ([Crawford et al. 2012](https://academic.oup.com/sysbio/article/61/5/717/1735316)).

**Experimental design — 36 parameter settings, 100 replicates each:**

| Factor | Values | Meaning |
|---|---|---|
| `DUP` / `LOS` | 0.0, 3 × 10⁻⁴, 4 × 10⁻⁴ | Gene duplication / loss rate (per gene, per generation). Duplication and loss rates are equal in every setting. `0.0` disables both, so no hidden paralogy can arise. |
| `RV` (rate variation) | N, G, L | None / Gene-specific / Lineage-specific substitution-rate variation |
| `N_ind` | 1, 2 | Individuals sampled per species |
| `SF` (scale factor) | 0.5, 1.0 | Effective-population-size scaling factor; higher = more incomplete lineage sorting |
| `genelen` | 1000 | Gene length, in base pairs (constant across all settings) |

Each setting is encoded in a string of the form
`DUP<d>-LOS<l>-RV<r>-N_ind<n>-SF<s>-genelen<g>` (e.g., `DUP0.0003-LOS0.0003-RVG-N_ind1-SF0.5-genelen1000`). This string is the row identifier
in every CSV in this deposit.

Per setting we ran 100 replicates of the full pipeline:
SimPhy → paralogy filter → Seq-Gen → IQ-TREE → ASTRAL-IV → SNaQ + find_graphs.

---

## 2. Files and variables

This deposit contains 6 CSV files. Three are the primary summary tables; three
are pre-computed cross-tabulations used to make the paper's combined figure.

All files have **one row per parameter setting (36 rows)**, *except*
`combined_hypothesis_acceptance_marginal.csv` which has 10 rows (one per
factor level). Each file uses comma separators and a single header row.

### 2.1 Missing-data convention

`NaN` (Not a Number) is used throughout to indicate **"this metric is not
defined under this parameter combination"**. It is *never* used to indicate
a measurement failure. There are two reasons it appears:

- **`*_dup_and_loss`, `*_false_HP`, `*_weak_HP`, `*_strong_HP` columns are NaN
  in the 12 rows where `DUP/LOS = 0.0`.** With no duplication or loss events,
  no gene tree can fall into the dup-and-loss bucket or any of the hidden-
  paralogy buckets. The category is structurally empty.
- **The six `*_loss_only` columns are NaN in *all 36 rows*.** These columns
  are legacy QA placeholders from an earlier design that allowed missing taxa
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
| 3 | `H0Accepted` | int | Replicates where the best find_graphs model chosen by the *standard* worst-residual (WR) threshold (WR ≤ 3.0) had k = 0 admixture events. |
| 4 | `H1Accepted` | int | As above, but model chosen had k = 1. |
| 5 | `BT1Accepted` | int | As above, but model chosen had k > 1 ("BT1" = "best k beyond 1"). `H0Accepted + H1Accepted + BT1Accepted = total_replicates`. |
| 6 | `pct_true_tree_H0` | float (%) | % of replicates whose best H = 0 tree equals the true species tree (Robinson-Foulds distance = 0). |
| 7 | `avg_num_blocks` | float | Mean number of SNP blocks used by find_graphs across replicates (target: 1000). |
| 8 | `avg_gamma1_H1` | float ∈ [0, 1] | Mean major edge weight (γ₁) of the best H = 1 graph across replicates. |
| 9 | `avg_gamma2_H1` | float ∈ [0, 1] | Mean minor edge weight (γ₂ = 1 − γ₁) of the best H = 1 graph. |
| 10 | `avg_best_gamma1` | float | Mean γ₁ of the single best-likelihood H = 1 graph per replicate. |
| 11 | `avg_best_gamma2` | float | Mean γ₂ of the single best-likelihood H = 1 graph per replicate. |
| 12 | `pct_true_tree_H1` | float (%) | % of replicates where any H = 1 graph displays the true species tree as its major or minor backbone. |
| 13 | `pct_true_tree_H1_noF` | float (%) | Same as column 12 but counted *without the WR-filtering step* applied to per-replicate graphs (the "noF" suffix = "no filter"). Provides an upper bound on tree recoverability before model selection. |
| 14 | `mean_H0_trees_found` | float | Mean number of distinct H = 0 trees found per replicate after deduplication. |
| 15 | `median_H0_trees_found` | float | Median (per replicate). |
| 16 | `sd_H0_trees_found` | float | Standard deviation (per replicate). |
| 17 | `min_H0_trees_found` | int | Minimum across replicates. |
| 18 | `max_H0_trees_found` | int | Maximum across replicates. |
| 19–23 | `mean/median/sd/min/max_H1_graphs_found` | numeric | Same statistics for H = 1 graphs. |
| 24 | `count_H0_best_is_true_tree` | int | Replicates where the highest-likelihood H = 0 tree equals the true species tree (RF = 0). |
| 25 | `count_H0_best_is_true_tree_noF` | int | Same, without WR filtering. |
| 26 | `count_H1_best_displays_true_tree` | int | Replicates where the highest-likelihood H = 1 graph displays the true species tree. |
| 27 | `count_H1_best_displays_true_tree_noF` | int | Same, without WR filtering. |
| 28 | `H0Accepted_WR_3.7` | int | Like column 3, but using the *relaxed* WR threshold (WR ≤ 3.7). |
| 29 | `H1Accepted_WR_3.7` | int | As column 4, relaxed threshold. |
| 30 | `BT1Accepted_WR_3.7` | int | As column 5, relaxed threshold. |

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
| 10 | `find_true_net0_noF` | int | Same, without filtering (upper bound). |
| 11 | `find_true_net1` | int | Replicates where the best H = 1 network displays the true tree as either its major or minor backbone. |
| 12 | `find_true_net1_noF` | int | Same, without filtering. |
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

- `loss_only`     — gene trees that experienced loss but **no** duplication. Always empty in this design (see §2.1).
- `dup_and_loss`  — gene trees that experienced both duplication and loss. Always empty when DUP/LOS = 0.
- `nothing`       — gene trees with neither duplication nor loss.
- `false_HP`      — duplicated copies present but all duplicates were sampled and labeled (paralogy is *not* hidden).
- `weak_HP`       — some duplicates lost, but enough copies remain that paralogy is detectable.
- `strong_HP`     — all but one paralog lost, producing classic *hidden* paralogy (cannot be detected from the gene tree alone).

| # | Column | Type | Description |
|--:|---|---|---|
| 1 | `parameter_setting` | string | Parameter-setting identifier. |
| 2 | `n_genes_mean` | float | Mean number of single-copy gene trees retained per replicate (target = 1000). |
| 3 | `n_genes_min` | int | Minimum per replicate. |
| 4 | `n_genes_max` | int | Maximum per replicate. |
| 5 | `n_iterations` | float | Mean number of SimPhy re-draws required to reach `n_genes_min` (averaged across replicates). |
| 6 | `n_repeated_taxa_removed` | float | Mean number of gene trees dropped per replicate due to within-species paralogs (the species had ≥ 2 copies of the gene). |
| 7 | `n_insufficient_taxa_removed` | float | Mean number of gene trees dropped per replicate for having too few remaining taxa to be phylogenetically informative. |
| 8 | `percentage_genes_meet_min` | float (%) | % of replicates that met the minimum-genes threshold without exhausting the SimPhy re-draw quota. |
| 9 | `pert_trees_experiencing_gene_loss_only` | float (%) | % of gene trees in the `loss_only` bucket (see HP categories above). Note: `pert_` is a typographical artifact in the column schema; read as "percent". |
| 10 | `pert_trees_experiencing_gene_duplication_and_loss` | float (%) | % in the `dup_and_loss` bucket. |
| 11 | `pert_trees_experiencing_nothing` | float (%) | % in the `nothing` bucket. |
| 12 | `pert_false_HP` | float (%) | % in the `false_HP` bucket. |
| 13 | `pert_weak_HP` | float (%) | % in the `weak_HP` bucket. |
| 14 | `pert_strong_HP` | float (%) | % in the `strong_HP` bucket. |
| 15 | `mean_RF_true_and_estimated_species_trees` | float | Mean Robinson-Foulds distance between the true species tree and the ASTRAL-IV estimate, across replicates. |
| 16 | `num_estimated_species_tree_diff_from_truth` | int | Replicates whose ASTRAL estimate differs from the true species tree (RF > 0). |
| 17 | `mean_num_taxa_all_genes` | float | Mean number of taxa per gene tree, across all gene trees in all replicates. |
| 18 | `mean_num_taxa_loss_only` | float | Same, restricted to `loss_only` trees. **All NaN by design — see §2.1.** |
| 19 | `mean_num_taxa_dup_and_loss` | float | Same, restricted to `dup_and_loss` trees. NaN when DUP/LOS = 0. |
| 20 | `mean_num_taxa_nothing` | float | Same, `nothing` trees. |
| 21 | `mean_num_taxa_false_HP` | float | Same, `false_HP`. NaN when DUP/LOS = 0. |
| 22 | `mean_num_taxa_weak_HP` | float | Same, `weak_HP`. NaN when DUP/LOS = 0. |
| 23 | `mean_num_taxa_strong_HP` | float | Same, `strong_HP`. NaN when DUP/LOS = 0. |
| 24 | `mean_internal_bl_locus_all` | float | Mean of average internal branch length on the *locus* tree (the true SimPhy gene tree before sequence simulation), across all genes. Units: SimPhy's default branch-length units (substitutions per site). |
| 25 | `mean_internal_bl_locus_loss_only` | float | Same, restricted to `loss_only`. **All NaN by design.** |
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
| 5 | `snaq_rec_all` | % of replicates where SNaQ's best network (any H) recovers the true species tree as its backbone. |
| 6 | `snaq_rec_noF` | Same, without filtering. |
| 7 | `fg_rec_alltrees` | % of replicates where find_graphs recovers the true tree among **all** returned trees (not just the best). |
| 8 | `fg_rec_alltrees_noF` | Same, without filtering. |
| 9 | `fg_rec_best` | % of replicates where find_graphs' **best** tree equals the true species tree. |
| 10 | `fg_rec_best_noF` | Same, without filtering. |

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

---

## 3. Code and software

All code is available on GitHub: **<URL of tagged release — fill in on Dryad submission>**.

The repository contains the full simulation pipeline, the per-setting and
cross-setting summary scripts that produced the CSVs in this deposit, and the
Quarto notebooks that produce the paper figures.

**Pipeline summary** (from `readme.md` in the repository):

1. `scripts/simulation.jl`         — SimPhy → paralogy filter → Seq-Gen → IQ-TREE → ASTRAL-IV (per replicate, parallel)
2. `scripts/snaq.jl`               — SNaQ at H = 0 and H = 1 (per replicate)
3. `scripts/findgraphs.jl`         — find_graphs / qpgraph at k = 0 and k = 1 (per replicate)
4. `scripts/run_postprocessing.jl` — aggregate per-replicate results
5. `scripts/summary_simulation.jl`, `scripts/summary_snaq.jl`, `scripts/summary_findgraph.jl` — produce `summary_concatenated.csv`, `SNaQ_summary.csv`, `findgraph_summary.csv` (the three primary files in this deposit)
6. `visualization_scripts/visual_combined.qmd` — produces the paper's combined figure from the `combined_*.csv` files

**Software versions** used to produce these data (Linux x86-64):

| Software | Version |
|---|---|
| Julia | ≥ 1.11 |
| SimPhy | 1.0.2 |
| Seq-Gen | 1.3.5 |
| IQ-TREE | 2.4.0 |
| ASTER (ASTRAL-IV) | 1.24.4.8 |
| snp-sites | 2.5.1 |
| R | system R |
| R `admixtools` | 2.0.8 (igraph pinned to 1.6.0 — see `executables/README.md`) |
| Python | 3.12 |

Per-pipeline-stage runtime, dependency versions, and the deterministic seed
scheme are documented in `scripts/readme.md` and `plots/seed_control.png`.

---

## 4. Access information

**License**: CC0 (per Dryad policy).

**Recommended citation**:
> Li, B., Solís-Lemus, C., & Ané, C. (2026). *Data from: Substitution rate
> variation, not hidden paralogy, drives false hybridization signal in
> phylogenetic network inference.* Dryad Digital Repository. <DOI on
> acceptance>.

**Corresponding author**: Cécile Ané (UW–Madison) — see manuscript title page.

**Code corresponding author** (for questions about the pipeline or CSV
schemas): Bing Li — bingli8899@gmail.com.

**Linked manuscript**: *Systematic Biology* USYB-2026-124, "Substitution rate
variation, not hidden paralogy, drives false hybridization signal in
phylogenetic network inference."
