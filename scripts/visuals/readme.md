# Figures and tables

R scripts for the figures and tables of the paper, with the reptile and the
angiosperm tree in the same figure. 
They read the per-replicate summaries
in `results/summary_reptile/<method>/` and `results/summary_plant/<method>/`
(written by scripts/run_postprocessing.jl) 
and write to
`visualization_results/`. 
The tables are then copied to `results/tables/`.
Run from the repo root:

```bash
Rscript scripts/visuals/visual_main_figures.R
Rscript scripts/visuals/visual_appendix_figures.R
Rscript scripts/visuals/visual_tables.R
```

- `visual_main_figures.R`: Fig. 1 (four panels, one point per setting) and
  Fig. 2 (mirrored histograms of the minor gamma), plus per_setting.csv
  with the numbers behind Fig. 1 and T_star.csv, 
  in `visualization_results/main_figures/`
- `visual_appendix_figures.R`: the supplementary figures, one folder each
  (A gene tree discordance, B worst residual and T of every replicate,
  D placement of the PhyNEST reticulation, E model choice, F true tree
  recovery, G minor gamma, H best find_graphs graph), 
  in `visualization_results/appendix_two_trees/`
- `visual_tables.R`: the supplementary tables, each with a reptile and a
  plant column under every statistic: 
  worst residual of the true tree (find_graphs), per-site T (PhyNEST), model choice, gene tree discordance
  and estimation error, minor gamma of each method, best find_graphs graph.
  csv and tex in `visualization_results/tables/`
- `visual_utilities.R`: sourced by the three others: colors and labels,
  dup/loss rate to hidden paralogy level, 
  the loader of the summaries,
  T* and the y axis rule of the mirrored histograms

In all figures x is the hidden paralogy level (none, medium, high: dup/loss
0, 3e-4, 4e-4 for the reptiles and 0, 1e-4, 
1.4e-4 for the plants), 
each level split into a reptile and a plant sub-column.

The two gene tree tables read `results/summary_<tree>/genetrees/`, the branch
counts of scripts/genetree_tiers.jl. T* of PhyNEST (95th percentile of the
per-site T over the settings without lineage rate variation, 
per tree and ILS level) is recomputed from the summaries at every run and printed.
