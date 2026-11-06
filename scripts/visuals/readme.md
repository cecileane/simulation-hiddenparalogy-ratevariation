# Visualization scripts

R scripts for the figures for both reptile and angiosperm
species tree in the same figure. 
They read the per-replicate summaries of in summary_reptile/<method>/ and
summary_plant/<method>/ (written by scripts/run_postprocessing). 

| Script | Content | Output |
|---|---|---|
| visual_main_figures.R | Fig. 1 (four panels, one point per setting) and Fig. 2 (mirrored histograms of the minor gamma); also per_setting.csv and T_star.csv | visualization_results/main_figures/ |
| visual_appendix_figures.R | supplementary figures A, B, D, E, F, G, one folder each (see the README there) | visualization_results/appendix_two_trees/ |
| visual_tables.R | supplementary tables with a reptile and a plant subcolumn under every statistic: worst residual of the true tree (find_graphs), per-site T (PhyNEST), model choice, gene tree discordance (locus category x ILS), gene tree estimation error (rate variation, one tabular per number of individuals) and the minor gamma under h = 1 | visualization_results/tables/ |
| visual_utilities.R | sourced by all three: colors and labels, dup/loss rate to hidden-paralogy level, and utility functions | – |

```bash
Rscript scripts/visuals/visual_main_figures.R
Rscript scripts/visuals/visual_appendix_figures.R
Rscript scripts/visuals/visual_tables.R
```

Figure design: x = hidden paralogy level (none / medium / high; dup/loss
0, 3e-4, 4e-4 for the reptiles and 0, 1e-4, 1.4e-4 for the plants), each level
split into a reptile and a plant sub-column. 

The two gene tree tables read summary_<tree>/genetrees/, the per-setting
branch counts of scripts/genetree_tiers.jl. 

T* of PhyNEST is the 95th percentile of the per-site T pooled over the
settings without lineage rate variation, per tree and ILS level. It is
recomputed from the data in use at every run and printed by both scripts.
