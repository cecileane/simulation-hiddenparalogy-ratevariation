# Visualization Scripts

Quarto notebooks (`.qmd`) that generate figures and summary tables from the aggregated CSV outputs produced by the summary scripts in `scripts/`. Run these after `summary_simulation.jl`, `summary_snaq.jl`, and `summary_findgraph.jl` have completed.

## Notebooks

| Notebook | Content |
|---|---|
| `visualization_simulation.qmd` | Gene tree discordance, RF distances, ASTRAL accuracy |
| `visualization_snaq.qmd` | SNaQ topology recovery and model selection rates |
| `visualization_finsgraphs.qmd` | find_graphs WR distributions, gamma estimates, type I error |
| `visual_combined.qmd` | Side-by-side comparison of both methods across all conditions |
| `visualize_baselinetree.qmd` | Baseline (no model violation) results |

## How to render

Open any `.qmd` in VS Code (Quarto extension) and click **Render**, or run from the repo root:

```bash
quarto render visualization_scripts/visual_combined.qmd
```
Or we can just run the codes chunk by chunk.  

Rendered figures go to `visualization_results/`. The main results figure used in the paper is [plots/combined_three_panel_figure.pdf](../plots/combined_three_panel_figure.pdf).

R is required (`ggplot2`, `dplyr`, and packages listed in `scripts/visual_utilities.R`).
