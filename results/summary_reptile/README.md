Per-replicate summaries of the reptile simulations, one csv per setting and
method, copied here from the run folders on the servers. These are the
files the figure and table scripts in scripts/visuals/ read.

- simulation/, snaq/, findgraphs/, phynest/: the summary written by
  scripts/<method>_postprocess.jl in each setting folder, one row per
  replicate (100 rows)
- findgraphs_bestgraph/: the two trees displayed by the best h = 1
  find_graphs graph of every replicate, with their RF distance to the true
  tree (scripts/findgraphs_bestgraph_trees.R)
- genetrees/: the branch counts of scripts/genetree_tiers.jl, one row per
  locus category, for the gene tree discordance and estimation error tables

Column dictionaries: ../readme-dryad.md
