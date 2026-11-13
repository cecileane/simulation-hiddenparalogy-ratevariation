
# For every replicate of one setting, 
# take the best find_graphs graph
# (lowest score) under h = 1 and write its two displayed trees: 
# the major tree (minor admixture edge removed) and the minor tree (major edge removed), 
# with their RF distance to the true species tree and the clades
# that the true tree does not have. 

# findgraphs_postprocess.jl only records whether one of the two displayed trees 
# is the true tree, so this script is supplied as additional testing. 

# Input: <output_dir>/<paramname>/rep*/findgraph/rep*_admix1_unique_graphs.rds
#   and rep*_admix1_summary_table.txt, from findgraphs_1rep.R

# Output: <output_dir>/<paramname>/findgraph-bestgraph-<paramname>.csv,
#   one row per replicate

# Usage: Rscript scripts/findgraphs_bestgraph_trees.R --tree plant \
#   --dup_rate 0.0 --loss_rate 0.0 --ratevar L \
# --n_inds 1 --SF 1.0


library(optparse)
library(igraph)
library(ape)

opts <- list(
  make_option("--tree", type = "character", default = "reptile",
              help = "reptile or plant"),
  make_option("--dup_rate", type = "character", default = "0.0"),
  make_option("--loss_rate", type = "character", default = "0.0"),
  make_option("--ratevar", type = "character", default = "N"),
  make_option("--n_inds", type = "integer", default = 1L),
  make_option("--SF", type = "character", default = "1.0"),
  make_option("--gene_len", type = "integer", default = 1000L),
  make_option("--n_reps", type = "integer", default = 100L),
  make_option("--output_dir", type = "character", default = NULL,
              help = "folder holding the setting folders (default: output, or
                output_plant for the plant tree)"))

opt <- parse_args(OptionParser(option_list = opts))

paramname <- sprintf("DUP%s-LOS%s-RV%s-N_ind%d-SF%s-genelen%d", opt$dup_rate,
                     opt$loss_rate, opt$ratevar, opt$n_inds, opt$SF,
                     opt$gene_len)
output_dir <- if (!is.null(opt$output_dir)) opt$output_dir else
  if (opt$tree == "reptile") "output" else paste0("output_", opt$tree)
setting_dir <- file.path(output_dir, paramname)

stopifnot(dir.exists(setting_dir))

true_newick <- if (startsWith(opt$tree, "plant"))
  "(A,(B,((C,D),(E,(H,(F,G))))));" 
  else "(A,((((B,C),(D,E)),F),(G,H)));"
true_tree <- read.tree(text = true_newick)

# clades of a rooted tree as sorted taxon strings, root clade excluded
clades <- function(tree) {
  tips <- tree$tip.label
  s <- sapply(prop.part(tree), function(p) paste(sort(tips[p]), collapse = ","))
  setdiff(s, paste(sort(tips), collapse = ","))
}
true_clades <- clades(true_tree)

# newick of a directed tree (igraph) from its root, leaves named by taxon
to_newick <- function(g) {
  rec <- function(v) {
    kids <- igraph::neighbors(g, v, mode = "out")$name
    if (length(kids) == 0) return(v)
    paste0("(", paste(sapply(kids, rec), collapse = ","), ")")
  }
  root <- names(which(igraph::degree(g, mode = "in") == 0))
  paste0(rec(root), ";")
}

# drop one admixture edge, then remove the dead ends it leaves 
# (internal nodes without children) and suppress the nodes left with a single child
# (vertices are handled by name; igraph:: 
# everywhere because ape masks degree)
displayed_tree <- function(edges, drop, taxa) {
  keep <- edges[-drop, c("from", "to")]
  g <- igraph::graph_from_data_frame(keep, directed = TRUE)
  repeat {
    deg <- igraph::degree(g, mode = "out")
    dead <- setdiff(names(deg)[deg == 0], taxa)
    single <- names(deg)[deg == 1]
    if (length(dead) > 0) {
      g <- igraph::delete_vertices(g, dead[1])
      next
    }
    if (length(single) == 0) break
    v <- single[1]
    kid <- igraph::neighbors(g, v, mode = "out")$name
    par <- igraph::neighbors(g, v, mode = "in")$name
    if (length(par) == 1) g <- igraph::add_edges(g, c(par, kid))
    g <- igraph::delete_vertices(g, v)
  }
  to_newick(g)
}

rows <- list()
for (r in seq_len(opt$n_reps)) {
  rep_id <- sprintf(paste0("%0", nchar(opt$n_reps), "d"), r)
  fg_dir <- file.path(setting_dir, paste0("rep", rep_id), "findgraph")
  rds <- file.path(fg_dir, paste0("rep", rep_id, "_admix1_unique_graphs.rds"))
  tab <- file.path(fg_dir, paste0("rep", rep_id, "_admix1_summary_table.txt"))
  if (!file.exists(rds) || !file.exists(tab)) next
  graphs <- readRDS(rds)
  best <- graphs[[which.min(sapply(graphs, function(g) g$score))]]
  summ <- read.table(tab, header = TRUE, stringsAsFactors = FALSE)
  brow <- summ[which.min(summ$score), ]
  edges <- as.data.frame(best$edges)
  adm <- which(edges$type == "admix")
  if (length(adm) != 2) next
  major_edge <- adm[which.max(edges$weight[adm])]
  minor_edge <- adm[which.min(edges$weight[adm])]
  major <- displayed_tree(edges, minor_edge, true_tree$tip.label)
  minor <- displayed_tree(edges, major_edge, true_tree$tip.label)
  tmaj <- read.tree(text = major)
  tmin <- read.tree(text = minor)
  rf <- function(t) dist.topo(unroot(t), unroot(true_tree))
  rows[[length(rows) + 1]] <- data.frame(
    rep_id = r, score = best$score, worst_residual = best$wr,
    gamma_major = edges$weight[major_edge],
    gamma_minor = edges$weight[minor_edge],
    hybrid_taxon = brow$hybrid_taxon, major_donor = brow$major_donor,
    minor_donor = brow$minor_donor,
    major_tree = major, minor_tree = minor,
    rf_major = rf(tmaj), rf_minor = rf(tmin),
    major_wrong_clades = paste(setdiff(clades(tmaj), true_clades),
                               collapse = ";"),
    minor_wrong_clades = paste(setdiff(clades(tmin), true_clades),
                               collapse = ";"),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
outfile <- file.path(setting_dir, paste0("findgraph-bestgraph-", paramname,
                                         ".csv"))
write.csv(out, outfile, row.names = FALSE)
cat(nrow(out), "replicates written to", outfile, "\n")
