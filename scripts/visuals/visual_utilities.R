# shared pieces of the two-tree figure scripts 
# (visual_main_figures.R, visual_appendix_figures.R, visual_tables.R): 
# colors and labels, dup/loss rate -> hidden paralogy level, 
# the loader of the per-replicate summaries
# (load_method), 
# the T* calibration and the y axis rule of the mirrored histograms

# source("scripts/visuals/visual_utilities.R") from the repo root
# load_method() reads results/summary_reptile/<method>/ and
# results/summary_plant/<method>/

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

results_root <- "visualization_results" # where the figures and tables go

# colors and labels for the substitution rate variation
rv_colors <- c("none" = "#898686", "gene" = "#039afe", "lineage" = "#fbc531")
rv_labels <- c("none" = "none", "gene" = "across genes",
  "lineage" = "across lineages")
tree_labels <- c("reptile" = "reptile", "plant" = "angiosperm")
tree_short <- c("reptile" = "reptile", "plant" = "plant")
# T* lines of PhyNEST: one color per tree, line type by ILS level
tstar_colors <- c("reptile" = "#7f93a8", "plant" = "aquamarine4")
tstar_types <- c("high" = "solid", "low" = "dashed")
hp_levels <- c("none", "medium", "high") # hidden paralogy levels

# the taxon dropped in the "noF" columns: F for reptiles, B for plants
pruned <- c("reptile" = "F", "plant" = "B")

# hidden paralogy level of a dup/loss rate
# the two trees use different rates for the same level: 0 / 3e-4 / 4e-4
# for the reptiles and 0 / 1e-4 / 1.4e-4 for the plants
rate_to_hp <- function(tree, rate) {
  r <- as.numeric(rate)
  case_when(r == 0 ~ "none",
      tree == "reptile" & abs(r - 3e-4) < 1e-9 ~ "medium",
      tree == "reptile" & abs(r - 4e-4) < 1e-9 ~ "high",
      tree == "plant"   & abs(r - 1e-4) < 1e-9 ~ "medium",
      tree == "plant"   & abs(r - 1.4e-4) < 1e-9 ~ "high",
      TRUE ~ NA_character_)
}

# the parameter setting from a summary file name
# (DUP*-LOS*-RV*-N_ind*-SF*-genelen*)
parse_setting <- function(file) {
  b <- basename(file)
  tibble(file = file,
    rate = str_match(b, "DUP([0-9.eE+-]+)-LOS")[, 2],
    rv = c(N = "none", G = "gene",
      L = "lineage")[str_match(b, "-RV([NGL])-")[, 2]],
    n = as.integer(str_match(b, "N_ind([12])")[, 2]),
    sf = as.numeric(str_match(b, "SF([0-9.]+)")[, 2]))
}

# read all the summary files of one folder: one row per setting,
# the per-replicate table in the list column "data"
read_rep_files <- function(dir, pattern, tree) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  parse_setting(files) %>%
    mutate(tree = tree, hp = rate_to_hp(tree, rate),
           ils = ifelse(sf == 1, "high", "low")) %>%
    filter(!is.na(hp)) %>%
    mutate(data = lapply(file, function(f) suppressMessages(
      read_csv(f, show_col_types = FALSE, guess_max = 200))))
}

# file name pattern and folders per method
method_spec <- list(
  findgraphs = list(pat = "^findgraph-DUP.*\\.csv$", sub = "findgraphs"),
  snaq = list(pat = "^SNaQ-DUP.*-summary\\.csv$", sub = "snaq"),
  phynest = list(pat = "^PhyNEST-DUP.*-summary\\.csv$", sub = "phynest"),
  simulation = list(pat = "^summary_DUP.*\\.csv$", sub = "simulation"),
  # displayed trees of the best h = 1 graph (findgraphs_bestgraph_trees.R)
  findgraphs_bestgraph = list(pat = "^findgraph-bestgraph-DUP.*\\.csv$",
                              sub = "findgraphs_bestgraph"),
  # per setting and locus category (genetree_tiers.jl), not per replicate
  genetrees = list(pat = "^genetree_tiers_DUP.*\\.csv$", sub = "genetrees"))

# load the per-replicate summaries of one method, for both trees
# method: "findgraphs", "snaq", "phynest", "simulation", 
# ... (see method_spec)
# returns one row per setting: tree, hp, rv, n, ils, file, data 
# (list of per-replicate tibbles) 
load_method <- function(method) {
  spec <- method_spec[[method]]
  out <- list()
  for (tree in c("reptile", "plant")) {
    tbl <- read_rep_files(file.path("results", paste0("summary_", tree),
                                    spec$sub),
                          spec$pat, tree)
    if (is.null(tbl) || !nrow(tbl)) next
    out[[tree]] <- tbl
  }
  bind_rows(out) %>% mutate(method = method) %>%
    mutate(hp = factor(hp, hp_levels), rv = factor(rv, names(rv_colors)),
           tree = factor(tree, names(tree_labels)),
           ils = factor(ils, c("high", "low")), n = factor(n, c(1, 2)))
}

# T* of PhyNEST: 95th percentile of the per-site T, pooled over the
# settings without lineage rate variation, one T* per tree and ILS level
t_star_from <- function(ph) ph %>% filter(rv != "lineage", n == 1) %>%
  select(tree, ils, data) %>% unnest(data) %>%
  group_by(tree, ils) %>%
  summarise(t_star = quantile(T_per_site, 0.95), n_pooled = n(),
            .groups = "drop")

# y axis of the mirrored histograms (Fig. 2 and appendix G, H)
# the limit of a row is the tallest bar rounded up to a round number (next 5
# for counts up to 100, next 10 up to 300, next 50 up to 1000 ...),
# plus one more unit when the bar would touch the top.

# the step is half the limit when that is a round number (30 gives 15,
# 50 gives 25), otherwise a third (45 gives 15, 60 gives 20),
# so there are at most three labels on each half of the axis
NICE_STEPS <- c(5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 150, 200, 250,
                300, 400, 500, 750, 1000)

nice_axis <- function(top) {
  unit <- if (top <= 100) 5 else if (top <= 300) 10 else
    if (top <= 1000) 50 else 100
  lim <- ceiling(top / unit) * unit
  if (lim - top < unit / 5) lim <- lim + unit
  step <- if ((lim / 2) %in% NICE_STEPS) lim / 2 else
    if ((lim / 3) %in% NICE_STEPS) lim / 3 else
      max(NICE_STEPS[NICE_STEPS <= lim / 2])
  list(limit = lim, step = step)
}

# limit and step for every level of the grouping column, from the signed
# bar heights in column h
nice_axes <- function(df, group, h) {
  m <- df %>% group_by(.data[[group]]) %>%
    summarise(top = max(abs(.data[[h]]), na.rm = TRUE), .groups = "drop")
  ax <- lapply(m$top, nice_axis)
  m$limit <- sapply(ax, function(a) a$limit)
  m$step <- sapply(ax, function(a) a$step)
  m
}

# breaks for scale_y_continuous(): the panel limits tell which row we are
# in, this needs expand = c(0, 0) so that the panel limit is the row limit
mirror_breaks <- function(ax) function(l) {
  top <- max(abs(l))
  i <- which.min(abs(ax$limit - top))
  b <- seq(0, ax$limit[i], ax$step[i])
  c(-rev(b[-1]), b)
}
