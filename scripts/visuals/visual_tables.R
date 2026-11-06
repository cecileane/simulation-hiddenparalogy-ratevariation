#!/usr/bin/env Rscript
# ============================================================================
# scripts/visuals/visual_tables.R
#
# Supplementary tables of the two-tree comparison, each with a reptile and a
# plant sub-column under every statistic:
#   wr_summary       worst residual (WR) of the true tree, find_graphs,
#                    pooled over the replicates of each factor level
#   T_summary        per-site score difference T of PhyNEST (one individual
#                    per taxon), same layout
#   model_choice     % of replicates selecting h = 1 and h > 1 (SNaQ,
#                    find_graphs at WR <= 3.0 and 3.7) or h >= 1 (PhyNEST,
#                    T > T*), per factor level; one CSV, two tabulars
#                    (_findgraphs, _snaq_phynest)
#   genetree_discordance      normalized RF between the true gene trees and
#                    the species tree and % of gene trees by number of
#                    discordant branches, per locus category and ILS level
#   genetree_estimation_error same between the true and the estimated gene
#                    trees, per rate variation, one table per number of
#                    individuals (_n1, _n2)
#   gamma_<method>   minor gamma of the h = 1 network (find_graphs: best
#                    graph; SNaQ; PhyNEST, one individual per taxon), per
#                    factor level: mean, median, SD, % below 0.05/0.10/0.25
#                    These two give the mean over the settings sharing the
#                    row's factor levels and the range across them.
# Inputs  : per-replicate summaries of the kappa-corrected rerun, read with
#           load_method() from visual_utilities.R; per-setting branch counts
#           of genetree_tiers.jl for the last two tables
# Outputs : <results root>/tables/<name>.csv  (long format, one row per
#           tree x factor level x statistic) and <name>.tex (the tabular
#           only, to \input in the supplement)
# Usage   : Rscript scripts/visuals/visual_tables.R   (repo root)
# ============================================================================

source("scripts/visuals/visual_utilities.R")

out_dir <- file.path(results_root, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fg <- load_method("findgraphs")
sq <- load_method("snaq")
ph <- load_method("phynest")
t_star <- t_star_from(ph)

# ---------------------------------------------------------------------------
# Factor levels: every replicate of a setting carries the setting's factors,
# so pooling within a level is a filter on the long table
# ---------------------------------------------------------------------------
rv_level <- c(none = "none", gene = "across genes", lineage = "across lineages")
factor_levels <- list(
  list(factor = "overall", level = "all",
       keep = function(d) rep(TRUE, nrow(d))),
  list(factor = "rate variation", level = "none",
       keep = function(d) d$rv == "none"),
  list(factor = "rate variation", level = "across genes",
       keep = function(d) d$rv == "gene"),
  list(factor = "rate variation", level = "across lineages",
       keep = function(d) d$rv == "lineage"),
  list(factor = "duplication/loss", level = "none",
       keep = function(d) d$hp == "none"),
  list(factor = "duplication/loss", level = "medium",
       keep = function(d) d$hp == "medium"),
  list(factor = "duplication/loss", level = "high",
       keep = function(d) d$hp == "high"),
  list(factor = "ILS", level = "high", keep = function(d) d$ils == "high"),
  list(factor = "ILS", level = "low", keep = function(d) d$ils == "low"),
  list(factor = "individuals/taxon", level = "1", keep = function(d) d$n == 1),
  list(factor = "individuals/taxon", level = "2", keep = function(d) d$n == 2))

# one row per replicate with the setting's factors and the chosen columns
replicates <- function(tbl, f) tbl %>%
  select(tree, hp, rv, n, ils, data) %>%
  mutate(data = lapply(data, f)) %>% unnest(data)

# apply the statistics to every tree x factor level; `stats` is a named list
# of functions of a numeric vector
pooled_table <- function(rep, value, stats, levels = factor_levels) {
  rows <- list()
  for (tree in c("reptile", "plant")) for (fl in levels) {
    x <- rep[[value]][rep$tree == tree & fl$keep(rep)]
    x <- x[!is.na(x)]
    if (!length(x)) next
    for (s in names(stats))
      rows[[length(rows) + 1]] <- tibble(
        tree = tree, factor = fl$factor, level = fl$level, n = length(x),
        statistic = s, value = stats[[s]](x))
  }
  bind_rows(rows)
}

# ---------------------------------------------------------------------------
# LaTeX writer: factor / level rows, a reptile and a plant column under
# every statistic (booktabs, multirow for the factor names). Rows keep the
# order of first appearance in `long`. With columns lo and hi in `long`,
# every cell is the value over its range, stacked with makecell.
# ---------------------------------------------------------------------------
tex_table <- function(long, stat_names, digits, file,
                      stat_labels = stat_names,
                      head = c("factor", "level")) {
  cols <- as.vector(t(outer(stat_names, c("reptile", "plant"), paste,
                            sep = "__")))
  ord <- unique(paste(long$factor, long$level))
  wide_of <- function(what) {
    w <- long %>% filter(statistic %in% stat_names) %>%
      mutate(col = paste(statistic, tree, sep = "__")) %>%
      select(factor, level, col, value = all_of(what)) %>%
      pivot_wider(names_from = col, values_from = value)
    for (cc in setdiff(cols, names(w))) w[[cc]] <- NA_real_
    w[order(match(paste(w$factor, w$level), ord)), ]
  }
  wide <- wide_of("value")
  has_range <- all(c("lo", "hi") %in% names(long))
  if (has_range) { wlo <- wide_of("lo"); whi <- wide_of("hi") }
  fmt <- function(v, d) formatC(v, format = "f", digits = d)
  cell <- function(k, i) {
    v <- wide[[cols[k]]][i]
    d <- digits[ceiling(k / 2)]
    if (is.na(v)) return("--")
    if (!has_range) return(fmt(v, d))
    sprintf("\\makecell{%s\\\\(%s--%s)}", fmt(v, d),
            fmt(wlo[[cols[k]]][i], d), fmt(whi[[cols[k]]][i], d))
  }
  lines <- c(
    sprintf("\\begin{tabular}{ll%s}", paste(rep("cc", length(stat_names)),
                                             collapse = "")),
    "\\toprule",
    paste0("& & ", paste(sprintf("\\multicolumn{2}{c}{%s}", stat_labels),
                         collapse = " & "), " \\\\"),
    paste(sprintf("\\cmidrule(l{3pt}r{3pt}){%d-%d}",
                  2 * seq_along(stat_names) + 1,
                  2 * seq_along(stat_names) + 2), collapse = " "),
    paste0(head[1], " & ", head[2], " & ",
           paste(rep("reptile & plant", length(stat_names)), collapse = " & "),
           " \\\\"),
    "\\midrule")
  prev <- ""
  for (i in seq_len(nrow(wide))) {
    f <- wide$factor[i]
    if (f != prev && i > 1) lines <- c(lines, "\\midrule")
    nlev <- sum(wide$factor == f)
    fcell <- if (f != prev) {
      if (nlev > 1) sprintf("\\multirow{%d}{*}{%s}", nlev, f) else f
    } else ""
    vals <- sapply(seq_along(cols), function(k) cell(k, i))
    lines <- c(lines, paste0(fcell, " & ", wide$level[i], " & ",
                             paste(vals, collapse = " & "), " \\\\"))
    prev <- f
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file)
}

# ---------------------------------------------------------------------------
# 1. worst residual of the true tree (find_graphs)
# ---------------------------------------------------------------------------
stats5 <- list(mean = mean, median = median, sd = sd,
               q95 = function(x) quantile(x, 0.95),
               q99 = function(x) quantile(x, 0.99))
stat5_labels <- c("mean", "median", "SD", "95th percentile",
                  "99th percentile")
wr_rep <- replicates(fg, function(d) tibble(wr = d$true_tree_wr))
wr_tab <- pooled_table(wr_rep, "wr", stats5)
write_csv(wr_tab, file.path(out_dir, "wr_summary.csv"))
tex_table(wr_tab, names(stats5), rep(2, 5),
          file.path(out_dir, "wr_summary.tex"), stat5_labels)

# ---------------------------------------------------------------------------
# 2. per-site score difference T (PhyNEST, one individual per taxon)
# ---------------------------------------------------------------------------
T_rep <- replicates(ph %>% filter(n == 1), function(d) tibble(T = d$T_per_site))
T_levels <- factor_levels[sapply(factor_levels,
                                 function(fl) fl$factor != "individuals/taxon")]
T_tab <- pooled_table(T_rep, "T", stats5, T_levels)
write_csv(T_tab, file.path(out_dir, "T_summary.csv"))
tex_table(T_tab, names(stats5), rep(4, 5),
          file.path(out_dir, "T_summary.tex"), stat5_labels)
write_csv(t_star, file.path(out_dir, "T_star.csv"))

# ---------------------------------------------------------------------------
# 3. model choice: % of replicates selecting h = 1 and h > 1
# ---------------------------------------------------------------------------
# selected h per replicate: 0, 1 or 2 (for h > 1)
k_of <- function(x) {
  x <- as.character(x); ifelse(x == "0", 0L, ifelse(x == "1", 1L, 2L))
}
mc_rep <- bind_rows(
  replicates(sq, function(d) tibble(
    h = ifelse(d$p_H0 > 0.05, 0L, ifelse(d$p_H1 > 0.05, 1L, 2L)))) %>%
    mutate(method = "SNaQ"),
  replicates(fg, function(d) tibble(h = k_of(d$best_k))) %>%
    mutate(method = "fg30"),
  replicates(fg, function(d) tibble(h = k_of(d$best_k_new_WR_3.7))) %>%
    mutate(method = "fg37"),
  replicates(ph %>% left_join(t_star, by = c("tree", "ils")) %>% rowwise() %>%
               mutate(data = list(mutate(data, ts = t_star))) %>% ungroup(),
             function(d) tibble(h = as.integer(d$T_per_site > d$ts))) %>%
    mutate(method = "phynest"))
# the two percentages of every method become statistics of one table
pct_h <- function(k) function(h) 100 * mean(h == k)
mc_tab <- bind_rows(lapply(c("SNaQ", "fg30", "fg37", "phynest"), function(m) {
  st <- if (m == "phynest") list(h1 = pct_h(1)) else
    list(h1 = pct_h(1), hgt1 = pct_h(2))
  pooled_table(mc_rep %>% filter(method == m), "h", st) %>%
    mutate(statistic = paste(m, statistic, sep = "_"))
}))
write_csv(mc_tab, file.path(out_dir, "model_choice.csv"))
# two tabulars: find_graphs at its two thresholds, then SNaQ and PhyNEST
tex_table(mc_tab, c("fg30_h1", "fg30_hgt1", "fg37_h1", "fg37_hgt1"),
          rep(1, 4), file.path(out_dir, "model_choice_findgraphs.tex"),
          c("$\\text{WR} \\leq 3.0$, $h=1$", "$\\text{WR} \\leq 3.0$, $h>1$",
            "$\\text{WR} \\leq 3.7$, $h=1$", "$\\text{WR} \\leq 3.7$, $h>1$"))
tex_table(mc_tab, c("SNaQ_h1", "SNaQ_hgt1", "phynest_h1"), rep(1, 3),
          file.path(out_dir, "model_choice_snaq_phynest.tex"),
          c("SNaQ $h=1$", "SNaQ $h>1$", "\\ph{} $h \\geq 1$"))

# ---------------------------------------------------------------------------
# 4. gene tree discordance and gene tree estimation error, in branches:
#    one value per setting (genetree_tiers.jl), then the mean over the
#    settings sharing the row's factor levels and the range across them
# ---------------------------------------------------------------------------
gt <- load_method("genetrees") %>% select(tree, hp, rv, n, ils, data) %>%
  unnest(data)
# `stats` are functions of the per-setting table giving one value per row
range_table <- function(tbl, rows, stats) {
  out <- list()
  for (tree in c("reptile", "plant")) for (r in rows) {
    d <- tbl[tbl$tree == tree & r$keep(tbl), ]
    if (!nrow(d)) next
    for (s in names(stats)) {
      x <- stats[[s]](d)
      out[[length(out) + 1]] <- tibble(
        tree = tree, factor = r$factor, level = r$level, n = nrow(d),
        statistic = s, value = mean(x), lo = min(x), hi = max(x))
    }
  }
  bind_rows(out)
}
pct <- function(d, cols) 100 * rowSums(d[, cols, drop = FALSE])

# discordance: locus category x ILS level, all settings with the category
categories <- c(ortholog = "orthologs", falseHP = "false hidden paralogs",
                weakHP = "weak hidden paralogs",
                strongHP = "strong hidden paralogs")
disc_rows <- list()
for (cat in names(categories)) for (ils in c("high", "low"))
  disc_rows[[length(disc_rows) + 1]] <- list(
    factor = categories[[cat]], level = ils,
    keep = local({ cat <- cat; ils <- ils
      function(d) d$category == cat & d$ils == ils }))
disc_stats <- list(
  rf = function(d) d$mean_RF_species,
  d0 = function(d) pct(d, "d0"), d1 = function(d) pct(d, "d1"),
  d2 = function(d) pct(d, "d2"), d3 = function(d) pct(d, c("d3", "d4", "d5")))
disc_tab <- range_table(gt, disc_rows, disc_stats)
write_csv(disc_tab, file.path(out_dir, "genetree_discordance.csv"))
tex_table(disc_tab, names(disc_stats), c(3, 1, 1, 1, 1),
          file.path(out_dir, "genetree_discordance.tex"),
          c("mean normalized RF", "0 (identical)", "1", "2", "$\\geq 3$"),
          c("locus category", "ILS"))

# estimation error: rate variation, all genes, one table per number of
# individuals because the tiers differ (5 and 13 internal branches)
err_rows <- lapply(names(rv_level), function(v) list(
  factor = "rate variation", level = rv_level[[v]],
  keep = function(d) d$category == "all" & d$rv == v))
err_stats <- list(
  "1" = list(rf = function(d) d$mean_RF_estimated,
             e0 = function(d) pct(d, "e0"), e1 = function(d) pct(d, "e1"),
             e2 = function(d) pct(d, "e2"),
             e3 = function(d) pct(d, c("e3", "e4", "e5")),
             all = function(d) pct(d, "e5")),
  "2" = list(rf = function(d) d$mean_RF_estimated,
             e0 = function(d) pct(d, "e0"),
             e12 = function(d) pct(d, c("e1", "e2")),
             e35 = function(d) pct(d, c("e3", "e4", "e5")),
             e6 = function(d) pct(d, paste0("e", 6:13)),
             all = function(d) pct(d, "e13")))
err_labels <- list(
  "1" = c("mean normalized RF", "0 (exact)", "1", "2", "$\\geq 3$", "all 5"),
  "2" = c("mean normalized RF", "0 (exact)", "1--2", "3--5", "$\\geq 6$",
          "all 13"))
err_tab <- bind_rows(lapply(c("1", "2"), function(k)
  range_table(gt %>% filter(n == k), err_rows, err_stats[[k]]) %>%
    mutate(n_inds = as.integer(k))))
write_csv(err_tab, file.path(out_dir, "genetree_estimation_error.csv"))
for (k in c("1", "2"))
  tex_table(err_tab %>% filter(n_inds == as.integer(k)),
            names(err_stats[[k]]), c(3, 1, 1, 1, 1, 1),
            file.path(out_dir, paste0("genetree_estimation_error_n", k,
                                      ".tex")), err_labels[[k]], c("", ""))

# ---------------------------------------------------------------------------
# 5. minor gamma of the h = 1 network: mean, median, SD and the % of
#    replicates below 0.05, 0.10 and 0.25, one table per method
# ---------------------------------------------------------------------------
below <- function(c) function(x) 100 * mean(x < c)
gamma_stats <- list(mean = mean, median = median, sd = sd,
                    lt05 = below(0.05), lt10 = below(0.10),
                    lt25 = below(0.25))
gamma_labels <- c("mean", "median", "SD", "$\\gamma < 0.05$ (\\%)",
                  "$\\gamma < 0.10$ (\\%)", "$\\gamma < 0.25$ (\\%)")
gamma_rep <- list(
  findgraphs = replicates(fg, function(d) tibble(
    gamma = as.numeric(d$best_graph_gamma2))),
  snaq = replicates(sq, function(d) tibble(gamma = as.numeric(d$gamma_2))),
  phynest = replicates(ph, function(d) tibble(gamma = as.numeric(d$gamma_2))))
for (m in names(gamma_rep)) {
  tab <- pooled_table(gamma_rep[[m]], "gamma", gamma_stats,
                      if (m == "phynest") T_levels else factor_levels)
  write_csv(tab, file.path(out_dir, paste0("gamma_", m, ".csv")))
  tex_table(tab, names(gamma_stats), c(3, 3, 3, 1, 1, 1),
            file.path(out_dir, paste0("gamma_", m, ".tex")), gamma_labels)
}

cat("tables written to", out_dir, "\n")
print(as.data.frame(t_star))
