
# appendix figures of the two-tree comparison
# input: the per-replicate summaries, with load_method()

# output, in visualization_results/appendix_two_trees/:
# A_genetree_rf: gene tree discordance and estimation error
# B_wr_T_distributions: A worst residual of the true tree for every
#   replicate (find_graphs), B per-site T for every replicate with T* (PhyNEST)
# D_phynest_T_placement: T of the settings without lineage rates,
#   shape = where the h = 1 reticulation was placed
# E_model_choice: stacked bars of the selected h
# F_true_tree_recovery: % of replicates displaying the true tree under
#   h = 0, h = 1, and h = 1 with the rogue taxon removed
# G_minor_gamma: mirrored histograms of the minor gamma under h = 1
# H_bestgraph: minor gamma of the best h = 1 graph by whether its major tree
#   is the true tree, and its reticulation

# Rscript scripts/visuals/visual_appendix_figures.R, from the repo root
# MODEL_CHOICE_PLANT=solid Rscript ... picks the color of the plant bars in
# the model choice figures: "solid" (same colors as the reptile bars),
# "light" (a lighter tint) or "both" (default: both versions, the light
# one with the suffix _lightplant)

# same encoding as Fig. 1: x = hidden paralogy level (none / medium / high)
# split into a reptile and a plant column, the plant column on a very light
# grey band, solid lines between the levels. color = substitution rate
# variation, shape = ILS (circle high, triangle low), fill = individuals per
# taxon (filled n = 1, hollow n = 2)

library(ggplot2)
library(cowplot)
library(scales)
source("scripts/visuals/visual_utilities.R")

MODEL_CHOICE_PLANT <- match.arg(Sys.getenv("MODEL_CHOICE_PLANT", "both"),
                                c("both", "solid", "light"))
out_root <- file.path(results_root, "appendix_two_trees")

odir <- function(name) {
  d <- file.path(out_root, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

save_both <- function(p, dir, stem, w, h) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 130, bg = "white")
}

hp_dl <- c(none = "no duplication & loss", medium = "medium duplication & loss",
           high = "high duplication & loss")
ils_levels <- c("ILS = high", "ILS = low")

fg <- load_method("findgraphs")
sq <- load_method("snaq")
ph <- load_method("phynest")
sim <- load_method("simulation")
t_star <- t_star_from(ph)
cat("settings: findgraphs", nrow(fg), " snaq", nrow(sq), " phynest", nrow(ph),
    " simulation", nrow(sim), "\n")
print(as.data.frame(t_star))


#------- shared pieces -------#
BASE <- 17 # base font size of the panels
PT <- 3.4 # point size
shape_ils <- c("high" = 21, "low" = 24) # circle, triangle
# x offsets: tree within a level, rate variation within a tree
tree_off <- c("reptile" = -0.2, "plant" = 0.2)
rv_off <- c("none" = -0.085, "gene" = 0, "lineage" = 0.085)
tree_x <- as.vector(outer(c(-0.2, 0.2), 1:3, "+")) # the six sub-columns

# x position of a setting (as in Fig. 1), a small jitter spreads the
# replicates in the per-replicate figures
x_nested <- function(hp, tree, rv, n = NULL, jitter = 0) {
  as.numeric(hp) + tree_off[as.character(tree)] + rv_off[as.character(rv)] +
    (if (is.null(n)) 0 else ifelse(n == 1, -0.035, 0.035)) +
    (if (jitter > 0) runif(length(hp), -jitter, jitter) else 0)
}

# the panel spans exactly the three levels, so the edges of the blocks fall
# on the panel border and on the separating lines
scale_x_nested <- function(top = TRUE, fixed = TRUE) scale_x_continuous(
  breaks = tree_x, labels = rep(tree_short, 3),
  limits = if (fixed) c(0.5, 3.5) else NULL,
  expand = expansion(add = 0), minor_breaks = NULL,
  sec.axis = if (top) dup_axis(breaks = 1:3, labels = hp_levels, name = NULL)
             else waiver())

# background: each level is the block [x - 0.5, x + 0.5], its plant half
# [x, x + 0.5] on a very light grey band, solid lines at x +- 0.5
# vlines can have a facet column, so that a panel only gets its own levels
block_guides <- function(ylo = -Inf, vlines = tibble(x = 1:3)) list(
  geom_rect(data = vlines %>% mutate(xmin = x, xmax = x + 0.5),
            aes(xmin = xmin, xmax = xmax, ymin = ylo, ymax = Inf),
            inherit.aes = FALSE, fill = "grey94", alpha = 0.5),
  geom_vline(data = vlines %>% mutate(x = x + 0.5), aes(xintercept = x),
             color = "grey55", linewidth = 0.7),
  geom_vline(data = vlines %>% mutate(x = x - 0.5), aes(xintercept = x),
             color = "grey55", linewidth = 0.7))
theme_tt <- function(base = BASE) theme_bw(base) + theme(
  # strip on top, the level labels between the strip and the panel
  strip.placement = "outside",
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "grey78", linewidth = 0.55),
  strip.background = element_rect(fill = "grey90"),
  strip.text = element_text(face = "bold", size = base - 1),
  axis.text = element_text(color = "black", size = base - 2),
  axis.title = element_text(size = base),
  axis.text.x.bottom = element_text(size = base - 2, angle = 45, hjust = 1,
                                    vjust = 1),
  axis.text.x.top = element_text(size = base - 1, face = "bold"),
  axis.ticks.x.top = element_blank(),
  legend.position = "bottom", legend.box = "horizontal",
  legend.text = element_text(size = base - 3),
  legend.title = element_text(size = base - 3, face = "bold"),
  legend.key.size = unit(1.1, "lines"), legend.spacing.x = unit(0.4, "lines"),
  legend.background = element_rect(color = "black"),
  plot.margin = margin(4, 6, 2, 4))

# points for per-setting values: color = rate variation, shape = ILS, fill = n
pts_setting <- function(df, y, size = PT, fixed_x = TRUE) list(
  geom_point(data = df %>% mutate(
               x = x_nested(hp, tree, rv, n),
               fillc = ifelse(n == 1, rv_colors[as.character(rv)], "white")),
             aes(x = x, y = .data[[y]], shape = ils, color = rv, fill = fillc),
             size = size, stroke = 1.1),
  scale_fill_identity(),
  scale_shape_manual(values = shape_ils,
                     labels = c(high = "high ILS", low = "low ILS"),
                     name = "ILS",
                     guide = guide_legend(order = 2,
                                          override.aes = list(fill = "black"))),
  scale_color_manual(values = rv_colors, labels = rv_labels,
                     name = "rate variation",
                     guide = guide_legend(order = 1, override.aes = list(
                       shape = 21, fill = rv_colors))),
  scale_x_nested(fixed = fixed_x), labs(x = NULL))

# the fill of the points is an identity scale (no legend), so the figures
# that need the n legend get one made from a dummy data set
legend_rv_ils_n <- function(with_n = TRUE) {
  d <- expand_grid(rv = factor(names(rv_colors), names(rv_colors)),
                   ils = factor(c("high", "low"), c("high", "low")),
                   n = factor(c("n = 1", "n = 2")))
  p <- ggplot(d, aes(1, 1)) +
    geom_point(aes(color = rv), shape = 21, size = 4.5, stroke = 1.1,
               fill = NA) +
    geom_point(aes(shape = ils), size = 4.5, stroke = 1.1, color = "black",
               fill = "black") +
    scale_color_manual(values = rv_colors, labels = rv_labels,
                       name = "rate variation",
                       guide = guide_legend(order = 1, override.aes = list(
                         shape = 21, fill = rv_colors, color = rv_colors))) +
    scale_shape_manual(values = shape_ils,
                       labels = c(high = "high ILS", low = "low ILS"),
                       name = "ILS", guide = guide_legend(order = 2))
  if (with_n) p <- p +
    geom_point(aes(fill = n), shape = 21, size = 4.5, stroke = 1.1,
               color = "black") +
    scale_fill_manual(values = c("n = 1" = "black", "n = 2" = "white"),
                      name = "individuals / taxon",
                      guide = guide_legend(order = 3))
  cowplot::get_legend(p + theme_tt())
}

with_legend <- function(p, leg, rel = 0.08)
  plot_grid(p + theme(legend.position = "none"), leg, ncol = 1,
            rel_heights = c(1, rel))

#--- A. gene tree discordance and estimation error ---#
dA <- odir("A_genetree_rf")

sim_long <- sim %>% rowwise() %>% mutate(
  disc_nothing = data$mean_RF_genetree_vs_sptree_nothing,
  disc_false = data$mean_RF_genetree_vs_sptree_false_HP,
  disc_weak = data$mean_RF_genetree_vs_sptree_weak_HP,
  disc_strong = data$mean_RF_genetree_vs_sptree_strong_HP,
  disc_all = data$mean_RF_genetree_vs_sptree_all,
  err_nothing = data$mean_RF_true_vs_est_genetree_nothing,
  err_false = data$mean_RF_true_vs_est_genetree_false_HP,
  err_weak = data$mean_RF_true_vs_est_genetree_weak_HP,
  err_strong = data$mean_RF_true_vs_est_genetree_strong_HP,
  err_all = data$mean_RF_true_vs_est_genetree_all) %>%
  ungroup() %>% select(-data) %>%
  # without dup/loss every locus is in the "nothing" class: overall mean
  mutate(disc_nothing = ifelse(hp == "none", disc_all, disc_nothing),
         err_nothing = ifelse(hp == "none", err_all, err_nothing)) %>%
  pivot_longer(c(starts_with("disc_"), starts_with("err_")),
               names_to = c("what", "locus"), names_sep = "_") %>%
  filter(locus != "all", !is.na(value), is.finite(value)) %>%
  mutate(what = factor(what, c("disc", "err"),
                       c("A  gene tree discordance:\ntrue gene vs species tree",
                         paste0("B  gene tree estimation error:\n",
                                "true vs estimated gene tree"))),
         locus = factor(locus, c("nothing", "false", "weak", "strong"),
                        c("no duplication & loss", "false paralogy",
                          "weak hidden paralogy", "strong hidden paralogy")))

# the hidden-paralogy classes have no "none" level
vl_A <- sim_long %>% distinct(locus, hp) %>%
  transmute(locus, x = as.numeric(hp))
pA <- ggplot() + block_guides(vlines = vl_A) +
  pts_setting(sim_long, "value", fixed_x = FALSE) +
  facet_grid(what ~ locus, scales = "free", switch = "y") +
  labs(y = "mean normalized RF distance") + theme_tt() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90))
save_both(with_legend(pA, legend_rv_ils_n()), dA, "genetree_rf_two_trees",
          14, 9)
write_csv(sim_long, file.path(dA, "genetree_rf_per_setting.csv"))

#--- B. every replicate: WR of the true tree (find_graphs, panel A) ---#
dB <- odir("B_wr_T_distributions")
#-------- and per-site T with T* (PhyNEST, panel B) ------#
set.seed(1) # same jitter every time

wr_rep <- fg %>% select(tree, hp, rv, n, ils, data) %>%
  mutate(data = lapply(data, function(d) tibble(wr = d$true_tree_wr))) %>%
  unnest(data) %>%
  mutate(x = x_nested(hp, tree, rv, jitter = 0.03),
         n = factor(paste0("n = ", n)),
         ils = factor(paste0("ILS = ", ils), ils_levels))
pB <- ggplot() + block_guides(ylo = min(wr_rep$wr, na.rm = TRUE) * 0.9) +
  geom_point(data = wr_rep, aes(x, wr, color = rv), alpha = 0.55, size = 1.3,
             stroke = 0) +
  geom_hline(yintercept = 3.0, linetype = "dashed", color = "red",
             linewidth = 0.6) +
  geom_hline(yintercept = 3.7, linetype = "dashed", color = "green4",
             linewidth = 0.6) +
  scale_color_manual(values = rv_colors, labels = rv_labels,
                     name = "rate variation",
                     guide = guide_legend(override.aes = list(size = 4,
                                                              alpha = 1))) +
  scale_y_log10(breaks = c(1, 2, 3, 5, 10, 20, 30)) + scale_x_nested() +
  facet_grid(n ~ ils) +
  labs(x = NULL, y = "worst residual of the true tree topology (log scale)",
       title = "A  find_graphs") +
  theme_tt() + theme(plot.title = element_text(face = "bold", size = BASE + 1))

T_rep <- ph %>% filter(n == 1) %>%
  select(tree, hp, rv, ils, data) %>%
  mutate(data = lapply(data, function(d) tibble(
    T = d$T_per_site, hybrid_taxon = d$hybrid_taxon,
    major_donor = d$major_donor))) %>%
  unnest(data) %>%
  mutate(x = x_nested(hp, tree, rv, jitter = 0.03),
         ils_lab = factor(paste0("ILS = ", ils), ils_levels))
# one T* line per tree and ILS level, as in Fig. 1: color by tree, line type
# by ILS. The values are given in the caption
tstar_line <- t_star %>%
  mutate(ils_lab = factor(paste0("ILS = ", ils), ils_levels))
tstar_lines <- function() list(
  ggnewscale::new_scale_color(),
  geom_hline(data = tstar_line,
             aes(yintercept = t_star, color = tree, linetype = ils),
             linewidth = 0.7),
  scale_color_manual(values = tstar_colors, labels = tree_short,
                     name = "T* of"),
  scale_linetype_manual(values = tstar_types, guide = "none"))
y_T <- "per-site T (log scale)"
# T is negative in a few replicates, those are not shown on the log axis

pC <- ggplot() +
  block_guides(ylo = min(T_rep$T[T_rep$T > 0], na.rm = TRUE) * 0.9) +
  geom_point(data = T_rep, aes(x, T, color = rv), alpha = 0.55, size = 1.3,
             stroke = 0) +
  scale_color_manual(values = rv_colors, labels = rv_labels,
                     name = "rate variation",
                     guide = guide_legend(override.aes = list(size = 4,
                                                              alpha = 1))) +
  tstar_lines() +
  scale_y_log10(breaks = c(0.001, 0.003, 0.01, 0.03, 0.1),
                labels = label_number(accuracy = 0.001)) +
  scale_x_nested() + facet_wrap(~ils_lab) +
  labs(x = NULL, y = y_T, title = "B  PhyNEST") +
  theme_tt() + theme(plot.title = element_text(face = "bold", size = BASE + 1))
# both panels share the legend of B (it also has the T* line colors)
pBC <- plot_grid(pB + theme(legend.position = "none"), pC, ncol = 1,
                 rel_heights = c(2, 1.35), align = "v", axis = "lr")
save_both(pBC, dB, "wr_and_T_distributions_two_trees", 9.5, 13.5)
# each panel on its own, for the supplement
save_both(pB + labs(title = NULL), dB, "wr_distributions_two_trees", 9.5, 8)
save_both(pC + labs(title = NULL), dB, "T_distributions_two_trees", 9.5, 6)


#--- D. T without lineage rates, shape = where the reticulation goes ---#
dD <- odir("D_phynest_T_placement")

# "deep": the reticulation joins the deep clade and the rest of the ingroup
# (G,H for the reptiles, C,D for the plants)
deep_placement <- function(tree, hybrid, major) {
  ifelse(tree == "reptile",
         hybrid == "B,C,D,E,F" & grepl("^(G,H|G|H)$", major),
         hybrid %in% c("B,E,F,G,H", "E,F,G,H") & grepl("^(C,D|C|D)$", major))
}

T_pl <- T_rep %>% filter(rv != "lineage") %>%
  mutate(deep = deep_placement(tree, hybrid_taxon, major_donor),
         placement = factor(ifelse(deep, "deep", "elsewhere"),
                            c("deep", "elsewhere")))
pD <- ggplot() +
  block_guides(ylo = min(T_pl$T[T_pl$T > 0], na.rm = TRUE) * 0.9) +
  geom_point(data = T_pl, aes(x, T, color = rv, shape = placement),
             alpha = 0.55, size = 1.8, stroke = 0.6) +
  scale_shape_manual(values = c(deep = 17, elsewhere = 1),
                     name = "reticulation placed",
                     labels = c(deep = paste0(
                       "between the deep clade and the rest of the ingroup\n",
                       "(reptile: (G,H); angiosperm: (C,D))"),
                       elsewhere = "elsewhere")) +
  scale_color_manual(values = rv_colors[c("none", "gene")],
                     labels = rv_labels, name = "rate variation",
                     guide = guide_legend(override.aes = list(
                       size = 4, alpha = 1, shape = 16))) +
  tstar_lines() +
  scale_y_log10(breaks = c(0.001, 0.003, 0.01, 0.03, 0.1),
                labels = label_number(accuracy = 0.001)) +
  scale_x_nested() + facet_wrap(~ils_lab) +
  labs(x = NULL, y = y_T) + theme_tt() +
  theme(legend.box = "vertical")
save_both(pD, dD, "phynest_T_placement_two_trees", 12, 9.5)
write_csv(T_pl %>% count(tree, ils, rv, hp, placement),
          file.path(dD, "placement_counts.csv"))


#--- E. model choice ---#
dE <- odir("E_model_choice")
model_cols <- c("h = 0 (tree)" = "#2b6cb0", "h = 1 (network)" = "#f6ad55",
                "h > 1 (network)" = "#c53030")

# % of replicates per selected h, for one method
# f(d) returns the selected h (0, 1, 2 for h > 1) of every replicate
mc_of <- function(tbl, f) tbl %>%
  select(tree, hp, rv, n, ils, data) %>%
  mutate(data = lapply(data, function(d) tibble(h = f(d)))) %>% unnest(data) %>%
  filter(!is.na(h)) %>%
  mutate(h = factor(pmin(h, 2), 0:2, names(model_cols))) %>%
  count(tree, hp, rv, n, ils, h, name = "k") %>%
  group_by(tree, hp, rv, n, ils) %>%
  mutate(pct = 100 * k / sum(k)) %>% ungroup()

# best_k is "0", "1" or ">1" in the find_graphs summaries
k_of <- function(x) {
  x <- as.character(x); ifelse(x == "0", 0L, ifelse(x == "1", 1L, 2L))
}
mc_fg30 <- mc_of(fg, function(d) k_of(d$best_k))
mc_fg37 <- mc_of(fg, function(d) k_of(d$best_k_new_WR_3.7))
mc_sq <- mc_of(sq, function(d)
  ifelse(d$p_H0 > 0.05, 0L, ifelse(d$p_H1 > 0.05, 1L, 2L)))
# phynest: T* of the tree and ILS level added to each replicate first
mc_ph <- mc_of(ph %>% left_join(t_star, by = c("tree", "ils")) %>%
    rowwise() %>%
    mutate(data = list(mutate(data, ts = t_star))) %>%
    ungroup(),
  function(d) as.integer(d$T_per_site > d$ts))
# x = rate variation, with the reptile and the plant bar side by side
rv_x <- c("none" = 1, "gene" = 2, "lineage" = 3)
bar_off <- c("reptile" = -0.16, "plant" = 0.16)

# stacked bars of the selected h, one panel per dup/loss level and row
# fill_labels: legend labels of the fitted models (PhyNEST fits only h = 0
# and h = 1, so its orange bar stands for h >= 1)
bars <- function(mc, title, rows = "n", light_plant = FALSE,
                 fill_labels = names(model_cols)) {
  mc <- mc %>% mutate(
    x = rv_x[as.character(rv)] + bar_off[as.character(tree)],
    hp_lab = factor(hp_dl[as.character(hp)], hp_dl),
    row = if (rows == "n") paste0("ILS = ", ils, "\nn = ", n)
          else paste0("ILS = ", ils),
    a = if (light_plant) ifelse(tree == "plant", 0.55, 1) else 1)
  ggplot(mc, aes(x, pct, fill = h)) +
    geom_col(aes(alpha = a), width = 0.30, color = "grey20", linewidth = 0.3) +
    scale_alpha_identity() +
    scale_fill_manual(values = model_cols, name = "model choice",
                      labels = fill_labels, drop = FALSE,
                      breaks = names(model_cols)[seq_along(fill_labels)]) +
    scale_x_continuous(breaks = as.vector(outer(c(-0.16, 0.16), 1:3, "+")),
                       labels = rep(tree_short, 3),
                       sec.axis = dup_axis(breaks = 1:3, name = NULL,
                                           labels = c("none", "genes",
                                                      "lineages"))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    facet_grid(row ~ hp_lab) +
    labs(x = NULL, y = "percentage of replicates (%)", title = title) +
    theme_tt() + theme(panel.grid.major.x = element_blank(),
                       plot.title = element_text(face = "bold", size = BASE))
}

no_leg <- theme(legend.position = "none")
for (lp in switch(MODEL_CHOICE_PLANT, both = c(FALSE, TRUE), solid = FALSE,
                  light = TRUE)) {
  sfx <- if (lp && MODEL_CHOICE_PLANT == "both") "_lightplant" else ""
  p_fg <- plot_grid(
    bars(mc_fg30, "A  find_graphs model choice when WR \u2264 3.0",
         light_plant = lp) + no_leg,
    bars(mc_fg37, "B  find_graphs model choice when WR \u2264 3.7",
         light_plant = lp),
    ncol = 1, rel_heights = c(1, 1.07))
  save_both(p_fg, dE, paste0("model_choice_findgraphs_two_trees", sfx),
            13, 14)
  # each panel keeps its own legend: PhyNEST only fits h = 0 and h = 1
  p_sp <- plot_grid(
    bars(mc_sq, "A  SNaQ", light_plant = lp),
    bars(mc_ph, "B  PhyNEST (T > T*)", rows = "ils", light_plant = lp,
         fill_labels = c("h = 0 (tree)", "h \u2265 1 (network)")),
    ncol = 1, rel_heights = c(1.5, 1.1))
  save_both(p_sp, dE, paste0("model_choice_snaq_phynest_two_trees", sfx),
            13, 14)
}
write_csv(bind_rows(mc_fg30 %>% mutate(m = "find_graphs WR 3.0"),
                    mc_fg37 %>% mutate(m = "find_graphs WR 3.7"),
                    mc_sq %>% mutate(m = "SNaQ"),
                    mc_ph %>% mutate(m = "PhyNEST")),
          file.path(dE, "model_choice_counts.csv"))

#--- F. true tree recovery: h = 0, h = 1, h = 1 without the rogue taxon ---#
dF <- odir("F_true_tree_recovery")
pct <- function(x) 100 * mean(x, na.rm = TRUE)
rec_cols <- c("h = 0", "h = 1",
              "h = 1, rogue taxon removed\n(turtle, Amborella)")
# f0, f1, f1p: functions giving the per-replicate yes/no of the three columns
rec_of <- function(tbl, f0, f1, f1p, panel) tbl %>% rowwise() %>%
  mutate(h0 = pct(f0(data)), h1 = pct(f1(data)), h1p = pct(f1p(data))) %>%
  ungroup() %>% select(-data) %>%
  pivot_longer(c(h0, h1, h1p), names_to = "h") %>%
  mutate(h = factor(h, c("h0", "h1", "h1p"), rec_cols), panel = panel)
lg <- function(x) as.logical(x)
rec_fg_best <- rec_of(
  fg, function(d) lg(d$H0_best_tree_is_true_tree),
  function(d) lg(d$H1_best_graph_displayed_true_tree),
  function(d) lg(d$H1_best_graph_displayed_true_tree_noF),
  "A  find_graphs:\nbest-scoring graph")
rec_fg_all <- rec_of(
  fg, function(d) lg(d$true_tree_found_in_H0),
  function(d) lg(d$True_tree_displayed_H1_major) |
    lg(d$True_tree_displayed_H1_minor),
  function(d) lg(d$True_tree_noF_displayed_H1_major) |
    lg(d$True_tree_noF_displayed_H1_minor),
  "B  find_graphs:\nany retained graphs")
rec_sq <- rec_of(
  sq, function(d) d$RF_net0_true == 0,
  function(d) d$RF_net1_1_true == 0 | d$RF_net1_2_true == 0,
  function(d) d$RF_net1_1_true_noF == 0 | d$RF_net1_2_true_noF == 0,
  "A  SNaQ")
rec_ph <- rec_of(
  ph, function(d) d$RF_net0_true == 0,
  function(d) d$RF_net1_1_true == 0 | d$RF_net1_2_true == 0,
  function(d) d$RF_net1_1_true_noF == 0 | d$RF_net1_2_true_noF == 0,
  "B  PhyNEST")
rec_plot <- function(df) ggplot() + block_guides() + pts_setting(df, "value") +
  facet_grid(panel ~ h, scales = "free_y") +
  labs(y = "% of displaying the true tree") + theme_tt()
save_both(with_legend(rec_plot(bind_rows(rec_fg_best, rec_fg_all)),
                      legend_rv_ils_n()),
          dF, "true_tree_recovery_findgraphs_two_trees", 15, 10)
save_both(with_legend(rec_plot(bind_rows(rec_sq, rec_ph)), legend_rv_ils_n()),
          dF, "true_tree_recovery_snaq_phynest_two_trees", 15, 10)
write_csv(bind_rows(rec_fg_best, rec_fg_all, rec_sq, rec_ph) %>% select(-file),
          file.path(dF, "true_tree_recovery_per_setting.csv"))

#--- G. minor gamma under h = 1, reptile up and plant down ---#
dG <- odir("G_minor_gamma")
gam_of <- function(tbl, f, m) tbl %>%
  select(tree, hp, rv, n, ils, data) %>%
  mutate(data = lapply(data, function(d) tibble(gamma = f(d)))) %>%
  unnest(data) %>% filter(!is.na(gamma)) %>% mutate(m = m)
gam <- bind_rows(gam_of(fg, function(d) d$best_graph_gamma2, "find_graphs"),
                 gam_of(sq, function(d) d$gamma_2, "SNaQ"),
                 gam_of(ph, function(d) d$gamma_2, "PhyNEST"))

# counts per bin of width 0.01, plant counts negative
gh <- gam %>% mutate(bin = floor(pmin(gamma, 0.4999) / 0.01) * 0.01 + 0.005) %>%
  count(m, tree, hp, rv, n, ils, bin, name = "k") %>%
  mutate(k = ifelse(tree == "reptile", k, -k),
         row = paste0("ILS = ", ils, ", n = ", n),
         hp_lab = factor(hp_dl[as.character(hp)], hp_dl))

# every row gets its own y range (rule in visual_utilities.R)
# rows_n = FALSE for phynest: one individual per taxon only
gamma_hist <- function(gm, title, rows_n = TRUE) {
  if (!rows_n) gm <- gm %>% mutate(row = paste0("ILS = ", ils))
  ax <- nice_axes(gm, "row", "k")
  lim_df <- ax %>% select(row, limit) %>% crossing(sgn = c(1, -1)) %>%
    mutate(y = sgn * limit, x = 0.25)
  lab <- gm %>% distinct(row, hp_lab) %>%
    left_join(ax %>% select(row, limit), by = "row") %>%
    crossing(tree = c("reptile", "plant")) %>%
    mutate(y = ifelse(tree == "reptile", limit * 0.88, -limit * 0.88))
  ggplot(gm) +
    geom_blank(data = lim_df, aes(x, y)) +
    geom_col(aes(bin, k, fill = rv), width = 0.01, position = "identity",
             alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(data = lab, aes(0.49, y, label = tree), hjust = 1, size = 5.4) +
    facet_grid(row ~ hp_lab, scales = "free_y") +
    scale_fill_manual(values = rv_colors, labels = rv_labels,
                      name = "rate variation") +
    scale_y_continuous(labels = function(v) abs(v), breaks = mirror_breaks(ax),
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.4, 0.1),
                       expand = c(0.01, 0)) +
    labs(x = "estimated gene flow proportion (minor \u03b3, h = 1)",
         y = "number of replicates", title = title) +
    theme_tt() + theme(
      panel.grid.major.x = element_line(),
      plot.title = element_text(face = "bold", size = BASE),
      axis.text.x.bottom = element_text(angle = 0, hjust = 0.5),
      panel.spacing.x = unit(1, "lines"))
}
save_both(gamma_hist(gh %>% filter(m == "find_graphs"), "find_graphs"), dG,
          "minor_gamma_findgraphs_two_trees", 12.5, 12)
save_both(gamma_hist(gh %>% filter(m == "SNaQ"), "SNaQ"), dG,
          "minor_gamma_snaq_two_trees", 12.5, 12)
save_both(gamma_hist(gh %>% filter(m == "PhyNEST"),
                     "PhyNEST (one individual per taxon)", rows_n = FALSE),
          dG, "minor_gamma_phynest_two_trees", 12.5, 7.5)
#--- H. best h = 1 graph of find_graphs ---#
# minor gamma by whether its major tree is the true tree, and where its
# reticulation goes
dH <- odir("H_bestgraph")
bg <- load_method("findgraphs_bestgraph") %>%
  select(tree, hp, rv, n, ils, data) %>% unnest(data) %>%
  mutate(shown = factor(ifelse(rf_major == 0, "yes", "no"), c("yes", "no")),
         retic = paste(hybrid_taxon, "from", minor_donor))

# panel A: counts per bin of width 0.01, reptile up and plant down (Fig. 2)
bh <- bg %>%
  mutate(bin = floor(pmin(gamma_minor, 0.4999) / 0.01) * 0.01 + 0.005) %>%
  count(tree, rv, shown, bin, name = "k") %>%
  mutate(k = ifelse(tree == "reptile", k, -k))
axH <- nice_axes(bh, "rv", "k")
limH <- axH %>% select(rv, limit) %>% crossing(sgn = c(1, -1)) %>%
  mutate(y = sgn * limit, x = 0.25)
labH <- axH %>% crossing(tree = c("reptile", "plant")) %>%
  mutate(y = ifelse(tree == "reptile", limit * 0.88, -limit * 0.88))
pHa <- ggplot(bh) +
  geom_blank(data = limH, aes(x, y)) +
  geom_col(aes(bin, k, fill = shown), width = 0.01, position = "identity",
           alpha = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(data = labH, aes(0.49, y, label = tree), hjust = 1, size = 5.4) +
  facet_grid(rv ~ ., labeller = labeller(rv = rv_labels), scales = "free_y") +
  scale_fill_manual(values = c(yes = "#228B22", no = "#DC143C"),
                    name = "major tree of the best graph = true tree?") +
  scale_y_continuous(labels = function(v) abs(v), breaks = mirror_breaks(axH),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     expand = c(0.01, 0)) +
  labs(x = "minor \u03b3 of the best graph (h = 1)", y = "number of replicates",
       title = "A  estimated gene flow") +
  theme_tt() + theme(
    legend.position = "bottom", panel.grid.major.x = element_line(),
    plot.title = element_text(face = "bold", size = BASE),
    axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))

# panel B: the reticulation of the best graph (hybrid taxon and its minor
# parent), the six most frequent ones named, the rest pooled as "other"
top_retic <- bg %>% count(retic, sort = TRUE) %>% slice_head(n = 6) %>%
  pull(retic)
bp <- bg %>%
  mutate(retic = factor(ifelse(retic %in% top_retic, retic, "other"),
                        c(top_retic, "other"))) %>%
  count(tree, rv, retic, name = "k") %>%
  group_by(tree, rv) %>% mutate(pct = 100 * k / sum(k)) %>% ungroup() %>%
  mutate(x = rv_x[as.character(rv)] + bar_off[as.character(tree)])
pHb <- ggplot(bp, aes(x, pct, fill = retic)) +
  geom_col(width = 0.30, color = "grey20", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2",
                    name = "reticulation: hybrid taxon, from minor parent") +
  scale_x_continuous(breaks = as.vector(outer(c(-0.16, 0.16), 1:3, "+")),
                     labels = rep(tree_short, 3),
                     sec.axis = dup_axis(breaks = 1:3, name = NULL,
                                         labels = c("none", "genes",
                                                    "lineages"))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  guides(fill = guide_legend(ncol = 3, title.position = "top")) +
  labs(x = NULL, y = "percentage of replicates (%)",
       title = "B  reticulation of the best graph") +
  theme_tt() + theme(legend.position = "bottom",
                     panel.grid.major.x = element_blank(),
                     plot.title = element_text(face = "bold", size = BASE))


                  
save_both(plot_grid(pHa, pHb, ncol = 2, rel_widths = c(1, 1)), dH,
          "bestgraph_two_trees", 16, 9.5)
cat("appendix figures written to", out_root, "\n")
