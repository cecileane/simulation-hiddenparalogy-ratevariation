#!/usr/bin/env Rscript
# ============================================================================
# scripts/visuals/visual_main_figures.R
#
# Main-text figures of the two-tree comparison (reptile and angiosperm
# species tree in the same figure):
#   Fig. 1  four panels, one point per setting, both trees and both
#           numbers of individuals per taxon
#   Fig. 2  mirrored histograms of the estimated minor gamma under h = 1
# Inputs  : per-replicate summaries of the kappa-corrected rerun, read with
#           load_method() from visual_utilities.R
# Outputs : visualization_results/main_figures/
#             fig1_two_trees.pdf/.png
#             fig2_gamma_mirrored_overlap.pdf/.png
#             per_setting.csv  (the numbers behind Fig. 1)
#             T_star.csv
# Usage   : Rscript scripts/visuals/visual_main_figures.R   (repo root)
#
# Fig. 1: x = hidden paralogy level (none / medium / high), each
# level split into a reptile and a plant sub-column; the plant sub-column
# sits on a very light grey band and solid lines separate the levels.
# color = substitution rate variation, shape = ILS (circle high, triangle
# low), fill = individuals per taxon (filled n = 1, hollow n = 2).
# Panels: A 95th percentile of the worst residual of the true tree
# (find_graphs), B 95th percentile of T relative to T* (PhyNEST), C % of
# replicates whose h = 1 network displays the true tree, D type I error.
# ============================================================================

library(ggplot2)
library(cowplot)
library(scales)
source("scripts/visuals/visual_utilities.R")

out_dir <- file.path(results_root, "main_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fg <- load_method("findgraphs")
sq <- load_method("snaq")
ph <- load_method("phynest")
t_star <- t_star_from(ph)

# --------------------------------------------------------------
# Per-setting summaries
# --------------------------------------------------------------
pct <- function(x) 100 * mean(x, na.rm = TRUE)

summ_fg <- function(tbl) tbl %>% rowwise() %>% mutate(
  wr95 = quantile(data$true_tree_wr, 0.95, na.rm = TRUE),
  typeI_wr30 = pct(as.character(data$best_k) != "0"),
  typeI_wr37 = pct(as.character(data$best_k_new_WR_3.7) != "0"),
  rec_h0 = pct(as.logical(data$H0_best_tree_is_true_tree)),
  rec_h1 = pct(as.logical(data$H1_best_graph_displayed_true_tree)),
  rec_h1_pruned = pct(as.logical(data$H1_best_graph_displayed_true_tree_noF)),
  n_reps = nrow(data)) %>% ungroup()

summ_sq <- function(tbl) tbl %>% rowwise() %>% mutate(
  typeI = pct(data$p_H0 <= 0.05),
  rec_h0 = pct(data$RF_net0_true == 0),
  rec_h1 = pct(data$RF_net1_1_true == 0 | data$RF_net1_2_true == 0),
  rec_h1_pruned = pct(data$RF_net1_1_true_noF == 0 |
                        data$RF_net1_2_true_noF == 0),
  n_reps = nrow(data)) %>% ungroup()

summ_ph <- function(tbl) tbl %>% left_join(t_star, by = c("tree", "ils")) %>%
  rowwise() %>% mutate(
    t95 = quantile(data$T_per_site, 0.95),
    typeI = pct(data$T_per_site > t_star),
    rec_h0 = pct(data$RF_net0_true == 0),
    rec_h1 = pct(data$RF_net1_1_true == 0 | data$RF_net1_2_true == 0),
    rec_h1_pruned = pct(data$RF_net1_1_true_noF == 0 |
                          data$RF_net1_2_true_noF == 0),
    n_reps = nrow(data)) %>% ungroup()

agg <- bind_rows(summ_fg(fg) %>% select(-data), summ_sq(sq) %>% select(-data),
                 summ_ph(ph) %>% select(-data))

write_csv(agg %>% select(-file), file.path(out_dir, "per_setting.csv"))
write_csv(t_star, file.path(out_dir, "T_star.csv"))

cat("settings loaded (method x tree):\n")
print(agg %>% count(method, tree) %>% as.data.frame())
cat("T* per tree and ILS:\n")
print(as.data.frame(t_star))

# -------------------------------------------------------------
# Fig. 1
# -------------------------------------------------------------
BASE <- 15    # base font size of the panels
PT   <- 3.1  # point size
LAB  <- 20  # panel letter size
shape_ils <- c("high" = 21, "low" = 24)
tree_x <- as.vector(outer(c(-0.2, 0.2), 1:3, "+"))   # the six sub-columns
hp_edges <- c(0.45, 1.5, 2.5, 3.55)

# x position of a setting: the level, then the tree, the rate variation and
# n shift the point so that no two settings share an x
x_of <- function(df) as.numeric(df$hp) +
  ifelse(df$tree == "reptile", -0.2, 0.2) +
  c(none = -0.085, gene = 0, lineage = 0.085)[as.character(df$rv)] +
  ifelse(df$n == 1, -0.035, 0.035)

# point layer plus its scales: 
# color = rate variation, shape = ILS, fill = n 
pts <- function(df, y) {
  df <- df %>% mutate(x = x_of(df), shp = shape_ils[as.character(ils)],
                      fillc = ifelse(n == 1, rv_colors[as.character(rv)],
                                     "white"),
                      alpha = ifelse(corrected, 1, 0.4))
  list(geom_point(data = df, aes(x = x, y = .data[[y]], shape = shp,
                                 color = rv, fill = fillc, alpha = alpha),
                  size = PT, stroke = 1.2),
       scale_shape_identity(), scale_fill_identity(), scale_alpha_identity(),
       scale_color_manual(values = rv_colors, labels = rv_labels,
                          name = "rate variation"),
       scale_x_continuous(breaks = tree_x, labels = rep(tree_short, 3),
                          limits = c(0.45, 3.55), minor_breaks = NULL,
                          sec.axis = dup_axis(breaks = 1:3, labels = hp_levels,
                                              name = NULL)),
       labs(x = NULL))
}

# background: the plant half of each level on a very light grey band, and
# a solid line between the levels. ylo must be finite on a log axis.
guides_x <- function(ylo = -Inf) list(
  geom_rect(data = tibble(xmin = 1:3, xmax = hp_edges[2:4]),
            aes(xmin = xmin, xmax = xmax, ymin = ylo, ymax = Inf),
            inherit.aes = FALSE, fill = "grey94", alpha = 0.5),
  geom_vline(xintercept = c(1.5, 2.5), color = "grey55", linewidth = 0.7))

theme_f1 <- function(base = BASE) theme_bw(base) + theme(
  # strip on top, the level labels between the strip and the panel
  strip.placement = "outside",
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "grey78", linewidth = 0.55),
  strip.background = element_rect(fill = "grey90"),
  strip.text = element_text(face = "bold", size = base),
  axis.text = element_text(color = "black", size = base - 2),
  axis.title.y = element_text(size = base - 1),
  legend.position = "none",
  axis.text.x.bottom = element_text(size = base - 2, angle = 45, hjust = 1,
                                    vjust = 1),
  axis.text.x.top = element_text(size = base - 2, face = "bold"),
  axis.ticks.x.top = element_blank(),
  plot.margin = margin(4, 6, 2, 4))

# The fill of the points is an identity scale, so the legend is built from
# a small dummy data set with the same encoding
legend_f1 <- function() {
  d <- expand_grid(rv = factor(names(rv_colors), names(rv_colors)),
                   ils = factor(c("high", "low"), c("high", "low")),
                   n = factor(c(1, 2)))
  th <- theme_bw(BASE) + theme(
    legend.position = "bottom", legend.box = "horizontal",
    legend.text = element_text(size = BASE - 1),
    legend.title = element_text(size = BASE - 1, face = "bold"),
    legend.background = element_rect(color = "black"),
    legend.key.size = unit(1.1, "lines"))
  p <- ggplot(d, aes(1, 1)) +
    geom_point(aes(color = rv), shape = 21, size = 4.5, stroke = 1.2,
               fill = NA) +
    scale_color_manual(values = rv_colors, labels = rv_labels,
                       name = "rate variation",
                       guide = guide_legend(order = 1, override.aes = list(
                         shape = 21, fill = rv_colors, color = rv_colors))) +
    geom_point(aes(shape = ils), size = 4.5, stroke = 1.2, color = "black",
               fill = "black") +
    scale_shape_manual(values = shape_ils,
                       labels = c(high = "high ILS", low = "low ILS"),
                       name = "ILS", guide = guide_legend(order = 2)) +
    geom_point(aes(fill = n), shape = 21, size = 4.5, stroke = 1.2,
               color = "black") +
    scale_fill_manual(values = c("1" = "black", "2" = "white"),
                      labels = c("n = 1", "n = 2"),
                      name = "individuals / taxon",
                      guide = guide_legend(order = 3))
  cowplot::get_legend(p + th)
}

afg <- agg %>% filter(method == "findgraphs")
asq <- agg %>% filter(method == "snaq")
aph <- agg %>% filter(method == "phynest") %>% mutate(t95_rel = t95 / t_star)

# A: 95th percentile of the worst residuals, with the two WR thresholds
pA <- ggplot() + guides_x(ylo = min(afg$wr95, na.rm = TRUE) * 0.9) +
  pts(afg, "wr95") +
  geom_hline(yintercept = c(3.0, 3.7), linetype = "dashed",
             color = c("red", "green4"), linewidth = 0.6) +
  scale_y_log10(breaks = c(2, 3, 5, 10, 20, 30)) +
  facet_wrap(~"find_graphs") +
  labs(y = "95th pct. of worst residuals (log)") + theme_f1()
# B: 95th percentile of the per-site T, with T* of each tree (solid at
# high ILS, dashed at low ILS), lines colored by tree
tsl <- aph %>% distinct(tree, ils, t_star)
pB <- ggplot() + guides_x(ylo = min(aph$t95) * 0.9) + pts(aph, "t95") +
  geom_hline(data = tsl, aes(yintercept = t_star, linetype = ils),
             color = tstar_colors[as.character(tsl$tree)], linewidth = 0.7,
             show.legend = FALSE) +
  scale_linetype_manual(values = tstar_types) +
  scale_y_log10(breaks = c(0.01, 0.02, 0.03, 0.05, 0.1, 0.15)) +
  facet_wrap(~"PhyNEST") +
  labs(y = "95th pct. of per-site T (log)") + theme_f1()
# C: true tree displayed by the h = 1 network
rec <- bind_rows(afg %>% mutate(m = "find_graphs"), asq %>% mutate(m = "SNaQ"),
                 aph %>% mutate(m = "PhyNEST")) %>%
  mutate(m = factor(m, c("find_graphs", "SNaQ", "PhyNEST")))
pC <- ggplot() + guides_x() + pts(rec, "rec_h1") +
  scale_y_continuous(limits = c(0, 102), expand = c(0, 0)) +
  facet_wrap(~m, nrow = 1) +
  labs(y = "% of true tree displayed") + theme_f1()
# D: type I error rate
keep <- c("hp", "rv", "ils", "tree", "n", "corrected")
m_t1 <- c("find_graphs (WR ≤ 3.0)", "find_graphs (WR ≤ 3.7)", "SNaQ",
          "PhyNEST (T > T*)")
t1 <- bind_rows(
  afg %>% select(all_of(keep), v = typeI_wr30) %>% mutate(m = m_t1[1]),
  afg %>% select(all_of(keep), v = typeI_wr37) %>% mutate(m = m_t1[2]),
  asq %>% select(all_of(keep), v = typeI) %>% mutate(m = m_t1[3]),
  aph %>% select(all_of(keep), v = typeI) %>% mutate(m = m_t1[4])) %>%
  mutate(m = factor(m, m_t1))
pD <- ggplot() + guides_x() + pts(t1, "v") +
  scale_y_continuous(limits = c(-3, 103), expand = c(0, 0)) +
  facet_wrap(~m, nrow = 1) +
  labs(y = "type I error rate (%)") + theme_f1()

top <- plot_grid(pA, pB, pC, nrow = 1, rel_widths = c(1, 1, 2.6),
                 labels = c("A", "B", "C"), label_size = LAB, align = "h",
                 axis = "tb")
fig1 <- plot_grid(plot_grid(top, pD, ncol = 1, labels = c("", "D"),
                            label_size = LAB, rel_heights = c(1, 1)),
                  legend_f1(), ncol = 1, rel_heights = c(1, 0.09))
ggsave(file.path(out_dir, "fig1_two_trees.pdf"), fig1, width = 14,
       height = 9.5, device = cairo_pdf)
ggsave(file.path(out_dir, "fig1_two_trees.png"), fig1, width = 14,
       height = 9.5, dpi = 130, bg = "white")

# -------------------------------------------------------------------
# Fig. 2: minor gamma under h = 1, reptile up and plant down
# -------------------------------------------------------------------
# one row per rep of the n = 1 settings: 
# the minor gamma of the h = 1 network and whether it displays the true tree
gamma_rows <- function(tbl, m, f) tbl %>% filter(n == 1) %>%
  select(tree, rv, corrected, data) %>%
  mutate(data = lapply(data, f)) %>% unnest(data) %>% mutate(m = m)
gam <- bind_rows(
  gamma_rows(fg, "find_graphs", function(d) tibble(
    gamma = d$best_graph_gamma2,
    shown = as.logical(d$H1_best_graph_displayed_true_tree))),
  gamma_rows(sq, "SNaQ", function(d) tibble(
    gamma = d$gamma_2,
    shown = d$RF_net1_1_true == 0 | d$RF_net1_2_true == 0)),
  gamma_rows(ph, "PhyNEST", function(d) tibble(
    gamma = d$gamma_2,
    shown = d$RF_net1_1_true == 0 | d$RF_net1_2_true == 0))) %>%
  mutate(m = factor(m, c("find_graphs", "SNaQ", "PhyNEST")),
         shown = factor(ifelse(shown, "yes", "no"), c("yes", "no")))

# counts per bin of width 0.01; the plant counts are negative so that they
# hang below the zero line. 
# The two trees do not have the same number of
# finished settings, which the caption should say.
gh <- gam %>% mutate(bin = floor(pmin(gamma, 0.4999) / 0.01) * 0.01 + 0.005) %>%
  count(m, rv, tree, shown, bin) %>%
  mutate(k = ifelse(tree == "reptile", n, -n))
rv_strip <- c(none = "no rate variation", gene = "variation across genes",
              lineage = "variation across lineages")

# y range of each row (rate variation) from the rule in visual_utilities.R;
# facet_grid cannot take per-row breaks, so the breaks function recognises
# the row by the panel limit (hence expand = 0 below)
ax <- nice_axes(gh, "rv", "k") %>%
  mutate(rv = factor(rv, names(rv_colors))) %>% arrange(rv)
lim_df <- ax %>% select(rv, limit) %>% crossing(sgn = c(1, -1)) %>%
  mutate(y = sgn * limit, x = 0.25)
lab_df <- lim_df %>%
  mutate(y = y * 0.85, label = ifelse(sgn > 0, "reptile", "plant"), x = 0.49)
cat("Fig. 2 row limits:", paste(ax$rv, ax$limit, ax$step, sep = "/"), "\n")
gh <- gh %>% left_join(ax %>% select(rv, limit), by = "rv") %>%
  mutate(ymin = pmin(0, pmax(k, -limit)), ymax = pmax(0, pmin(k, limit)))

# the "yes" and "no" bars both start at zero and overlap with some alpha
fig2 <- ggplot(gh) +
  geom_blank(data = lim_df, aes(x, y)) +
  geom_rect(aes(xmin = bin - 0.005, xmax = bin + 0.005, ymin = ymin,
                ymax = ymax, fill = shown), alpha = 0.75) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_text(data = lab_df, aes(x, y, label = label), hjust = 1, size = 4.5) +
  facet_grid(rv ~ m, labeller = labeller(rv = rv_strip), scales = "free_y") +
  scale_fill_manual(
    values = c(yes = "#228B22", no = "#DC143C"),
    name = "is the true tree displayed in the estimated network?") +
  scale_y_continuous(labels = function(v) abs(v), expand = c(0, 0),
                     breaks = mirror_breaks(ax)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     expand = c(0.01, 0)) +
  labs(x = "estimated gene flow proportion", y = "number of replicates") +
  theme_bw(15) + theme(
    legend.position = "bottom", panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 13),
    strip.text.y = element_text(size = 12),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.box.spacing = unit(0.2, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(1.4, "lines"),
    panel.spacing.y = unit(0.7, "lines"))

  
ggsave(file.path(out_dir, "fig2_gamma_mirrored_overlap.pdf"), fig2,
       width = 8.5, height = 8.5, device = cairo_pdf)

ggsave(file.path(out_dir, "fig2_gamma_mirrored_overlap.png"), fig2,
       width = 8.5, height = 8.5, dpi = 200, bg = "white")

cat("main figures written to", out_dir, "\n")
