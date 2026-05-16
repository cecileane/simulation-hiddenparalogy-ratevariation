#!/usr/bin/env Rscript
# ============================================================================
# scripts/visual_utilities.R
#
# Purpose : Shared R plotting helpers used by the summary_*.jl scripts and by
#           the Quarto notebooks in visualization_scripts/. Defines color
#           palettes, parameter parsers (DUP*-LOS*-RV*-N_ind*-SF*-genelen*),
#           and per-quantity plotting functions (RF distance distributions,
#           worst-residual histograms, gamma scatter plots, etc.).
# Inputs  : Sourced (`source("scripts/visual_utilities.R")`) from other code;
#           individual plot functions read CSV files from results/, snaq_summary/
#           and findgraph_summary/.
# Outputs : PNG / PDF figures written under visualization_results/.
# Usage   : Not run directly. Source from another R / Julia script.
# Note    : `plot_taxon_recovery_heatmaps()` is LEGACY and gated by stop()
#           (see banner above that function).
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(readr)
library(tibble)

# Global color palette for substitution rate variation
rv_colors <- c("none" = "#898686",
               "gene" = "#039afe",
               "lineage" = "#fbc531")

# Global color palette for tree display status
# (H=1 displays vs not displays true tree)
color_displayed_tree <- c("True tree displayed" = "#228B22",
                          "True tree NOT displayed" = "#DC143C",
                          "H=1 displays true tree" = "#228B22",
                          "H=1 does NOT display true tree" = "#DC143C")

#' Plot overlapping ratevar distributions using facet_grid
#' 
#' Creates a 12-panel plot with facet_grid where:
#' - Rows represent (n_inds, ILS) combinations
#' - Columns represent duplication/loss rates
#' - Colors represent substitution rate variations (gene, lineage, none)
#'
#' @param input_dir Directory containing CSV files
#' @param column_name_1 First gamma column name (e.g., "best_graph_gamma1")
#' @param column_name_2 Second gamma column name (e.g., "best_graph_gamma2")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param max_y Maximum y-axis value for all plots (default: 100)
plot_overlapping_ratevar_by_n_inds_sf <- function(input_dir, 
                                                   column_name_1, 
                                                   column_name_2,
                                                   output_dir, 
                                                   output_filename,
                                                   max_y = 70,
                                                   alpha = 0.6, 
                                                   rv_colors = get(
                                                     "rv_colors",
                                                     envir = globalenv()),
                                                   plot_title = NULL) {
  
  # Print max_y for debugging
  cat("Using max_y =", max_y, "\n")
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Read and combine all data
  all_data <- data.frame()
  
  for (csv_file in csv_files) {
    # Extract parameters from filename
    filename <- basename(csv_file)
    
    n_inds_match <- str_match(filename, "N_ind(\\d+)")
    sf_match <- str_match(filename, "SF([\\d.]+)")
    rate_match <- str_match(filename, "DUP([\\d.e-]+)-LOS([\\d.e-]+)")
    rv_match <- str_match(filename, "-(RVG|RVL|RVN|RVGL)-")
    
    if (any(is.na(c(n_inds_match, sf_match, rate_match, rv_match)))) {
      next
    }
    
    # Read the CSV file
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      if (!(column_name_1 %in% names(df) && column_name_2 %in% names(df))) {
        next
      }
      
      # Calculate minor gamma
      df$minor_gamma <- pmin(df[[column_name_1]], df[[column_name_2]])
      df <- df[!is.na(df$minor_gamma), ]
      
      if (nrow(df) == 0) next
      
      # Add grouping variables
      df$n_inds <- n_inds_match[2]
      df$sf <- sf_match[2]
      df$rate <- rate_match[2]
      df$ratevar <- rv_match[2]
      
      all_data <- rbind(all_data,
        df[, c("minor_gamma", "n_inds", "sf", "rate", "ratevar")])
    }, error = function(e) {
      cat("Error reading", csv_file, ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(all_data) == 0) {
    cat("No valid data found\n")
    return(invisible(NULL))
  }
  
  # Map codes to readable names
  all_data$ratevar_name <- recode(all_data$ratevar,
                                  "RVG" = "gene",
                                  "RVL" = "lineage", 
                                  "RVN" = "none")
  
  all_data$ils_level <- recode(all_data$sf,
                               "1.0" = "high",
                               "0.5" = "low")
  
  # Create faceting variables
  all_data$row_facet <- paste0(
    "individuals / taxon = ", all_data$n_inds,
    "\nILS = ", all_data$ils_level)
  
  # Format rate in scientific notation for column facets
  all_data$rate_sci <- sapply(all_data$rate, function(r) {
    rate_num <- as.numeric(r)
    if (rate_num == 0) {
      "0.0"
    } else {
      # Convert to scientific notation (e.g., 0.0003 -> 3e-4)
      formatted <- format(rate_num, scientific = TRUE, digits = 1)
      # Clean up the format (remove + sign, simplify)
      formatted <- gsub("e\\+00", "", formatted)
      formatted <- gsub("e-0", "e-", formatted)
      formatted
    }
  })
  all_data$col_facet <- paste0("dup/loss rate = ", all_data$rate_sci)
  
  # Order facets properly
  all_data$row_facet <- factor(all_data$row_facet,
    levels = c("individuals / taxon = 1\nILS = high",
               "individuals / taxon = 1\nILS = low",
               "individuals / taxon = 2\nILS = high",
               "individuals / taxon = 2\nILS = low"))
  
  all_data$col_facet <- factor(all_data$col_facet,
                              levels = c("dup/loss rate = 0.0",
                                       "dup/loss rate = 3e-4",
                                       "dup/loss rate = 4e-4"))
  
  all_data$ratevar_name <- factor(all_data$ratevar_name,
                                 levels = c("gene", "lineage", "none"))
  
  # Define colors (more distinctive: blue, orange, dark gray)
  colors <- rv_colors
  
  # Create the plot with facet_grid
  p <- ggplot(all_data, aes(x = minor_gamma, fill = ratevar_name)) +
    geom_histogram(binwidth = 0.01, alpha = alpha, position = "identity") +
    facet_grid(row_facet ~ col_facet, scales = "fixed") +
    scale_fill_manual(values = colors,
                     name = "rate variation:",
                     labels = c("across genes",
                              "across lineages",
                              "none")) +
    scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, max_y)) +
    labs(x = "estimated gene flow proportion",
         y = "number of replicates",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 21, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 21, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 20, color = "black", face = "bold", angle = 45, hjust = 1),
      axis.text.y       = element_text(
        size = 20, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 18, face = "bold"),
      strip.text.x      = element_text(color = "black"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "bottom",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 18, face = "bold"),
      legend.title       = element_text(size = 18, face = "bold"),
      legend.key.size    = unit(1.5, "cm"),
      panel.spacing.x    = unit(1.0, "lines"),
      panel.spacing.y    = unit(1.0, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, "_combined.png"))
  ggsave(output_path, plot = p,
    width = 14, height = 16, dpi = 300, units = "in")
  
  cat("Saved combined 12-panel figure:", output_path, "\n")
  
  invisible(p)
}

#' Plot hypothesis acceptance using facet_grid
#' 
#' Creates a multi-panel plot using facet_grid where:
#' - Rows represent n_ind values (1, 2)  
#' - Columns represent ILS levels (low SF=0.5, high SF=1.0) x dup/loss rates
#' - X-axis shows rate variation (RV: N, G, L, GL)
#' - Stacked bars show hypothesis acceptance (H=0, H=1, H>1)
#' - Works for both SNaQ and findgraphs data
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind, 
#'   and hypothesis columns (pct_H0/H0Accepted, pct_H1/H1Accepted,
#'   pct_BT1/BT1Accepted or pct_H_greater/H>1Accepted)
#' @param value_col_prefix Prefix for value columns: "pct_" for percentages
#'   or "" for counts
#' @param h_greater_col Name of the H>1 column: "pct_H_greater", "pct_BT1",
#'   "H>1Accepted", or "BT1Accepted"
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @param y_label Y-axis label (default: "Percentage of replicates")
#' @return ggplot object
plot_hypothesis_acceptance_grid <- function(df,
                                 value_col_prefix = "pct_",
                                 h0_col = NULL,
                                 h1_col = NULL,
                                 h_greater_col = "pct_H_greater",
                                 output_dir,
                                 output_filename,
                                 plot_title = NULL,
                                 y_label = "Percentage of replicates (%)") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Standardize column names
  df_work <- df %>%
    mutate(
      n_ind = if("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L", "GL"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF))
  
  # Determine column names based on prefix if not explicitly provided
  if (is.null(h0_col)) {
    h0_col <- paste0(value_col_prefix, "H0")
  }
  if (is.null(h1_col)) {
    h1_col <- paste0(value_col_prefix, "H1")
  }
  
  # Handle alternative column names for H0 and H1
  if (!h0_col %in% names(df_work)) {
    if ("H0Accepted" %in% names(df_work)) h0_col <- "H0Accepted"
    else if (paste0(value_col_prefix, "H=0") %in% names(df_work))
      h0_col <- paste0(value_col_prefix, "H=0")
    else if ("H=0Accepted" %in% names(df_work)) h0_col <- "H=0Accepted"
  }
  if (!h1_col %in% names(df_work)) {
    if ("H1Accepted" %in% names(df_work)) h1_col <- "H1Accepted"
    else if (paste0(value_col_prefix, "H=1") %in% names(df_work))
      h1_col <- paste0(value_col_prefix, "H=1")
    else if ("H=1Accepted" %in% names(df_work)) h1_col <- "H=1Accepted"
  }
  
  # Prepare plot data
  plot_data <- df_work %>%
    select(dup_loss_rate, RV, SF, n_ind, 
           H0 = all_of(h0_col), 
           H1 = all_of(h1_col),
           H_greater = all_of(h_greater_col)) %>%
    pivot_longer(cols = c(H0, H1, H_greater),
                 names_to = "hypothesis",
                 values_to = "value") %>%
    mutate(
      hypothesis = factor(hypothesis,
                         levels = c("H0", "H1", "H_greater"),
                         labels = c("h=0 (tree)", "h=1 (network)",
                                   "h>1 (network)")),

      RV = droplevels(RV),
      ILS = ifelse(SF == 1.0, "high", "low"),
      n_ind_label = paste0("individuals / taxon = ", n_ind),
      dup_rate_label = paste0("dup/loss rate = ", dup_loss_rate),
      facet_col = paste0("ILS = ", ILS, "\n", dup_rate_label)
    )
  
  # Order facets properly
  plot_data <- plot_data %>%
    mutate(
      n_ind_label = factor(
        n_ind_label, levels = unique(n_ind_label[order(n_ind)])),
      facet_col = factor(facet_col)
    )

  # Create plot title if not provided
  if (is.null(plot_title)) {
    if (grepl("../visualization_results/findgraph", output_dir, fixed = TRUE)) {
      plot_title <- "find_graphs model choice by parameter settings"
    } else if (grepl(
        "../visualization_results/snaq", output_dir, fixed = TRUE)) {
      plot_title <- "SNaQ model choice by parameter settings"
    } else {
      plot_title <- "model choice by parameter settings"
    }
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = RV, y = value, fill = hypothesis)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_grid(n_ind_label ~ facet_col) +
    scale_fill_manual(values = c("h=0 (tree)" = "#2166ac",
                                 "h=1 (network)" = "#fdae61",
                                 "h>1 (network)" = "#d73027")) +
    scale_x_discrete(labels = c("N" = "none", "G" = "gene",
                                "L" = "lineage", "GL" = "gene+lineage")) +
    labs(x = "substitution rate variation",
         y = y_label,
         fill = "model choice :",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 16, color = "black", face = "bold", angle = 45, hjust = 1),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 14, face = "bold"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey60", linewidth = 0.7, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey82", linewidth = 0.4, linetype = "dotted"),
      legend.position    = "bottom",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.title       = element_text(size = 16, face = "bold"),
      legend.text        = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(0.5, "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, output_filename)
  ggsave(output_path, plot = p, width = 14.2, height = 8, units = "in")
  
  cat("Saved hypothesis acceptance grid plot:", output_path, "\n")
  
  invisible(p)
}

#' Plot true graph recovery rates by duplication/loss rate
#' 
#' Creates a multi-panel line plot using facet_grid where:
#' - Rows represent n_ind values (1, 2)
#' - Columns represent RV (N, G, L, GL)
#' - X-axis shows duplication/loss rate
#' - Y-axis shows percentage of replicates finding true graph
#' - Lines for H=0 (find_true_net0) and H=1
#'   (find_true_net1 and find_true_net1_noF)
#' - Works for both SNaQ and findgraphs data
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind
#' @param h0_col Column name for H=0 true graph counts
#'   (default: "find_true_net0")
#' @param h1_col Column name for H=1 true graph counts
#'   (default: "find_true_net1")
#' @param h1_noF_col Column name for H=1 true graph without F counts
#'   (default: "find_true_net1_noF")
#' @param n_ind_filter Filter to specific n_ind value (optional)
#' @param sf_filter Filter to specific SF value (optional)
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @return ggplot object
plot_true_graph_recovery_lines <- function(df,
                                          h0_col = "find_true_net0",
                                          h1_col = "find_true_net1",
                                          h1_noF_col = "find_true_net1_noF",
                                          n_ind_filter = NULL,
                                          sf_filter = NULL,
                                          output_dir,
                                          output_filename,
                                          plot_title = NULL) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Standardize column names
  df_work <- df %>%
    mutate(
      n_ind = if("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF))
  
  # Apply filters if provided
  if (!is.null(n_ind_filter)) {
    df_work <- df_work %>% filter(n_ind == n_ind_filter)
  }
  if (!is.null(sf_filter)) {
    df_work <- df_work %>% filter(SF == sf_filter)
  }
  
  # Prepare plot data
  plot_data <- df_work %>%
    select(dup_loss_rate, RV, SF, n_ind, 
           H0 = all_of(h0_col), 
           H1 = all_of(h1_col),
           H1_noF = all_of(h1_noF_col)) %>%
    pivot_longer(cols = c(H0, H1, H1_noF),
                 names_to = "metric",
                 values_to = "percentage") %>%
    mutate(
      metric = factor(metric,
                     levels = c("H0", "H1", "H1_noF"),
                     labels = c("H=0 finds true graph", 
                               "H=1 finds true graph", 
                               "H=1 finds true graph (no F)")),
      dup_loss_rate_cat = factor(dup_loss_rate),
      RV = droplevels(RV),
      RV_label = recode(RV,
                       "N" = "substitution variation = none",
                       "G" = "substitution variation = gene",
                       "L" = "substitution variation = lineage"),
      ILS = ifelse(SF == 1.0, "high", "low"),
      n_ind_label = paste0("individuals / taxon = ", n_ind),
      facet_col = paste0("ILS = ", ILS)
    )
  
  # Order facets properly
  plot_data <- plot_data %>%
    mutate(
      n_ind_label = factor(
        n_ind_label, levels = unique(n_ind_label[order(n_ind)])),
      facet_col = factor(facet_col, levels = c("ILS = low", "ILS = high"))
    )
  
  # Create plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- "True Graph Recovery Rate by Duplication/Loss Rate"
  }
  
  # Define colors and linetypes
  colors <- c("H=0 finds true graph" = "#2C7FB8",
             "H=1 finds true graph" = "#C51B7D",
             "H=1 finds true graph (no F)" = "#999999")
  
  linetypes <- c("H=0 finds true graph" = "solid",
                "H=1 finds true graph" = "solid",
                "H=1 finds true graph (no F)" = "dashed")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = dup_loss_rate_cat, y = percentage,
                            color = metric, linetype = metric,
                            group = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_grid(n_ind_label ~ RV_label + facet_col, scales = "free_x") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_y_continuous(limits = c(-1, 105), expand = c(0, 0)) +
    labs(x = "duplication/loss rate",
         y = "number of replicates",
         color = NULL,
         linetype = NULL,
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 14, face = "bold"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey60", linewidth = 0.7, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey82", linewidth = 0.4, linetype = "dotted"),
      legend.position    = "bottom",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(0.5, "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, output_filename)
  ggsave(output_path, plot = p, width = 16, height = 8, units = "in")
  
  cat("Saved true graph recovery line plot:", output_path, "\n")
  
  invisible(p)
}

#' Plot true graph recovery rates aggregated across RV
#'
#' Single-panel jitter plot (per metric facet) where:
#' - Columns represent metric (H=0, H=1, H=1 no F)
#' - Color encodes substitution rate variation (none/gene/lineage)
#' - Shape + fill encode ILS × n_ind (4 combos):
#'     n_ind=1, high ILS → filled circle   (shape 21)
#'     n_ind=1, low  ILS → filled triangle (shape 24)
#'     n_ind=2, high ILS → hollow circle   (shape 21)
#'     n_ind=2, low  ILS → hollow triangle (shape 24)
#' - No connecting lines; points jittered horizontally
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind
#' @param h0_col   Column for H=0 true graph counts (default: "find_true_net0")
#' @param h1_col   Column for H=1 true graph counts (default: "find_true_net1")
#' @param h1_noF_col Column for H=1 no-F counts (default: "find_true_net1_noF")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @param y_range  Numeric vector of length 2 for y-axis limits
#'   (default: c(75, 105))
#' @param jitter_width Horizontal jitter width (default: 0.2)
#' @param rv_colors Named color vector for rate variation levels
#' @return ggplot object invisibly
plot_true_graph_recovery_lines_aggregated <- function(df,
                    h0_col = "find_true_net0",
                    h1_col = "find_true_net1",
                    h1_noF_col = "find_true_net1_noF",
                    output_dir,
                    output_filename,
                    plot_title = NULL,
                    y_range = c(75, 105),
                    jitter_width = 0.2,
                    rv_colors = get("rv_colors", envir = globalenv()),
                    return_parts = FALSE) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combo_levels <- c(
    "n_ind=1, high ILS",
    "n_ind=2, high ILS",
    "n_ind=1, low ILS",
    "n_ind=2, low ILS"
  )

  fmt_rate <- function(x) {
    v <- as.numeric(x)
    ifelse(v == 0, "0", sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
  }

  plot_data <- df %>%
    mutate(
      n_ind = if ("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L", "GL"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF), n_ind %in% c(1, 2)) %>%
    select(dup_loss_rate, RV, SF, n_ind,
           H0 = all_of(h0_col),
           H1 = all_of(h1_col),
           H1_noF = all_of(h1_noF_col)) %>%
    pivot_longer(cols = c(H0, H1, H1_noF),
                 names_to = "metric", values_to = "percentage") %>%
    mutate(
      metric = factor(metric,
                      levels = c("H0", "H1", "H1_noF"),
                      labels = c("h=0", "h=1", "h=1, taxon F pruned")),
      ratevar_name = recode(as.character(droplevels(RV)),
                            "N" = "none", "G" = "gene", "L" = "lineage"),
      ILS      = ifelse(SF == 1.0, "high", "low"),
      fill_var = ifelse(n_ind == 1, ratevar_name, "hollow"),
      combo    = case_when(
        n_ind == 1 & ILS == "high" ~ combo_levels[1],
        n_ind == 2 & ILS == "high" ~ combo_levels[2],
        n_ind == 1 & ILS == "low"  ~ combo_levels[3],
        n_ind == 2 & ILS == "low"  ~ combo_levels[4]
      ),
      dup_loss_rate_cat = factor(
        fmt_rate(dup_loss_rate),
        levels = fmt_rate(sort(unique(dup_loss_rate)))
      )
    ) %>%
    mutate(
      ratevar_name = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var     = factor(
        fill_var, levels = c("none", "gene", "lineage", "hollow")),
      combo        = factor(combo,         levels = combo_levels)
    )

  if (is.null(plot_title)) {
    plot_title <- plot_title 
  }

  fill_values <- c(rv_colors, "hollow" = "white")
  shape_values <- setNames(c(21, 21, 24, 24), combo_levels)
  legend_fill_override <- c("grey55", "white", "grey55", "white")

  base_theme <- theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 21, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 20, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 20, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 20, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 20, face = "bold"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(1.8, "lines"),
      legend.text        = element_text(size = 20, face = "bold"),
      legend.title       = element_text(size = 20, face = "bold"),
      panel.spacing      = unit(1.0, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    )

  base_scales <- list(
    scale_color_manual(values = rv_colors,
                       name = "rate variation:",
                       labels = c("none" = "none", "gene" = "across genes",
                                  "lineage" = "across lineages")),
    scale_fill_manual(values = fill_values, guide = "none"),
    scale_shape_manual(values = shape_values,
                       # trailing spaces in name are intentional for layout
                       name = "ILS & individuals / taxon (n):    ",
                       labels = c("n_ind=1, high ILS" = "ILS = high, n = 1",
                                  "n_ind=2, high ILS" = "ILS = high, n = 2",
                                  "n_ind=1, low ILS"  = "ILS = low, n = 1",
                                  "n_ind=2, low ILS"  = "ILS = low, n = 2")),
    scale_y_continuous(limits = y_range, expand = c(0, 0))
  )

  base_plot <- ggplot(plot_data,
                      aes(x = dup_loss_rate_cat, y = percentage,
                          color = ratevar_name, fill = fill_var,
                          shape = combo)) +
    geom_jitter(size   = 5.5, stroke = 1.2,
                position = position_jitter(
                  width = jitter_width, height = 0, seed = 42)) +
    facet_wrap(~ metric, nrow = 1) +
    base_scales +
    labs(x = "duplication/loss rate",
         y = "prob. of displaying the true tree (%)",
         title = plot_title) +
    base_theme

  # Extract right legend (color only)
  p_right_leg <- base_plot +
    theme(legend.position = "right",
    legend.box.spacing = unit(0.0, "lines"),
          legend.spacing.x   = unit(0.0, "cm"),
          legend.background  = element_rect(color = NA),) +
    guides(
      color = guide_legend(override.aes = list(size = 7, stroke = 1.6),
                           ncol = 1, reverse = TRUE),
      shape = "none"
    )

  p_bottom_leg <- base_plot +
  theme(
    legend.position  = "bottom",
    legend.box.spacing = unit(0.0, "lines"),
    legend.spacing.x   = unit(0.0, "cm"),
    legend.margin      = margin(10, 10, 10, 10),
    legend.key.width   = unit(2.5, "cm"),
    legend.title       = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 6)
    )
  ) +
  guides(
    color = "none",
    shape = guide_legend(
      override.aes = list(
        size = 7,
        stroke = 1.6,
        fill = legend_fill_override,
        color = "black"
      ),
      nrow = 2
    )
  )

  legend_right  <- cowplot::get_legend(p_right_leg)
  legend_bottom <- cowplot::get_legend(p_bottom_leg)

  # Wrap bottom legend as a ggdraw object so we can apply negative top margin
  # to eliminate the gap that plot_grid otherwise inserts between rows
  legend_bottom_plot <- cowplot::ggdraw(legend_bottom) +
    theme(plot.margin = unit(c(-1.5, -1.0, -1.0, -1.0), "cm"))

  # Main plot with no legend
  p_no_legend <- base_plot + theme(
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "lines"))

  # If caller wants the component parts, return them without assembling
  if (return_parts) {
    return(list(
      p_no_legend   = p_no_legend,
      base_plot     = base_plot,
      legend_right  = legend_right,
      legend_bottom = legend_bottom
    ))
  }

  # Assemble: main plot + right legend side by side, then bottom legend below
  top_row <- cowplot::plot_grid(p_no_legend, legend_right,
                                nrow = 1, rel_widths = c(1, 0.25))
  p <- cowplot::plot_grid(top_row, legend_bottom_plot,
                          ncol = 1, rel_heights = c(1, 0.30))

  output_path <- file.path(output_dir, output_filename)
  ggsave(output_path, plot = p, width = 17.3, height = 8, units = "in")

  cat("Saved aggregated true graph recovery jitter plot:", output_path, "\n")

  invisible(p)
}

#' Plot H=1 true graph recovery rates by n_inds
#' 
#' Plot H=1 true graph recovery rates (single jitter panel)
#'
#' Single-panel jitter plot for H=1 recovery only:
#' - Color encodes substitution rate variation (none/gene/lineage)
#' - Shape + fill encode ILS × n_ind (4 combos):
#'     n_ind=1, high ILS → filled circle   (shape 21)
#'     n_ind=1, low  ILS → filled triangle (shape 24)
#'     n_ind=2, high ILS → hollow circle   (shape 21)
#'     n_ind=2, low  ILS → hollow triangle (shape 24)
#' - No connecting lines; points jittered horizontally
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind
#' @param h1_col Column name for H=1 true graph counts
#'   (default: "find_true_net1")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param y_range Numeric vector for y-axis limits (default: c(75, 105))
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @param jitter_width Horizontal jitter width (default: 0.2)
#' @param rv_colors Named color vector for rate variation levels
#' @return ggplot object invisibly
plot_true_graph_recovery_H1_only <- function(df,
                                             h1_col = "find_true_net1",
                                             output_dir,
                                             output_filename,
                                             y_range = c(75, 105),
                                             plot_title = NULL,
                                             jitter_width = 0.2,
                                             rv_colors = get(
                                               "rv_colors",
                                               envir = globalenv())) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combo_levels <- c(
    "n_ind=1, high ILS",
    "n_ind=2, high ILS",
    "n_ind=1, low ILS",
    "n_ind=2, low ILS"
  )

  fmt_rate <- function(x) {
    v <- as.numeric(x)
    ifelse(v == 0, "0", sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
  }

  plot_data <- df %>%
    mutate(
      n_ind = if ("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF), n_ind %in% c(1, 2)) %>%
    select(dup_loss_rate, RV, SF, n_ind, H1 = all_of(h1_col)) %>%
    mutate(
      ratevar_name = recode(as.character(droplevels(RV)),
                            "N" = "none", "G" = "gene", "L" = "lineage"),
      ILS      = ifelse(SF == 1.0, "high", "low"),
      fill_var = ifelse(n_ind == 1, ratevar_name, "hollow"),
      combo    = case_when(
        n_ind == 1 & ILS == "high" ~ combo_levels[1],
        n_ind == 2 & ILS == "high" ~ combo_levels[2],
        n_ind == 1 & ILS == "low"  ~ combo_levels[3],
        n_ind == 2 & ILS == "low"  ~ combo_levels[4]
      ),
      dup_loss_rate_cat = factor(
        fmt_rate(dup_loss_rate),
        levels = fmt_rate(sort(unique(dup_loss_rate)))
      )
    ) %>%
    mutate(
      ratevar_name = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var     = factor(
        fill_var, levels = c("none", "gene", "lineage", "hollow")),
      combo        = factor(combo,         levels = combo_levels)
    )

  if (is.null(plot_title)) {
    plot_title <- "H=1 true graph recovery rate by duplication/loss rate"
  }

  fill_values <- c(rv_colors, "hollow" = "white")
  shape_values <- setNames(c(21, 21, 24, 24), combo_levels)
  legend_fill_override <- c("grey55", "white", "grey55", "white")

  p <- ggplot(plot_data,
              aes(x = dup_loss_rate_cat, y = H1,
                  color = ratevar_name, fill = fill_var, shape = combo)) +
    geom_jitter(size   = 5.5, stroke = 1.2,
                position = position_jitter(
                  width = jitter_width, height = 0, seed = 42)) +
    scale_color_manual(values = rv_colors,    name = "rate variation") +
    scale_fill_manual(values  = fill_values,  guide = "none") +
    scale_shape_manual(values = shape_values, name = "ILS \u00d7 n_ind") +
    scale_y_continuous(limits = y_range, expand = c(0, 0)) +
    labs(x = "duplication/loss rate",
         y = "number of replicates finding true graph under H=1",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.box         = "horizontal",
      legend.box.spacing = unit(0.4, "lines"),
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(1.8, "lines"),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 7, stroke = 1.6), nrow = 1),
      shape = guide_legend(override.aes = list(size = 7, stroke = 1.6,
                                               fill = legend_fill_override,
                                               color = "black"), nrow = 2)
    )

  output_path <- file.path(output_dir, output_filename)
  ggsave(output_path, plot = p, width = 11, height = 9, units = "in")

  cat("Saved H=1 only true graph recovery jitter plot:", output_path, "\n")

  invisible(p)
}

#' Plot combined SNaQ vs findgraphs H=1 true graph recovery (two-panel jitter)
#'
#' Creates a two-panel jitter plot (SNaQ | findgraphs) for H=1 recovery rates:
#' - Color encodes substitution rate variation (none/gene/lineage)
#' - Shape + fill encode ILS × n_ind (4 combos):
#'     n_ind=1, high ILS → filled circle   (shape 21)
#'     n_ind=1, low  ILS → filled triangle (shape 24)
#'     n_ind=2, high ILS → hollow circle   (shape 21)
#'     n_ind=2, low  ILS → hollow triangle (shape 24)
#' - No connecting lines; points jittered horizontally
#'
#' @param snaq_df      Data frame for SNaQ (must contain dup_loss_rate, RV,
#'   SF, n_ind/N_ind, and h1_col)
#' @param findgraph_df Data frame for findgraphs (same structure)
#' @param h1_col_snaq      Column name for H=1 counts in snaq_df
#'   (default: "find_true_net1")
#' @param h1_col_findgraph Column name for H=1 counts in findgraph_df
#'   (default: "find_true_net1")
#' @param output_dir   Directory to save output figure
#' @param output_filename Output filename (without path, no extension)
#' @param y_range      Numeric vector for y-axis limits (default: c(75, 105))
#' @param plot_title   Plot title (optional, auto-generated if NULL)
#' @param jitter_width Horizontal jitter width (default: 0.2)
#' @param rv_colors    Named color vector for rate variation levels
#' @return ggplot object invisibly
plot_combined_true_graph_recovery_H1_only <- function(snaq_df,
                                           findgraph_df,
                                           h1_col_snaq      = "find_true_net1",
                                           h1_col_findgraph = "find_true_net1",
                                           output_dir,
                                           output_filename,
                                           y_range = c(75, 105),
                                           plot_title = NULL,
                                           jitter_width = 0.2,
                                           rv_colors = get(
                                             "rv_colors",
                                             envir = globalenv())) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combo_levels <- c(
    "n_ind=1, high ILS",
    "n_ind=2, high ILS",
    "n_ind=1, low ILS",
    "n_ind=2, low ILS"
  )

  fmt_rate <- function(x) {
    v <- as.numeric(x)
    ifelse(v == 0, "0", sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
  }

  prep <- function(df, h1_col, source_label) {
    df %>%
      mutate(
        n_ind = if ("n_ind" %in% names(.)) n_ind else N_ind,
        SF = as.numeric(SF),
        dup_loss_rate = as.numeric(dup_loss_rate),
        RV = factor(RV, levels = c("N", "G", "L"))
      ) %>%
      filter(!is.na(RV), !is.na(n_ind), !is.na(SF), n_ind %in% c(1, 2)) %>%
      select(dup_loss_rate, RV, SF, n_ind, H1 = all_of(h1_col)) %>%
      mutate(source = source_label)
  }

  combined <- bind_rows(
    prep(snaq_df,      h1_col_snaq,      "SNaQ"),
    prep(findgraph_df, h1_col_findgraph, "find_graphs")
  ) %>%
    mutate(
      ratevar_name = recode(as.character(droplevels(RV)),
                            "N" = "none", "G" = "gene", "L" = "lineage"),
      ILS      = ifelse(SF == 1.0, "high", "low"),
      fill_var = ifelse(n_ind == 1, ratevar_name, "hollow"),
      combo    = case_when(
        n_ind == 1 & ILS == "high" ~ combo_levels[1],
        n_ind == 2 & ILS == "high" ~ combo_levels[2],
        n_ind == 1 & ILS == "low"  ~ combo_levels[3],
        n_ind == 2 & ILS == "low"  ~ combo_levels[4]
      ),
      dup_loss_rate_cat = factor(
        fmt_rate(dup_loss_rate),
        levels = fmt_rate(sort(unique(dup_loss_rate)))
      ),
      source = factor(source, levels = c("find_graphs", "SNaQ"))
    ) %>%
    mutate(
      ratevar_name = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var     = factor(
        fill_var, levels = c("none", "gene", "lineage", "hollow")),
      combo        = factor(combo,         levels = combo_levels)
    )

  fill_values <- c(rv_colors, "hollow" = "white")
  shape_values <- setNames(c(21, 21, 24, 24), combo_levels)
  legend_fill_override <- c("grey55", "white", "grey55", "white")

  p <- ggplot(combined,
              aes(x = dup_loss_rate_cat, y = H1,
                  color = ratevar_name, fill = fill_var, shape = combo)) +
    geom_jitter(size   = 5.5, stroke = 1.2,
                position = position_jitter(
                  width = jitter_width, height = 0, seed = 42)) +
    facet_wrap(~ source, nrow = 1) +
    scale_color_manual(values = rv_colors,
                       name   = "rate variation:",
                       breaks = c("lineage", "gene", "none"),
                       labels = c("lineage" = "across lineages",
                                  "gene"    = "across genes",
                                  "none"    = "none")) +
    scale_fill_manual(values  = fill_values,  guide = "none") +
    scale_shape_manual(values = shape_values,
                       name = "ILS & individuals / taxon (n)") +
    scale_y_continuous(limits = y_range, expand = c(0, 0)) +
    labs(x = "duplication and loss rate",
         y = "probability of displaying the true tree (%)",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 20, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 20, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 18, angle = 45, hjust = 1, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 18, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 16, face = "bold"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "right",
      legend.box         = "vertical",
      legend.box.spacing = unit(0.4, "lines"),
      legend.background  = element_blank(),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(1.8, "lines"),
      legend.text        = element_text(size = 18, face = "bold"),
      legend.title       = element_text(size = 18, face = "bold"),
      panel.spacing      = unit(1.5, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 7, stroke = 1.6),
                           nrow = 3, position = "right", order = 1),
      shape = guide_legend(
        override.aes = list(size = 7, stroke = 1.6,
                            fill = legend_fill_override,
                            color = "black"),
        nrow = 2, position = "bottom",
        theme = theme(legend.background = element_rect(
          color = "black", linewidth = 1.0, fill = "white")))
    )

  output_path <- file.path(output_dir, paste0(output_filename, ".pdf"))
  ggsave(output_path, plot = p, width = 12, height = 8, units = "in")

  cat("Saved combined SNaQ/findgraphs H=1 recovery jitter plot:",
      output_path, "\n")

  invisible(p)
}

#' Plot statistics by true tree display status (generic for SNaQ and findgraphs)
#' 
#' Creates a histogram with overlapping distributions showing statistics
#' separated by true tree display status:
#' - Green: True tree is displayed
#' - Red: True tree is NOT displayed
#' Works with any numeric column from the input data.
#'
#' @param input_dir Directory containing CSV summary files
#' @param column_to_plot Column name to plot (e.g., "best_graph_gamma2",
#'   "H1_best_graph_WR", "minor_gamma")
#' @param display_status_col Column name indicating display status
#'   (e.g., "H1_network_displays_true_tree",
#'   "H1_best_graph_displayed_true_tree")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Title for the plot
#' @param x_label X-axis label (optional, defaults to column_to_plot)
#' @param y_label Y-axis label (default: "Frequency")
plot_statistics_by_tree_display <- function(input_dir,
                                           column_to_plot,
                                           display_status_col,
                                           output_dir,
                                           output_filename,
                                           plot_title,
                                           x_label = NULL,
                                           y_label = "Frequency",
                                           alpha = 0.6,
                                           color_displayed_tree = get(
                                             "color_displayed_tree",
                                             envir = globalenv())) {

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Default x_label to column_to_plot if not provided
  if (is.null(x_label)) {
    x_label <- column_to_plot
  }
  
  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Initialize vectors to store data
  values_displays <- numeric()  # Displays true tree (green)
  values_not_displays <- numeric()  # Does NOT display true tree (red)
  
  # Process each CSV file
  for (csv_file in csv_files) {
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Check that required columns exist
      if (!column_to_plot %in% names(df)) {
        warning("Column '", column_to_plot,
                "' not found in: ", basename(csv_file))
        next
      }
      if (!display_status_col %in% names(df)) {
        warning("Column '", display_status_col,
                "' not found in: ", basename(csv_file))
        next
      }
      
      # Extract columns and remove rows with NA
      values_col <- df[[column_to_plot]]
      status_col <- df[[display_status_col]]
      
      # Handle different truth value representations (TRUE, "True", 1, etc.)
      valid_rows <- !is.na(values_col) & !is.na(status_col)
      values_col <- values_col[valid_rows]
      status_col <- status_col[valid_rows]
      
      # Separate data by display status
      # Handle different representations of TRUE
      displays <- values_col[status_col == TRUE | 
                            status_col == "True" | 
                            status_col == "TRUE" |
                            status_col == 1]
      not_displays <- values_col[status_col == FALSE | 
                                status_col == "False" | 
                                status_col == "FALSE" |
                                status_col == 0]
      
      # Add to vectors (remove NAs)
      values_displays <- c(values_displays, na.omit(as.numeric(displays)))
      values_not_displays <- c(values_not_displays,
        na.omit(as.numeric(not_displays)))
      
    }, error = function(e) {
      cat("Error reading", basename(csv_file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (length(values_displays) == 0 && length(values_not_displays) == 0) {
    cat("No valid data found\n")
    return(invisible(NULL))
  }
  
  # Prepare data for ggplot
  plot_data <- data.frame(
    value = c(values_displays, values_not_displays),
    tree_display = c(rep("True tree displayed", length(values_displays)),
                    rep("True tree NOT displayed", length(values_not_displays)))
  )
  
  plot_data$tree_display <- factor(plot_data$tree_display,
                                   levels = c("True tree displayed",
                                             "True tree NOT displayed"))
  
  # Define colors (green for displays, red for not displays)
  colors <- color_displayed_tree[
    c("True tree displayed", "True tree NOT displayed")]
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = value, fill = tree_display)) +
    geom_histogram(alpha = alpha, bins = 30, position = "identity") +
    scale_fill_manual(values = colors,
                     name = "Network status") +
    labs(x = x_label,
         y = y_label,
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(1.0, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, units = "in")
  
  cat("Saved statistics by tree display plot:", output_path, "\n")
  cat("True tree displayed: ", length(values_displays), " records\n")
  cat("True tree NOT displayed: ", length(values_not_displays), " records\n")
  
  invisible(p)
}

#' Plot SNaQ minor gamma distributions by true tree display status
#' 
#' Wrapper around plot_statistics_by_tree_display for SNaQ data.
#' Creates a histogram with overlapping distributions showing:
#' - Green: H=1 network displays the true tree
#' - Red: H=1 network does NOT display the true tree
#'
#' @param snaq_input_dir Directory containing SNaQ summary CSV files
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
plot_snaq_minor_gamma_by_tree_display <- function(snaq_input_dir,
                                                   output_dir,
                                                   output_filename,
                                                   alpha = 0.6,
                                                   color_displayed_tree = get(
                                                     "color_displayed_tree",
                                                     envir = globalenv())) {

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Find all CSV files
  csv_files <- list.files(
    snaq_input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", snaq_input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Initialize vectors to store data
  minor_gamma_displays <- numeric()  # H=1 displays true tree (green)
  minor_gamma_not_displays <- numeric()  # H=1 does NOT display true tree (red)
  
  # Process each CSV file
  for (csv_file in csv_files) {
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Check required columns exist for SNaQ
      if (!all(c("gamma_1", "gamma_2",
                 "RF_net1_1_true", "RF_net1_2_true") %in% names(df))) {
        next
      }
      
      # Calculate minor gamma
      df$minor_gamma <- pmin(df$gamma_1, df$gamma_2)
      df <- df[!is.na(df$minor_gamma), ]
      
      if (nrow(df) == 0) next
      
      # Categorize rows based on whether H=1 displays true tree
      # H=1 displays true tree if RF_net1_1_true==0 OR RF_net1_2_true==0
      displays_true <- (df$RF_net1_1_true == 0.0) | (df$RF_net1_2_true == 0.0)
      
      minor_gamma_displays <- c(
        minor_gamma_displays, df$minor_gamma[displays_true])
      minor_gamma_not_displays <- c(
        minor_gamma_not_displays, df$minor_gamma[!displays_true])
      
    }, error = function(e) {
      cat("Error reading", basename(csv_file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (length(minor_gamma_displays) == 0 &&
      length(minor_gamma_not_displays) == 0) {
    cat("No valid data found\n")
    return(invisible(NULL))
  }
  
  # Prepare data for ggplot
  plot_data <- data.frame(
    minor_gamma = c(minor_gamma_displays, minor_gamma_not_displays),
    tree_display = c(
      rep("H=1 displays true tree", length(minor_gamma_displays)),
      rep("H=1 does NOT display true tree",
          length(minor_gamma_not_displays)))
  )
  
  plot_data$tree_display <- factor(plot_data$tree_display,
                                   levels = c("H=1 displays true tree",
                                             "H=1 does NOT display true tree"))
  
  # Define colors (green for displays, red for not displays)
  colors <- color_displayed_tree[c(
    "H=1 displays true tree", "H=1 does NOT display true tree")]
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = minor_gamma, fill = tree_display)) +
    geom_histogram(binwidth = 0.01, alpha = alpha, position = "identity") +
    scale_fill_manual(values = colors,
                     name = "Network status") +
    scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    labs(x = "minor gamma value (estimated gene flow proportion)",
         y = "frequency",
         title = paste0("SNaQ minor gamma distribution:\n",
           "H=1 network displays vs. does not display true tree")) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(1.0, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p, width = 7, height = 7, dpi = 300, units = "in")
  
  cat("Saved SNaQ minor gamma distribution plot:", output_path, "\n")
  cat("H=1 displays true tree: ", length(minor_gamma_displays), " records\n")
  cat("H=1 does NOT display true tree: ",
      length(minor_gamma_not_displays), " records\n")
}
                                                  

#' Plot combined overlapping minor gamma distributions for SNaQ and findgraphs
#' 
#' Creates a 2-panel plot aggregating minor gamma distributions across all
#' parameter settings:
#' - Left panel: SNaQ data
#' - Right panel: findgraphs data
#' - Colors represent substitution rate variations (gene, lineage, none)
#' - Data is aggregated across all n_inds, SF, and dup/loss rate combinations
#'
#' @param snaq_input_dir Directory containing SNaQ summary CSV files
#' @param findgraph_input_dir Directory containing findgraph summary CSV files
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param max_y Maximum y-axis value for all plots (optional, auto-calculated
#'   if NULL)
plot_combined_snaq_findgraph_aggregated <- function(snaq_input_dir,
                                                     findgraph_input_dir,
                                                     output_dir,
                                                     output_filename,
                                                     max_y = NULL,
                                                     alpha = 0.6,
                                                     rv_colors = get(
                                               "rv_colors",
                                               envir = globalenv())) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize combined dataframe
  combined_data <- data.frame()
  
  # Process SNaQ files
  cat("Processing SNaQ files...\n")
  snaq_files <- list.files(
    snaq_input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in snaq_files) {
    filename <- basename(csv_file)
    rv_match <- str_match(filename, "-(RVG|RVL|RVN)-")
    
    if (is.na(rv_match[2])) {
      next
    }
    
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      if (!("gamma_1" %in% names(df) && "gamma_2" %in% names(df))) {
        next
      }
      
      # Calculate minor gamma
      df$minor_gamma <- pmin(df$gamma_1, df$gamma_2)
      df <- df[!is.na(df$minor_gamma), ]
      
      if (nrow(df) == 0) next
      
      # Add source and ratevar
      df$source <- "SNaQ"
      df$ratevar <- rv_match[2]
      
      combined_data <- rbind(combined_data, 
        df[, c("minor_gamma", "source", "ratevar")])
    }, error = function(e) {
      cat("Error reading SNaQ file", filename, ":", conditionMessage(e), "\n")
    })
  }
  
  # Process findgraph files
  cat("Processing findgraph files...\n")
  findgraph_files <- list.files(
    findgraph_input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in findgraph_files) {
    filename <- basename(csv_file)
    rv_match <- str_match(filename, "-(RVG|RVL|RVN)-")
    
    if (is.na(rv_match[2])) {
      next
    }
    
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      if (!("best_graph_gamma1" %in% names(df) &&
            "best_graph_gamma2" %in% names(df))) {
        next
      }
      
      # Calculate minor gamma
      df$minor_gamma <- pmin(df$best_graph_gamma1, df$best_graph_gamma2)
      df <- df[!is.na(df$minor_gamma), ]
      
      if (nrow(df) == 0) next
      
      # Add source and ratevar
      df$source <- "findgraph"
      df$ratevar <- rv_match[2]
      
      combined_data <- rbind(combined_data,
        df[, c("minor_gamma", "source", "ratevar")])
    }, error = function(e) {
      cat("Error reading findgraph file", filename, ":",
          conditionMessage(e), "\n")
    })
  }
  
  if (nrow(combined_data) == 0) {
    cat("No valid data found\n")
    return(invisible(NULL))
  }
  
  cat("Total rows in combined data:", nrow(combined_data), "\n")
  
  # Map codes to readable names
  combined_data$ratevar_name <- recode(combined_data$ratevar,
                                       "RVG" = "gene",
                                       "RVL" = "lineage",
                                       "RVN" = "none")
  
  combined_data$source_label <- recode(combined_data$source,
                                       "SNaQ" = "SNaQ",
                                       "findgraph" = "findgraph")
  
  # Factor levels
  combined_data$ratevar_name <- factor(combined_data$ratevar_name,
                                       levels = c("gene", "lineage", "none"))
  
  combined_data$source_label <- factor(combined_data$source_label,
                                       levels = c("SNaQ", "findgraph"))
  
  # Define colors (more distinctive: blue, orange, dark gray)
  colors <- rv_colors
  
  # Auto-calculate max_y if not provided
  if (is.null(max_y)) {
    # Use hist with same breaks as the actual plot (binwidth = 0.01)
    hist_data <- hist(combined_data$minor_gamma,
      breaks = seq(0, 0.5, 0.01), plot = FALSE)
    max_y <- ceiling(max(hist_data$counts)) * 1.15
    if (is.na(max_y) || max_y == 0) {
      max_y <- 100  # Fallback default
    }
  }
  
  cat("Using max_y =", max_y, "\n")
  
  # Create the plot with facet_wrap for source
  p <- ggplot(combined_data, aes(x = minor_gamma, fill = ratevar_name)) +
    geom_histogram(binwidth = 0.01, alpha = alpha, position = "identity") +
    facet_wrap(~source_label, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = colors,
                     name = NULL,
                     labels = c("substitution rate: gene",
                              "substitution rate: lineage",
                              "substitution rate: none")) +
    scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, max_y)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "estimated gene flow proportion",
         y = "number of replicates",
         title = paste0("estimated gene flow proportion distributions:",
           " SNaQ vs findgraph")) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 14, face = "bold"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      legend.key.size    = unit(1.2, "cm"),
      panel.spacing      = unit(1.5, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 16, height = 10, dpi = 300, units = "in")
  
  cat("Saved combined 2-panel figure:", output_path, "\n")
  
  # Return the combined data and plot for further use
  invisible(list(data = combined_data, plot = p))
}

#' Plot WR (worst residual) distributions by n_inds and SF
#' 
#' Creates 12-panel plots for H0 and H1 WR distributions with facet_grid where:
#' - Rows represent (n_inds, ILS) combinations
#' - Columns represent duplication/loss rates
#' - Colors represent substitution rate variations (gene, lineage, none)
#'
#' @param input_dir Directory containing CSV files
#' @param wr_column Column name to plot (e.g., "H0_best_tree_WR" or
#'   "H1_best_graph_WR")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param max_y Maximum y-axis value for all plots (optional, auto-calculated
#'   if NULL)
plot_WR_distributions <- function(input_dir, 
                                 wr_column,
                                 output_dir, 
                                 output_filename,
                                 plot_title = NULL, 
                                 max_y = NULL,
                                 alpha = 0.6) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Read and combine all data
  all_data <- data.frame()
  
  for (csv_file in csv_files) {
    # Extract parameters from filename
    filename <- basename(csv_file)
    
    n_inds_match <- str_match(filename, "N_ind(\\d+)")
    sf_match <- str_match(filename, "SF([\\d.]+)")
    rate_match <- str_match(filename, "DUP([\\d.e-]+)-LOS([\\d.e-]+)")
    rv_match <- str_match(filename, "-(RVG|RVL|RVN|RVGL)-")
    
    if (any(is.na(c(n_inds_match, sf_match, rate_match, rv_match)))) {
      next
    }
    
    # Read the CSV file
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      if (!(wr_column %in% names(df))) {
        cat("Warning: Column", wr_column, "not found in", filename, "\n")
        next
      }
      
      # Get WR values, removing NA and non-numeric values
      df[[wr_column]] <- suppressWarnings(as.numeric(df[[wr_column]]))
      df_clean <- df[!is.na(df[[wr_column]]), ]
      
      if (nrow(df_clean) == 0) next
      
      # Add grouping variables
      df_clean$n_inds <- n_inds_match[2]
      df_clean$sf <- sf_match[2]
      df_clean$rate <- rate_match[2]
      df_clean$ratevar <- rv_match[2]
      df_clean$wr_value <- df_clean[[wr_column]]
      
      all_data <- rbind(all_data, 
                df_clean[, c("wr_value", "n_inds", "sf", "rate", "ratevar")])
    }, error = function(e) {
      cat("Error reading", csv_file, ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(all_data) == 0) {
    cat("No valid data found for column", wr_column, "\n")
    return(invisible(NULL))
  }
  
  # Map codes to readable names
  all_data$ratevar_name <- recode(all_data$ratevar,
                                  "RVG" = "gene",
                                  "RVL" = "lineage", 
                                  "RVN" = "none")
  
  all_data$ils_level <- recode(all_data$sf,
                               "1.0" = "high",
                               "0.5" = "low")
  
  # Create faceting variables
  all_data$row_facet <- paste0("individuals / taxon = ", 
                        all_data$n_inds, "\nILS=", 
                        all_data$ils_level)
  
  # Format rate in scientific notation for column facets
  all_data$rate_sci <- sapply(all_data$rate, function(r) {
    rate_num <- as.numeric(r)
    if (rate_num == 0) {
      "0.0"
    } else {
      formatted <- format(rate_num, scientific = TRUE, digits = 1)
      formatted <- gsub("e\\+00", "", formatted)
      formatted <- gsub("e-0", "e-", formatted)
      formatted
    }
  })
  all_data$col_facet <- paste0("dup/loss rate = ", all_data$rate_sci)
  
  # Order facets properly
  all_data$row_facet <- factor(all_data$row_facet, 
                               levels = c("individuals / taxon = 1\nILS=high", 
                                        "individuals / taxon = 1\nILS=low",
                                        "individuals / taxon = 2\nILS=high",
                                        "individuals / taxon = 2\nILS=low"))
  
  all_data$col_facet <- factor(all_data$col_facet,
                              levels = c("dup/loss rate = 0.0",
                                       "dup/loss rate = 3e-4",
                                       "dup/loss rate = 4e-4"))
  
  all_data$ratevar_name <- factor(all_data$ratevar_name,
                                 levels = c("gene", "lineage", "none"))
  
  # Define colors (more distinctive: blue, orange, dark gray)
  colors <- rv_colors
  
  # Auto-calculate max_y if not provided
  if (is.null(max_y)) {
    max_y <- ceiling(max(table(cut(all_data$wr_value, breaks = 50)))) * 1.1
  }
  
  # Determine appropriate binwidth based on data range
  data_range <- max(all_data$wr_value, na.rm = TRUE) -
    min(all_data$wr_value, na.rm = TRUE)
  binwidth <- data_range / 50  # Aim for ~50 bins

  # set up the graph title based on wr_column
  if (wr_column == "H0_best_tree_WR") {
    wr_column_name <- "WR of best tree under H=0"
  } else if (wr_column == "H1_best_graph_WR") {
    wr_column_name <- "WR of best tree under H=1"
  } else {
    wr_column_name <- "WR of true tree"
  }
  
  # Create the plot with facet_grid
  p <- ggplot(all_data, aes(x = wr_value, fill = ratevar_name)) +
    geom_histogram(binwidth = binwidth, alpha = alpha, position = "identity") +
    geom_vline(aes(xintercept = 3.0, color = "x = 3.0"),
               linetype = "dashed", linewidth = 1.2) +
    geom_vline(aes(xintercept = 3.7, color = "x = 3.7"),
               linetype = "dashed", linewidth = 1.2) +
    facet_grid(row_facet ~ col_facet, scales = "free_x") +
    scale_fill_manual(values = colors,
                     name = "rate variation:",
                     labels = c("across genes", "across lineage", "none")) +
    scale_color_manual(values = c("x = 3.0" = "red", "x = 3.7" = "green4"),
                       name = "WR threshold:") +
    scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
    labs(x = "worst residual on the true tree topology",
         y = "number of replicates",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 21, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 21, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 20, color = "black", face = "bold", angle = 45, hjust = 1),
      axis.text.y       = element_text(
        size = 20, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 18, face = "bold"),
      strip.text.x      = element_text(color = "black"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "bottom",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 18, face = "bold"),
      legend.title       = element_text(size = 18, face = "bold"),
      legend.key.size    = unit(1.5, "cm"),
      panel.spacing.x    = unit(1.0, "lines"),
      panel.spacing.y    = unit(1.0, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 14, height = 16, dpi = 300, units = "in")
  
  cat("Saved WR distribution plot:", output_path, "\n")
  
  invisible(p)
}

#' Plot combined minor gamma distributions by tree display status
#'
#' Creates a 2-panel facet_grid plot with SNaQ and findgraphs side-by-side:
#' - Columns represent data source (SNaQ, findgraphs)
#' - Colors represent tree display status (green = displays, red = doesn't)
#' - Overlapping histograms show distribution of minor_gamma values
#' - Y-axis shows frequency
#'
#' @param snaq_input_dir Directory containing SNaQ summary CSV files
#' @param findgraph_input_dir Directory containing findgraph summary CSV files
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param max_y Maximum y-axis value for all plots (optional, auto-calculated
#'   if NULL)
plot_combined_snaq_findgraph_by_tree_display <- function(
    snaq_input_dir,
    findgraph_input_dir,
    output_dir,
    output_filename,
    max_y = NULL,
    max_y_by_ratevar = NULL,
    alpha = 0.6,
    plot_title = NULL,
    color_displayed_tree = get(
      "color_displayed_tree", envir = globalenv())) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize combined dataframe
  combined_data <- data.frame()
  
  # Process SNaQ files
  cat("Processing SNaQ files for tree display analysis...\n")
  snaq_files <- list.files(
    snaq_input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in snaq_files) {
    tryCatch({
      df <- read_csv(csv_file, show_col_types = FALSE)
      
      # Extract rate variation from filename
      ratevar_code  <- str_extract(basename(csv_file), "RV[GLN]")
      ratevar_label <- case_when(
        ratevar_code == "RVG" ~ "gene",
        ratevar_code == "RVL" ~ "lineage",
        ratevar_code == "RVN" ~ "none",
        TRUE ~ NA_character_
      )
      
      # SNaQ has gamma_1 and gamma_2 columns; minor gamma is the minimum
      # Tree display: true if RF_net1_1_true==0 OR RF_net1_2_true==0
      if ("gamma_1" %in% names(df) && "gamma_2" %in% names(df) &&
          "RF_net1_1_true" %in% names(df) && "RF_net1_2_true" %in% names(df)) {
        
        temp_data <- df %>%
          mutate(minor_gamma = pmin(gamma_1, gamma_2),
                 tree_displays = (RF_net1_1_true == 0 | RF_net1_2_true == 0),
                 source = "SNaQ",
                 ratevar = ratevar_label,
                 tree_display_text = ifelse(tree_displays, 
                                           "True tree displayed",
                                           "True tree NOT displayed")) %>%
          select(minor_gamma, tree_display_text, source, ratevar) %>%
          filter(!is.na(minor_gamma),
                 !is.na(tree_display_text), !is.na(ratevar))

        combined_data <- rbind(combined_data, temp_data)
        cat("  Found", nrow(temp_data), "SNaQ replicates\n")
      } else {
        cat("  Missing required columns in", basename(csv_file), "\n")
      }
    }, error = function(e) {
      cat("Error reading SNaQ file", basename(csv_file), ":",
          conditionMessage(e), "\n")
    })
  }
  
  # Process findgraph files
  cat("Processing findgraph files for tree display analysis...\n")
  findgraph_files <- list.files(
    findgraph_input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in findgraph_files) {
    tryCatch({
      df <- read_csv(csv_file, show_col_types = FALSE)
      
      # Extract rate variation from filename
      ratevar_code  <- str_extract(basename(csv_file), "RV[GLN]")
      ratevar_label <- case_when(
        ratevar_code == "RVG" ~ "gene",
        ratevar_code == "RVL" ~ "lineage",
        ratevar_code == "RVN" ~ "none",
        TRUE ~ NA_character_
      )
      
      # findgraph: minor gamma = best_graph_gamma2;
      # display col = H1_best_graph_displayed_true_tree
      if ("best_graph_gamma2" %in% names(df) &&
          "H1_best_graph_displayed_true_tree" %in% names(df)) {

        temp_data <- df %>%
          mutate(source = "find_graphs",
                 ratevar = ratevar_label,
                 tree_display_text = ifelse(
                   H1_best_graph_displayed_true_tree == 1 |
                   H1_best_graph_displayed_true_tree == "True",
                   "True tree displayed",
                   "True tree NOT displayed")) %>%
          select(minor_gamma = "best_graph_gamma2",
                 tree_display_text, source, ratevar) %>%
          filter(!is.na(minor_gamma),
                 !is.na(tree_display_text), !is.na(ratevar))
        
        combined_data <- rbind(combined_data, temp_data)
        cat("  Found", nrow(temp_data), "findgraph replicates\n")
      } else {
        cat("  Missing required columns in", basename(csv_file), "\n")
      }
    }, error = function(e) {
      cat("Error reading findgraph file", basename(csv_file), ":",
          conditionMessage(e), "\n")
    })
  }

  if (nrow(combined_data) == 0) {
    cat("No valid data found for tree display analysis\n")
    return(invisible(NULL))
  }
  
  cat("Total records in combined data:", nrow(combined_data), "\n")
  
  # Factor levels
  combined_data$source <- factor(combined_data$source,
                                 levels = c("find_graphs", "SNaQ"))
  
  combined_data$ratevar <- factor(combined_data$ratevar,
                                  levels = c("none", "gene", "lineage"),
                                  labels = c("no rate variation", 
                                          "rate variation across genes", 
                                          "rate variation across lineages"))
  
  combined_data$tree_display_text <- factor(combined_data$tree_display_text,
                                           levels = c("True tree displayed", 
                                                     "True tree NOT displayed"))
  
  # Define colors
  colors <- color_displayed_tree[
    c("True tree displayed", "True tree NOT displayed")]
  
  # Auto-calculate global max_y (fallback for per-ratevar limits)
  if (is.null(max_y)) {
    max_y <- max(
      ggplot_build(ggplot(combined_data, aes(x = minor_gamma)) +
        geom_histogram(binwidth = 0.005))$data[[1]]$count,
      na.rm = TRUE) * 1.1
  }

  # Per-ratevar y-limit lookup (keys: "none", "gene", "lineage")
  ratevar_levels <- c("no rate variation",
    "rate variation across genes", "rate variation across lineages")
  ratevar_keys   <- c("none", "gene", "lineage")
  ylim_vals <- setNames(
    sapply(ratevar_keys, function(k) {
      if (!is.null(max_y_by_ratevar) &&
          !is.null(max_y_by_ratevar[[k]])) max_y_by_ratevar[[k]] else max_y
    }),
    ratevar_levels
  )
  cat("Using y limits by rate variation:\n")
  print(ylim_vals)
  
  # Dummy data to enforce per-row y limits with scales = "free_y"
  dummy_limits <- data.frame(
    minor_gamma       = 0,
    ratevar           = factor(ratevar_levels, levels = ratevar_levels),
    source            = factor("find_graphs",
      levels = c("find_graphs", "SNaQ")),
    ylim_max          = unname(ylim_vals),
    tree_display_text = factor("True tree displayed",
      levels = c("True tree displayed", "True tree NOT displayed"))
  )
  
  # Create the plot: rows = rate variation, cols = source (SNaQ / findgraphs)
  p <- ggplot(combined_data, aes(x = minor_gamma, fill = tree_display_text)) +
    geom_blank(data = dummy_limits, aes(y = ylim_max)) +
    geom_histogram(binwidth = 0.005, alpha = alpha, position = "identity") +
    facet_grid(ratevar ~ source, scales = "free_y") +
    scale_fill_manual(values = colors,
                     name = "is the true tree displayed in estimated networks?",
                     labels = c("yes",
                              "no")) +
    scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "estimated gene flow proportion",
         y = "number of replicates",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 21, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 21, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 20, color = "black", face = "bold", angle = 45, hjust = 1),
      axis.text.y       = element_text(
        size = 20, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 18, face = "bold"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "bottom",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 19.5, face = "bold"),
      legend.title       = element_text(size = 19.5, face = "bold"),
      legend.key.size    = unit(1.2, "cm"),
      panel.spacing      = unit(1.5, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 10.3, height = 14.1, dpi = 300, units = "in")
  
  cat("Saved combined SNaQ/findgraphs tree display plot:", output_path, "\n")
  
  # Print summary stats
  summary_data <- combined_data %>%
    group_by(source, tree_display_text) %>%
    summarise(count = n(), .groups = "drop")
  
  cat("\nSummary of records by source and tree display status:\n")
  print(summary_data)
  
  invisible(p)
}

#' Plot hypothesis acceptance percentage for H>=1
#' 
#' Creates two side-by-side line plots showing H>=1 (combined H=1 and H>1)
#' acceptance percentage by duplication/loss rate:
#' - Left panel: n_ind = 1
#' - Right panel: n_ind = 2
#' - Different line types represent ILS levels (solid=high ILS, dashed=low ILS)
#' - Different colors represent RV types (none, gene, lineage)
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind,
#'           and columns for H=0, H=1, and H>1 acceptance percentages
#' @param h0_col Column name for H=0 acceptance percentage (default: "pct_H0")
#' @param h1_col Column name for H=1 acceptance percentage (default: "pct_H1")
#' @param h_greater_col Column for H>1 acceptance percentage
#'   (default: "pct_H_greater")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @return ggplot object
plot_hypothesis_acceptance_by_rate_variation <- function(df,
                                                h0_col = "pct_H0",
                                                h1_col = "pct_H1",
                                                h_greater_col = "pct_H_greater",
                                                output_dir,
                                                output_filename,
                                                plot_title = NULL,
                                                rv_colors = get(
                                                  "rv_colors",
                                                  envir = globalenv())) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Standardize column names
  df_work <- df %>%
    mutate(
      n_ind = if("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF))
  
  # Handle alternative column names for H0, H1, and H_greater
  if (!h0_col %in% names(df_work)) {
    if ("H0Accepted" %in% names(df_work)) h0_col <- "H0Accepted"
    else if ("H=0Accepted" %in% names(df_work)) h0_col <- "H=0Accepted"
  }
  if (!h1_col %in% names(df_work)) {
    if ("H1Accepted" %in% names(df_work)) h1_col <- "H1Accepted"
    else if ("H=1Accepted" %in% names(df_work)) h1_col <- "H=1Accepted"
  }
  if (!h_greater_col %in% names(df_work)) {
    if ("H>1Accepted" %in% names(df_work)) h_greater_col <- "H>1Accepted"
  }
  
  # Prepare plot data - calculate H>=1 as H=1 + H>1
  plot_data <- df_work %>%
    select(dup_loss_rate, RV, SF, n_ind,
           H0 = all_of(h0_col),
           H1 = all_of(h1_col),
           H_greater = all_of(h_greater_col)) %>%
    mutate(
      # Calculate H>=1 as sum of H=1 and H>1
      H_greater_equal_1 = H1 + H_greater,
      dup_loss_rate_cat = factor(dup_loss_rate),
      RV = droplevels(RV),
      RV_label = recode(RV,
                       "N" = "none",
                       "G" = "gene",
                       "L" = "lineage"),
      ILS = ifelse(SF == 1.0, "high", "low"),
      n_ind_label = paste0("individuals / taxon = ", n_ind)
    ) %>%
    filter(n_ind %in% c(1, 2))
  
  # Order factors properly
  plot_data <- plot_data %>%
    mutate(
      RV_label = factor(RV_label, levels = c("none", "gene", "lineage")),
      n_ind_label = factor(n_ind_label,
        levels = c("individuals / taxon = 1",
                   "individuals / taxon = 2")),
      ILS = factor(ILS, levels = c("high", "low"))
    )
  
  # Create plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste0("H>=1 (Network) Hypothesis Acceptance",
                         " by Duplication/Loss Rate")
  }
  
  # Define colors for RV and line types for ILS
  ils_linetypes <- c("high" = "solid",
                     "low" = "dashed")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = dup_loss_rate_cat, y = H_greater_equal_1,
                            color = RV_label, linetype = ILS,
                            group = interaction(RV_label, ILS))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_grid(. ~ n_ind_label) +
    scale_color_manual("Substitution\nRate Variation", values = rv_colors) +
    scale_linetype_manual("ILS Level", values = ils_linetypes) +
    scale_y_continuous(limits = c(-1, 105), expand = c(0, 0)) +
    labs(x = "duplication/loss rate",
         y = "percentage of replicates",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 14, face = "bold"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey60", linewidth = 0.7, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey82", linewidth = 0.4, linetype = "dotted"),
      legend.position    = "right",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(0.8, "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, output_filename)
  ggsave(output_path, plot = p, width = 12, height = 6, units = "in")
  
  cat("Saved hypothesis acceptance (H>=1) by rate variation plot:",
      output_path, "\n")
  
  invisible(p)
}

#' Plot Type I error rate (H>=1 acceptance) as a single jitter panel
#'
#' For each parameter combination, computes \code{100 - pct_H0} as the Type I
#' error rate (percentage of replicates selecting H=1 or H>1).  All n_inds and
#' ILS combinations are shown in a single panel using fillable point shapes:
#' \itemize{
#'   \item n_inds=1, high ILS → filled circle  (shape 21, interior = RV color)
#'   \item n_inds=1, low  ILS → filled triangle (shape 24, interior = RV color)
#'   \item n_inds=2, high ILS → hollow circle  (shape 21, interior = white)
#'   \item n_inds=2, low  ILS → hollow triangle (shape 24, interior = white)
#' }
#' Point color (border) encodes substitution rate variation.  Points are
#' jittered horizontally.  No connecting lines or reference lines are drawn.
#'
#' @param df Data frame with columns: dup_loss_rate, RV, SF, n_ind or N_ind,
#'           and a column for H=0 acceptance percentage
#' @param h0_col Column name for H=0 acceptance percentage (default: "pct_H0")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path, no extension)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @param y_range Numeric vector for y-axis limits (default: c(-5, 105))
#' @param jitter_width Horizontal jitter width (default: 0.2)
#' @param rv_colors Named color vector for rate variation levels
#' @return ggplot object invisibly
plot_hypothesis_acceptance_H0_by_rate_variation <- function(
    df,
    h0_col = "pct_H0",
    output_dir,
    output_filename,
    plot_title = NULL,
    y_range = c(-1, 105),
    jitter_width = 0.2,
    rv_colors = get("rv_colors", envir = globalenv())) {

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Standardize column names
  df_work <- df %>%
    mutate(
      n_ind = if ("n_ind" %in% names(.)) n_ind else N_ind,
      SF = as.numeric(SF),
      dup_loss_rate = as.numeric(dup_loss_rate),
      RV = factor(RV, levels = c("N", "G", "L"))
    ) %>%
    filter(!is.na(RV), !is.na(n_ind), !is.na(SF))

  # Handle alternative column names for H0
  if (!h0_col %in% names(df_work)) {
    if ("H0Accepted" %in% names(df_work))    h0_col <- "H0Accepted"
    else if ("H=0Accepted" %in% names(df_work)) h0_col <- "H=0Accepted"
  }

  # 4-level combo labels (matches plot_WR_percentiles_jitter_combined)
  combo_levels <- c(
    "n_inds=1, high ILS",
    "n_inds=2, high ILS",
    "n_inds=1, low ILS",
    "n_inds=2, low ILS"
  )

  # Prepare plot data
  plot_data <- df_work %>%
    filter(n_ind %in% c(1, 2)) %>%
    select(dup_loss_rate, RV, SF, n_ind,
           H0_pct = all_of(h0_col)) %>%
    mutate(
      # Type I error = replicates NOT choosing H=0
      type1_error   = 100 - H0_pct,
      dup_loss_rate_cat = factor(dup_loss_rate),
      RV            = droplevels(RV),
      ratevar_name  = recode(as.character(RV),
                             "N" = "none",
                             "G" = "gene",
                             "L" = "lineage"),
      ILS           = ifelse(SF == 1.0, "high", "low"),
      # fill: RV color for n_ind=1, white for n_ind=2
      fill_var      = ifelse(n_ind == 1, ratevar_name, "hollow"),
      combo         = case_when(
        n_ind == 1 & ILS == "high" ~ combo_levels[1],
        n_ind == 2 & ILS == "high" ~ combo_levels[2],
        n_ind == 1 & ILS == "low"  ~ combo_levels[3],
        n_ind == 2 & ILS == "low"  ~ combo_levels[4]
      )
    ) %>%
    mutate(
      ratevar_name      = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var          = factor(fill_var,
        levels = c("none", "gene", "lineage", "hollow")),
      combo             = factor(combo,         levels = combo_levels),
      dup_loss_rate_cat = factor(dup_loss_rate_cat)
    )

  # Create plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste0(
      "Type I Error Rate (H\u22651 Accepted) by Duplication/Loss Rate")
  }

  # Scale definitions (mirroring plot_WR_percentiles_jitter_combined)
  fill_values <- c(rv_colors, "hollow" = "white")

  shape_values <- setNames(
    c(21, 21, 24, 24),
    combo_levels
  )

  legend_fill_override <- c("grey55", "white", "grey55", "white")

  # Build plot
  p <- ggplot(plot_data,
              aes(x     = dup_loss_rate_cat,
                  y     = type1_error,
                  color = ratevar_name,
                  fill  = fill_var,
                  shape = combo)) +
    geom_jitter(size     = 6.5,
                stroke   = 1.2,
                position = position_jitter(width = jitter_width, height = 0,
                                           seed  = 42)) +
    scale_color_manual(values = rv_colors,
                       name   = "Rate Variation") +
    scale_fill_manual(values  = fill_values,
                      guide   = "none") +
    scale_shape_manual(values = shape_values,
                       name   = "ILS \u00d7 n_inds") +
    scale_y_continuous(limits = y_range,
                       expand = c(0, 0)) +
    labs(x     = "duplication/loss rate",
         y     = "Type I error rate (%)",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.box         = "horizontal",
      legend.box.spacing = unit(0.4, "lines"),
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(2.0, "lines"),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 8.0, stroke = 1.6),
                           nrow = 1),
      shape = guide_legend(override.aes = list(size   = 8.0,
                                               stroke = 1.6,
                                               fill   = legend_fill_override,
                                               color  = "black"),
                           nrow = 2)
    )

  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 13, height = 10.5, dpi = 300, units = "in")

  cat("Saved Type I error rate jitter plot:", output_path, "\n")

  invisible(p)
}

#' Plot combined H0 acceptance (type I error) for SNaQ and findgraphs
#'
#' Creates a three-panel faceted jitter plot showing type I error rate
#' (H>=1 accepted):
#'   Panel 1: SNaQ
#'   Panel 2: findgraphs under WR threshold = 3.0
#'   Panel 3: findgraphs under WR threshold = 3.7
#'
#' @param snaq_df              Data frame of SNaQ summary results
#' @param findgraph_df         Data frame of findgraphs summary results
#' @param h0_col_snaq          Column for H0 rate in snaq_df (default "pct_H0")
#' @param h0_col_findgraph_wr30 Column for H0 rate (WR<=3.0) in findgraph_df
#'   (default "pct_H0")
#' @param h0_col_findgraph_wr37 Column for H0 rate (WR<=3.7) in findgraph_df
#'   (default "pct_H0_wr_3.7")
#' @param output_dir           Directory to save output figure
#' @param output_filename      Output filename (without extension)
#' @param plot_title           Optional plot title
#' @param y_range              Y-axis limits (default c(-1, 105))
#' @param jitter_width         Width for jitter (default 0.2)
#' @param rv_colors            Named vector of colors for rate variation levels
#' @return ggplot object invisibly
plot_combined_hypothesis_acceptance_H0 <- function(
    snaq_df,
    findgraph_df,
    h0_col_snaq           = "pct_H0",
    h0_col_findgraph_wr30 = "pct_H0",
    h0_col_findgraph_wr37 = "pct_H0_wr_3.7",
    output_dir,
    output_filename,
    plot_title    = NULL,
    y_range       = c(-1, 105),
    jitter_width  = 0.2,
    rv_colors = get("rv_colors", envir = globalenv())) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combo_levels <- c(
    "ILS = high, n = 1",
    "ILS = high, n = 2",
    "ILS = low, n = 1",
    "ILS = low, n = 2"
  )

  fmt_rate <- function(x) {
    v <- as.numeric(x)
    ifelse(v == 0, "0", sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
  }

  # Standardise one data frame and resolve a H0 column name
  standardise_df <- function(df) {
    df %>%
      mutate(
        n_ind = if ("n_ind" %in% names(.)) n_ind else N_ind,
        SF    = as.numeric(SF),
        dup_loss_rate = as.numeric(dup_loss_rate),
        RV    = factor(RV, levels = c("N", "G", "L"))
      ) %>%
      filter(!is.na(RV), !is.na(n_ind), !is.na(SF))
  }

  resolve_col <- function(df, h0_col) {
    if (!h0_col %in% names(df)) {
      if ("H0Accepted"    %in% names(df)) return("H0Accepted")
      if ("H=0Accepted"  %in% names(df)) return("H=0Accepted")
    }
    h0_col
  }

  snaq_std      <- standardise_df(snaq_df)
  findgraph_std <- standardise_df(findgraph_df)

  prep <- function(df, h0_col, source_label) {
    col <- resolve_col(df, h0_col)
    df %>%
      filter(n_ind %in% c(1, 2)) %>%
      select(dup_loss_rate, RV, SF, n_ind, H0_pct = all_of(col)) %>%
      mutate(
        type1_error       = 100 - H0_pct,
        dup_loss_rate_cat = factor(
          fmt_rate(dup_loss_rate),
          levels = fmt_rate(sort(unique(dup_loss_rate)))),
        RV                = droplevels(RV),
        ratevar_name      = recode(as.character(RV),
                                   "N" = "none", "G" = "gene", "L" = "lineage"),
        ILS               = ifelse(SF == 1.0, "high", "low"),
        fill_var          = ifelse(n_ind == 1, ratevar_name, "hollow"),
        combo             = case_when(
          n_ind == 1 & ILS == "high" ~ combo_levels[1],
          n_ind == 2 & ILS == "high" ~ combo_levels[2],
          n_ind == 1 & ILS == "low"  ~ combo_levels[3],
          n_ind == 2 & ILS == "low"  ~ combo_levels[4]
        ),
        source = source_label
      )
  }

  source_levels <- c("find_graphs (WR ≤ 3.0)",
                     "find_graphs (WR ≤ 3.7)",
                     "SNaQ")

  combined <- bind_rows(
    prep(findgraph_std, h0_col_findgraph_wr30,  source_levels[1]),
    prep(findgraph_std, h0_col_findgraph_wr37,  source_levels[2]),
    prep(snaq_std,      h0_col_snaq,            source_levels[3])
  ) %>%
    mutate(
      ratevar_name = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var     = factor(
        fill_var, levels = c("none", "gene", "lineage", "hollow")),
      combo        = factor(combo,        levels = combo_levels),
      source       = factor(source,       levels = source_levels)
    )

  fill_values          <- c(rv_colors, "hollow" = "white")
  shape_values         <- setNames(c(21, 21, 24, 24), combo_levels)
  legend_fill_override <- c("grey55", "white", "grey55", "white")

  p <- ggplot(combined,
              aes(x     = dup_loss_rate_cat,
                  y     = type1_error,
                  color = ratevar_name,
                  fill  = fill_var,
                  shape = combo)) +
    geom_jitter(size     = 5.5,
                stroke   = 1.2,
                position = position_jitter(
                  width = jitter_width, height = 0, seed = 42)) +
    facet_wrap(~ source, nrow = 1) +
    scale_color_manual(values = rv_colors,
                       name   = "rate variation:",
                       breaks = c("lineage", "gene", "none"),
                       labels = c("lineage" = "across lineages",
                                  "gene"    = "across genes",
                                  "none"    = "none")) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_shape_manual(values = shape_values,
                       name   = "ILS & individuals / taxon (n)") +
    scale_y_continuous(limits = y_range, expand = c(0, 0)) +
    labs(x     = "duplication and loss rate",
         y     = "type I error rate (%)",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 20, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 20, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 19, angle = 45, hjust = 1, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 19, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 19, face = "bold"),
      strip.background  = element_rect(fill = "grey90"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "right",
      legend.box         = "vertical",
      legend.box.spacing = unit(0.4, "lines"),
      legend.background  = element_blank(),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(1.8, "lines"),
      legend.text        = element_text(size = 19, face = "bold"),
      legend.title       = element_text(size = 19, face = "bold"),
      panel.spacing      = unit(1.5, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 7, stroke = 1.6),
                           nrow = 3, position = "right", order = 1),
      shape = guide_legend(override.aes = list(size   = 7,
                                               stroke = 1.6,
                                               fill   = legend_fill_override,
                                               color  = "black"),
                           nrow = 2, position = "bottom",
                           theme = theme(legend.background = element_rect(
                             color = "black", linewidth = 1.0,
                             fill = "white")))
    )

  output_path <- file.path(output_dir, paste0(output_filename, ".pdf"))
  ggsave(output_path, plot = p,
    width = 16, height = 8, dpi = 300, units = "in")

  cat("Saved combined SNaQ/findgraphs type I error jitter plot:",
    output_path, "\n")

  invisible(p)
}

#' Plot WR (worst residual) distribution by rate variation
#'
#' Creates an overlapping histogram showing WR distributions by rate variation.
#' Aggregates data across all parameter settings (n_inds, SF, dup/loss rate).
#' Similar to plot_combined_snaq_findgraph_aggregated but for WR data.
#'
#' @param input_dir Directory containing findgraph summary CSV files
#' @param wr_column Column name to plot (e.g., "H0_best_tree_WR",
#'   "H1_best_graph_WR", "true_tree_wr")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param max_y Maximum y-axis value (optional, auto-calculated if NULL)
#' @param binwidth Bin width for histogram (default: 0.01)
#' @return ggplot object invisibly
plot_WR_by_rate_variation <- function(
    input_dir,
    wr_column,
    output_dir,
    output_filename,
    max_y = NULL,
    rv_colors = get("rv_colors", envir = globalenv())) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Initialize data frame
  all_data <- data.frame()
  
  # Read and combine all data
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    rv_match <- str_match(filename, "-(RVG|RVL|RVN)-")
    
    if (is.na(rv_match[2])) {
      next
    }
    
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Check if WR column exists
      if (!(wr_column %in% names(df))) {
        cat("Warning: Column", wr_column, "not found in", filename, "\n")
        next
      }
      
      # Extract WR values and remove NA
      wr_values <- df[[wr_column]]
      wr_values <- wr_values[!is.na(wr_values)]
      
      if (length(wr_values) == 0) next
      
      # Add rate variation
      df_subset <- data.frame(
        wr_value = wr_values,
        ratevar = rv_match[2]
      )
      
      all_data <- rbind(all_data, df_subset)
    }, error = function(e) {
      cat("Error reading file", filename, ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(all_data) == 0) {
    cat("No valid WR data found\n")
    return(invisible(NULL))
  }
  
  cat("Total WR values in combined data:", nrow(all_data), "\n")
  
  # Map codes to readable names
  all_data$ratevar_name <- recode(all_data$ratevar,
                                   "RVG" = "gene",
                                   "RVL" = "lineage",
                                  "RVN" = "none")
  
  # Define colors
  colors <- rv_colors
  
  # Determine binwidth based on data range (matching plot_WR_distributions)
  data_range <- max(all_data$wr_value, na.rm = TRUE) -
    min(all_data$wr_value, na.rm = TRUE)
  binwidth <- data_range / 50  # Aim for ~50 bins
  
  # Auto-calculate max_y if not provided
  if (is.null(max_y)) {
    # Calculate breaks based on actual data range
    data_min <- min(all_data$wr_value, na.rm = TRUE)
    data_max <- max(all_data$wr_value, na.rm = TRUE)
    # Create breaks that span the actual data range
    breaks <- seq(floor(data_min / binwidth) * binwidth, 
                  ceiling(data_max / binwidth) * binwidth, 
                  by = binwidth)
    
    # Use hist with data-driven breaks
    hist_data <- hist(all_data$wr_value, breaks = breaks, plot = FALSE)
    max_y <- ceiling(max(hist_data$counts)) * 1.15
    if (is.na(max_y) || max_y == 0) {
      max_y <- 100  # Fallback default
    }
  }
  
  cat("Using max_y =", max_y, "\n")
  
  # Determine plot title based on wr_column
  if (wr_column == "H0_best_tree_WR") {
    plot_title <- paste0(
      "H0 best tree WR distribution by rate variation")
    x_label <- "worst residual (WR)"
  } else if (wr_column == "H1_best_graph_WR") {
    plot_title <- paste0(
      "H1 best graph WR distribution by rate variation")
    x_label <- "worst residual (WR)"
  } else if (wr_column == "true_tree_wr") {
    plot_title <- paste0(
      "true tree WR distribution by rate variation")
    x_label <- "worst residual (WR)"
  } else {
    plot_title <- paste(wr_column, "Distribution by Rate Variation")
    x_label <- wr_column
  }
  
  # Create the plot
  p <- ggplot(all_data, aes(x = wr_value, fill = ratevar_name)) +
    geom_histogram(binwidth = binwidth, alpha = 0.5, position = "identity") +
    geom_vline(xintercept = 3.0, color = "red",
      linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = 3.7, color = "green4",
      linetype = "dashed", linewidth = 1.2) +
    scale_fill_manual(values = colors,
                     name = NULL,
                     labels = c("substitution rate: gene",
                              "substitution rate: lineage",
                              "substitution rate: none")) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, max_y)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x_label,
         y = "frequency",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      legend.position    = "top",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      legend.key.size    = unit(1.2, "cm"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 16, height = 10, dpi = 300, units = "in")
  
  cat("Saved WR distribution plot:", output_path, "\n")
  
  invisible(p)
}

#' Plot WR percentiles (95th and 99th) by rate variation and n_inds
#'
#' Creates a 2-panel line plot showing WR percentiles (95th and 99th) where:
#' - Panels represent n_ind values (1, 2)
#' - X-axis shows duplication/loss rate
#' - Y-axis shows percentile values
#' - Colors represent RV (N, G, L)
#' - Line types represent ILS level (solid = high SF, dashed = low SF)
#' - Two lines per RV: one for 95th percentile, one for 99th percentile
#'
#' @param input_dir Directory containing findgraph CSV summary files
#' @param wr_column Column name to calculate percentiles for
#'   (e.g., "H0_best_tree_WR", "H1_best_graph_WR", "true_tree_wr")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @return ggplot object invisibly
plot_WR_percentiles_by_rate_variation <- function(
    input_dir,
    wr_column,
    output_dir,
    output_filename,
    plot_title = NULL,
    ymax = 10,
    rv_colors = get("rv_colors", envir = globalenv())) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }
  
  # Initialize data frame to store percentile results
  all_percentiles <- data.frame()
  
  # Process each file
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    
    # Extract parameters from filename
    # Example: "findgraph-DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5.csv"
    n_inds_match <- str_match(filename, "N_ind(\\d+)")
    sf_match <- str_match(filename, "SF([\\d.]+)")
    # Extract dup/loss rates (DUP and LOS followed by rate value)
    dup_match <- str_match(filename, "DUP([0-9eE.+-]+)")
    los_match <- str_match(filename, "LOS([0-9eE.+-]+)")
    rv_match <- str_match(filename, "-(RVG|RVL|RVN)-")
    
    if (any(is.na(
      c(n_inds_match, sf_match, dup_match, los_match, rv_match)))) {
      next
    }
    
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Check if WR column exists
      if (!(wr_column %in% names(df))) {
        cat("Warning: Column", wr_column, "not found in", filename, "\n")
        next
      }
      
      # Extract WR values and calculate percentiles
      wr_values <- df[[wr_column]]
      wr_values <- wr_values[!is.na(wr_values)]
      
      if (length(wr_values) == 0) next
      
      # Calculate percentiles
      pct_95 <- quantile(wr_values, 0.95, na.rm = TRUE)
      pct_99 <- quantile(wr_values, 0.99, na.rm = TRUE)
      
      # Get the rate - create a categorical label using the DUP rate
      dup_rate_str <- dup_match[2]
      
      # Create data frame row
      df_row <- data.frame(
        n_inds = as.numeric(n_inds_match[2]),
        SF = as.numeric(sf_match[2]),
        dup_loss_rate = dup_rate_str,
        ratevar = rv_match[2],
        pct_95 = as.numeric(pct_95),
        pct_99 = as.numeric(pct_99),
        stringsAsFactors = FALSE
      )
      
      all_percentiles <- rbind(all_percentiles, df_row)
      
    }, error = function(e) {
      cat("Error reading", filename, ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(all_percentiles) == 0) {
    cat("No valid WR data found\n")
    return(invisible(NULL))
  }
  
  cat("Total records in percentile data:", nrow(all_percentiles), "\n")
  
  # Map codes to readable names; prepare plot data (95th percentile only)
  plot_data <- all_percentiles %>%
    mutate(
      ratevar_name = recode(ratevar,
                           "RVG" = "gene",
                           "RVL" = "lineage",
                           "RVN" = "none"),
      ILS = ifelse(SF == 1.0, "high", "low"),
      n_ind_label = paste0("individuals / taxon = ", n_inds),
      dup_loss_rate_cat = factor(dup_loss_rate),
      value = pct_95
    ) %>%
    select(n_inds, n_ind_label, SF, ILS, dup_loss_rate, dup_loss_rate_cat, 
           ratevar, ratevar_name, value)
  
  # Order factors
  plot_data <- plot_data %>%
    mutate(
      n_ind_label = factor(n_ind_label, 
                          levels = c("individuals / taxon = 1",
                                    "individuals / taxon = 2")),
      ratevar_name = factor(ratevar_name, 
                           levels = c("none", "gene", "lineage")),
      dup_loss_rate_cat = factor(dup_loss_rate_cat),
      ILS = factor(ILS, levels = c("high", "low"))
    )
  
  # Create plot title if not provided
  if (is.null(plot_title)) {
    if (wr_column == "H0_best_tree_WR") {
      plot_title <- "H0 Best Tree WR 95-Percentiles by Rate Variation"
    } else if (wr_column == "H1_best_graph_WR") {
      plot_title <- "H1 Best Graph WR 95-Percentiles by Rate Variation"
    } else if (wr_column == "true_tree_wr") {
      plot_title <- "True Tree WR 95-Percentiles by Rate Variation"
    } else {
      plot_title <- paste(wr_column, "Percentiles by Rate Variation")
    }
  }
  
  # Define linetypes for ILS
  ils_linetypes <- c("high" = "solid",
                     "low" = "dashed")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = dup_loss_rate_cat, y = value,
                            color = ratevar_name,
                            linetype = ILS,
                            group = interaction(ratevar_name, ILS))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ n_ind_label, ncol = 2) +
    scale_color_manual(values = rv_colors,
                      name = "Substitution\nRate Variation") +
    scale_linetype_manual(values = ils_linetypes,
                         name = "ILS Level") +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(x = "duplication/loss rate",
         y = "95th percentile of worst f4 residual Z-score (WR)",
         title = plot_title) +
    theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x       = element_text(
        angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 16, color = "black", face = "bold"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text        = element_text(size = 14, face = "bold"),
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey60", linewidth = 0.7, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey82", linewidth = 0.4, linetype = "dotted"),
      legend.position    = "top",
      legend.background  = element_rect(
        color = "black", linewidth = 1.0, fill = "white"),
      legend.margin      = margin(5, 10, 5, 10),
      legend.text        = element_text(size = 16, face = "bold"),
      legend.title       = element_text(size = 16, face = "bold"),
      panel.spacing      = unit(1.5, "lines"),
      plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
  
  # Save the plot
  output_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(output_path, plot = p,
    width = 14, height = 10, dpi = 300, units = "in")

  cat("Saved WR percentiles plot:", output_path, "\n")
  
  invisible(p)
}

#' Plot WR (worst residual) 95th percentile – n_inds combined in one panel
#'
#' Creates a single-panel jitter scatter plot showing the 95th percentile WR
#' across all parameter combinations (n_inds = 1 and 2 together):
#' - X-axis shows duplication/loss rate
#' - Y-axis shows 95th percentile of worst f4 residual Z-score (WR) value
#' - Colors represent RV (none, gene, lineage)
#' - Shapes represent ILS level: circle (shape 21) = high SF=1.0,
#'                                triangle (shape 24) = low SF=0.5
#' - Fill represents n_inds: filled (interior = RV color) = n_inds=1,
#'                            hollow (interior = white)    = n_inds=2
#' - Points at the same dup/loss rate are jittered horizontally
#' - No connecting lines
#'
#' @param input_dir Directory containing findgraph CSV summary files
#' @param wr_column Column name to calculate percentiles for
#'   (e.g., "H0_best_tree_WR", "H1_best_graph_WR", "true_tree_wr")
#' @param output_dir Directory to save output figure
#' @param output_filename Output filename (without path, no extension)
#' @param plot_title Plot title (optional, auto-generated if NULL)
#' @param ymax Maximum y value for the axis (default: 105)
#' @param jitter_width Horizontal jitter width (default: 0.3)
#' @param rv_colors Named color vector for rate variation levels
#' @return ggplot object invisibly
plot_WR_percentiles_jitter_combined <- function(input_dir,
                          wr_column,
                          output_dir,
                          output_filename,
                          plot_title = NULL,
                          ymax = 11.5,
                          jitter_width = 0.2,
                          facet_label = NULL,
                          rv_colors = get("rv_colors", envir = globalenv())) {

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Find all CSV files
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    cat("No CSV files found in", input_dir, "\n")
    return(invisible(NULL))
  }

  # Initialize data frame to collect per-file percentiles
  all_percentiles <- data.frame()

  for (csv_file in csv_files) {
    filename <- basename(csv_file)

    # Extract parameters from filename
    # e.g. "findgraph-DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5-genelen1000.csv"
    n_inds_match <- str_match(filename, "N_ind(\\d+)")
    sf_match     <- str_match(filename, "SF([\\d.]+)")
    dup_match    <- str_match(filename, "DUP([0-9.eE+]+)")
    los_match    <- str_match(filename, "LOS([0-9eE.+-]+)")
    rv_match     <- str_match(filename, "-(RVG|RVL|RVN)-")

    if (any(is.na(
      c(n_inds_match, sf_match, dup_match, los_match, rv_match)))) next

    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE)

      if (!(wr_column %in% names(df))) {
        cat("Warning: Column", wr_column, "not found in", filename, "\n")
        next
      }

      wr_values <- df[[wr_column]]
      wr_values <- wr_values[!is.na(wr_values)]
      if (length(wr_values) == 0) next

      pct_95 <- as.numeric(quantile(wr_values, 0.95, na.rm = TRUE))

      all_percentiles <- rbind(all_percentiles, data.frame(
        n_inds       = as.numeric(n_inds_match[2]),
        SF           = as.numeric(sf_match[2]),
        dup_loss_rate = dup_match[2],
        ratevar      = rv_match[2],
        pct_95       = pct_95,
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      cat("Error reading", filename, ":", conditionMessage(e), "\n")
    })
  }

  if (nrow(all_percentiles) == 0) {
    cat("No valid WR data found\n")
    return(invisible(NULL))
  }

  cat("Total records in percentile data:", nrow(all_percentiles), "\n")

  # ── Prepare plot data ──────────────────────────────────────────────────────
  # 4-level combo: ILS (circle/triangle) × n_inds (filled/hollow)
  combo_levels <- c(
    "ILS = high, n = 1",
    "ILS = high, n = 2",
    "ILS = low, n = 1",
    "ILS = low, n = 2"
  )

  # Helper: format dup/loss rate as display label
  fmt_rate <- function(x) {
    v <- as.numeric(x)
    ifelse(v == 0, "0", sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
  }

  plot_data <- all_percentiles %>%
    mutate(
      ratevar_name = recode(ratevar,
                            "RVG" = "gene",
                            "RVL" = "lineage",
                            "RVN" = "none"),
      ILS = ifelse(SF == 1.0, "high", "low"),
      fill_var = ifelse(n_inds == 1,
                        recode(ratevar,
                               "RVG" = "gene",
                               "RVL" = "lineage",
                               "RVN" = "none"),
                        "hollow"),
      combo = case_when(
        n_inds == 1 & ILS == "high" ~ combo_levels[1],
        n_inds == 2 & ILS == "high" ~ combo_levels[2],
        n_inds == 1 & ILS == "low"  ~ combo_levels[3],
        n_inds == 2 & ILS == "low"  ~ combo_levels[4]
      ),
      dup_loss_rate_cat = factor(
        fmt_rate(dup_loss_rate),
        levels = fmt_rate(sort(unique(as.numeric(dup_loss_rate))))
      ),
      value = pct_95
    ) %>%
    mutate(
      ratevar_name = factor(ratevar_name,
        levels = c("none", "gene", "lineage")),
      fill_var     = factor(
        fill_var, levels = c("none", "gene", "lineage", "hollow")),
      combo        = factor(combo,         levels = combo_levels)
    )

  if (!is.null(facet_label)) {
    plot_data <- plot_data %>% mutate(panel_label = facet_label)
  }

  # ── Scale definitions ──────────────────────────────────────────────────────
  # Fill: RV colors for filled points (n_inds=1), white for hollow (n_inds=2)
  fill_values <- c(rv_colors, "hollow" = "white")

  # Shape: 4-level combo → fillable circle (21) or triangle (24)
  shape_values <- setNames(
    c(21, 21, 24, 24),
    combo_levels
  )

  # Fill overrides for the shape legend keys:
  # odd entries (filled) get a neutral grey; even entries (hollow) get white
  legend_fill_override <- c("grey55", "white", "grey55", "white")

  # ── Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot(plot_data,
              aes(x     = dup_loss_rate_cat,
                  y     = value,
                  color = ratevar_name,
                  fill  = fill_var,
                  shape = combo)) +
    geom_hline(data = data.frame(threshold = "y=3.0", yintercept = 3.0),
               aes(yintercept = yintercept, linetype = threshold),
               color = "red", linewidth = 1.5) +
    geom_hline(data = data.frame(threshold = "y=3.7", yintercept = 3.7),
               aes(yintercept = yintercept, linetype = threshold),
               color = "green4", linewidth = 1.5) +
    geom_jitter(size   = 6.5,
                stroke = 1.2,
                position = position_jitter(width = jitter_width, height = 0,
                                           seed  = 42)) +
    scale_color_manual(values = rv_colors,
                       name   = "rate variation:",
                       breaks = c("lineage", "gene", "none"),
                       labels = c("lineage" = "across lineages",
                                  "gene"    = "across genes",
                                  "none"    = "none")) +
    scale_fill_manual(values  = fill_values,
                      guide   = "none") +
    scale_shape_manual(values = shape_values,
                       name   = "ILS & individuals / taxon (n)") +
    scale_linetype_manual(
      name   = "WR threshold:",
      values = c("y=3.0" = "dashed", "y=3.7" = "dashed"),
      breaks = c("y=3.7", "y=3.0")
    ) +
    scale_y_continuous(
      limits = c(0, ymax),
      expand = c(0, 0),
      breaks = sort(unique(c(seq(0, ymax, by = 2), 3.0, 3.7))),
      labels = function(x) {
        ifelse(abs(x - 3.7) < 1e-8 | abs(x - 3.0) < 1e-8,
              "",
              as.character(x))
      }
    )  +               
    labs(
      x     = "duplication and loss rate",
      y     = "95th percentile of worst f4 residual Z-score (WR)",
      title = plot_title
    ) +
    theme_bw() +
    theme(
      # Title & caption
      plot.title        = element_text(hjust = 0.5, size = 22, face = "bold"),
      plot.caption      = element_text(
        size = 11, hjust = 0.5, color = "grey25", face = "italic"),
      # Axis titles – bold, large, with breathing room
      axis.title.x      = element_text(
        size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(
        size = 16, face = "bold", margin = margin(r = 10)),
      # Axis tick labels
      axis.text.x       = element_text(size = 15, angle = 45, 
                        hjust = 1, color = "black", face = "bold"),
      axis.text.y       = element_text(
        size = 15, color = "black", face = "bold"),
      # Prominent axis ticks
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      # Panel border – thick and clearly black
      panel.border      = element_rect(
        color = "black", linewidth = 1.8, fill = NA),
      panel.grid.major  = element_line(
        color = "grey87", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor  = element_line(
        color = "grey93", linewidth = 0.3, linetype = "dotted"),
      # Legend – prominent box, larger keys and text
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.box.spacing = unit(0.4, "lines"),
      legend.background  = element_blank(),
      legend.margin      = margin(5, 10, 5, 10),
      legend.key.size    = unit(2.0, "lines"),
      legend.text        = element_text(size = 14, face = "bold"),
      legend.title       = element_text(size = 14, face = "bold"),
      plot.margin        = unit(c(0.5, 0.5, 0.8, 0.5), "lines")
    ) +
    guides(
      color    = guide_legend(override.aes = list(size = 8.0, stroke = 1.6),
                              nrow = 3, position = "right", order = 1),
      shape    = guide_legend(override.aes = list(size   = 8.0,
                                                  stroke = 1.6,
                                                  fill   = legend_fill_override,
                                                  color  = "black"),
                              nrow = 2, position = "bottom",
                              theme = theme(legend.background = element_rect(
                                color = "black", linewidth = 1.0,
                                fill = "white"))),
      linetype = guide_legend(
        override.aes = list(color = c("green4", "red"), linewidth = 1.5),
        nrow = 2, position = "right", order = 2)
    )

  if (!is.null(facet_label)) {
    p <- p +
      facet_wrap(~ panel_label) +
      theme(strip.text       = element_text(size = 19, face = "bold"),
            strip.background = element_rect(fill = "grey90"))
  }

  # ── Save ───────────────────────────────────────────────────────────────────
  output_path <- file.path(output_dir, paste0(output_filename, ".pdf"))
  ggsave(output_path, plot = p, width = 8, height = 8, dpi = 300, units = "in")

  cat("Saved WR percentiles jitter combined plot:", output_path, "\n")

  invisible(p)
}


# ─────────────────────────────────────────────────────────────────────────
# LEGACY (disabled for Dryad / paper submission, 2026-05).
# plot_taxon_recovery_heatmaps() consumes results/*_taxon_recovery.csv,
# which are no longer produced by the pipeline (the upstream code in
# scripts/summary_findgraph.jl and scripts/summary_snaq.jl is commented
# out). This function is retained as inert reference and will stop() if
# called. To re-enable: remove the stop() at the top of the function body,
# then re-enable the upstream taxon-recovery code.
# ─────────────────────────────────────────────────────────────────────────
#' Plot taxon recovery heatmaps by role and parameter setting
#'
#' Generates two figures showing taxon recovery rates:
#' - Figure 1: Per-parameter-setting heatmaps in a grid (3 columns)
#' - Figure 2: Aggregated heatmaps by rate variation (3 panels side-by-side)
#'
#' Heatmaps show normalized recovery rates (count/n_replicates) for each taxon
#' (rows: A, B, C, D, F, G, H) across roles
#' (columns: hybrid_taxon, major_donor, minor_donor).
#'
#' @param input_long Path to long-format CSV with columns:
#'   param_setting, taxon, role, count_recovered, count_not_recovered,
#'   total, pct_correct
#' @param input_dir Directory containing original summary CSV files
#'   (to determine n_replicates)
#' @param output_dir Directory to save output PDF figures
#' @param colormap Name of colormap (default: "Blues")
#' @param annotation_digits Decimal places for cell annotations (default: 2)
#'
#' @details
#' Extracts rate variation condition from param_setting:
#'   - contains "RVL" → "Lineage (RVL)"
#'   - contains "RVG" → "Gene-specific (RVG)"
#'   - otherwise     → "No variation (RVN)"
#'
#' Normalized value = count_recovered / n_replicates (can exceed 1.0)
#'
#' @return Invisibly returns NULL. Saves two PDF figures:
#'   - heatmap_per_param.pdf   : Individual heatmaps per parameter setting
#'   - heatmap_aggregated.pdf  : Aggregated by rate variation
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import stringr
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' \dontrun{
#' plot_taxon_recovery_heatmaps(
#'   input_long = "results/findgraph_taxon_recovery.csv",
#'   input_dir = "results",
#'   output_dir = "visualization_results/taxon_recovery"
#' )
#' }
#'
plot_taxon_recovery_heatmaps <- function(input_long,
                                        input_dir,
                                        output_dir,
                                        colormap = "Blues",
                                        annotation_digits = 2) {

  stop("LEGACY: plot_taxon_recovery_heatmaps() is disabled for Dryad/paper ",
       "submission (2026-05). See the banner comment above the function for ",
       "how to re-enable.")

  # Fixed taxon and role orders
  taxon_order <- c("A", "B", "C", "D", "F", "G", "H")
  role_order <- c("hybrid_taxon", "major_donor", "minor_donor")
  
  # Create output directory if needed
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Loading long-format taxon recovery CSV...\n")
  df_long <- read_csv(input_long, show_col_types = FALSE) %>%
    mutate(
      taxon = as.character(taxon),
      role = as.character(role)
    )
  
  # ──────────────────────────────────────────────────────────────────────────
  # Extract parameters and normalize
  # ──────────────────────────────────────────────────────────────────────────
  cat("Extracting parameters and normalizing recovery rates...\n")
  
  # Count n_replicates per param_setting
  replica_counts <- df_long %>%
    group_by(param_setting) %>%
    summarise(n_replicates = max(total, na.rm = TRUE), .groups = "drop")
  
  df_prep <- df_long %>%
    left_join(replica_counts, by = "param_setting") %>%
    mutate(
      # Extract dup/loss rate (supports 0.0 and scientific notation like 3e-4)
      dup_loss_rate = as.numeric(str_extract(
        param_setting, "(?<=DUP)[0-9.]+(?:[eE][+-]?[0-9]+)?")),
      # Extract rate variation type
      rate_var = case_when(
        str_detect(param_setting, "RVL") ~ "Lineage (RVL)",
        str_detect(param_setting, "RVG") ~ "Gene-specific (RVG)",
        TRUE                             ~ "No variation (RVN)"
      ),
      # Extract SF (ILS level)
      SF = as.numeric(str_extract(param_setting, "(?<=SF)[0-9.]+")),
      # Extract N_ind
      N_ind = as.integer(str_extract(param_setting, "(?<=N_ind)[0-9]+")),
      # ILS category
      ILS_cat = if_else(SF < 1.0, "low ILS", "high ILS"),
      # ILS row label for faceting
      ILS_row = factor(paste0(ILS_cat, ", n_ind=", N_ind),
                       levels = c("low ILS, n_ind=1", "high ILS, n_ind=1", 
                                 "low ILS, n_ind=2", "high ILS, n_ind=2")),
      # Normalize: count / n_replicates
      n_replicates = if_else(
        is.na(n_replicates), as.numeric(total), n_replicates),
      norm_value = count_recovered / n_replicates
    ) %>%
    filter(!is.na(dup_loss_rate), !is.na(N_ind), !is.na(SF)) %>%
    select(param_setting, taxon, role, norm_value, rate_var,
      dup_loss_rate, N_ind, ILS_cat, ILS_row)
  
  cat("Data prepared successfully.\n")
  
  # ──────────────────────────────────────────────────────────────────────────
  # FIGURE 1: Per-parameter-setting heatmaps using facet_grid
  # ──────────────────────────────────────────────────────────────────────────
  cat("Creating Figure 1: Per-parameter-setting heatmaps using facet_grid...\n")
  
  # Compute global min/max
  vmin <- 0.0
  vmax <- max(df_prep$norm_value, na.rm = TRUE)
  
  # Prepare heatmap data: complete all taxon x role combos per param_setting
  df_fig1 <- df_prep %>%
    # Keep metadata for each param_setting
    distinct(param_setting, rate_var, dup_loss_rate, ILS_row) %>%
    # Expand to all taxon x role combinations
    expand_grid(taxon = factor(taxon_order, levels = taxon_order),
                role = factor(role_order, levels = role_order)) %>%
    # Left join with actual data
    left_join(
      df_prep %>% select(param_setting, taxon, role, norm_value),
      by = c("param_setting", "taxon", "role")
    ) %>%
    # Fill missing values with 0
    mutate(norm_value = replace_na(norm_value, 0)) %>%
    mutate(
      label = sprintf(paste0("%.", annotation_digits, "f"), norm_value),
      dup_label = case_when(
        dup_loss_rate < 0.001 ~ sprintf("DUP/LOS = %.1e", dup_loss_rate),
        TRUE ~ sprintf("DUP/LOS = %.4f", dup_loss_rate)
      )
    )
  
  # Create the faceted heatmap
  p1 <- ggplot(df_fig1, aes(x = role, y = taxon, fill = norm_value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = 2, color = "black") +
    scale_fill_distiller(
      palette = colormap, direction = 1, limits = c(vmin, vmax),
      guide = "none") +
    scale_y_discrete(limits = rev(taxon_order)) +
    facet_grid(ILS_row ~ dup_label + rate_var,
      scales = "free_x", space = "free_x",
      labeller = labeller(.rows = label_value, .cols = label_value)) +
    labs(
      title = "Taxon recovery by parameter setting", x = "Role", y = "Taxon") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y = element_text(size = 10, margin = margin(r = 10)),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.3, fill = NA),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(angle = 0),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
    )
  
  # Save Figure 1
  ggsave(file.path(output_dir, "heatmap_per_param.pdf"),
         plot = p1, width = 16, height = 10, dpi = 300, units = "in")
  cat("Saved Figure 1: heatmap_per_param.pdf\n")
  
  # ──────────────────────────────────────────────────────────────────────────
  # FIGURE 2: Aggregated by rate variation (3 panels side-by-side)
  # ──────────────────────────────────────────────────────────────────────────
  cat("Creating Figure 2: Aggregated heatmaps by rate variation...\n")
  
  # Define rate variation order with custom titles
  rate_var_order <- c(
    "No variation (RVN)", "Gene-specific (RVG)", "Lineage (RVL)")
  rate_var_titles <- c(
    "No variation (RVN)" = "No rate variation",
    "Gene-specific (RVG)" = "Variation across genes",
    "Lineage (RVL)" = "Variation across lineage"
  )
  
  # Aggregate: mean and std across all param_settings
  df_agg <- df_prep %>%
    group_by(rate_var, taxon, role) %>%
    summarise(
      mean_val = mean(norm_value, na.rm = TRUE),
      std_val = sd(norm_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(
      rate_var = factor(rate_var_order, levels = rate_var_order),
      taxon = factor(taxon_order, levels = taxon_order),
      role = factor(role_order, levels = role_order),
      fill = list(mean_val = 0, std_val = 0)
    ) %>%
    mutate(
      rate_var = factor(rate_var, levels = rate_var_order),
      label = case_when(
        std_val > 0.05 ~ sprintf(
          paste0("%.", annotation_digits, "f\n±%.", annotation_digits, "f"),
          mean_val, std_val),
        TRUE ~ sprintf(paste0("%.", annotation_digits, "f"), mean_val)
      )
    )
  
  # Create faceted heatmap
  p2 <- ggplot(df_agg, aes(x = role, y = taxon, fill = mean_val)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_distiller(
      palette = colormap, direction = 1, limits = c(vmin, vmax),
      name = "Mean norm\ncount") +
    scale_y_discrete(limits = rev(taxon_order)) +
    facet_wrap(~rate_var, nrow = 1,
      labeller = labeller(rate_var = rate_var_titles)) +
    labs(
      title = "Taxon recovery aggregated by rate variation",
      x = "Role", y = "Taxon") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, margin = margin(r = 10)),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
      strip.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      aspect.ratio = 1,
      legend.position = "right"
    )
  
  # Save Figure 2
  ggsave(file.path(output_dir, "heatmap_aggregated.pdf"),
         plot = p2, width = 14, height = 5.5, dpi = 300, units = "in")
  cat("Saved Figure 2: heatmap_aggregated.pdf\n")
  
  cat("Done! Figures saved to:", output_dir, "\n")
  
  invisible(NULL)
}
