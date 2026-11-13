
# find_graphs for one replicate, called by findgraphs.jl (not run directly)
# f2 statistics from the eigenstrat files, then find_graphs() of admixtools
# with num_admix (h) admixture events, several independent runs. The graphs
# are filtered by log likelihood, duplicates removed (same hash), then
# qpgraph gives the worst residual of each graph. For h = 0 we also check
# if the true species tree is among the graphs found.
# Input: <input_dir>/<prefix>.geno .snp .ind, and the seed file
# Output: <output_dir>/rep<id>_admix<k>_unique_graphs.rds,
#   rep<id>_admix<k>_summary_table.txt and rep<id>_f2.rds (h = 0 only)
# main steps:
# 1. compute f2 statistics from input data (h = 0), or load them (h > 0)
# 2. run find_graphs() several times (runs), one seed per run
# 3. select graphs by top-n ranking or ll threshold (first filter per run)
# 4. deduplicate topologies via hashing (remove duplicates across all runs)
# 5. apply second filter using selection method on unique graphs
# 6. evaluate graphs with qpgraph (likelihood, residuals)
# 7. compare against the true tree when h=0
# the true tree is included only when h=0, for h>0 we only look at the fit
# example:
# Rscript findgraphs_1rep.R --input_dir data --output_dir results \
#   --prefix prefix --num_admix 0 --runs 100 --rep_id 1 \
#   --selection_method keep_five_lowest_ll_graphs --seed_file_path seeds.txt


library(admixtools)
library(optparse)
library(dplyr)
library(tidyverse)
library(Rcpp)
library(igraph)
# needs igraph 1.x (tested with 1.6.0): in igraph 2.x get.edge.ids() does not
# take igraph.vs objects any more, and admixtools:::graph_to_pwts passes one,
# so igraph >= 2.0 fails at run time.
# to install: remotes::install_version('igraph', version='1.6.0')

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Input directory"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
              help = "Input file prefix"),
  # Parameters for find_graphs 
  # num_admix: number of admixture events
  make_option(c("-k", "--num_admix"), type = "integer", default = 0,
              help = "Number of admixture events"),
  # stop_gen: number of generations to stop find_graphs
  make_option("--stop_gen", type = "integer", default = 100,
              help = "Number of generations to stop find_graphs"),
  # outgroup: outgroup population name          
  make_option("--outgroup", type = "character", default = "A",
              help = "Outgrosup (Default Homo \"A\")"),
  # runs: number of independent runs
  make_option("--runs", type = "integer", default = 100,
              help = "Number of times to run find_graphs"),
  make_option(c("-b", "--blgsize"), type = "integer", default = 1000,
              help = "Block sizes used in find_graphs"),
  make_option("--maxmiss", type = "numeric", default = 1.0,
              help = "Max fraction of missing data per SNP (default: 1.0)"),
  make_option(c("-r", "--rep_id"), type = "character",
              help = "Replication ID for printing and result saving"),
  make_option(c("-s", "--seed_file_path"), type = "character",
              help = "A path to the seed file to be used"),
  make_option("--output_graph_suffix", type = "character",
              default = "_graphs.rds",
              help = "Suffix to the output graphs files"),
  make_option("--output_f2_suffix", type = "character",
              default = "_f2.rds",
              help = "Suffix to the output graphs files"),
  make_option("--output_summary_table_suffix", type = "character",
              default = "_summary_table.txt",
              help = "Suffix to the table outputs"),
  make_option("--true_tree_newick", type = "character",
              default = "(A,((((B,C),(D,E)),F),(G,H)));",
              help = "True species tree in Newick format"),
  make_option("--rootfolder", type = "character",
              default = "__use_cwd__",
              help = "rootfolder to save any warnings"),
  make_option("--selection_method", type = "character",
              default = NULL,
              help = "Method to select graphs"),
  make_option("--threshold", type = "integer", 
              default = 2,
              help = "Threshold for graph selection, default: 2"),
  make_option("--debug_mode", action = "store_true", default = FALSE,
              help = "Enable debug mode to save individual run results")
              ) 

#-------------- Parse options --------------#
opt <- parse_args(OptionParser(option_list = option_list))

# ----------- set up variables --------------#
input_dir <- opt$input_dir
output_dir <- opt$output_dir
prefix <- opt$prefix
outgroup <- opt$outgroup
num_admix <- opt$num_admix
stop_gen <- opt$stop_gen
runs <- opt$runs
rep_id <- opt$rep_id
blgsize <- opt$blgsize
maxmiss <- opt$maxmiss
seed_file_path <- opt$seed_file_path
output_graph_suffix <- opt$output_graph_suffix
output_f2_suffix <- opt$output_f2_suffix
output_summary_table_suffix <- opt$output_summary_table_suffix
true_tree_newick <- opt$true_tree_newick
selection_method <- opt$selection_method 
threshold <- opt$threshold 
debug_mode <- opt$debug_mode

if (is.null(threshold) || is.na(threshold) || !is.numeric(threshold)) {
  stop("threshold parameter is NULL, NA, or not numeric")
}
if (is.null(num_admix) || is.na(num_admix)) {
  stop("num_admix parameter is NULL or NA")
} 

# selection_method: how we want to select graphs generated from find_graphs 
# select_threshold_within_ll: 
#      -> select graphs within a threshold of the lowest ll score from each run 
# keep_five_lowest_ll_graphs:
#      -> keep the five lowest ll graphs from each run 
if (selection_method == "select_threshold_within_ll" &&
    (is.null(threshold) || threshold <= 0)) {
  stop("Please provide a valid positive threshold value for graph selection.")
} else if (selection_method == "select_threshold_within_ll" &&
    !is.null(threshold)) {
  message(paste0("Select graphs within ", threshold,
  " of the best graph with lowest ll score."))
} else if (selection_method == "keep_five_lowest_ll_graphs") {
  message("Keep the five lowest ll graphs from each run.") 
} else if (selection_method == "both") {
  message(paste0("Select graphs within ", threshold,
    " of best graph and keep five lowest ll graphs per run."))
} else {
  stop(paste0("Invalid selection method. Choose ",
    "'select_threshold_within_ll', 'keep_five_lowest_ll_graphs', or 'both'."))
}

# Debug mode message
if (debug_mode) {
  message("DEBUG MODE: Individual run results will be saved to .rds files")
} else if (selection_method == "both") {
  message("Individual run results saved to .rds files per selection method")
}

rootfolder <- opt$rootfolder # set up rootfolder 
if (rootfolder == "__use_cwd__") {
  rootfolder <- getwd()
}
# Set up paths:
# input_prefix: the input prefix for f2_from_geno 
input_prefix <- file.path(input_dir, prefix) 

# output_dir: the output directory to save results 
output_graph_file <- file.path(output_dir, 
  paste0("rep", rep_id, "_admix", num_admix, output_graph_suffix)) # nolint

# output_graph_f2: the output file for f2 statistics
# This is the same for both num_admix = 0 and num_admix > 0 
output_graph_f2 <- file.path(output_dir, 
  paste0("rep", rep_id, output_f2_suffix))  # nolint

# output_summary_table: the output file for summary table 
output_summary_table <- file.path(output_dir,
  paste0("rep", rep_id, "_admix", num_admix, # nolint
         output_summary_table_suffix))

#------------ set up seeds --------------#
# seed matrix: one row per replicate, one column per run
seed_matrix <- read.table(seed_file_path, header = TRUE)

rep_index <- as.integer(rep_id)
if (rep_index > nrow(seed_matrix)) {
  stop(paste("Replicate ID", rep_id, 
            "exceeds number of available seed rows:", nrow(seed_matrix)))
}

replicate_seeds <- as.numeric(seed_matrix[rep_index, ])
if (length(replicate_seeds) < runs) {
  stop(paste("Not enough seeds for", runs, 
            "runs. Available seeds:", length(replicate_seeds)))
}

cat("Using", length(replicate_seeds), "seeds for replicate", rep_id, "\n")
cat("First few seeds:", head(replicate_seeds, 5), "\n")

#-------------- Print analysis parameters --------------#
cat("\n=== Analysis Parameters ===\n")
cat("Replicate ID:", rep_id, "\n")
cat("Number of admixture events:", num_admix, "\n")
cat("Number of runs:", runs, "\n")
cat("Block size:", blgsize, "\n")
cat("Max missing data per SNP:", maxmiss, "\n")
cat("Selection method:", selection_method, "\n")
if (!is.null(threshold)) cat("Threshold:", threshold, "\n")
cat("===========================\n\n")

# snp_file <- file.path(input_dir, paste0(prefix, ".snp"))
# snp_ids <- read.table(snp_file, stringsAsFactors = FALSE)[[1]]

#-------------- Create f2 statistics --------------#
# f2 are computed when num_admix = 0 and saved, then loaded for num_admix > 0
if (num_admix == 0) {
  # adjust_pseudohaploid: TRUE because the data are haploid
  # blgsize: block size for the jackknife
  f2 <- f2_from_geno(pref = input_prefix,
                    adjust_pseudohaploid = TRUE,
                    blgsize = blgsize,
                    remove_na = FALSE # safety check: no missing data
                    )
                    
  # number of blocks left after filtering
  num_blocks <- dim(f2)[3] 
  cat("Number of blocks in f2 statistics:", num_blocks, "\n")
  saveRDS(f2, file = output_graph_f2)
} else {
  f2_file_for_admix <- file.path(output_dir, 
    paste0("rep", rep_id, output_f2_suffix))  # nolint
  if (!file.exists(f2_file_for_admix)) {
    stop(paste("F2 statistics file for replicate", 
                rep_id, 
                "with admix0 not found at:", 
                f2_file_for_admix))
  }
  f2 <- readRDS(f2_file_for_admix)
  cat("Loaded pre-computed f2 statistics from:", f2_file_for_admix, "\n")
  num_blocks <- dim(f2)[3]
  cat("Number of blocks in loaded f2 statistics:", num_blocks, "\n")
}

#-------------- Create true species tree graph --------------#
cat("Tree topology:", true_tree_newick, "\n")

if (is.null(true_tree_newick) || length(true_tree_newick) == 0 ||
    true_tree_newick == "" || nchar(true_tree_newick) == 0) {
  stop("true_tree_newick is null, empty, or invalid")
}

# newick -> edge matrix (admixtools) -> igraph, then its hash
edges <- newick_to_edges(true_tree_newick)
true_tree_graph <- graph_from_edgelist(edges, directed = TRUE)
cat("Change a newick into a graph:\n")

true_tree_hash <- graph_hash(true_tree_graph)
cat("True species tree hash:", true_tree_hash, "\n")

#-------------- Run find_graphs --------------#
# summary table file name with the method name before the extension,
# for selection_method = "both"
get_summary_table_filename <- function(output_dir, rep_id, num_admix,
    method_name, output_summary_table_suffix) {
  if (grepl("\\.", output_summary_table_suffix)) {
    parts <- strsplit(output_summary_table_suffix, "\\.")[[1]]
    name_part <- paste(parts[-length(parts)], collapse = ".")
    extension <- parts[length(parts)]
    filename <- paste0("rep", rep_id, "_admix", num_admix, 
                        name_part, "_", method_name, ".", extension)
  } else {
    filename <- paste0("rep", rep_id, "_admix", num_admix, 
                        output_summary_table_suffix, "_", method_name)
  }
  return(file.path(output_dir, filename))
}

# run find_graphs() "runs" times on f2, one seed per run, and keep the best
# graphs of each run (first filter)
# Input:
#   f2: f2 object from f2_from_geno()
#   num_admix: number of admixture events (h)
#   stop_gen, outgroup: passed to find_graphs()
#   runs, replicate_seeds: number of runs and their seeds
#   selection_method, threshold: how to pick the graphs of each run
#   debug_mode, output_dir: save the graphs of each run as rds
# Output: list by run ("run1", "run2", ...), each a list of graphs with
#   $graph (igraph), $edges, $score (ll), $rank (in the run), $run_id
#   for selection_method = "both": list(main, keep_five_lowest_ll_graphs,
#   select_threshold_within_ll)
# e.g. graph_lst[["run1"]][[1]]$graph is the best graph of run 1
run_find_graphs_replicate <- function(f2,
                                      num_admix,
                                      stop_gen,
                                      outgroup,
                                      runs,
                                      rep_id,
                                      replicate_seeds,
                                      selection_method,
                                      threshold = NULL,
                                      debug_mode = FALSE,
                                      output_dir = NULL) {

  # the 5 graphs with the lowest score of one run
  select_top5_graphs <- function(opt_results, run_id) {
    Top5_graphs_per_run <- opt_results %>% 
                  dplyr::slice_min(score, n = 5, with_ties = TRUE) 
    top5_graphs_list_per_run <- lapply(
        1:nrow(Top5_graphs_per_run), function(j) {
          list(graph = Top5_graphs_per_run$graph[[j]], 
              edges = Top5_graphs_per_run$edges[[j]],
              score = Top5_graphs_per_run$score[j],
              rank = j, # rank: rank of the graph based on score in each run
              run_id = run_id)})
    return(top5_graphs_list_per_run)
  }
  
  # the graphs within threshold of the lowest score of one run
  select_threshold_graphs <- function(opt_results, threshold, run_id) {
    min_score <- min(opt_results$score)
    threshold_score <- min_score + threshold
    selected_graphs <- opt_results %>%
      dplyr::filter(score <= threshold_score) %>%
      dplyr::arrange(score) # Sort by score
    selected_graphs_list <- lapply(1:nrow(selected_graphs), function(j) {
      list(graph = selected_graphs$graph[[j]],
          edges = selected_graphs$edges[[j]], 
          score = selected_graphs$score[j],
          rank = j, # rank: rank of the graph based on score in each run
          run_id = run_id)
    })
    return(selected_graphs_list)
  }
  
  # save the graphs of one run as rds (debug mode)
  save_debug_results <- function(graph_list, 
                                  method_name, 
                                  run_id, rep_id, 
                                  num_admix, output_dir) {
    if (!is.null(output_dir)) {
      debug_filename <- file.path(output_dir, 
        paste0("rep", rep_id, "_admix", num_admix, "_", 
                    method_name, "_run", run_id, ".rds"))
      saveRDS(graph_list, file = debug_filename)
      cat("DEBUG: Saved run", run_id, method_name, "to:", debug_filename, "\n")
    }
  }

  graph_lst <- list()
  
  for (i in 1:runs) {
    current_seed <- replicate_seeds[i]

    set.seed(current_seed)
    message(paste0("Running find_graphs for rep ", rep_id, " at run ", 
                    i, " with seed ", current_seed)) # nolint

    opt_results <- find_graphs(f2, 
                              numadmix = num_admix,
                              stop_gen = stop_gen,
                              outpop = outgroup
                              )

    if (selection_method == "keep_five_lowest_ll_graphs") {
      top5_graphs_list_per_run <- select_top5_graphs(opt_results, i)
      graph_lst[[paste0("run", i)]] <- top5_graphs_list_per_run
      if (debug_mode) {
        save_debug_results(top5_graphs_list_per_run,
          "keep_five_lowest_ll_graphs", i, rep_id, num_admix, output_dir)
      }

    } else if (selection_method == "select_threshold_within_ll") {
      selected_graphs_list <- select_threshold_graphs(opt_results, threshold, i)
      graph_lst[[paste0("run", i)]] <- selected_graphs_list
      if (debug_mode) {
        save_debug_results(selected_graphs_list, "select_threshold_within_ll", 
                            i, rep_id, num_admix, output_dir)
      }

    } else if (selection_method == "both") {
      top5_graphs_list_per_run <- select_top5_graphs(opt_results, i)
      selected_graphs_list <- select_threshold_graphs(opt_results, threshold, i)
      if (!exists("graph_lst_method1")) graph_lst_method1 <- list()
      if (!exists("graph_lst_method2")) graph_lst_method2 <- list()
      graph_lst_method1[[paste0("run", i)]] <- top5_graphs_list_per_run
      graph_lst_method2[[paste0("run", i)]] <- selected_graphs_list
      graph_lst[[paste0("run", i)]] <- top5_graphs_list_per_run # the main list
      # with "both", the graphs of each run are always saved
      save_debug_results(top5_graphs_list_per_run,
        "keep_five_lowest_ll_graphs", i, rep_id, num_admix, output_dir)
      save_debug_results(selected_graphs_list, "select_threshold_within_ll", 
                        i, rep_id, num_admix, output_dir)

    } else { # This is checked before but just in case 
      stop(paste0("Invalid selection method. Choose ",
        "'select_threshold_within_ll',",
        " 'keep_five_lowest_ll_graphs', or 'both'."))
    }
  }
  
  cat("Finishing findgraphs for rep", rep_id, "with", runs, "runs.\n") # nolint
  
  if (selection_method == "both") {
    return(list(
      main = graph_lst,
      keep_five_lowest_ll_graphs = graph_lst_method1,
      select_threshold_within_ll = graph_lst_method2
    ))
  } else {
    return(graph_lst)
  }
}

graph_lst <- run_find_graphs_replicate(f2 = f2,
                                      num_admix = num_admix,
                                      stop_gen = stop_gen,
                                      outgroup = outgroup,
                                      runs = runs,
                                      rep_id = rep_id,
                                      replicate_seeds = replicate_seeds,
                                      selection_method = selection_method,
                                      threshold = threshold,
                                      debug_mode = debug_mode,
                                      output_dir = output_dir
                                      )

#-------------- Remove duplicated graphs --------------#
# add $hash (graph_hash of $graph) to every graph of graph_lst
# graphs with the same hash have the same topology
add_hashes_to_graph_lst <- function(graph_lst) {
  for (run_name in names(graph_lst)) {
    for (i in seq_along(graph_lst[[run_name]])) {
      graph_obj <- graph_lst[[run_name]][[i]]$graph
      graph_hash_val <- graph_hash(graph_obj)
      graph_lst[[run_name]][[i]]$hash <- graph_hash_val
    }
  }
  return(graph_lst)
}

# second filter, after removing the duplicates: keep the 5 lowest scores,
# or the graphs within threshold of the lowest score, over all runs together
# Input: unique_graphs (list), selection_method, threshold
# Output: the filtered list
apply_second_filter <- function(unique_graphs, selection_method,
    threshold = NULL) {
  if (length(unique_graphs) == 0) {
    cat("No graphs to filter - returning empty list\n")
    return(unique_graphs)
  }
  
  cat("Applying second filter:", selection_method,
      "to", length(unique_graphs), "unique graphs\n")
  
  scores <- sapply(unique_graphs, function(g) g$score)
  sorted_indices <- order(scores)
  sorted_graphs <- unique_graphs[sorted_indices]
  
  if (selection_method == "keep_five_lowest_ll_graphs") {
    n_keep <- min(5, length(sorted_graphs))
    cat("Keeping", n_keep, "lowest scoring graphs\n")
    return(sorted_graphs[1:n_keep])
    
  } else if (selection_method == "select_threshold_within_ll") {
    if (is.null(threshold)) {
      stop("Threshold must be provided for 'select_threshold_within_ll' method")
    }
    min_score <- min(scores)
    threshold_score <- min_score + threshold
    selected_indices <- which(scores <= threshold_score)
    cat("Keeping", length(selected_indices), "graphs within threshold",
        threshold, "of best score", min_score,
        "(threshold score:", threshold_score, ")\n")
    return(unique_graphs[selected_indices])
    
  } else { # "both": nothing more to do here
    cat("No additional filtering - keeping all",
        length(unique_graphs), "graphs\n")
    return(unique_graphs)
  }
}

# remove the graphs with the same hash (same topology) across runs
# two graphs with the same hash should have the same score, if they differ
# by more than score_tol a warning is written to findgraph_warnings.txt
# Input: graph_lst from run_find_graphs_replicate, rootfolder (for the
#   warning file), score_tol (default 1e-6)
# Output: flat list of the unique graphs
remove_duplicate_graphs <- function(graph_lst, rootfolder, score_tol = 1e-6) {
  
  graph_lst <- add_hashes_to_graph_lst(graph_lst)

  hash_seen <- list()
  unique_graph_lst <- list()

  warning_log_path <- file.path(rootfolder, "findgraph_warnings.txt")
  warning_file <- NULL # opened only if there is a warning

  for (run_name in names(graph_lst)) {
    for (graph_info in graph_lst[[run_name]]) {
      h <- graph_info$hash
      score <- graph_info$score
      rank <- graph_info$rank

      if (!h %in% names(hash_seen)) {
        # graph_info$run_id <- run_name
        hash_seen[[h]] <- list(score = score, graph_info = graph_info)
        unique_graph_lst[[length(unique_graph_lst) + 1]] <- graph_info
      } else {
        prev_score <- hash_seen[[h]]$score
        if (abs(prev_score - score) > score_tol) {
          warning_msg <- paste0(
            "WARNING: Graph from rep ", rep_id, "run ", run_name,
            " with rank ", rank,
            " has identical hash with different scores: ",
            prev_score, " vs ", score,
            " (difference = ", abs(prev_score - score), ")."
          )
          if (is.null(warning_file)) {
            warning_file <- file(warning_log_path, open = "at")
          }
          writeLines(warning_msg, con = warning_file)
          warning(warning_msg)
        } else {
          msg <- paste0("For ", run_name, " rank ", rank,
                        " graph is duplicated and removed.")
          message(msg)
        }
      }
    }
  }
  # Only save the warning_log if there is some warning else no save
  if (!is.null(warning_file)) close(warning_file) 
  return(unique_graph_lst)
}

#-------------- Run qpgraphs & Organize Outputs --------------#
# qpgraph on each unique graph (and on the true tree when given): score,
# worst residual, is it the true tree, gammas when num_admix = 1
# Input: unique_graphs (list), f2, output_file_path, num_admix,
#   true_tree_graph and true_tree_hash (NULL for num_admix > 0), num_blocks
# Output: the summary table is written to output_file_path (true tree first,
#   then the graphs ranked by score), returns unique_graphs with $wr added
run_qpgraph <- function(unique_graphs, 
                        f2, 
                        output_file_path, 
                        num_admix = 0, 
                        true_tree_graph = NULL, 
                        true_tree_hash = NULL,
                        num_blocks = NA){

  true_tree_result <- NULL
  if (!is.null(true_tree_graph)) {
    cat("Processing true species tree...\n")
    true_tree_qpgraph_result <- qpgraph(f2, true_tree_graph,
      return_fstats = TRUE)
    true_tree_wr <- true_tree_qpgraph_result$worst_residual
    true_tree_WR_smaller_than_3 <- true_tree_wr <= 3
    score <- true_tree_qpgraph_result$score

    cat("===== Do we find the true tree =====\n")
    cat("True tree score:", score, "\n")
    cat("True tree worst residual:", true_tree_wr, "\n")
    cat("True tree qpgraph result structure:\n")
    str(true_tree_qpgraph_result)
    if (!is.null(true_tree_qpgraph_result$f2)) {
      f2fit <- range(true_tree_qpgraph_result$f2$fit, na.rm = TRUE)
      cat("True tree F2 fitted range:", f2fit, "\n")
    }
    if (!is.null(true_tree_qpgraph_result$f4)) {
      f4fit <- range(true_tree_qpgraph_result$f4$fit, na.rm = TRUE)
      f4diff <- range(true_tree_qpgraph_result$f4$diff, na.rm = TRUE)
      cat("True tree F4 fitted range:", f4fit, "\n")
      cat("True tree F4 diff range:", f4diff, "\n")
    }
    cat("=====================================\n")

    true_tree_result <- data.frame(
      run_id = "true_tree",
      rank = 0,
      score = score,
      worst_residual = true_tree_wr,
      WR_smaller_than_3 = true_tree_WR_smaller_than_3,
      hash = if (!is.null(true_tree_hash)) true_tree_hash
             else graph_hash(true_tree_graph),
      true_tree_or_not = TRUE,
      gamma1 = NA, # no admixture in the true tree
      gamma2 = NA,
      num_blocks = num_blocks,
      stringsAsFactors = FALSE
    )
    
    cat("True tree worst residual:", true_tree_wr, "\n")
  }

  results_list <- vector("list", length(unique_graphs))

  if (length(unique_graphs) == 0) { # empty table with the right columns
    results_df <- data.frame(
      run_id = character(0),
      rank = integer(0),
      score = numeric(0),
      worst_residual = numeric(0),
      WR_smaller_than_3 = logical(0),
      hash = character(0),
      true_tree_or_not = logical(0),
      gamma1 = numeric(0),
      gamma2 = numeric(0),
      num_blocks = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    for (i in seq_along(unique_graphs)) {

      graph_info <- unique_graphs[[i]]
      graph <- graph_info$graph
      score <- graph_info$score
      rank <- graph_info$rank
      run_id <- graph_info$run_id
      hash_val <- if ("hash" %in% names(graph_info)) graph_info$hash
                 else graph_hash(graph)

      result <- qpgraph(f2, graph, return_fstats = TRUE)
      wr <- result$worst_residual # returns abs(wr)
      graph_info$wr <- wr # store this value as well
      unique_graphs[[i]] <- graph_info # update unique graphs

      cat("Graph", i, "- WR:", wr, "- Score:", score, "- Hash:", hash_val, "\n")
      if (i <= 3) { # more details for the first graphs

        cat("  Full qpgraph result structure:\n")
        cat("    Score:", result$score, "\n")
        cat("    Worst residual:", result$worst_residual, "\n")
        cat("    F-stats available:", !is.null(result$f2), "\n")
        if (!is.null(result$f2)) {
          cat("    F2 fitted range:", range(result$f2$fit, na.rm = TRUE), "\n")
        }
        if (!is.null(result$f4)) {
          cat("    F4 fitted range:", range(result$f4$fit, na.rm = TRUE), "\n")
        }
      }

      WR_smaller_than_3 <- wr <= 3

      is_true_tree <- if (!is.null(true_tree_hash)) {
        hash_val == true_tree_hash
      } else {
        FALSE
      }

      # gammas, only when num_admix = 1 and the edges are there
      gamma1 <- NA
      gamma2 <- NA
      if (num_admix == 1 && "edges" %in% names(graph_info)) {
        edges_df <- graph_info$edges
        if (!is.null(edges_df) && "type" %in% names(edges_df)) {
          admix_edges <- edges_df[edges_df$type == "admix", ]
          if (nrow(admix_edges) >= 1) {
            # gammas are in the high and low columns, without the NA
            if ("high" %in% names(admix_edges) &&
                "low" %in% names(admix_edges)) {
              high_values <- admix_edges$high[!is.na(admix_edges$high)]
              low_values <- admix_edges$low[!is.na(admix_edges$low)]
              all_gamma_values <- c(high_values, low_values)
              unique_gamma_values <- unique(all_gamma_values)
              gamma1 <- max(unique_gamma_values) # gamma1 = the major one
              gamma2 <- min(unique_gamma_values)
              # they should sum to 1
              if (abs(gamma1 + gamma2 - 1.0) > 1e-6) {
                stop(paste("ERROR: Gamma values do not sum to 1.0!",
                          "gamma1 =", gamma1, 
                          "gamma2 =", gamma2,
                          "sum =", gamma1 + gamma2,
                          "for graph with hash:", hash_val))
              }
            }
          }
        }
      }

      results_list[[i]] <- data.frame(
        run_id = run_id, 
        rank = rank,
        score = score,
        worst_residual = wr,
        WR_smaller_than_3 = WR_smaller_than_3,
        hash = hash_val,
        true_tree_or_not = is_true_tree,
        gamma1 = gamma1,
        gamma2 = gamma2,
        num_blocks = num_blocks,
        stringsAsFactors = FALSE
      )
    }

  results_df <- do.call(rbind, results_list)
  }

  # true tree first, same column order
  if (!is.null(true_tree_result)) {
    if (nrow(results_df) > 0) {
      results_df <- results_df[, names(true_tree_result)]
    }
    results_df <- rbind(true_tree_result, results_df)
  }

  # rank over all runs by score, true tree not counted
  if (nrow(results_df) > 0) {
    if (!is.null(true_tree_result) && nrow(results_df) > 1) {
      true_tree_row <- results_df[1, ]
      other_results <- results_df[-1, ]
      other_results <- other_results[order(other_results$score), ]
      other_results$rank <- 1:nrow(other_results)
      results_df <- rbind(true_tree_row, other_results)
    } else if (is.null(true_tree_result)) {
      results_df <- results_df[order(results_df$score), ]
      results_df$rank <- 1:nrow(results_df)
    }
  }

  cat("Writing summary table to:", output_file_path, "\n")
  write.table(results_df, file = output_file_path, row.names = FALSE)
  return(unique_graphs)
}

# run_qpgraph with the true tree when num_admix = 0, without otherwise
run_qpgraph_for_graphs <- function(unique_graphs, method_suffix = "") {
  output_file <- if (method_suffix == "") {
    output_summary_table
  } else {
    get_summary_table_filename(output_dir, 
                              rep_id, 
                              num_admix, 
                              method_suffix, 
                              output_summary_table_suffix)
  }
  
  if (num_admix == 0) {
    return(run_qpgraph(
      unique_graphs = unique_graphs,
      f2 = f2,
      output_file_path = output_file,
      num_admix = num_admix,
      true_tree_graph = true_tree_graph,
      true_tree_hash = true_tree_hash,
      num_blocks = num_blocks
    ))
  } else {
    return(run_qpgraph(
      unique_graphs = unique_graphs,
      f2 = f2,
      output_file_path = output_file,
      num_admix = num_admix,
      true_tree_graph = NULL,
      true_tree_hash = NULL,
      num_blocks = num_blocks
    ))
  }
}

#-------------- Two-stage graph filtering --------------#
# 1. remove the duplicates (same hash)  2. second filter on the unique graphs
if (selection_method == "both") {
  combined_graph_lst <- graph_lst$main
  unique_graphs_all <- remove_duplicate_graphs(combined_graph_lst, rootfolder)
  unique_graphs_method1 <- apply_second_filter(
    unique_graphs_all, "keep_five_lowest_ll_graphs")
  unique_graphs_method2 <- apply_second_filter(
    unique_graphs_all, "select_threshold_within_ll", threshold)

  unique_graph_with_residual_method1 <- run_qpgraph_for_graphs(
    unique_graphs_method1, "keep_five_lowest_ll_graphs")
  unique_graph_with_residual_method2 <- run_qpgraph_for_graphs(
    unique_graphs_method2, "select_threshold_within_ll")
  
  # method 1 is the one used for the rest
  unique_graphs <- unique_graphs_method1
  unique_graph_with_residual <- unique_graph_with_residual_method1
  
} else {
  unique_graphs_all <- remove_duplicate_graphs(graph_lst, rootfolder)
  unique_graphs <- apply_second_filter(unique_graphs_all,
    selection_method, threshold)
  unique_graph_with_residual <- run_qpgraph_for_graphs(unique_graphs)
}

#-------------- Check if true tree was found when h = 0 --------------#
if (num_admix == 0) {
  cat("\n=== Checking if true species tree was found ===\n")
  cat("True species tree hash:", true_tree_hash, "\n")
  
  discovered_hashes <- sapply(unique_graphs, function(graph_info) {
    if ("hash" %in% names(graph_info)) {
      return(graph_info$hash)
    } else {
      return(graph_hash(graph_info$graph))
    }
  })
  
  cat("Number of unique graphs discovered:", length(discovered_hashes), "\n")
  cat("Discovered graph hashes:\n")
  for (i in seq_along(discovered_hashes)) {
    cat("  Graph", i, "hash:", discovered_hashes[i], "\n")
  }
  
  true_tree_found <- true_tree_hash %in% discovered_hashes
  
  if (true_tree_found) {
    matching_indices <- which(discovered_hashes == true_tree_hash)
    cat("\n*** Oh yay! True species tree was found! ***\n")
    cat("True tree hash", true_tree_hash, "matches graph(s):",
        paste(matching_indices, collapse = ", "), "\n")
    
    for (idx in matching_indices) {
      graph_info <- unique_graphs[[idx]]
      cat("  Matching graph", idx, "details:\n")
      cat("    Run ID:", graph_info$run_id, "\n")
      cat("    Rank:", graph_info$rank, "\n")
      cat("    Score:", graph_info$score, "\n")
      if ("wr" %in% names(graph_info)) {
        cat("    Worst residual:", graph_info$wr, "\n")
      }
    }
  } else {
    cat("\n*** Oh no! True species tree was NOT found ***\n")
    cat("True tree hash", true_tree_hash,
        "does not match any discovered graphs\n")
  }
}

#--- Save deduplicated graphs & f2 statistics for this rep ---#
saveRDS(unique_graph_with_residual, file = output_graph_file)

#-------------- Summary output for true tree finding analysis --------------#
cat("\n=== SUMMARY FOR REPLICATE", rep_id, "===\n")
cat("Number of admixture events (h):", num_admix, "\n")
cat("Number of runs:", runs, "\n")
cat("Selection method:", selection_method, "\n")
if (selection_method == "select_threshold_within_ll") {
  cat("Threshold:", threshold, "\n")
}
cat("True species tree hash:", true_tree_hash, "\n")
cat("Number of unique graphs found:", length(unique_graphs), "\n")

if (num_admix == 0) {
  discovered_hashes <- sapply(unique_graphs, function(graph_info) {
    if ("hash" %in% names(graph_info)) {
      return(graph_info$hash)
    } else {
      return(graph_hash(graph_info$graph))
    }
  })
  true_tree_found <- true_tree_hash %in% discovered_hashes
  cat("True species tree found:", ifelse(true_tree_found, "YES", "NO"), "\n")
}
cat("=====================================\n")
