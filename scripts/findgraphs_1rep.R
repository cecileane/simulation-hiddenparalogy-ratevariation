# R script to run each individual replicate with n runs 

library(admixtools)
library(optparse)
library(dplyr)
library(tidyverse)
library(Rcpp)

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Path to the folder containing Eigenstrat input files"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Path to output folder for results"),
  make_option(c("-p", "--prefix"), type = "character", default = "",
              help = "Prefix to find the eigenstrat files"),
  make_option(c("-k", "--num_admix"), type = "integer", default = 0,
              help = "Number of admixture events"),
  make_option("--stop_gen", type = "integer", default = 100,
              help = "Number of generations to stop find_graphs"),
  make_option("--outgroup", type = "character", default = "A",
              help = "Outgroup (Default Homo \"A\")"),
  make_option("--runs", type = "integer", default = 100,
              help = "Number of times to run find_graphs"),
  make_option(c("-b", "--blgsize"), type = "integer", default = 300,
              help = "Block sizes used in find_graphs"),
  make_option(c("-r", "--rep_id"), type = "character",
              help = "Replication ID for printing and result saving"),
  make_option(c("-s", "--seed_file_path", type = "string",
                help = "A path to the seed file to be used")),
  make_option("--output_graph_suffix", type = "character",
              default = "_graphs.rds",
              help = "Suffix to the output graphs files"),
  make_option("--output_f2_suffix", type = "character",
              default = "_f2.rds",
              help = "Suffix to the output graphs files"),
  make_option("--output_summary_table_suffix", type = "character",
              default = "_summary_table.txt",
              help = "Suffix to the table outputs"),
  make_option("--rootfolder", type = "character",
              default = "__use_cwd__",
              help = "rootfolder to save any warnings")
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
seed_file_path <- opt$seed_file_path
output_graph_suffix <- opt$output_graph_suffix
output_f2_suffix <- opt$output_f2_suffix
output_summary_table_suffix <- opt$output_summary_table_suffix
rootfolder <- opt$rootfolder
if (rootfolder == "__use_cwd__") {
  rootfolder <- getwd()
}

input_prefix <- file.path(input_dir, prefix)
output_graph_file <- file.path(output_dir, paste0("rep", rep_id, "_admix", num_admix, output_graph_suffix)) # nolint
output_graph_f2 <- file.path(output_dir, paste0("rep", rep_id, "_admix", num_admix, output_f2_suffix))  # nolint
output_summary_table <- file.path(output_dir, paste0("rep", rep_id, "_admix", num_admix, output_summary_table_suffix)) # nolint

#------------ set up seeds --------------#
seed_array_two_columns <- read.table(seed_file_path, header = TRUE)

if (num_admix == 0) { # if num_admix = 0, then use the first column
  seed_array <- seed_array_two_columns[, 1]
} else { # if num_admix = 1, then use the second column
  seed_array <- seed_array_two_columns[, 2]
}

#-------------- Create f2 statistics --------------#
f2 <- f2_from_geno(pref = input_prefix,
                   blgsize = blgsize,
                   outpop = outgroup)

#-------------- Run find_graphs --------------#
#' Run find_graphs repeatedly and extract top-scoring graphs per run
#'
#' This function performs multiple runs of "find_graphs()" on a given f2 matrix,
#' extracts the top N graphs from each run, and returns a structured list
#' storing the graphs and their scores by run.
#'
#' @param f2 An f2 object,typically from "f2_from_geno()".
#' @param num_admix Integer. Number of admixture events (K).
#' @param stop_gen Integer. Number of generations to stop graph search.
#' @param outgroup Character. Name of the outgroup population.
#' @param runs Integer. Number of independent "find_graphs()" runs to perform.
#' @param rep Integer or string. Replicate ID (used for printing/logging).
#'
#' @return A named list ("graph_lst") where each element corresponds to a run
#'         and contains a list of top-scoring graphs. Each graph is represented
#'         as a list with two elements: "$graphs" (the igraph object),and
#'         "$score" (the log-likelihood score).
#'
#' @examples
#' result <- run_find_graphs_replicate(f2,num_admix = 1,stop_gen = 100,
#'                                     initgraph = g,outgroup = "A",
#'                                     runs = 100,rep = 1)
#' result[["run1"]][[1]]$graph  # Access the first graph from run 1
#' result[["run1"]][[1]]$score  # Get its score
run_find_graphs_replicate <- function(f2, num_admix, stop_gen, outgroup, runs, rep_id, seed_array) { # nolint

  graph_lst <- list()

  for (i in 1:runs) {

    # rep_id is used to specify the output file name (string)
    # rep is used to index seed array (int)
    rep <- as.integer(rep_id)
    set.seed(seed_array[rep]) # set seed for specific

    message(paste0("Running find_graphs for rep ", rep_id, " at run ", i)) # nolint

    opt_results <- find_graphs(f2, numadmix = num_admix,
                              stop_gen = stop_gen, # nolint
                              outpop = outgroup)

    # Select top 5 graphs with lowest LL score
    Top5_graphs_per_run <- opt_results %>% dplyr::slice_min(score, n = 5, with_ties = TRUE) # nolint
    # create a list for top 5 graphs per run -> examples:
    # To extract the 1st graph top5_graphs_list_per_run[[1]]$graphs
    # To extract the 1st graph's score top5_graphs_list_per_run[[1]]$scores
    top5_graphs_list_per_run <- lapply(1:nrow(Top5_graphs_per_run), function(j) { # nolint
          list(graph = Top5_graphs_per_run$graph[[j]], # nolint
               score = Top5_graphs_per_run$score[j],
               rank = j, # rank: rank of the graph based on score in each run
               run_id = i)})

    graph_lst[[paste0("run", i)]] <- top5_graphs_list_per_run
  cat("Finishing findgraphs for rep", rep_id, "with", runs, "runs.\n") # nolint
  }
  return(graph_lst)
}

graph_lst <- run_find_graphs_replicate(f2 = f2,
                                      num_admix = num_admix,
                                      stop_gen = stop_gen,
                                      outgroup = outgroup,
                                      runs = runs,
                                      rep_id = rep_id,
                                      seed_array = seed_array)

#-------------- Remove duplicated graphs --------------#
#' Add hash values to each graph in a list of graph search results
#'
#' This function loops over a nested list of results from find_graphs,
#' It finds a hash for each graph using `graph_hash()`, and appends it
#' as a $hash to each graph entry.
#' Graphs with the same hash are believed to be topologically similar
#'
#' @param graph_lst A nested list: each element corresponds to a run
#' (e.g., "run1") and contains a list of graph results,
#' each with `$graph`, `$score`, and `$rank`.
#'
#' @return The same `graph_lst` structure,
#' but with an added `$hash` field for each graph.

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

#' Remove duplicate graphs based on hash and likelihoood scores
#'
#' This removes graphs with duplicate topology (same has from graph_hash)
#' across runs, while checking for consistency in likelihood scores.
#' If two graphs have the same hash (they are believed to be topologically
#' identitcal) but scores that differ more than a tolerance, an error will
#' be raised.
#'
#' @param graph_lst A list of graphs (output from 'run_find_graphs_replicate'),
#' potentially with duplicates.
#' @param score_tol Numeric. Allowed difference in log-likelihood scores for
#' graphs with identical topology. Default = 1e-6.
#'
#' @return A new graph_lst with duplicates removed
#' @notes This function is called within ''
remove_duplicate_graphs <- function(graph_lst, rootfolder, score_tol = 1e-6) {
  
  graph_lst <- add_hashes_to_graph_lst(graph_lst)

  hash_seen <- list()
  unique_graph_lst <- list()

  # Save warnings to file
  warning_log_path <- file.path(rootfolder, "findgraph_warnings.txt")
  warning_file <- file(warning_log_path, open = "at")

  for (run_name in names(graph_lst)) {
    for (graph_info in graph_lst[[run_name]]) {
      h <- graph_info$hash
      score <- graph_info$score
      rank <- graph_info$rank

      if (!h %in% names(hash_seen)) {
        graph_info$run_id <- run_name
        hash_seen[[h]] <- list(score = score, graph_info = graph_info)
        unique_graph_lst[[length(unique_graph_lst) + 1]] <- graph_info
      } else {
        prev_score <- hash_seen[[h]]$score
        if (abs(prev_score - score) > score_tol) {
          warning_msg <- paste0(
            "WARNING: Graph from rep ", rep_id, "mrun ", run_name,
            " with rank ", rank,
            " has identical hash with different scores: ",
            prev_score, " vs ", score,
            " (difference = ", abs(prev_score - score), ")."
          )
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
#' This processes a list of graphs, calculates the worst residuals for each,
#' and outputs the results to a specified file.
#'
#' @param graph_lst A list of graph objects to be processed.
#' @param f2 The f2 statistics used in the analysis.
#' @param output_file_path The file path where the output will be saved.
#'
#' @return A list of unique graphs after processing.
#'        wirte the summary table into a summary table
run_qpgraph <- function(unique_graphs, f2, output_file_path){

  results_list <- vector("list", length(unique_graphs))

  for (i in seq_along(unique_graphs)) {

    graph_info <- unique_graphs[[i]]
    graph <- graph_info$graph
    ll_score <- graph_info$score
    rank <- graph_info$rank
    run_id <- graph_info$run_id

    # de-bugging:
    # cat("---- Processing graph", i, "----\n")
    # cat("run_id:", graph_info$run_id, "\n")
    # cat("rank:", graph_info$rank, "\n")
    # cat("hash:", graph_info$hash, "\n")
    # cat("score:", graph_info$score, "\n")
    # cat("graph class:", class(graph_info$graph), "\n\n")

    result <- qpgraph(f2, graph, return_fstats = TRUE)

    wr <- result$worst_residual # returns abs(wr)
    graph_info$wr <- wr # store this value as well
    unique_graphs[[i]] <- graph_info # update unique graphs

    WR_smaller_than_3 <- wr <= 3

    results_list[[i]] <- data.frame(
      run_id = run_id,
      rank = rank,
      ll_score = ll_score,
      worst_residual = wr,
      WR_smaller_than_3 = WR_smaller_than_3,
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)

  write.table(results_df, file = output_file_path, row.names = FALSE)
  
  return(unique_graphs)
}

unique_graphs <- remove_duplicate_graphs(graph_lst, rootfolder)
unique_graph_with_residual <- run_qpgraph(unique_graphs = unique_graphs,
                    f2 = f2, # nolint
                    output_file_path = output_summary_table)

#--- Save deduplicated graphs & f2 statistics for this rep ---#
saveRDS(unique_graphs, file = output_graph_file)
saveRDS(f2, file = output_graph_f2)