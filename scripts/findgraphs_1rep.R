# Script to run findgraphs for one replicate

library(admixtools)
library(optparse)
library(dplyr)
library(tidyverse)

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
              help = "Block sizes used in findgraphs"),
  make_option(c("-r", "--simulation_rep"), type = "integer",
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
                help = "Suffix to the table outputs")
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
simulation_rep <- opt$simulation_rep
blgsize <- opt$blgsize
seed_file_path <- opt$seed_file_path
output_graph_suffix <- opt$output_graph_suffix
output_f2_suffix <- opt$output_f2_suffix
output_summary_table_suffix <- opt$output_summary_table_suffix

input_prefix <- file.path(input_dir, prefix)
output_graph_file <- file.path(output_dir, paste0("rep", simulation_rep, "_admix", num_admix, output_graph_suffix)) # nolint
output_graph_f2 <- file.path(output_dir, paste0("rep", simulation_rep, "_admix", num_admix, output_f2_suffix))  # nolint
output_summary_table <- file.path(output_dir, paste0("rep", simulation_rep, "_admix", num_admix, output_summary_table_suffix)) # nolint 

#------------ set up seeds --------------#
seed_array_two_columns <- read.table(seed_file_path, header = TRUE)

if (num_admix == 0) { # if num_admix = 0, then use the first column
  seed_array <- seed_array_two_columns[, 1]
} else { # if num_admix = 1, then use the second column
  seed_array <- seed_array_two_columns[, 2]
}

# # Determine block length from last line in .snp file
# last_line_in_snp <- tail(read.csv(paste0(prefix,".snp"),
#                           sep = "\t",header = FALSE),n = 1)
# blgsize <- as.integer(strsplit(last_line_in_snp[[1]],
#                                   split = " +")[[1]][4]) + 1)

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
#'         "$scores" (the log-likelihood score).
#'
#' @examples
#' result <- run_find_graphs_replicate(f2,num_admix = 1,stop_gen = 100,
#'                                     initgraph = g,outgroup = "A",
#'                                     runs = 100,rep = 1)
#' result[["run1"]][[1]]$graphs  # Access the first graph from run 1
#' result[["run1"]][[1]]$scores  # Get its score
run_find_graphs_replicate <- function(f2, num_admix, stop_gen, outgroup, runs, simulation_rep, seed_array) { # nolint

  graph_lst <- list()

  for (i in 1:runs) {

    set.seed(seed_array[simulation_rep]) # set seed for specific

    message(paste0("Running find_graphs for rep ", simulation_rep, " at run ", i)) # nolint

    opt_results <- find_graphs(f2, numadmix = num_admix,
                              stop_gen = stop_gen, # nolint
                              outpop = outgroup)
    # Select top 5 graphs with lowest LL score
    Top5_graphs_per_run <- opt_results %>% dplyr::slice_min(score, n = 5, with_ties = TRUE) # nolint
    # create a list for top 5 graphs per run -> examples:
    # To extract the 1st graph top5_graphs_list_per_run[[1]]$graphs
    # To extract the 1st graph's score top5_graphs_list_per_run[[1]]$scores
    top5_graphs_list_per_run <- lapply(1:nrow(Top5_graphs_per_run), function(j) { # nolint
          list(graphs = Top5_graphs_per_run$graph[[j]], # nolint
               scores = Top5_graphs_per_run$score[j],
               rank = j)}) # rank: rank of the graph based on score in each run

    graph_lst[[paste0("run", i)]] <- top5_graphs_list_per_run
  cat("Finishing findgraphs for rep", simulation_rep, "with", runs, "runs.\n") # nolint
  }
  return(graph_lst)
}

graph_lst <- run_find_graphs_replicate(f2 = f2,
                                      num_admix = num_admix,
                                      stop_gen = stop_gen,
                                      outgroup = outgroup,
                                      runs = runs,
                                      simulation_rep = simulation_rep,
                                      seed_array = seed_array)

#-------------- Remove duplicated graphs --------------#
#' Remove topologically redundant graphs from a list of runs
#'
#' @param graph_lst A list of runs,each containing a list of top-scoring graphs
#'        with components "$graphs" and "$scores".
#'
#' @return List of unique graphs with "$graph","$score","$run", and "$topo_key"
#'
#' @examples
#' dedup_graphs <- deduplicate_graphs(graph_lst)
#' dedup_graphs[[1]]$graph  # igraph object
#' dedup_graphs[[1]]$score  # score
#' dedup_graphs[[1]]$run    # run number
deduplicate_graphs <- function(graph_lst) {
  all_graphs <- list()
  index <- 1

  for (run_name in names(graph_lst)) { # get repID-runID
    run_graphs <- graph_lst[[run_name]]
    run_number <- as.integer(gsub("run", "", run_name))

    for (g in run_graphs) {
      # graph_lst[[runs]] --> Each run saves a list of 5(+) graph object
      # all_graph --> all_graphs[[i]] stores a list of the graph object
      all_graphs[[index]] <- list(
        graph = g$graphs,
        score = g$scores,
        rank = g$rank,
        run = run_number
      )
      index <- index + 1 # number of graphs across all runs in this rep
    }
  }

  get_topology_key <- function(graph) {
    # create an unique string for each graph topology
    edges <- igraph::as_edgelist(graph, names = TRUE)
    edge_str <- apply(edges, 1, function(x) paste(x[1], x[2], sep = "->"))
    paste(sort(edge_str), collapse = ";") # sort edges so we can compare
  }

  topo_keys <- sapply(all_graphs, function(x) get_topology_key(x$graph))
  unique_indices <- !duplicated(topo_keys)
  unique_graphs <- all_graphs[unique_indices]

  for (i in seq_along(unique_graphs)) {
    # store topology key as metadata
    unique_graphs[[i]]$topo_key <- topo_keys[unique_indices][i]
  }

  return(unique_graphs)
}

#' Organize Worst Residual LL-scores Outputs
#'
#' This processes a list of graphs, calculates the worst residuals for each,
#' and outputs the results to a specified file.
#'
#' @param graph_lst A list of graph objects to be processed.
#' @param f2 The f2 statistics used in the analysis.
#' @param output_file_path The file path where the output will be saved.
#'
#' @return A list of unique graphs after processing.
#'        wirte the summary table into a summary table
#'
#' @examples
#' \dontrun{
#' unique_graphs <- organize_wr_outputs(graph_lst, f2, "output_summary.txt")
#' }
organize_wr_outputs <- function(graph_lst, f2, output_file_path){
  unique_graphs <- deduplicate_graphs(graph_lst)
  results_list <- vector("list", length(unique_graphs))

  for (i in seq_along(unique_graphs)) {
    graph_info <- unique_graphs[[i]]
    graph <- graph_info$graph
    ll_score <- graph_info$score
    rank <- graph_info$rank
    run_id <- graph_info$run

    result <- qpgraph(f2, graph, return_fstats = TRUE)
    wr <- abs(result$worst_residual) 
    # This value doesn't have the sign, so use result$f4 %>% slice_max(abs(z), with_ties = F) to get the sign 

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

unique_graphs <- organize_wr_outputs(graph_lst = graph_lst,
                    f2 = f2,
                    output_file_path = output_summary_table)

#--- Save deduplicated graphs & f2 statistics for this rep ---#
saveRDS(unique_graphs, file = output_graph_file)
saveRDS(f2, file = output_graph_f2)