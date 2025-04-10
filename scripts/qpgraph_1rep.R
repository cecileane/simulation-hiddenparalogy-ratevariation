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
                help = "A path to the seed file to be used"))
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
input_prefix <- file.path(input_dir, prefix)
output_graph_file <- file.path(output_dir, paste0("rep", simulation_rep, "_unique_graphs.rds")) # nolint
output_graph_f2 <- file.path(output_dir, paste0("rep", simulation_rep, "_f2.rds"))  # nolint


# ----------- Read input files + Remove duplicated across --------------# 

for (file in rds_files) {
  graph_list <- readRDS(file)
  for (g in graph_list) {
    all_graphs[[index]] <- g
    index <- index + 1
  }
}