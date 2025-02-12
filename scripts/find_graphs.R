# R script to run through admixture tools 

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
  devtools::install_github("uqrmaie1/admixtools")
}

library(admixtools)
library(optparse) # parse arguments 

option_list <- list(
  make_option(c("--prefix"), type = "character", default = NULL, 
              help = "Prefix of the eigenstrat files", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Get f2 statistics 
f2_results <- f2_from_geno(
  pref = opt$prefix
)

output_file <- paste0(opt$prefix, "_f2_results.RData")
save(f2_results, file = output_file)
 



