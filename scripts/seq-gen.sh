#!/bin/bash

# In scripts/simulation.jl, $input_path and $output_path cannot be pre-specified. 
rep_number_string="$1"
input_file="$2"
output_file="$3" 
seed="$4"

# seqgen specific parameters
alpha="$5"  # shape parameter for Gamma distribution of rates across sites
kappa="$6"  # transition/transversion ratio
A="$7"      # base frequency of A
C="$8"      # base frequency of C
G="$9"      # base frequency of G
T="${10}"   # base frequency of T 

# The length of the gene sequence to be simulated  
gene_len="${11}" 

# Tree files not statisfying the requirements are removed by modify_newick from utilities. 
# This caused missing tree file, 
# so the below part is crucial to deal with missing tree file: 
if [ ! -s "$input_file" ]; then # Check if input file exists and is non-empty
    echo "Input file $input_file is missing or empty. Skipping seq-gen for this gene."
    exit 0  # Exit without error
fi

echo "Running seq-gen for input file: $input_file" 
echo "Sequence length: $gene_len, alpha: $alpha, kappa: $kappa"

# Run seq-gen with specified parameters 
executables/seq-gen \
    -l"$gene_len" \
    -mHKY \
    -a"$alpha" \
    -t"$kappa" \
    -f"$A","$C","$G","$T" \
    -z "$seed" \
    -on < "$input_file" > "$output_file"

# Run seq-gen only if the input file is not empty
# - to simulate all genes with the same substitution model, use:
#   * HKY (-m option) with transition/transversion ratio kappa = 4.143 (option -t)
#   * base frequencies 0.316,0.182,0.183,0.319 (-f option)
#   * shape alpha = 0.356 (-a option) for the Gamma distribution of rates across sites
# Under all same model, all genes will use the same parameters above. 
# executables/seq-gen \
#         -l"$seq_len" \
#         -mHKY \
#         -a0.356 \
#         -t4.143 \
#         -f0.316,0.182,0.183,0.319 \
#         -z "$seed" \
#         -on < "$input_file" > "$output_file" 
# Under all different model, each gene will have its own parameters sampled from empirical distributions. 
# to simulate each gene with its own substitution model, use HKY with:
#   * kappa from LogNormal(μ=1.4215, σ=0.2798)
#   * frequencies from Dirichlet(66.59, 38.41, 38.61, 67.12)
#   * alpha from Gamma(α=3.267, θ=0.109).
