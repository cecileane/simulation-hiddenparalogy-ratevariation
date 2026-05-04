#!/bin/bash

# simulation.jl cannot pre-specify $input_path/$output_path, so pass as args
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

# modify_newick removes invalid tree files; guard for missing input:
if [ ! -s "$input_file" ]; then # Check if input file exists and is non-empty
    echo "Input $input_file is missing or empty. Skipping seq-gen."
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
#   * HKY (-m option) with kappa = 4.143 (-t option)
#   * base frequencies 0.316,0.182,0.183,0.319 (-f option)
#   * Gamma shape alpha = 0.356 (-a option)
# Under all same model, all genes will use the same parameters above. 
# executables/seq-gen \
#         -l"$seq_len" \
#         -mHKY \
#         -a0.356 \
#         -t4.143 \
#         -f0.316,0.182,0.183,0.319 \
#         -z "$seed" \
#         -on < "$input_file" > "$output_file" 
# Per-gene model: parameters sampled from empirical distributions.
# to simulate each gene with its own substitution model, use HKY with:
#   * kappa from LogNormal(μ=1.4215, σ=0.2798)
#   * frequencies from Dirichlet(66.59, 38.41, 38.61, 67.12)
#   * alpha from Gamma(α=3.267, θ=0.109).
