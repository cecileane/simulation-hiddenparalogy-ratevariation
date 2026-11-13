#!/bin/bash
# Seq-Gen for one gene: simulates sequences along one gene tree under
# HKY + Gamma, with the rate, kappa and base frequencies given as arguments.
# Does nothing if the gene tree file is empty (gene dropped by the
# paralogy filter). Called by simulation.jl for each gene of each replicate

# simulation.jl cannot pre-specify $input_path/$output_path, so pass as args
rep_number_string="$1"
input_file="$2"
output_file="$3" 
seed="$4"

# seqgen specific parameters
alpha="$5"  # shape parameter for Gamma distribution of rates across sites
tstv="$6"   # expected transition/transversion ratio (seq-gen -t)
A="$7"      # base frequency of A
C="$8"      # base frequency of C
G="$9"      # base frequency of G
T="${10}"   # base frequency of T 

# The length of the gene sequence to be simulated  
gene_len="${11}" 
tree="${12}"

# modify_newick removes invalid tree files; guard for missing input:
if [ ! -s "$input_file" ]; then # Check if input file exists and is non-empty
    echo "Input $input_file is missing or empty. Skipping seq-gen."
    exit 0  # Exit without error
fi

case "$tree" in plant*) treelabel=plant;; *) treelabel=reptile;; esac
# print if this is a reptile or plant tree
echo "Simulating molecular sequences for a $treelabel tree" 
echo "Running seq-gen for input file: $input_file" 
echo "Sequence length: $gene_len, alpha: $alpha, ts/tv: $tstv"

# Run seq-gen with specified parameters 
executables/seq-gen \
    -l"$gene_len" \
    -mHKY \
    -a"$alpha" \
    -t"$tstv" \
    -f"$A","$C","$G","$T" \
    -z "$seed" \
    -on < "$input_file" > "$output_file"
