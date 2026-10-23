#!/bin/bash
# ============================================================================
# scripts/seq-gen.sh
#
# Purpose : Per-gene Seq-Gen wrapper. Simulates DNA sequences along a single
#           gene tree under HKY+Gamma with per-gene rate, kappa, and base
#           frequencies passed in by the caller. Skips silently if the input
#           gene-tree file is empty (removed during paralogy filtering).
# Inputs  : $2 = gene tree file (Newick)
# Outputs : $3 = Seq-Gen sequence file (relaxed-PHYLIP, -on)
# Usage   : Called from simulation.jl per gene per replicate. see simulation.jl
# ============================================================================

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
