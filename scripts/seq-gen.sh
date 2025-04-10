#!/bin/bash

# In scripts/simulation.jl, $input_path and $output_path cannot be pre-specified. 
rep_number_string="$1"
input_file="$2"
output_file="$3" 
seed="$4"

# Tree files not statisfying the requirements are removed by modify_newick from utilities. 
# This caused missing tree file, and the below part is crucial to deal with missing tree file: 
if [ ! -s "$input_file" ]; then # Check if input file exists and is non-empty
    echo "Input file $input_file is missing or empty. Skipping seq-gen for this gene."
    exit 0  # Exit without error
fi

# Run seq-gen only if the input file is not empty
executables/seq-gen  -l1000 -mHKY -a0.356 -t4.143 -f0.316,0.182,0.183,0.319 -z "$seed" -on < "$input_file" > "$output_file" 
 