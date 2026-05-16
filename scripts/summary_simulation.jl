#!/usr/bin/env julia
# ============================================================================
# scripts/summary_simulation.jl
#
# Purpose : Concatenate per-parameter simulation summary CSVs into one
#           cross-setting table covering gene-tree counts, hidden-paralogy
#           categories, RF distances, and branch-length statistics.
# Inputs  : simulation_summary/summary_<paramname>.csv     (one per setting)
# Outputs : results/summary_concatenated.csv               (cross-setting table)
# Usage   : julia --project=. scripts/summary_simulation.jl
#           (optional overrides: --input_dir, --output_file)
# Note    : Run after simulation_postprocess.jl has produced per-setting CSVs
#           and run_postprocessing.jl has staged them under simulation_summary/.
# ============================================================================

using CSV
using DataFrames
using ArgParse

"""
function to parse command line arguments
"""
function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input_dir"
            help = "Directory containing summary CSV files"
            default = "./simulation_summary"
        "--output_file"
            help = "Output file path for concatenated results"
            default = "./results/summary_concatenated.csv"
    end
    return ArgParse.parse_args(s)
end

"""
function to find all summary CSV files in the input directory
"""
function find_summary_files(input_dir)
    files = readdir(input_dir, join=true)
    summary_files = []
    
    println("Found $(length(files)) files in the input directory")
    
    for file in files
        if occursin(r"^summary_.*\.csv$", basename(file))
            push!(summary_files, file)
        end
    end
    
    println("Found $(length(summary_files)) summary CSV files to concatenate")
    return summary_files
end

"""
Main function to concatenate all summary CSV files
""" 
function main()

    mkpath("./results") 

    args = parse_arguments()
    
    input_dir = get(args, "input_dir", ".") 
    output_path = get(args, "output_file", "./summary_concatenated.csv") 
    
    summary_files = find_summary_files(input_dir)
    
    if isempty(summary_files)
        println("No summary CSV files found to concatenate.")
        return
    end
    
    # Read and concatenate all CSV files
    all_data = DataFrame()
    
    for (i, file) in enumerate(summary_files)
        println("Processing $(i)/$(length(summary_files)): $(basename(file))")
        try
            df = CSV.read(file, DataFrame)
            if i == 1
                all_data = df
            else
                all_data = vcat(all_data, df, cols=:union)
            end
        catch e
            println("Warning: Could not read file $(basename(file)): $e")
        end
    end
    
    if nrow(all_data) > 0
        CSV.write(output_path, all_data)
        println("Concatenated $(length(summary_files)) files → $output_path")
        println("Total rows in output: $(nrow(all_data))")
    else
        println("No data to write to output file.")
    end
end

main()