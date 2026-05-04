#!/usr/bin/env julia

using CSV
using DataFrames
using Statistics
using ArgParse

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir"
            help = "Directory containing findgraph summary files"
            required = true
        "--output_file"
            help = "Output CSV file path"
            default = "findgraph_summary_stats.csv"
    end
    return parse_args(s)
end

"""
Extract parameter name root from filename
Example: "findgraph-DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5.csv" -> "DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5"
"""
function extract_paramname(filename::String)
    # Remove "findgraph-" prefix and ".csv" suffix
    paramname = replace(filename, r"^findgraph-" => "")
    paramname = replace(paramname, r"\.csv$" => "")
    return paramname
end

"""
Process a single findgraph file and return summary statistics
"""
function process_file(filepath::String, paramname::String)
    # Read the CSV file
    df = CSV.read(filepath, DataFrame)
    
    # Total number of replicates
    total_reps = nrow(df)
    
    # Count best_k categories (handle both string and numeric values)
    h0_accepted = count(x -> (x == 0 || x == "0"), df.best_k)
    h1_accepted = count(x -> (x == 1 || x == "1"), df.best_k)
    bt1_accepted = count(x -> (x == ">1"), df.best_k)
    
    # Percentage of true_tree_found_in_H0 = True
    true_tree_h0_count = count(df.true_tree_found_in_H0 .== true)
    pct_true_tree_h0 = true_tree_h0_count / total_reps * 100
    
    # Average num_blocks (handling missing values)
    avg_num_blocks = mean(skipmissing(df.num_blocks))
    
    # Average gamma values
    avg_gamma1_h1 = mean(skipmissing(df.gamma1_H1))
    avg_gamma2_h1 = mean(skipmissing(df.gamma2_H1))
    
    # Percentage of true tree displayed in H1 (either major or minor)
    true_tree_h1_count = count((df.True_tree_displayed_H1_major .== true) .| 
                               (df.True_tree_displayed_H1_minor .== true))
    pct_true_tree_h1 = true_tree_h1_count / total_reps * 100
    
    return (
        paramname_root = paramname,
        total_replicates = total_reps,
        H0Accepted = h0_accepted,
        H1Accepted = h1_accepted,
        BT1Accepted = bt1_accepted,
        pct_true_tree_H0 = pct_true_tree_h0,
        avg_num_blocks = avg_num_blocks,
        avg_gamma1_H1 = avg_gamma1_h1,
        avg_gamma2_H1 = avg_gamma2_h1,
        pct_true_tree_H1 = pct_true_tree_h1
    )
end

"""
Main function
"""
function main()
    args = parse_commandline()
    
    input_dir = args["input_dir"]
    output_file = args["output_file"]
    
    # Check if input directory exists
    if !isdir(input_dir)
        error("Input directory does not exist: $input_dir")
    end
    
    # Find all files starting with "findgraph"
    all_files = readdir(input_dir)
    findgraph_files = filter(f -> startswith(f, "findgraph") && endswith(f, ".csv"), all_files)
    
    if isempty(findgraph_files)
        error("No findgraph files found in $input_dir")
    end
    
    println("Found $(length(findgraph_files)) findgraph files")
    
    # Process each file
    results = []
    for filename in findgraph_files
        println("Processing: $filename")
        filepath = joinpath(input_dir, filename)
        paramname = extract_paramname(filename)
        
        try
            stats = process_file(filepath, paramname)
            push!(results, stats)
        catch e
            @warn "Error processing $filename: $e"
        end
    end
    
    # Convert to DataFrame
    results_df = DataFrame(results)
    
    # Sort by paramname_root for consistent output
    sort!(results_df, :paramname_root)
    
    # Write output
    CSV.write(output_file, results_df)
    println("\nSummary written to: $output_file")
    println("Processed $(nrow(results_df)) files successfully")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
