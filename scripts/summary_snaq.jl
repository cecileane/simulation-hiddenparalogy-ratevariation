#!/usr/bin/env julia

using CSV
using DataFrames
using Statistics
using ArgParse

include("visual_utilities.jl")

using Plots

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir"
            help = "Directory containing findgraph summary files"
            default = "snaq_summary" 
        "--output_file"
            help = "Output CSV file path"
            default = "results/SNaQ_summary.csv"
        "--visualization_output_dir" 
            help = "Output directory for visualization files"
            default = "visualization_results/snaq" 
    end
    return parse_args(s)
end


"""
Script to summarize goodness-of-fit results from SNaQ analysis.

This script processes all CSV files in the snaq_summary folder and creates
a summary table with hypothesis testing results and mean statistics.
"""

function extract_parameter_setting(filename::String)
    """Extract parameter setting from filename by removing the suffix."""
    # Remove the -_summary.csv suffix
    setting = replace(filename, r"-_summary\.csv$" => "")
    return setting
end

function classify_hypothesis(p_H0::Float64, p_H1::Float64)
    """
    Classify which hypothesis is supported based on p-values.
    
    Returns:
    - H0: true if p_H0 > 0.05
    - H1: true if p_H0 <= 0.05 and p_H1 > 0.05  
    - H_gt_1: true if both p_H0 <= 0.05 and p_H1 <= 0.05
    """
    if p_H0 > 0.05
        return (true, false, false)  # H0, H1, H>1
    elseif p_H1 > 0.05
        return (false, true, false)  # H0, H1, H>1
    else
        return (false, false, true)  # H0, H1, H>1
    end
end

function process_csv_file(filepath::String)
    """Process a single CSV file and return summary statistics."""
    try
        # Read the CSV file
        df = CSV.read(filepath, DataFrame)
        
        # Get parameter setting from filename
        filename = basename(filepath)
        parameter_setting = extract_parameter_setting(filename)
        
        # Count how many replicates support each hypothesis
        h0_count = 0
        h1_count = 0
        h_gt_1_count = 0
        
        for row in eachrow(df)
            h0_support, h1_support, h_gt_1_support = classify_hypothesis(row.p_H0, row.p_H1)
            if h0_support
                h0_count += 1
            elseif h1_support
                h1_count += 1
            elseif h_gt_1_support
                h_gt_1_count += 1
            end
        end
        
        # Calculate mean scores
        mean_score_H0 = mean(df.score_H0)
        mean_score_H1 = mean(df.score_H1)
        
        # Calculate mean gamma values
        mean_gamma_1 = mean(df.gamma_1)
        mean_gamma_2 = mean(df.gamma_2)
        
        # Calculate mean p-values for reference
        mean_p_H0 = mean(df.p_H0)
        mean_p_H1 = mean(df.p_H1)
        
        # Count how many replicates found the true network topology
        # find_true_net0: count 0.0 values in RF_net0_true
        find_true_net0 = sum(df.RF_net0_true .== 0.0)
        find_true_net0_noF = sum(df.RF_net0_true_noF .== 0.0)
        
        # find_true_net1: count replicates where EITHER major or minor tree matched the true network
        # Use logical OR to avoid double-counting replicates where both trees match
        find_true_net1 = sum((df.RF_net1_1_true .== 0.0) .| (df.RF_net1_2_true .== 0.0))
        find_true_net1_noF = sum((df.RF_net1_1_true_noF .== 0.0) .| (df.RF_net1_2_true_noF .== 0.0))

        # Check if net0's RF distance to alter1, alter2, alter3 together
        # Only count rows where either RF_net1_1_true != 0 AND RF_net1_2_true != 0
        # Skip rows where both RF_net1_1_true == 0 OR RF_net1_2_true == 0
        valid_rows_net0 = df.RF_net1_1_true .!= 0.0  
        valid_rows_net1 = (df.RF_net1_1_true .!= 0.0) .&& (df.RF_net1_2_true .!= 0.0)
        
        find_alter1_net0 = sum((df.RF_net0_alter1 .== 0.0) .& valid_rows_net0)
        find_alter2_net0 = sum((df.RF_net0_alter2 .== 0.0) .& valid_rows_net0)
        find_alter3_net0 = sum((df.RF_net0_alter3 .== 0.0) .& valid_rows_net0)
        find_alter_net0 = find_alter1_net0 + find_alter2_net0 + find_alter3_net0

        # check if net1's major hybrid's RF distance to alter1, alter2, alter3
        find_alter1_net1_major = sum((df.RF_net1_1_alter1 .== 0.0) .& valid_rows_net1)
        find_alter2_net1_major = sum((df.RF_net1_1_alter2 .== 0.0) .& valid_rows_net1)
        find_alter3_net1_major = sum((df.RF_net1_1_alter3 .== 0.0) .& valid_rows_net1)
        
        # check if net1's minor hybrid's RF distance to alter1, alter2, alter3
        find_alter1_net1_minor = sum((df.RF_net1_2_alter1 .== 0.0) .& valid_rows_net1)
        find_alter2_net1_minor = sum((df.RF_net1_2_alter2 .== 0.0) .& valid_rows_net1)
        find_alter3_net1_minor = sum((df.RF_net1_2_alter3 .== 0.0) .& valid_rows_net1)
        
        # Return results as a named tuple
        return (
            parameter_setting = parameter_setting,
            H_eq_0 = h0_count,
            H_eq_1 = h1_count, 
            H_gt_1 = h_gt_1_count,
            mean_score_H0 = mean_score_H0,
            mean_score_H1 = mean_score_H1,
            mean_gamma_1 = mean_gamma_1,
            mean_gamma_2 = mean_gamma_2,
            # whether net0 and net1 matched true network: 
            find_true_net0 = find_true_net0,
            find_true_net0_noF = find_true_net0_noF, 
            find_true_net1 = find_true_net1, # include both major and minor tree 
            find_true_net1_noF = find_true_net1_noF,
            mean_p_H0 = mean_p_H0,
            mean_p_H1 = mean_p_H1,
            find_alter_net0 = find_alter_net0,
            # if net1's two displated trees matches alter1, alter2, alter3
            find_alter1_net1_major = find_alter1_net1_major, 
            find_alter2_net1_major = find_alter2_net1_major,
            find_alter3_net1_major = find_alter3_net1_major,
            find_alter1_net1_minor = find_alter1_net1_minor,
            find_alter2_net1_minor = find_alter2_net1_minor,
            find_alter3_net1_minor = find_alter3_net1_minor,
        )
        
    catch e
        println("Error processing file $filepath: $e")
        return nothing
    end
end

"""
Filter SNaQ results where both RF_net1_1_true != 0.0 AND RF_net1_2_true != 0.0
Extract minor gamma values and output to CSV with histogram visualization.
"""
function filter_and_extract_minor_gamma(snaq_summary_dir::String, output_csv::String, output_plot::String)
    """
    Process all SNaQ summary files to extract rows where both RF_net1_1_true != 0.0 
    and RF_net1_2_true != 0.0. Extract the minor gamma value (smaller of gamma_1, gamma_2).
    
    Output CSV columns: parameter_setting, repID, RF_net1_1_true, RF_net1_2_true, minor_gamma
    Also generates a histogram of minor gamma values.
    """
    
    if !isdir(snaq_summary_dir)
        println("Error: Directory '$snaq_summary_dir' not found!")
        return
    end
    
    # Find all CSV files
    csv_files = filter(x -> endswith(x, ".csv"), readdir(snaq_summary_dir))
    
    if isempty(csv_files)
        println("No CSV files found in '$snaq_summary_dir' directory!")
        return
    end
    
    # Collect results from all files
    results = []
    all_minor_gammas = Float64[]
    
    for csv_file in csv_files
        filepath = joinpath(snaq_summary_dir, csv_file)
        
        try
            # Read the CSV file
            df = CSV.read(filepath, DataFrame)
            
            # Extract parameter setting from filename
            parameter_setting = extract_parameter_setting(csv_file)
            
            # Filter rows where both RF_net1_1_true != 0.0 AND RF_net1_2_true != 0.0
            filtered_rows = df[(df.RF_net1_1_true .!= 0.0) .& (df.RF_net1_2_true .!= 0.0), :]
            
            # Process filtered rows
            for row in eachrow(filtered_rows)
                # Calculate minor gamma (smaller of gamma_1 and gamma_2)
                minor_gamma = min(row.gamma_1, row.gamma_2)
                
                # Store the result
                push!(results, (
                    parameter_setting = parameter_setting,
                    repID = row.repID,
                    RF_net1_1_true = row.RF_net1_1_true,
                    RF_net1_2_true = row.RF_net1_2_true,
                    minor_gamma = minor_gamma
                ))
                
                push!(all_minor_gammas, minor_gamma)
            end
            
            println("Processed: $csv_file - Found $(nrow(filtered_rows)) matching replicates")
            
        catch e
            println("Error processing file $filepath: $e")
        end
    end
    
    if isempty(results)
        println("No replicates found matching the filter criteria!")
        return
    end
    
    # Convert results to DataFrame
    result_df = DataFrame(results)
    
    # Write to CSV
    CSV.write(output_csv, result_df)
    println("\nResults written to: $output_csv")
    println("Total matching replicates: $(nrow(result_df))")
    
    # Display summary statistics
    println("\nSummary statistics of minor_gamma values:")
    println("  Mean: $(mean(all_minor_gammas))")
    println("  Median: $(median(all_minor_gammas))")
    println("  Min: $(minimum(all_minor_gammas))")
    println("  Max: $(maximum(all_minor_gammas))")
    println("  Std Dev: $(std(all_minor_gammas))")
    
    return result_df
end

function main()
    """Main function to process all CSV files and create summary."""

    args = parse_commandline()
    
    # Define input and output paths
    snaq_summary_dir = args["input_dir"] 
    output_file = args["output_file"] 
    snaq_visualization_dir = args["visualization_output_dir"] 
    
    # Check if the snaq_summary directory exists
    if !isdir(snaq_summary_dir)
        println("Error: Directory '$snaq_summary_dir' not found!")
        return
    end
    
    # Find all CSV files in the snaq_summary directory
    csv_files = filter(x -> endswith(x, ".csv"), readdir(snaq_summary_dir))
    
    if isempty(csv_files)
        println("No CSV files found in '$snaq_summary_dir' directory!")
        return
    end
    
    println("Found $(length(csv_files)) CSV files to process:")
    for file in csv_files
        println("  - $file")
    end
    
    # Process each CSV file
    results = []
    for csv_file in csv_files
        filepath = joinpath(snaq_summary_dir, csv_file)
        println("Processing: $csv_file")
        
        result = process_csv_file(filepath)
        if result !== nothing
            push!(results, result)
        end
    end
    
    if isempty(results)
        println("No results to summarize!")
        return
    end
    
    # Convert results to DataFrame
    summary_df = DataFrame(results)
    
    # Reorder columns as specified
    column_order = [
        :parameter_setting,
        :H_eq_0,
        :H_eq_1, 
        :H_gt_1,
        :mean_score_H0,
        :mean_score_H1,
        :mean_gamma_1,
        :mean_gamma_2,
        :find_true_net0,
        :find_true_net0_noF,
        :find_true_net1,
        :find_true_net1_noF,
        :mean_p_H0,
        :mean_p_H1,
        :find_alter_net0,
        :find_alter1_net1_major,
        :find_alter2_net1_major,
        :find_alter3_net1_major,
        :find_alter1_net1_minor,
        :find_alter2_net1_minor,
        :find_alter3_net1_minor
    ]
    
    summary_df = select(summary_df, column_order)
    
    # Rename columns to match specification
    rename!(summary_df, 
        :parameter_setting => :parameter_setting,
        :H_eq_0 => Symbol("H=0Accepted"),
        :H_eq_1 => Symbol("H=1Accepted"),
        :H_gt_1 => Symbol("H>1Accepted")
    )
    
    # Write summary to CSV
    CSV.write(output_file, summary_df)
    
    # Print summary
    println("\nSummary completed!")
    println("Results written to: $output_file")
    println("\nSummary statistics:")
    println(summary_df)
    
    # Print some basic statistics
    h0_count = sum(summary_df[!, Symbol("H=0Accepted")])
    h1_count = sum(summary_df[!, Symbol("H=1Accepted")])
    h_gt_1_count = sum(summary_df[!, Symbol("H>1Accepted")])

    # Generate distribution plots
    mkpath(snaq_visualization_dir)
    # plot_column_distributions(snaq_summary_dir, :gamma_1, :gamma_2, 
    #                             snaq_visualization_dir, 
    #                             "SNaQ_minor_gamma_distributions.png")
    # plot_by_ratevar(snaq_summary_dir, :gamma_1, :gamma_2,
    #                snaq_visualization_dir,
    #                "SNaQ_minor_gamma_by_ratevar.png")
    # plot_by_duploss_rate(snaq_summary_dir, :gamma_1, :gamma_2,
    #                     snaq_visualization_dir,
    #                     "SNaQ_minor_gamma_by_duploss_rate.png")

    # plot_overlapping_categories(snaq_summary_dir, :gamma_1, :gamma_2,
    #                            snaq_visualization_dir,
    #                            "SNaQ_overlapping_ratevar.png",
    #                            "ratevar")
    # plot_overlapping_categories(snaq_summary_dir, :gamma_1, :gamma_2,
    #                            snaq_visualization_dir,
    #                            "SNaQ_overlapping_duploss.png",
    #                            "duploss")
    
    # plot_ratevar_overlapping_by_duploss(snaq_summary_dir, :gamma_1, :gamma_2,
    #                                    snaq_visualization_dir,
    #                                    "SNaQ_ratevar_overlapping_by_duploss.png")

    # plot_overlapping_ratevar_by_n_inds_sf(snaq_summary_dir, :gamma_1, :gamma_2, 
    #                                         snaq_visualization_dir, 
    #                                         "SNAQ_minorG_duploss_SF_ninds")
    
    # Generate 12-panel combined plot using R with facet_grid
    println("Generating combined 12-panel plot using R facet_grid...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                 plot_overlapping_ratevar_by_n_inds_sf('$snaq_summary_dir', 
                                     'gamma_1', 
                                     'gamma_2', 
                                     '$snaq_visualization_dir', 
                                     'SNAQ_minorG', 
                                     70)"`)
        println("R plot generated successfully")
    catch e
        @warn "Could not generate R plot: $e"
    end
    
    # Generate gamma threshold summary
    summarize_gamma_by_threshold(snaq_summary_dir, :gamma_1, :gamma_2,
                                snaq_visualization_dir,
                                "snaq_gamma_summary.txt")
                               
    println("Finished generating visualizations in: $snaq_visualization_dir")
    
    println("\nHypothesis support counts across all parameter settings:")
    println("  H=0 (no hybridization): $h0_count replicates")
    println("  H=1 (one hybridization): $h1_count replicates") 
    println("  H>1 (multiple hybridizations): $h_gt_1_count replicates")
    
    # Filter for replicates with both net1 trees != 0 and extract minor gamma
    println("\n" * "="^80)
    println("FILTERING FOR REPLICATES WITH BOTH RF_net1_1_true ≠ 0 AND RF_net1_2_true ≠ 0")
    println("="^80)
    filtered_output_csv = joinpath(snaq_visualization_dir,
                                "SNaQ_minor_gamma_filtered.csv")
    filter_and_extract_minor_gamma(snaq_summary_dir, 
                                    filtered_output_csv, "")
    
    # Generate visualization using R infrastructure
    println("\nGenerating SNaQ minor gamma by tree display visualization...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                 plot_snaq_minor_gamma_by_tree_display('$snaq_summary_dir', 
                    '$snaq_visualization_dir', 
                    'SNaQ_minor_gamma_by_tree_display')"`)
        println("SNaQ minor gamma by tree display visualization generated successfully")
    catch e
        @warn "Could not generate SNaQ visualization: $e"
    end
    
end

# Run the main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end