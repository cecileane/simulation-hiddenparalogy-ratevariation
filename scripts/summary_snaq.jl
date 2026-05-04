#!/usr/bin/env julia

using CSV
using DataFrames
using Statistics
using ArgParse

using Plots
using PrettyTables

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
        "--taxon_recovery_output"
            help = "Output CSV for long-format taxon recovery summary"
            default = "results/SNaQ_taxon_recovery.csv"
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
            h0_support, h1_support, h_gt_1_support =
                classify_hypothesis(row.p_H0, row.p_H1)
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
        
        # find_true_net1: count reps where major OR minor tree matched
        find_true_net1 = sum(
            (df.RF_net1_1_true .== 0.0) .| (df.RF_net1_2_true .== 0.0))
        find_true_net1_noF = sum(
            (df.RF_net1_1_true_noF .== 0.0) .| (df.RF_net1_2_true_noF .== 0.0))

        # Check if net0's RF distance to alter1, alter2, alter3 together
        # for alters: only count rows where the true tree was NOT recovered
        valid_rows_net0 = df.RF_net1_1_true .!= 0.0
        valid_rows_net1 = (df.RF_net1_1_true .!= 0.0) .&&
            (df.RF_net1_2_true .!= 0.0)
        
        find_alter1_net0 = sum((df.RF_net0_alter1 .== 0.0) .& valid_rows_net0)
        find_alter2_net0 = sum((df.RF_net0_alter2 .== 0.0) .& valid_rows_net0)
        find_alter3_net0 = sum((df.RF_net0_alter3 .== 0.0) .& valid_rows_net0)
        find_alter_net0 = find_alter1_net0 + find_alter2_net0 + find_alter3_net0

        # check if net1's major hybrid's RF distance to alter1, alter2, alter3
        find_alter1_net1_major = sum(
            (df.RF_net1_1_alter1 .== 0.0) .& valid_rows_net1)
        find_alter2_net1_major = sum(
            (df.RF_net1_1_alter2 .== 0.0) .& valid_rows_net1)
        find_alter3_net1_major = sum(
            (df.RF_net1_1_alter3 .== 0.0) .& valid_rows_net1)
        find_alter1_net1_minor = sum(
            (df.RF_net1_2_alter1 .== 0.0) .& valid_rows_net1)
        find_alter2_net1_minor = sum(
            (df.RF_net1_2_alter2 .== 0.0) .& valid_rows_net1)
        find_alter3_net1_minor = sum(
            (df.RF_net1_2_alter3 .== 0.0) .& valid_rows_net1)
        
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
            find_true_net1 = find_true_net1,
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
function filter_and_extract_minor_gamma(snaq_summary_dir::String,
                                        output_csv::String, output_plot::String)
    """
    filter_and_extract_minor_gamma(snaq_summary_dir, output_csv, output_plot)
    Extract rows where both RF_net1_1_true != 0.0 and RF_net1_2_true != 0.0,
    then write minor gamma (min of gamma_1, gamma_2) and a histogram to disk.
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
            
            filtered_rows = df[(df.RF_net1_1_true .!= 0.0) .&
                (df.RF_net1_2_true .!= 0.0), :]
            
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
            
            println("$csv_file: $(nrow(filtered_rows)) matching reps")
            
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

"""
Explode a comma-separated role column into one row per taxon.
Carries param_setting, repID, and true_tree_recovered.
"""
function explode_role_column(df::DataFrame, col::Symbol, role_label::String)
    rows = NamedTuple[]
    for row in eachrow(df)
        for t in strip.(split(string(row[col]), ","))
            isempty(t) && continue
            push!(rows, (param_setting        = row.param_setting,
                         repID               = row.repID,
                         taxon               = t,
                         role                = role_label,
                         true_tree_recovered = row.true_tree_recovered))
        end
    end
    return DataFrame(rows)
end

"""
Long-format stratified taxon recovery summary from all CSV files in input_dir.
Steps: read+tag → recovery flag → explode roles → group+summarize → print+save.
"""
function compute_taxon_recovery_summary(input_dir::String, output_file::String)
    required_cols = [:repID, :hybrid_taxon, :major_donor, :minor_donor,
                     :RF_net1_1_true, :RF_net1_2_true]

    # STEP 1: Read and tag each CSV file; skip files missing required columns
    csv_files = filter(f -> endswith(f, ".csv"), readdir(input_dir))
    tagged_dfs = DataFrame[]
    for f in csv_files
        df = CSV.read(joinpath(input_dir, f), DataFrame)
        missing_cols = [c for c in required_cols if !(c in propertynames(df))]
        if !isempty(missing_cols)
            @warn "Skipping $f — missing: $(join(string.(missing_cols), ", "))"
            continue
        end
        df[!, :param_setting] = fill(replace(f, r"\.csv$" => ""), nrow(df))
        push!(tagged_dfs, select(df, vcat(required_cols, [:param_setting])))
    end
    if isempty(tagged_dfs)
        @warn "No usable files found for taxon recovery analysis — skipping."
        return
    end

    # Stack all DataFrames; keep only shared columns
    combined = vcat(tagged_dfs...; cols=:intersect)

    # Filter out replicates where H1 never ran (both RF values are NaN);
    # NaN == 0.0 is false in Julia, so keeping these would silently count
    # missing replicates as "not recovered" rather than as missing data.
    combined = filter(
        row -> !isnan(row.RF_net1_1_true) || !isnan(row.RF_net1_2_true),
        combined)

    # Also drop rows where taxon columns were never filled ("NA" default)
    combined = filter(row -> !(string(row.hybrid_taxon) == "NA" &&
                               string(row.major_donor)  == "NA" &&
                               string(row.minor_donor)  == "NA"), combined)

    # STEP 2: Recovery flag — true if either RF orientation equals 0
    # NaN RF values (only one tree displayed) are treated as not matching.
    combined[!, :true_tree_recovered] =
        (combined.RF_net1_1_true .== 0.0) .| (combined.RF_net1_2_true .== 0.0)

    # STEP 3: Explode hybrid_taxon/major_donor/minor_donor (one row per taxon);
    # explode_role_column already skips empty strings — "NA" entries are
    # filtered above so they will not appear as a taxon name.
    long_all = vcat(
        explode_role_column(combined, :hybrid_taxon, "hybrid_taxon"),
        explode_role_column(combined, :major_donor,  "major_donor"),
        explode_role_column(combined, :minor_donor,  "minor_donor")
    )
    # Deduplicate: each (param_setting, repID, taxon, role) counted at most once
    long_dedup = unique(long_all, [:param_setting, :repID, :taxon, :role])

    # STEP 4: Stratified summary grouped by (param_setting, taxon, role)
    summary = combine(
        groupby(long_dedup, [:param_setting, :taxon, :role]),
        :true_tree_recovered => sum             => :count_recovered,
        :true_tree_recovered => (x -> sum(.!x)) => :count_not_recovered,
        :true_tree_recovered => length          => :total
    )
    summary[!, :pct_correct] = map(
        (r, t) -> t == 0 ? NaN : round(100.0 * r / t, digits=1),
        summary.count_recovered, summary.total
    )
    sort!(summary, [:param_setting, order(:total, rev=true)])

    # STEP 5: Print preview of first 20 rows and save full table to CSV
    n_preview = min(20, nrow(summary))
    println("\nTaxon recovery summary — first $n_preview rows:")
    pretty_table(summary[1:n_preview, :];
        header=names(summary), tf=tf_unicode_rounded)
    mkpath(dirname(output_file) == "" ? "." : dirname(output_file))
    CSV.write(output_file, summary)
    println("Taxon recovery ($(nrow(summary)) rows) → $output_file")
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
    println("FILTERING: RF_net1_1_true ≠ 0 AND RF_net1_2_true ≠ 0")
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
        println("SNaQ minor gamma by tree display visualization generated")
    catch e
        @warn "Could not generate SNaQ visualization: $e"
    end

    # Taxon recovery summary — long-format table by (param_setting, taxon, role)
    compute_taxon_recovery_summary(snaq_summary_dir,
        args["taxon_recovery_output"])
    
end

# Run the main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end