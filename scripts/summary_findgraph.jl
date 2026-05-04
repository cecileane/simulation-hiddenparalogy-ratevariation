#!/usr/bin/env julia

using CSV
using DataFrames
using Statistics
using ArgParse

using PrettyTables

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir"
            help = "Directory containing findgraph summary files"
            default = "findgraph_summary"
        "--output_file"
            help = "Output CSV file path"
            default = "./results/findgraph_summary.csv"
        "--visualization_output_dir"
            help = "Output directory for visualization files"
            default = "visualization_results/findgraph" 
        "--taxon_recovery_output"
            help = "Output CSV for long-format taxon recovery summary"
            default = "results/findgraph_taxon_recovery.csv"
    end
    return parse_args(s)
end

"""
Extract parameter name root from filename
Example: "findgraph-DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5.csv"
→ "DUP0.0004-LOS0.0004-RVN-N_ind1-SF0.5"
"""
function extract_paramname(filename::String)
    # Remove "findgraph-" prefix and ".csv" suffix
    paramname = replace(filename, r"^findgraph-" => "")
    paramname = replace(paramname, r"\.csv$" => "")
    return paramname
end

"""
Extract WR threshold from column name
Example: "best_k_new_WR_3.7" -> "3.7"
"""
function extract_wr_threshold(colname::String)
    pattern = r"best_k_new_WR_([0-9.]+)"
    m = Base.match(pattern, colname)
    if m !== nothing
        return m.captures[1]
    end
    return nothing
end

"""
Find all best_k_new_WR_{threshold} columns in the given dataframe
"""
function find_wr_threshold_columns(df::DataFrame)
    wr_columns = Dict{String, String}()  # threshold -> column_name
    for colname in names(df)
        threshold = extract_wr_threshold(colname)
        if threshold !== nothing
            wr_columns[threshold] = colname
        end
    end
    return wr_columns
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

    # Process WR threshold columns
    wr_thresholds = find_wr_threshold_columns(df)
    wr_statistics = NamedTuple()
    
    # Count and store statistics for each WR threshold
    for (threshold, colname) in wr_thresholds
        h0_wr = count(x -> (x == 0 || x == "0"), df[!, colname])
        h1_wr = count(x -> (x == 1 || x == "1"), df[!, colname])
        bt1_wr = count(x -> (x == ">1" || string(x) == ">1"), df[!, colname])
        
        wr_statistics = merge(wr_statistics, (
            Symbol("H0Accepted_WR_$(threshold)") => h0_wr,
            Symbol("H1Accepted_WR_$(threshold)") => h1_wr,
            Symbol("BT1Accepted_WR_$(threshold)") => bt1_wr,
        ))
    end
    
    # Percentage of true_tree_found_in_H0 = True
    true_tree_h0_count = count(df.true_tree_found_in_H0 .== true)
    pct_true_tree_h0 = true_tree_h0_count / total_reps * 100
    
    # Average num_blocks (handling missing values)
    avg_num_blocks = mean(skipmissing(df.num_blocks))
    
    # Average gamma values (mean of means across replicates)
    avg_gamma1_H1 = mean(skipmissing(df.avg_gamma1_H1))
    avg_gamma2_H1 = mean(skipmissing(df.avg_gamma2_H1))
    
    # Average of best graph gamma values across replicates
    avg_best_gamma1 = mean(skipmissing(df.best_graph_gamma1))
    avg_best_gamma2 = mean(skipmissing(df.best_graph_gamma2))
    
    # Percentage of true tree displayed in H1 (either major or minor)
    true_tree_h1_count = count((df.True_tree_displayed_H1_major .== true) .| 
                               (df.True_tree_displayed_H1_minor .== true))
    pct_true_tree_h1 = true_tree_h1_count / total_reps * 100

    # Same but without F: pct reps finding true tree in H1
    true_tree_h1_count_noF = count(
        (df.True_tree_noF_displayed_H1_major .== true) .|
        (df.True_tree_noF_displayed_H1_minor .== true))
    pct_true_tree_h1_noF = true_tree_h1_count_noF / total_reps * 100 
    
    # Summary statistics for H0_trees_found
    mean_H0_trees = mean(df.H0_trees_found)
    median_H0_trees = median(df.H0_trees_found)
    sd_H0_trees = std(df.H0_trees_found)
    min_H0_trees = minimum(df.H0_trees_found)
    max_H0_trees = maximum(df.H0_trees_found)
    
    # Summary statistics for H1_graphs_found
    mean_H1_graphs = mean(df.H1_graphs_found)
    median_H1_graphs = median(df.H1_graphs_found)
    sd_H1_graphs = std(df.H1_graphs_found)
    min_H1_graphs = minimum(df.H1_graphs_found)
    max_H1_graphs = maximum(df.H1_graphs_found)
    
    istrue(x) = x == "True" || x == true
    count_H0_best_is_true = count(istrue, df.H0_best_tree_is_true_tree)
    count_H0_best_is_true_noF = count(istrue, df.H0_best_tree_is_true_tree_noF)
    count_H1_best_displays_true = count(
        istrue, df.H1_best_graph_displayed_true_tree)
    count_H1_best_displays_true_noF = count(
        istrue, df.H1_best_graph_displayed_true_tree_noF)
    
    return merge((
        paramname_root = paramname,
        total_replicates = total_reps,
        H0Accepted = h0_accepted,
        H1Accepted = h1_accepted,
        BT1Accepted = bt1_accepted,
        pct_true_tree_H0 = pct_true_tree_h0,
        avg_num_blocks = avg_num_blocks,
        avg_gamma1_H1 = avg_gamma1_H1,
        avg_gamma2_H1 = avg_gamma2_H1,
        avg_best_gamma1 = avg_best_gamma1,
        avg_best_gamma2 = avg_best_gamma2,
        pct_true_tree_H1 = pct_true_tree_h1,
        pct_true_tree_H1_noF = pct_true_tree_h1_noF,
        # H0_trees_found statistics
        mean_H0_trees_found = mean_H0_trees,
        median_H0_trees_found = median_H0_trees,
        sd_H0_trees_found = sd_H0_trees,
        min_H0_trees_found = min_H0_trees,
        max_H0_trees_found = max_H0_trees,
        # H1_graphs_found statistics
        mean_H1_graphs_found = mean_H1_graphs,
        median_H1_graphs_found = median_H1_graphs,
        sd_H1_graphs_found = sd_H1_graphs,
        min_H1_graphs_found = min_H1_graphs,
        max_H1_graphs_found = max_H1_graphs,
        # Counts of replicates with true tree/graph found
        count_H0_best_is_true_tree = count_H0_best_is_true,
        count_H0_best_is_true_tree_noF = count_H0_best_is_true_noF,
        count_H1_best_displays_true_tree = count_H1_best_displays_true,
        count_H1_best_displays_true_tree_noF = count_H1_best_displays_true_noF
    ), wr_statistics)
end

"""
Validate that all files have consistent WR thresholds
"""
function validate_wr_thresholds(input_dir::String)
    # Find all files starting with "findgraph"
    all_files = readdir(input_dir)
    findgraph_files = filter(
        f -> startswith(f, "findgraph") && endswith(f, ".csv"), all_files)
    
    all_thresholds = Set()
    
    for filename in findgraph_files
        filepath = joinpath(input_dir, filename)
        df = CSV.read(filepath, DataFrame)
        thresholds = keys(find_wr_threshold_columns(df))
        
        if isempty(thresholds)
            @warn "File $filename has no best_k_new_WR columns"
        else
            union!(all_thresholds, thresholds)
        end
    end
    
    # Check consistency
    if length(all_thresholds) > 1
        error("Inconsistent WR thresholds in $input_dir: " *
            "$(sort(collect(all_thresholds))). " *
            "All files must have the same threshold.")
    elseif length(all_thresholds) == 1
        threshold = collect(all_thresholds)[1]
        println("Using WR threshold: $threshold")
    else
        @warn "No best_k_new_WR columns found in any files"
    end
    
    return all_thresholds
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
    required_cols = [:repID, :best_graph_hybrid_taxon,
        :best_graph_major_donor, :H1_best_graph_displayed_true_tree]

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

    # Drop rows where no H1 graph was found ("NA" in all taxon columns).
    # Keeping these would conflate "not computed" with "computed False".
    combined = filter(row ->
        !(string(row.best_graph_hybrid_taxon) == "NA" &&
          string(row.best_graph_major_donor)  == "NA" &&
          string(row.best_graph_minor_donor)  == "NA"), combined)

    # STEP 2: Recovery flag from H1_best_graph_displayed_true_tree
    # (always set unconditionally in findgraphs_postprocess.jl)
    combined[!, :true_tree_recovered] = map(
        x -> x == true || x == "True",
        combined.H1_best_graph_displayed_true_tree)

    # STEP 3: Explode taxon role columns (one row per taxon);
    # "NA" rows are filtered above so they will not appear as a taxon name.
    long_all = vcat(
        explode_role_column(combined, :best_graph_hybrid_taxon, "hybrid_taxon"),
        explode_role_column(combined, :best_graph_major_donor,  "major_donor"),
        explode_role_column(combined, :best_graph_minor_donor,  "minor_donor")
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

"""
Main function
"""
function main()

    mkpath("./results")

    args = parse_commandline()
    
    input_dir = args["input_dir"]
    output_file = args["output_file"]
    visualization_output_dir = args["visualization_output_dir"] 
    
    # Check if input directory exists
    if !isdir(input_dir)
        error("Input directory does not exist: $input_dir")
    end
    
    # Validate that all files have consistent WR thresholds
    wr_thresholds = validate_wr_thresholds(input_dir)
    
    # Find all files starting with "findgraph"
    all_files = readdir(input_dir)
    findgraph_files = filter(
        f -> startswith(f, "findgraph") && endswith(f, ".csv"), all_files)
    
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

    #--------- Generate distribution plots ---------# 
    # Plot gamma distribution 
    mkpath(visualization_output_dir)
    
    # Generate 12-panel combined plot using R with facet_grid
    println("Generating combined 12-panel plot using R facet_grid...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
        plot_overlapping_ratevar_by_n_inds_sf('$input_dir', 
                                'best_graph_gamma1', 
                                'best_graph_gamma2', 
                                '$visualization_output_dir', 
                                'findgraph_minorG', 
                                50)"`
                                )
        println("R plot generated successfully")
    catch e
        @warn "Could not generate R plot: $e"
    end
    
    println("Generating the true graph WR distribution plot...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                plot_WR_distributions('$input_dir', 
                                    'true_tree_wr', 
                                    '$visualization_output_dir',
                                    'findgraph_true_tree_WR_distributions',
                                    NULL, 
                                    25)"` ) 
        println("True tree WR distribution plot generated successfully")
    catch e
        @warn "Could not generate true tree WR plot: $e"
    end
    
    println("Generating true tree WR percentiles by rate variation...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                plot_WR_percentiles_by_rate_variation('$input_dir', 
                            'true_tree_wr', 
                            '$visualization_output_dir', 
                            'findgraph_true_tree_WR_percentiles_by_ratevar',
                            rv_colors = rv_colors)"`) 
        println("True tree WR percentiles plot generated successfully")
    catch e
        @warn "Could not generate true tree WR percentiles plot: $e"
    end

    println("Generating true tree WR percentiles jitter-combined plot...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R');
                plot_WR_percentiles_jitter_combined('$input_dir',
                    'true_tree_wr', '$visualization_output_dir',
                    'findgraph_true_tree_WR_95-percentiles_jitter_combined',
                    ymax = 9, rv_colors = rv_colors)"`)
        println("True tree WR percentiles jitter-combined plot generated")
    catch e
        @warn "Could not generate true tree WR jitter-combined plot: $e"
    end

    # Generate gamma threshold summary
    summarize_gamma_by_threshold(input_dir, :best_graph_gamma1, 
                                :best_graph_gamma2,
                                visualization_output_dir,
                                "findgraph_gamma_summary.txt")
    
    # Generate H0 best tree worst residual summary
    summarize_WR_by_threshold(input_dir, :H0_best_tree_WR,
                             visualization_output_dir,
                             "findgraph_H0_WR_summary.txt",
                             "H0 Best Tree Worst Residual",
                             [0.95, 0.99])
    
    # Generate H1 best graph worst residual summary
    summarize_WR_by_threshold(input_dir, :H1_best_graph_WR,
                             visualization_output_dir,
                             "findgraph_H1_WR_summary.txt",
                             "H1 Best Graph Worst Residual",
                             [0.95, 0.99])
    summarize_WR_by_threshold(input_dir, :true_tree_wr,
                             visualization_output_dir,
                             "findgraph_true_tree_WR_summary.txt",
                             "True Tree Worst Residual",
                             [0.95, 0.99])

    # Generate findgraph statistics by tree display visualization using R
    println("\nGenerating findgraph minor gamma (gamma2) by tree display...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                 plot_statistics_by_tree_display('$input_dir', 
                    'best_graph_gamma2', 
                    'H1_best_graph_displayed_true_tree',
                    '$visualization_output_dir', 
                    'findgraph_gamma2_by_tree_display',
                    'Findgraph H1 Minor Gamma (gamma2) by True Tree Display',
                    'Minor Gamma (gamma2)')"`)
        println("Findgraph gamma2 by tree display visualization generated")
    catch e
        @warn "Could not generate findgraph gamma2 visualization: $e"
    end
    
    println("Generating findgraph H1 best graph WR by tree display...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R'); 
                 plot_statistics_by_tree_display('$input_dir', 
                    'H1_best_graph_WR', 
                    'H1_best_graph_displayed_true_tree',
                    '$visualization_output_dir', 
                    'findgraph_H1_WR_by_tree_display',
                    'Findgraph H1 Best Graph WR by True Tree Display',
                    'Worst Residual (WR)')"`)
        println("Findgraph H1 WR by tree display visualization generated")
    catch e
        @warn "Could not generate findgraph H1 WR visualization: $e"
    end
                            
    println("Finished generating visualizations in: $visualization_output_dir")
    #--------- End of distribution plots ---------# 

    println("\nSummary written to: $output_file")
    println("Processed $(nrow(results_df)) files successfully")

    # Taxon recovery summary — long-format table by (param_setting, taxon, role)
    compute_taxon_recovery_summary(input_dir, args["taxon_recovery_output"])

end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
