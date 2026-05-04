using CSV
using DataFrames
using Plots, StatsPlots
using KernelDensity
using StatsBase
using Statistics

# Functions here are often used for quick plotting for sanity check for the results 
# Visualiztaion functions used in the manuscriot is in visual_utilities.R files 

"""
    plot_column_distributions(input_dir::String, 
                            column_name::Symbol, 
                            output_dir::String, 
                            output_filename::String)

Plot distributions of a specified column from all CSV files in input_dir.
Creates a multi-panel figure with one panel per parameter setting 
# Arguments
- `input_dir`: Directory containing CSV files
- `column_name`: Column to plot (as Symbol, e.g., :gamma_1)
- `output_dir`: Directory to save output figure
- `output_filename`: Name of output file (without path)
"""
function plot_column_distributions(input_dir::String, 
                                column_name_1::Symbol, 
                                column_name_2::Symbol,  
                                output_dir::String, 
                                output_filename::String)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    isempty(csv_files) && (println("No CSV files found in $input_dir"); return)
    
    mkpath(output_dir)
    
    # Combine all data into single DataFrame with faceting column
    combined_df = DataFrame(minor_gamma=Float64[], parameter=String[])
    
    for csv_file in csv_files
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        (!hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)) && continue
        
        param = replace(csv_file, r"^(SNaQ-|findgraph-)" => "", r"(-summary)?\.csv$" => "")
        minor_gamma = filter(!isnan, min.(df[!, column_name_1], df[!, column_name_2]))
        isempty(minor_gamma) && continue
        
        append!(combined_df, DataFrame(minor_gamma=minor_gamma, parameter=param))
    end
    
    isempty(combined_df) && (println("No valid data found"); return)
    
    # Create faceted plot with shared axes and single legend
    params = unique(combined_df.parameter)
    n_params = length(params)
    n_cols = min(3, n_params)
    n_rows = ceil(Int, n_params / n_cols)
    
    plots_array = []
    for param in sort(params)
        data = filter(row -> row.parameter == param, combined_df).minor_gamma
        p = histogram(data, bins=0.0:0.01:0.5, title=param, 
                     legend=false, titlefontsize=8, guidefontsize=9, tickfontsize=7,
                     xlim=(0.0, 0.5), margin=3Plots.mm)
        push!(plots_array, p)
    end
    
    final_plot = plot(plots_array..., layout=(n_rows, n_cols), 
                     size=(300*n_cols, 250*n_rows),
                     plot_title="Minor Gamma Distributions",
                     xlabel="Minor Gamma Values", ylabel="Num of replicates")
    
    savefig(final_plot, joinpath(output_dir, output_filename))
    println("Distribution plot saved to: $(joinpath(output_dir, output_filename))")
end


"""
    plot_by_ratevar(input_dir::String, column_name_1::Symbol, 
        column_name_2::Symbol, output_dir::String, output_filename::String)
Plot minor gamma distributions grouped by ratevar level (RVG, RVL, RVN, RVGL).
Each panel shows combined data for all files with the same ratevar.
"""
function plot_by_ratevar(input_dir::String, 
                        column_name_1::Symbol, column_name_2::Symbol, 
                        output_dir::String, output_filename::String)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    isempty(csv_files) && (println("No CSV files found in $input_dir"); return)
    mkpath(output_dir)
    
    # Combine data with ratevar grouping
    combined_df = DataFrame(minor_gamma=Float64[], ratevar=String[])
    
    for csv_file in csv_files
        m = match(r"-(RVG|RVL|RVN|RVGL)-", csv_file)
        m === nothing && continue
        
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        (!hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)) && continue
        
        minor_gamma = filter(!isnan, min.(df[!, column_name_1], df[!, column_name_2]))
        isempty(minor_gamma) && continue
        
        append!(combined_df, DataFrame(minor_gamma=minor_gamma, ratevar=m.captures[1]))
    end
    
    isempty(combined_df) && (println("No valid data found"); return)
    
    # Create faceted plot
    ratevars = sort(unique(combined_df.ratevar))
    plots_array = [histogram(filter(row -> row.ratevar == rv, combined_df).minor_gamma,
                            bins=0.0:0.01:0.5, title=rv, legend=false,
                            titlefontsize=10, xlim=(0.0, 0.5), margin=3Plots.mm)
                   for rv in ratevars]
    
    n_cols = min(3, length(plots_array))
    n_rows = ceil(Int, length(plots_array) / n_cols)
    final_plot = plot(plots_array..., layout=(n_rows, n_cols), 
                     size=(350*n_cols, 300*n_rows),
                     plot_title="Minor Gamma by Rate Variation",
                     xlabel="Minor Gamma Values", ylabel="Num of replicates")
    
    savefig(final_plot, joinpath(output_dir, output_filename))
    println("Ratevar plot saved to: $(joinpath(output_dir, output_filename))")
end


"""
    plot_by_duploss_rate(input_dir::String, column_name_1::Symbol, 
        column_name_2::Symbol, output_dir::String, output_filename::String)

Plot minor gamma distributions grouped by duplication/loss rate level.
Each panel shows combined data for all files with the same dup/loss rate.
"""
function plot_by_duploss_rate(input_dir::String, column_name_1::Symbol, column_name_2::Symbol,
                             output_dir::String, output_filename::String)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    isempty(csv_files) && (println("No CSV files found in $input_dir"); return)
    mkpath(output_dir)
    
    # Combine data with rate grouping
    combined_df = DataFrame(minor_gamma=Float64[], rate=String[])
    
    for csv_file in csv_files
        m = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", csv_file)
        m === nothing && continue
        
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        (!hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)) && continue
        
        minor_gamma = filter(!isnan, min.(df[!, column_name_1], df[!, column_name_2]))
        isempty(minor_gamma) && continue
        
        append!(combined_df, DataFrame(minor_gamma=minor_gamma, rate=m.captures[1]))
    end
    
    isempty(combined_df) && (println("No valid data found"); return)
    
    # Create faceted plot
    rates = sort(unique(combined_df.rate))
    plots_array = [histogram(filter(row -> row.rate == r, combined_df).minor_gamma,
                            bins=0.0:0.01:0.5, title="Rate = $r", legend=false,
                            titlefontsize=10, xlim=(0.0, 0.5), margin=3Plots.mm)
                   for r in rates]
    
    n_cols = min(3, length(plots_array))
    n_rows = ceil(Int, length(plots_array) / n_cols)
    final_plot = plot(plots_array..., layout=(n_rows, n_cols),
                     size=(350*n_cols, 300*n_rows),
                     plot_title="Minor Gamma by Dup/Loss Rate",
                     xlabel="Minor Gamma Values", ylabel="Num of replicates")
    
    savefig(final_plot, joinpath(output_dir, output_filename))
    println("Dup/loss rate plot saved to: $(joinpath(output_dir, output_filename))")
end


"""
    plot_overlapping_categories(input_dir::String, 
                               column_name_1::Symbol, 
                               column_name_2::Symbol, 
                               output_dir::String, 
                               output_filename::String,
                               category_type::String)

Plot overlapping distributions for different parameter categories on a single plot.
Uses different colors and transparency (alpha) to show overlaps between distributions.

# Arguments
- `input_dir`: Directory containing CSV files
- `column_name_1`: First gamma column (as Symbol, e.g., :best_graph_gamma1)
- `column_name_2`: Second gamma column (as Symbol, e.g., :best_graph_gamma2)
- `output_dir`: Directory to save output figure
- `output_filename`: Name of output file (without path)
- `category_type`: Either "ratevar" or "duploss" to specify grouping
"""
function plot_overlapping_categories(input_dir::String, 
                                    column_name_1::Symbol, 
                                    column_name_2::Symbol,  
                                    output_dir::String, 
                                    output_filename::String,
                                    category_type::String)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    
    if isempty(csv_files)
        println("No CSV files found in $input_dir")
        return
    end
    
    mkpath(output_dir)
    
    # Dictionary to store data by category
    category_data = Dict{String, Vector{Float64}}()
    
    for csv_file in csv_files
        # Extract category based on type
        if category_type == "ratevar"
            m = match(r"-(RVG|RVL|RVN|RVGL)-", csv_file)
            if m === nothing
                continue
            end
            category = m.captures[1]
        elseif category_type == "duploss"
            m = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", csv_file)
            if m === nothing
                continue
            end
            rate_str = m.captures[1]
            category = "Rate = $rate_str"
        else
            error("category_type must be either 'ratevar' or 'duploss'")
        end
        
        filepath = joinpath(input_dir, csv_file)
        df = CSV.read(filepath, DataFrame)
        
        if !hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)
            continue
        end
        
        minor_gamma = min.(df[!, column_name_1], df[!, column_name_2])
        minor_gamma = filter(!isnan, minor_gamma)
        
        if haskey(category_data, category)
            append!(category_data[category], minor_gamma)
        else
            category_data[category] = minor_gamma
        end
    end
    
    categories = sort(collect(keys(category_data)))
    
    if isempty(categories)
        println("No valid categories found")
        return
    end
    
    # Define color palette with transparency (blue, orange, dark gray, purple)
    colors = [:dodgerblue, :orangered, :darkgray, :purple]
    alpha_val = 0.5
    
    # Create the overlapping plot
    p = plot(xlabel="Minor Gamma Values", 
             ylabel="Num of replicates", 
             title="Overlapping Distributions by $(uppercasefirst(category_type))",
             legend=:topright,
             legendfontsize=9,
             guidefontsize=10,
             tickfontsize=9,
             titlefontsize=12,
             xlim=(0.0, 0.5),
             size=(800, 600),
             margin = 5Plots.mm)
    
    for (idx, category) in enumerate(categories)
        data = category_data[category]
        color = colors[mod1(idx, length(colors))]
        
        # Plot histogram with transparency
        histogram!(p, data, 
                  bins=0.0:0.01:0.5,
                  alpha=alpha_val,
                  color=color,
                  label=category,
                  linewidth=0)
    end
    
    output_path = joinpath(output_dir, output_filename)
    savefig(p, output_path)
    println("Overlapping $(category_type) plot saved to: $output_path")
end


"""
    plot_ratevar_overlapping_by_duploss(input_dir::String, 
                                       column_name_1::Symbol, 
                                       column_name_2::Symbol,  
                                       output_dir::String, 
                                       output_filename::String)

Plot overlapping ratevar distributions (color-coded with alpha) on separate plots per duploss rate.
Each subplot shows all ratevar levels with different colors and transparency.
"""
function plot_ratevar_overlapping_by_duploss(input_dir::String, 
                                            column_name_1::Symbol, 
                                            column_name_2::Symbol,  
                                            output_dir::String, 
                                            output_filename::String,
                                            y_max::Int=400)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    isempty(csv_files) && (println("No CSV files found"); return)
    mkpath(output_dir)
    
    # Combine all data with grouping columns
    combined_df = DataFrame(minor_gamma=Float64[], rate=String[], ratevar=String[])
    
    for csv_file in csv_files
        m_rate = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", csv_file)
        m_rv = match(r"-(RVG|RVL|RVN|RVGL)-", csv_file)
        (m_rate === nothing || m_rv === nothing) && continue
        
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        (!hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)) && continue
        
        minor_gamma = filter(!isnan, min.(df[!, column_name_1], df[!, column_name_2]))
        isempty(minor_gamma) && continue
        
        append!(combined_df, DataFrame(minor_gamma=minor_gamma, 
                                       rate=m_rate.captures[1], 
                                       ratevar=m_rv.captures[1]))
    end
    
    isempty(combined_df) && (println("No valid data found"); return)
    
    # Create faceted plot with overlapping distributions
    rates = sort(unique(combined_df.rate))
    colors = [:dodgerblue, :orangered, :darkgray, :purple]
    
    plots_array = []
    for rate in rates
        rate_df = filter(row -> row.rate == rate, combined_df)
        p = plot(title="Rate = $rate", legend=:topright, legendfontsize=7,
                xlim=(0.0, 0.5), margin=3Plots.mm)
        for (idx, rv) in enumerate(sort(unique(rate_df.ratevar)))
            data = filter(row -> row.ratevar == rv, rate_df).minor_gamma
            histogram!(p, data, bins=0.0:0.01:0.5, alpha=0.5, 
                      color=colors[mod1(idx, length(colors))], 
                      label=rv, linewidth=0)
        end
        push!(plots_array, p)
    end
    
    n_cols = min(3, length(plots_array))
    n_rows = ceil(Int, length(plots_array) / n_cols)
    final_plot = plot(plots_array..., layout=(n_rows, n_cols),
                     size=(350*n_cols, 300*n_rows),
                     plot_title="Rate Variation by Dup/Loss Rate",
                     xlabel="Minor Gamma Values", ylabel="Num of replicates")
    
    savefig(final_plot, joinpath(output_dir, output_filename))
    println("Ratevar overlapping plot saved to: $(joinpath(output_dir, output_filename))")
end


"""
    plot_overlapping_ratevar_by_n_inds_sf(input_dir::String, 
                                         column_name_1::Symbol, 
                                         column_name_2::Symbol,  
                                         output_dir::String, 
                                         output_filename::String)

Plot overlapping ratevar distributions by (n_inds, SF) combination.
Each figure shows separate subplots per duploss rate with overlapped ratevar levels.
Separate PNG for each (n_inds, SF) pair. Empty subplots are skipped.
"""
function plot_overlapping_ratevar_by_n_inds_sf(input_dir::String, 
                                              column_name_1::Symbol, 
                                              column_name_2::Symbol,  
                                              output_dir::String, 
                                              output_filename::String,
                                              y_max::Int = 400)
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    isempty(csv_files) && (println("No CSV files found"); return)
    mkpath(output_dir)
    
    # Combine all data with grouping columns
    combined_df = DataFrame(minor_gamma=Float64[], n_inds=String[], sf=String[], 
                           rate=String[], ratevar=String[])
    
    for csv_file in csv_files
        m_ninds = match(r"N_ind(\d+)", csv_file)
        m_sf = match(r"SF([\d.]+)", csv_file)
        m_rate = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", csv_file)
        m_rv = match(r"-(RVG|RVL|RVN|RVGL)-", csv_file)
        (m_ninds === nothing || m_sf === nothing || m_rate === nothing || m_rv === nothing) && continue
        
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        (!hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)) && continue
        
        minor_gamma = filter(!isnan, min.(df[!, column_name_1], df[!, column_name_2]))
        isempty(minor_gamma) && continue
        
        append!(combined_df, DataFrame(minor_gamma=minor_gamma,
                                       n_inds=m_ninds.captures[1],
                                       sf=m_sf.captures[1],
                                       rate=m_rate.captures[1],
                                       ratevar=m_rv.captures[1]))
    end
    
    isempty(combined_df) && (println("No valid data found"); return)
    
    # Map ratevar codes to readable names
    ratevar_map = Dict("RVG" => "gene", "RVL" => "lineage", "RVN" => "none")
    ils_map = Dict("1.0" => "high", "0.5" => "low")
    # Colors: blue, orange, dark gray (more distinctive palette)
    colors = Dict("gene" => :dodgerblue, "lineage" => :orangered, "none" => :darkgray)
    
    # Define facet grid structure: rows = (n_inds, ILS), columns = duplication/loss rate
    param_order = [("1", "1.0"), ("1", "0.5"), ("2", "1.0"), ("2", "0.5")]
    rate_order = ["0.0", "0.0003", "0.0004"]
    
    # Create legend plot separately
    legend_plot = plot(framestyle=:none, legend=:inside, legendfontsize=10, 
                      size=(300, 50), margin=0Plots.mm, foreground_color_legend=nothing,
                      background_color_legend=nothing)
    for (rv_name, color) in [("gene", :dodgerblue), ("lineage", :orangered), ("none", :forestgreen)]
        plot!(legend_plot, [], [], label="substitution rate: $(rv_name)", 
              color=color, fillalpha=0.5, linewidth=0, fillrange=0)
    end
    
    # Create all faceted panels with minimal spacing
    all_plots = []
    
    for (row_idx, (n_inds, sf)) in enumerate(param_order)
        param_df = filter(row -> row.n_inds == n_inds && row.sf == sf, combined_df)
        ils_level = get(ils_map, sf, sf)
        
        for (col_idx, rate) in enumerate(rate_order)
            rate_df = filter(row -> row.rate == rate, param_df)
            
            # Show y-axis label only on leftmost column
            show_ylabel = (col_idx == 1)
            # Show x-axis label only on bottom row
            show_xlabel = (row_idx == length(param_order))
            # Show row label (n_inds, ILS) in middle column
            show_row_label = (col_idx == 2)
            
            p = plot(xlim=(0.0, 0.5), ylim=(0, y_max),
                    xlabel=show_xlabel ? "estimated gene flow proportion" : "",
                    ylabel=show_ylabel ? "number of replicates" : "",
                    title=show_row_label ? "n_inds=$(n_inds), ILS=$(ils_level)" : "",
                    legend=false,
                    titlefontsize=10,
                    guidefontsize=10,
                    tickfontsize=9,
                    left_margin=show_ylabel ? 5Plots.mm : 1Plots.mm,
                    right_margin=1Plots.mm,
                    top_margin=show_row_label ? 3Plots.mm : 1Plots.mm,
                    bottom_margin=show_xlabel ? 5Plots.mm : 1Plots.mm,
                    grid=true, gridalpha=0.3)
            
            if !isempty(rate_df)
                for rv_code in ["RVG", "RVL", "RVN"]
                    rv_data = filter(row -> row.ratevar == rv_code, rate_df).minor_gamma
                    isempty(rv_data) && continue
                    
                    rv_name = get(ratevar_map, rv_code, rv_code)
                    color = get(colors, rv_name, :gray)
                    
                    histogram!(p, rv_data, bins=0.0:0.01:0.5, alpha=0.5,
                              color=color, linewidth=0)
                end
            end
            
            push!(all_plots, p)
        end
    end
    
    # Create layout with column headers for duplication/loss rates
    # Custom layout: title row + legend row + column header row + 4 data rows with 3 columns each
    l = @layout [
        a{0.001h}
        b{0.001h}
        c{0.001h}
        grid(4, 3)
    ]
    
    # Title plot - minimal spacing
    title_plot = plot(framestyle=:none, title="variation in substitution rate across: gene, lineage, and none",
                     titlefontsize=14, titlelocation=:center, 
                     top_margin=0Plots.mm, bottom_margin=-2Plots.mm, margin=0Plots.mm)
    
    # Column header plot for duplication/loss rates - closer to plots
    col_header = plot(framestyle=:none, showaxis=false, xlim=(0,1), ylim=(0,1),
                     top_margin=-2Plots.mm, bottom_margin=-3Plots.mm, margin=0Plots.mm)
    for (idx, rate) in enumerate(rate_order)
        x_pos = (idx - 0.5) / 3
        annotate!(col_header, x_pos, 0.5, text("duplication/loss rate = $rate", 11, :center, :bold, :black))
    end
    
    # Combine all plots
    final_plot = plot(title_plot, legend_plot, col_header, all_plots...,
                     layout=l,
                     size=(1200, 1400),
                     margin=0Plots.mm)
    
    out_file = "$(output_filename)_combined.png"
    savefig(final_plot, joinpath(output_dir, out_file))
    println("Saved combined 12-panel figure: $(joinpath(output_dir, out_file))")
end

"""
    parse_parameter_setting(parameter_setting::String)

Parse parameter setting string to extract key parameters.
Returns a dictionary with keys: ratevar, duploss_rate, SF, n_inds

# Example
- Input: "DUP0.0003-LOS0.0003-RVN-N_ind2-SF1.0-genelen1000"
- Output: Dict("ratevar" => "RVN", "duploss_rate" => "0.0003", "SF" => "1.0", "n_inds" => "2")
"""
function parse_parameter_setting(parameter_setting::String)
    params = Dict{String, String}()
    
    # Extract ratevar
    m = match(r"-(RVG|RVL|RVN|RVGL)-", parameter_setting)
    params["ratevar"] = m !== nothing ? m.captures[1] : "unknown"
    
    # Extract duplication/loss rate
    m = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", parameter_setting)
    params["duploss_rate"] = m !== nothing ? m.captures[1] : "unknown"
    
    # Extract SF (scaling factor)
    m = match(r"SF([\d\.]+)", parameter_setting)
    params["SF"] = m !== nothing ? m.captures[1] : "unknown"
    
    # Extract n_inds (number of individuals)
    m = match(r"N_ind(\d+)", parameter_setting)
    params["n_inds"] = m !== nothing ? m.captures[1] : "unknown"
    
    return params
end


"""
    summarize_gamma_by_threshold(input_dir::String, 
                                 column_name_1::Symbol, 
                                 column_name_2::Symbol,
                                 output_dir::String, 
                                 output_filename::String,
                                 thresholds::Vector{Float64}=[0.05, 0.1, 0.25])

Generate summary statistics for gamma values below specified thresholds, grouped by parameters.

# Arguments
- `input_dir`: Directory containing CSV files
- `column_name_1`: First gamma column (e.g., :best_graph_gamma1 or :gamma_1)
- `column_name_2`: Second gamma column (e.g., :best_graph_gamma2 or :gamma_2)
- `output_dir`: Directory to save output file
- `output_filename`: Name of output file (e.g., "findgraph_gamma_summary.txt")
- `thresholds`: Vector of threshold values (default: [0.05, 0.1, 0.25])
"""
function summarize_gamma_by_threshold(input_dir::String, 
                                     column_name_1::Symbol, 
                                     column_name_2::Symbol,
                                     output_dir::String, 
                                     output_filename::String,
                                     thresholds::Vector{Float64}=[0.05, 0.1, 0.25])
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    
    if isempty(csv_files)
        println("No CSV files found in $input_dir")
        return
    end
    
    mkpath(output_dir)
    
    # Store data for each category
    category_stats = Dict{String, Dict{String, Any}}()
    categories = ["ratevar", "duploss_rate", "SF", "n_inds"]
    
    for category in categories
        category_stats[category] = Dict{String, Vector{Float64}}()
    end
    
    # Overall statistics
    all_gamma_values = Float64[]
    
    # Process each CSV file
    for csv_file in csv_files
        filepath = joinpath(input_dir, csv_file)
        df = CSV.read(filepath, DataFrame)
        
        if !hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)
            println("Warning: Columns not found in $csv_file")
            continue
        end
        
        # Extract parameter setting (remove prefix like "SNaQ-" or "findgraph-")
        parameter_setting = replace(csv_file, r"^(SNaQ-|findgraph-)" => "")
        parameter_setting = replace(parameter_setting, r"(-summary)?\.csv$" => "")
        
        # Parse parameters
        params = parse_parameter_setting(parameter_setting)
        
        # Calculate minor gamma
        minor_gamma = min.(df[!, column_name_1], df[!, column_name_2])
        minor_gamma = filter(!isnan, collect(skipmissing(minor_gamma)))
        
        if isempty(minor_gamma)
            continue
        end
        
        # Add to overall statistics
        append!(all_gamma_values, minor_gamma)
        
        # Add to category-specific statistics
        for category in categories
            key = params[category]
            if haskey(category_stats[category], key)
                append!(category_stats[category][key], minor_gamma)
            else
                category_stats[category][key] = copy(minor_gamma)
            end
        end
    end
    
    # Write summary to file
    output_path = joinpath(output_dir, output_filename)
    open(output_path, "w") do io
        println(io, "Gamma summary statistics")
        println(io, "Thresholds: ", join(thresholds, ", "))
        println(io, "="^70)
        println(io)
        
        # Overall statistics
        n_total = length(all_gamma_values)
        
        println(io, "Overall statistics")
        println(io, "-"^70)
        println(io, "Total gamma values: $n_total")
        println(io, "Mean: $(round(mean(all_gamma_values), digits=3))")
        println(io, "Median: $(round(median(all_gamma_values), digits=3))")
        println(io, "Std dev: $(round(std(all_gamma_values), digits=3))")
        println(io)
        
        for threshold in thresholds
            n_below = count(x -> x < threshold, all_gamma_values)
            pct_below = n_total > 0 ? round(n_below / n_total * 100, digits=1) : 0.0
            println(io, "  Below $threshold: $n_below ($pct_below%)")
        end
        println(io)
        println(io)
        
        # Category-specific statistics
        for category in categories
            println(io, "Statistics by $category")
            println(io, "-"^70)
            
            category_data = category_stats[category]
            sorted_keys = sort(collect(keys(category_data)))
            
            for key in sorted_keys
                data = category_data[key]
                n_total_cat = length(data)
                mean_gamma = round(mean(data), digits=3)
                median_gamma = round(median(data), digits=3)
                std_gamma = round(std(data), digits=3)
                
                println(io, "  $key:")
                println(io, "    n=$n_total_cat, mean=$mean_gamma, median=$median_gamma, std=$std_gamma")
                
                for threshold in thresholds
                    n_below_cat = count(x -> x < threshold, data)
                    pct_below_cat = n_total_cat > 0 ? round(n_below_cat / n_total_cat * 100, digits=2) : 0.0
                    println(io, "    <$threshold: $n_below_cat ($pct_below_cat%)")
                end
                println(io)
            end
            
            println(io)
        end
        
        println(io, "="^70)
        println(io, "Source: $input_dir")
    end
    
    println("Gamma summary saved to: $output_path")
    
    # Create CSV table with all statistics
    csv_rows = []
    
    # Add overall statistics row
    n_total = length(all_gamma_values)
    mean_gamma = round(mean(all_gamma_values), digits=4)
    median_gamma = round(median(all_gamma_values), digits=4)
    sd_gamma = round(std(all_gamma_values), digits=4)
    
    overall_row = Dict(
        "category" => "Overall",
        "value" => "all",
        "n" => n_total,
        "mean_gamma" => mean_gamma,
        "median_gamma" => median_gamma,
        "sd_gamma" => sd_gamma
    )
    
    for threshold in thresholds
        n_below = count(x -> x < threshold, all_gamma_values)
        pct_below = n_total > 0 ? round(n_below / n_total * 100, digits=1) : 0.0
        col_name = "gamma_below_$(threshold)"
        overall_row[col_name] = pct_below
    end
    
    push!(csv_rows, overall_row)
    
    # Add category-specific statistics
    for category in categories
        category_data = category_stats[category]
        sorted_keys = sort(collect(keys(category_data)))
        
        for key in sorted_keys
            data = category_data[key]
            n_total_cat = length(data)
            mean_gamma_cat = round(mean(data), digits=3)
            median_gamma_cat = round(median(data), digits=3)
            sd_gamma_cat = round(std(data), digits=3)
            
            row = Dict(
                "category" => category,
                "value" => key,
                "n" => n_total_cat,
                "mean_gamma" => mean_gamma_cat,
                "median_gamma" => median_gamma_cat,
                "sd_gamma" => sd_gamma_cat
            )
            
            for threshold in thresholds
                n_below_cat = count(x -> x < threshold, data)
                pct_below_cat = n_total_cat > 0 ? round(n_below_cat / n_total_cat * 100, digits=1) : 0.0
                col_name = "gamma_below_$(threshold)"
                row[col_name] = pct_below_cat
            end
            
            push!(csv_rows, row)
        end
    end
    
    # Convert to DataFrame and save as CSV
    summary_df = DataFrame(csv_rows)
    
    # Reorder columns for better readability
    col_order = ["category", "value", "n", "mean_gamma", "median_gamma", "sd_gamma"]
    for threshold in thresholds
        push!(col_order, "gamma_below_$(threshold)")
    end
    select!(summary_df, col_order)
    
    # Generate CSV filename
    csv_filename = replace(output_filename, r"\.(txt|log)$" => ".csv")
    csv_output_path = joinpath(output_dir, csv_filename)
    
    CSV.write(csv_output_path, summary_df)
    println("Gamma summary CSV saved to: $csv_output_path")
end


"""
    summarize_WR_by_threshold(input_dir::String, 
                              column_name::Symbol,
                              output_dir::String, 
                              output_filename::String,
                              metric_name::String,
                              percentiles::Vector{Float64}=[0.95, 0.99])

Generate summary statistics for worst residual (WR) values with percentile thresholds, 
grouped by parameters. Works for single-column metrics like H0_best_tree_WR and H1_best_graph_WR.

# Arguments
- `input_dir`: Directory containing CSV files
- `column_name`: Column to analyze (e.g., :H0_best_tree_WR or :H1_best_graph_WR)
- `output_dir`: Directory to save output file
- `output_filename`: Name of output file (e.g., "findgraph_H0_WR_summary.txt")
- `metric_name`: Display name for the metric (e.g., "H0 Best Tree WR")
- `percentiles`: Vector of percentiles to calculate (default: [0.95, 0.99] for 95th and 99th percentiles)
"""
function summarize_WR_by_threshold(input_dir::String, 
                                   column_name::Symbol,
                                   output_dir::String, 
                                   output_filename::String,
                                   metric_name::String,
                                   percentiles::Vector{Float64}=[0.95, 0.99])
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))
    
    if isempty(csv_files)
        println("No CSV files found in $input_dir")
        return
    end
    
    mkpath(output_dir)
    
    # Store data for each category
    category_stats = Dict{String, Dict{String, Any}}()
    categories = ["ratevar", "duploss_rate", "SF", "n_inds"]
    
    for category in categories
        category_stats[category] = Dict{String, Vector{Float64}}()
    end
    
    # Overall statistics
    all_wr_values = Float64[]
    
    # Process each CSV file
    for csv_file in csv_files
        filepath = joinpath(input_dir, csv_file)
        df = CSV.read(filepath, DataFrame)
        
        if !hasproperty(df, column_name)
            println("Warning: Column $(column_name) not found in $csv_file")
            continue
        end
        
        # Extract parameter setting (remove prefix like "SNaQ-" or "findgraph-")
        parameter_setting = replace(csv_file, r"^(SNaQ-|findgraph-)" => "")
        parameter_setting = replace(parameter_setting, r"(-summary)?\.csv$" => "")
        
        # Parse parameters
        params = parse_parameter_setting(parameter_setting)
        
        # Extract WR values and filter NaN/missing
        wr_values = filter(!isnan, collect(skipmissing(df[!, column_name])))
        
        if isempty(wr_values)
            continue
        end
        
        # Add to overall statistics
        append!(all_wr_values, wr_values)
        
        # Add to category-specific statistics
        for category in categories
            key = params[category]
            if haskey(category_stats[category], key)
                append!(category_stats[category][key], wr_values)
            else
                category_stats[category][key] = copy(wr_values)
            end
        end
    end
    
    # Write summary to file
    output_path = joinpath(output_dir, output_filename)
    open(output_path, "w") do io
        println(io, "$metric_name summary statistics")
        println(io, "Percentiles: ", join(map(p -> "$(round(p*100, digits=0))%", percentiles), ", "))
        println(io, "="^70)
        println(io)
        
        # Overall statistics
        n_total = length(all_wr_values)
        
        println(io, "Overall statistics")
        println(io, "-"^70)
        println(io, "Total WR values: $n_total")
        if !isempty(all_wr_values)
            println(io, "Mean: $(round(mean(all_wr_values), digits=4))")
            println(io, "Median: $(round(median(all_wr_values), digits=4))")
            println(io, "Std dev: $(round(std(all_wr_values), digits=4))")
        else
            println(io, "No valid data found for this metric")
        end
        println(io)
        
        # Calculate percentile values for overall data
        for pct in percentiles
            if !isempty(all_wr_values)
                q_value = round(quantile(all_wr_values, pct), digits=4)
                pct_label = round(Int, pct * 100)
                println(io, "  $pct_label-th percentile: $q_value ($(pct_label)% of values below this)")
            end
        end
        println(io)
        println(io)
        
        # Category-specific statistics
        for category in categories
            println(io, "Statistics by $category")
            println(io, "-"^70)
            
            category_data = category_stats[category]
            sorted_keys = sort(collect(keys(category_data)))
            
            for key in sorted_keys
                data = category_data[key]
                n_total_cat = length(data)
                if !isempty(data)
                    mean_wr = round(mean(data), digits=4)
                    median_wr = round(median(data), digits=4)
                    std_wr = round(std(data), digits=4)
                    
                    println(io, "  $key:")
                    println(io, "    n=$n_total_cat, mean=$mean_wr, median=$median_wr, std=$std_wr")
                    
                    # Calculate percentile values for this category
                    for pct in percentiles
                        if !isempty(data)
                            q_value = round(quantile(data, pct), digits=4)
                            pct_label = round(Int, pct * 100)
                            println(io, "    $pct_label-th percentile: $q_value")
                        end
                    end
                end
                println(io)
            end
            
            println(io)
        end
        
        println(io, "="^70)
        println(io, "Source: $input_dir")
    end
    
    println("$metric_name summary saved to: $output_path")
    
    # Create CSV table with all statistics
    csv_rows = []
    
    # Add overall statistics row
    n_total = length(all_wr_values)
    overall_row = Dict(
        "category" => "Overall",
        "value" => "all",
        "n" => n_total
    )
    
    if !isempty(all_wr_values)
        mean_wr = round(mean(all_wr_values), digits=4)
        median_wr = round(median(all_wr_values), digits=4)
        sd_wr = round(std(all_wr_values), digits=4)
        overall_row["mean_WR"] = mean_wr
        overall_row["median_WR"] = median_wr
        overall_row["sd_WR"] = sd_wr
        
        # Calculate percentile values for overall data
        for pct in percentiles
            q_value = round(quantile(all_wr_values, pct), digits=4)
            pct_label = round(Int, pct * 100)
            col_name = "percentile_$(pct_label)"
            overall_row[col_name] = q_value
        end
    end
    
    push!(csv_rows, overall_row)
    
    # Add category-specific statistics
    for category in categories
        category_data = category_stats[category]
        sorted_keys = sort(collect(keys(category_data)))
        
        for key in sorted_keys
            data = category_data[key]
            n_total_cat = length(data)
            
            row = Dict(
                "category" => category,
                "value" => key,
                "n" => n_total_cat
            )
            
            if !isempty(data)
                mean_wr_cat = round(mean(data), digits=4)
                median_wr_cat = round(median(data), digits=4)
                sd_wr_cat = round(std(data), digits=4)
                
                row["mean_WR"] = mean_wr_cat
                row["median_WR"] = median_wr_cat
                row["sd_WR"] = sd_wr_cat
                
                # Calculate percentile values for this category
                for pct in percentiles
                    q_value = round(quantile(data, pct), digits=4)
                    pct_label = round(Int, pct * 100)
                    col_name = "percentile_$(pct_label)"
                    row[col_name] = q_value
                end
            end
            
            push!(csv_rows, row)
        end
    end
    
    # Convert to DataFrame and save as CSV
    summary_df = DataFrame(csv_rows)
    
    # Reorder columns for better readability (only include columns that exist)
    col_order = ["category", "value", "n"]
    if hasproperty(summary_df, :mean_WR)
        push!(col_order, "mean_WR")
    end
    if hasproperty(summary_df, :median_WR)
        push!(col_order, "median_WR")
    end
    if hasproperty(summary_df, :sd_WR)
        push!(col_order, "sd_WR")
    end
    for pct in percentiles
        pct_label = round(Int, pct * 100)
        col_name = "percentile_$(pct_label)"
        if hasproperty(summary_df, Symbol(col_name))
            push!(col_order, col_name)
        end
    end
    if length(col_order) > 3  # Only reorder if there are extra columns beyond the base three
        select!(summary_df, col_order)
    end
    
    # Generate CSV filename
    csv_filename = replace(output_filename, r"\.(txt|log)$" => ".csv")
    csv_output_path = joinpath(output_dir, csv_filename)
    
    CSV.write(csv_output_path, summary_df)
    println("$metric_name CSV summary saved to: $csv_output_path")
end

