#!/usr/bin/env julia
# ============================================================================
# scripts/summary_phynest.jl
#
# Purpose : Aggregate per-parameter PhyNEST CSVs into one cross-setting summary
#           table (model choice, scores, gamma estimates, tree-recovery counts,
#           search stability) and produce diagnostic plots via R.
#           Model choice between h=0 and h=1 uses the per-site score
#           difference T = (score_H0 - score_H1) / M, with M the number of
#           sites of the replicate's alignment (the raw difference grows
#           linearly with M, and M is smaller when fewer loci were
#           retained; columns T_per_site and T of the input tables), in two
#           ways: the naive rule T > 0, and the calibrated rule T > T*. As
#           for the WR <= 3.7 threshold of find_graphs, T* is the 95th
#           percentile of T pooled over the settings without lineage rate
#           variation (no and gene-specific rates, all duplication/loss
#           rates);
#           because T differs between ILS levels, one T* is computed per
#           ILS level (SF) and number of individuals. T* is computed here,
#           after all PhyNEST runs, and saved to its own file. It is
#           specific to this design: 8 taxa (70 quartets), 100 search runs
#           per model.
# Inputs  : phynest_summary/PhyNEST-<paramname>-summary.csv  (one per setting)
# Outputs : results/PhyNEST_summary.csv             (cross-setting table)
#           results/PhyNEST_T_threshold.csv         (T* per ILS level)
#           results/PhyNEST_T_summary_by_factor.csv     (T by factor level)
#           results/PhyNEST_gamma_summary_by_factor.csv (gamma_2 by level)
#           results/PhyNEST_model_choice_by_factor.csv  (% h=1 by level)
#           visualization_results/phynest/*             (diagnostic plots)
# Usage   : julia --project=. scripts/summary_phynest.jl
#           (optional overrides: --input_dir, --output_file, --threshold_file,
#                                --factor_summary_file, --gamma_summary_file,
#                                --model_choice_file,
#                                --visualization_output_dir, --quantile)
# Note    : Run after phynest_postprocess.jl has populated phynest_summary/.
# ============================================================================

using CSV
using DataFrames
using Statistics
using Distributions
using ArgParse

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir"
            help = "Directory containing PhyNEST summary files"
            default = "phynest_summary"
        "--output_file"
            help = "Output CSV file path"
            default = "results/PhyNEST_summary.csv"
        "--threshold_file"
            help = "Output CSV with the calibrated thresholds T*"
            default = "results/PhyNEST_T_threshold.csv"
        "--factor_summary_file"
            help = "Output CSV with summary statistics of T by factor level"
            default = "results/PhyNEST_T_summary_by_factor.csv"
        "--gamma_summary_file"
            help = "Output CSV with summary statistics of gamma_2 by level"
            default = "results/PhyNEST_gamma_summary_by_factor.csv"
        "--model_choice_file"
            help = "Output CSV with the percentage of h=1 by factor level"
            default = "results/PhyNEST_model_choice_by_factor.csv"
        "--visualization_output_dir"
            help = "Output directory for visualization files"
            default = "visualization_results/phynest"
        "--quantile"
            help = "Quantile of the pooled distribution of T used as T*"
            arg_type = Float64
            default = 0.95
    end
    return parse_args(s)
end

"""
Extract the parameter setting from a file name such as
PhyNEST-DUP0.0-LOS0.0-RVN-N_ind1-SF0.5-genelen1000-summary.csv
"""
function parse_parameter_setting(filename::String)
    m = match(r"DUP([\d.e-]+)-LOS([\d.e-]+)-RV([A-Z]+)-N_ind(\d+)" *
              r"-SF([\d.]+)-genelen(\d+)", filename)
    m === nothing && return nothing
    return (parameter_setting = m.match,
            dup_rate = parse(Float64, m[1]),
            loss_rate = parse(Float64, m[2]),
            ratevar = String(m[3]),
            n_inds = parse(Int, m[4]),
            SF = parse(Float64, m[5]))
end

# null setting: no duplication or loss, and no rate variation
is_null_setting(p) = p.dup_rate == 0 && p.loss_rate == 0 && p.ratevar == "N"
# settings pooled to calibrate T*: everything without lineage rate variation
is_pooled_setting(p) = p.ratevar in ("N", "G")

# median ignoring NaN (e.g. rank of the best run with the true tree when no
# run had it)
nanmedian(x) = (v = filter(!isnan, x); isempty(v) ? NaN : median(v))

"""
Clopper-Pearson 95% confidence interval (in %) for k successes out of n.
"""
function clopper_pearson(k::Int, n::Int)
    lo = k == 0 ? 0.0 : quantile(Beta(k, n - k + 1), 0.025)
    hi = k == n ? 1.0 : quantile(Beta(k + 1, n - k), 0.975)
    return (100 * lo, 100 * hi)
end

"""
Calibrate T*: pool the per-site statistic T over the settings without
lineage rate variation (no and gene-specific rates, all duplication/loss
rates), one threshold per (SF, n_inds) cell. The per-setting quantiles,
the quantile of the null setting alone (no duplication/loss, no rate
variation) and the quantile of the unscaled T are kept for reference.
"""
function calibrate_thresholds(settings, q::Float64)
    pooled = [s for s in settings if is_pooled_setting(s.params)]
    isempty(pooled) && error("No setting without lineage rate variation found")
    cells = sort(unique([(s.params.n_inds, s.params.SF) for s in pooled]))

    rows = []
    for (n_inds, SF) in cells
        cell = [s for s in pooled
                if s.params.SF == SF && s.params.n_inds == n_inds]
        T = vcat([s.df.T_per_site for s in cell]...)
        setting_quantiles = [quantile(s.df.T_per_site, q) for s in cell]
        null = [s for s in cell if is_null_setting(s.params)]
        push!(rows, (cell = "SF$SF-N_ind$n_inds", SF = SF, n_inds = n_inds,
                     n_settings = length(cell), n_reps = length(T),
                     quantile = q,
                     T_star = quantile(T, q),
                     median_T_per_site = median(T),
                     min_setting_quantile = minimum(setting_quantiles),
                     max_setting_quantile = maximum(setting_quantiles),
                     T_star_null_only = isempty(null) ? NaN :
                                        quantile(null[1].df.T_per_site, q),
                     T_star_raw = quantile(vcat([s.df.T for s in cell]...), q)))
    end
    return DataFrame(rows)
end

"""
Levels of the simulation factors over which replicates are pooled in the
supplementary tables: (factor name, level name, predicate on a setting).
"""
function factor_levels(settings)
    groups = [("overall", "all", s -> true),
              ("rate variation", "across genes", s -> s.params.ratevar == "G"),
              ("rate variation", "across lineages",
               s -> s.params.ratevar == "L"),
              ("rate variation", "none", s -> s.params.ratevar == "N")]
    for d in sort(unique([s.params.dup_rate for s in settings]))
        push!(groups, ("duplication/loss rate", string(d),
                       s -> s.params.dup_rate == d))
    end
    for (level, SF) in (("low", 0.5), ("high", 1.0))
        push!(groups, ("ILS level", level, s -> s.params.SF == SF))
    end
    return groups
end

# values of `column` pooled over the settings selected by `keep`
pooled(settings, keep, column) =
    vcat([s.df[!, column] for s in settings if keep(s)]...)

"""
Summary statistics of T pooled within each level of each simulation factor
(rate variation, duplication/loss rate, ILS level), as reported for the worst
residual of find_graphs in the paper's supplement.
"""
function summarize_T_by_factor(settings)
    stats(x) = (n = length(x), mean = mean(x), median = median(x), sd = std(x),
                q95 = quantile(x, 0.95), q99 = quantile(x, 0.99))
    rows = [merge((factor = f, level = l),
                  stats(pooled(settings, keep, :T_per_site)))
            for (f, l, keep) in factor_levels(settings)]
    return DataFrame(rows)
end

"""
Summary statistics of the minor inheritance probability gamma_2 pooled
within each factor level, with the percentage of replicates below each
threshold, as in the supplementary table for find_graphs and SNaQ.
"""
function summarize_gamma_by_factor(settings; thresholds = [0.05, 0.1, 0.25])
    below(x) = NamedTuple{Tuple(Symbol("pct_below_$t") for t in thresholds)}(
        Tuple(100 * count(<(t), x) / length(x) for t in thresholds))
    stats(x) = merge((n = length(x), mean = mean(x), median = median(x),
                      sd = std(x)), below(x))
    rows = [merge((factor = f, level = l),
                  stats(pooled(settings, keep, :gamma_2)))
            for (f, l, keep) in factor_levels(settings)]
    return DataFrame(rows)
end

"""
Percentage of replicates selecting h=1 under the calibrated rule, averaged
over the settings of each factor level (`pct_h1` maps a parameter setting
to its percentage), as in the supplementary table of model choice.
"""
function summarize_model_choice_by_factor(settings, pct_h1::Dict)
    rows = [(factor = f, level = l,
             n_settings = count(keep, settings),
             pct_h1 = mean(pct_h1[s.params.parameter_setting]
                           for s in settings if keep(s)))
            for (f, l, keep) in factor_levels(settings)]
    return DataFrame(rows)
end

"""
Summarize one parameter setting given its matched threshold T*.
Column names follow results/SNaQ_summary.csv where the quantity is the same.
"""
function summarize_setting(params, df::DataFrame, T_star::Float64)
    n = nrow(df)
    h1_cal = count(t -> t > T_star, df.T_per_site)
    h1_naive = count(t -> t > 0, df.T)
    ci_low, ci_high = clopper_pearson(h1_cal, n)

    true0 = df.RF_net0_true .== 0.0
    true1_major = df.RF_net1_1_true .== 0.0
    true1 = true1_major .| (df.RF_net1_2_true .== 0.0)
    # alternative placements of F, counted only when the true tree is missed
    alt_net0 = .!true0 .& ((df.RF_net0_alter1 .== 0.0) .|
        (df.RF_net0_alter2 .== 0.0) .| (df.RF_net0_alter3 .== 0.0))
    alt1(col) = sum((df[!, col] .== 0.0) .& .!true1)

    return (
        parameter_setting = params.parameter_setting,
        n_reps = n,
        H_eq_0 = n - h1_cal,
        H_eq_1 = h1_cal,
        H_gt_1 = 0, # only h=0 and h=1 are fitted
        H_eq_1_CI_low = ci_low,
        H_eq_1_CI_high = ci_high,
        H_eq_0_naive = n - h1_naive,
        H_eq_1_naive = h1_naive,
        T_star = T_star,
        mean_n_sites = mean(df.n_sites),
        min_n_sites = minimum(df.n_sites),
        median_T_per_site = median(df.T_per_site),
        q95_T_per_site = quantile(df.T_per_site, 0.95),
        mean_T = mean(df.T),
        median_T = median(df.T),
        q95_T = quantile(df.T, 0.95),
        min_T = minimum(df.T),
        mean_score_H0 = mean(df.score_H0),
        mean_score_H1 = mean(df.score_H1),
        mean_score_truetree = mean(df.score_truetree),
        mean_gamma_1 = mean(df.gamma_1),
        mean_gamma_2 = mean(df.gamma_2),
        find_true_net0 = sum(true0),
        find_true_net0_noF = sum(df.RF_net0_true_noF .== 0.0),
        find_true_net1 = sum(true1),
        find_true_net1_noF = sum((df.RF_net1_1_true_noF .== 0.0) .|
                                 (df.RF_net1_2_true_noF .== 0.0)),
        find_true_net1_major = sum(true1_major),
        find_alter_net0 = sum(alt_net0),
        find_alter1_net1_major = alt1(:RF_net1_1_alter1),
        find_alter2_net1_major = alt1(:RF_net1_1_alter2),
        find_alter3_net1_major = alt1(:RF_net1_1_alter3),
        find_alter1_net1_minor = alt1(:RF_net1_2_alter1),
        find_alter2_net1_minor = alt1(:RF_net1_2_alter2),
        find_alter3_net1_minor = alt1(:RF_net1_2_alter3),
        median_runs_at_best_H0 = median(df.n_runs_at_best_H0),
        median_runs_at_best_H1 = median(df.n_runs_at_best_H1),
        # where the true tree sits among the runs of each search
        median_runs_true_H0 = median(df.n_runs_true_H0),
        median_pct_runs_better_than_truetree_H0 =
            median(df.pct_runs_better_than_truetree_H0),
        find_true_any_run_H1 = count(>(0), df.n_runs_true_H1),
        median_runs_true_H1 = median(df.n_runs_true_H1),
        median_rank_best_true_H1 = nanmedian(df.rank_best_true_H1),
        median_gap_best_true_H1 = nanmedian(df.gap_best_true_H1),
        median_pct_runs_better_than_truetree_H1 =
            median(df.pct_runs_better_than_truetree_H1),
    )
end

function main()
    args = parse_commandline()
    input_dir = args["input_dir"]
    output_file = args["output_file"]
    threshold_file = args["threshold_file"]
    factor_file = args["factor_summary_file"]
    gamma_file = args["gamma_summary_file"]
    model_choice_file = args["model_choice_file"]
    visualization_dir = args["visualization_output_dir"]
    q = args["quantile"]

    if !isdir(input_dir)
        println("Error: Directory '$input_dir' not found!")
        return
    end

    # Read every per-setting file
    settings = []
    for csv_file in sort(filter(x -> endswith(x, ".csv"), readdir(input_dir)))
        params = parse_parameter_setting(csv_file)
        if params === nothing
            println("Skipping $csv_file: cannot parse the parameter setting")
            continue
        end
        df = CSV.read(joinpath(input_dir, csv_file), DataFrame)
        push!(settings, (params = params, df = df))
        println("Read $csv_file: $(nrow(df)) replicates")
    end
    if isempty(settings)
        println("No CSV files found in '$input_dir' directory!")
        return
    end

    #-----------------------------------------------#
    #   Calibrate T* from the null settings
    #-----------------------------------------------#
    thresholds = calibrate_thresholds(settings, q)
    mkpath(dirname(threshold_file))
    CSV.write(threshold_file, thresholds)
    println("\nCalibrated thresholds (quantile $q of T pooled over the " *
            "settings without lineage rate variation):")
    println(thresholds)
    println("Thresholds written to: $threshold_file")

    # threshold of the setting's own ILS level and number of individuals
    function matched_T_star(params)
        row = findfirst((thresholds.SF .== params.SF) .&
                        (thresholds.n_inds .== params.n_inds))
        row === nothing && error("No threshold for SF=$(params.SF), " *
                                 "n_inds=$(params.n_inds)")
        return thresholds.T_star[row]
    end

    #-----------------------------------------------#
    #   Summarize each setting
    #-----------------------------------------------#
    results = [summarize_setting(s.params, s.df, matched_T_star(s.params))
               for s in settings]
    summary_df = DataFrame(results)
    rename!(summary_df,
        :H_eq_0 => Symbol("H=0Accepted"),
        :H_eq_1 => Symbol("H=1Accepted"),
        :H_gt_1 => Symbol("H>1Accepted"),
        :H_eq_1_CI_low => Symbol("H=1Accepted_CI_low"),
        :H_eq_1_CI_high => Symbol("H=1Accepted_CI_high"),
        :H_eq_0_naive => Symbol("H=0Accepted_naive"),
        :H_eq_1_naive => Symbol("H=1Accepted_naive"))

    mkpath(dirname(output_file))
    CSV.write(output_file, summary_df)
    println("\nSummary completed!")
    println("Results written to: $output_file")

    factor_df = summarize_T_by_factor(settings)
    CSV.write(factor_file, factor_df)
    println("T by factor level written to: $factor_file")
    println(factor_df)

    gamma_df = summarize_gamma_by_factor(settings)
    CSV.write(gamma_file, gamma_df)
    println("gamma_2 by factor level written to: $gamma_file")
    println(gamma_df)

    pct_h1 = Dict(r.parameter_setting => 100 * r.H_eq_1 / r.n_reps
                  for r in results)
    choice_df = summarize_model_choice_by_factor(settings, pct_h1)
    CSV.write(model_choice_file, choice_df)
    println("% h=1 by factor level written to: $model_choice_file")
    println(choice_df)
    println(summary_df[:, [:parameter_setting, Symbol("H=1Accepted"),
                          Symbol("H=1Accepted_naive"), :T_star,
                          :find_true_net0, :find_true_net1,
                          :find_true_any_run_H1, :median_rank_best_true_H1,
                          :median_pct_runs_better_than_truetree_H1]])

    #-----------------------------------------------#
    #   Diagnostic plots (R)
    #-----------------------------------------------#
    mkpath(visualization_dir)
    low = thresholds.T_star[findfirst(thresholds.cell .== "SF0.5-N_ind1")]
    high = thresholds.T_star[findfirst(thresholds.cell .== "SF1.0-N_ind1")]

    println("\nGenerating T percentile plot with T* lines...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R');
                 plot_T_percentiles_jitter('$input_dir',
                     '$visualization_dir', 'PhyNEST_T_95-percentiles_jitter',
                     t_star = c(low = $low, high = $high),
                     facet_label = 'PhyNEST')"`)
    catch e
        @warn "Could not generate T percentile plot: $e"
    end

    println("Generating T distributions by setting...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R');
                 plot_T_distributions('$input_dir',
                     '$visualization_dir', 'PhyNEST_T_distributions',
                     t_star = c(low = $low, high = $high))"`)
    catch e
        @warn "Could not generate T distribution plot: $e"
    end

    println("Generating minor gamma plots...")
    try
        run(`Rscript -e "source('scripts/visual_utilities.R');
                 plot_overlapping_ratevar_by_n_inds_sf('$input_dir',
                     'gamma_1', 'gamma_2', '$visualization_dir',
                     'PhyNEST_minorG', 70)"`)
        run(`Rscript -e "source('scripts/visual_utilities.R');
                 plot_snaq_minor_gamma_by_tree_display('$input_dir',
                     '$visualization_dir',
                     'PhyNEST_minor_gamma_by_tree_display',
                     method_label = 'PhyNEST', split_by_placement = TRUE)"`)
    catch e
        @warn "Could not generate minor gamma plots: $e"
    end

    println("Finished generating visualizations in: $visualization_dir")
end

# Run the main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
