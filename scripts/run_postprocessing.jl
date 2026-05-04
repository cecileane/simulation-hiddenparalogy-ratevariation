#!/usr/bin/env julia
#=
Coordinator script to run postprocessing across all parameter sets in output/.

Usage:
    julia -p N scripts/run_postprocessing.jl \
        --mode <simulation|snaq|findgraphs> [options]

Options:
    --mode          Postprocessing mode: simulation, snaq, or findgraphs
    --n_reps        Number of replicates per parameter set (default: 100)
    --output_dir    Directory with parameter folders (default: "output")
    --saved_path    Directory to copy summary CSVs into after each run
                    (default: <project_root>/<mode>_summary)
    --n_procs       Workers to pass to each sub-script (default: nworkers())

Example:
    julia -p 10 scripts/run_postprocessing.jl --mode snaq
    julia -p 4  scripts/run_postprocessing.jl --mode findgraphs --n_reps 10
    julia       scripts/run_postprocessing.jl --mode simulation --n_procs 8
=#

using Distributed
using ArgParse

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

function parse_commandline()
    s = ArgParseSettings(
        description = "Run postprocessing for every parameter set in output/."
    )

    @add_arg_table s begin
        "--mode"
            help    = "Postprocessing mode: simulation, snaq, or findgraphs. " *
                      "snaq and findgraphs also copy consensus results to " *
                      "<project_root>/snaq_consensus_tree or " *
                      "findgraphs_consensus_tree."
            arg_type = String
            required = true

        "--n_reps"
            help    = "Number of replicates per parameter set"
            arg_type = Int
            default = 100

        "--output_dir"
            help    = "Directory with output parameter folders"
            arg_type = String
            default = "output"

        "--saved_path"
            help    = "Directory to copy the summary CSV file(s) into. " *
                      "Defaults to <project_root>/<mode>_summary."
            arg_type = String
            default = ""

        "--n_procs"
            help    = "Workers to pass via -p to each sub-script. " *
                      "Defaults to nworkers() (-p N given to this script)."
            arg_type = Int
            default = -1   # sentinel → auto-detect
    end

    return parse_args(s)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
Parse a parameter-folder name like
    DUP0.0-LOS0.0-RVG-N_ind1-SF0.5-genelen1000
into a NamedTuple.  Returns `nothing` when the name does not match.
"""
function parse_paramname(dirname::AbstractString)
    m = match(
        r"^DUP([\d.]+)-LOS([\d.]+)-RV([A-Z]+)-N_ind(\d+)" *
        r"-SF([\d.]+)-genelen(\d+)$",
        dirname
    )
    m === nothing && return nothing
    return (
        dup_rate  = parse(Float64, m[1]),
        loss_rate = parse(Float64, m[2]),
        ratevar   = String(m[3]),
        n_inds    = parse(Int,     m[4]),
        SF        = parse(Float64, m[5]),
        gene_len  = parse(Int,     m[6]),
    )
end

"""
Return the path to the summary CSV that the given postprocessing
script is expected to produce for `paramname`.
"""
function summary_file_path(output_dir::String, paramname::String, mode::String)
    if mode == "simulation"
        return joinpath(output_dir, paramname, "summary_$(paramname).csv")
    elseif mode == "snaq"
        return joinpath(output_dir, paramname, "SNaQ-$(paramname)-summary.csv")
    elseif mode == "findgraphs"
        return joinpath(output_dir, paramname, "findgraph-$(paramname).csv")
    else
        error("Unknown mode: $mode")
    end
end

"""
Build the shell command that runs the appropriate sub-script.
"""
function build_cmd(mode::String, script_dir::String, params,
    n_reps::Int, n_procs::Int)
    script = joinpath(script_dir, "$(mode)_postprocess.jl")

    julia_flags = n_procs > 0 ? ["-p", string(n_procs)] : String[]

    args = [
        "--dup_rate",  string(params.dup_rate),
        "--loss_rate", string(params.loss_rate),
        "--ratevar",   params.ratevar,
        "--n_reps",    string(n_reps),
        "--n_inds",    string(params.n_inds),
        "--SF",        string(params.SF),
        "--gene_len",  string(params.gene_len),
    ]

    return Cmd(vcat(["julia"], julia_flags, [script], args))
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    parsed_args = parse_commandline()

    mode       = parsed_args["mode"]
    n_reps     = parsed_args["n_reps"]
    output_dir = parsed_args["output_dir"]
    saved_path = parsed_args["saved_path"]
    n_procs    = parsed_args["n_procs"]

    # Validate mode
    valid_modes = ("simulation", "snaq", "findgraphs")
    if !(mode in valid_modes)
        error("Invalid --mode '$mode'. Options: $(join(valid_modes, ", "))")
    end

    base_mode = mode

    # Resolve worker count: prefer explicit --n_procs, then fall back to
    # the workers already spawned via -p N on this process.
    if n_procs == -1
        n_procs = nworkers()   # equals N when launched as  julia -p N ...
        # nworkers()==1 means no -p flag; pass 0 (no -p to child)
        n_procs = n_procs == 1 ? 0 : n_procs
    end

    # Determine project root (one level above this script)
    script_dir   = dirname(abspath(PROGRAM_FILE))
    project_root = normpath(joinpath(script_dir, ".."))

    # Resolve output_dir relative to project root if not absolute
    if !isabspath(output_dir)
        output_dir = joinpath(project_root, output_dir)
    end

    # Resolve saved_path
    if isempty(saved_path)
        saved_path = joinpath(project_root, "$(base_mode)_summary")
    elseif !isabspath(saved_path)
        saved_path = joinpath(project_root, saved_path)
    end

    # Ensure saved_path exists
    mkpath(saved_path)

    # Consensus collection dir (snaq and findgraphs collect consensus results)
    consensus_collect_dir = (mode in ("snaq", "findgraphs")) ?
        joinpath(project_root, "$(mode)_consensus_tree") : ""
    if !isempty(consensus_collect_dir)
        mkpath(consensus_collect_dir)
    end

    # Collect all matching parameter-set directories
    if !isdir(output_dir)
        error("output_dir does not exist: $output_dir")
    end

    param_dirs = sort(filter(readdir(output_dir)) do d
        isdir(joinpath(output_dir, d)) && parse_paramname(d) !== nothing
    end)

    if isempty(param_dirs)
        @warn "No parameter-set directories found in $output_dir"
        return
    end

    println("Found $(length(param_dirs)) parameter set(s):")
    for d in param_dirs
        println("  $d")
    end
    println()

    # -----------------------------------------------------------------------
    # Run the postprocessing script for each parameter set
    # -----------------------------------------------------------------------
    successes             = String[]
    failures              = String[]
    consensus_plots_copied = String[]   # populated for snaq / findgraphs modes

    for paramname in param_dirs
        params = parse_paramname(paramname)   # never nothing here

        println("-" ^ 60)
        println("Processing: $paramname")

        cmd = build_cmd(base_mode, script_dir, params, n_reps, n_procs)
        println("Running: $cmd\n")

        # Run from the project root so relative paths inside the sub-script
        # (e.g. "output/...") resolve correctly.
        success = true
        try
            run(setenv(cmd, dir = project_root))
        catch e
            @warn "Sub-script failed for $paramname: $e"
            success = false
        end

        if !success
            push!(failures, paramname)
            println("FAILED: $paramname\n")
            continue
        end

        # -------------------------------------------------------------------
        # Copy summary CSV to saved_path (force-overwrite)
        # -------------------------------------------------------------------
        src = summary_file_path(output_dir, paramname, base_mode)

        if isfile(src)
            dst = joinpath(saved_path, basename(src))
            cp(src, dst; force = true)
            println("Copied summary → $dst")
            push!(successes, paramname)
        else
            @warn "Expected summary file not found after run: $src"
            push!(failures, paramname)
        end

        # -------------------------------------------------------------------
        # Copy consensus plot PDFs (flat) to central collection dir
        # Only the _plot.pdf files are copied; no subdirectory is created.
        # -------------------------------------------------------------------
        if !isempty(consensus_collect_dir)
            prefix      = base_mode == "snaq" ? "SNaQ" : "findgraphs"
            file_prefix = base_mode == "snaq" ? "SNaQ" : "findgraph"
            cons_folder = joinpath(output_dir, paramname,
                                   "$(prefix)-$(paramname)-consensus_nets")
            if isdir(cons_folder)
                for h in ("H0", "H1")
                    pdf_name = "$(file_prefix)-$(paramname)" *
                        "-consensus_$(h)_plot.pdf"
                    src_pdf  = joinpath(cons_folder, pdf_name)
                    if isfile(src_pdf)
                        dst_pdf = joinpath(consensus_collect_dir, pdf_name)
                        cp(src_pdf, dst_pdf; force = true)
                        push!(consensus_plots_copied, dst_pdf)
                        println("Copied $(h) plot → $dst_pdf")
                    else
                        @warn "Consensus plot not found (skipping): $src_pdf"
                    end
                end
            else
                @warn "Consensus folder not found (skipping): $cons_folder"
            end
        end

        println()
    end

    # -----------------------------------------------------------------------
    # Final report
    # -----------------------------------------------------------------------
    println("=" ^ 60)
    println("POST-PROCESSING COMPLETE  (mode = $mode)")
    println("  Succeeded : $(length(successes)) / $(length(param_dirs))")
    println("  Failed    : $(length(failures))")
    if !isempty(failures)
        println("  Failed parameter sets:")
        for f in failures
            println("    $f")
        end
    end
    println("  Summary files saved to: $saved_path")
    if !isempty(consensus_collect_dir)
        println("  Consensus plot dir: $consensus_collect_dir")
        if isempty(consensus_plots_copied)
            println("    (no consensus plots were copied)")
        else
            println("  Consensus plots copied " *
                "($(length(consensus_plots_copied))):")
            for p in consensus_plots_copied
                println("    $p")
            end
        end
    end
    println("=" ^ 60)
end

main()
