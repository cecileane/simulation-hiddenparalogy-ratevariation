
# Runs one postprocessing script (<mode>_postprocess.jl) 
# on every parameter set in output/
# and copies each summary csv to <mode>_summary/ 
# (for snaq and findgraphs, also the consensus plots
# to <mode>_consensus_tree/)
# Usage: julia -p N scripts/run_postprocessing.jl 
# e.g. julia -p 10 scripts/run_postprocessing.jl --mode snaq
#      julia scripts/run_postprocessing.jl --mode phynest 
#      --output_dir Phynest_results

using Distributed
using ArgParse

#--------------- argument parsing ---------------#
function parse_commandline()
    s = ArgParseSettings(
        description = "Run postprocessing for every parameter set in output/."
    )

    @add_arg_table s begin # see ArgParse package tutorial 
        "--mode"
            help    = "Postprocessing mode: simulation, snaq, findgraphs " *
                      "or phynest. " *
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
            default = -1   # sentinel -> auto-detect

        "--tree"
            help    = "Species tree of the runs: reptile or plant " *
                      "(passed to every sub-script; plant results are kept " *
                      "apart with _plant suffixes on the default folders)"
            arg_type = String
            default = "reptile"
    end

    return parse_args(s)
end

#--------------- functions ---------------#
"""
Input: dirname (String): name of a parameter folder, 
e.g. DUP0.0-LOS0.0-RVG-N_ind1-SF0.5-genelen1000
Output: NamedTuple (dup_rate, loss_rate, ratevar, n_inds, SF, gene_len)
nothing if the name does not match the pattern 
    (then it is not a parameter folder)
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
Input: output_dir, paramname, mode (Strings)
Output: path of the summary csv that <mode>_postprocess.jl 
    writes for this parameter set
(each postprocess script names its csv differently)
"""
function summary_file_path(output_dir::String, paramname::String, mode::String)
    if mode == "simulation"
        return joinpath(output_dir, paramname, "summary_$(paramname).csv")
    elseif mode == "snaq"
        return joinpath(output_dir, paramname, "SNaQ-$(paramname)-summary.csv")
    elseif mode == "findgraphs"
        return joinpath(output_dir, paramname, "findgraph-$(paramname).csv")
    elseif mode == "phynest"
        return joinpath(output_dir, paramname,
                        "PhyNEST-$(paramname)-summary.csv")
    else
        error("Unknown mode: $mode")
    end
end

"""
Input: mode, script_dir, params (from parse_paramname), 
n_reps, n_procs, output_dir, tree
Output: Cmd to run <mode>_postprocess.jl for this parameter set, 
 with -p n_procs if n_procs > 0
"""
function build_cmd(mode::String, script_dir::String, params,
    n_reps::Int, n_procs::Int, output_dir::String, tree::String)
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
    # only phynest_postprocess.jl takes --output_dir, 
    # the others read output/ or output_plant/
    if mode == "phynest"
        append!(args, ["--output_dir", output_dir])
    end
    # all of them need --tree: picks the reference trees and the output folder
    append!(args, ["--tree", tree])

    return Cmd(vcat(["julia"], julia_flags, [script], args))
end

#--------------- main function here --------------#
function main()
    parsed_args = parse_commandline()

    mode       = parsed_args["mode"]
    n_reps     = parsed_args["n_reps"]
    output_dir = parsed_args["output_dir"]
    saved_path = parsed_args["saved_path"]
    n_procs    = parsed_args["n_procs"]
    tree       = parsed_args["tree"]

    # check mode
    valid_modes = ("simulation", "snaq", "findgraphs", "phynest")
    if !(mode in valid_modes)
        error("Invalid --mode '$mode'. Options: $(join(valid_modes, ", "))")
    end

    base_mode = mode

    # workers: --n_procs if given, otherwise the -p N of this script
    if n_procs == -1
        n_procs = nworkers()   # equals N when launched as  julia -p N ...
        # nworkers() == 1 means no -p flag, then pass 0 (no -p for the sub-script)
        n_procs = n_procs == 1 ? 0 : n_procs
    end

    # project root = one level above scripts/
    script_dir   = dirname(abspath(PROGRAM_FILE))
    project_root = normpath(joinpath(script_dir, ".."))

    # plant runs: output_plant/ and a _plant suffix on the summary folder, unless given
    suffix = (tree == "reptile" ? "" : "_$tree")
    if tree != "reptile" && output_dir == "output"
        output_dir = "output$suffix"
    end

    # relative paths are relative to the project root
    if !isabspath(output_dir)
        output_dir = joinpath(project_root, output_dir)
    end

    # same for saved_path
    if isempty(saved_path)
        saved_path = joinpath(project_root, "$(base_mode)_summary$suffix")
    elseif !isabspath(saved_path)
        saved_path = joinpath(project_root, saved_path)
    end

    # make sure saved_path exists
    mkpath(saved_path)

    # snaq and findgraphs also have consensus plots to collect
    consensus_collect_dir = (mode in ("snaq", "findgraphs")) ?
        joinpath(project_root, "$(mode)_consensus_tree$suffix") : ""
    if !isempty(consensus_collect_dir)
        mkpath(consensus_collect_dir)
    end

    # all parameter folders in output_dir
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

    # run the postprocess script for each parameter set
    successes = String[]
    failures = String[]
    consensus_plots_copied = String[] # populated for snaq / findgraphs modes

    for paramname in param_dirs
        params = parse_paramname(paramname)# never nothing here

        println("-" ^ 60)
        println("Processing: $paramname")

        cmd = build_cmd(base_mode, script_dir, params, n_reps, n_procs,
                        output_dir, tree)
        println("Running: $cmd\n")

        # run from the project root, so that relative paths (output/...) in the sub-script work
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

        # copy the summary csv to saved_path, overwrite if there
        src = summary_file_path(output_dir, paramname, base_mode)

        if isfile(src)
            dst = joinpath(saved_path, basename(src))
            cp(src, dst; force = true)
            println("Copied summary -> $dst")
            push!(successes, paramname)
        else
            @warn "Expected summary file not found after run: $src"
            push!(failures, paramname)
        end

        # copy the consensus plot pdfs (only the _plot.pdf files, no subfolder) to one folder
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
                        println("Copied $(h) plot -> $dst_pdf")
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

    # final report
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
