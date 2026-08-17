#= rerun_gof.jl
This script reloads net0, net1, and CFs from existing SNaQ outputs
and reruns quarnetGoFtest! for all replicates of a given parameter setting.
It takes the same parameter arguments as snaq.jl.
The purpose of this script is to rerun GoF tests with the original seeds.
This saves as a quality assurance check for debugging and 
    was used during the development of SNaQ.jl sript. 
This script needs to be kept so that in the future, it 
    is easier to rerun GoF tests with the original seeds if needed.  
=# 

using ArgParse
using Distributed
using TimerOutputs
@everywhere using Dates
@everywhere using TimeZones
@everywhere using DataFrames
@everywhere using CSV 
@everywhere using Printf
@everywhere using PhyloNetworks
@everywhere using SNaQ
@everywhere using QuartetNetworkGoodnessFit
@everywhere include("utilities.jl")

const to = TimerOutput()
tz = TimeZone("America/Chicago")
current_time_tz = ZonedDateTime(now(), tz)
time_str = Dates.format(current_time_tz, "yyyy-mm-dd HH:MM:SS zzz")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--dup_rate"
            help = "Parameter setting (duplication rate)"
            arg_type = Float64
            required = true
        "--loss_rate"
            help = "Parameter setting (gene loss rate)"
            arg_type = Float64
            required = true
        "--ratevar"
            help = "Parameter setting (variation rate)"
            arg_type = String
            required = true
        "--n_reps"
            help = "How many reps this parameter set was run in total"
            arg_type = Int
            default = 100
        "--rep_start"
            help = "Start of the range of replicates to rerun GoF"
            arg_type = Int
            default = 1
        "--rep_end"
            help = "End of the range of replicates to rerun GoF"
            arg_type = Int
            default = -1
        "--n_inds"
            help = "Number of individuals per species"
            arg_type = Int
            default = 1
        "--SF"
            help = "Scaling factor for effective population Ne (Default = 1.0)"
            arg_type = Float64
            default = 1.0
        "--gene_len"
            help = "Length of simulated gene sequences (Default = 1000 bp)"
            arg_type = Int
            default = 1000

    end

    parsed_args = parse_args(s)
    if parsed_args["rep_end"] == -1
        parsed_args["rep_end"] = parsed_args["n_reps"]
    end
    return parsed_args
end

parsed_args = parse_commandline()

dup_rate  = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar   = parsed_args["ratevar"]
n_reps    = parsed_args["n_reps"]
n_inds    = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end   = parsed_args["rep_end"]
SF        = parsed_args["SF"]
gene_len  = parsed_args["gene_len"]

# Reconstruct folder names (same logic as snaq.jl)
paramname_root = set_up_paramname_root(dup_rate, 
            loss_rate, ratevar, n_inds, 
            SF, gene_len)
outfolder = "output/$paramname_root"

index_length = rep_end - rep_start + 1

#-----------------------------------------------#
#    Reconstruct seed array (identical to snaq.jl)
#-----------------------------------------------#
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
master_seed = generate_master_seed(params_dict_for_seed_setting)

software_names = ["snaq"]
seed_dic = generate_software_seeds(master_seed, software_names)
seed_snaq = seed_dic["snaq"]

# n_reps x 6: columns 3 and 4 are seed_qgof0 and seed_qgof1
seed_array = seed_generator(seed_snaq, n_reps, 6, 
                    outfolder, "random_seed_snaq.txt")

#-----------------------------------------------#
#    Build folder path lists
#-----------------------------------------------#
folder_path_list = []
for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
    push!(folder_path_list, rep_folder_path)
end

snaqfolder_list = []
for ind in 1:index_length
    snaqfolder = setup_rep_output_folders(folder_path_list, 
                                        ind, "snaqfolder")
    push!(snaqfolder_list, snaqfolder)
end

#-----------------------------------------------#
#    Broadcast shared variables to workers
#-----------------------------------------------#
@everywhere seed_array  = $seed_array
@everywhere snaqfolder_list = $snaqfolder_list
@everywhere folder_path_list = $folder_path_list
@everywhere outfolder = $outfolder
@everywhere rep_start = $rep_start

@everywhere begin
"""
    rerun_gof_for_replicate(ind::Int)

Load net0, net1, and CFs from existing SNaQ output 
    for replicate `ind` and rerun quarnetGoFtest! 
    using the original seeds.
Overwrites the existing snaq_gof_results_H0.csv 
    and snaq_gof_results_H1.csv files.
"""
function rerun_gof_for_replicate(ind::Int)

    snaqfolder = snaqfolder_list[ind]
    H0folder   = joinpath(snaqfolder, "H0_output")
    H1folder   = joinpath(snaqfolder, "H1_output")

    # --- Load networks ---
    net0_path = joinpath(H0folder, "H0.out")
    net1_path = joinpath(H1folder, "H1.out")

    if !isfile(net0_path)
        println("Worker $(myid()): MISSING $net0_path — skipping rep $ind")
        return
    end
    if !isfile(net1_path)
        println("Worker $(myid()): MISSING $net1_path — skipping rep $ind")
        return
    end

    net0 = readnewick(net0_path)
    net1 = readnewick(net1_path)

    # --- Extract gamma values from net1 ---
    gammas = [e.gamma for e in net1.edge if e.hybrid]
    println("Worker $(myid()): rep $ind — " *
        "$(length(gammas)) hybrid edge(s), gammas=$gammas")

    if length(gammas) >= 2
        gamma_1, gamma_2 = gammas[1], gammas[2]
        if gamma_1 < gamma_2
            gamma_1, gamma_2 = gamma_2, gamma_1
        end
    else
        println("Worker $(myid()): WARNING — <2 hybrid edges for rep $ind;" *
            " gammas will be missing.")
        gamma_1 = length(gammas) >= 1 ? gammas[1] : NaN
        gamma_2 = NaN
    end

    # --- Load CFs ---
    cffile = joinpath(snaqfolder, "CF_results.csv")
    if !isfile(cffile)
        println("Worker $(myid()): MISSING $cffile — skipping replicate $ind")
        return
    end
    qCF = CSV.read(cffile, DataFrame)

    # --- Seeds (same columns as in snaq.jl) ---
    # seed_generator returns Int32; quarnetGoFtest! requires Int64, 
    # so convert explicitly.
    ind_in_seed_array = ind + rep_start - 1
    seed_qgof0 = Int64(seed_array[ind_in_seed_array, 3])
    seed_qgof1 = Int64(seed_array[ind_in_seed_array, 4])

    # ---- Paths for output ----
    gof0_path = joinpath(H0folder, "snaq_gof_results_H0.csv")
    gof1_path = joinpath(H1folder, "snaq_gof_results_H1.csv")

    # --- Check for zero-gene-tree quartets ---
    zero_gene_quartets = filter(row -> row.ngenes == 0, qCF)

    if nrow(zero_gene_quartets) > 0
        println("\n" * "="^70)
        println("Worker $(myid()): WARNING — " *
            "$(nrow(zero_gene_quartets)) quartet(s) with ngenes=0 in rep $ind")
        println("Skipping quarnetGoFtest! for this replicate.")
        for row in eachrow(zero_gene_quartets)
            println("  $(row.t1), $(row.t2), $(row.t3), $(row.t4)")
        end
        println("="^70 * "\n")

        score_net0 = loglik(net0)
        score_net1 = loglik(net1)

        df_0 = DataFrame(type  = ["score", "status"],
            value = [score_net0, "SKIPPED: missing quartets"])
        CSV.write(gof0_path, df_0)

        df_1 = DataFrame(type  = ["score", "gamma_1", "gamma_2", "status"],
            value = [score_net1, gamma_1, gamma_2, "SKIPPED: missing quartets"])
        CSV.write(gof1_path, df_1)

        println("Worker $(myid()): Saved minimal results (no GoF) to " *
            "$gof0_path and $gof1_path")
        return
    end

    # --- Run quarnetGoFtest! (same two-pass pipeline as snaq_1rep.jl) ---
    score_net0 = loglik(net0)
    score_net1 = loglik(net1)

    println("Worker $(myid()): rep $ind — GoFtest! H0 warm-up")
    res0 = quarnetGoFtest!(net0, qCF, true; seed=201, nsim=10)
    println("Worker $(myid()): rep $ind — GoFtest! H1 warm-up")
    res1 = quarnetGoFtest!(net1, qCF, true; seed=202, nsim=50)

    net0 = res0[5]
    net1 = res1[5]

    println("Worker $(myid()): rep $ind — GoFtest! H0 (seed=$seed_qgof0)")
    res0 = quarnetGoFtest!(net0, qCF, false; seed=seed_qgof0, nsim=1000)
    println("Worker $(myid()): rep $ind — GoFtest! H1 (seed=$seed_qgof1)")
    res1 = quarnetGoFtest!(net1, qCF, false; seed=seed_qgof1, nsim=1000)

    # --- Serialise results (same format as snaq_1rep.jl) ---
    result_names = ["p", "z_uncorrected", "sigma", "bootstrap_values",
                    "network", "z_bootstrap"]

    function serialise_res(res)
        rows = []
        for (i, x) in enumerate(res)
            name = i <= length(result_names) ? result_names[i] : "unknown_$i"
            if typeof(x) <: Number
                push!(rows, (type=name, value=x))
            elseif typeof(x) <: AbstractArray
                push!(rows, (type=name, value=join(x, ";")))
            elseif typeof(x) == HybridNetwork
                push!(rows, (type=name, value=string(x)))
            else
                push!(rows, (type=name, value=string(x)))
            end
        end
        return DataFrame(rows)
    end

    # Build new GoF-specific rows from the freshly computed results.
    # Convert value to String for consistent vcat with rows loaded from CSV.
    df_0_gof = serialise_res(res0)
    df_0_gof.value = string.(df_0_gof.value)

    df_1_gof = serialise_res(res1)
    df_1_gof.value = string.(df_1_gof.value)

    # Load existing CSVs; preserve non-GoF rows (score, gamma_*, status, etc.).
    # If the file does not yet exist, start with an empty frame.
    gof_row_names = Set(result_names)

    df_0_existing = isfile(gof0_path) ?
        CSV.read(gof0_path, DataFrame, types=Dict("value" => String)) :
        DataFrame(type=String[], value=String[])
    df_0_preserved = filter(row -> !(row.type in gof_row_names), df_0_existing)

    df_1_existing = isfile(gof1_path) ?
        CSV.read(gof1_path, DataFrame, types=Dict("value" => String)) :
        DataFrame(type=String[], value=String[])
    df_1_preserved = filter(row -> !(row.type in gof_row_names), df_1_existing)

    # Merge: new GoF rows first (matching snaq_1rep.jl column order),
    # then preserved non-GoF rows (score, gamma_1, gamma_2, …).
    df_0 = vcat(df_0_gof, df_0_preserved)
    df_1 = vcat(df_1_gof, df_1_preserved)

    println(df_0)
    println(df_1)

    CSV.write(gof0_path, df_0)
    CSV.write(gof1_path, df_1)

    println("Worker $(myid()): rep $ind — saved GoF to " *
        "$gof0_path and $gof1_path")
end
end # @everywhere begin

#-----------------------------------------------#
#    Run in parallel over replicates
#-----------------------------------------------#
@timeit to "Rerunning GoF from rep$rep_start to rep$rep_end" begin
    pmap(ind -> begin
        println("Worker $(myid()): Starting GoF rerun for replicate index $ind")
        rerun_gof_for_replicate(ind)
    end, 1:index_length)
end

#-----------------------------------------------#
#    Write timing log
#-----------------------------------------------#
host_name = gethostname()
processors = nprocs()

log_text = """
  #=====================================================#
  #-------------- Rerun quarnetGoFtest! ----------------#
  #=====================================================#
  paramname_root  = $paramname_root
  rep_start       = $rep_start
  rep_end         = $rep_end
  processors      = $processors
  server          = $host_name
  time            = $time_str
  --- Running time ---
  """

log_file = joinpath(outfolder, "arguments-rerun-gof-$paramname_root.log")
open(log_file, "a") do io
    println(io, log_text)
    show(io, to)
end

println("=============================================")
println("GoF rerun completed!")
println("Log saved to $log_file")
