# SNaQ postprocessing: re-analyzes failed reps after snaq.jl finishes.
# Usage: julia snaq_postprocessing.jl --dup_rate X --loss_rate Y
#          --ratevar Z --n_reps N [--rep_start S --rep_end E]

using ArgParse
using Distributed
@everywhere using Printf
@everywhere using PhyloNetworks
@everywhere using PhyloSummaries
@everywhere using CSV
@everywhere using DataFrames
@everywhere include("utilities.jl")
using PhyloPlots
using RCall

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        # Parameter settings to identify the correct output folder
        "--dup_rate"
            help = "Parameter setting (duplication rate) used in snaq"
            arg_type = Float64
            required = true
        "--loss_rate"
            help = "Parameter setting (gene loss rate) used in snaq"
            arg_type = Float64
            required = true
        "--ratevar"
            help = "Parameter setting (variation rate) used in snaq"
            arg_type = String
            required = true
        "--n_reps"
            help = "Total number of reps for this parameter set"
            arg_type = Int
            default = 100
        "--rep_start"
            help = "Start of the range of replicates to process"
            arg_type = Int
            default = 1
        "--rep_end"
            help = "End of the range of replicates to process"
            arg_type = Int
            default = -1
        "--n_inds"
            help = "Number of individuals per species"
            arg_type = Int
            default = 1
        "--SF"
            help = "Scaling factor to scale effective population Ne
                (Default = 1.0, no scaling)"
            arg_type = Float64
            default = 1.0
        "--gene_len"
            help = "Gene length used in simulation (default = 1000)"
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

# Parse arguments
dup_rate = parsed_args["dup_rate"]
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"]
n_reps = parsed_args["n_reps"]
n_inds = parsed_args["n_inds"]
rep_start = parsed_args["rep_start"]
rep_end = parsed_args["rep_end"]
SF = parsed_args["SF"] 
gene_len = parsed_args["gene_len"]  

# Set up folders
paramname_root = set_up_paramname_root(dup_rate, loss_rate, 
                                    ratevar, n_inds, SF,
                                    gene_len) 
outfolder = "output/$paramname_root"

# Check if the output folder exists
if !isdir(outfolder)
    error("Output folder does not exist: $outfolder")
end

println("SNaQ Postprocessing Script")
println("Parameter set: $paramname_root")
println("Processing replicates $rep_start to $rep_end")
println("Output folder: $outfolder")

#-----------------------------------------------#       
#     Check if summary file exists, create if missing
#-----------------------------------------------#

# Define the goodness-of-fit summary file path
gof_summary_file = joinpath(outfolder, "SNaQ-$paramname_root-summary.csv")

# Check if the summary file already exists
if !isfile(gof_summary_file)
    println("Goodness-of-fit summary file not found: $gof_summary_file")
    println("Creating goodness-of-fit summary from individual result files...")
    
    # This flag indicates we need to create the summary file
    create_summary = true
else
    println("Goodness-of-fit summary file exists: $gof_summary_file")
    println("Proceeding with RF recovery and updating existing summary...")
    
    # This flag indicates we're updating an existing file
    create_summary = false
end

#-----------------------------------------------#       
#   Calculate RF distances to true species tree
#-----------------------------------------------#

# True species tree used in calculations
true_species_tree = readnewick("(A,((((B,C),(D,E)),F),(G,H)));")

#=
 If the estimated tree under H=0 or displayed trees under H=1 are 
 not the true species tree (A,((G,H),(F,((B,C),(D,E))))); 
 
 Then, we want to understand if it could be the following alternatives 
 1. alter1 (A,((G,H),(((B,C),F),(D,E)))); F is clustered with (B,C) 
 2. alter2 (A,((G,H),((B,C),(F,(D,E))))); F is clustered with (D,E) 
 3. alter3 (A,(((G,H),((D,E),(B,C))),F)); F is outside of (B,C),(D,E),(G,H)
Therefore, not only we calculate RF distance to true species tree,
    but also to these three alternatives. This will help us understand 
    why the estimated tree is different from the true species tree. 
=# 
alter1 = readnewick("(A,((G,H),(((B,C),F),(D,E))));") # F with (B,C)
alter2 = readnewick("(A,((G,H),((B,C),(F,(D,E)))));") # F with (D,E)
alter3 = readnewick("(A,(((G,H),((D,E),(B,C))),F));") # F outside all
species_tree_noF = readnewick("(A,((G,H),((B,C),(D,E))));") # no F
species_tree_noG = readnewick("(A,((((B,C),(D,E)),F),H));") # no G
species_tree_noH = readnewick("(A,((((B,C),(D,E)),F),G));") # no H

"""
Calculate RF distances between estimated tree and alternative trees
"""
function check_alternative_trees(est_tree::PhyloNetworks.HybridNetwork)
    rf_alt1 = mudistance_semidirected(est_tree, alter1, preorder=true)
    rf_alt2 = mudistance_semidirected(est_tree, alter2, preorder=true)
    rf_alt3 = mudistance_semidirected(est_tree, alter3, preorder=true) 
    return (rf_alt1, rf_alt2, rf_alt3) 
end 

for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
    
    net0_file = joinpath(rep_folder_path, "snaqfolder", "H0_output", "H0.out")
    results0_file = joinpath(rep_folder_path,
        "snaqfolder", "H0_output", "snaq_gof_results_H0.csv")
    net1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", "H1.out")
    results1_file = joinpath(rep_folder_path,
        "snaqfolder", "H1_output", "snaq_gof_results_H1.csv")
    
    # Calculate RF distance for H0
    if isfile(net0_file) && isfile(results0_file)
        df0 = CSV.read(results0_file, DataFrame; types=Dict(:value => String))
        net0 = readnewick(net0_file) 
        rf0 = mudistance_semidirected(net0, true_species_tree, preorder=true) 

        # Check if estimated tree matches species tree with no F
        rf0_noF = NaN
        try
            net0_noF = deepcopy(net0)
            deleteleaf!(net0_noF, "F") 
            rf0_noF = mudistance_semidirected(
                net0_noF, species_tree_noF, preorder=true)
        catch e
            println(
                "Warning: RF_net0_noF failed for rep $rep_number_string: $e")
        end 

        (rf_alt1, rf_alt2, rf_alt3) = check_alternative_trees(net0) 
        println("Rep $rep_number_string: RF=$rf0 alt1=$rf_alt1 alt2=$rf_alt2" *
            " alt3=$rf_alt3")
        
        rf_row_idx = findfirst(row -> row.type == "RF_net0_true", eachrow(df0))
        rf_alt1_idx = findfirst(
            row -> row.type == "RF_net0_alter1", eachrow(df0))
        rf_alt2_idx = findfirst(
            row -> row.type == "RF_net0_alter2", eachrow(df0))
        rf_alt3_idx = findfirst(
            row -> row.type == "RF_net0_alter3", eachrow(df0))
        rf_true_noF_idx = findfirst(
            row -> row.type == "RF_net0_true_noF", eachrow(df0))

        if rf_row_idx === nothing || rf_alt1_idx === nothing ||
                rf_alt2_idx === nothing || rf_alt3_idx === nothing ||
                rf_true_noF_idx === nothing
            # Add any missing rows
            if rf_row_idx === nothing
                push!(df0, (type="RF_net0_true", value=string(rf0)))
            else
                df0[rf_row_idx, :value] = string(rf0)
            end
            if rf_alt1_idx === nothing
                push!(df0, (type="RF_net0_alter1", value=string(rf_alt1)))
            else
                df0[rf_alt1_idx, :value] = string(rf_alt1)
            end
            if rf_alt2_idx === nothing
                push!(df0, (type="RF_net0_alter2", value=string(rf_alt2)))
            else
                df0[rf_alt2_idx, :value] = string(rf_alt2)
            end
            if rf_alt3_idx === nothing
                push!(df0, (type="RF_net0_alter3", value=string(rf_alt3)))
            else
                df0[rf_alt3_idx, :value] = string(rf_alt3)
            end
            if rf_true_noF_idx === nothing
                push!(df0, (type="RF_net0_true_noF", value=string(rf0_noF)))
            else
                df0[rf_true_noF_idx, :value] = string(rf0_noF)
            end
        else
            # All rows exist, just update them
            df0[rf_row_idx, :value] = string(rf0)
            df0[rf_alt1_idx, :value] = string(rf_alt1)
            df0[rf_alt2_idx, :value] = string(rf_alt2)
            df0[rf_alt3_idx, :value] = string(rf_alt3)
            df0[rf_true_noF_idx, :value] = string(rf0_noF)
        end
        CSV.write(results0_file, df0)
    end
    
    # Calculate RF distances for H1 displayed trees
    if isfile(net1_file) && isfile(results1_file)
        df1 = CSV.read(results1_file, DataFrame; types=Dict(:value => String))
        net1 = readnewick(net1_file)
        # deepcopy to avoid mutating net1 before RF calculations
        hybrid_info = get_hybrid_info(deepcopy(net1), "A")
        displayed_trees = displayedtrees(net1, 0.0)
        
        if length(displayed_trees) != 2
            error("Expected 2 displayed trees for $net1_file," *
                " found $(length(displayed_trees))")
        end
        
        rf1 = mudistance_semidirected(
            displayed_trees[1], true_species_tree, preorder=true)
        rf2 = mudistance_semidirected(
            displayed_trees[2], true_species_tree, preorder=true)

        # Also check if either displayed tree matches species tree with no F
        rf1_noF = NaN
        rf2_noF = NaN
        rf1_noG = NaN
        rf2_noG = NaN 
        rf1_noH = NaN 
        rf2_noH = NaN 
        try
            displayed_tree1_noF = deepcopy(displayed_trees[1])
            displayed_tree1_noG = deepcopy(displayed_trees[1]) 
            displayed_tree1_noH = deepcopy(displayed_trees[1]) 
            deleteleaf!(displayed_tree1_noF, "F") 
            rf1_noF = mudistance_semidirected(
                displayed_tree1_noF, species_tree_noF, preorder=true)
            deleteleaf!(displayed_tree1_noG, "G")
            rf1_noG = mudistance_semidirected(
                displayed_tree1_noG, species_tree_noG, preorder=true)
            deleteleaf!(displayed_tree1_noH, "H")
            rf1_noH = mudistance_semidirected(
                displayed_tree1_noH, species_tree_noH, preorder=true)
        catch e
            println(
                "Warning: RF_net1_1_noF/G/H failed rep $rep_number_string: $e")
        end
        try
            displayed_tree2_noF = deepcopy(displayed_trees[2])
            deleteleaf!(displayed_tree2_noF, "F")
            rf2_noF = mudistance_semidirected(
                displayed_tree2_noF, species_tree_noF, preorder=true)
        catch e
            println("Warning: RF_net1_2_noF failed rep $rep_number_string: $e")
        end
        try
            displayed_tree2_noG = deepcopy(displayed_trees[2])
            deleteleaf!(displayed_tree2_noG, "G")
            rf2_noG = mudistance_semidirected(
                displayed_tree2_noG, species_tree_noG, preorder=true)
        catch e
            println("Warning: RF_net1_2_noG failed rep $rep_number_string: $e")
        end
        try
            displayed_tree2_noH = deepcopy(displayed_trees[2])
            deleteleaf!(displayed_tree2_noH, "H")
            rf2_noH = mudistance_semidirected(
                displayed_tree2_noH, species_tree_noH, preorder=true)
        catch e
            println("Warning: RF_net1_2_noH failed rep $rep_number_string: $e")
        end

        # two idx lookups done first so re-runs don't clobber existing RF values
        rf1_idx = findfirst(row -> row.type == "RF_net1_1_true", eachrow(df1))
        rf2_idx = findfirst(row -> row.type == "RF_net1_2_true", eachrow(df1))

        rf_alt1_tree1, rf_alt2_tree1, rf_alt3_tree1 =
            check_alternative_trees(displayed_trees[1])
        rf_alt1_tree2, rf_alt2_tree2, rf_alt3_tree2 =
            check_alternative_trees(displayed_trees[2])
        
        findr(t) = findfirst(row -> row.type == t, eachrow(df1))
        rf_alt1_tree1_idx = findr("RF_net1_1_alter1")
        rf_alt2_tree1_idx = findr("RF_net1_1_alter2")
        rf_alt3_tree1_idx = findr("RF_net1_1_alter3")
        rf_alt1_tree2_idx = findr("RF_net1_2_alter1")
        rf_alt2_tree2_idx = findr("RF_net1_2_alter2")
        rf_alt3_tree2_idx = findr("RF_net1_2_alter3")
        rf_true_noF_idx1  = findr("RF_net1_1_true_noF")
        rf_true_noF_idx2  = findr("RF_net1_2_true_noF")
        rf_true_noG_idx1  = findr("RF_net1_1_true_noG")
        rf_true_noG_idx2  = findr("RF_net1_2_true_noG")
        rf_true_noH_idx1  = findr("RF_net1_1_true_noH")
        rf_true_noH_idx2  = findr("RF_net1_2_true_noH")

        # Check if any rows are missing and add/update accordingly
        if rf1_idx === nothing || rf2_idx === nothing || 
            rf_alt1_tree1_idx === nothing || rf_alt2_tree1_idx === nothing || 
            rf_alt3_tree1_idx === nothing || rf_alt1_tree2_idx === nothing || 
            rf_alt2_tree2_idx === nothing || rf_alt3_tree2_idx === nothing || 
            rf_true_noF_idx1 === nothing || rf_true_noF_idx2 === nothing || 
            rf_true_noG_idx1 === nothing || rf_true_noG_idx2 === nothing || 
            rf_true_noH_idx1 === nothing || rf_true_noH_idx2 === nothing
            # Add or update each row individually
            if rf1_idx === nothing
                push!(df1, (type="RF_net1_1_true", value=string(rf1)))
            else
                df1[rf1_idx, :value] = string(rf1)
            end
            if rf2_idx === nothing
                push!(df1, (type="RF_net1_2_true", value=string(rf2)))
            else
                df1[rf2_idx, :value] = string(rf2)
            end
            function upsert!(df, idx, t, v)
                if idx === nothing
                    push!(df, (type=t, value=string(v)))
                else
                    df[idx, :value] = string(v)
                end
            end
            upsert!(df1, rf_alt1_tree1_idx, "RF_net1_1_alter1", rf_alt1_tree1)
            upsert!(df1, rf_alt2_tree1_idx, "RF_net1_1_alter2", rf_alt2_tree1)
            upsert!(df1, rf_alt3_tree1_idx, "RF_net1_1_alter3", rf_alt3_tree1)
            upsert!(df1, rf_alt1_tree2_idx, "RF_net1_2_alter1", rf_alt1_tree2)
            upsert!(df1, rf_alt2_tree2_idx, "RF_net1_2_alter2", rf_alt2_tree2)
            upsert!(df1, rf_alt3_tree2_idx, "RF_net1_2_alter3", rf_alt3_tree2)
            if rf_true_noF_idx1 === nothing
                push!(df1, (type="RF_net1_1_true_noF", value=string(rf1_noF)))
            else
                df1[rf_true_noF_idx1, :value] = string(rf1_noF)
            end
            if rf_true_noF_idx2 === nothing
                push!(df1, (type="RF_net1_2_true_noF", value=string(rf2_noF)))
            else
                df1[rf_true_noF_idx2, :value] = string(rf2_noF)
            end
            if rf_true_noG_idx1 === nothing
                push!(df1, (type="RF_net1_1_true_noG", value=string(rf1_noG)))
            else
                df1[rf_true_noG_idx1, :value] = string(rf1_noG)
            end
            if rf_true_noG_idx2 === nothing
                push!(df1, (type="RF_net1_2_true_noG", value=string(rf2_noG)))
            else
                df1[rf_true_noG_idx2, :value] = string(rf2_noG)
            end
            if rf_true_noH_idx1 === nothing
                push!(df1, (type="RF_net1_1_true_noH", value=string(rf1_noH)))
            else
                df1[rf_true_noH_idx1, :value] = string(rf1_noH)
            end
            if rf_true_noH_idx2 === nothing
                push!(df1, (type="RF_net1_2_true_noH", value=string(rf2_noH)))
            else
                df1[rf_true_noH_idx2, :value] = string(rf2_noH)
            end
        else
            # All rows exist, just update them
            df1[rf1_idx, :value] = string(rf1)
            df1[rf2_idx, :value] = string(rf2)
            df1[rf_alt1_tree1_idx, :value] = string(rf_alt1_tree1)
            df1[rf_alt2_tree1_idx, :value] = string(rf_alt2_tree1)
            df1[rf_alt3_tree1_idx, :value] = string(rf_alt3_tree1)
            df1[rf_alt1_tree2_idx, :value] = string(rf_alt1_tree2)
            df1[rf_alt2_tree2_idx, :value] = string(rf_alt2_tree2)
            df1[rf_alt3_tree2_idx, :value] = string(rf_alt3_tree2)
            df1[rf_true_noF_idx1, :value] = string(rf1_noF)
            df1[rf_true_noF_idx2, :value] = string(rf2_noF)
            df1[rf_true_noG_idx1, :value] = string(rf1_noG)
            df1[rf_true_noG_idx2, :value] = string(rf2_noG)
            df1[rf_true_noH_idx1, :value] = string(rf1_noH)
            df1[rf_true_noH_idx2, :value] = string(rf2_noH)
        end
        # Upsert hybrid_taxon / major_donor / minor_donor rows
        for (rtype, rval) in [("hybrid_taxon", hybrid_info.hybrid_taxon),
                               ("major_donor",  hybrid_info.major_donor),
                               ("minor_donor",  hybrid_info.minor_donor)]
            row_idx = findfirst(r -> r.type == rtype, eachrow(df1))
            if row_idx === nothing
                push!(df1, (type=rtype, value=rval))
            else
                df1[row_idx, :value] = rval
            end
        end
        CSV.write(results1_file, df1)
    end
end

#-----------------------------------------------#       
#     Analyze and update goodness-of-fit summary
#-----------------------------------------------#

println("Analyzing goodness-of-fit results from all replicates...")

# Initialize arrays to store results
rep_ids = String[]
p_H0_values = Float64[]
p_H1_values = Float64[]
z_uncorrected_H0_values = Float64[]
z_uncorrected_H1_values = Float64[]
sigma_H0_values = Float64[]
sigma_H1_values = Float64[]
score_H0_values = Float64[]
score_H1_values = Float64[]
gamma_1_values = Float64[]
gamma_2_values = Float64[]
rf_net0_true_values = Float64[]
rf_net1_1_true_values = Float64[]
rf_net1_2_true_values = Float64[]
RF_net0_alter1_values = Float64[]
RF_net0_alter2_values = Float64[]
RF_net0_alter3_values = Float64[]
RF_net1_1_alter1_values = Float64[]
RF_net1_1_alter2_values = Float64[]
RF_net1_1_alter3_values = Float64[]
RF_net1_2_alter1_values = Float64[]
RF_net1_2_alter2_values = Float64[]
RF_net1_2_alter3_values = Float64[] 
RF_net0_true_noF_values = Float64[] 
RF_net1_1_true_noF_values = Float64[]
RF_net1_2_true_noF_values = Float64[]
RF_net1_1_true_noG_values = Float64[]
RF_net1_2_true_noG_values = Float64[]
RF_net1_1_true_noH_values = Float64[]
RF_net1_2_true_noH_values = Float64[]
hybrid_taxon_values = String[]
major_donor_values = String[]
minor_donor_values = String[]

# Process each replicate
for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
    
    # Construct paths to the goodness-of-fit result files
    gof_H0_file = joinpath(rep_folder_path, "snaqfolder", "H0_output", 
                            "snaq_gof_results_H0.csv")
    gof_H1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", 
                            "snaq_gof_results_H1.csv")
    
    # Initialize values
    p_H0, p_H1 = NaN, NaN
    z_uncorrected_H0, z_uncorrected_H1 = NaN, NaN
    sigma_H0, sigma_H1 = NaN, NaN
    score_H0, score_H1 = NaN, NaN
    gamma_1, gamma_2 = NaN, NaN
    rf_net0_true = NaN
    rf_net0_true_noF = NaN 
    rf_net1_1_true, rf_net1_2_true = NaN, NaN
    rf_net1_1_true_noF, rf_net1_2_true_noF = NaN, NaN 
    rf_net1_1_true_noG, rf_net1_2_true_noG = NaN, NaN 
    rf_net1_1_true_noH, rf_net1_2_true_noH = NaN, NaN  
    hybrid_taxon, major_donor, minor_donor = "NA", "NA", "NA"
    rf_net0_alter1_values, rf_net0_alter2_values,
        rf_net0_alter3_values = NaN, NaN, NaN
    rf_net1_1_alter1_values, rf_net1_1_alter2_values,
        rf_net1_1_alter3_values = NaN, NaN, NaN
    rf_net1_2_alter1_values, rf_net1_2_alter2_values,
        rf_net1_2_alter3_values = NaN, NaN, NaN
    
    # Parse H0 results
    if isfile(gof_H0_file)
        try
            df_H0 = CSV.read(gof_H0_file, DataFrame;
                types=Dict(:value => String))
            
            for row in eachrow(df_H0)
                if row.type == "p"
                    p_H0 = parse(Float64, string(row.value))
                elseif row.type == "z_uncorrected"
                    z_uncorrected_H0 = parse(Float64, string(row.value))
                elseif row.type == "sigma"
                    sigma_H0 = parse(Float64, string(row.value))
                elseif row.type == "score"
                    score_H0 = parse(Float64, string(row.value))
                elseif row.type == "RF_net0_true"
                    rf_net0_true = parse(Float64, string(row.value))
                elseif row.type == "RF_net0_true_noF"
                    rf_net0_true_noF = parse(Float64, string(row.value))
                elseif row.type == "RF_net0_alter1"
                    rf_net0_alter1_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net0_alter2"
                    rf_net0_alter2_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net0_alter3"
                    rf_net0_alter3_values = parse(Float64, string(row.value))
                end
            end
        catch e
            println("Warning: Failed to parse $gof_H0_file - $e")
        end
    end
    
    # Parse H1 results
    if isfile(gof_H1_file)
        try
            df_H1 = CSV.read(gof_H1_file, DataFrame;
                types=Dict(:value => String))
            
            for row in eachrow(df_H1)
                if row.type == "p"
                    p_H1 = parse(Float64, string(row.value))
                elseif row.type == "z_uncorrected"
                    z_uncorrected_H1 = parse(Float64, string(row.value))
                elseif row.type == "sigma"
                    sigma_H1 = parse(Float64, string(row.value))
                elseif row.type == "score"
                    score_H1 = parse(Float64, string(row.value))
                elseif row.type == "gamma_1"
                    gamma_1 = parse(Float64, string(row.value))
                elseif row.type == "gamma_2"
                    gamma_2 = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_1_true"
                    rf_net1_1_true = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_true"
                    rf_net1_2_true = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_1_true_noF"
                    rf_net1_1_true_noF = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_true_noF"
                    rf_net1_2_true_noF = parse(Float64, string(row.value)) 
                elseif row.type == "RF_net1_1_true_noG"
                    rf_net1_1_true_noG = parse(Float64, string(row.value)) 
                elseif row.type == "RF_net1_2_true_noG"
                    rf_net1_2_true_noG = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_1_true_noH"
                    rf_net1_1_true_noH = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_true_noH"
                    rf_net1_2_true_noH = parse(Float64, string(row.value))  
                elseif row.type == "RF_net1_1_alter1"
                    rf_net1_1_alter1_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_1_alter2"
                    rf_net1_1_alter2_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_1_alter3"
                    rf_net1_1_alter3_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_alter1"
                    rf_net1_2_alter1_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_alter2"
                    rf_net1_2_alter2_values = parse(Float64, string(row.value))
                elseif row.type == "RF_net1_2_alter3"
                    rf_net1_2_alter3_values = parse(Float64, string(row.value))
                elseif row.type == "hybrid_taxon"
                    hybrid_taxon = string(row.value)
                elseif row.type == "major_donor"
                    major_donor = string(row.value)
                elseif row.type == "minor_donor"
                    minor_donor = string(row.value)
                end
            end
            
            #= Ensure gamma_1 is always bigger than gamma_2 
            This makes sure gamma1 is the major and gamma2 is the minor =# 
            if !isnan(gamma_1) && !isnan(gamma_2) && gamma_1 < gamma_2
                gamma_1, gamma_2 = gamma_2, gamma_1
            end
            
        catch e

            println("Warning: Failed to parse $gof_H1_file - $e")
        end
    end
    
    # Store results
    push!(rep_ids, rep_number_string)
    push!(p_H0_values, p_H0)
    push!(p_H1_values, p_H1)
    push!(z_uncorrected_H0_values, z_uncorrected_H0)
    push!(z_uncorrected_H1_values, z_uncorrected_H1)
    push!(sigma_H0_values, sigma_H0)
    push!(sigma_H1_values, sigma_H1)
    push!(score_H0_values, score_H0)
    push!(score_H1_values, score_H1)
    push!(gamma_1_values, gamma_1)
    push!(gamma_2_values, gamma_2)
    push!(rf_net0_true_values, rf_net0_true)
    push!(rf_net1_1_true_values, rf_net1_1_true)
    push!(rf_net1_2_true_values, rf_net1_2_true)
    push!(RF_net0_alter1_values, rf_net0_alter1_values)
    push!(RF_net0_alter2_values, rf_net0_alter2_values)
    push!(RF_net0_alter3_values, rf_net0_alter3_values)
    push!(RF_net1_1_alter1_values, rf_net1_1_alter1_values)
    push!(RF_net1_1_alter2_values, rf_net1_1_alter2_values)
    push!(RF_net1_1_alter3_values, rf_net1_1_alter3_values)
    push!(RF_net1_2_alter1_values, rf_net1_2_alter1_values)
    push!(RF_net1_2_alter2_values, rf_net1_2_alter2_values)
    push!(RF_net1_2_alter3_values, rf_net1_2_alter3_values)
    push!(RF_net0_true_noF_values, rf_net0_true_noF)
    push!(RF_net1_1_true_noF_values, rf_net1_1_true_noF)
    push!(RF_net1_2_true_noF_values, rf_net1_2_true_noF) 
    push!(RF_net1_1_true_noG_values, rf_net1_1_true_noG)
    push!(RF_net1_2_true_noG_values, rf_net1_2_true_noG) 
    push!(RF_net1_1_true_noH_values, rf_net1_1_true_noH)
    push!(RF_net1_2_true_noH_values, rf_net1_2_true_noH)
    push!(hybrid_taxon_values, hybrid_taxon)
    push!(major_donor_values, major_donor)
    push!(minor_donor_values, minor_donor)
    
    println("Processed rep$rep_number_string: " *
        "p_H0=$p_H0, p_H1=$p_H1, RF_H0=$rf_net0_true")
end

# Create comprehensive summary DataFrame
summary_df = DataFrame(
    repID = rep_ids,
    p_H0 = p_H0_values,
    p_H1 = p_H1_values,
    z_uncorrected_H0 = z_uncorrected_H0_values,
    z_uncorrected_H1 = z_uncorrected_H1_values,
    sigma_H0 = sigma_H0_values,
    sigma_H1 = sigma_H1_values,
    score_H0 = score_H0_values,
    score_H1 = score_H1_values,
    gamma_1 = gamma_1_values,
    gamma_2 = gamma_2_values,
    RF_net0_true = rf_net0_true_values,
    RF_net0_true_noF = RF_net0_true_noF_values,
    RF_net1_1_true = rf_net1_1_true_values,
    RF_net1_2_true = rf_net1_2_true_values,
    RF_net1_1_true_noF = RF_net1_1_true_noF_values,
    RF_net1_2_true_noF = RF_net1_2_true_noF_values,
    RF_net1_1_true_noG = RF_net1_1_true_noG_values,
    RF_net1_2_true_noG = RF_net1_2_true_noG_values,
    RF_net1_1_true_noH = RF_net1_1_true_noH_values,
    RF_net1_2_true_noH = RF_net1_2_true_noH_values, 
    RF_net0_alter1 = RF_net0_alter1_values,
    RF_net0_alter2 = RF_net0_alter2_values,
    RF_net0_alter3 = RF_net0_alter3_values,
    RF_net1_1_alter1 = RF_net1_1_alter1_values,
    RF_net1_1_alter2 = RF_net1_1_alter2_values,
    RF_net1_1_alter3 = RF_net1_1_alter3_values,
    RF_net1_2_alter1 = RF_net1_2_alter1_values,
    RF_net1_2_alter2 = RF_net1_2_alter2_values,
    RF_net1_2_alter3 = RF_net1_2_alter3_values,
    hybrid_taxon = hybrid_taxon_values,
    major_donor = major_donor_values,
    minor_donor = minor_donor_values
)

# Save updated summary results
CSV.write(gof_summary_file, summary_df)

if create_summary
    println("Created new goodness-of-fit summary: $gof_summary_file")
else
    println("Updated existing goodness-of-fit summary: $gof_summary_file")
end
println("Summary statistics:")
println("  Total replicates processed: $(length(rep_ids))")
println("  H0 files found: $(count(x -> !isnan(x), p_H0_values))")
println("  H1 files found: $(count(x -> !isnan(x), p_H1_values))")
println("  RF calculations available:")
println("    RF_net0_true: $(count(x -> !isnan(x), rf_net0_true_values))")

#-----------------------------------------------#
#   Build consensus trees/networks across replicates
#   (treat each replicate as one bootstrap sample)
#-----------------------------------------------#

println("\nBuilding consensus trees/networks across replicates " *
    "(treating each rep as a bootstrap)...")

h0_consensus_trees = HybridNetwork[]
h1_consensus_nets  = HybridNetwork[]

for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
    net0_file = joinpath(rep_folder_path, "snaqfolder", "H0_output", "H0.out")
    net1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", "H1.out")

    if isfile(net0_file)
        try
            push!(h0_consensus_trees, readnewick(net0_file))
        catch e
            println("Warning: could not read H0 tree for " *
                "rep $rep_number_string: $e")
        end
    end

    if isfile(net1_file)
        try
            push!(h1_consensus_nets, readnewick(net1_file))
        catch e
            println("Warning: could not read H1 network for " *
                "rep $rep_number_string: $e")
        end
    end
end

println("Collected $(length(h0_consensus_trees)) H=0 trees and " *
    "$(length(h1_consensus_nets)) H=1 networks for consensus")

consensus_dir  = joinpath(outfolder, "SNaQ-$(paramname_root)-consensus_nets")
mkpath(consensus_dir)
println("Consensus output directory: $consensus_dir")
con_h0_prefix  = joinpath(consensus_dir, "SNaQ-$(paramname_root)-consensus_H0")
con_h0_nwk     = "$(con_h0_prefix).nwk"
con_h0_sup_csv = "$(con_h0_prefix)_edge_support.csv"
con_h1_prefix  = joinpath(consensus_dir, "SNaQ-$(paramname_root)-consensus_H1")

# H=0: greedy consensus tree (unrooted, no outgroup specified)
if length(h0_consensus_trees) >= 2
    println("Computing H=0 consensus tree from " *
        "$(length(h0_consensus_trees)) trees...")
    con_h0 = consensustree(h0_consensus_trees)
    rootatnode!(con_h0, "A")
    writenewick(con_h0, con_h0_nwk; support=true)
    esup_h0 = DataFrame(
        edge_number = [e.number for e in con_h0.edge if !isexternal(e)],
        support = [round(e.y, digits=6) for e in con_h0.edge if !isexternal(e)]
    )
    CSV.write(con_h0_sup_csv, esup_h0)

    # Plot H=0 consensus tree while con_h0 is live
    # (e.y is not preserved on file round-trip)
    con_h0_plot = "$(con_h0_prefix)_plot.pdf"
    esup_h0_plot = DataFrame(
        number  = [e.number for e in con_h0.edge if !isexternal(e)],
        support = [round(100 * e.y, digits=1)
            for e in con_h0.edge if !isexternal(e)]
    )
    try
        R"pdf($con_h0_plot, width=8, height=6)"
        R"par(mar=c(1,1,3,1))"
        plot(con_h0;
             edgelabel     = esup_h0_plot,
             edgewidth     = 3,
             tipcex        = 1.5,
             edgecex       = 0.9,
             edgelabeladj  = [0.5, -0.3])
        title_h0 = "H=0 consensus tree — SNaQ: $(paramname_root)"
        R"title(main=$title_h0, cex.main=0.85)"
        R"dev.off()"
        println("H=0 plot saved: $con_h0_plot")
    catch err
        @warn "H=0 plot failed: $err"
        try; R"dev.off()"; catch; end
    end
println("H=0 consensus tree computed.")
else
    println("Warning: Fewer than 2 H=0 trees available – " *
        "skipping H=0 consensus (found $(length(h0_consensus_trees))).")
    con_h0_nwk = nothing
    con_h0_sup_csv = nothing
end

# H=1: consensus level-1 network (unrooted, no outgroup specified)
if length(h1_consensus_nets) >= 2
    println("Computing H=1 consensus level-1 network from " *
        "$(length(h1_consensus_nets))...")
    res_h1 = consensus_level1network(h1_consensus_nets, 
            outgroup = "A", 
            suppressinfo=true)
    consensus_level1network_save(res_h1, con_h1_prefix)

    # Plot H=1 consensus network while res_h1 is live
    # (blob/hybrid tables have :edge column)
    con_h1_plot = "$(con_h1_prefix)_plot.pdf"
    blb_df = DataFrame(res_h1[:blob_table],   copycols=false)
    hyb_df = DataFrame(res_h1[:hybrid_table], copycols=false)
    try
        R"pdf($con_h1_plot, width=8, height=6)"
        R"par(mar=c(1,1,3,1))"
        plot(res_h1[:net];
             edgewidth      = 3,
             tipcex         = 1.5,
             nodelabeladj   = -0.1,
             edgelabeladj   = [0.5, -0.3],
             nodelabelcolor = "orangered",
             edgelabelcolor = "deepskyblue",
             nodelabel = select(blb_df, [:node, :support_partition]),
             edgelabel = select(hyb_df, [:edge, :support_hybrid]))
        title_h1 = "H=1 consensus network — SNaQ: $(paramname_root)"
        R"title(main=$title_h1, cex.main=0.85)"
        R"dev.off()"
        println("H=1 plot saved: $con_h1_plot")
    catch err
        @warn "H=1 plot failed: $err"
        try; R"dev.off()"; catch; end
    end
    println("H=1 consensus network computed.")
else
    println("Warning: Fewer than 2 H=1 networks available – " *
        "skipping H=1 consensus (found $(length(h1_consensus_nets))).")
end

println()
println("========================================")
println("SNaQ Postprocessing: Output Files Summary")
println("========================================")
println("  GOF summary CSV       : $gof_summary_file")
if con_h0_nwk !== nothing 
    println("Saved to consensus directory: $consensus_dir") 
end
println("========================================")

println()
println("SNaQ Postprocessing completed successfully!")
