# SNaQ Postprocessing Script
# This script runs after snaq.jl to perform additional analyses and fix failed calculations
# Usage: julia snaq_postprocessing.jl --dup_rate X --loss_rate Y --ratevar Z --n_reps N [--rep_start S --rep_end E]

using ArgParse
using Distributed
@everywhere using Printf
@everywhere using PhyloNetworks
@everywhere using CSV
@everywhere using DataFrames
@everywhere include("utilities.jl")

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
alter1 = readnewick("(A,((G,H),(((B,C),F),(D,E))));") # F is clustered with (B,C) 
alter2 = readnewick("(A,((G,H),((B,C),(F,(D,E)))));") # F is clustered with (D,E)  
alter3 = readnewick("(A,(((G,H),((D,E),(B,C))),F));") # F is outside of (B,C),(D,E),(G,H)
species_tree_noF = readnewick("(A,((G,H),((B,C),(D,E))));") # species tree without F 

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
    results0_file = joinpath(rep_folder_path, "snaqfolder", "H0_output", "snaq_gof_results_H0.csv")
    net1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", "H1.networks")
    results1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", "snaq_gof_results_H1.csv")
    
    # Calculate RF distance for H0
    if isfile(net0_file) && isfile(results0_file)
        df0 = CSV.read(results0_file, DataFrame)
        net0 = readnewick(net0_file) 
        rf0 = mudistance_semidirected(net0, true_species_tree, preorder=true) 

        # Check if estimated tree matches species tree with no F
        rf0_noF = NaN
        try
            net0_noF = deepcopy(net0)
            deleteleaf!(net0_noF, "F") 
            rf0_noF = mudistance_semidirected(net0_noF, species_tree_noF, preorder=true)
        catch e
            println("Warning: Could not calculate RF_net0_noF for rep $rep_number_string: $e")
        end 

        (rf_alt1, rf_alt2, rf_alt3) = check_alternative_trees(net0) 
        println("Rep $rep_number_string: RF to true species tree = $rf0;" *
        "RF to alter1 = $rf_alt1; RF to alter2 = $rf_alt2; RF to alter3 = $rf_alt3") 
        
        rf_row_idx = findfirst(row -> row.type == "RF_net0_true", eachrow(df0))
        rf_alt1_idx = findfirst(row -> row.type == "RF_net0_alter1", eachrow(df0))
        rf_alt2_idx = findfirst(row -> row.type == "RF_net0_alter2", eachrow(df0))
        rf_alt3_idx = findfirst(row -> row.type == "RF_net0_alter3", eachrow(df0))
        rf_true_noF_idx = findfirst(row -> row.type == "RF_net0_true_noF", eachrow(df0))
        
        if rf_row_idx === nothing || rf_alt1_idx === nothing || rf_alt2_idx === nothing || rf_alt3_idx === nothing || rf_true_noF_idx === nothing
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
        df1 = CSV.read(results1_file, DataFrame)
        net1 = readnewick(net1_file)
        displayed_trees = displayedtrees(net1, 0.0)
        
        if length(displayed_trees) != 2
            error("Expected 2 displayed trees for $net1_file, but found $(length(displayed_trees))")
        end
        
        rf1 = mudistance_semidirected(displayed_trees[1], true_species_tree, preorder=true)
        rf2 = mudistance_semidirected(displayed_trees[2], true_species_tree, preorder=true)

        # Also check if either displayed tree matches species tree with no F
        rf1_noF = NaN
        rf2_noF = NaN
        try
            displayed_tree1_noF = deepcopy(displayed_trees[1])
            deleteleaf!(displayed_tree1_noF, "F") 
            rf1_noF = mudistance_semidirected(displayed_tree1_noF, species_tree_noF, preorder=true)
        catch e
            println("Warning: Could not calculate RF_net1_1_noF for rep $rep_number_string: $e")
        end
        try
            displayed_tree2_noF = deepcopy(displayed_trees[2])
            deleteleaf!(displayed_tree2_noF, "F")
            rf2_noF = mudistance_semidirected(displayed_tree2_noF, species_tree_noF, preorder=true)
        catch e
            println("Warning: Could not calculate RF_net1_2_noF for rep $rep_number_string: $e")
        end 
        
        # Keep those two lines below because this keeps us to re-run snaq-postprocess.jl 
        # without messing up existing RF values 
        rf1_idx = findfirst(row -> row.type == "RF_net1_1_true", eachrow(df1))
        rf2_idx = findfirst(row -> row.type == "RF_net1_2_true", eachrow(df1))

        rf_alt1_tree1, rf_alt2_tree1, rf_alt3_tree1 = check_alternative_trees(displayed_trees[1])
        rf_alt1_tree2, rf_alt2_tree2, rf_alt3_tree2 = check_alternative_trees(displayed_trees[2])
        
        rf_alt1_tree1_idx = findfirst(row -> row.type == "RF_net1_1_alter1", eachrow(df1))
        rf_alt2_tree1_idx = findfirst(row -> row.type == "RF_net1_1_alter2", eachrow(df1))
        rf_alt3_tree1_idx = findfirst(row -> row.type == "RF_net1_1_alter3", eachrow(df1))
        rf_alt1_tree2_idx = findfirst(row -> row.type == "RF_net1_2_alter1", eachrow(df1))
        rf_alt2_tree2_idx = findfirst(row -> row.type == "RF_net1_2_alter2", eachrow(df1))
        rf_alt3_tree2_idx = findfirst(row -> row.type == "RF_net1_2_alter3", eachrow(df1))
        rf_true_noF_idx1 = findfirst(row -> row.type == "RF_net1_1_true_noF", eachrow(df1))
        rf_true_noF_idx2 = findfirst(row -> row.type == "RF_net1_2_true_noF", eachrow(df1))

        # Check if any rows are missing and add/update accordingly
        if rf1_idx === nothing || rf2_idx === nothing || rf_alt1_tree1_idx === nothing || rf_alt2_tree1_idx === nothing || rf_alt3_tree1_idx === nothing || rf_alt1_tree2_idx === nothing || rf_alt2_tree2_idx === nothing || rf_alt3_tree2_idx === nothing || rf_true_noF_idx1 === nothing || rf_true_noF_idx2 === nothing
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
            if rf_alt1_tree1_idx === nothing
                push!(df1, (type="RF_net1_1_alter1", value=string(rf_alt1_tree1)))
            else
                df1[rf_alt1_tree1_idx, :value] = string(rf_alt1_tree1)
            end
            if rf_alt2_tree1_idx === nothing
                push!(df1, (type="RF_net1_1_alter2", value=string(rf_alt2_tree1)))
            else
                df1[rf_alt2_tree1_idx, :value] = string(rf_alt2_tree1)
            end
            if rf_alt3_tree1_idx === nothing
                push!(df1, (type="RF_net1_1_alter3", value=string(rf_alt3_tree1)))
            else
                df1[rf_alt3_tree1_idx, :value] = string(rf_alt3_tree1)
            end
            if rf_alt1_tree2_idx === nothing
                push!(df1, (type="RF_net1_2_alter1", value=string(rf_alt1_tree2)))
            else
                df1[rf_alt1_tree2_idx, :value] = string(rf_alt1_tree2)
            end
            if rf_alt2_tree2_idx === nothing
                push!(df1, (type="RF_net1_2_alter2", value=string(rf_alt2_tree2)))
            else
                df1[rf_alt2_tree2_idx, :value] = string(rf_alt2_tree2)
            end
            if rf_alt3_tree2_idx === nothing
                push!(df1, (type="RF_net1_2_alter3", value=string(rf_alt3_tree2)))
            else
                df1[rf_alt3_tree2_idx, :value] = string(rf_alt3_tree2)
            end
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

# Process each replicate
for simulation_rep in rep_start:rep_end
    rep_number_string = pad_number(simulation_rep, n_reps)
    rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
    
    # Construct paths to the goodness-of-fit result files
    gof_H0_file = joinpath(rep_folder_path, "snaqfolder", "H0_output", "snaq_gof_results_H0.csv")
    gof_H1_file = joinpath(rep_folder_path, "snaqfolder", "H1_output", "snaq_gof_results_H1.csv")
    
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
    rf_net0_alter1_values, rf_net0_alter2_values, rf_net0_alter3_values = NaN, NaN, NaN
    rf_net1_1_alter1_values, rf_net1_1_alter2_values, rf_net1_1_alter3_values = NaN, NaN, NaN
    rf_net1_2_alter1_values, rf_net1_2_alter2_values, rf_net1_2_alter3_values = NaN, NaN, NaN
    
    # Parse H0 results
    if isfile(gof_H0_file)
        try
            df_H0 = CSV.read(gof_H0_file, DataFrame)
            
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
            df_H1 = CSV.read(gof_H1_file, DataFrame)
            
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
    
    println("Processed rep$rep_number_string: p_H0=$p_H0, p_H1=$p_H1, RF_H0=$rf_net0_true")
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
    RF_net0_alter1 = RF_net0_alter1_values,
    RF_net0_alter2 = RF_net0_alter2_values,
    RF_net0_alter3 = RF_net0_alter3_values,
    RF_net1_1_alter1 = RF_net1_1_alter1_values,
    RF_net1_1_alter2 = RF_net1_1_alter2_values,
    RF_net1_1_alter3 = RF_net1_1_alter3_values,
    RF_net1_2_alter1 = RF_net1_2_alter1_values,
    RF_net1_2_alter2 = RF_net1_2_alter2_values,
    RF_net1_2_alter3 = RF_net1_2_alter3_values 
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

println()
println("SNaQ Postprocessing completed successfully!")
