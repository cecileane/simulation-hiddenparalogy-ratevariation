#=
This script postprocesses simulation results by calculating RF distances
between true and estimated trees.

Usage:
julia scripts/simulation_postprocess.jl --dup_rate 0.0 --loss_rate 0.0 --ratevar N --n_reps 100

Output: Creates CSV files with RF distance comparisons in the output folder.
=#

using Distributed
using ArgParse
using Statistics

using DataFrames
using CSV
using PhyloNetworks
using Printf

# Ensure utilities are included from the script directory so paths are stable
const _UTILS_PATH = joinpath(@__DIR__, "utilities.jl")
include(_UTILS_PATH)                 # include on master
@everywhere include($(_UTILS_PATH)) # include on workers (absolute path)

function parse_commandline()
  s = ArgParseSettings() 
  @add_arg_table s begin 

    "--dup_rate"
      help = """Gene duplication rate: specify the gene duplication rate.
                                      If 0, then no duplication"""
      arg_type = Float64
      required = true

    "--loss_rate" 
      help = """Gene loss rate: specify the gene loss rate, 
                                if 0, then no gene loss"""
      arg_type = Float64
      required = true

    "--ratevar" 
      help = """'N': No rate variation';
        'G: Gene specific rate variation'; 
        'L': Lineage specific rate variation;
        'GL' or 'G*L': genexlineage rate variation"""
      arg_type = String
      required = true

    "--n_reps"
      help = "Number of replicates" 
      arg_type = Int
      required = true

    "--n_inds" 
      help = "Number of individuals/accessions per species (Default = 1)"  
      arg_type = Int
      default = 1

    "--SF"
      help = "Scaling factor to scale the branch lengths of the species tree
        (Default = 1.0, no scaling)"
      arg_type = Float64
      default = 1.0
    
    "--gene_len" 
      help = "The length of simulated gene sequences (Default = 1000 bp)" 
      arg_type = Int 
      default = 1000

  end 
  return parse_args(s) 
end

#--------------- Parse arguments ---------------# 
parsed_args = parse_commandline()
dup_rate = parsed_args["dup_rate"] 
loss_rate = parsed_args["loss_rate"]
ratevar = parsed_args["ratevar"] 
n_reps = parsed_args["n_reps"] 
n_inds = parsed_args["n_inds"] 
scaling_factor_branch_length = parsed_args["SF"]
scaling_factor_branch_length = Float64(scaling_factor_branch_length)  
gene_len = parsed_args["gene_len"] 

#--------------- set up folders and paths ---------------# 
# Use the repository (script) directory as the rootfolder so paths are stable
rootfolder = normpath(joinpath(@__DIR__, ".."))

# set up paramname_root based on given parameter values -> see utilities.jl for details 
paramname_root = set_up_paramname_root(dup_rate, 
                                      loss_rate, 
                                      ratevar, 
                                      n_inds, 
                                      scaling_factor_branch_length,
                                      gene_len)
outfolder = joinpath(rootfolder, "output", paramname_root) 

if !isdir(outfolder)
    error("Output folder does not exist: $outfolder. Please run the simulation first.")
end

folder_path_list = [] 
for n in 1:n_reps 
  rep_number_string = pad_number(n, n_reps)
  rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
  push!(folder_path_list, rep_folder_path)
end 

# Set global variables for @everywhere functions
@everywhere global dup_rate = $dup_rate 
@everywhere global n_inds = $n_inds 
@everywhere global paramname_root = $paramname_root
@everywhere global outfolder = $outfolder
@everywhere global n_reps = $n_reps
@everywhere global rootfolder = $rootfolder
@everywhere global folder_path_list = $folder_path_list


#-----------------------------------------------#
# Create SimPhy vs species tree RF summary
# Calculate RF distance between
# SimPhy (true) gene trees vs true species tree 
# Only for n_inds = 1
#-----------------------------------------------#
function create_simphy_species_rf_summary(gene_duplication_and_loss_trees::Set{String}, 
                                        gene_loss_only_trees::Set{String})
    """
    Collect RF data between SimPhy gene trees and true species tree 
        for individual replicates.
    
    This function:
    1. Collects RF data from individual replicate files
    2. Classifies gene trees based on events (hidden paralogy, gene loss, nothing)
    3. Returns the collected data for further processing
    
    Inputs:
        gene_duplication_and_loss_trees::Set{String}: 
            Set of gene tree file paths experiencing gene duplication and loss
        gene_loss_only_trees::Set{String}: 
            Set of gene tree file paths experiencing gene loss only

    Returns:
        Vector{Tuple}: Collected RF data as tuples (gene_file_path, rf_score, rf_score_noF,
            gene_duplication_and_loss_flag, gene_loss_only_flag, nothing_flag)
    """
    
    # Collect all RF data from individual rep files
    simphy_species_csv_data = []
    
    println("Looking for RF files in $(n_reps) replicates...")
    for simulation_rep in 1:n_reps
        rep_output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")  
        rf_simphy_true_file = joinpath(rep_output_folder, 
                                      "RF_btw_simphy_vs_true_species_tree_$simulation_rep.txt")
        rf_simphy_true_file_noF = joinpath(rep_output_folder, 
                                      "RF_btw_simphy_vs_true_species_tree_noF_$simulation_rep.txt")
        
        println("Checking RF file: $rf_simphy_true_file")
        
        # Read noF RF scores if available
        rf_noF_dict = Dict{String, Float64}()
        if isfile(rf_simphy_true_file_noF)
            rf_lines_noF = readlines(rf_simphy_true_file_noF)
            for line in rf_lines_noF
                if !isempty(strip(line))
                    parts = split(strip(line), " ")
                    if length(parts) >= 2
                        rf_noF_dict[parts[1]] = parse(Float64, parts[2])
                    end
                end
            end
        end
        
        if isfile(rf_simphy_true_file)
            rf_lines = readlines(rf_simphy_true_file)
            println("Found $(length(rf_lines)) lines in RF file for rep $simulation_rep")
            
            for line in rf_lines
                if !isempty(strip(line))
                    # Parse line format: "filepath rf_score"
                    parts = split(strip(line), " ")
                    if length(parts) >= 2
                        gene_file_path = parts[1]
                        rf_score = parse(Float64, parts[2])
                        rf_score_noF = get(rf_noF_dict, gene_file_path, NaN)  # Get noF score or NaN
                        
                        # Determine classification based on file path
                        gene_duplication_and_loss_flag = gene_file_path in gene_duplication_and_loss_trees ? "Y" : "N"
                        gene_loss_only_flag = gene_file_path in gene_loss_only_trees ? "Y" : "N" 
                        nothing_flag = (gene_duplication_and_loss_flag == "N" && gene_loss_only_flag == "N") ? "Y" : "N"
                        
                        push!(simphy_species_csv_data, 
                            (gene_file_path, rf_score, rf_score_noF,
                            gene_duplication_and_loss_flag, 
                            gene_loss_only_flag, 
                            nothing_flag))
                    end
                end
            end
        else
            println("Warning: RF file not found: $rf_simphy_true_file")
        end
    end
    
    println("Total collected RF data points: $(length(simphy_species_csv_data))")
    return simphy_species_csv_data
end 

#-----------------------------------------------#
# Summarize RF distance 
# There are two parts for this: 
# summarize rep-level statistics from genetrees_stats_$paramname_root.csv 
# and save it into simulation_$paramname_root_summary.csv 
# summarize parameter-level statistics from genetrees_stats_$paramname_root.csv 
# and simulation_$paramname_root_summary.csv 
# and save it into summary_$paramname_root.csv 
#-----------------------------------------------# 


# First summarize statistics at the replicate level: 

#= In this df_rf csv file, there are many columns include: 
1. gene_tree_file_path: Save all genetree file path across all reps in this parameter set
2. num_taxa_gene_tree: Number of taxa in each gene tree  
3. avg_internal_branch_length_locus_tree: Average internal branch length for each locus tree 
4. avg_internal_branch_length_gene_tree: Average internal branch length for each gene tree
5. RF_true_genetrees_vs_species_tree: RF distance between true gene trees and species tree
6. RF_true_vs_estimated_genetrees: RF distance between true gene trees and estimated gene trees
7. gene_duplication_and_loss_or_not: Whether the gene tree underwent duplication and loss
8. gene_loss_only_or_not: Whether the gene tree underwent only loss
9. nothing_or_not: Whether the gene tree underwent no changes
10. false_HP: Whether the gene tree underwent false hidden paralogy
11. weak_HP: Whether the gene tree underwent weak hidden paralogy
12. strong_HP: Whether the gene tree underwent strong hidden paralogy

What are want to compute in the resulting summary csv file: 
This will be append to the simulation_$paramname_root_summary.csv file 

# Anything related to number of taxa:    
1. mean_num_taxa_all_gene_trees: Mean number of taxa for all gene trees in each rep 
2. mean_num_taxa_loss_only_gene_trees: Mean number of taxa for gene loss only genes
3. mean_num_taxa_gene_duplication_and_loss_gene_trees: Mean number of taxa for gene duplication and loss genes
4. mean_num_taxa_nothing_gene_trees: Mean number of taxa for genes with no changes
5. mean_num_taxa_false_HP: Mean number of taxa for false hidden paralogy genes
6. mean_num_taxa_weak_HP: Mean number of taxa for weak hidden paralogy genes
7. mean_num_taxa_strong_HP: Mean number of taxa for strong hidden paralogy genes 

# Anything about internal branch length: 
8. mean_avg_internal_branch_length_locus_trees_all_genes:
        -> Mean average internal branch length for locus trees in each rep
9. mean_avg_internal_branch_length_gene_trees_all_genes:
        -> Mean average internal branch length for gene trees in each rep
10. mean_avg_internal_branch_length_locus_trees_loss_only_gene_trees:
        -> Mean average internal branch length for locus trees for gene loss only genes
11. mean_avg_internal_branch_length_gene_trees_loss_only_gene_trees:
        -> Mean average internal branch length for gene trees for gene loss only genes
12. mean_avg_internal_branch_length_locus_trees_gene_duplication_and_loss_gene_trees:
    -> Mean average internal branch length for locus trees for gene duplication and loss genes
13. mean_avg_internal_branch_length_gene_trees_gene_duplication_and_loss_gene_trees:
    -> Mean average internal branch length for gene trees for gene duplication and loss genes
14. mean_avg_internal_branch_length_locus_trees_nothing_gene_trees:
    -> Mean average internal branch length for locus trees for genes with no changes
15. mean_avg_internal_branch_length_gene_trees_nothing_gene_trees:
    -> Mean average internal branch length for gene trees for genes with no changes
16. mean_avg_internal_branch_length_locus_trees_false_HP:
    -> Mean average internal branch length for locus trees for false hidden paralogy genes
17. mean_avg_internal_branch_length_gene_trees_false_HP:
    -> Mean average internal branch length for gene trees for false hidden paralogy genes
18. mean_avg_internal_branch_length_locus_trees_weak_HP:
    -> Mean average internal branch length for locus trees for weak hidden paralogy genes
19. mean_avg_internal_branch_length_gene_trees_weak_HP:
    -> Mean average internal branch length for gene trees for weak hidden paralogy genes
20. mean_avg_internal_branch_length_locus_trees_strong_HP:
    -> Mean average internal branch length for locus trees for strong hidden paralogy genes
21. mean_avg_internal_branch_length_gene_trees_strong_HP:
    -> Mean average internal branch length for gene trees for strong hidden paralogy genes

# RF distance between true genetrees vs species tree 
22. mean_RF_true_genetrees_vs_species_tree_all_genes: 
        -> Mean RF distance between true gene trees and species tree across all genes in each rep
23. mean_RF_true_genetrees_vs_species_tree_loss_only_gene_trees: 
        -> Mean RF distance between true gene trees and species tree for gene loss only genes
24. mean_RF_true_genetrees_vs_species_tree_gene_duplication_and_loss_gene_trees: 
    -> Mean RF distance between true gene trees and species tree for gene duplication and loss genes
25. mean_RF_true_genetrees_vs_species_tree_nothing_gene_trees:
    -> Mean RF distance between true gene trees and species tree for genes with no changes
26. mean_RF_true_genetrees_vs_species_tree_false_HP:
    -> Mean RF distance between true gene trees and species tree for false hidden paralogy genes
27. mean_RF_true_genetrees_vs_species_tree_weak_HP:
    -> Mean RF distance between true gene trees and species tree for weak hidden paralogy genes
28. mean_RF_true_genetrees_vs_species_tree_strong_HP:
    -> Mean RF distance between true gene trees and species tree for strong hidden paralogy genes

# RF distance between true genetrees vs estimated genetrees
29. mean_RF_true_vs_estimated_genetrees_all_genes:
        -> Mean RF distance between true gene trees and estimated genetrees for all genes in each rep
30. mean_RF_true_vs_estimated_genetrees_loss_only_gene_trees:
    -> Mean RF distance between true gene trees and estimated genetrees for gene loss only genes
31. mean_RF_true_vs_estimated_genetrees_gene_duplication_and_loss_gene_trees:
    -> Mean RF distance between true gene trees and estimated genetrees for gene duplication and loss genes
32. mean_RF_true_vs_estimated_genetrees_nothing_gene_trees:
    -> Mean RF distance between true gene trees and estimated genetrees for genes with no changes
33. mean_RF_true_vs_estimated_genetrees_false_HP:
    -> Mean RF distance between true gene trees and estimated genetrees for false hidden paralogy genes
34. mean_RF_true_vs_estimated_genetrees_weak_HP:
    -> Mean RF distance between true gene trees and estimated genetrees for weak hidden paralogy genes
35. mean_RF_true_vs_estimated_genetrees_strong_HP:
    -> Mean RF distance between true gene trees and estimated genetrees for strong hidden paralogy genes

Those new columns will be appended to the simulation_$paramname_root_summary.csv file 
=#
"""
summarize_rf_distances:
    Summarize information from genetrees_stats_$paramname_root.csv. 
Inputs:
    input_csv_file::String: Path to the genetrees_stats CSV file.
    summary_csv_file::String: Path to the existing simulation summary CSV file.
    paramname_root::String: Parameter set identifier. 
"""
function summarize_rf_distances(input_csv_file::String, 
                                summary_csv_file::String,
                                paramname_root::String)
    if !isfile(input_csv_file) || !isfile(summary_csv_file)
        error("CSV files not found: $input_csv_file or $summary_csv_file")
        return
    end

    df_rf = CSV.read(input_csv_file, DataFrame)
    df_summary = CSV.read(summary_csv_file, DataFrame)
    
    # Extract replicate number from gene tree file paths
    df_rf[!, :RepID] = [parse(Int, match(r"rep(\d+)", row.gene_tree_file_path).captures[1]) 
                        for row in eachrow(df_rf)]
    
    # Initialize new columns in summary dataframe
    new_columns = initialize_summary_columns()
    
    # Create all new columns with NaN values
    for col in new_columns
        df_summary[!, col] = fill(NaN, nrow(df_summary))
    end
    
    # Calculate statistics for each replicate
    for rep_id in df_summary.RepID
        rep_data = df_rf[df_rf.RepID .== rep_id, :]
        
        if nrow(rep_data) == 0
            fill_rep_with_missing!(df_summary, rep_id, new_columns)
            continue
        end
        
        # Calculate and store all statistics for this replicate
        calculate_rep_statistics!(df_summary, rep_data, rep_id, new_columns)
    end
    
    # Save updated summary file
    output_file = summary_csv_file
    println("Writing updated summary to: $output_file")
    CSV.write(output_file, df_summary)
    println("Summary statistics successfully appended to $output_file")
end

"""
initialize_summary_columns:
    Create list of new column names to be added to summary dataframe.
"""
function initialize_summary_columns()
    return [
        # Number of taxa columns
        :mean_num_taxa_all_gene_trees,
        :mean_num_taxa_loss_only_gene_trees,
        :mean_num_taxa_gene_duplication_and_loss_gene_trees,
        :mean_num_taxa_nothing_gene_trees,
        :mean_num_taxa_false_HP,
        :mean_num_taxa_weak_HP,
        :mean_num_taxa_strong_HP,
        
        # Internal branch length - locus trees
        :mean_avg_internal_branch_length_locus_trees_all_genes,
        :mean_avg_internal_branch_length_locus_trees_loss_only_gene_trees,
        :mean_avg_internal_branch_length_locus_trees_gene_duplication_and_loss_gene_trees,
        :mean_avg_internal_branch_length_locus_trees_nothing_gene_trees,
        :mean_avg_internal_branch_length_locus_trees_false_HP,
        :mean_avg_internal_branch_length_locus_trees_weak_HP,
        :mean_avg_internal_branch_length_locus_trees_strong_HP,
        
        # Internal branch length - gene trees
        :mean_avg_internal_branch_length_gene_trees_all_genes,
        :mean_avg_internal_branch_length_gene_trees_loss_only_gene_trees,
        :mean_avg_internal_branch_length_gene_trees_gene_duplication_and_loss_gene_trees,
        :mean_avg_internal_branch_length_gene_trees_nothing_gene_trees,
        :mean_avg_internal_branch_length_gene_trees_false_HP,
        :mean_avg_internal_branch_length_gene_trees_weak_HP,
        :mean_avg_internal_branch_length_gene_trees_strong_HP,
        
        # RF distance: true gene trees vs species tree
        :mean_RF_true_genetrees_vs_species_tree_all_genes,
        :mean_RF_true_genetrees_vs_species_tree_loss_only_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_gene_duplication_and_loss_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_nothing_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_false_HP,
        :mean_RF_true_genetrees_vs_species_tree_weak_HP,
        :mean_RF_true_genetrees_vs_species_tree_strong_HP,
        
        # RF distance: true gene trees vs species tree (noF)
        :mean_RF_true_genetrees_vs_species_tree_noF_all_genes,
        :mean_RF_true_genetrees_vs_species_tree_noF_loss_only_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_noF_gene_duplication_and_loss_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_noF_nothing_gene_trees,
        :mean_RF_true_genetrees_vs_species_tree_noF_false_HP,
        :mean_RF_true_genetrees_vs_species_tree_noF_weak_HP,
        :mean_RF_true_genetrees_vs_species_tree_noF_strong_HP,
        
        # RF distance: true vs estimated gene trees
        :mean_RF_true_vs_estimated_genetrees_all_genes,
        :mean_RF_true_vs_estimated_genetrees_loss_only_gene_trees,
        :mean_RF_true_vs_estimated_genetrees_gene_duplication_and_loss_gene_trees,
        :mean_RF_true_vs_estimated_genetrees_nothing_gene_trees,
        :mean_RF_true_vs_estimated_genetrees_false_HP,
        :mean_RF_true_vs_estimated_genetrees_weak_HP,
        :mean_RF_true_vs_estimated_genetrees_strong_HP
    ]
end

"""
fill_rep_with_missing!:
    Fill a replicate row with NaN values for all new columns.
"""
function fill_rep_with_missing!(df::DataFrame, rep_id::Int, columns::Vector{Symbol})
    row_idx = findfirst(df.RepID .== rep_id)
    for col in columns
        df[row_idx, col] = NaN
    end
end

"""
calculate_rep_statistics!:
    Calculate all summary statistics for a single replicate and update the summary dataframe.
"""
function calculate_rep_statistics!(df_summary::DataFrame, 
                                   rep_data::DataFrame, 
                                   rep_id::Int,
                                   columns::Vector{Symbol})
    row_idx = findfirst(df_summary.RepID .== rep_id)
    
    # Filter data by gene tree categories
    loss_only = rep_data[rep_data.gene_loss_only_or_not .== "Y", :]
    dup_and_loss = rep_data[rep_data.gene_duplication_and_loss_or_not .== "Y", :]
    nothing_trees = rep_data[rep_data.nothing_or_not .== "Y", :]
    false_hp = rep_data[rep_data.false_HP .== "Y", :]
    weak_hp = rep_data[rep_data.weak_HP .== "Y", :]
    strong_hp = rep_data[rep_data.strong_HP .== "Y", :]
    
    # Calculate mean number of taxa
    df_summary[row_idx, :mean_num_taxa_all_gene_trees] = mean(rep_data.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_loss_only_gene_trees] = safe_mean(loss_only.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_nothing_gene_trees] = safe_mean(nothing_trees.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_false_HP] = safe_mean(false_hp.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_weak_HP] = safe_mean(weak_hp.num_taxa_gene_tree)
    df_summary[row_idx, :mean_num_taxa_strong_HP] = safe_mean(strong_hp.num_taxa_gene_tree)
    
    # Calculate mean internal branch lengths - locus trees
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_all_genes] = mean(rep_data.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_loss_only_gene_trees] = safe_mean(loss_only.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_nothing_gene_trees] = safe_mean(nothing_trees.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_false_HP] = safe_mean(false_hp.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_weak_HP] = safe_mean(weak_hp.avg_internal_branch_length_locus_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_locus_trees_strong_HP] = safe_mean(strong_hp.avg_internal_branch_length_locus_tree)
    
    # Calculate mean internal branch lengths - gene trees
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_all_genes] = mean(rep_data.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_loss_only_gene_trees] = safe_mean(loss_only.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_nothing_gene_trees] = safe_mean(nothing_trees.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_false_HP] = safe_mean(false_hp.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_weak_HP] = safe_mean(weak_hp.avg_internal_branch_length_gene_tree)
    df_summary[row_idx, :mean_avg_internal_branch_length_gene_trees_strong_HP] = safe_mean(strong_hp.avg_internal_branch_length_gene_tree)
    
    # Calculate mean RF distance: true gene trees vs species tree
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_all_genes] = mean(rep_data.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_loss_only_gene_trees] = safe_mean(loss_only.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_nothing_gene_trees] = safe_mean(nothing_trees.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_false_HP] = safe_mean(false_hp.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_weak_HP] = safe_mean(weak_hp.RF_true_genetrees_vs_species_tree)
    df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_strong_HP] = safe_mean(strong_hp.RF_true_genetrees_vs_species_tree)
    
    # Calculate mean RF distance: true gene trees vs species tree (noF) - if column exists
    if "RF_true_genetrees_vs_species_tree_noF" in names(rep_data)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_all_genes] = mean(skipmissing(rep_data.RF_true_genetrees_vs_species_tree_noF))
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_loss_only_gene_trees] = safe_mean(loss_only.RF_true_genetrees_vs_species_tree_noF)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.RF_true_genetrees_vs_species_tree_noF)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_nothing_gene_trees] = safe_mean(nothing_trees.RF_true_genetrees_vs_species_tree_noF)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_false_HP] = safe_mean(false_hp.RF_true_genetrees_vs_species_tree_noF)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_weak_HP] = safe_mean(weak_hp.RF_true_genetrees_vs_species_tree_noF)
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_strong_HP] = safe_mean(strong_hp.RF_true_genetrees_vs_species_tree_noF)
    else
        # Fill with NaN if column doesn't exist
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_all_genes] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_loss_only_gene_trees] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_gene_duplication_and_loss_gene_trees] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_nothing_gene_trees] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_false_HP] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_weak_HP] = NaN
        df_summary[row_idx, :mean_RF_true_genetrees_vs_species_tree_noF_strong_HP] = NaN
    end
    
    # Calculate mean RF distance: true vs estimated gene trees
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_all_genes] = mean(rep_data.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_loss_only_gene_trees] = safe_mean(loss_only.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_gene_duplication_and_loss_gene_trees] = safe_mean(dup_and_loss.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_nothing_gene_trees] = safe_mean(nothing_trees.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_false_HP] = safe_mean(false_hp.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_weak_HP] = safe_mean(weak_hp.RF_true_vs_estimated_genetrees)
    df_summary[row_idx, :mean_RF_true_vs_estimated_genetrees_strong_HP] = safe_mean(strong_hp.RF_true_vs_estimated_genetrees)
end

"""
safe_mean:
    Calculate mean of a vector, returning NaN if the vector is empty.
"""
function safe_mean(values::AbstractVector)
    if isempty(values)
        return NaN
    end
    return mean(values)
end


#------------------------------------------------#
# second, summarize statistics at the parameter level: 

"""
summarize_parameter_level_RF_statistics:
    Calculate mean statistics across all genes (ignoring replication level) 
    from genetrees_stats CSV file.
    
Inputs:
    paramname_root::String: Parameter set identifier
    outfolder::String: Output folder containing the genetrees_stats CSV
    output_csv_file_path::String: Path to save the parameter-level summary CSV
"""
function summarize_parameter_level_RF_statistics(paramname_root::String, 
                                        outfolder::String, 
                                        output_csv_file_path::String)

    genetree_stats_csv_file = joinpath(outfolder, 
                             "genetrees_stats_$(paramname_root).csv") 

    if !isfile(genetree_stats_csv_file)
        error("Genetree stats CSV file not found: $genetree_stats_csv_file")
        return
    end

    df = CSV.read(genetree_stats_csv_file, DataFrame)
    
    # Read n_genes per replicate from simulation CSV
    simulation_csv_file = joinpath(outfolder, "simulation_$(paramname_root).csv")
    n_genes = CSV.read(simulation_csv_file, DataFrame)[!, :n_genes]

    # Filter data by gene tree categories (across all genes/reps)
    all_genes = df
    loss_only = df[df.gene_loss_only_or_not .== "Y", :]
    dup_and_loss = df[df.gene_duplication_and_loss_or_not .== "Y", :]
    nothing_trees = df[df.nothing_or_not .== "Y", :]
    false_hp = df[df.false_HP .== "Y", :]
    weak_hp = df[df.weak_HP .== "Y", :]
    strong_hp = df[df.strong_HP .== "Y", :]
    
    # Create summary DataFrame with means across all genes
    df_param_summary = DataFrame(
        parameter_setting = [paramname_root],
        n_genes_mean = [mean(n_genes)],
        n_genes_min = [minimum(n_genes)],
        n_genes_max = [maximum(n_genes)],
        
        # Mean number of taxa
        mean_num_taxa_all_genes = [mean(all_genes.num_taxa_gene_tree)],
        mean_num_taxa_loss_only = [safe_mean(loss_only.num_taxa_gene_tree)],
        mean_num_taxa_dup_and_loss = [safe_mean(dup_and_loss.num_taxa_gene_tree)],
        mean_num_taxa_nothing = [safe_mean(nothing_trees.num_taxa_gene_tree)],
        mean_num_taxa_false_HP = [safe_mean(false_hp.num_taxa_gene_tree)],
        mean_num_taxa_weak_HP = [safe_mean(weak_hp.num_taxa_gene_tree)],
        mean_num_taxa_strong_HP = [safe_mean(strong_hp.num_taxa_gene_tree)],
        
        # Mean internal branch lengths - locus trees
        mean_internal_bl_locus_all = [mean(all_genes.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_loss_only = [safe_mean(loss_only.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_dup_and_loss = [safe_mean(dup_and_loss.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_nothing = [safe_mean(nothing_trees.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_false_HP = [safe_mean(false_hp.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_weak_HP = [safe_mean(weak_hp.avg_internal_branch_length_locus_tree)],
        mean_internal_bl_locus_strong_HP = [safe_mean(strong_hp.avg_internal_branch_length_locus_tree)],
        
        # Mean internal branch lengths - gene trees
        mean_internal_bl_gene_all = [mean(all_genes.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_loss_only = [safe_mean(loss_only.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_dup_and_loss = [safe_mean(dup_and_loss.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_nothing = [safe_mean(nothing_trees.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_false_HP = [safe_mean(false_hp.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_weak_HP = [safe_mean(weak_hp.avg_internal_branch_length_gene_tree)],
        mean_internal_bl_gene_strong_HP = [safe_mean(strong_hp.avg_internal_branch_length_gene_tree)],
        
        # Mean RF: true gene trees vs species tree
        mean_RF_genetree_vs_sptree_all = [mean(all_genes.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_loss_only = [safe_mean(loss_only.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_dup_and_loss = [safe_mean(dup_and_loss.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_nothing = [safe_mean(nothing_trees.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_false_HP = [safe_mean(false_hp.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_weak_HP = [safe_mean(weak_hp.RF_true_genetrees_vs_species_tree)],
        mean_RF_genetree_vs_sptree_strong_HP = [safe_mean(strong_hp.RF_true_genetrees_vs_species_tree)],
        
        # Mean RF: true gene trees vs species tree (noF) - if column exists
        mean_RF_genetree_vs_sptree_noF_all = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [mean(skipmissing(all_genes.RF_true_genetrees_vs_species_tree_noF))] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_loss_only = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(loss_only.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_dup_and_loss = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(dup_and_loss.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_nothing = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(nothing_trees.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_false_HP = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(false_hp.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_weak_HP = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(weak_hp.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        mean_RF_genetree_vs_sptree_noF_strong_HP = ("RF_true_genetrees_vs_species_tree_noF" in names(all_genes) ? [safe_mean(strong_hp.RF_true_genetrees_vs_species_tree_noF)] : [NaN]),
        
        # Mean RF: true vs estimated gene trees
        mean_RF_true_vs_est_genetree_all = [mean(all_genes.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_loss_only = [safe_mean(loss_only.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_dup_and_loss = [safe_mean(dup_and_loss.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_nothing = [safe_mean(nothing_trees.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_false_HP = [safe_mean(false_hp.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_weak_HP = [safe_mean(weak_hp.RF_true_vs_estimated_genetrees)],
        mean_RF_true_vs_est_genetree_strong_HP = [safe_mean(strong_hp.RF_true_vs_estimated_genetrees)]
    )
    
    # Read the file created by extract_summarize_simphy_summary and add RF columns
    df = CSV.read(output_csv_file_path, DataFrame)
    row_idx = findfirst(df.parameter_setting .== paramname_root)
    
    if row_idx !== nothing
        # Update the row with new RF statistics
        for col in names(df_param_summary)
            if col != "parameter_setting"
                if col in names(df)
                    df[row_idx, col] = df_param_summary[1, col]
                else
                    df[!, col] = fill(NaN, nrow(df))
                    df[row_idx, col] = df_param_summary[1, col]
                end
            end
        end
    end
    
    CSV.write(output_csv_file_path, df)
    println("Parameter-level RF statistics appended to: $output_csv_file_path") 
end



"""
extract_simphy_summary:
    Extract and summarize statistics from simulation CSV file for a parameter setting.
    Calculates means across replicates and percentages for gene tree categories.

Inputs:
    paramname_root::String: Parameter set identifier (e.g., "DUP0.0001-LOS0.0001-RVL-N_ind1-SF1.0")
    input_dir::String: Directory containing the simulation CSV file
    output_csv_file_path::String: Path where the summary should be written

Outputs:
    Named tuple containing summary statistics including:
    - Mean values for n_genes, n_iterations, taxa removed
    - Percentages of genes meeting minimum requirements
    - Percentages of trees experiencing different events (loss only, duplication+loss, nothing)
    - Percentages of hidden paralogy types (false, weak, strong)
    - Mean RF distance between true and estimated species trees
    - Count of replicates where estimated species tree differs from truth
"""
function extract_summarize_simphy_summary(paramname_root::String,
                                outfolder::String, 
                                output_csv_file_path::String)

    csv_path = joinpath(outfolder, "simulation_$(paramname_root).csv")
    df = CSV.read(csv_path, DataFrame)
    
    # Extract column data
    n_genes = df[!,:n_genes] 
    n_iterations = df[!,:n_iterations]
    n_repeated_taxa_removed = df[!,:n_repeated_taxa_removed] 
    n_insufficient_taxa_removed = df[!,:n_insufficient_taxa_removed]

    # Gene tree event categories
    num_trees_experiencing_gene_loss_only = df[!,:num_trees_experiencing_gene_loss_only] 
    num_trees_experiencing_gene_duplication_and_loss = df[!,:num_trees_experiencing_gene_duplication_and_loss]
    num_trees_experiencing_nothing = df[!,:num_trees_experiencing_nothing] 
    num_false_HP = "num_false_HP" in names(df) ? df[!, :num_false_HP] : zeros(nrow(df))
    num_weak_HP = "num_weak_HP" in names(df) ? df[!, :num_weak_HP] : zeros(nrow(df))
    num_strong_HP = "num_strong_HP" in names(df) ? df[!, :num_strong_HP] : zeros(nrow(df))

    # Count genes meeting minimum requirements
    bool_min_genes = df[!, Symbol(">=min_num_genes")]
    num_T = count(x -> x == 1, bool_min_genes)
    num_F = count(x -> x == 0, bool_min_genes)

    # RF distance statistics
    RF = df[!, :RF_astral_true]
    num_estimated_species_tree_diff_from_truth = count(x -> !isnan(x) && x != 0, RF)

    # Calculate total genes across all reps
    total_n_genes = sum(n_genes) 

    # Create summary statistics as DataFrame
    df_summary = DataFrame(
        parameter_setting = [paramname_root],
        n_genes_mean = [mean(n_genes)],
        n_genes_min = [minimum(n_genes)], # find the min n_genes across reps in the setting 
        n_genes_max = [maximum(n_genes)], # find the max n_genes across reps in the setting 
        n_iterations = [mean(n_iterations)],
        n_repeated_taxa_removed = [mean(n_repeated_taxa_removed)],
        n_insufficient_taxa_removed = [mean(n_insufficient_taxa_removed)],
        percentage_genes_meet_min = [num_T / (num_T + num_F) * 100],
        
        # Percentages of trees experiencing different events
        pert_trees_experiencing_gene_loss_only = [sum(num_trees_experiencing_gene_loss_only) / total_n_genes * 100],
        pert_trees_experiencing_gene_duplication_and_loss = [sum(num_trees_experiencing_gene_duplication_and_loss) / total_n_genes * 100],
        pert_trees_experiencing_nothing = [sum(num_trees_experiencing_nothing) / total_n_genes * 100],
        
        # Percentages of hidden paralogy types
        pert_false_HP = [sum(num_false_HP) / total_n_genes * 100],
        pert_weak_HP = [sum(num_weak_HP) / total_n_genes * 100],
        pert_strong_HP = [sum(num_strong_HP) / total_n_genes * 100],
        
        # RF distance metrics
        mean_RF_true_and_estimated_species_trees = [mean(RF)],
        num_estimated_species_tree_diff_from_truth = [num_estimated_species_tree_diff_from_truth]
    )
    
    # Write to CSV file
    CSV.write(output_csv_file_path, df_summary)
    println("Summary statistics written to: $output_csv_file_path")
end 

#------------------------------------------------#
# Categorize gene trees with hidden paralogy
#------------------------------------------------# 
#=
There are a few senarioes for gene trees with hidden paralogy. 
  Our deinition for generating classifying that a gene going through hidden paralogy is:
  simphy this gene goes through duplication and loss, and there are no multiple genes 
  copies from the same accession is observed in the gene tree. 
  However, there are three other senarioes for hidden paralogy: 
    
  False hidden paralogy: 
  There are two senarioes for false HP: 
    a. A gene got duplicated, but the duplicated gene copy quickly got lost. 
    In this case, no hidden paralogy is happened. However, technically gene loss and 
    duplication did happen. For example, A_0 got copied into A_0 and A_1, and then 
    A_1 got lost so that A_0 is left. In this case, I will call it "false" hidden paralogy.
    In this case, the suffixes of all the tips in the locus tree are the same, 
    such as A_0, B_0, C_0, D_0, E_0, F_0, G_0, H_0. 
    b. The suffixes of all the tips in the locus tree are not the same but the gene loss and 
    duplication process doesn't change the topology or the branch length of the locus tree. 
    In this case, we used distance method to figure out if the branch length of the locus tree 
    change compared to the true species tree. 
  
  Weak hidden paralogy:
    2. A gene got duplicated, and then the another copy of gene copy got lost. For example, 
    A_0 got copied and became A_0 and A_1, and then the original gene copy A_0 got lost. However, 
    the topology of the locus tree and species tree are the same. 
    In this case, I will call it "weak" hidden paralogy.
  
  Strong HP: 
    3. After the previous two senarioes, a gene tree has different topology from the species tree.
    In this case, I will call it "strong" hidden paralogy. 
This requires to examine the locus tree and species tree. 

=# 
# Individual worker function for processing a single gene tree path
@everywhere begin 
  """
  process_single_HP_line:
  Process a single gene tree file path to classify hidden paralogy.

  Inputs:
    line::String: File path to the gene tree.
    n_reps::Int: Total number of replicates.
    folder_path_list::Vector: List of paths to each replicate folder.
    species_tree::String: Newick string of the true species tree.
  Outputs:
    Tuple{String, String, Float64}: A tuple containing:
      - The original file path.
      - The hidden paralogy category ("false", "weak", "strong", or "skip").
      - The RF distance
  """
    function process_single_HP_line(line::String, 
                                   n_reps::Int, 
                                   folder_path_list::Vector, 
                                   delta = 1e-4
                                   )
        gene_tree_name = split(line, "/")[end] # find the last item in the path 
        
        # Extract rep_id from the path (e.g., "rep01" from the path)
        path_parts = split(line, "/")
        rep_part = ""
        for part in path_parts
            if startswith(part, "rep")
                rep_part = part
                break
            end
        end
        simulation_rep = parse(Int, replace(rep_part, "rep" => "")) # Extract number from "rep01" -> 1
        
        # Skip if this replicate is not in our current processing range
        if simulation_rep > n_reps
            return (line, "skip", NaN)
        end
        rep_id = pad_number(simulation_rep, n_reps) # pad with leading zeros if necessary

        # find gene and int numbers from the filename 
        # example gene tree name: g_trees_noLocusID_Gene4_Int6.trees 
        gene_part = split(gene_tree_name, "_")[4] # "Gene4"
        gene_id = parse(Int, replace(gene_part, "Gene" => "")) # Extract number: 4
        int_part = split(split(gene_tree_name, "_")[5], ".")[1] # "Int6" -> "Int6"
        int_id = parse(Int, replace(int_part, "Int" => "")) # Extract number: 6
        # This identifies to the correct replicate folder, e.g., "rep01"
        # gene id, such as Gene4 
        # int id, such as Int6 

        # Find the corresponding locus tree file 
        simphy_raw_gene_tree_folder = setup_rep_output_folders(folder_path_list, 
                                                              simulation_rep, 
                                                              "genetrees_simphy") 

        # Our legacy code structure include "1" subfolder for locus trees 
        # New structure may not have it, so we check both 
        temp_dir = joinpath(simphy_raw_gene_tree_folder, "Int$(int_id)", "1") 
        if isdir(temp_dir)
            locus_tree_file = joinpath(temp_dir, "l_trees.trees") 
        else 
            locus_tree_file = joinpath(simphy_raw_gene_tree_folder, 
                                    "Int$(int_id)", "l_trees.trees") 
        end 
        
        locus_trees = readlines(locus_tree_file) 
        locus_tree = readnewick(locus_trees[gene_id]) # read the corresponding locus tree 
        
        # senario 1: false hidden paralogy 
        # check if the locus tree has all tips ending with the same suffix 
        tip_labels = tipLabels(locus_tree) 
        suffixes = Set{String}()

        for tip in tip_labels
            if startswith(tip, "Lost-") 
                # This means this tip is lost, so we skip it 
                deleteleaf!(locus_tree, tip) # prune the lost tip from the locus tree 
                continue 
            else 
                suffix = split(tip, "_")[end] # get the suffix after the last underscore
                push!(suffixes, suffix)
            end
        end

        # The first senarios for false hidden paralogy: check if all suffixes are the same
        if length(suffixes) == 1 # false hidden paralogy 
            return (line, "false", 0.0)

        else 
            true_species_tree = readnewick("(A:3440,((((B:880,C:880):1710,(D:930,E:930):1660):170,F:2760):180,(G:500,H:500):2440):500);") 
            # prune tips on locus tree so that it first prunes all tips starting with "lost-"
            # This has been done above
            # trim anything after "_" in the tip names
            modify_locus_tree_labels!(locus_tree) # see utilities.jl for details

            pruned_locus_tree, pruned_species_tree = 
                prune_nonoverlapping_tips(locus_tree, true_species_tree) 

            # have a de-bug checking file 
            tiplabels_locus = tipLabels(pruned_locus_tree) # will be used later 
            sorted_tiplabels_locus = sort(tiplabels_locus)
            tiplabels_species = tipLabels(pruned_species_tree) # will be used later 
            sorted_tiplabels_species = sort(tiplabels_species) 

            if sorted_tiplabels_locus != sorted_tiplabels_species
                error("Tip labels differ between locus and species trees for $gene_tree_name in rep $simulation_rep")
            end 

            rf_distance = mudistance_semidirected(pruned_locus_tree, pruned_species_tree, preorder=true)

            if rf_distance == 0 
                # weak hidden paralogy or false hidden paralogy? 
                # Even if the rf_distance = 0, it could be false or weak hidden paralogy 
                # The difference is whether the branch length of the locus tree is changed compared to the species tree 
                # In changed, then it is weak hidden paralogy
                # If not, it is still false hidden paralogy: 
                
                # First, our Lineage variation is applied to scale the species tree branch 
                # but it does not change the branch length of the species tree. 
                m_locus = pairwisetaxondistancematrix(pruned_locus_tree)
                m_species = pairwisetaxondistancematrix(pruned_species_tree) 
                taxalist_locus = sortperm(tiplabels_locus) 
                taxalist_species = sortperm(tiplabels_species) 
                sorted_m_locus = m_locus[taxalist_locus, taxalist_locus]
                sorted_m_species = m_species[taxalist_species, taxalist_species] 

                if isapprox(sorted_m_locus, sorted_m_species; atol=delta, rtol=0.0)
                    return (line, "false", rf_distance)
                else
                    return (line, "weak", rf_distance)
                end

            else 
                # strong hidden paralogy 
                # If rf distance > 0 this must be strong hidden paralogy 
                return (line, "strong", rf_distance)
            end 
        end
    end
end 

# Main function that coordinates the distributed processing
"""
Post-examine hidden paralogy based on locus tree information.
Inputs:
    outfolder::String: Output folder containing simulation results.
    n_reps::Int: Number of replicates.
    paramname_root::String: Parameter set identifier.
    folder_path_list::Vector: List of paths to each replicate folder.
    species_tree::String: Newick string of the true species tree (default is "(A,((((B,C),(D,E)),F),(G,H)));").
Outputs:
    Tuple{Vector{String}, Vector{String}, Vector{String}}: Three vectors containing file paths for:
        - False hidden paralogy gene trees.
        - Weak hidden paralogy gene trees.
        - Strong hidden paralogy gene trees.
"""
function postexamine_HP_based_on_locus_tree(outfolder::String,
                                    n_reps::Int, 
                                    paramname_root::String, 
                                    folder_path_list::Vector, 
                                    species_tree::String = "(A,((((B,C),(D,E)),F),(G,H)));"
                                    ) 
  HP_file_path = joinpath(outfolder, "Simphy_gene_duplication_and_loss_$paramname_root.csv")
  
  # Read the CSV file properly to extract just the gene tree file paths
  if isfile(HP_file_path)
      hp_df = CSV.read(HP_file_path, DataFrame)
      data_lines = hp_df[!, "gene_tree_file_path"]  # Extract just the file paths
  else
      error("HP file not found: $HP_file_path")
  end
  
  # Process all lines in parallel using pmap
  results = pmap(line -> process_single_HP_line(line, 
                          n_reps, 
                          folder_path_list), 
                          data_lines)
  
  # Collect results by category
  false_HP_paths = String[]
  weak_HP_paths = String[]
  strong_HP_paths = String[]
  
  for (line, category, rf_dist) in results
    if category == "false"
      push!(false_HP_paths, line)
    elseif category == "weak"
      push!(weak_HP_paths, line)
    elseif category == "strong"
      push!(strong_HP_paths, line)
    end
  end
  
  # Return the results for further processing if needed
  return false_HP_paths, weak_HP_paths, strong_HP_paths
end 

#-----------------------------------------------#
# Calculate RF distance between 
# estimated species tree from astral vs true species tree
#-----------------------------------------------# 
@everywhere begin
  """
  compare_astral_truespecies_1rep:
  Compare the estimated species tree from ASTRAL with the true species 
  tree for a given replicate.

  Inputs:
    simulation_rep::Int: Replication ID.
    folder_path_list::Vector: List of paths to each replicate folder.
    rootfolder::String: Root folder of the simulation pipeline.

  Outputs:
    Returns the Robinson-Foulds (RF) distance between the true species tree 
    and the estimated species tree.

  Notes:
    - The true species tree is hardcoded as "(A,((((B,C),(D,E)),F),(G,H)));".
    - If n_inds = 1, the gene trees are treated as single-copy.
        - Both trees are rooted at node "A" before comparison.
    - If n_inds > 1, the gene trees are treated as multi-copy. 
        - Both trees are rooted at node "A" before comparison. 
    - The RF distance is calculated using the `mudistance_semidirected` 
      function from PhyloNetworks.
  """
  function compare_astral_simphy_truespecies_1rep(
                    simulation_rep::Int, 
                    folder_path_list::Vector, 
                    rootfolder::String,
                    n_inds::Int
                    )
    output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")

    # The true species tree VS. estimated species tree from astral 
    astralfolder = joinpath(output_folder, "astralfolder")
    astralfile = joinpath(astralfolder, "astral.tre") 
    
    # Read and clean the ASTRAL tree file to handle nan values
    local astraltre

    # Error out immediately if ASTRAL tree file is missing
    if !isfile(astralfile)
        error("ASTRAL tree file not found for replicate $simulation_rep: $astralfile")
    end

    astral_tree_string = read(astralfile, String)
    # Replace nan values with 0.0 to make the tree parseable
    cleaned_tree_string = replace(astral_tree_string, r":nan\b" => ":0.0")
    cleaned_tree_string = replace(cleaned_tree_string, r"\bnan:" => "0.0:")
    cleaned_tree_string = replace(cleaned_tree_string, r"\bnan\b" => "0.0")
    astraltre = readnewick(cleaned_tree_string) # estimated astral species tree
    true_species = readnewick("(A,((((B,C),(D,E)),F),(G,H)));") # true species tree 

    # The astral tree might not have matching tips with the true species tree
    # although this should be very rare because with n_genes = 1000, 
    # unless the duplication/loss rates are extremely high 
    # it is hard to completely loss a species in the estimated species tree 
    # However just to double check, we prune the non-overlapping tips 
    pruned_true_species, pruned_astraltre = 
        prune_nonoverlapping_tips(true_species, astraltre) 
    rf_btw_estimated_and_true_tree = mudistance_semidirected(pruned_true_species, 
                                                            pruned_astraltre, 
                                                            preorder=true) # check if the function works
    # For trees, this is equivalent to the Robinson-Foulds distance 

    # Normalize the RF distance 
    n_true_species_tree = length(tipLabels(pruned_true_species)) 
    n_astral_tree = length(tipLabels(pruned_astraltre))
    # Again, after pruning, the two trees should have the same number of leaves 
    # Let's double check 
    if n_true_species_tree != n_astral_tree
        error("Number of leaves differ between true species ($n_true_species_tree) and ASTRAL ($n_astral_tree) in replicate $simulation_rep")
    end 

    rf_btw_estimated_and_true_tree /= (2 * (n_true_species_tree - 3)) # normalized RF distance

    # The true species tree VS. gene trees from simphy (true gene trees)
    simphyfolder = setup_rep_output_folders(folder_path_list, simulation_rep, "genetrees_singlecopy") 
    astralfolder = setup_rep_output_folders(folder_path_list, simulation_rep, "astralfolder") 
    mappingfile = joinpath(astralfolder, "astral_mapping.txt")
    
    # Error out immediately if mapping file or simphy folder is missing
    if !isfile(mappingfile)
        error("Mapping file not found for replicate $simulation_rep: $mappingfile")
    end
    if !isdir(simphyfolder)
        error("SimPhy folder not found for replicate $simulation_rep: $simphyfolder")
    end
    
    # Use the new mapping-aware RF distance calculation
    rf_result = calculate_rf_distance_simphy_vs_species(
        "(A,((((B,C),(D,E)),F),(G,H)));",  # true species tree
        simphyfolder,
        mappingfile,
        n_inds,
        true  # remove_F_analysis = true
    )
    
    # Unpack results based on whether noF analysis was done
    if length(rf_result) == 3
        gene_tree_files_with_paths, rf_btw_simphy_genetree_and_true_tree, rf_btw_simphy_genetree_and_true_tree_noF = rf_result
    else
        gene_tree_files_with_paths, rf_btw_simphy_genetree_and_true_tree = rf_result
        rf_btw_simphy_genetree_and_true_tree_noF = Float64[]  # Empty if not calculated
    end

    return rf_btw_estimated_and_true_tree, gene_tree_files_with_paths, rf_btw_simphy_genetree_and_true_tree, rf_btw_simphy_genetree_and_true_tree_noF
  end
end


#-----------------------------------------------#
# Calculate RF distance between 
# gene trees from simphy vs gene trees estimated by IQ-TREE  
# This gave an estimation of how well IQ-TREE recovers gene trees 
# from simulated molecular sequences
#-----------------------------------------------# 
@everywhere begin
  """
  Compare gene trees generated by SimPhy (true gene trees) and 
  gene trees estimated by IQ-TREE for a given replicate.

  This function uses the mapping.csv file in the IQ-TREE output folder 
  to match each estimated gene tree to its corresponding true gene tree 
  from SimPhy.
  Inputs:
    - simulation_rep::Int: Replication ID (1-based index for the replicate).
    - folder_path_list::Vector: List of output folder paths for each replicate.
    - rootfolder::String: Root directory of the simulation pipeline.
  Output:
    - Tuple{Vector{Float64}, Vector{String}}: A tuple containing:
      - RF distances between each true gene tree (SimPhy) and its corresponding estimated gene tree (IQ-TREE) for the replicate
      - File paths to the corresponding SimPhy gene tree files
  """
  function compare_gene_trees_iqtree_simphy_1rep(
                    simulation_rep::Int, 
                    folder_path_list::Vector, 
                    rootfolder::String
                    )

    output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")
    iqtreefolder = joinpath(output_folder, "iqtreefolder")

    # Check if IQ-TREE folder and mapping file exist
    # Read the mapping file to get correspondence between alignment files and IQ-TREE trees
    mapping_file = joinpath(iqtreefolder, "mapping.csv") 
    mapping_lines = readlines(mapping_file)[2:end] # skip header
    # Extract gene names and IQ-TREE newick strings from mapping file
    gene_names = String[]
    iqtree_newick_strings = String[]
    
    for line in mapping_lines
        # Parse each line: phylip_file : "newick_string"
        if occursin(".phy", line) && occursin(":", line)
            # Extract gene name from path (e.g., "g_trees_noLocusID_Gene01_Int0.phy")
            gene_match = match(r"(g_trees_noLocusID_Gene\d+_Int\d+)\.phy", line)
            # Extract newick string (everything after " : " and within quotes)
            newick_match = match(r":\s*\"([^\"]+)\"", line)
            
            if gene_match !== nothing && newick_match !== nothing
                gene_name = gene_match.captures[1]
                newick_string = newick_match.captures[1]
                push!(gene_names, gene_name)
                push!(iqtree_newick_strings, newick_string)
            end
        end
    end

    # Read corresponding SimPhy gene trees and calculate RF distances
    singlecopy_folder = setup_rep_output_folders(folder_path_list, 
                                                simulation_rep, 
                                                "genetrees_singlecopy")
    
    rf_distances = Float64[]
    gene_file_paths = String[]
    
    for (i, gene_name) in enumerate(gene_names)
        # Read the corresponding SimPhy gene tree
        simphy_tree_file = joinpath(singlecopy_folder, gene_name * ".trees")
        
        if isfile(simphy_tree_file)
            try
                simphy_tree = readnewick(simphy_tree_file)
                # Parse IQ-TREE newick string from mapping file
                iqtree_tree = readnewick(iqtree_newick_strings[i])
                
                # Calculate RF distance between SimPhy and IQ-TREE gene trees
                rf_distance = mudistance_semidirected(simphy_tree, iqtree_tree, preorder=true)
                
                n_simphytree = length(tipLabels(simphy_tree)) # number of leaves in the gene tree
                n_iqtreetree = length(tipLabels(iqtree_tree)) # number of leaves in the gene tree  

                # We don't need to prune the non-overlapping tips here 
                # because both are gene trees so their tips should match, 
                # But the below function just tries to double check 
                if n_simphytree != n_iqtreetree
                    error("Warning: Number of leaves differ between SimPhy ($n_simphytree) and IQ-TREE ($n_iqtreetree) for gene $gene_name in replicate $simulation_rep")
                end

                # Calculate normalized RF distance
                # Normalized RF distance = RF distance / (2 * (n - 3)) 
                normalized_rf_distance = rf_distance / (2 * (n_simphytree - 3))

                println("Normalized RF distance for gene $gene_name:\n$normalized_rf_distance")

                push!(rf_distances, normalized_rf_distance)
                push!(gene_file_paths, simphy_tree_file)
            catch e
                println("Warning: Failed to parse tree files for gene $gene_name: $e")
                push!(rf_distances, NaN)
                push!(gene_file_paths, simphy_tree_file)
            end
        else
            println("Warning: SimPhy tree file not found: $simphy_tree_file")
            push!(rf_distances, NaN)
            push!(gene_file_paths, simphy_tree_file)
        end
    end

    return rf_distances, gene_file_paths
  end 
end

#=
The below code calculates the RF distance between simphy gene tree
and the species tree, which is the RF distance between the true gene 
tree and the true species tree. 
=# 
@everywhere begin
    """
        calculate_rf_distance_simphy_vs_species(species_tree_string::String, 
                                            simphy_folder::String, 
                                            mapping_file::String)

    Calculate RF distances between SimPhy-generated gene trees and the true species tree
    using individual-to-species mapping.

    This function reads SimPhy gene trees, maps individuals to species using the astral
    mapping file format, and calculates Robinson-Foulds (RF) distances between each 
    mapped gene tree and the true species tree.

    # Arguments
    - `species_tree_string::String`: Newick string of the true species tree
    - `simphy_folder::String`: Path to folder containing SimPhy gene tree files  
    - `mapping_file::String`: Path to astral_mapping.txt file with 
        individual-to-species mapping

    # Returns
    - `Vector{Float64}`: RF distances between each gene tree and the species tree

    # Notes
    - Expects gene tree files with names matching pattern `g_trees_noLocusID_Gene*.trees`
    - Mapping file should have format: `individual_name species_name` (space-separated)
    - Uses PhyloNetworks.mapindividuals to perform the mapping
    - Skips gene trees that cannot be successfully mapped or processed

    # Examples
    ```julia-repl
    julia> rf_distances = calculate_rf_distance_simphy_vs_species(
            "(A,((((B,C),(D,E)),F),(G,H)));",
            "output/rep1/genetrees_singlecopy/",
            "output/rep1/astralfolder/astral_mapping.txt"
        )
    [0.0, 2.0, 0.0, 4.0, ...]
    ```
    """
    function calculate_rf_distance_simphy_vs_species(species_tree_string::String,
                                                simphy_folder::String,
                                                mapping_file::String, 
                                                n_inds::Int,
                                                remove_F_analysis::Bool = True)
        #= if remove_F_analysis = True (default), 
            Then, we want to remove "F" from the gene tree and the species tree 
            and then re-calculate RF distance. 
            This will be included in the summary_csv as 
            mean_RF_genetreeNoF_vs_sptreeNoF_[category]
        =# 
        # Read the true species tree
        species_tree = readnewick(species_tree_string)

        if remove_F_analysis
            species_tree_noF = deepcopy(species_tree)
            deleteleaf!(species_tree_noF, "F")
        end 

        individual_to_species = create_mapping_file_based_on_gene_tree(mapping_file) # utilities.jl 

        rf_distances = Float64[]
        rf_distances_noF = Float64[]  # Initialize noF distances list
        gene_tree_files_with_paths = Tuple{String, String}[]  # (filename, full_path)
        
        # Get all gene tree files (preserve original order)
        gene_tree_files = filter(f -> startswith(f, "g_trees_noLocusID_Gene") && endswith(f, ".trees"), readdir(simphy_folder))
        
        for gene_file in gene_tree_files
            
            # Read the gene tree
            gene_tree_path = joinpath(simphy_folder, gene_file)
            gene_tree_content = read(gene_tree_path, String) 

            if n_inds > 1 && !isempty(gene_tree_content) 
                
                # This list stores the mean RF distance across all subtrees
                # extracted from one gene tree file 
                subtree_normalized_rf_distances = Float64[] 
                subtree_normalized_rf_distances_noF = Float64[]  # For noF analysis

                gene_tree = readnewick(gene_tree_content)  
                sisterdict = get_sisters_for_all_leafs(gene_tree) # get sisters for all leaves 
                prune_related_leaves!(gene_tree, sisterdict) # prune related leaves 
            
                # pointer --> too tired. Will work on this later 
                # can use simplify_tip_labels(pruned_tree) to simplify tip labels 

                # Now for any trees with duplicated species, for example, "A_0" and "A_1", 
                # both occur in the pruned tree, 
                # we want to extract two sub-trees: one with "A_0" and the other with "A_1". 
                # This needs to be done for all combinations of duplicated species. 
                # In theory, with 8 taxa, we gonna have 2^8 = 256 combinations at most
                # but we pruned sister taxa within the same species first, 
                # so this number should be smaller. 

                species_dict = group_tips_by_species(gene_tree)
                # Above creates a dictionary where keys are species names (e.g., "A") 
                # and values are vectors of tip names belonging to that species (e.g., ["A_0", "A_1"]) =# 
                
                all_comb = all_combinations(species_dict) 
                # Above gives all the combinations of tips that will be kept 

                # Iterate through all combinations and extract sub-trees 
                for combo in all_comb
                    tips_to_keep = collect(combo)
                    subtree = subtree_for_combination(gene_tree, tips_to_keep)
                    # for each subtree, we modify the tip so that "A_0" become "A" for example
                    mapped_subgene_tree = map_gene_tree_based_on_mapping_file(subtree,
                                                        individual_to_species)

                    # If remove_F_analysis, create a copy without F
                    if remove_F_analysis
                        mapped_subgene_tree_noF = deepcopy(mapped_subgene_tree)
                        # Check if F exists before trying to delete
                        if "F" in tipLabels(mapped_subgene_tree_noF)
                            deleteleaf!(mapped_subgene_tree_noF, "F")
                        end
                    end

                    # Calculate RF distance
                    pruned_species_tree, pruned_mapped_subgene_tree =
                        prune_nonoverlapping_tips(species_tree, mapped_subgene_tree)

                    rf_dist = mudistance_semidirected(pruned_species_tree,
                                                    pruned_mapped_subgene_tree,
                                                    preorder=true)

                    # If remove_F_analysis, also calculate RF distance without F
                    if remove_F_analysis
                        pruned_species_tree_noF, pruned_mapped_subgene_tree_noF =
                            prune_nonoverlapping_tips(species_tree_noF, mapped_subgene_tree_noF)
                        rf_dist_noF = mudistance_semidirected(pruned_species_tree_noF,
                                                        pruned_mapped_subgene_tree_noF,
                                                        preorder=true)
                    end

                    # Normalized RF distance -> RF distance / (2 * (n - 3)) 
                    n = length(tipLabels(pruned_species_tree)) 
                    # This calculates RF distance for each subtree 
                    normalized_rf_dist_per_subtree = rf_dist / (2 * (n - 3))
                    
                    if remove_F_analysis
                        n_noF = length(tipLabels(pruned_species_tree_noF))
                        normalized_rf_dist_per_subtree_noF = rf_dist_noF / (2 * (n_noF - 3))
                    end
                    
                    # For each gene_tree file, we may have multiple subtrees 
                    # so we want to have mean RF distance across all subtrees
                    # as the RF distance for this gene tree file
                    push!(subtree_normalized_rf_distances, normalized_rf_dist_per_subtree) 
                    
                    if remove_F_analysis
                        push!(subtree_normalized_rf_distances_noF, normalized_rf_dist_per_subtree_noF)
                    end
                    # The list stores all normalized RF distances for all subtrees 
                    # for one particular gene tree file 
                    
                end 

                mean_rf_dist_for_one_gene_tree = mean(subtree_normalized_rf_distances)
                push!(rf_distances, mean_rf_dist_for_one_gene_tree)
                push!(gene_tree_files_with_paths, (gene_file, gene_tree_path))
                
                if remove_F_analysis
                    mean_rf_dist_for_one_gene_tree_noF = mean(subtree_normalized_rf_distances_noF)
                    push!(rf_distances_noF, mean_rf_dist_for_one_gene_tree_noF)
                end

            end

            if n_inds == 1 && !isempty(gene_tree_content)
                gene_tree = readnewick(gene_tree_content)  

                #= when one individual per species, 
                There is no need to getsisters and then loop through all possible sub-trees: 
                just map the individuals to species directly =# 
                mapped_gene_tree = map_gene_tree_based_on_mapping_file(gene_tree, 
                                                        individual_to_species) 

                if remove_F_analysis # need to be run after mapping 
                    mapped_gene_tree_noF = deepcopy(mapped_gene_tree)
                    # Check if F exists before trying to delete
                    if "F" in tipLabels(mapped_gene_tree_noF)
                        deleteleaf!(mapped_gene_tree_noF, "F")
                    end
                end 

                # Calculate RF distance 
                pruned_species_tree, pruned_mapped_gene_tree = 
                    prune_nonoverlapping_tips(species_tree, mapped_gene_tree) 
                rf_dist = mudistance_semidirected(pruned_species_tree, 
                                                pruned_mapped_gene_tree, 
                                                preorder=true) 
                if remove_F_analysis
                    pruned_species_tree_noF, pruned_mapped_gene_tree_noF = 
                        prune_nonoverlapping_tips(species_tree_noF, mapped_gene_tree_noF) 
                    rf_dist_noF = mudistance_semidirected(pruned_species_tree_noF, 
                                                    pruned_mapped_gene_tree_noF, 
                                                    preorder=true) 
                end
                

                # Normalized RF distance -> RF distance / (2 * (n - 3)) 
                n = length(tipLabels(pruned_species_tree)) 
                normalized_rf_dist = rf_dist / (2 * (n - 3)) 
                
                if remove_F_analysis
                    n_noF = length(tipLabels(pruned_species_tree_noF))
                    normalized_rf_dist_noF = rf_dist_noF / (2 * (n_noF - 3))
                end 

                push!(rf_distances, normalized_rf_dist) 
                push!(gene_tree_files_with_paths, (gene_file, gene_tree_path)) 

                if remove_F_analysis 
                    push!(rf_distances_noF, normalized_rf_dist_noF)
                end 
            end 

        end
        
        if remove_F_analysis
            return gene_tree_files_with_paths, rf_distances, rf_distances_noF
        else
            return gene_tree_files_with_paths, rf_distances
        end
    end
end  # end @everywhere




#------------------------------------------------# 
# Counts the number of taxa in the gene trees 
# Plus calculate average internal branch lengths for locus trees 
#-------------------------------------------------#
#=
Make a csv recording every single locus tree file and its path 
Record the number of taxa in each gene/locus tree. The below function does: 
    1. count the number of taxa in each gene tree and locus tree 
        -> record in a csv file 
        -> error out if the number of taxa in gene tree and locus tree do not match 
    2. calculate the average internal branch length for each locus tree 
        -> record in the same csv file later
=#
@everywhere begin 
    """
    gene_and_locus_tree_taxa_and_internal_branch_stats:
    For a given replicate, this function compares gene trees and their
    corresponding locus trees to count the number of taxa and calculate
    average internal branch lengths.
    return: 
        1. genetree_filepath_list_per_rep: list of gene tree file paths 
        2. num_taxa_gene_tree_list_per_rep: list of number of taxa in each gene tree
            -> the same as num_taxa_locus_tree_list_per_rep
        3. avg_internal_branch_length_locus_tree_list_per_rep: 
            list of average internal branch lengths for each locus tree
        4. avg_internal_branch_length_gene_tree_list_per_rep: 
            list of average internal branch lengths for each gene tree
    Those returns will then be concated together for all reps in main() 
    """
    function gene_and_locus_tree_taxa_and_internal_branch_stats(
                        simulation_rep::Int,
                        folder_path_list::Vector,
                        rootfolder::String, 
                        n_inds::Int
                        )
        
        # initialization: 
        genetree_filepath_list_per_rep = String[]
        num_taxa_gene_tree_list_per_rep = Int[]
        avg_internal_branch_length_gene_tree_list_per_rep = Float64[]
        avg_internal_branch_length_locus_tree_list_per_rep = Float64[] 

        output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")

        # To identify the locus tree corresponding to each gene tree
        simphy_raw_gene_tree_folder = setup_rep_output_folders(folder_path_list,
                                                            simulation_rep, 
                                                            "genetrees_simphy")
        # To identify the single-copy gene trees 
        simphy_singlecopy_gene_tree_folder = setup_rep_output_folders(folder_path_list, 
                                                            simulation_rep, 
                                                            "genetrees_singlecopy")  
        genetrees_list = readdir(simphy_singlecopy_gene_tree_folder) # stores all gene tree files
        
        # Based on gene tree file name, find the corresponding locus tree file 
        for genetree_filename in genetrees_list
            # find gene tree path and update the list 
            gene_tree_filepath = joinpath(simphy_singlecopy_gene_tree_folder, genetree_filename)
            push!(genetree_filepath_list_per_rep, gene_tree_filepath)

            # Extract int_id from the gene tree file name
            int_part = replace(split(split(genetree_filename, "_")[5], ".")[1], "Int" => "") # "Int06" -> "06" 
            gene_part = replace(split(split(genetree_filename, "_")[4], ".")[1], "Gene" => "") # "Gene06" -> "06" 

            int_id = parse(Int, int_part) # Change this to integer "06" -> 6 
            gene_id = parse(Int, gene_part) # Change this to integer "06" -> 6 
            # parse the n_inds part: 

            temp_dir = joinpath(simphy_raw_gene_tree_folder, "Int$(int_id)", "1") 
            if isdir(temp_dir)
                locus_tree_file = joinpath(temp_dir, "l_trees.trees")
            else
                locus_tree_file = joinpath(simphy_raw_gene_tree_folder,
                                    "Int$(int_id)", "l_trees.trees")
            end 

            gene_tree = readnewick(gene_tree_filepath) # read the gene tree into network 
            locus_tree = readnewick(readlines(locus_tree_file)[gene_id]) # read the locus tree into network 

            to_remove = [t for t in tiplabels(locus_tree) if startswith(t, "Lost")] 
            for lost_tip in to_remove # only locus tree needs to remove lost tips 
                deleteleaf!(locus_tree, lost_tip)
            end

            # count the number of taxa in the gene tree
            num_taxa_gene_tree = length(tipLabels(gene_tree))  
            num_taxa_locus_tree = length(tipLabels(locus_tree))
            # The below check is a double-safty to make sure the pipeline is correct
            if num_taxa_gene_tree != num_taxa_locus_tree * n_inds # double check 
                error("Number of taxa differ between gene tree for gene $gene_tree_filepath")
            end 
            push!(num_taxa_gene_tree_list_per_rep, num_taxa_gene_tree)  

            # caculate average internal branch length for this specific locus tree 
            internal_branch_lengths_per_locus_tree = Float64[] 
            # calculate internal branch lengths for this specific gene tree 
            internal_branch_lengths_per_gene_tree = Float64[] 
            # remove external branch starts with "Lost-" first
            removedegree2nodes!(locus_tree) # remove degree 2 nodes first  
            removedegree2nodes!(gene_tree) # remove degree 2 nodes first 

            for edge in locus_tree.edge # internal branch lengths for locus tree 
                if !isexternal(edge) 
                    push!(internal_branch_lengths_per_locus_tree, edge.length)
                end     
            end
            for edge in gene_tree.edge # internal branch lengths for gene tree 
                if !isexternal(edge) 
                    push!(internal_branch_lengths_per_gene_tree, edge.length)
                end
            end

            # calculate average internal branch length
            avg_internal_branch_length_locus_tree = mean(internal_branch_lengths_per_locus_tree)
            avg_internal_branch_length_gene_tree = mean(internal_branch_lengths_per_gene_tree)
            push!(avg_internal_branch_length_locus_tree_list_per_rep, 
                  avg_internal_branch_length_locus_tree)
            push!(avg_internal_branch_length_gene_tree_list_per_rep, 
                  avg_internal_branch_length_gene_tree) 
        end

        return genetree_filepath_list_per_rep, 
               num_taxa_gene_tree_list_per_rep, 
               avg_internal_branch_length_locus_tree_list_per_rep,
               avg_internal_branch_length_gene_tree_list_per_rep 
    end
end         

#= Summarzation: 
The below function summarize the information from the above function 
and it adds more information regarding duplication/loss status 
and hidden paralogy classification 
It could be better integrated into the pipeline, but for now 
I think having a modularized function is better for future maintenance 

This function summarizes gene tree statistics into a CSV file.
It reads gene tree statistics from a CSV file and classifies each gene tree
based on duplication/loss status and hidden paralogy classification.
The function adds new columns to the statistics DataFrame to indicate:
- Whether the gene tree underwent duplication and loss.
- Whether the gene tree underwent loss only.
- Whether the gene tree underwent neither duplication nor loss.
- The hidden paralogy classification (false, weak, strong, or NaN).  
=#
"""
summarize_genetree_stats:
Summarizes gene tree statistics and classifies them based on duplication/loss status
and hidden paralogy classification.
"""
function summarize_genetree_stats(
        outfolder::String,
        paramname_root::String
    )
    #-- set up paths to read required csv files --# 
    df_internal_branch_length_stats = joinpath(outfolder, 
        "genetrees_taxa_and_internal_branch_$paramname_root.csv") 
    rf_true_genetrees_and_species_tree_path = joinpath(outfolder, 
        "RF_true_genetrees_and_species_tree_$paramname_root.csv")  
    rf_true_estimated_genetrees_stats_csv_path = joinpath(outfolder, 
        "RF_true_and_estimated_genetrees_$paramname_root.csv")

    if !isfile(df_internal_branch_length_stats) 
        error("Required CSV files not found in $df_internal_branch_length_stats")
    elseif !isfile(rf_true_genetrees_and_species_tree_path)
        error("Required CSV file not found: $rf_true_genetrees_and_species_tree_path") 
    elseif !isfile(rf_true_estimated_genetrees_stats_csv_path) 
        error("Required CSV file not found: $rf_true_estimated_genetrees_stats_csv_path")
    else
        df_stats = CSV.read(df_internal_branch_length_stats, DataFrame)
        df_rf_true_genetrees_and_species_tree = CSV.read(rf_true_genetrees_and_species_tree_path, DataFrame)
        df_rf_true_and_estimated_genetrees = CSV.read(rf_true_estimated_genetrees_stats_csv_path, DataFrame)  
    end  

    #-------------------------------#
    # Add information about the RF distance between 
    # true gene trees vs species tree 
    # genetrees vs species tree 
    #-------------------------------# 
    # change the columns name 
    rename!(df_rf_true_genetrees_and_species_tree, 
            "RF_score" => "RF_true_genetrees_vs_species_tree") 
    # Rename noF column if it exists
    if "RF_score_noF" in names(df_rf_true_genetrees_and_species_tree)
        rename!(df_rf_true_genetrees_and_species_tree, 
                "RF_score_noF" => "RF_true_genetrees_vs_species_tree_noF")
    end
    rename!(df_rf_true_and_estimated_genetrees, 
            "RF_score" => "RF_true_vs_estimated_genetrees") 
    # Select columns to keep - include noF if it exists
    cols_to_keep = ["gene_tree_file_path", "RF_true_genetrees_vs_species_tree"]
    if "RF_true_genetrees_vs_species_tree_noF" in names(df_rf_true_genetrees_and_species_tree)
        push!(cols_to_keep, "RF_true_genetrees_vs_species_tree_noF")
    end
    select!(df_rf_true_genetrees_and_species_tree, cols_to_keep)
    select!(df_rf_true_and_estimated_genetrees, 
            ["gene_tree_file_path", "RF_true_vs_estimated_genetrees"]) 
    # merge the two dataframes based on gene_tree_file_path
    df_stats = leftjoin(df_stats, df_rf_true_genetrees_and_species_tree, 
                        on = "gene_tree_file_path") 
    df_stats = leftjoin(df_stats, df_rf_true_and_estimated_genetrees, 
                        on = "gene_tree_file_path") 

    #-------------------------------#
    # For the final stats csv file 
    # get the gene duplication and loss classification 
    # and the gene loss or nothing classification
    #-------------------------------# 
    # get gene duplication and loss classification 
    gene_duplication_loss_csv_path = joinpath(outfolder, 
        "Simphy_gene_duplication_and_loss_$paramname_root.csv")  
    # read csv for gene duplication and loss classification 
    if !isfile(gene_duplication_loss_csv_path)
        error("Gene duplication and loss CSV file not found: $gene_duplication_loss_csv_path")
    else 
        gene_duplication_loss_csv = CSV.read(gene_duplication_loss_csv_path, DataFrame)
        gene_duplication_and_loss_trees = Set{String}(gene_duplication_loss_csv[!, "gene_tree_file_path"]) 

        # find HP classifications: 
        # Keep only rows with false_HP == "Y" (if the column exists)
        filtered_false_HP = gene_duplication_loss_csv[gene_duplication_loss_csv[!, "false_HP"] .== "Y", :]
        false_HP_trees = Set{String}(filtered_false_HP[!, "gene_tree_file_path"])
        filtered_weak_HP = gene_duplication_loss_csv[gene_duplication_loss_csv[!, "weak_HP"] .== "Y", :]
        weak_HP_trees = Set{String}(filtered_weak_HP[!, "gene_tree_file_path"])
        filtered_strong_HP = gene_duplication_loss_csv[gene_duplication_loss_csv[!, "strong_HP"] .== "Y", :]
        strong_HP_trees = Set{String}(filtered_strong_HP[!, "gene_tree_file_path"])
    end 

    # read csv for gene loss only classification 
    gene_loss_only_csv_path = joinpath(outfolder, 
        "Simphy_gene_loss_only_$paramname_root.csv") 
    if !isfile(gene_loss_only_csv_path)
        error("Gene loss only CSV file not found: $gene_loss_only_csv_path")
    else
        gene_loss_only_csv = CSV.read(gene_loss_only_csv_path, DataFrame)
        gene_loss_only_trees = Set{String}(gene_loss_only_csv[!, "gene_tree_file_path"]) 
    end

    # Initialize new columns
    df_stats[!, "gene_duplication_and_loss_or_not"] = fill("N", nrow(df_stats))
    df_stats[!, "gene_loss_only_or_not"] = fill("N", nrow(df_stats))
    df_stats[!, "nothing_or_not"] = fill("N", nrow(df_stats))
    df_stats[!, "false_HP"] = fill("NaN", nrow(df_stats))
    df_stats[!, "weak_HP"] = fill("NaN", nrow(df_stats))
    df_stats[!, "strong_HP"] = fill("NaN", nrow(df_stats))

    # Populate new columns
    for i in 1:nrow(df_stats)
        gene_path = df_stats[i, "gene_tree_file_path"]
        gene_path_str = String(gene_path)  # Convert to String to avoid SubString issues

        # Duplication and loss classification
        if gene_path_str in gene_duplication_and_loss_trees
            df_stats[i, "gene_duplication_and_loss_or_not"] = "Y"
            # only if the gene experiences both duplication and loss  
            # then we first assign N to HP types otherwise should be NaN
            df_stats[i, "false_HP"] = "N" 
            df_stats[i, "weak_HP"] = "N"
            df_stats[i, "strong_HP"] = "N" 
        elseif gene_path_str in gene_loss_only_trees
            df_stats[i, "gene_loss_only_or_not"] = "Y"
        else
            df_stats[i, "nothing_or_not"] = "Y"
        end

        # Hidden paralogy classification
        if gene_path_str in false_HP_trees
            df_stats[i, "false_HP"] = "Y"
        elseif gene_path_str in weak_HP_trees
            df_stats[i, "weak_HP"] = "Y"
        elseif gene_path_str in strong_HP_trees
            df_stats[i, "strong_HP"] = "Y"
        else
            continue # remains "NaN" 
        end 
    end

    # write back to csv 
    stats_csv_path = joinpath(outfolder, 
        "genetrees_stats_$paramname_root.csv") 
    CSV.write(stats_csv_path, df_stats) 
end 




#-----------------------------------------------#
# Main postprocessing function  
#-----------------------------------------------#
function main()
    println("Starting postprocessing for parameter set: $paramname_root")
    #-----------------------------------------------#
    # Run ASTRAL-species tree comparisons 
    # Compare estimated species tree from ASTRAL vs true species tree 
    #-----------------------------------------------# 
    println("="^30)
    println("Calculating RF distances between ASTRAL and true species trees...")
    results = pmap(simulation_rep -> compare_astral_simphy_truespecies_1rep(
                    simulation_rep, 
                    folder_path_list, 
                    rootfolder,
                    n_inds
                    ), 1:n_reps)
    RF_astral_true_score = [result[1] for result in results]
    RF_astral_simphy_gene_files = [result[2] for result in results]
    RF_astral_simphy_score = [result[3] for result in results]
    RF_astral_simphy_score_noF = [result[4] for result in results]

    # Add RF score into the simphy csv file 
    csv_file = joinpath(outfolder, "simulation_$paramname_root.csv")
    if isfile(csv_file)
        df = CSV.read(csv_file, DataFrame)
        # Only update rows for the replicates we're processing
        if "RF_astral_true" ∉ names(df)
            df[!, "RF_astral_true"] = fill(NaN, nrow(df))
        end
        for (i, rep_id) in enumerate(1:n_reps)
            if rep_id <= nrow(df)
                df[rep_id, "RF_astral_true"] = RF_astral_true_score[i]
            end
        end
        CSV.write(csv_file, df)
        println("Updated Simphy CSV file with RF scores")
    else
        println("Warning: Simphy CSV file not found: $csv_file")
        println("Creating a basic CSV file with RF scores...")
        # Create a basic DataFrame with the RF scores
        df = DataFrame(
            replicate = 1:n_reps,
            RF_astral_true = RF_astral_true_score
        )
        CSV.write(csv_file, df)
        println("Created basic Simphy CSV file with RF scores")
    end

    # Save the RF score between gene trees from simphy and ground truth species tree 
    println("Saving RF scores between SimPhy gene trees and true species tree...")
    for simulation_rep in 1:n_reps
      rep_output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")  
      rf_simphy_true_file = joinpath(rep_output_folder, 
                                    "RF_btw_simphy_vs_true_species_tree_$simulation_rep.txt")
      rf_simphy_true_file_noF = joinpath(rep_output_folder, 
                                    "RF_btw_simphy_vs_true_species_tree_noF_$simulation_rep.txt")
      
      # Save regular RF scores
      open(rf_simphy_true_file, "w") do f
          gene_files = RF_astral_simphy_gene_files[simulation_rep]
          rf_scores = RF_astral_simphy_score[simulation_rep]
          
          for (i, (filename, filepath)) in enumerate(gene_files)
              if i <= length(rf_scores)
                  write(f, "$filepath $(rf_scores[i])\n")
              end
          end
      end
      
      # Save noF RF scores if available
      if !isempty(RF_astral_simphy_score_noF[simulation_rep])
          open(rf_simphy_true_file_noF, "w") do f
              gene_files = RF_astral_simphy_gene_files[simulation_rep]
              rf_scores_noF = RF_astral_simphy_score_noF[simulation_rep]
              
              for (i, (filename, filepath)) in enumerate(gene_files)
                  if i <= length(rf_scores_noF)
                      write(f, "$filepath $(rf_scores_noF[i])\n")
                  end
              end
          end
      end
    end 

    #-----------------------------------------------# 
    # Compare gene trees from SimPhy and IQ-TREE
    # Calculate RF distance between
    # gene trees from simphy vs gene trees estimated by IQ-TREE 
    #-----------------------------------------------# 
    println("="^30)
    println("Calculating RF distances between SimPhy and IQ-TREE gene trees...")
    gene_tree_rf_results = pmap(simulation_rep -> compare_gene_trees_iqtree_simphy_1rep(
                          simulation_rep, 
                          folder_path_list, 
                          rootfolder
                          ), 1:n_reps)
    rf_csv_file = joinpath(outfolder, 
              "RF_true_and_estimated_genetrees_$paramname_root.csv")
    
    # Prepare data for CSV - flatten results from all replicates
    all_gene_paths = vcat([result[2] for result in gene_tree_rf_results]...)
    all_rf_scores = vcat([result[1] for result in gene_tree_rf_results]...)
    df_rf = DataFrame(
        gene_tree_file_path = all_gene_paths,
        RF_score = all_rf_scores
    )
    CSV.write(rf_csv_file, df_rf)
    println("Created RF distance CSV: $rf_csv_file")

    #-----------------------------------------------#
    # Read gene duplication and loss classification files
    # make HP classification based on locus trees 
    #-----------------------------------------------#
    println("="^30)
    println("Making HP classification based on locus trees...")

    # Read the SimPhy classification files
    gene_duplication_and_loss_file = joinpath(outfolder, "Simphy_gene_duplication_and_loss_$paramname_root.csv")
    gene_loss_only_file = joinpath(outfolder, "Simphy_gene_loss_only_$paramname_root.csv")
    
    # Read the gene tree lists (if files exist and are not empty)
    gene_duplication_and_loss_trees = Set{String}()
    gene_loss_only_trees = Set{String}()

    if isfile(gene_duplication_and_loss_file) 
        # For dup != 0 and loss != 0 the file should not be empty 
        gene_duplication_and_loss_df = CSV.read(gene_duplication_and_loss_file, DataFrame)
        data_lines = gene_duplication_and_loss_df[!, "gene_tree_file_path"]  # Extract just the file paths
        # ensure it's a Set of Strings -> Thus, it can be handled for empty when dup_rate == 0 
        gene_duplication_and_loss_trees = Set{String}(filter(!isempty, data_lines)) 
        # If the file is empty, this will be an empty set -> gene_duplication_and_loss_trees = Set{String}() 
    else 
        error("Warning: Gene trees with duplication and loss file not found: $gene_duplication_and_loss_file")
    end 
    
    if isfile(gene_loss_only_file)  
        gene_loss_only_df = CSV.read(gene_loss_only_file, DataFrame)
        data_lines = gene_loss_only_df[!, "gene_tree_file_path"]  # Extract just the file paths
        gene_loss_only_trees  = Set{String}(filter(!isempty, data_lines)) 
        # if the file is empty, this will be an empty set -> gene_loss_only_trees = Set{String}()
    else 
        error("Warning: Gene loss only file not found: $gene_loss_only_file")
    end 

    #-----------------------------------------------# 
    # Classify gene trees with hidden paralogy
    # Add hidden paralogy classification columns to the duplication and loss CSV
    false_HP_paths, weak_HP_paths, strong_HP_paths = postexamine_HP_based_on_locus_tree(
        outfolder, n_reps, paramname_root, folder_path_list)
    
    # Read the existing CSV file and add the new columns
    if isfile(gene_duplication_and_loss_file)
        # Read existing CSV
        df_existing = CSV.read(gene_duplication_and_loss_file, DataFrame)
        
        # Initialize new columns with "N"
        df_existing[!, "false_HP"] = fill("N", nrow(df_existing))
        df_existing[!, "weak_HP"] = fill("N", nrow(df_existing))
        df_existing[!, "strong_HP"] = fill("N", nrow(df_existing))
        
        # Set "Y" for paths that match each category
        for i in 1:nrow(df_existing)
            gene_path = df_existing[i, "gene_tree_file_path"]
            if gene_path in false_HP_paths
                df_existing[i, "false_HP"] = "Y"
            elseif gene_path in weak_HP_paths
                df_existing[i, "weak_HP"] = "Y"
            elseif gene_path in strong_HP_paths
                df_existing[i, "strong_HP"] = "Y"
            else
                println("  -> No match found")
            end
        end
        # Write back to the same CSV file
        CSV.write(gene_duplication_and_loss_file, df_existing)
        println("Updated gene duplication and loss CSV with hidden paralogy classification")
    else
        println("Warning: Could not update CSV file - file not found: $gene_duplication_and_loss_file")
    end

    # -----------------------------------------------# 
    # Calculate RF distance between 
    #= SimPhy gene trees vs true species tree 
    For n_inds = 1 and n_inds > 1 cases
    Those two cases are different =# 
    # -----------------------------------------------# 
    println("="^30)
    println("Calculating RF dist between SimPhy gene trees and true species tree...") 
    println("Directly compare trees for n_inds == $n_inds.") 
    
    # Get RF data from the function
    simphy_species_csv_data = create_simphy_species_rf_summary(
      gene_duplication_and_loss_trees,
      gene_loss_only_trees
    )
    # pointer to point back 
    # Write SimPhy vs species tree RF summary to CSV (gene_tree_file_path, RF_score, and RF_score_noF)
    rf_simphy_species_csv = joinpath(outfolder, 
                            "RF_true_genetrees_and_species_tree_$paramname_root.csv")
    if !isempty(simphy_species_csv_data)
        df_simphy_species = DataFrame(
            gene_tree_file_path = [row[1] for row in simphy_species_csv_data],
            RF_score = [row[2] for row in simphy_species_csv_data],
            RF_score_noF = [row[3] for row in simphy_species_csv_data]
        )
        CSV.write(rf_simphy_species_csv, df_simphy_species)
        println("Created RF distance CSV: $rf_simphy_species_csv")
    end
    
    #-----------------------------------------------#
    # Add HP counts to simulation CSV file
    #-----------------------------------------------# 
    # Count HP categories by replicate using the RF CSV files
    println("="^30)
    println("Counting hidden paralogy categories by replicate...")
    if isfile(csv_file)
        df_sim = CSV.read(csv_file, DataFrame)
        
        # Initialize HP count columns if they don't exist
        if "num_false_HP" ∉ names(df_sim)
            df_sim[!, "num_false_HP"] = fill(0, nrow(df_sim))
        end
        if "num_weak_HP" ∉ names(df_sim)
            df_sim[!, "num_weak_HP"] = fill(0, nrow(df_sim))
        end
        if "num_strong_HP" ∉ names(df_sim)
            df_sim[!, "num_strong_HP"] = fill(0, nrow(df_sim))
        end
        
        # Count HP categories from the RF CSV files for each replicate
        gene_duplication_loss_csv_file = joinpath(outfolder, "Simphy_gene_duplication_and_loss_$paramname_root.csv")
        df_rf = CSV.read(gene_duplication_loss_csv_file, DataFrame)
        for rep_id in 1:n_reps
            if rep_id <= nrow(df_sim)
                # Filter RF data for this replicate
                rep_string = pad_number(rep_id, n_reps)
                # Build a boolean mask robust to Missing values in the path column
                rep_mask = [(!ismissing(path) && occursin("/rep$rep_string/", path)) for path in df_rf.gene_tree_file_path]
                rep_data = df_rf[rep_mask, :]

                # Count HP categories for this replicate, guarding against Missing
                num_false = count(x -> (!ismissing(x) && x == "Y"), rep_data.false_HP)
                num_weak = count(x -> (!ismissing(x) && x == "Y"), rep_data.weak_HP)
                num_strong = count(x -> (!ismissing(x) && x == "Y"), rep_data.strong_HP)
                nrows = length(rep_data.gene_tree_file_path)
                hp_sum = num_false + num_weak + num_strong
                if hp_sum != nrows # validation check 
                    error("HP validation failed")
                end
                
                # Update simulation CSV
                df_sim[rep_id, "num_false_HP"] = num_false
                df_sim[rep_id, "num_weak_HP"] = num_weak
                df_sim[rep_id, "num_strong_HP"] = num_strong
            end
        end
        
        # Write updated simulation CSV
        CSV.write(csv_file, df_sim)
        println("Updated simulation CSV with HP counts")
    end

    #-----------------------------------------------#
    # Gene and locus tree taxa and internal branch stats
    #-----------------------------------------------#
    println("="^30)
    println("Calculating gene and locus tree taxa and internal branch stats...")
    all_genetree_filepaths = String[]
    all_num_taxa_gene_trees = Int[]
    all_avg_internal_branch_length_locus_trees = Float64[]
    all_avg_internal_branch_length_gene_trees = Float64[]
    for simulation_rep in 1:n_reps
        gene_tree_filepaths, 
        num_taxa_gene_trees, 
        avg_internal_branch_length_locus_trees,
        avg_internal_branch_length_gene_trees = gene_and_locus_tree_taxa_and_internal_branch_stats(
            simulation_rep,
            folder_path_list,
            rootfolder,
            n_inds
        )
        append!(all_genetree_filepaths, gene_tree_filepaths)
        append!(all_num_taxa_gene_trees, num_taxa_gene_trees)
        append!(all_avg_internal_branch_length_locus_trees, avg_internal_branch_length_locus_trees)
        append!(all_avg_internal_branch_length_gene_trees, avg_internal_branch_length_gene_trees)
    end
    # Write to CSV
    taxa_internal_branch_csv = joinpath(outfolder, 
        "genetrees_taxa_and_internal_branch_$paramname_root.csv")
    df_taxa_internal = DataFrame(
        gene_tree_file_path = all_genetree_filepaths,
        num_taxa_gene_tree = all_num_taxa_gene_trees,
        avg_internal_branch_length_locus_tree = all_avg_internal_branch_length_locus_trees,
        avg_internal_branch_length_gene_tree = all_avg_internal_branch_length_gene_trees
    )
    # save to the genetrees_taxa_and_internal_branch_*.csv file  
    CSV.write(taxa_internal_branch_csv, df_taxa_internal) 

    
    println("="^30)
    println("Summarizing gene tree statistics with classifications...") 

    summarize_genetree_stats(outfolder, paramname_root) 

    println("Information about gene tree category is added to")
    println("=> $taxa_internal_branch_csv\n")  

    # -----------------------------------------------#
    # major summarization 
    # -----------------------------------------------# 
    # summarize rep-level information 
    summarize_rf_distances(
        joinpath(outfolder, "genetrees_stats_$paramname_root.csv"), 
        joinpath(outfolder, "simulation_$paramname_root.csv"),
        paramname_root
    )

    # Extract parameter set level information 
    output_csv_file_path = joinpath(outfolder, "summary_$paramname_root.csv")
    extract_summarize_simphy_summary(paramname_root, 
                                    outfolder, 
                                    output_csv_file_path)
    
    # Add parameter-level RF statistics to the same file
    summarize_parameter_level_RF_statistics(paramname_root, 
                                        outfolder, 
                                        output_csv_file_path)

    #-----------------------------------------------#
    # Finished postprocessing --> oh yay! 
    #-----------------------------------------------# 
    println("="^30)
    println("Postprocessing completed successfully for parameter set: $paramname_root")

end


# Run the main function
main()
