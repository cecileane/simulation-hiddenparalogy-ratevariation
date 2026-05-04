#=
================================================================================
  simulation.jl — Full Phylogenomic Simulation Pipeline
================================================================================

OVERVIEW
--------
This script runs an end-to-end phylogenomic simulation pipeline across multiple
replicates in parallel. Starting from a fixed species tree, it produces estimated
gene trees and a coalescent-based species tree for each parameter configuration.

PIPELINE STAGES
---------------
1. GENE TREE SIMULATION (SimPhy)
   - Simulates gene trees along the species tree under the multispecies coalescent,
     optionally with gene duplication and loss rates (--dup_rate, --loss_rate).
   - Supports rate variation across genes (G), lineages (L), or both (GL/G*L).
   - Multiple individuals per taxon are supported (--n_inds).

2. HIDDEN PARALOGY FILTERING (newick string modification)
   - After SimPhy, gene trees are post-processed to mimic hidden paralogy:
       * Trees with >1 gene copy per individual (detectable paralogy) → discarded
       * Trees with ≤3 taxa after gene loss (uninformative) → discarded
       * Trees with exactly 1 copy per individual and ≥4 taxa → retained;
         gene copy IDs are stripped from taxon names
   - If too few trees pass filtering, SimPhy is re-run (up to --max_iteration_simphy
     times) until the target count is reached or the minimum threshold (--min_gene_proportion
     × --n_genes) is met.

3. SEQUENCE SIMULATION (Seq-Gen)
   - Simulates DNA alignments along each filtered gene tree using the HKY model.
   - Substitution parameters (κ, base frequencies, Γ-shape α) are either fixed
     across all genes (all_genes_same_HKY) or drawn independently per gene from
     empirical distributions (all_genes_diff_HKY; default).

4. ALIGNMENT CONCATENATION (Python)
   - Individual per-gene NEXUS alignments are converted and concatenated into a
     single FASTA file per replicate (via concatenate_seq.py).

5. GENE TREE ESTIMATION (IQ-Tree)
   - IQ-Tree infers a maximum-likelihood gene tree for each simulated alignment.
   - Trees are collected into a single besttrees.tre file per replicate.

6. SPECIES TREE INFERENCE (ASTRAL)
   - A taxon mapping file is generated automatically from the gene trees.

REPRODUCIBILITY
---------------
All random seeds are deterministically derived from a single master seed, which
is itself generated from the parameter configuration name. This ensures that
results for any parameter setting are fully reproducible.

OUTPUT STRUCTURE
----------------
output/
└── <paramname_root>/          e.g. DUP0.0003-LOS0.0003-RVG-N_ind1-SF1.0-genelen1000
    ├── rep001/
    │   ├── genetrees_simphy/     Raw SimPhy output (kept for inspection)
    │   ├── genetrees_singlecopy/ Filtered & renamed gene trees (.trees)
    │   ├── seqgenfolder/         Simulated alignments (.fasta)
    │   ├── iqtreefolder/         IQ-Tree output (besttrees.tre, …)
    │   └── astralfolder/         ASTRAL output (astral.tre, mapping file)
    ├── rep002/ …
    ├── simulation_<paramname>.csv     Per-replicate SimPhy run statistics
    ├── Simphy_gene_duplication_and_loss_<paramname>.csv
    ├── Simphy_gene_loss_only_<paramname>.csv
    ├── random_seed_simphy/seqgen/iqtree.txt
    ├── arguments-<paramname>.log      All arguments + runtime summary
    └── screen_<paramname>.log         Full console output (if --log true)

USAGE
-----
Run from the repository root, e.g.:

  julia -p <nprocs> scripts/simulation.jl \
      --dup_rate 0.0003 --loss_rate 0.0003 \
      --ratevar G --n_reps 100 --n_genes 500 \
      --n_inds 1 --SF 0.5 --gene_len 1000

  On the darwin cluster:  use /nobackup  and run inside tmux
  On the franklin cluster: use /nobackup2 and run inside tmux

REQUIRED ARGUMENTS
------------------
  --dup_rate   FLOAT   Gene duplication rate (0 = no duplication)
  --loss_rate  FLOAT   Gene loss rate (0 = no loss); must be ≤ dup_rate
  --ratevar    STR     Rate variation: N | G | L | GL | G*L
  --n_reps     INT     Number of independent replicates
  --n_genes    INT     Target number of gene trees per replicate

OPTIONAL ARGUMENTS
------------------
  --n_inds              INT    Individuals per taxon (default: 1)
  --Ne                  INT    Effective population size (default: 1000)
  --SF                  FLOAT  Ne scaling factor; SF<1 → more ILS (default: 1.0)
  --gene_len            INT    Alignment length in bp (default: 1000)
  --seqgen_model        STR    all_genes_same_HKY | all_genes_diff_HKY (default)
  --max_iteration_simphy INT   Max SimPhy re-runs per replicate (default: 20000)
  --min_gene_proportion FLOAT  Min fraction of target genes required (default: 0.8)
  --max_taxa_missing    INT    Max taxa allowed missing after gene loss (default: 0)
  --min_locus_tree_tips INT    Min tips required in a locus tree (default: 8)
  --debug_mode          BOOL   Keep all intermediate folders (default: false)
  --log                 BOOL   Save console output to log file (default: true)
=#

using Distributed  
using ArgParse
using TimerOutputs
using Dates
using TimeZones 
using Statistics

@everywhere using DataFrames
@everywhere using CSV
@everywhere using PhyloNetworks
@everywhere using Printf 
@everywhere include("utilities.jl")

const to = TimerOutput() 
tz = TimeZone("America/Chicago") 
current_time_tz = ZonedDateTime(now(), tz) 
time = Dates.format(current_time_tz, "yyyy-mm-dd HH:MM:SS zzz") 

function parse_commandline()
  s = ArgParseSettings() 
  @add_arg_table s begin 
    "--dup_rate"
      help = """Gene duplication and loss rate: specify the gene duplication rate.
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
    "--n_genes"
      help = "Number of genes" 
      arg_type = Int
      required = true
    "--n_inds" 
      help = "Number of individuals/accessions per species (Default = 1)"  
      arg_type = Int
      default = 1
    "--Ne"
      help = "Effective population size (Default = 1000)"
      arg_type = Int
      default = 1000
    "--generation_time"
      help = "Time unit per generation (Default = 1)"
      arg_type = Int
      default = 1
    "--max_iteration_simphy"
      help = "Maximum iteration to re-run simphy to get enough gene trees"
      arg_type = Int
      default = 20000 
    "--min_gene_porportion"
      help = "Num of simulated gene trees in a rep >= min_gene_proportion * n_genes"
      arg_type = Float64
      default = 0.8 
    "--SF"
      help = "Scaling factor to scale effective population Ne
        (Default = 1.0, no scaling)"
      arg_type = Float64
      default = 1.0
    "--max_taxa_missing"
      help = "Maximum number of taxa allowed to be missing in a gene tree after gene loss"
      arg_type = Int
      default = 0
    "--min_locus_tree_tips"
      help = "Minimum number of tips required in locus tree after gene dup/loss"
      arg_type = Int
      default = 8
    "--seqgen_model"
      help = "Substitution model for seq-gen (Default = all_genes_diff_HKY)" * 
              "Options: all_genes_same_HKY, all_genes_diff_HKY" 
      arg_type = String
      default = "all_genes_diff_HKY"  
    "--gene_len" 
      help = "The length of simulated gene sequences (Default = 1000 bp)" 
      arg_type = Int 
      default = 1000                

    # logging and debugging 
    "--debug_mode" 
      help = "Debugging mode: Keep all temporary folders"
      default = false
      arg_type = Bool
    "--log"
      help = "Enable logging: Save all screen output to screen_<paramname>.log"
      default = true
      arg_type = Bool

  end 
  return parse_args(s) 
end

#-----------------------------------------------#       
#         Initialization
#-----------------------------------------------#
#--------------- Parse arguments ---------------# 
parsed_args = parse_commandline()
# A number indicates rate of gene duplication rates (0 = no dup):
dup_rate = parsed_args["dup_rate"] 
# A number indicates rate of gene loss rates (0 = no loss):
loss_rate = parsed_args["loss_rate"]
# "N", "G", "L" or "GL" or "G*L" to include genexlineage rate variation: 
ratevar = parsed_args["ratevar"] 
n_reps = parsed_args["n_reps"] # number of replicates 
n_genes = parsed_args["n_genes"] # number of genes 
n_inds = parsed_args["n_inds"] # Number of individuals per taxa -- default = 1
Ne = parsed_args["Ne"] # Effective population size, default = 1000 
generation_time = parsed_args["generation_time"] # Generation time, default = 1 
# Maximum number of iteration of re-running simphy: 
max_iteration_simphy = parsed_args["max_iteration_simphy"] 
# Percentage * n_genes =  n_genes_min, the mimum number of gene trees in each rep: 
min_gene_porportion = parsed_args["min_gene_porportion"] 
# After reaching max_iteration_simphy, -only if num of genes > n_genes_min --> proceed: 
n_genes_min = n_genes * min_gene_porportion  
# Scaling factor to scale the effective population Ne
SF = parsed_args["SF"] # SF stands for scaling factor 
SF = Float64(SF) 
# debug mode - if true, keep all temporary folders, else, remove them to save space 
debug_mode = parsed_args["debug_mode"]
# log mode - if true, save all screen output to a log file
log_output = parsed_args["log"] 
max_taxa_missing = parsed_args["max_taxa_missing"] 
seqgen_model = parsed_args["seqgen_model"]  # substitution model for seq-gen 
min_locus_tree_tips = parsed_args["min_locus_tree_tips"]   
gene_len = parsed_args["gene_len"]  # length of simulated gene sequences  

#--------------- set up folders and Path ---------------# 
# set all configuration parameters here 
rootfolder = pwd()
paramname_root = set_up_paramname_root(dup_rate, loss_rate, ratevar, n_inds, SF, gene_len) 
outfolder = joinpath(rootfolder, "output", paramname_root) 

# The below function checks if a path exists:
# If yes, the function will ask user's input tp decide if to remove the path 
check_existing_dir([outfolder]) # see utilities.jl --> input needs to be a list
mkpath(outfolder) # create a new one if there is no pre-existing outfolder

folder_path_list = [] # create folder path 
for n in 1:n_reps # All outputs goes to each rep file 
  rep_number_string = pad_number(n, n_reps)
  rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
  push!(folder_path_list, rep_folder_path)
end 

#--------------- set up seeds ---------------# 
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
master_seed = generate_master_seed(params_dict_for_seed_setting) 

software_names = ["simphy", "seqgen", "iqtree"] 
# set up seeds for all softwares 
seed_dic = generate_software_seeds(master_seed, software_names) # see utility.jl 
seed_simphy = seed_dic["simphy"] 
seed_seqgen = seed_dic["seqgen"]
seed_iqtree = seed_dic["iqtree"]
# astral is a deterministic program with no need to specify a seed 

#--------------- Set up and check general parameters ---------------# 
# Check if ratevar belongs to the following: 
valid_ratevars = ["N", "G", "L", "GL", "GxL", "G*L"] 
if !(ratevar in valid_ratevars) 
  error("Invalid value for --ratevar: $ratevar. Valid options are: $(join(valid_ratevars, ", "))")
end

# SimPhy sets dup_rate > loss_rate, or otherwise there will be error message 
if dup_rate < loss_rate
  error("Invalid loss rate $loss_rate -- have it lower than dup_rate $dup_rate")
end

#= Second, find the species tree with lineage variation, which is used for L, GL and G*L. 
  See speciestree.jl: 
    1) each branch in tree_CU has a specific length in CU 
    2) each branch in tree_sub has a specific length in substitution rate per site 
  To calculate substitution rate per site per generation: 
  substitution rate per site per generation = substitution rate per site / CU / 2 * Ne for each specific branch 
=# 

#-----------------------------------------------#       
#         set up SimPhy parameters
#-----------------------------------------------#
# Modify the master simphy config file based on arguments and then save it to new config file
master_conf = joinpath(rootfolder,"simphy-configs/", "simphysim-conf-master")
conf_content = read(master_conf, String) # read the master config file into a string 

#Set up the parameters
parameters = ""

# set up effective population size 
scaled_Ne = Ne * SF 
parameters *= "-sp f:$scaled_Ne //Effective population size\n" 
parameters *= "-sg f:1 //Generation time\n" 
parameters *= "-ll $min_locus_tree_tips // Number of minimum locus tree tips\n"

# To simulate different gene duplication and loss rates
if dup_rate != 0 # if dup_rate is 0 then no -lb parameter 
  parameters *= "-lb f:$dup_rate // gene duplication rate\n" 
end
if loss_rate != 0 # if loss_rate is 0 then no -ld paramater
  parameters *= "-ld f:$loss_rate // gene loss rate\n" 
end

# To simulate substitution rate variation
if occursin("G", ratevar) # gene-family-speciic rate heterogenity : "G" or "GL" or "G*L"
  parameters *= "-hl ln:-0.19,0.6164414002968976 //log-normal distribution of gene rates\n"
  # hl is Gene-family-specific rate heterogeneity modifiers 
end

# To simulate variation across lineages (ratevar = "L" or "GL" or "GxL") 
if occursin("L", ratevar) # add tree with variation across lineages
  # For lineage-specific rate variation, we need to modify the species tree
  # see speciestree.jl for details
  # For instance, if scaling factor = 0.5, 
  # then effective population size is multiplied by 0.5 
  species_tree = "(A:3440*1.7784334944675035,((((B:880*0.18954057365099278,C:880*0.165743416346785):1710*0.6875429065797596,(D:930*1.1237159594044575,E:930*0.951588827019626):1660*0.679379951364618):170*0.05993108469820257,F:2760*0.9500336728638588):180*0.09042274190364852,(G:500*2.043344482162159,H:500*4.601350105819277):2440*2.380352926246597):500*1.778619857472514);"
  parameters *= "-s $species_tree\n" 
  parameters *= "-su f:0.000019526049565237014 //substitution rate\n" #see speciestree.jl
else 
  # set up the scaling factor for branch lengths: 
  species_tree = "(A:3440,((((B:880,C:880):1710,(D:930,E:930):1660):170,F:2760):180,(G:500,H:500):2440):500);"  
  parameters *= "-s $species_tree\n"
  # If ratevar doesn't contain L, then add tree without variations across lineages 
  parameters *= "-su f:0.000019526049565237014 //substitution rate\n" #see speciestree.jl
end

# To simulate multiple individuals per species 
if n_inds > 1 # if n_inds == 1, no "-si" argument & use the default 
  parameters *= "-si f:$n_inds // number of individuals per tree tip\n" 
end 

# Store all parameters into one string, which will be updated in the below loop
simphy_conf_content = parameters * conf_content # combine parameters with master config

print(simphy_conf_content) # print the content to check if everything is correct

#-----------------------------------------------#       
#  Push global params to all processors
#-----------------------------------------------#
@everywhere global folder_path_list = $folder_path_list
@everywhere global seed_simphy = $seed_simphy
@everywhere global seed_seqgen = $seed_seqgen
@everywhere global seed_iqtree = $seed_iqtree
@everywhere global max_iteration_simphy = $max_iteration_simphy
@everywhere global n_genes = $n_genes
@everywhere global simphy_conf_content = $simphy_conf_content
@everywhere global rootfolder = $rootfolder
@everywhere global n_genes_min = $n_genes_min
@everywhere global n_reps = $n_reps
@everywhere global Ne = $Ne 
@everywhere global paramname_root = $paramname_root
@everywhere global debug_mode = $debug_mode
@everywhere global dup_rate = $dup_rate 
@everywhere global n_inds = $n_inds 
@everywhere global SF = $SF
@everywhere global time = $time 
@everywhere global species_tree = $species_tree
@everywhere global log_output = $log_output
@everywhere global outfolder = $outfolder
@everywhere global max_taxa_missing = $max_taxa_missing
@everywhere global seqgen_model = $seqgen_model 
@everywhere global gene_len = $gene_len 

# Initialize worker-specific logging variable
@everywhere global worker_logfile = nothing

#-----------------------------------------------#       
#  run SimPhy + modify newick strings 
#-----------------------------------------------#
#= The below function has two goals: 
    1. simulate gene trees using Simphy (rep = 1 and loop through n_reps)
    2. modify newick strings to mimic hidden paralogy 
  If the gene tree has
  * more than one repeated gene copy: 
    Then, it is obvious that gene duplication happened, so paralogy is not hidden. 
    There will be *no* associated gene tree file passed to downstream analysis.
  * Some gene loss events and there are <= 3 taxa left: 
    Not phylogenetically informative. 
    There will be *no* associated gene tree file passed to downstream analysis.
  * 0 or 1 gene copy per individual, and >= 4 taxa left: 
    If there is paralogy, then it's hidden. 
    The associated gene tree file is the same as the original, except for taxon names.
    Gene copy id is removed, the taxon name still contains the species + individual number.
  * If dup_rate = 0, then all gene copies are orthologous: 
    There cannot be any hidden paralogy. 
    The below code will only change the tip names by removing gene copy id. 
  * This newick string modification might remove some trees from the pipeline. 
    a. If the resulting num of gene trees <= n_genes (target), 
      we re-run simphy with <= max_iteration_simphy times to generate enough gene trees. 
    b. If iteration > max_iteration_simphy, num of gene trees > n_genes_min 
      -> pipeline continues 
    c. If iteration > max_iteration_simphy, num of gene tress < n_genes_min 
      -> pipeline continues but will write a warning.log
    Info about num of gene trees generated by a rep is saved in simphy_simulation.csv
=#

@everywhere begin 
  """
  rerun_simphy_1rep_1int: 
    Run simphy and modify string for one rep and one interation 
  Inputs: 
    batch: Estimated batch size (num of gene trees) to re-run simphy 
    simulation_rep: Each replicate ID 
    seed: seed for current run (seed_array[n,m] see utilities.jl)
    iteration: current itereation 
    output_dir: dir to store gene trees simulated by re-running simphy 
    simphy_conf_cont: simphy configuration content (=global simphy_conf_cont)
    matching_pattern: mattching pattern to search modified trees 
                    (here = "g_trees_noLocusID_Gene") -> see utilities.jl
    final_output_per_rep: final dir to store modified newick strings 
  Output: 
    -> Run simphy for one replicate and output simphy files into output_dir
    -> Modify simphy-generated trees based on the criteria above and output 
      to final_output_per_rep 
  """
  function rerun_simphy_1rep_1int(
                        batch::Int, 
                        simulation_rep::Int, 
                        seed::Int32, 
                        iteration::Int, 
                        output_dir::String, 
                        simphy_conf_content::String, 
                        final_modified_trees_output::String, 
                        rootfolder::String, 
                        max_taxa_missing::Int 
                        )

    # Goal 1: run simphy with n_rep = 1
    updated_parameters = """
    -cs $seed  // seed # Update the seed
    -rl f:$batch  // Number of loci (genes) per replicate - f means a fixed value
    """
    updated_content = updated_parameters * simphy_conf_content
    new_conf_file = joinpath(output_dir, "simphysim-conf-Int$iteration-Rep$simulation_rep")
    write(new_conf_file, updated_content)
    
    # Run SimPhy with error handling - throw exception if it fails
    try
      run(`$rootfolder/executables/simphy -i $new_conf_file -o $output_dir`)
    catch e
      # Diagnose the failure by examining locus trees
      diagnostic_msg = diagnose_simphy_locus_trees(output_dir, simulation_rep, iteration)
      # Re-throw with detailed diagnostic information
      error("SimPhy execution failed for Rep$simulation_rep Int$iteration: $e\n" *
            "Diagnostic analysis:\n$diagnostic_msg")
    end

    # modify newick string within one rep
    temp_path = joinpath(output_dir, "1") 
    if isdir(temp_path) 
      simphy_output = temp_path 
    else
      simphy_output = output_dir 
    end

    # check utility.jl for the function modify_newick_for_n_genes 
    tree_repeated_taxa_1int, 
    tree_insufficient_taxa_1int, 
    num_trees_experiencing_gene_loss_only, 
    num_trees_experiencing_gene_duplication_and_loss, 
    num_trees_experiencing_nothing, 
    avg_num_leaf_left_in_trees,
    gene_trees_gene_loss_only,
    gene_trees_gene_duplication_and_loss = modify_newick_for_n_genes(
                                simphy_output, 
                                final_modified_trees_output,
                                batch, 
                                iteration,
                                dup_rate, # indicator helping with file searching 
                                n_inds, # indicator helping with file searching 
                                max_taxa_missing,
                                debug_mode, 
                                ) # see utilities.jl

    return tree_repeated_taxa_1int, 
            tree_insufficient_taxa_1int, 
            num_trees_experiencing_gene_loss_only, 
            num_trees_experiencing_gene_duplication_and_loss, 
            num_trees_experiencing_nothing, 
            avg_num_leaf_left_in_trees,
            gene_trees_gene_loss_only,
            gene_trees_gene_duplication_and_loss
  end 
end

@everywhere begin 
  """
  run_simulation_1rep: 
  Goal: Run simulation for one replicate
    -> Simulate gene trees using SimPhy, modify strings, set seed_array[n,m], and write information
       for one replicate when iteration <= max_iteration_simphy. 
    -> This function calls "rerun_simphy_1rep_1int" to re-run SimPhy until either:
       a) max_iteration_simphy is reached, or
       b) the number of target genes is achieved.

Inputs: 
    max_iteration_simphy: Maximum number of iterations for SimPhy.
    simulation_rep: Replication ID.
    n_genes: Total number of genes to be simulated (target).
    simphy_conf_cont: SimPhy configuration content (global simphy_conf_cont).
    rootfolder: Root folder (the GitHub repository).
    n_genes_min: The minimum number of gene trees to be retained. 
        -> If the number of gene trees simulated < n_genes_min, the loop breaks and raises an error.
    seed_array: The seed array generated by the master random seed.
    folder_path_list: A list to store paths to each replicate folder.

Outputs:
    csv_info: Information about the simulation for this replicate, including:
        - Total number of gene trees generated.
        - Number of iterations run.
        - Number of trees removed due to repeated taxa or insufficient taxa.
        - Average number of leaves left in trees.
        - Counts of trees experiencing gene loss only, gene duplication and gene loss, 
          or no events.
    log_info: Log information for debugging and tracking.
    gene_trees_gene_duplication_and_loss_1rep: 
      List of gene trees experiencing gene duplication and loss.
    gene_trees_gene_loss_only_1rep: 
      List of gene trees experiencing gene loss only.

Notes:
    - The function dynamically removes intermediate SimPhy output folders 
      to save space unless debug_mode is enabled.
    - If no valid trees are generated in an iteration, 
      NaN values are handled gracefully to avoid pipeline errors.
  """
  function run_simulation_1rep(
                      max_iteration_simphy::Int, 
                      simulation_rep::Int, 
                      n_genes::Int, 
                      simphy_conf_content::String, 
                      rootfolder::String, 
                      n_genes_min::Float64, 
                      seed_array::AbstractArray, # both vector or matrix
                      folder_path_list::Vector,
                      max_taxa_missing::Int
                      )

    # Set up dir paths 
    simphyfolder = setup_rep_output_folders(
                  folder_path_list, 
                  simulation_rep, 
                  "genetrees_simphy"
                  ) # each rep has a genetrees_simphy dir to store all outputs from simphy. 
    modified_genetree_folder = setup_rep_output_folders(
                                folder_path_list, 
                                simulation_rep, 
                                "genetrees_singlecopy") 
    # Above: each rep has a genetrees_singlecopy dir to store outputs 

    # Initializing
    mkpath(simphyfolder)
    mkpath(modified_genetree_folder)

    csv_info = ""
    log_info = ""
    iteration = 0 # Inter starts from 0 
    total_trees = 0 # Track total number of gene trees in this rep  
    num_tree_insufficient_taxa = 0 
    num_tree_repeated_taxa = 0 
    enough_genes_status = false

    #= collate information about (gene loss only) and (gene duplication and loss)
    Those information will be used to calculate: 
    1. the percentage of gene loss in each rep
    2. the percentage of gene duplication and loss in each rep
    3. the percentage of nothing happened in each rep
    4. the average number of leaves left in each gene tree
    =#
    tol_num_trees_experiencing_gene_loss_only = 0
    tol_num_trees_experiencing_gene_duplication_and_loss = 0
    tol_num_trees_experiencing_nothing = 0
    avg_num_leaf_left_in_trees_one_rep = []  

    # list to store gene tree name for both categories: 
    gene_trees_gene_duplication_and_loss_1rep = String[] # all gene trees experiencing gene duplication and loss in this rep
    gene_trees_gene_loss_only_1rep = String[] # all gene trees experiencing gene loss only in this rep 

    while iteration < max_iteration_simphy

      # this rep has enough genes already, skip and break the while-loop 
      if enough_genes_status
        break
      end

      # seed for this rep and this int 
      # Below: Interation + 1 because it starts from 0 
      current_seed = seed_array[simulation_rep, iteration + 1] # see utilities.jl
      
      # count num of gene trees in the output dir 
      # Add sync to ensure directory operations are complete
      if isdir(modified_genetree_folder)
        files = readdir(modified_genetree_folder) 
        num = count(file -> occursin(r"g_trees_noLocusID_.*", file), files) 
      else
        num = 0
      end
      # Above: for this rep, count the number of gene trees 

      total_trees = num # Update the total number of trees 

      if num < n_genes 
        # re-running simphy for each rep 
        simphyfolder_int = joinpath(simphyfolder, "Int$iteration")
        mkpath(simphyfolder_int) # each int has its own folder inside simphyfolder
        batch = calculate_batch(num, n_genes) # see utilities.jl 

        # Initialize variables before try-catch to ensure they're in scope
        local tree_repeated_taxa_1int, tree_insufficient_taxa_1int
        local num_trees_experiencing_gene_loss_only, num_trees_experiencing_gene_duplication_and_loss
        local num_trees_experiencing_nothing, avg_num_leaf_left_in_trees
        local gene_trees_gene_loss_only, gene_trees_gene_duplication_and_loss
        
        # Both output below are stats about one iteration
        # Wrap in try-catch to handle SimPhy failures gracefully
        try
          tree_repeated_taxa_1int, 
          tree_insufficient_taxa_1int, 
          num_trees_experiencing_gene_loss_only, 
          num_trees_experiencing_gene_duplication_and_loss, 
          num_trees_experiencing_nothing, 
          avg_num_leaf_left_in_trees,
          gene_trees_gene_loss_only,
          gene_trees_gene_duplication_and_loss = rerun_simphy_1rep_1int(
                                                        batch, 
                                                        simulation_rep, 
                                                        current_seed, 
                                                        iteration, 
                                                        simphyfolder_int, 
                                                        simphy_conf_content, 
                                                        modified_genetree_folder, 
                                                        rootfolder,
                                                        max_taxa_missing
                                                        )
        catch e
          # SimPhy failed for this iteration - log detailed diagnostic information
          error_msg = sprint(showerror, e)
          worker_println("=" ^ 80)
          worker_println("WARNING: Rep$simulation_rep Iter$iteration: SimPhy failed")
          worker_println("=" ^ 80)
          
          # Parse the error message to extract diagnostic information
          if occursin("Diagnostic analysis:", error_msg)
            # Extract and display the diagnostic section
            parts = split(error_msg, "Diagnostic analysis:")
            if length(parts) >= 2
              worker_println("Error: $(strip(parts[1]))")
              worker_println("")
              worker_println("Diagnostic Analysis:")
              worker_println("-" ^ 80)
              worker_println(strip(parts[2]))
            else
              worker_println(error_msg)
            end
          else
            worker_println(error_msg)
          end
          
          worker_println("=" ^ 80)
          worker_println("Skipping this iteration and continuing...")
          worker_println("")
          
          log_info *= "Rep$simulation_rep Iter$iteration: SimPhy crashed - skipped (see worker log for details)\n"
          
          # Set default values for this failed iteration
          tree_repeated_taxa_1int = 0
          tree_insufficient_taxa_1int = 0
          num_trees_experiencing_gene_loss_only = 0
          num_trees_experiencing_gene_duplication_and_loss = 0
          num_trees_experiencing_nothing = 0
          avg_num_leaf_left_in_trees = 0.0  # Use 0.0 instead of NaN for failed iterations
          gene_trees_gene_loss_only = String[]
          gene_trees_gene_duplication_and_loss = String[]
          
          # Don't push NaN/0.0 to the average list for failed iterations
          # This prevents affecting the overall average
          
          # Increment iteration and continue to next iteration
          iteration += 1
          continue
    
        end

        # Collect gene tree names for both categories 
        # Those two list will be uploaded later 
        # This is because we want to remove the gene trees 
        # if the total number of gene trees > n_genes 
        if gene_trees_gene_duplication_and_loss != []
          append!(gene_trees_gene_duplication_and_loss_1rep, gene_trees_gene_duplication_and_loss)
        end                                         
        if gene_trees_gene_loss_only != []
          append!(gene_trees_gene_loss_only_1rep, gene_trees_gene_loss_only)
        end
        
        # Summarize statistics over all interations in one rep 
        # Those are information about this rep collected over interations 
        num_tree_repeated_taxa += tree_repeated_taxa_1int
        num_tree_insufficient_taxa += tree_insufficient_taxa_1int

        # Those below lists will be updated later 
        # Some notes: because we want to remove gene trees if 
        # the total number of gene trees > n_genes 
        # This means that we need to keep track of  the total number of trees for now 
        # The list will be updated later
        tol_num_trees_experiencing_gene_loss_only += num_trees_experiencing_gene_loss_only
        tol_num_trees_experiencing_gene_duplication_and_loss += num_trees_experiencing_gene_duplication_and_loss
        tol_num_trees_experiencing_nothing += num_trees_experiencing_nothing
        
        # Handle NaN in avg_num_leaf_left_in_trees for this iteration
        avg_display = isnan(avg_num_leaf_left_in_trees) ? "N/A (no valid trees)" : 
                      @sprintf("%.2f", avg_num_leaf_left_in_trees)
        
        worker_println("Rep$simulation_rep Iter$iteration: " * 
                "Trees with repeated taxa removed: $tree_repeated_taxa_1int, " * 
                "Trees with insufficient taxa removed: $tree_insufficient_taxa_1int, " * 
                "Trees experiencing gene loss only: $num_trees_experiencing_gene_loss_only, " * 
                "Trees experiencing gene duplication and loss: $num_trees_experiencing_gene_duplication_and_loss, " * 
                "Trees experiencing nothing: $num_trees_experiencing_nothing, " * 
                "Total trees so far: $total_trees / $n_genes, " *
                "Avg leaves: $avg_display")
        
        # Only push valid (non-NaN) values to maintain accurate overall average
        if !isnan(avg_num_leaf_left_in_trees)
          push!(avg_num_leaf_left_in_trees_one_rep, avg_num_leaf_left_in_trees)
        end 

        iteration += 1 # increase one iteration before setting up seed 

      else
        enough_genes_status = true 
      end
    end
    
    # Final check after loop exits - catch silent failures
    files = readdir(modified_genetree_folder)
    final_count = count(file -> occursin(r"g_trees_noLocusID_.*", file), files)
    if final_count == 0 && iteration > 0
      worker_println("ERROR: Rep$simulation_rep completed $iteration iterations but generated 0 trees!")
      worker_println("This may indicate a silent SimPhy failure or file I/O issue.")
    end
    total_trees = final_count  # Update with actual final count

    # Limit gene trees to exactly n_genes if more were generated
    # If the final number of simulated gene trees > n_genes, 
    # Then, we only keep the first n_genes trees after sorting the files 
    files = readdir(modified_genetree_folder) 
    gene_tree_files = filter(file -> occursin(r"g_trees_noLocusID_.*", file), files)
    if length(gene_tree_files) > n_genes
      # Sort files to ensure consistent selection of first n_genes
      sort!(gene_tree_files)
      # Keep only the first n_genes files
      files_to_keep = gene_tree_files[1:n_genes]
      files_to_remove = gene_tree_files[(n_genes+1):end]
      
      # Remove excess gene tree files
      for file_to_remove in files_to_remove

        file_path = joinpath(modified_genetree_folder, file_to_remove)
        rm(file_path; force=true)

        # Remove excess gene trees from genetrees_simphy folder as well... 
        int_id = match(r"Int(\d+)", file_to_remove).captures[1]
        gene_id = match(r"g_trees_noLocusID_Gene(\d+)", file_to_remove).captures[1] 

        temp_path = joinpath(simphyfolder, "Int$int_id", "1") 
        if isdir(temp_path) 
          simphy_tree_file_path = temp_path 
        else
          simphy_tree_file_path = joinpath(simphyfolder, "Int$int_id") 
        end 

        mapsl_file_to_remove = joinpath(simphy_tree_file_path, "$(gene_id)l1g.mapsl") 
        raw_genetree_file_to_remove = joinpath(simphy_tree_file_path, "g_trees$gene_id.trees")
        maplg_file_to_remove = joinpath(simphy_tree_file_path, "$(gene_id)l1g.maplg") 
        # permanently remove those files only if they exist
        if isfile(mapsl_file_to_remove)
            rm(mapsl_file_to_remove; force=true)
        end
        if isfile(raw_genetree_file_to_remove)
            rm(raw_genetree_file_to_remove; force=true)
        end
        if isfile(maplg_file_to_remove)
            rm(maplg_file_to_remove; force=true)
        end

        # IMPORTANT: remove those gene trees from the lists 
        # recording for gene trees with gene loss only + gene duplication and loss
        full_path_to_remove = joinpath(modified_genetree_folder, file_to_remove)
        
        # Check if the file exists in gene_loss_only list and remove it
        if full_path_to_remove in gene_trees_gene_loss_only_1rep
          filter!(x -> x != full_path_to_remove, gene_trees_gene_loss_only_1rep)
          tol_num_trees_experiencing_gene_loss_only -= 1
        # Check if the file exists in gene_duplication_and_loss list and remove it
        elseif full_path_to_remove in gene_trees_gene_duplication_and_loss_1rep
          filter!(x -> x != full_path_to_remove, gene_trees_gene_duplication_and_loss_1rep)
          tol_num_trees_experiencing_gene_duplication_and_loss -= 1
        # If not in either list, it must be a "nothing" tree
        else
          tol_num_trees_experiencing_nothing -= 1
        end 

      end

      # Update total_trees count to reflect the actual number kept
      total_trees = n_genes
      # Have a check function here to make sure the counts are correct 
      total_num_check = tol_num_trees_experiencing_gene_loss_only + 
                        tol_num_trees_experiencing_gene_duplication_and_loss + 
                        tol_num_trees_experiencing_nothing 
      if total_num_check != total_trees # Important check!
        error("Count mismatch after removing excess trees: " * 
              "Expected $total_trees, but got $total_num_check") 
      end

      worker_println("Rep$simulation_rep: Limited to first $n_genes gene trees") 
      worker_println("Removed $(length(files_to_remove)) excess trees")
    end

    # Calculate average number of leaves left in trees
    # Note: NaN values are already filtered out when pushing to the array above
    # This ensures we only average over iterations that produced valid trees
    average_num_leaf_left_in_trees_for_1rep = isempty(avg_num_leaf_left_in_trees_one_rep) ? 
                                               0.0 : mean(avg_num_leaf_left_in_trees_one_rep)
    
    # Print summary for this replicate
    worker_println("Rep$simulation_rep completed: Average leaves per tree = " * 
            (average_num_leaf_left_in_trees_for_1rep == 0.0 ? "N/A (no valid trees)" : 
             @sprintf("%.2f", average_num_leaf_left_in_trees_for_1rep)))

    # Check if n_genes_min is met and write the tracking info 
    if total_trees >= n_genes
      csv_info = join([
                  simulation_rep,
                  total_trees,
                  iteration,
                  num_tree_repeated_taxa,
                  num_tree_insufficient_taxa,
                  "T",
                  average_num_leaf_left_in_trees_for_1rep,
                  tol_num_trees_experiencing_gene_loss_only,
                  tol_num_trees_experiencing_gene_duplication_and_loss,
                  tol_num_trees_experiencing_nothing
                  ], ",") * "\n"

    elseif total_trees >= n_genes_min
      csv_info = join([
                  simulation_rep,
                  total_trees,
                  iteration,
                  num_tree_repeated_taxa,
                  num_tree_insufficient_taxa,
                  "T",
                  average_num_leaf_left_in_trees_for_1rep,
                  tol_num_trees_experiencing_gene_loss_only,
                  tol_num_trees_experiencing_gene_duplication_and_loss,
                  tol_num_trees_experiencing_nothing
                  ], ",") * "\n"

    else # total_trees < n_genes_min 
      num_less = n_genes_min - total_trees 
      info = """
        Warning: Rep$simulation_rep: Insufficient trees ($num_less less than min threshold)\n
        """ 
      log_info *= info 
      # error("Number of genes doesn't meet expectation." * 
      #  "Check simulation_info.log and re-run simulation with updated parameters") 
      csv_info = join([
                  simulation_rep,
                  total_trees,
                  iteration,
                  num_tree_repeated_taxa,
                  num_tree_insufficient_taxa,
                  "F",
                  average_num_leaf_left_in_trees_for_1rep,
                  tol_num_trees_experiencing_gene_loss_only,
                  tol_num_trees_experiencing_gene_duplication_and_loss,
                  tol_num_trees_experiencing_nothing
                  ], ",") * "\n"
    end
  
    return log_info, csv_info, 
          gene_trees_gene_duplication_and_loss_1rep, gene_trees_gene_loss_only_1rep 
  end
end 

"""
run_simulation: 
Goal: -> 1) run simulation with n_rep and < max_iteration_simphy, 
      -> 2) Individual simulation_rep is procssed by different processor
Input: 
  max_iteration_simphy: Maximum iteration for simphy (default = 100)
  n_reps: Total number of replicates
  n_genes: Total number of genes to be simulated (target)
  seed_simphy: The master seed to generate seed_array (see utility.jl)
  simphy_conf_cont: simphy configuration content (=global simphy_conf_cont)
  rootfolder: rootfolder (the github repo)
  n_genes_min: The minimum number of genes trees to be retained. 
    -> If num of gene trees simulated < n_genes_min, the loop breaks and raise an error
  folder_path_list: A list to store path to each rep  
"""

function run_simulation(
          max_iteration_simphy::Int, 
          n_reps::Int, 
          n_genes::Int, 
          seed_array_simphy::Array, 
          simphy_conf_content::String, 
          rootfolder::String, 
          n_genes_min::Float64, 
          folder_path_list::Vector, 
          paramname_root::String
          )

  # Initialize worker logs if logging is enabled
  if log_output
    @everywhere open_worker_log()
  end

  # For simphy simulation, each rep and each iter gets a unique seed
  results = pmap(simulation_rep -> run_simulation_1rep(
                                              max_iteration_simphy, 
                                              simulation_rep, 
                                              n_genes, 
                                              simphy_conf_content, 
                                              rootfolder, 
                                              n_genes_min, 
                                              seed_array_simphy, 
                                              folder_path_list,
                                              max_taxa_missing
                                              ), 1:n_reps)


  # Write info to file: 
  log_info_lst = getindex.(results, 1)  # all first elements
  csv_info_lst = getindex.(results, 2)  # all second elements

  # gene trees experiencing gene duplication and loss
  gene_trees_gene_duplication_and_loss_lst = getindex.(results, 3)
  # gene trees experiencing gene loss only
  gene_trees_gene_loss_only_lst = getindex.(results, 4)
  # Close worker logs if logging is enabled
  if log_output
    @everywhere close_worker_log()
  end

  # Save gene tree list to two csv files: 
  gene_duplication_and_loss_csv = joinpath(outfolder, "Simphy_gene_duplication_and_loss_$paramname_root.csv")
  open(gene_duplication_and_loss_csv, "w") do f
      write(f, "gene_tree_file_path\n")
      for gene_trees in gene_trees_gene_duplication_and_loss_lst
          for tree in gene_trees
              write(f, "$(tree)\n")
          end
      end 
  end

  gene_loss_only_csv = joinpath(outfolder, "Simphy_gene_loss_only_$paramname_root.csv")
  open(gene_loss_only_csv, "w") do f
    write(f, "gene_tree_file_path\n")
      for gene_trees in gene_trees_gene_loss_only_lst
          for tree in gene_trees
              write(f, "$(tree)\n")
          end
      end
  end 

  log_info = vcat(log_info_lst...)
  csv_info = vcat(csv_info_lst...) 

  # A general log file is printed for this parameter settings
  # The log file contains which rep has insufficient num genes (< min threshold)
  # Other scripts (eg. snaq.jl and findgraph.jl) will append warnings into the same log 
  combined_log_info = join(log_info, "")

  #= A csv file is written to record the info about simphy simulation, including: 
    a. RepID;
    b. Number of total gene trees generated per rep; 
    c. Number of iteration run per rep; 
    d. How many gene trees got removed because of duplicated taxa 
    e. How many gene trees got removed because of insufficient taxa (<= 3)  
    This information is saved to keep track of how dup_rate and loss_rate affect simulation
  =# 
  simulation_csv = joinpath(outfolder, "simulation_$paramname_root.csv")

  # If file doesn't exist, create it and write the header
  if !isfile(simulation_csv)
      open(simulation_csv, "w") do f
          write(f, "RepID,n_genes,n_iterations,n_repeated_taxa_removed,"*
          "n_insufficient_taxa_removed,>=min_num_genes"*
          ",avg_num_leaves_left_in_trees,"*
          "num_trees_experiencing_gene_loss_only,"*
          "num_trees_experiencing_gene_duplication_and_loss,"*
          "num_trees_experiencing_nothing\n")
      end
  end

  # Append new csv_info (must be a single string)
  open(simulation_csv, "a") do f
      write(f, join(csv_info, "")) 
      # above ensures to change csv_info from a Vector to String 
  end

  # Return the simphy warnings for integration with concatenation warnings
  return combined_log_info   

end 

#-----------------------------------------------#       
#  simulate molecular sequences using seq-gen
#-----------------------------------------------#
#= Goal: run seq-gen on each replicate and each gene
  Inputs: 
  -Simulation_rep: Replication ID
  -folder_path_list: A list to store path to each rep 
  Outputs: 
  -store concate_alignment_rep{rep_id}.fasta into repID/seqgenfolder/
=# 

@everywhere function run_seqgen_1rep(simulation_rep::Int, 
                                    folder_path_list::Vector, 
                                    seed_array_seqgen::Array, 
                                    seqgen_model::String,
                                    gene_length::Int)  

  rep_number_string = pad_number(simulation_rep, n_reps) # see utilities.jl 
  input_dir = setup_rep_output_folders(folder_path_list, 
                                      simulation_rep, 
                                      "genetrees_singlecopy")  

  output_dir = setup_rep_output_folders(folder_path_list, 
                                        simulation_rep,
                                        joinpath("seqgenfolder", "nexus_folder")
                                        ) 
  # Above: need a separate dir for all nexus because in iqtree.jl, 
  # iqtree is running on all listed files in seqgendir

  mkpath(output_dir)

  # Check if the genetrees_singlecopy directory exists
  # If not, it means the SimPhy stage failed or produced no valid trees
  if !isdir(input_dir)
    worker_println("WARNING: Rep$simulation_rep: genetrees_singlecopy directory not found")
    worker_println("This likely means SimPhy produced no valid gene trees for this replicate")
    return
  end

  genetreefiles = readdir(input_dir)

  seed_idx = 1 # starting from the second seed in the array 
  for input_file_name in genetreefiles 

    # current_seed is updated for each gene tree 
    current_seed = seed_array_seqgen[simulation_rep, seed_idx]
    seed_idx += 1 

    output_file_name = replace(input_file_name, r"\.trees$" => ".nex") 
    input_file_path = joinpath(input_dir, input_file_name)
    output_file_path = joinpath(output_dir, output_file_name)

    if seqgen_model == "all_genes_same_HKY"
      #= to simulate all genes with the same substitution model, use:
      * HKY (-m option) with transition/transversion ratio kappa = 4.143 (option -t)
      * base frequencies 0.316,0.182,0.183,0.319 (-f option)
      * shape alpha = 0.356 (-a option) for the Gamma distribution of rates across sites
      =# 
      kappa, basefreq, alpha = 4.143, [0.316,0.182,0.183,0.319], 0.356 

    elseif seqgen_model == "all_genes_diff_HKY"
      #= to simulate each gene with its own substitution model, use HKY with:
      * kappa from LogNormal(μ=1.4215, σ=0.2798)
      * frequencies from Dirichlet(66.59, 38.41, 38.61, 67.12)
      * alpha from Gamma(α=3.267, θ=0.109).
      =# 
      kappa, basefreq, alpha = sample_substitution_params(current_seed)

    else
      error("Unsupported seqgen_model: $seqgen_model")
    end
    
    A_freq = basefreq[1]
    C_freq = basefreq[2]
    G_freq = basefreq[3]
    T_freq = basefreq[4]

    run(`bash scripts/seq-gen.sh \
      $rep_number_string \
      $input_file_path \
      $output_file_path \
      $current_seed \
      $alpha $kappa $A_freq $C_freq $G_freq $T_freq \
      $gene_length`) 
  end  

end

#-----------------------------------------------#       
#  Convert and concatenate nexus into fasta
#-----------------------------------------------#
#=  Goal 1: concatenate .nex output from seq-gen into a .fasta 
    Goal 2: Store intermediate individual gene (.fasta) in the input folder  
  Inputs: 
  -Simulation_rep: Replication ID
  -folder_path_list: A list to store path to each rep 
  Outputs: 
  -store concate_alignment_rep{rep_id}.fasta into repID/seqgenfolder/
=# 
@everywhere function concatenate_nexus_1rep(
                    simulation_rep::Int, 
                    folder_path_list::Vector)

  rep_number_string = pad_number(simulation_rep, n_reps)
  seqgenfolder = setup_rep_output_folders(
                folder_path_list, 
                simulation_rep, 
                "seqgenfolder"
                )

  input_nexus_folder = joinpath(seqgenfolder, "nexus_folder")

  # Check if the nexus_folder exists and contains files
  if !isdir(input_nexus_folder) || isempty(readdir(input_nexus_folder))
    worker_println("WARNING: Rep$simulation_rep: nexus_folder not found or empty")
    worker_println("Skipping seq-gen concatenation for this replicate")
    return ""
  end

  # Run the Python script and capture output
  result = read(`python3 scripts/concatenate_seq.py \
      $input_nexus_folder \
      $seqgenfolder \
      $rep_number_string`, String)
  
  # Return the message from the Python script
  return result
end

#-----------------------------------------------#       
#  Infer gene trees using IQ-tree
#-----------------------------------------------#
# Specify iqtree and astral temp folder path: 
iqtreefolder_tmp_root = "iqtreefolder"
astralfolder_tmp_root = "astralfolder"

@everywhere global astralfolder_tmp_root = $astralfolder_tmp_root
@everywhere global iqtreefolder_tmp_root = $iqtreefolder_tmp_root 

# temporary folders that need to be 1 level from the repo root folder, for iqtree.pl,
# later moved into their proper place down the folder hierarchy
# Important to have $params-$rep to identify the exact rep and parameter set  
@everywhere function tmp_iqtreeastral_folders(params, rep)  
  return ("$iqtreefolder_tmp_root-$params-$rep", "$astralfolder_tmp_root-$params-$rep") 
end 

# run iqtree on each gene of each rep: use iqtree.pl on each rep
@everywhere function run_iqtree_perl_1rep(
                    simulation_rep::Int, 
                    folder_path_list::Vector, 
                    paramname_root::String, 
                    seed_array_iqtree::Array
                    )
  
  current_seed = seed_array_iqtree[simulation_rep, 1] # get the current seed for this replicate 

  rep_number_string = pad_number(simulation_rep, n_reps) # utilities.jl
  seqgenfolder = setup_rep_output_folders(folder_path_list, 
                simulation_rep, 
                joinpath("seqgenfolder", "nexus_folder"))
  
  # Check if seqgen output exists
  if !isdir(seqgenfolder) || isempty(readdir(seqgenfolder))
    worker_println("WARNING: Rep$simulation_rep: seqgen nexus_folder not found or empty")
    worker_println("Skipping IQ-tree inference for this replicate")
    return
  end
  
  tmpiqtreedir, tmpastraldir =  tmp_iqtreeastral_folders(paramname_root, rep_number_string)

  run(`perl ./scripts/iqtree.pl \
        --seqdir=$seqgenfolder \
        --iqtreedir=$tmpiqtreedir \
        --astraldir=$tmpastraldir \
        --seed_iqtree=$current_seed`)

end

# iqtree.pl requires folders that are 1 level from where the script is run.
# below: move these folders to their proper place
@everywhere function mv_iqtree_folder_1rep(
                      simulation_rep::Int, 
                      folder_path_list::Vector, 
                      paramname_root::String
                      )
  output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")
  iqtreefolder = joinpath(output_folder, "iqtreefolder") 

  rep_number_string = pad_number(simulation_rep, n_reps) # utilities.jl  
  tmpiqtreedir, tmpastraldir =  tmp_iqtreeastral_folders(
                                paramname_root, 
                                rep_number_string)
  
  # Check if temp iqtree directory exists
  if !isdir(tmpiqtreedir)
    worker_println("WARNING: Rep$simulation_rep: temp iqtree folder not found")
    worker_println("Skipping IQ-tree folder move for this replicate")
    return
  end
  
  run(`mv $tmpiqtreedir $iqtreefolder`) 
  # Above: move temp iqtree folder to repID and rename it as iqtreefolder 

end 

#-----------------------------------------------# 
# Reconsurct species tree using astral 
#-----------------------------------------------#  
#= Goal: 
  -> Build a mappingfile which could be used for astral 
  -> Run astral for each rep 
=# 

@everywhere function run_astral_with_mapping_1rep(
                    simulation_rep::Int, 
                    folder_path_list::Vector
                    )
  output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")
  astralfolder = joinpath(output_folder, "astralfolder")
  iqtreefolder = joinpath(output_folder, "iqtreefolder")
  mkpath(astralfolder)
  iqtreefile = joinpath(iqtreefolder, "besttrees.tre")
  astralfile = joinpath(astralfolder, "astral.tre")

  # Check if iqtree output exists
  if !isfile(iqtreefile)
    worker_println("WARNING: Rep$simulation_rep: IQ-tree output file not found")
    worker_println("Skipping ASTRAL species tree inference for this replicate")
    return
  end

  # write mapping file 
  gene_trees = readmultinewick(iqtreefile) # read it into a network obj first 
  # Above: readmultinewick is compatible under PhyloNetworks 1.1.0 
  _,mappingfile = map_accessions_to_species_dict(gene_trees, astralfolder, "astral") 

  # run astral 
  # run(`executables/astral-pro3 -i $iqtreefile -a $mappingfile -o $astralfile`)
  run(`executables/astral -i $iqtreefile -a $mappingfile -o $astralfile`)
end 

#-----------------------------------------------# 
# Clean up the pipeline
#-----------------------------------------------# 
#= goal:
Remove temporary simphy folder to save space 
  -> only if not in debugging mode 
=# 
# The below folder is useful for checking the simphy output 
# so the code is commented out 
# @everywhere function remove_tmp_simphy_folder(
#                     simulation_rep:: Int,
#                     folder_path_list::Vector)
#   output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "") 
#   temp_simphyfolder = joinpath(output_folder, "genetrees_simphy")
#   rm(temp_simphyfolder; force=true, recursive=true)
# end

#= Goal: 
Remove temporary seqgen nexus and phylip folders 
and iqtree gene.treefile to save space
=# 

@everywhere function remove_temp_seqgen_nexus_iqtree_folder(
                    simulation_rep:: Int,
                    folder_path_list::Vector)
  output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "") 
  temp_seqgen_nexus_folder = joinpath(output_folder, "seqgenfolder", "nexus_folder")
  temp_seqgen_phylip_folder = joinpath(output_folder, "seqgenfolder", "phylip")

  # The gene.treefile is important for checking the iqtree output 
  # temp_iqtree_file = joinpath(output_folder, "iqtreefolder", "gene.treefile") # important but besttrees.tre is the same as this file 
  
  # Only remove if they exist
  if isdir(temp_seqgen_nexus_folder)
    rm(temp_seqgen_nexus_folder; force=true, recursive=true)
  end
  if isdir(temp_seqgen_phylip_folder)
    rm(temp_seqgen_phylip_folder; force=true, recursive=true)
  end
  # rm(temp_iqtree_file; force=true)
end

#-----------------------------------------------#       
#  Run the pipeline and time the process
#-----------------------------------------------# 
function main() 
  # Run major parts and time the outputs 
  
  # Set up logging if enabled
  logfile = nothing
  original_stdout = stdout
  
  if log_output
    log_filename = joinpath(outfolder, "screen_$paramname_root.log")
    logfile = open(log_filename, "w")
    println("Logging enabled: Output will be saved to $log_filename")
    println(logfile, "=== Simulation Log Started: $time ===")
    println(logfile, "Parameter set: $paramname_root")
    println(logfile, "")
    flush(logfile)
  end

  @timeit to "Running Simphy and modify newick strings" begin
    seed_array_simphy = seed_generator(
                        seed_simphy, 
                        n_reps, 
                        max_iteration_simphy, 
                        outfolder, 
                        "random_seed_simphy.txt")
    simphy_warnings = run_simulation(
              max_iteration_simphy,
              n_reps, n_genes, 
              seed_array_simphy, 
              simphy_conf_content, 
              rootfolder, 
              n_genes_min, 
              folder_path_list,
              paramname_root)
  end 

  @timeit to "Running Seq-Gen" begin 
    # for seed used in seq-gen, 
    # the each replicate get its own seed array of n_genes seeds 
    # the output is a 2D array: n_reps x n_genes) 
    # when using seq-gen, for each gene tree generated from simphy, 
    # we use [rep_id, gene_id] to get a seed for running this specific seqgen 
    # This seed is also used in sample_substitution_params() function 
    # to get random substitution parameters for this gene 
    # when seqgen_model = "all_genes_diff_HKY" 
    seed_array_seqgen = seed_generator(
                        seed_seqgen, 
                        n_reps, 
                        n_genes, 
                        outfolder, 
                        "random_seed_seqgen.txt")
    pmap(simulation_rep -> run_seqgen_1rep(
                            simulation_rep,
                            folder_path_list, 
                            seed_array_seqgen,
                            seqgen_model,
                            gene_len
                            ),1:n_reps)
  end 

  @timeit to "Modify and concatecate seq-gen sequences to fasta" begin
    concatenate_messages = pmap(simulation_rep -> concatenate_nexus_1rep(
                          simulation_rep, 
                          folder_path_list
                          ), 1:n_reps)
  end 

  @timeit to "Running iqtree" begin 
    seed_array_iqtree = seed_generator(
                        seed_iqtree, 
                        n_reps, 
                        1, 
                        outfolder, 
                        "random_seed_iqtree.txt")
    pmap(simulation_rep -> run_iqtree_perl_1rep(
                            simulation_rep, 
                            folder_path_list, 
                            paramname_root, 
                            seed_array_iqtree
                            ), 1:n_reps)
    pmap(simulation_rep -> mv_iqtree_folder_1rep(
                            simulation_rep, 
                            folder_path_list, 
                            paramname_root
                            ), 1:n_reps) 
  end

  @timeit to "Running astral" begin 
    pmap(simulation_rep -> run_astral_with_mapping_1rep(
                    simulation_rep, 
                    folder_path_list
                    ), 1:n_reps)
  end

  if !debug_mode
      msg = "Debugging mode is OFF."
      println(msg)
      if logfile !== nothing; println(logfile, msg); flush(logfile); end
      
      msg = "Removing temporary seqgen nexus and phylip folders and iqtree gene.treefile to save space..."
      println(msg)
      if logfile !== nothing; println(logfile, msg); flush(logfile); end
      
      pmap(simulation_rep -> remove_temp_seqgen_nexus_iqtree_folder(
                              simulation_rep, 
                              folder_path_list
                              ), 1:n_reps)
  end

  msg = "Simulation for $paramname_root is DONE!"
  println(msg)
  if logfile !== nothing; println(logfile, msg); flush(logfile); end
  
  # Close log file if logging was enabled
  if log_output && logfile !== nothing
    println(logfile, "")
    println(logfile, "=== Simulation Log Ended ===")
    close(logfile)
    println("Log file closed: screen_$paramname_root.log")
  end
  
  # Return both simphy warnings and concatenate messages for logging
  return simphy_warnings, concatenate_messages

end 

if abspath(PROGRAM_FILE) == @__FILE__

  simphy_warnings, concatenate_messages = main()

  # println("----debugging info----")
  # println("simphy_warnings = $simphy_warnings")
  # println("concatenate_messages = $concatenate_messages")

  n_processors = nprocs() # Include this information 
  # The number of processors is -p specified in the script + 1
  # For example julia -p 4 simulation_iqtree.jl --args ..., this will show 5 processes 
  host_name = gethostname() # host name for server 

  # Process concatenation and simphy warnings by rep
  rep_warnings = ""

  for rep_id in 1:n_reps
    rep_warnings_list = String[]
    
    # Add simphy warnings for this rep
    if !isempty(simphy_warnings)
      rep_simphy_warnings = filter(line -> contains(line, "Rep$rep_id:"), split(simphy_warnings, '\n'))
      if !isempty(rep_simphy_warnings)
        append!(rep_warnings_list, rep_simphy_warnings)
      end
    end
    
    # Add concatenation warnings for this rep
    if rep_id <= length(concatenate_messages)
      concat_msg = strip(concatenate_messages[rep_id])
      if !isempty(concat_msg)
        push!(rep_warnings_list, concat_msg)
      end
    end
    
    # Only add rep section if there are warnings
    if !isempty(rep_warnings_list)
      global rep_warnings 
      rep_warnings *= "Rep $rep_id: \n"
      for warning in rep_warnings_list
        rep_warnings *= warning * "\n"
      end
      rep_warnings *= "\n"
    end
  end

  # Write the general log file with appropriate content
  # pointer
  if !isempty(simphy_warnings) 
    simphy_log_warnings = "simphy_warnings\n"
  else
    simphy_log_warnings = "All genes meet min target number of genes\n"
  end

  # save all parameters and write into the output dir
  arguments = """
  #==============================================================#
  #--------------SimPhy, SeqGen, IQtree, Astral------------------#
  #==============================================================#
  ---Arguments used for simulating this output dataset---
  duplication rate = $dup_rate;
  loss rate = $loss_rate; 
  rate variation = $ratevar; 
  number of replicates = $n_reps; 
  number of genes = $n_genes; 
  number of individuals per taxon = $n_inds; 
  scaling factor for effective population size = $SF, Ne = $scaled_Ne;
  seqgen substitution model = $seqgen_model;
  gene length (num of bases per gene simulated by seq-gen) = $gene_len;
  species tree used:
  $species_tree;  

  ---Seed Information---
  master seed = $master_seed, used to generate seeds for simphy, seqgen, and iqtree;
  seed_simphy = $seed_simphy, used to generate seed array (rep x iter) used for Simphy; 
  seed_seqgen = $seed_seqgen, used to generate seed array (rep x 1) used for seqgen; 
  seed_iqtree = $seed_iqtree, used to generate seed array (rep x 1) used for iqtree; 

  ---Rerunning SimPhy--- 
  max_iteration_simphy = $max_iteration_simphy, maximum iteration to re-run simphy; 
  n_genes_min = $n_genes_min, num of simulated gene trees >= n_genes_min; 

  ---Other information---
  Number of processors used = $n_processors;
  Server for running the script = $host_name.stat.wisc.edu

  Time of running the script = $time;
  ---Important Warnings to double check before processing---

  $rep_warnings
  $simphy_log_warnings
  ---Running time---
  """
  argument_file = joinpath(outfolder, "arguments-$paramname_root.log")
  write(argument_file, arguments)
  open(argument_file, "a") do io
    show(io, to)  # Write timer output
  end

end

