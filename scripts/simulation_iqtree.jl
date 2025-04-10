#=
this script does:

1. Simulate gene trees based on the grand species tree using SimPhy 
    * After simulation, newick strings are modified to mimic hidden paralogy. 
      Tree not meeting the criteria (see below) are not used for downstream pipeline 
        Thus, simphy are re-run (<= max_iteration_simphy times) to get enough gene_trees. 
        Error raises if not enough trees (< n_genes_min) simulated 
2. Simulate molecular sequences using Seq-Gen based on gene trees output from SimPhy
3. Estmate gene trees using iqtree and infer species trees using Astral 

example: run this below
 on darwin cluster use the `/nobackup` dir with tmux
 on franlin use `/nobackup2`
 on UW-Madison Botany pink server 
 from the root of the repository:
julia scripts/simulation.jl

output: it will create folder 'output' in the root of the repo with inside:
-rep ID: Each folder is labeled as repID 
  Inside the repID folder: 
  - sim-phy-outfiles: Simulated gene trees and configuration files 
  - seq-gen-outfiles: Simulated molecular sequences 
  - iqtree-outfiles: Estimated gene trees 
  - astral-outfiles: Inferred species tree
=#

using Distributed  
using ArgParse
using TimerOutputs
@everywhere include("utilities.jl")
const to = TimerOutput() 

function parse_commandline()
  s = ArgParseSettings() 
  @add_arg_table s begin 

    # Specify parameters for SimPhy: 
    "--dup_rate"
      help = "Gene duplication rate: specify the gene duplication rate, if 0, then no duplication"
      arg_type = Float64
      required = true
    "--loss_rate" 
      help = "Gene loss rate: specify the gene loss rate, if 0, then no gene loss"
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
    "--max_iteration_simphy"
      help = "Maximum iteration to re-run simphy to get enough gene trees"
      arg_type = Int
      default = 100 
    "--min_gene_porportion"
      help = "Number of simulated gene trees in each rep >= min_gene_proportion (0 to 1) * n_genes"
      arg_type = Float64
      default = 0.8 # let set the default as 0.8 hard coded 
  end 
  return parse_args(s) 
end

#-----------------------------------------------#       
#         Initialization
#-----------------------------------------------#
#--------------- Parse arguments ---------------# 
parsed_args = parse_commandline()
dup_rate = parsed_args["dup_rate"] # A number indicates rate of gene duplication rates (0 = no dup)
loss_rate = parsed_args["loss_rate"] # A number indicates rate of gene loss rates (0 = no loss)
ratevar = parsed_args["ratevar"] # "N", "G", "L" or "GL" or "G*L" to include genexlineage rate variation 
n_reps = parsed_args["n_reps"] # number of replicates 
n_genes = parsed_args["n_genes"] # number of genes 
n_inds = parsed_args["n_inds"] # Number of individuals per taxa -- default = 1
max_iteration_simphy = parsed_args["max_iteration_simphy"] # Maximum number of iteration of re-running simphy 
min_gene_porportion = parsed_args["min_gene_porportion"] # Percentage * n_genes =  n_genes_min, the mimum number of gene trees in each rep 
n_genes_min = n_genes * min_gene_porportion  # After reaching max_iteration_simphy, -only if num of genes > n_genes_min --> proceed  

#--------------- Set up configuration ---------------# 
# Check if ratevar belongs to the following: 
valid_ratevars = ["N", "G", "L", "GL", "GxL", "G*L"] 
if !(ratevar in valid_ratevars) 
  error("Invalid value for --ratevar: $ratevar. Valid options are: $(join(valid_ratevars, ", "))")
end

# SimPhy sets dup_rate > loss_rate, or otherwise there will be error message 
if dup_rate < loss_rate
  error("Invalid loss rate $loss_rate -- have it lower than dup_rate $dup_rate")
end

# set all configuration parameters here 
rootfolder = pwd()
paramname_root = "DUP$dup_rate-LOS$loss_rate-RV$ratevar-N_ind$n_inds"
outfolder = joinpath(rootfolder, "output", paramname_root) 

#--------------- set up seeds ---------------# 
params_dict_for_seed_setting = get_dict_for_seed_setting(paramname_root)
master_seed = generate_master_seed(params_dict_for_seed_setting) # master seed used for this script

software_names = ["simphy", "seqgen", "iqtree"] 
# set up seeds for all softwares 
seed_dic = generate_software_seeds(master_seed, software_names) # see utility.jl 
seed_simphy = seed_dic["simphy"] 
seed_seqgen = seed_dic["seqgen"]
seed_iqtree = seed_dic["iqtree"]
# astral is a deterministic program with no need to specify a seed 

#--------------- set up folders ---------------# 
# The below function checks if a path exists, if yes, the function will ask user's input tp decide if to remove the path 
check_existing_dir([outfolder]) # see utilities.jl --> input needs to be a list
mkpath(outfolder) # create a new one if there is no pre-existing outfolder

folder_path_list = [] # create folder path 
for n in 1:n_reps # All outputs goes to each rep file 
  rep_number_string = pad_number(n, n_reps)
  rep_folder_path = joinpath(outfolder, "rep$rep_number_string")
  push!(folder_path_list, rep_folder_path)
end 

#-----------------------------------------------#       
#         set up SimPhy parameters
#-----------------------------------------------#
# Modify the master simphy config file based on arguments and then save it to new config file
master_conf = joinpath(rootfolder,"simphy-configs/", "simphysim-conf-master")
conf_content = read(master_conf, String) # read the master config file into a string 

#Set up the parameters
parameters = ""

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
end

# To simulate variation across lineages (ratevar = "L" or "GL" or "GxL") 
if occursin("L", ratevar) # add tree with variation across lineages
  parameters *= "-s (A:3.44*0.0100947,((((B:0.88*0.0042057,C:0.88*0.0036776):1.71*0.0078509,(D:0.93*0.0235933,E:0.93*0.0199793):1.66*0.0079913):0.17*0.0068836,F:2.76*0.0067212):0.18*0.0098089,(G:0.5*0.0797969,H:0.5*0.1796924):2.44*0.0190487):0.5*0.0694588); // tree with lineage variation\n"
else 
  parameters *= "-s (A:3.44,((((B:0.88,C:0.88):1.71,(D:0.93,E:0.93):1.66):0.17,F:2.76):0.18,(G:0.5,H:0.5):2.44):0.5); // tree without lineage variation\n"
  # If ratevar doesn't contain L (G and N), then add tree without variations across lineages 
end

# To simulate multiple individuals per species 
if n_inds > 1 # if n_inds == 1, no "-si" argument & use the default 
  parameters *= "-si f:$n_inds // number of individuals per tree tip\n" 
end 

# Store all parameters into one string, which will be updated in the below loop
simphy_conf_content = parameters * conf_content # combine parameters with master config

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

#-----------------------------------------------#       
#  run SimPhy + modify newick strings 
#-----------------------------------------------#
#= The below function has two goals: 
    1. simulate gene trees using Simphy (rep = 1 and loop through n_reps)
    2. modify newick strings to mimic hidden paralogy 
  If the gene tree has
  * more than one repeated gene copy: it is obvious that gene duplication happened
    so paralogy is not hidden. There will be *no* associated gene tree file.
  * some gene loss events and there are <= 3 taxa left: not phylogenetically
    informative. There will be *no* associated gene tree file.
  * 0 or 1 gene copy per individual, and >= 4 taxa left: if there is paralogy,
    then it's hidden. The associated gene tree file is the same as the original
    except for taxon names (gene copy number is removed, the taxon name still
    contains the species + individual number).
  If dup_rate = 0, then all gene copies are orthologous: there cannot be any
  hidden paralogy. The below code will only change the tip names in the
  new gene tree files. 
  * This newick string modification might remove some trees from the pipeline. 
    a. If the resulting num of gene trees <= n_genes (target), we re-run simphy with <= max_iteration_simphy times to generate enough gene trees. 
    b. If iteration > max_iteration_simphy, num of gene trees > n_genes_min -> pipeline continues 
    c. If iteration > max_iteration_simphy, num of gene tress < n_genes_min -> pipeline breaks 
=#

"""
rerun_simphy_1rep_1int: Goal: Run simphy and modify string for one replicate and one particular interation 
Inputs: 
  batch: Estimated batch size (num of gene trees) to re-run simphy 
  simulation_rep: Each replicate ID 
  seed: seed for current run (seed_array[n,m] see utilities.jl)
  iteration: current itereation 
  output_dir: dir to store gene trees simulated by re-running simphy 
  simphy_conf_cont: simphy configuration content (=global simphy_conf_cont)
  matching_pattern: mattching pattern to search modified trees (here = "g_trees_noLocusID_Gene") -> see utilities.jl
  final_output_per_rep: final dir to store modified newick strings 
Output: 
  -> Run simphy for one replicate and output simphy files into output_dir
  -> Modify simphy-generated trees based on the criteria above and output to final_output_per_rep 
"""

@everywhere function rerun_simphy_1rep_1int(batch::Int, simulation_rep::Int, seed::Int, iteration::Int, output_dir::String, simphy_conf_content::String, final_modified_trees_output::String, rootfolder::String) 
  # Goal 1: run simphy with n_rep = 1
  updated_parameters = """
  -cs $seed  // seed # Update the seed
  -rs 1  // Number of replicates
  -rl f:$batch  // Number of loci (genes) per replicate - f means a fixed value
  """
  updated_content = updated_parameters * simphy_conf_content
  new_conf_file = joinpath(output_dir, "simphysim-conf-Int$iteration-Rep$simulation_rep")
  write(new_conf_file, updated_content)
  run(`$rootfolder/executables/simphy -i $new_conf_file -o $output_dir`)

  # modify newick string within one rep 
  simphy_output = joinpath(output_dir, "1")
  modify_newick_for_n_genes(simphy_output, final_modified_trees_output, batch, iteration)  # see utilities.jl 
end 

"""
run_simulation_1rep: 
Goal: Ren simulation (running simphy, modify strings, set seed_array[n,m], write information) for one replicate when iteration <= max_iteration_simphy. 
  This function calls the previous function to re-run until max_iteration_simphy reaches or the number of target genes reaches. 
Inputs: 
  max_iteration_simphy: Maximum iteration for simphy (default = 100)
  simulation_rep: Replication ID
  n_genes: Total number of genes to be simulated (target)
  simphy_conf_cont: simphy configuration content (=global simphy_conf_cont)
  rootfolder: rootfolder (the github repo)
  n_genes_min: The minimum number of genes trees to be retained. 
    -> If num of gene trees simulated < n_genes_min, the loop breaks and raise an error
  seed_array: the seed array generated by the master random seed
  folder_path_list: A list to store path to each rep  
"""

@everywhere function run_simulation_1rep(max_iteration_simphy::Int, simulation_rep::Int, n_genes::Int, simphy_conf_content::String, rootfolder::String, n_genes_min::Float64, seed_array, folder_path_list)

  # Set up dir paths 
  simphyfolder = setup_rep_output_folders(folder_path_list, simulation_rep, "genetrees_simphy") # each rep has a genetrees_simphy dir to store all outputs from simphy. 
  modified_genetree_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "genetrees_singlecopy") # each rep has a genetrees_singlecopy dir to store all outputs from modify_newick (utilities.jl)
  mkpath(simphyfolder)
  mkpath(modified_genetree_folder)

  # Initializing 
  info = "" # info for rep will be appened to tracking_info 
  iteration = 1 # Interation starts from 0 
  total_trees = 0 # Track total number of gene trees in this rep  
  enough_genes_status = false

  while iteration <= max_iteration_simphy

    # seed for this rep and this int 
    current_seed = seed_array[simulation_rep, iteration] # see utilities.jl 
    
    # this rep has enough genes already, skip and break the while-loop 
    if enough_genes_status
      info = String("This rep gets enough genes! Loop breaks at iteration $iteration. ")
      break
    end
    
    # count num of gene trees in the output dir 
    files = readdir(modified_genetree_folder) 
    num = count(file -> occursin(r"g_trees_noLocusID_.*", file), files) # For this rep, count num of gene trees 
    total_trees = num # Update the total number of trees 

    if num < n_genes 
      # re-running simphy for each rep 
      simphyfolder_int = joinpath(simphyfolder, "Int$iteration")
      mkpath(simphyfolder_int) # each int has its own folder inside simphyfolder
      batch = calculate_batch(num, n_genes) # see utilities.jl 
      rerun_simphy_1rep_1int(batch, simulation_rep, current_seed, iteration, simphyfolder_int, simphy_conf_content, modified_genetree_folder, rootfolder)

      iteration += 1 # increase one iter
    else
      enough_genes_status = true
    end
  
    # Check if n_genes_min is met and write the tracking info 
    if total_trees >= n_genes
      num_more_trees = total_trees - n_genes
      if enough_genes_status 
        str = ("""
        --------------- Rep$simulation_rep --------------------
        Number of gene trees reaches the target. 
        There are $num_more_trees more compared to the target $n_genes. 
        Total number of gene trees: $total_trees. 
        """)
        info = str * info 
      else 
        info = """
        --------------- Rep$simulation_rep --------------------
        Maximum iteration $max_iteration_simphy reaches and the bumber of gene trees reaches the target. 
        There are $num_more_trees more compared to the target $n_genes. 
        Total number of gene trees: $total_trees. """
      end

      elseif total_trees >= n_genes_min
        num_less_trees = n_genes - total_trees 
        num_more_than_lower = total_trees - n_genes_min
        info = """
          --------------- Rep$simulation_rep --------------------
          Maximum iteration $max_iteration_simphy reaches BUT there are $num_less_trees less trees than the target $n_genes. 
          $num_more_than_lower more gene trees HIGHER than n_genes_min.
          Total number of gene trees: $total_trees. """
      else
        num_less_than_lower = n_genes_min - total_trees 
        info = """
          --------------- Rep$simulation_rep --------------------
          Maximum iteration $max_iteration_simphy reaches BUT $num_less_than_lower less than n_genes_min. Simulation stopd. 
          Check parameters and re-run the simulation. Total number of gene trees: $total_trees. """ 
        # error("Number of genes doesn't meet expectation. Check simulation_info.log and re-run simulation with updated parameters") # The script stops if # gene trees in i rep doesn't reach the n_genes_min
      end
    end

  return info
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
 
function run_simulation(max_iteration_simphy::Int, n_reps::Int, n_genes::Int, seed_array_simphy::Array, simphy_conf_content::String, rootfolder::String, n_genes_min::Float64, folder_path_list::Vector)

  # For simphy simulation, each replicate and each iteration gets a unique seed based on seed_array_simphy generated from seed_simphy

  tracking_info_list= pmap(simulation_rep -> run_simulation_1rep(max_iteration_simphy, simulation_rep, n_genes, simphy_conf_content, rootfolder, n_genes_min, seed_array_simphy, folder_path_list), 1:n_reps)
  tracking_info = vcat(tracking_info_list...)

  # Write tracking info to file
  tracking_file_path = joinpath(outfolder, "simulation_tracking_file.log")
  open(tracking_file_path, "a") do io
      write(io, join(tracking_info, "\n"))
  end 
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

@everywhere function run_seqgen_1rep(simulation_rep::Int, folder_path_list::Vector, seed_array_seqgen::Array) 

  current_seed = seed_array_seqgen[simulation_rep, 1]

  rep_number_string = pad_number(simulation_rep, n_reps) # see utilities.jl 
  input_dir = setup_rep_output_folders(folder_path_list, simulation_rep, "genetrees_singlecopy") 
  output_dir = setup_rep_output_folders(folder_path_list, simulation_rep,  joinpath("seqgenfolder", "nexus_folder")) # Need a separate folder for all nexus because in iqtree.jl, iqtree is running on all listed files in seqgendir
  mkpath(output_dir)

  genetreefiles = readdir(input_dir)
  for input_file_name in genetreefiles 
    output_file_name = replace(input_file_name, r"\.trees$" => ".nex") 
    input_file_path = joinpath(input_dir, input_file_name)
    output_file_path = joinpath(output_dir, output_file_name)
    run(`bash scripts/seq-gen.sh $rep_number_string $input_file_path $output_file_path $current_seed`) 
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
@everywhere function concatenate_nexus_1rep(simulation_rep::Int, folder_path_list::Vector)
  rep_number_string = pad_number(simulation_rep, n_reps)
  seqgenfolder = setup_rep_output_folders(folder_path_list, simulation_rep, "seqgenfolder") # output folder. Save the concated fasta one layer above nexus_folder 
  input_nexus_folder = joinpath(seqgenfolder, "nexus_folder")
  run(`python3 scripts/concatenate_seq.py $input_nexus_folder $seqgenfolder $rep_number_string`) 
end

#-----------------------------------------------#       
#  simulate iqtree and astral using seq-gen
#-----------------------------------------------#
# Specify iqtree and astral temp folder path: 
iqtreefolder_tmp_root = "iqtreefolder"
iqtreefolder  = joinpath(outfolder, iqtreefolder_tmp_root)
astralfolder_tmp_root = "astralfolder"
astralfolder = joinpath(outfolder, astralfolder_tmp_root)

@everywhere global astralfolder_tmp_root = $astralfolder_tmp_root
@everywhere global iqtreefolder_tmp_root = $iqtreefolder_tmp_root 

# temporary folders that need to be 1 level from the repo root folder, for iqtree.pl,
# later moved into their proper place down the folder hierarchy
# Important to have $params-$rep to identify the exact rep and parameter set  
@everywhere function tmp_iqtreeastral_folders(params, rep)  
  return ("$iqtreefolder_tmp_root-$params-$rep", "$astralfolder_tmp_root-$params-$rep") 
end 

# run iqtree on each gene of each rep: use iqtree.pl on each rep
@everywhere function run_iqtree_perl_1rep(simulation_rep::Int, folder_path_list::Vector, paramname_root::String, seed_array_iqtree::Array)
  
  current_seed = seed_array_iqtree[simulation_rep, 1] # get the current seed for this replicate 

  rep_number_string = pad_number(simulation_rep, n_reps) # utilities.jl
  seqgenfolder = setup_rep_output_folders(folder_path_list, simulation_rep, joinpath("seqgenfolder", "nexus_folder"))
  tmpiqtreedir, tmpastraldir =  tmp_iqtreeastral_folders(paramname_root, rep_number_string)

  run(`perl ./scripts/iqtree.pl --seqdir=$seqgenfolder --iqtreedir=$tmpiqtreedir --astraldir=$tmpastraldir --seed_iqtree=$current_seed`)

end

# iqtree.pl requires folders that are 1 level from where the script is run.
# below: move these folders to their proper place
@everywhere function mv_iqtree_astral_folder_1rep(simulation_rep::Int, folder_path_list::Vector, paramname_root::String)
  output_folder = setup_rep_output_folders(folder_path_list, simulation_rep, "")
  rep_number_string = pad_number(simulation_rep, n_reps) # utilities.jl  
  tmpiqtreedir, tmpastraldir =  tmp_iqtreeastral_folders(paramname_root, rep_number_string)
  run(`mv $tmpiqtreedir $output_folder/iqtreefolder`) # move temp iqtree folder to repID and rename it as iqtreefolder 
  run(`mv $tmpastraldir $output_folder/astralfolder`) # move temp astral folder to repID and rename it as astralfolder 
end 

#-----------------------------------------------#       
#  Run the pipeline and time the process
#-----------------------------------------------# 
function main() 
  # Runn major parts and time the outputs 

  @timeit to "Running Simphy and modify newick strings" begin
    seed_array_simphy = seed_generator(seed_simphy, n_reps, max_iteration_simphy, outfolder, "random_seed_simphy.txt")
    run_simulation(max_iteration_simphy, n_reps, n_genes, seed_array_simphy, simphy_conf_content, rootfolder, n_genes_min, folder_path_list)
  end 

  @timeit to "Running Seq-Gen" begin 
    seed_array_seqgen = seed_generator(seed_seqgen, n_reps, 1, outfolder, "random_seed_seqgen.txt")
    pmap(simulation_rep -> run_seqgen_1rep(simulation_rep,folder_path_list, seed_array_seqgen), 1:n_reps) 
  end 

  @timeit to "Concatecate seq-gen sequences to fasta" begin
    pmap(simulation_rep -> concatenate_nexus_1rep(simulation_rep, folder_path_list), 1:n_reps)
  end 

  @timeit to "Running iqtree + astral" begin 
    seed_array_iqtree = seed_generator(seed_iqtree, n_reps, 1, outfolder, "random_seed_iqtree.txt")
    pmap(simulation_rep -> run_iqtree_perl_1rep(simulation_rep, folder_path_list, paramname_root, seed_array_iqtree), 1:n_reps)
    pmap(simulation_rep -> mv_iqtree_astral_folder_1rep(simulation_rep, folder_path_list, paramname_root), 1:n_reps) 
  end

  println("Simulation for $paramname_root is DONE!")

end 

if abspath(PROGRAM_FILE) == @__FILE__
  main()

  n_processors = nprocs() # Include this information 
  # The number of processors is -p specified in the script + 1
  # For example julia -p 4 simulation_iqtree.jl --args ..., this will show 5 processes 
  host_name = gethostname() # host name for server 

  # save all parameters and write into the output dir
  arguments = """
  #--------------SimPhy, SeqGen, IQtree, Astral------------------#
  Arguments used for simulating this output dataset: 
  duplication rate = $dup_rate;
  loss rate = $loss_rate; 
  rate variation = $ratevar; 
  number of replicates = $n_reps; 
  number of genes = $n_genes; 
  master seed = $master_seed, this will be further used to generate three seeds for simphy, seqgen, and iqtree;
  seed_simphy = $seed_simphy, used to generate seed array (replicate x iteration) used for Simphy; 
  seed_seqgen = $seed_seqgen, used to generate seed array (replicate x 1) used for seqgen; 
  seed_iqtree = $seed_iqtree, used to generate seed array (replicate x 1) used for iqtree; 
  number of individuals per taxon = $n_inds; 
  max_iteration_simphy = $max_iteration_simphy, maximum iteration to re-run simphy; 
  n_genes_min = $n_genes_min, num of simulated gene trees >= n_genes_min else pipeline stops; 
  Number of processors used = $n_processors;
  Server for running the script = $host_name.stat.wisc.edu

  Running time shown below: 
  """
  argument_file = joinpath(outfolder, "arguments-$paramname_root.log")
  write(argument_file, arguments)
  open(argument_file, "a") do io
    show(io, to)  # Write timer output
  end

end



