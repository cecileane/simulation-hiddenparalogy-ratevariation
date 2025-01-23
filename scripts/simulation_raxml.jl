#=
this script does:

1. Simulate gene trees based on the grand species tree using SimPhy 
2. Simulate molecular sequences using Seq-Gen based on gene trees output from SimPhy
3. Estmate gene trees using RAxML and infer species trees using Astral 

example: run this below
 on darwin cluster use the `/nobackup` dir with tmux
 on franlin use `/nobackup2`
 on UW-Madison Botany pink server 
 from the root of the repository:
julia scripts/simulation.jl

output: it will create folder 'output' in the root of the repo with inside:
- sim-phy-outfiles: Simulated gene trees and configuration files 
- seq-gen-outfiles: Simulated molecular sequences 
- raxml-outfiles: Estimated gene trees 
- astral-outfiles: Inferred species tree
=#

using ArgParse 
include("utilities.jl")

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
      help = "'N': No rate variation';\n 'G: Gene specific rate variation';\n 'L': Lineage specific rate variation;\n 'GL' or 'G*L': genexlineage rate variation"
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
    "--seed_simphy"
      help = "Seed for simphy"
      arg_type = Int
      required = true
    "--n_inds" 
      help = "Number of individuals/accessions per species (Default = 1)"  
      arg_type = Int
      default = 1
  end 
  return parse_args(s)
end

parsed_args = parse_commandline()

# Parameter settings for SimPhy
dup_rate = parsed_args["dup_rate"] # A number indicates rate of gene duplication rates (0 = no dup)
loss_rate = parsed_args["loss_rate"] # A number indicates rate of gene loss rates (0 = no loss)
ratevar = parsed_args["ratevar"] # "N", "G", "L" or "GL" or "G*L" to include genexlineage rate variation 
n_reps = parsed_args["n_reps"] # number of replicates 
n_genes = parsed_args["n_genes"] # number of genes 
seed_simphy = parsed_args["seed_simphy"] # seeds for SimPhy 
n_inds = parsed_args["n_inds"] # Number of individuals per taxa -- default = 1

# Check if ratevar belongs to the following: 
valid_ratevars = ["N", "G", "L", "GL", "GxL"] 
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
outfolder = "output/$paramname_root"
mkpath(outfolder)
simphyfolder = joinpath(outfolder, "sim-phy-outfiles")
seqgenfolder = joinpath(outfolder, "seq-gen-outfiles")
raxmlfolder_tmp_root = "raxml-outfiles"
raxmlfolder  = joinpath(outfolder, raxmlfolder_tmp_root)
astralfolder_tmp_root = "astral-outfiles"
astralfolder = joinpath(outfolder, astralfolder_tmp_root)
# temporary folders that need to be 1 level from the repo root folder, for raxml.pl,
# later moved into their proper place down the folder hierarchy
function tmp_raxmlastral_folders(params, rep)  
  return ("$raxmlfolder_tmp_root-$params-$rep", "$astralfolder_tmp_root-$params-$rep")
end

mkdir(simphyfolder)
mkdir(seqgenfolder)
mkdir(raxmlfolder)
mkdir(astralfolder)

# For future references, save all arguments used for this dataset and write it to the output folder
arguments = """
Arguments used for this output dataset: 
duplication rate = $dup_rate, 
loss rate = $loss_rate, 
rate variation = $ratevar, 
number of replicates = $n_reps, 
nnumber of genes = $n_genes, 
seed for SimPhy= $seed_simphy, 
number of individuals per taxon = $n_inds
"""
argument_files = joinpath(outfolder, "arguments-$paramname_root")
write(argument_files, arguments)

# Modify the master simphy config file based on arguments and then save it to new config file
master_conf = joinpath(rootfolder,"simphy-configs/", "simphysim-conf-master")
conf_content = read(master_conf, String) # read the master config file into a string 

#Set up the parameters (# replicates, # genes, and seeds)
parameters = """
# Parameters: DUP$dup_rate-LOS$loss_rate-RV:$ratevar, with replicate = $n_reps genes = $n_genes and seeds = $seed_simphy:
-rs $n_reps  // Number of replicates
-rl f:$n_genes  // Number of loci (genes) per replicate - f means a fixed value 
-cs $seed_simphy  // seed
"""

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
if n_inds > 1 # if n_inds == 1, nothing is padded so simphy simulated 1 ind per species
  parameters *= "-si f:$n_inds // number of individuals per tree tip\n" 
end 

# change new_conf_file 
combined_content = parameters * conf_content # combine parameters with master config
new_conf_file = joinpath(simphyfolder, "simphysim-conf-$paramname_root")
write(new_conf_file, combined_content) # write the combined config into SimPhy output folder

#-----------------------------------------------#       
#  simulate gene trees using SimPhy 
#-----------------------------------------------#
run(`$rootfolder/executables/simphy -i $new_conf_file -o $simphyfolder`)

#-----------------------------------------------#       
#  Modify the gene trees generated by SimPhy 
#-----------------------------------------------#
# The goal of this section is to modify the gene trees generated by simphy to mimic hidden paralogy 
# If the gene tree has more than one repeated gene copy, it is obvious that gene duplication happened so not hidden. We will not process those treees. 
# If the gene tree has many gene loss events and there are <= 3 taxa left, not possibly to process those trees. We will not process those trees. 
# When we don't want to simulate hidden paraology, then dup_rate = 0 and loss_rate = 0, then the below code will only change the tree tips. 
for simulation_rep in 1:n_reps
  rep_number_string = lpad(string(simulation_rep), ceil(Int, log10(n_reps+1)), '0')
  simphy_rep_path = joinpath(simphyfolder, rep_number_string)

  for gene_tree in 1:n_genes
    genenum_string = lpad(string(gene_tree), ceil(Int, log10(n_genes+1)), '0')
    input = joinpath(simphy_rep_path, "g_trees$genenum_string.trees")
    output = joinpath(simphy_rep_path, "g_trees_noLocusID_$genenum_string.trees")
    modify_newicks(input, output)
  end
end

#-----------------------------------------------#       
#  simulate molecular sequences using seq-gen
#-----------------------------------------------#
# run seq-gen on each replicate and each gene
for simulation_rep in 1:n_reps
  rep_number_string = lpad(string(simulation_rep), ceil(Int, log10(n_reps+1)), '0')
  run(`mkdir $seqgenfolder/rep$rep_number_string`)
  for gene_tree in 1:n_genes
    # account for sim phy naming convention
    genenum_string = lpad(string(gene_tree), ceil(Int, log10(n_genes+1)), '0')
    # seq-gen.sh knows about the directory structure & names
    run(`bash scripts/seq-gen.sh $rep_number_string $genenum_string $outfolder`)
  end
end

#-----------------------------------------------#       
#  Convert and concatenate nexus into fasta
#-----------------------------------------------#
# After simulating individual genes, the below code call concatenate_seq.py in root/scripts
# First converts the nexus output from seq-gen and then concatenate the sequences into a concatenated fasta 
# Intermediate .fasta files for each indivusal genes which will be stored in the same input folder
# The below code is used for the codes I wrote using concatenate_seq.py 
for simulation_rep in 1:n_reps
  rep_number_string = lpad(string(simulation_rep), ceil(Int, log10(n_reps+1)), '0') 
  input_nexus_dir = joinpath(seqgenfolder,"rep$rep_number_string")
  output_fasta_dir = joinpath(seqgenfolder)
  run(`python scripts/concatenate_seq.py $input_nexus_dir $output_fasta_dir`) 
end

#-----------------------------------------------#       
#  simulate raxml and astral using seq-gen
#-----------------------------------------------#
# run raxml on each gene of each rep: use raxml.pl on each rep
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  tmpraxmldir, tmpastraldir =  tmp_raxmlastral_folders(paramname_root, simulation_rep)
  run(`perl ./scripts/raxml.pl --seqdir=$seqgenfolder/rep$simulation_rep --raxmldir=$tmpraxmldir --astraldir=$tmpastraldir`)
end

# raxml.pl requires folders that are 1 level from where the script is run.
# below: move these folders to their proper place
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  tmpraxmldir, tmpastraldir =  tmp_raxmlastral_folders(paramname_root, simulation_rep)
  run(`mv $tmpraxmldir  $raxmlfolder`)
  run(`mv $tmpastraldir $astralfolder`)
end