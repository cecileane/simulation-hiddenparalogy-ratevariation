#=
this script does:

1. Simulate gene trees based on the grand species tree using SimPhy 
2. Simulate molecular sequences using Seq-Gen based on gene trees output from SimPhy
3. Estmate gene trees using RAxML and infer species trees using Astral 

example: run this below
 on darwin cluster use the `/nobackup` dir with tmux
 on franlin use `/nobackup2`
 from the root of the repository:
julia scripts/simulation.jl

output: it will create folder 'output' in the root of the repo with inside:
- sim-phy-outfiles: Simulated gene trees and configuration files 
- seq-gen-outfiles: Simulated molecular sequences 
- raxml-outfiles: Estimated gene trees 
- astral-outfiles: Inferred species tree
=#

# todo for week Oct 7: add lines of code to create specific configuation file:
# - modify the "new" simphysim-conf configure file to have less, named "master": Done 
# - here: read the master: Done
# - add stuff to the string to add n_reps, n_genes & seed: Done 
# - write specific configuration file in simphyfolder: Done 

# A few notes: 
# 1) The bottom loop cannot be removed (see notes below)
# 2) Add the duploss (yes and no)

# set all configuration parameters here
duploss = "no"
ratevar = "yes"
n_reps = 2
n_genes = 10
# number of sites / gene is always 1000: hard-coded into seq-gen.sh

seed = 100 

# Set up seeds 
#=
seed_set= 1:5 # start with 5 seeds first but this could be changed to a larger number 
for i in 1:5 
  seed = seed_set[i]
  println("Running simulation with seed = $seed")
end
=# 

rootfolder = pwd()
outfolder = "output/DL$duploss-RV$ratevar"
mkpath(outfolder)
simphyfolder = joinpath(outfolder, "sim-phy-outfiles")
seqgenfolder = joinpath(outfolder, "seq-gen-outfiles")
raxmlfolder  = joinpath(outfolder, "raxml-outfiles")
astralfolder = joinpath(outfolder, "astral-outfiles")
mkdir(simphyfolder)
mkdir(seqgenfolder)
mkdir(raxmlfolder)
mkdir(astralfolder)

master_conf = joinpath(rootfolder,"simphy-configs/", "simphysim-conf-master")
conf_content = read(master_conf, String) # read the master config file into a string 

#Set up the parameters (# replications, # genes, and seeds)
parameters = """
# Parameters: DL:$duploss-RV:$ratevar, with replicate = $n_reps genes = $n_genes and seeds = $seed: 
-rs $n_reps  // Number of replicates
-rl f:$n_genes  // Number of loci (genes) per replicate - f means a fixed value 
-cs $seed  // seed
"""
combined_content = parameters * conf_content # combine parameters with master config
new_conf_file = joinpath(simphyfolder, "simphysim-conf-DL$duploss-RV$ratevar")
write(new_conf_file, combined_content) # write the combined config into SimPhy output folder

#-----------------------------------------------#       
#  simulate gene trees using SimPhy 
#-----------------------------------------------#

run(`$rootfolder/executables/simphy -i $new_conf_file -o $simphyfolder`)

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
#  simulate raxml and astral using seq-gen
#-----------------------------------------------#
# run raxml on each gene of each rep: use raxml.pl on each rep
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`perl ./scripts/raxml.pl --seqdir=$seqgenfolder/rep$simulation_rep --raxmldir=raxml-outfiles$simulation_rep --astraldir=astral-outfiles$simulation_rep`)
end

#removing this could cause issues with setting up raxmlfolder 
# The issue came from the line in raxml.pl: die ("raxmldir should be only 1 level up\n") if ($raxmldir =~ /\//); 
# Commenting out this line will cause further errors when setting up philp file directory 
# This needs to be solved 
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`mv raxml-outfiles$simulation_rep $raxmlfolder`)
  run(`mv astral-outfiles$simulation_rep $astralfolder`)
end


