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

#= Previously, the duploss and ratevar are fixed 
duploss = "no"
ratevar = "yes"
n_reps = 2
n_genes = 10 # number of sites / gene is always 1000: hard-coded into seq-gen.sh
=# 

#= Notes: 
1) The very last part, which mv the output of raxml and astral to a folder, cannot be changed 
the issue is called by line 75 in raxml.pl: die ("raxmldir should be only 1 level up\n") if ($raxmldir =~ /\//); If comment out this line, then error occurs when setting up the phylipdir 
This could be resolved later or remain as it is 
2) Parameters for rate variations 
Based on README.md, it seems that only gene variation is examined. 
How about lineage variation or gene x lienage variations? 
3) I changed the last part of raxml.pl so that the final astral has one output -- need to double chec: Why do the original script generates boostrap trees? 
This could be changed to a newer version of unweighted astral later 
4) After changing the tips to letters, it seems that Seq-Gen
=# 

duploss = ARGS[1] # Yes or No
ratevar = ARGS[2] # "No", "Ge", "Li" or "GL” 
# In meeting, rate var is "yes" or "no" 
n_reps = parse(Int, ARGS[3])
n_genes = parse(Int, ARGS[4]) 
seed = parse(Int, ARGS[5])

# Print out the arguments 
println("Running simulation.pl for:\n")
println("duploss = $duploss, ratevar = $ratevar, n_reps = $n_reps, n_genes = $n_genes, seed = $seed")

# set all configuration parameters here 
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

# To simulate hidden paralogy, need to adjust duplication rate 
if duploss == "Yes" # Here, not sure about the duplication rate, so choose a very arbitrary number 
  parameters *= "-lb f:0.0001 // Duplication rate\n" # This should be changed 
  # Think: should I use -lb loss rate as well? 
  # What would be the dup rate to use and loss rate to use? 
end 

# To simulate substitution rate variation
if ratevar == "Ge" # gene-family-speciic rate heterogenity 
  parameters *= "-hl ln:-0.19,0.6164414002968976 //log-normal distribution of gene rates"
elseif ratevar == "Ln" # To simulate variation across lineages (ratevar = "Ln") 
  parameters *= "-hs f:0.01" # An arbitrary number for now -- need to change 
elseif ratevar == "GL"  # To simulate gene-by-lineage variation  
  parameters *= "-hh f:0.01" # An arbitrary number for now -- need to change 
end # else ratevar == "No" -- no additional para to be specified 


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


