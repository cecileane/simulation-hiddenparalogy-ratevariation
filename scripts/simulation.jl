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

## Now, start to add in ArgParse package 

using ArgParse 

function parse_commandline()
  s = ArgParseSettings() 
  @add_arg_table s begin 
    "--duploss"
      help = "Gene duplication rate: specify the gene duplication rate, if 0, then no duplication"
      arg_type = Number
      required = true
    "--ratevar" 
      help = "'N': No rate variation';\n 'G: Gene specific rate variation';\n 'L': Lineage specific rate variation;\n 'GL' or 'G*L': genexlineage rate variation"
      arg_type = String
      choices =  ['N', 'G', 'L', 'G*L']
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
  end 
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  println("Parsed args:")
  for (arg,val) in parsed_args
      println("  $arg  =>  $val")
  end

  duploss = parsed_args["duploss"] # A number indicates rate of gene duplication rates (0 = no dup)
  ratevar = parsed_args["ratevar"] # "N", "G", "L" or "GL" or "G*L" to include genexlineage rate variation 
  n_reps = parsed_args["n_reps"] # number of replicates 
  n_genes = parsed_args["n_genes"] # number of genes 
  seed_simphy = parsed_args["seed_simphy"] # seeds for SimPhy 

  println("duploss = $duploss, ratevar = $ratevar, n_reps = $n_reps, n_genes = $n_genes, seed = $seed_simphy")
end

main() 

# set all configuration parameters here 
rootfolder = pwd()
paramname_root = "DL$duploss-RV$ratevar"
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

master_conf = joinpath(rootfolder,"simphy-configs/", "simphysim-conf-master")
conf_content = read(master_conf, String) # read the master config file into a string 

#Set up the parameters (# replications, # genes, and seeds)
parameters = """
# Parameters: DL:$duploss-RV:$ratevar, with replicate = $n_reps genes = $n_genes and seeds = $seed_simphy: 
-rs $n_reps  // Number of replicates
-rl f:$n_genes  // Number of loci (genes) per replicate - f means a fixed value 
-cs $seed_simphy  // seed
"""

# To simulate hidden paralogy, need to adjust duplication rate 
if duploss == "Y" # Here, not sure about the duplication rate, so choose a very arbitrary number 
  parameters *= "-lb f:0.0001 // Duplication rate\n" # This should be changed 
  # Think: should I use -lb loss rate as well? 
  # What would be the dup rate to use and loss rate to use? 
end 

# To simulate substitution rate variation
if occursin("G", ratevar) # gene-family-speciic rate heterogenity : "G" or "GL"
  parameters *= "-hl ln:-0.19,0.6164414002968976 //log-normal distribution of gene rates"
end
if occursin("L", ratevar) # To simulate variation across lineages (ratevar = "L" or "GL" or "G*L") 
  parameters *= "-hs f:0.01" # An arbitrary number for now -- need to change
end
if occursin("G*L", ratevar)  # To simulate gene-by-lineage variation  
  parameters *= "-hh f:0.01" # An arbitrary number for now -- need to change 
end # else ratevar == "No" -- no additional para to be specified 

# change new_conf_file 
combined_content = parameters * conf_content # combine parameters with master config
new_conf_file = joinpath(simphyfolder, "simphysim-conf-$paramname_root")
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


