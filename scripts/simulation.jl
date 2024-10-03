#=
this script does:

1.
2.
3.

example: run this below
 on darwin cluster use the `/nobackup` dir with tmux
 on franlin use `/nobackup2`
 from the root of the repository:

julia scripts/simulation.jl

output: it will create folder xxx with inside:
- 
- 
=#

# set all configuration parameters here
duploss = "no"
ratevar = "yes"
n_reps = 2
n_genes = 10
# number of sites / gene is always 1000: hard-coded into seq-gen.sh

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

# todo: add lines of code to create specific configuation file:
# - modify the "new" simphysim-conf configure file to have less, named "master"
# - here: read the master
# - add stuff to the string to add n_reps, n_genes & seed
# - write specific configuration file in simphyfolder

run(`executables/SimPhy -i $(simphyfolder)/simphysim-conf-new -o $(simphyfolder)/sim-phy-outfiles`)

# run seq-gen on each replicate and each gene
for simulation_rep in 1:n_reps
  rep_number_string = lpad(string(simulation_rep), ceil(Int, log10(n_reps+1)), '0')
  run(`mkdir ../seq-gen-outfiles/simphy$rep_number_string`)
  for gene_tree in 1:n_genes
    # account for sim phy naming convention
    genenum_string = lpad(string(gene_tree), ceil(Int, log10(n_genes+1)), '0')
    # seq-gen.sh knows about the directory structure & names
    run(`bash executable/seq-gen.sh $rep_number_string $tree_string $outfolder`)
  end
end

# run raxml on each gene of each rep: use raxml.pl on each rep
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`perl raxml.pl --seqdir=$seqgenfolder/simphy$simulation_rep --raxmldir=$raxmlfolder/raxml-outfiles$simulation_rep --astraldir=$astralfolder/astral-outfiles$simulation_rep`)
end
#= remove this below if all works well
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`mv raxml-outfiles$simulation_rep $raxmlfolder`)
  run(`mv astral-outfiles$simulation_rep $astralfolder`)
end
=#

