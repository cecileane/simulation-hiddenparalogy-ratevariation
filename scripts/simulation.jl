n_reps=2
n_genes=10
#on darwin cluster use the /nobackup dir with tmux
mkdir("../sim-phy-outfiles")
mkdir("../seq-gen-outfiles")
mkdir("../raxml-outfiles")
mkdir("../astral-outfiles")


#do I need a for loop here or is the simphy config sufficent?
run(`../executables/SimPhy -i ../simphy-configs/simphysim-conf-new -o ../sim-phy-outfiles/sim_out`)

for simulation_rep in 1:n_reps #1000 on final number of replicates is hardcoded in config file
  repition_string = string(simulation_rep)
  repition_string = lpad(repition_string, ceil(Int, log10(n_reps+1)), '0')
  run(`mkdir ../seq-gen-outfiles/simphy$repition_string`)
  for gene_tree in 1:n_genes # 1000 on final
    
    tree_string = string(gene_tree)
    #account for sim phy naming convention
    tree_string = lpad(tree_string, ceil(Int, log10(n_genes+1)), '0')

    #runs seq gen for each tree in this rep of the simulation
    #change to a pipeline
    run(`bash seq-gen.sh $repition_string $tree_string`)
  
  end
end
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`perl raxml.pl --seqdir=../seq-gen-outfiles/simphy$simulation_rep --raxmldir=raxml-outfiles$simulation_rep --astraldir=astral-outfiles$simulation_rep`)
end
for simulation_rep in 1:n_reps
  simulation_rep = lpad(simulation_rep, ceil(Int, log10(n_reps+1)), '0')
  run(`mv raxml-outfiles$simulation_rep ../raxml-outfiles`)
  run(`mv astral-outfiles$simulation_rep ../astral-outfiles`)
end

