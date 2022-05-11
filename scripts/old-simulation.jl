#on darwin cluster use the /nobackup dir with tmux
run(`mkdir ../sim-phy-outfiles`)
run(`mkdir ../seq-gen-outfiles`)
run(`mkdir ../raxml-outfiles`)
run(`mkdir ../astral-outfiles`)


#do I need a for loop here or is the simphy config sufficent?
run(`../executables/SimPhy -i ../simphy-configs/simphysim-conf-new -o ../sim-phy-outfiles/sim_out`)

for simulation_rep in 1:2 #1000 on final
  repition_string = string(simulation_rep)
  run(`mkdir ../seq-gen-outfiles/simphy$repition_string`)
  for gene_tree in 1:2 # 1000 on final
    
    tree_string = string(gene_tree)
    #account for sim phy naming convention
    tree_string = lpad(tree_string, 2, '0')

    #runs seq gen for each tree in this rep of the simulation
    #change to a pipeline
    run(`bash seq-gen.sh $repition_string $tree_string`)
  
  end
end