
run(`mkdir ../sim-phy-outfiles`)
run(`mkdir ../seq-gen-outfiles`)

run(`../Sim-Phy-Source/SimPhy -i ../simphy-configs/simphysim-conf -o ../sim-phy-outfiles/sim_out`)

for simulation_rep in 1:2
  repition_string = string(simulation_rep)
  run(`mkdir ../seq-gen-outfiles/simphy$repition_string`)
  for gene_tree in 1:10
    
    tree_string = string(gene_tree)
    if length(tree_string) == 1
      tree_string = "0" * tree_string
    elseif length(tree_string) == 2
      tree_string = tree_string
    end

    run(`bash seq-gen.sh $repition_string $tree_string`)
  
  end
  #run IQ alignment for later
end
#run astral
#run snaq
#run PhylonetMPL
