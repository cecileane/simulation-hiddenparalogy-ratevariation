
run(`mkdir ../sim-phy-outfiles`)
run(`mkdir ../seq-gen-outfiles`)

run(`../Sim-Phy-Source/SimPhy -i ../simphy-configs/simphysim-conf -o ../sim-phy-outfiles/sim_out`)

for simulation_rep in 1:2 #1000 on final
  repition_string = string(simulation_rep)
  run(`mkdir ../seq-gen-outfiles/simphy$repition_string`)
  for gene_tree in 1:10 # 1000 on final
    
    tree_string = string(gene_tree)
    #account for sim phy naming convention
    if length(tree_string) == 1
      tree_string = "000" * tree_string
    elseif length(tree_string) == 2
      tree_string = "00" * tree_string
    elseif length(tree_string) == 3
      tree_string = "0" * tree_string
    elseif length(tree_string) == 4
      tree_string = tree_string
    end

    run(`bash seq-gen.sh $repition_string $tree_string`)
  
  end
  #run IQ alignment for later
end
#run astral
#run snaq
#run PhylonetMPL
