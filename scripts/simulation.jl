for simulation_rep in 1:2

  #run simphy and generate files
  string_rep = string(simulation_rep)
  run(`../Sim-Phy-Source/SimPhy -i ../simphy-configs/simphysim-conf -o ../outfiles/sim_out$string_rep`)

  for gene in 1:1 #start with 1 get to 1000
    

    #run seq-gen

    #run IQ alignment for later
  end
  #run astral
  #run snaq
  #run PhylonetMPL
end