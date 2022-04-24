for simulation_rep 1:1

  #run simphy and generate files
  simphy_command = './SimPhy -i configuration_files/simulation -o outfiles/sim_out'
  run(simphy_command)

  for gene in 1:1 #start with 1 get to 1000

    ## PLAN ON USING THE L_TREES FILE FOR THE SIMULATION SO I ONLY HAVE TO READ ONE FILE
    #run seq-gen

    #run IQ alignment for later
  end
  #run astral
  #run snaq
  #run PhylonetMPL