# simphy-configs

`simphysim-conf-master` is the SimPhy configuration that simulation.jl
starts from. The script reads it once, adds the parameters of the setting
(population size, duplication and loss rates, rate variation, species
tree) in front of it, then for every replicate and iteration adds the seed
and the number of loci and writes the result to
`simphysim-conf-Int<iter>-Rep<rep>`, which is passed to SimPhy with `-i`.
