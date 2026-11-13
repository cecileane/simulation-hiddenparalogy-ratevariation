# Executables

The external programs of the pipeline, as symlinks to where they are
installed (`~/private/software/` on our servers), so that the scripts can
call them as `executables/<name>`. The symlinks are not tracked by git.
The commands we used to install everything are in
[software_installation.sh](software_installation.sh).

- `simphy`: [SimPhy](https://github.com/adamallo/SimPhy) 1.0.2, gene trees
  under the coalescent with duplication and loss
- `seq-gen`: [Seq-Gen](https://github.com/rambaut/Seq-Gen) 1.3.5, sequences
  along the gene trees
- `iqtree2`: [IQ-TREE 2](https://github.com/iqtree/iqtree2) 2.4.0, gene tree
  estimation
- `astral`: [ASTER](https://github.com/chaoszhang/ASTER) 1.24.4.8, species
  tree estimation (ASTRAL-IV). astral-pro3 and wastral from the same
  package were only used for testing
- `snp-sites`: [snp-sites](https://github.com/sanger-pathogens/snp-sites)
  2.5.1, SNPs from the concatenated alignment, for find_graphs

To add one: install it somewhere, then from the repo root
`ln -s <path to the binary> executables/<name>` and check with
`executables/<name> -h`.

Also needed: Julia 1.11 (juliaup) with the packages of Project.toml
(`julia --project=. -e 'using Pkg; Pkg.instantiate()'` from the repo root).
PhyNEST needs an older PhyloNetworks than the rest, so it has its own
environment in `envs/phynest/` (Project.toml and Manifest.toml); the
PhyNEST scripts activate and instantiate it themselves, nothing to install
by hand. Python 3 with biopython (concatenate_seq.py, fasta2phylip.py), and R with
admixtools 2.0.8 (igraph 1.6, see scripts/findgraphs_1rep.R), optparse,
Rcpp, ape and the packages of the figure scripts (ggplot2, cowplot,
scales, dplyr, tidyr, readr, stringr, ggnewscale).
All of this was run on Linux x86-64.
