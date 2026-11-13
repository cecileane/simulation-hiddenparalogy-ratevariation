# Third-party code

Code from other projects, not written for this repository.

- `vcf2eigenstrat_biallelic.py` and `vcf2eigenstrat_multiallelic.py`:
  adapted from [mathii/gdc](https://github.com/mathii/gdc) by Iain
  Mathieson (Apache License 2.0, see LICENSE-APACHE-2.0.txt). Ported to
  Python 3 and adjusted for this project. The pipeline uses the biallelic
  version; the multiallelic one, which handles missing alleles, is kept
  for reference.
- `interop_admixtools.jl`: from [PhyloUtilities](https://github.com/JuliaPhylo/PhyloUtilities) of JuliaPhylo (MIT license,
  see mit-license.txt). Converts admixtools graphs to PhyloNetworks
  networks, used by findgraphs_postprocess.jl.
