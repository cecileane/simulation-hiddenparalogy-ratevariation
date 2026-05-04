Summary: SimPhy and Seq-Gen Interpretation for Lineage-Specific Simulations 

# Simphy effective population size (Ne): 
SimPhy defines the input Ne as the **haploid effective population size** — the number of gene copies in the population (see this link: https://groups.google.com/g/simphy/c/2QE0QeG2D_s). When we specify Ne in SimPhy, it is treated as-is (no internal multiplication by 2). This means for a diploid population of 1000 individuals, we should input `Ne = 2000` (since 1000 diploids = 2000 haploid gene copies).

# Simphy Substitution Rate Multiplier (*) 
The substitution rate multiplier in SimPhy is dimensionless. It scales the baseline substitution rate. It is expressed as $multiplier = branch specific substitution rate per site per generation / baseline substitution rate$

# Baseline substitution rate 
The baseline substitution rate is set using the `-su` parameter (or `-SU`) in the SimPhy configuration or command-line. The value represents **substitutions per site per generation**. SimPhy uses this to compute branch lengths on gene trees in expected substitutions per site, before passing to the sequence simulator.


# SimPhy outputs
Change gene tree with branch length in sub rate 







 



