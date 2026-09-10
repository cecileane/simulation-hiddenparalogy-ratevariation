# PhyNEST n_inds=2 design: one individual per species (settings later dropped)

Status: FINAL DECISION (2026-08-26): the n_inds=2 settings are dropped from
the PhyNEST analysis; only the 18 n_inds=1 settings are used. The design
below stays in the pipeline (fasta2phylip.py) and is kept here as the
record of what was considered.

## Current design

PhyNEST cannot map multiple individuals to one species
(https://github.com/sungsik-kong/PhyNEST.jl/issues/15). For n_inds >= 2
settings, the pipeline subsamples to one individual per species BEFORE
PhyNEST runs:

- `scripts/fasta2phylip.py` (called by `scripts/phynest_1rep.jl`) keeps
  only the first individual (`_0`) of each species and renames tips
  from e.g. "A_0" to "A".
- So PhyNEST always analyzes an 8-taxon species-level alignment, for
  all 36 settings. Verified in real output: networks in
  `output/DUP0.0-LOS0.0-RVL-N_ind2-SF0.5-genelen1000/` have exactly the
  8 tips A-H, and score_truetree is computed against the 8-taxon true
  tree.

## Basis for the decision

Email chain with Cecile (2026-08-17/18, re: USYB-2026-124):

- Bing proposed either (a) assume single lineage per species, or
  (b) a SNaQ-style analog: average quartet site-pattern counts over the
  2x2x2x2 individual choices per species quartet (not supported by
  PhyNEST; would require developing it ourselves).
- Cecile: "using a single individual would make sense" and finally
  "let's go with 1 individual, because using 2 individuals would be
  extra work and it's not supported by PhyNEST, so doesn't quite
  represent what a typical user would do." She agreed the site-pattern
  averaging idea "would make total sense" conceptually, but said not to
  do it.
- Bing's reading at the time: "Single individual it is, and that's what
  the pipeline already does."

Open question (why this is flagged): an alternative reading of the
emails is to keep all 16 sequences and let PhyNEST treat each
individual as its own tip ("each tip = separate lineage", 16-taxon
analysis). A short check-in email to Cecile was drafted to confirm
which she meant, before production runs.

## Which individual to keep, and why `_0`

- The two individuals are exchangeable: SimPhy samples both from the
  same population, and under dup/loss the probability that a given
  individual carries a hidden paralog at a gene is the same for `_0`
  and `_1`. Keeping `_0` has the same distribution as keeping `_1` or
  choosing at random.
- The choice must be pre-data (never based on which individual "looks
  cleaner"), or it would bias toward orthologs and understate
  confounding. A fixed label is pre-data and reproducible.
- Per-gene mixing of individuals was rejected: does not represent a
  typical user who sequenced one individual (Cecile's rationale).

## Monophyly

- Gene-tree non-monophyly of conspecific individuals (ILS, paralogy)
  is expected and handled where both individuals are used: ASTRAL via
  the mapping file (-a), SNaQ via CF averaging over individual quartets.
- After subsampling, PhyNEST has one tip per species, so within-species
  monophyly cannot arise in its input or output; PhyNEST's model (one
  lineage per species) is exactly satisfied.
- Paralogy in the retained `_0` sequences shows up as species-level
  topology error or spurious hybrid edges, which the postprocessing
  records (RF to true tree, alter1/2/3 diagnostics, hybrid
  recipient/donors).

## Known caveat for the paper

At n_inds=2 the three methods use the individuals differently: SNaQ
averages CFs over both individuals, find_graphs pools both into the
population, PhyNEST analyzes a one-individual subsample. One Methods
sentence should state this. The matched-null T* calibration is
unaffected: null and experimental settings at n_inds=2 go through the
same subsampling.

## If the design is ever switched to 16 tips (each individual a lineage)

Would require: fasta2phylip.py keeps all tips; worker needs an
n_inds-aware true tree (each species expanded into a 2-leaf cherry) and
outgroup "A_0"; postprocessing gains collapsing rules (incl.
non-monophyletic conspecific pairs and hybrids between them). Costs:
quartet count grows 70 -> 1820 (~26x per likelihood evaluation) plus a
much larger search space; T becomes scale-incomparable across n_inds
cells (matched-null calibration survives, pooling across n_inds does
not). Benefit: uses all data; paralogy breaking conspecific monophyly
becomes directly observable. No n_inds=2 production runs existed as of
2026-08-25, so switching would not discard compute.
