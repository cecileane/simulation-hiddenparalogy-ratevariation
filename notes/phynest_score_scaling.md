# What the model-selection scores depend on (find_graphs, SNaQ, PhyNEST)

Purpose: state, for the manuscript, what each method's score is a function
of, so that the PhyNEST threshold T* can be reported with the quantities it
depends on (number of taxa/quartets, alignment length) and readers know when
it must be recalibrated. Sources: the PhyNEST paper (Kong, Swofford &
Kubatko 2025, Syst Biol 74:53) and source code (PhyNEST v0.1.12,
`Theta.jl`, `Quartets.jl`, `readPhylip.jl`), SNaQ source (`pseudolik.jl`),
the ADMIXTOOLS 2 documentation (qpgraph "graphs" article).

## 1. find_graphs / qpgraph (ADMIXTOOLS 2)

Score:  S = 1/2 (g - f)' Q^-1 (g - f)
- g = observed f-statistics (f3 basis; 406 distinct f2/f3/f4 contrasts for
  8 populations as counted in the paper), f = values fitted by the graph,
  Q = block-jackknife covariance of g.
- Depends on the NUMBER OF POPULATIONS (number of f-statistics grows
  polynomially with n) and on the AMOUNT OF DATA: Q shrinks as the number of
  SNPs/blocks grows, so for a fixed misfit the score grows with the number
  of SNPs. Equal to minus the log-likelihood up to a constant (Maier et al.
  2023), hence the "score within 2 units" retention rule (Delta AIC = 4).
- Worst residual WR = max |(g_i - f_i)/SE_i| over all residuals: a
  STANDARDIZED quantity, so its null distribution does not grow with the
  number of SNPs, but it grows with the number of residuals (max over more
  standard-normal-like values), i.e. with the number of taxa. This is why
  WR <= 3.0 / 3.7 is portable across data sizes but, as the Discussion
  says, must depend on the number of taxa and the true tree.

## 2. SNaQ (PhyloNetworks / SNaQ.jl)

Score (`logPseudoLik`):
  -sum over 4-taxon sets q  sum over the 3 quartet topologies i of
   100 * obsCF_qi * log(expCF_qi / obsCF_qi)
- The CFs are PROPORTIONS of genes, so the score does not grow with the
  number of genes (only its variance shrinks); it is 100 x the summed
  KL-type divergence between observed and expected CFs over the C(n,4)
  quartets (70 for 8 taxa) and equals 0 at a perfect fit. Typical values in
  our runs: 3-11 for h=0, i.e. 0.05-0.15 per quartet.
- Not used for model selection in the paper (the quartet goodness-of-fit
  test with 1000 simulated data sets is used instead), because comparing
  pseudolikelihoods across h has no calibrated reference distribution.

## 3. PhyNEST (composite likelihood on site patterns)

Score = negative log composite likelihood
  score = sum over quartets q (all C(n,4)) sum over the 15 site-pattern
          classes k of  n_qk * (-log p_qk(tau, theta, gamma))
- n_qk = number of alignment sites showing pattern class k in quartet q.
  `readPhylip` counts EVERY site (constant sites are pattern class 1), and
  every site is counted once in every quartet, so the data enter as
  L x C(n,4) "observations" (10^6 sites x 70 quartets = 7 x 10^7 for us);
  p_qk = probability of the class under JC69 + strict clock + constant
  theta on the quartet network (Chifman & Kubatko 2015), averaged over the
  displayed trees with weights gamma.
- Magnitude: our scores are about 5 x 10^7 = 0.65-0.8 nats per site per
  quartet (the entropy of the site-pattern distribution), i.e. the score
  is EXTENSIVE in both alignment length and number of quartets.
- T = score(h=0) - score(h=1) = L x sum_q sum_k f_qk log(p_qk^h1/p_qk^h0),
  with f_qk the observed pattern frequencies. Therefore:
  (a) for a fixed degree of misfit, T grows LINEARLY with the alignment
      length L;
  (b) T sums over the quartets whose site-pattern probabilities change
      when the reticulation is added (for the (G,H)-vs-rest edge, the 55 of
      70 quartets containing G or H), so it grows with the number of taxa;
  (c) T is a maximum over search runs, so it also depends on the number of
      runs (median null T doubled from 10 to 100 runs).
- Why T is so large (about 2 x 10^4, or 0.04% of the score): under a
  correctly specified full likelihood, twice the log-likelihood ratio of
  two nested models would be O(1) (chi-square with a few degrees of
  freedom). The observed null T is four orders of magnitude larger because
  (i) the composite likelihood treats 7 x 10^7 dependent site-quartet
  observations as independent, and (ii) the data are HKY+Gamma with
  per-gene parameters while the model is JC69 without rate variation, so
  the extra edge absorbs systematic misfit whose contribution scales with
  L rather than with sqrt(L). The per-replicate consistency of the gain
  (the (G,H)-vs-rest reticulation always gains about 20,000 units) is the
  signature of misfit, not of sampling noise.

## 4. Reporting T* in the paper

T* was calibrated with n = 8 taxa (including the outgroup, 70 quartets),
1,000 loci x 1,000 bp = 10^6 sites, 100 runs per model, data simulated
under HKY+Gamma and analysed under PhyNEST's JC69:

| ILS level | T*     | per site (T*/L) | per site per quartet | per quartet |
|-----------|--------|-----------------|----------------------|-------------|
| low       | 23,501 | 0.0235          | 3.4 x 10^-4          | 336         |
| high      | 18,489 | 0.0185          | 2.6 x 10^-4          | 264         |

Recommendation for the text: give T* together with the design it applies
to, and state that, because the composite likelihood counts every site in
every quartet, the threshold scales with alignment length and with the
number of quartets affected, and depends on the substitution model and on
the true tree; it must therefore be recalibrated for other data sets, with
the same simulation-based procedure recommended for WR (simulate under the
best tree with the study's own design and take the 95th percentile of T).
Rough empirical check from the 100-gene pilot (10^5 sites, 10 runs, high
ILS null): 95th percentile of T about 880, versus about 4,800 for the
best-of-10-runs T at 10^6 sites, i.e. roughly 5-fold for a 10-fold longer
alignment (between sqrt(L) and L; single setting, indicative only).

## 5. Verification: is the "size" of each score the same in every setting?

Checked on the actual outputs (2026-09-05):

| method | unit summed over | count for our 8 taxa | same in every setting? |
|---|---|---|---|
| PhyNEST | quartets, C(n,4) | 70 (all 1,800 H0 logs: "Number of sequences: 8") | yes: depends only on the number of taxa |
| SNaQ | 4-taxon sets, C(n,4) | 70 rows in every `CF_results.csv` | yes |
| find_graphs | f-statistics | 406 distinct contrasts (28 f2 + 168 f3 + 210 f4) | yes: depends only on the number of populations |

What is NOT identical across settings is the alignment length L that PhyNEST
multiplies into every quartet. From the PhyNEST logs ("Sequence length"):
- dup/loss = 0 and 3e-4: 1,000 genes (10^6 sites) in every replicate
  (two replicates with 999);
- dup/loss = 4e-4: 940-976 genes (mean 957-958) in all 6 settings, i.e.
  about 4% fewer sites, because fewer genes survive the paralogy filter.
Since T is linear in L, T in the 4e-4 settings is expected to be ~4% lower
for the same misfit; this is within replicate noise but should be stated
(and the Methods sentence "at least min(1000, K) valid trees were passed to
Seq-Gen" describes the gene trees, not the number of genes that reach the
concatenated alignment).

## 6. PhyNEST's "equal rate" assumptions and which settings violate them

PhyNEST's site-pattern probabilities are the closed-form JC69 formulas of
Chifman & Kubatko (2015, J Theor Biol 374:35), which require (Kong,
Swofford & Kubatko 2025, Syst Biol 74:53, Methods and Discussion):
"coalescent independent sites that evolve according to the Jukes-Cantor
(JC69) substitution model", "a strict molecular clock on which [the] site
pattern probabilities depend", and "theta (= 4 Ne mu) is constant
throughout the tree". Four separate equal-rate assumptions follow:

1. Equal substitution rates ACROSS LINEAGES (strict molecular clock).
   Our data: violated only under lineage-specific rate variation, where
   the SimPhy species tree carries per-branch multipliers m_i = r_i / mean r
   (scripts/simulation.jl): tips A 0.52, B 0.22, C 0.19, D 1.21, E 1.02,
   F 0.34, G 4.09, H 9.20; internal branches 0.35-3.56, i.e. up to ~48-fold
   rate differences between lineages (H vs C). Pang, Liu & Zhang (2025, Mol
   Biol Evol 42:msaf216) show that site-pattern methods (D-statistic, HyDe)
   built on the same clock assumption reach false-positive rates of 35% and
   100% for rate differences of only 17% and 33%; a relaxed-clock version of
   the site-pattern probabilities exists (Richards & Kubatko 2022, J Theor
   Biol, PMID 35278472) but is not used by PhyNEST. This is the assumption
   whose violation produces the collapse under lineage rate variation
   (F misplaced in 58-80% of h=0 trees, gamma_minor ~0.45, T x7, 70-98%
   false detections).
2. Equal rates ACROSS SITES (no Gamma) and equal exchange rates / base
   frequencies (JC69). Our data: violated in EVERY setting: HKY with
   per-gene kappa ~ LogNormal(1.42, 0.28), Dirichlet base frequencies, and
   Gamma shape alpha ~ Gamma(3.267, 0.109), i.e. mean alpha = 0.36 (strong
   among-site rate heterogeneity). The PhyNEST paper tested robustness to
   HKY (kappa = 3, unequal base frequencies) but not to among-site rate
   variation; it notes that under JC69 "unequal base frequencies and
   unequal rates of substitution between different pairs of nucleotides
   cause underestimation of the amount of divergence", and lists "more
   complex substitution models than JC69" as future work. This chronic
   misfit is the likely source of the behaviour seen even under no rate
   variation: the (G,H) split estimated at age 0, the stereotyped
   (G,H)-vs-rest reticulation with gamma ~0.17 gaining a near-constant
   ~20,000 units, T four orders of magnitude above what a correctly
   specified likelihood-ratio would give, and the true tree displayed in
   only ~50% of h=1 networks.
3. Equal rates ACROSS GENES. Our data: violated under gene-specific rate
   variation (SimPhy -hl ln:-0.19,0.62: log-normal gene multipliers with
   mean 1, about 0.3x-3.3x). For PhyNEST, which concatenates all genes and
   treats sites as exchangeable, this is a block form of among-site rate
   heterogeneity; its effect was small (gene-rate settings behave like the
   no-rate-variation settings).
4. Equal population size (theta) ACROSS BRANCHES. Our data: NOT violated:
   SimPhy received one fixed Ne (-sp f:Ne x SF), so this assumption cannot
   explain any of the results.

In short: the advisor's "equal rate" remark covers (1) and (2). The
lineage-rate collapse is a strict-clock violation of a magnitude far beyond
what the site-pattern literature already flags as fatal; the poor null
behaviour is JC69-without-Gamma misfit on HKY+Gamma data, amplified by the
composite likelihood treating 7 x 10^7 site-quartet observations as
independent.

## 7. Exact wording from the PhyNEST paper (Kong, Swofford & Kubatko 2025)

From the Methods, section "Site Pattern Probability Distributions":
- "Define a site pattern arising from tree T with n tips as an assignment
  of nucleotide states s1 s2 ... sn where sy in {A,C,G,T}, y = 1,2,...,n,
  to the tips of the tree."
- "Under the JC69 model, there are 15 distinct site pattern probabilities
  on a quartet tree under the MSC" (their Eq. 1).
- "We assume that sites are unlinked and that each site is an independent
  and identically distributed observation from the species tree under the
  MSC ... a data type referred to as coalescent independent sites (CIS)."

From the section "Computing the Composite Likelihood of a Network":
- "There are R = (n choose 4) quartets in a network containing n tips."
- "The site pattern probabilities for Q^N are the weighted sum of the site
  pattern probabilities of each D_t for that subset of 4 taxa:
  p^N = p_{s1 s2 s3 s4 | (Q^N, tau, gamma, theta)}
      = sum_{t=1}^{2^h} Gamma_t p_{s1 s2 s3 s4 | (D_t, tau, theta)}"
  (D_t = displayed trees, Gamma_t = product of inheritance probabilities).
- "For each site m in a data matrix containing M sites, and for each
  quartet q, q = 1,2,...,R, we define I_q(m) to be a random vector of
  length 15 ..." (the indicator of which of the 15 patterns site m shows
  in quartet q; W_q = sum over sites = the 15 pattern counts).
- "the likelihood l_q for a quartet is expressed as a function of tau,
  gamma, and theta as: l_q(tau, gamma, theta | W_q) proportional to
  prod_{j=1}^{15} ((p^N)_j)^{(W_q)_j}"
- "the function L(tau, gamma, theta | N, d) = prod_{q=1}^{R} l_q(tau,
  gamma, theta | W_q), is the composite likelihood of N".

Read together: the alignment (M sites) is reduced, for each of the
R = C(n,4) quartets, to a vector W_q of 15 site-pattern counts summed over
ALL M sites; each quartet contributes a multinomial-type likelihood with
the counts as exponents; the network's composite likelihood is the
product over the R quartets (i.e. quartets are multiplied as if
independent, and every site is used R times). PhyNEST reports
-log L; T = -log L(h=0) + log L(h=1). This matches the code
(`Theta.jl::get_negative_log_clikelihood`: sum over quartets of
-count x log p over 15 classes; `readPhylip.jl::sitePatternCounts` counts
every site). For our data: M = 10^6 (0.94-1.0 x 10^6 at dup/loss 4e-4),
R = 70.

## 8. Should T be rescaled when alignment lengths differ (as WR is)?

What WR does: each f-statistic residual is divided by its own block-
jackknife standard error, so WR is a per-statistic z-score. The PhyNEST
analog would be T / SE(T) with SE from a block bootstrap over loci; PhyNEST
does not provide it and re-fitting 100 x 100 runs per bootstrap replicate
is not feasible. Dividing T by the standard deviation of T across
replicates would only rescale the threshold, it does not remove the
dependence on alignment length.

What does address length: T is linear in M (section 3), so T / M ("nats
per site") is the natural scale-free version. Checked on the actual data
(per-replicate M from the PhyNEST logs; pooled rule per ILS level):

| | raw T | T per site |
|---|---|---|
| T* low ILS (SF 0.5) | 23,501 | 0.02393 |
| T* high ILS (SF 1.0) | 18,489 | 0.01897 |
| replicates whose call changes | - | 22 of 1,800 |
| dup/loss 4e-4 vs 0 (no/gene rates, low ILS): median T ratio | 0.971 | 1.013 (M ratio 0.958) |

So per-site scaling behaves as expected (the 4% shorter alignments at
dup/loss 4e-4 give ~3% smaller raw T, and per-site T removes it), but the
decision changes for only 22/1800 replicates (at most +5 calls in one
setting) and none of the conclusions. Recommendation for the paper: keep
the calibration as implemented (raw T, pooled per ILS level) but REPORT
T* also per site (0.024 and 0.019 nats per site) and state that, because
T is a sum over sites and quartets, the per-site value is the one to carry
to other data sets, together with the number of quartets. If the advisor
prefers the fully scale-free version, switching the statistic to T / M is
a two-line change (record M per replicate in phynest_postprocess.jl from
the log line "Sequence length", divide in summary_phynest.jl) and shifts
the calls by the 22 replicates above.

## 9. Candidate manuscript sentences

Status (2026-09-08): the per-site statistic was adopted after discussion
with Cecile. Each replicate's T is divided by its own alignment length M
(read from the PhyNEST log, "Sequence length"), T* is the 95th percentile
of T/M pooled as before (per ILS level), and the manuscript states that
T* is for 8 taxa (R = 70 quartets) and 100 search runs; no division by R
(R is constant across our replicates, and T is not proportional to R
because a reticulation affects only the quartets that straddle it).
Implemented in phynest_postprocess.jl (columns n_sites, T_per_site),
summary_phynest.jl (T* on T/M, T_star_raw kept for reference) and the
R plotting functions (t_column = "T_per_site" by default).

Methods (PhyNEST paragraph, after the search description):
"PhyNEST reduces the concatenated alignment of M sites to, for each of
the R = C(8,4) = 70 quartets of taxa, the counts of the 15 site patterns
that are distinguishable under JC69, and maximizes the composite
likelihood, i.e. the product over quartets of the multinomial likelihood
of these counts (Kong et al. 2025). Because every site enters every
quartet, the score is a sum over M x R site-quartet terms, so T scales
with alignment length and with the number of quartets whose site-pattern
probabilities are affected by the reticulation."

Methods (threshold):
"T* was 23,501 (low ILS) and 18,489 (high ILS) for our design of 8 taxa
(70 quartets) and 10^6 sites, i.e. 0.024 and 0.019 per site; alignments
at the highest duplication and loss rate were 4% shorter (940-976 loci
retained), which shifts T by the same proportion and changed the
selection in 22 of 1,800 replicates when T was expressed per site."

Discussion (portability, parallel to the WR paragraph):
"As for the worst-residual threshold, no single value of T* is broadly
appropriate: the composite likelihood counts every site in every quartet,
so T grows with alignment length and with the number of taxa, and its
null distribution depends on the substitution model and on the true tree.
T* should be recalibrated for each study by simulating tree-like data
under the best tree with the study's own design."
