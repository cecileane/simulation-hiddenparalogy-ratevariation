# parameters for the substitution model used by seq-gen

look at log files from gene trees, to use parameters that fit the reptile data.

in `reptiles/crawford/iqtree`:
```sh
grep --no-filename -o "G[0-9]{[^}]*}" IQ_*/loci.best_model.nex | grep -oE "[0-9]\.[0-9]+" > gammashape.txt
grep --no-filename -o "HKY{[^}]*}" IQ_*/loci.best_model.nex | grep -oE "[0-9]\.[0-9]+" > kappa.txt
grep --no-filename -o "F{[^}]*}" IQ_*/loci.best_model.nex | grep -oE "[0-9,.]+" > basefrequency.txt
```

next, summarize the long list of parameters (1/gene) in julia:

```julia
using Distributions, CSV, DataFrames
# shape alpha for the Gamma distribution of rates across sites
alpha = CSV.read("gammashape.txt", DataFrame; header=false)[!,1]
for family in [Normal, LogNormal, Gamma]
  d = fit(family, alpha)
  loglik = sum(logpdf.(d, alpha))
  println("loglik=$(round(loglik,digits=1)), $d")
end
d = fit(Gamma, alpha)
mean(d) # αθ = 0.3564457282430214. 3.267*0.109 = 0.356103
median(d) # 0.32081130713038486
d = truncated(Gamma(3.267, 0.109); lower=0.10)
rand(d)    # to simulate 1 alpha value
rand(d, 5) # to simulate 5 values at once
# transition-transversion ratio kappa for HKY model
kappa = CSV.read("kappa.txt", DataFrame; header=false)[!,1]
for family in [Normal, LogNormal, Gamma]
  d = fit(family, kappa)
  loglik = sum(logpdf.(d, kappa))
  println("loglik=$(round(loglik,digits=1)), $d")
end
d = fit(LogNormal, kappa)
mean(d) # 4.309019014364004
median(d) # 4.143552395123885
# base frequencies, for HKY model
basefreq = CSV.read("basefrequency.txt", DataFrame; header=false) # 4 columns
basefreq = transpose(Matrix(basefreq)) # 4 rows: "fits" wants 1 column / sample
d = fit(Dirichlet, basefreq)
Dirichlet{Float64, Vector{Float64}, Float64}(alpha=[66.58969803920938, 38.414804512094776, 38.61218853342359, 67.11544454548475])
print(mean(d)) # [0.31599213779173824, 0.1822921046057452, 0.18322876298838114, 0.31848699461413543]
```

## conclusions

shape alpha for the Gamma distribution of rates across sites:

    loglik=110.9, Normal{Float64}(μ=0.3564457282430214, σ=0.20167590731757154)
    loglik=177.2, LogNormal{Float64}(μ=-1.1923578305291007, σ=0.596517100863769)
    loglik=191.6, Gamma{Float64}(α=3.266971881169338, θ=0.10910584517043342)

simulate using alpha=0.356 for all genes (-a option in seq-gen),
else simulate alpha for each gene before sending it to seq-gen
using Gamma(α=3.267, θ=0.109)

For HKY's kappa: LogNormal is the best fit (see below).
simulate using transition/transversion ratio kappa=4.143 (-t in seq-gen),
else simulate kappa for each gene from LogNormal(μ=1.4215, σ=0.2798)

    loglik=-864.9, Normal{Float64}(μ=4.313136127167631, σ=1.2807857687020354)
    loglik=-812.8, LogNormal{Float64}(μ=1.4215534863637036, σ=0.2798456180449342)
    loglik=-823.0, Gamma{Float64}(α=12.6295564985486, θ=0.34151128962155564)

For base frequencies:
simulate with a fixed vector of: 0.316,0.182,0.183,0.319
else simulate from
Dirichlet(66.59, 38.41, 38.61, 67.12).
