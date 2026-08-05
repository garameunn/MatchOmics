# MatchOmics

**Propensity score matching with weighted GEE for heterogeneity-aware omics association analysis.**

## Overview

MatchOmics implements a two-round propensity score matching algorithm that
improves statistical power for omics association studies in heterogeneous
cohorts. The key innovation is a rescue round that recovers unmatched treated
subjects without sacrificing PS-caliper quality.

### Algorithm

1. Estimate propensity scores from heterogeneity indices (PCA-based for
   proteomics/metabolomics; diversity-based for metagenomics).
2. **Round 1**: nearest-neighbour matching without replacement under a
   PS caliper.
3. **Round 2 (rescue)**: unmatched treated subjects are re-matched against
   all controls — controls remain re-eligible.
4. Build matched clusters across both rounds (union-find).
5. Assign mean-normalised weights: `w = m / mean(m)`, ensuring `sum(w) = N`
   for GEE.
6. Test marker–outcome association via weighted GEE:
   `outcome ~ marker_class + ps [+ outcome_covariates]` (propensity score
   included by default via `adjust_ps = TRUE`, matching the manuscript).

A standard single-round without-replacement matching is also available via
`method = "standard"`.

> **Terminology note**: the manuscript refers to `method = "two_round"` as
> **with-replacement matching** (a control may be reused across rounds 1 and
> 2) and to `method = "standard"` as **without-replacement matching**
> (round 1 only). The mechanics are identical; only the name differs.

## Installation

```r
install.packages("remotes")
remotes::install_github("garameunn/MatchOmics")
```

## Quick start

```r
library(MatchOmics)

# Load example data
data(toy_omics)           # 50 x 200 proteomics matrix
data(toy_outcome)         # length-200 binary outcome
data(toy_heterogeneity)   # 200 x 2 PCA-based heterogeneity

# Compute heterogeneity from omics matrix (if not pre-computed)
het <- compute_heterogeneity(toy_omics, omics_type = "proteomics")

# Test one marker — two-round matching (proposed method)
# Outcome model defaults to marker_class + ps (adjust_ps = TRUE)
fit <- MatchOmics(
  marker        = toy_omics[1, ],
  outcome       = toy_outcome,
  heterogeneity = het,
  caliper       = 0.3,
  method        = "two_round",    # or "standard"
  corstr        = "independence"  # or "exchangeable"
)

print(fit)
summary(fit)
coef(fit)

# Add nuisance covariates (e.g. age, sex) to the outcome model, or drop
# the propensity-score adjustment:
# MatchOmics(..., outcome_covariates = data.frame(age = age_vec))
# MatchOmics(..., adjust_ps = FALSE)
```

### Multi-marker loop

```r
results <- lapply(seq_len(nrow(toy_omics)), function(i) {
  tryCatch(
    MatchOmics(toy_omics[i, ], toy_outcome, het),
    error = function(e) NULL
  )
})
```

## Diagnostics

```r
# SMD before / after matching
mdata <- fit$matched_data

pre_smd  <- compute_smd(
  data.frame(toy_heterogeneity, marker_class = ifelse(toy_omics[1,] >= median(toy_omics[1,]), "up","down")),
  covnames = c("het1","het2")
)
post_smd <- compute_smd(mdata, covnames = c("het1","het2"))

sdf <- smd_summary(pre_smd, post_smd)
love_plot(sdf)

# Effective sample size
fit$neff
```

## Simulation code

Full simulation scripts used to benchmark methods (S1/S2/S3 scenarios,
N = 200/500/1000) are in `inst/simulation/`. See
`inst/simulation/README.md` for instructions.

## Citation

Kim, N.-E. et al. (2026). *MatchOmics: a propensity score matching framework
for heterogeneity-robust multi-omics association analysis.* (manuscript submitted)

## License

MIT
