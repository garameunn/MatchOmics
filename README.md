# MatchOmics

MatchOmics is an R package for heterogeneity-aware multi-omics association analysis using propensity score matching and weighted generalized estimating equations (GEE).

## Installation

Install from GitHub:


## Example
```r
install.packages("remotes")
remotes::install_github("garameunn/MatchOmics")

library(MatchOmics)

data("toy_marker")
data("toy_outcome")
data("toy_omics")

het <- prepare_heterogeneity(
  omics_type = "proteome",
  data = toy_omics
)

fit <- MatchOmics(
  marker = toy_marker,
  outcome = toy_outcome,
  heterogeneity = het,
  replace = TRUE,
  caliper = 0.2,
  corstr = "exchangeable",
  sandwich = TRUE
)

summary(fit)
```

## Description

MatchOmics performs heterogeneity-aware association testing for omics data by combining propensity score matching with weighted GEE.

The framework:
1. Dichotomizes the target marker using the median
2. Estimates propensity scores based on heterogeneity-related covariates
3. Performs nearest-neighbor matching
4. Constructs matching clusters
5. Conducts association testing using weighted GEE

## License

MIT
