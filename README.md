`````markdown

# MatchOmics

MatchOmics is an R package for heterogeneity-aware multi-omics association analysis using propensity score matching and weighted generalized estimating equations (GEE).

## Installation

Install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("garameunn/MatchOmics")
```

## Main features

- propensity score matching for omics association analysis  
- matching with replacement as the default option  
- optional caliper restriction  
- weighted GEE for valid inference after matching  
- omics-specific heterogeneity handling  

## Example

```r
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

fit
summary(fit)
coef(fit)
```

## Method overview

MatchOmics performs the following steps:

1. dichotomizes the target marker using the median  
2. estimates propensity scores based on heterogeneity-related covariates  
3. performs nearest-neighbor matching  
4. constructs matching clusters  
5. conducts association testing using weighted GEE  

## Output

`MatchOmics()` returns an object containing:

- fitted propensity score model  
- matching object  
- matched dataset  
- weighted GEE fit  
- coefficient table  

## Availability

The source code is available at:

- GitHub repository: https://github.com/garameunn/MatchOmics  

## License

MIT

```