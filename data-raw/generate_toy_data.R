## Run this script once from the package root to (re)generate toy data.
## Requires the package itself to be loaded (or source R/heterogeneity.R first).
##
##   source("data-raw/generate_toy_data.R")
##
## Output: data/toy_omics.rda, data/toy_outcome.rda, data/toy_heterogeneity.rda

set.seed(42)

n <- 200   # subjects
p <- 50    # proteins

## ── proteomics-like matrix (features x subjects) ──────────────────────────────
# Two latent batches to create genuine heterogeneity
batch <- rep(c(0, 1), each = n / 2)
toy_omics <- matrix(
  rnorm(p * n, mean = 10 + rep(batch, each = p) * 0.8, sd = 1.5),
  nrow = p, ncol = n
)
colnames(toy_omics) <- paste0("S", seq_len(n))
rownames(toy_omics) <- paste0("Protein", seq_len(p))

## ── heterogeneity: PCA ────────────────────────────────────────────────────────
pp <- prcomp(t(toy_omics), scale. = TRUE)
toy_heterogeneity <- data.frame(
  het1 = pp$x[, 1],
  het2 = pp$x[, 2],
  row.names = colnames(toy_omics)
)

## ── binary outcome ────────────────────────────────────────────────────────────
# Outcome depends on target marker (Protein1) + heterogeneity (confounding)
marker_true  <- toy_omics[1, ]
lp <- 0.4 * scale(marker_true)[, 1] +
      0.3 * scale(toy_heterogeneity$het1)[, 1] +
      rnorm(n, sd = 0.5)
prob         <- plogis(lp - quantile(lp, 0.70))   # ~30 % prevalence
toy_outcome  <- as.integer(rbinom(n, 1, prob))

## ── save ──────────────────────────────────────────────────────────────────────
save(toy_omics,         file = "data/toy_omics.rda",         compress = "xz")
save(toy_outcome,       file = "data/toy_outcome.rda",       compress = "xz")
save(toy_heterogeneity, file = "data/toy_heterogeneity.rda", compress = "xz")

message("Toy data saved to data/")
