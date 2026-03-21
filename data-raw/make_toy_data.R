set.seed(123)

n <- 120
p <- 20

toy_marker <- rnorm(n)

linpred <- 0.3 * scale(toy_marker)[, 1]
prob <- 1 / (1 + exp(-(-0.7 + linpred)))
toy_outcome <- rbinom(n, 1, prob)

toy_omics <- matrix(rexp(n * p), nrow = n, ncol = p)
colnames(toy_omics) <- paste0("feature_", seq_len(p))

usethis::use_data(
  toy_marker,
  toy_outcome,
  toy_omics,
  overwrite = TRUE
)