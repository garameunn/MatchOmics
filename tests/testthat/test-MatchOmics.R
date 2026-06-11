test_that("MatchOmics runs with two_round method", {
  data(toy_omics, envir = environment())
  data(toy_outcome, envir = environment())
  data(toy_heterogeneity, envir = environment())

  fit <- MatchOmics(
    marker        = toy_omics[1, ],
    outcome       = toy_outcome,
    heterogeneity = toy_heterogeneity,
    caliper       = 0.3,
    method        = "two_round",
    corstr        = "independence"
  )

  expect_s3_class(fit, "MatchOmics")
  expect_true(fit$n_matched > 0)
  expect_true(fit$n_matched <= fit$n_input)
  expect_true(is.finite(fit$neff))
  expect_true(nrow(fit$coef_table) >= 1)
})

test_that("MatchOmics runs with standard method", {
  data(toy_omics, envir = environment())
  data(toy_outcome, envir = environment())
  data(toy_heterogeneity, envir = environment())

  fit <- MatchOmics(
    marker        = toy_omics[2, ],
    outcome       = toy_outcome,
    heterogeneity = toy_heterogeneity,
    caliper       = 0.3,
    method        = "standard",
    corstr        = "independence"
  )

  expect_s3_class(fit, "MatchOmics")
  expect_equal(fit$method, "standard")
})

test_that("compute_heterogeneity returns correct dimensions", {
  data(toy_omics, envir = environment())
  het <- compute_heterogeneity(toy_omics, omics_type = "proteomics")
  expect_equal(nrow(het), ncol(toy_omics))
  expect_equal(ncol(het), 2)
  expect_named(het, c("het1", "het2"))
})

test_that("compute_neff returns value between 0 and N", {
  w <- rep(1, 100)
  expect_equal(compute_neff(w), 100)
  w2 <- c(rep(2, 50), rep(0.1, 50))
  neff <- compute_neff(w2 / mean(w2))
  expect_true(neff > 0 && neff <= 100)
})

test_that("input validation errors are informative", {
  data(toy_omics, envir = environment())
  data(toy_outcome, envir = environment())
  data(toy_heterogeneity, envir = environment())

  expect_error(
    MatchOmics(toy_omics[1, ], toy_outcome[-1], toy_heterogeneity),
    "same length"
  )
  expect_error(
    MatchOmics(toy_omics[1, ], toy_outcome, toy_heterogeneity[-1, ]),
    "one row per subject"
  )
})
