# MatchOmics 0.3.1

* The outcome GEE now uses the continuous marker (`outcome ~ marker + ps`),
  as described in the manuscript and used in its simulation and real-data
  analyses. The median split (`marker_class`) only defines the groups to be
  matched. Previously the dichotomised marker entered the outcome model, so
  marker coefficients from 0.3.0 are on a different scale (up vs down) and
  are not comparable with 0.3.1 (per unit of the marker).
* `compute_smd()` gains a `weight_col` argument for the weighted
  standardised mean difference (e.g. `weight_col = "weight"` on
  `fit$matched_data`), matching the post-matching balance criterion of the
  manuscript. The default remains unweighted.

# MatchOmics 0.3.0

Breaking changes to the default outcome model, to align the package with
the model used throughout the manuscript.

* `MatchOmics()` now fits `outcome ~ marker_class + ps` by default
  (`adjust_ps = TRUE`), instead of `outcome ~ marker_class` alone. Set
  `adjust_ps = FALSE` to recover the previous default.
* New `outcome_covariates` argument: an optional data frame of additional
  nuisance covariates (e.g. age, sex, batch) added to the outcome GEE,
  mirroring the manuscript's real-data analysis
  (`outcome ~ marker + PS + covariates`).
* The outcome GEE is now fit with `geeM::geem()` instead of
  `geepack::geeglm()`, matching the engine used in the manuscript's
  simulation and real-data analysis scripts. `coef_table` and `coef()`
  output are unchanged in shape; `gee_fit` is now a `geem` object rather
  than a `geeglm` object.

Numerical results from 0.2.0 (marker-only, geeglm) will differ from 0.3.0
defaults (marker + PS, geem).

# MatchOmics 0.2.0

Initial release.
