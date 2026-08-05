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
