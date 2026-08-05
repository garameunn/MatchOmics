#' Propensity score matching with weighted GEE for omics association analysis
#'
#' Performs two-round (or standard single-round) nearest-neighbour propensity
#' score matching on a dichotomised omics marker, then tests for association
#' with the binary outcome via weighted GEE.
#'
#' @param marker   Numeric vector of length n. A single omics feature (e.g. one
#'   protein or metabolite). Dichotomised at its median internally.
#' @param outcome  Integer/numeric 0/1 vector of length n. Binary outcome.
#' @param heterogeneity  Data frame (n rows). Heterogeneity indices used as PS
#'   covariates (e.g. output of \code{\link{compute_heterogeneity}}).
#' @param caliper  Numeric. Caliper width on the PS scale. \code{NULL} removes
#'   the caliper. Default \code{0.3}.
#' @param method   \code{"two_round"} (default) uses two-round rescue matching
#'   (proposed method; termed "with-replacement matching" in the manuscript,
#'   since a control may be reused across rounds 1 and 2); \code{"standard"}
#'   uses single-round without-replacement matching (round 1 only).
#' @param corstr   GEE correlation structure passed to \code{geeM::geem}.
#'   One of \code{"independence"} (default) or \code{"exchangeable"}.
#' @param adjust_ps  Logical. If \code{TRUE} (default), the propensity score
#'   is included as a covariate in the outcome GEE, matching the model used
#'   throughout the manuscript (\code{outcome ~ marker_class + ps}). Set to
#'   \code{FALSE} to fit \code{outcome ~ marker_class} only.
#' @param outcome_covariates  Optional data frame (n rows, one row per
#'   subject, same row order as \code{marker}/\code{outcome}) of additional
#'   nuisance covariates to include in the outcome GEE (e.g. age, sex, batch).
#'   \code{NULL} (default) adds none.
#'
#' @return A \code{MatchOmics} object (list) with components:
#'   \describe{
#'     \item{ps_model}{Fitted PS model (MatchIt object).}
#'     \item{match_result}{Raw output of the matching step.}
#'     \item{matched_data}{Data frame of matched subjects with weights.}
#'     \item{gee_fit}{Fitted \code{geem} object.}
#'     \item{coef_table}{Coefficient table from GEE summary.}
#'     \item{neff}{Kish effective sample size after matching.}
#'     \item{n_input}{Number of input subjects.}
#'     \item{n_matched}{Number of subjects retained after matching.}
#'   }
#'
#' @examples
#' data(toy_omics)
#' data(toy_outcome)
#' data(toy_heterogeneity)
#'
#' fit <- MatchOmics(
#'   marker        = toy_omics[1, ],
#'   outcome       = toy_outcome,
#'   heterogeneity = toy_heterogeneity,
#'   caliper       = 0.3,
#'   method        = "two_round",
#'   corstr        = "independence"
#' )
#' summary(fit)
#'
#' @export
MatchOmics <- function(marker,
                       outcome,
                       heterogeneity,
                       caliper  = 0.3,
                       method   = c("two_round", "standard"),
                       corstr   = c("independence", "exchangeable"),
                       adjust_ps = TRUE,
                       outcome_covariates = NULL) {

  method <- match.arg(method)
  corstr <- match.arg(corstr)

  n <- length(marker)
  if (length(outcome) != n)
    stop("'marker' and 'outcome' must have the same length.")
  if (!is.data.frame(heterogeneity) || nrow(heterogeneity) != n)
    stop("'heterogeneity' must be a data frame with one row per subject.")
  if (!all(outcome %in% c(0, 1, NA)))
    stop("'outcome' must be binary (0/1).")
  if (!is.null(outcome_covariates) &&
      (!is.data.frame(outcome_covariates) || nrow(outcome_covariates) != n))
    stop("'outcome_covariates' must be a data frame with one row per subject.")

  # Assign integer rownames for consistent subsetting
  outcome <- as.integer(outcome)
  rnames  <- as.character(seq_len(n))
  names(marker)  <- rnames
  names(outcome) <- rnames
  rownames(heterogeneity) <- rnames
  if (!is.null(outcome_covariates)) rownames(outcome_covariates) <- rnames

  # Dichotomise marker at median
  marker_class <- ifelse(marker >= stats::median(marker, na.rm = TRUE), "up", "down")
  names(marker_class) <- rnames

  # Build PS model data
  ps_data   <- data.frame(marker_class = as.factor(marker_class),
                          heterogeneity,
                          row.names = rnames)
  ps_form   <- stats::as.formula(paste("marker_class ~",
                                paste(names(heterogeneity), collapse = " + ")))

  # Match
  tmr <- if (method == "two_round")
    two_round_match(ps_data, ps_form, caliper = caliper)
  else
    standard_match(ps_data, ps_form, caliper = caliper)

  pairs <- rbind(tmr$pairs1,
                 if (!is.null(tmr$pairs2)) tmr$pairs2 else NULL)

  if (is.null(pairs) || nrow(pairs) == 0)
    stop("No matched pairs found. Consider increasing or removing 'caliper'.")

  # Build clusters and weights
  all_ids   <- unique(c(pairs$treated, pairs$control))
  clusters  <- build_clusters_from_pairs(pairs)
  w_vec     <- build_weights(all_ids, pairs)
  ps_vec    <- tmr$m1$distance
  round2_ids <- if (!is.null(tmr$pairs2))
    unique(c(tmr$pairs2$treated, tmr$pairs2$control)) else character(0)

  # Assemble matched data frame
  matched_df <- data.frame(
    id           = all_ids,
    outcome      = outcome[all_ids],
    marker       = marker[all_ids],
    marker_class = marker_class[all_ids],
    cluster      = as.factor(clusters[all_ids]),
    weight       = w_vec[all_ids],
    ps           = ps_vec[match(all_ids, names(ps_vec))],
    round2       = as.integer(all_ids %in% round2_ids),
    stringsAsFactors = FALSE,
    row.names    = all_ids
  )
  matched_df <- matched_df[order(as.integer(matched_df$cluster)), ]

  # Append heterogeneity columns (needed for post-matching SMD diagnostics)
  het_sub <- heterogeneity[match(rownames(matched_df), rownames(heterogeneity)), , drop = FALSE]
  rownames(het_sub) <- rownames(matched_df)
  matched_df <- cbind(matched_df, het_sub)

  # Append user-supplied outcome-model covariates
  outcome_terms <- character(0)
  if (!is.null(outcome_covariates)) {
    cov_sub <- outcome_covariates[match(rownames(matched_df), rownames(outcome_covariates)), ,
                                  drop = FALSE]
    rownames(cov_sub) <- rownames(matched_df)
    matched_df <- cbind(matched_df, cov_sub)
    outcome_terms <- names(outcome_covariates)
  }

  # Mean-normalise weights (guarantees sum(w) = N for GEE)
  matched_df$weight <- matched_df$weight / mean(matched_df$weight)

  # Effective sample size
  neff_val <- compute_neff(matched_df$weight)

  # Outcome model: marker_class [+ ps] [+ outcome_covariates], matching the
  # model used throughout the manuscript (outcome ~ marker + PS + covariates)
  rhs <- c("marker_class", if (adjust_ps) "ps", outcome_terms)
  outcome_form <- stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))

  # Weighted GEE (geeM::geem — matches the engine used in the manuscript
  # simulation and real-data analysis scripts, unlike geepack::geeglm)
  # suppressWarnings: weighted binomial GEE always triggers "non-integer #successes"
  # because weights shift the effective trial count off integers — results are unaffected.
  gee_fit <- suppressWarnings(
    geeM::geem(
      formula  = outcome_form,
      data     = matched_df,
      id       = matched_df$cluster,
      weights  = matched_df$weight,
      family   = stats::binomial(),
      corstr   = corstr,
      sandwich = TRUE,
      useP     = TRUE
    )
  )

  gee_summary <- summary(gee_fit)
  coef_tbl <- data.frame(
    Estimate    = gee_summary$beta,
    `Std.err`   = gee_summary$se.robust,
    Wald        = gee_summary$wald.test,
    `Pr(>|W|)`  = gee_summary$p,
    row.names   = gee_summary$coefnames,
    check.names = FALSE
  )

  structure(
    list(
      ps_model     = tmr$m1,
      match_result = tmr,
      matched_data = matched_df,
      gee_fit      = gee_fit,
      coef_table   = coef_tbl,
      neff         = neff_val,
      n_input      = n,
      n_matched    = nrow(matched_df),
      caliper      = caliper,
      method       = method,
      corstr       = corstr,
      adjust_ps    = adjust_ps
    ),
    class = "MatchOmics"
  )
}

# 'smd'/'covariate'/'timing' are resolved via NSE against `data =` in aes()
# at runtime — not undefined globals.
utils::globalVariables(c("smd", "covariate", "timing"))
