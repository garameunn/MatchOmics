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
#'   (proposed method); \code{"standard"} uses single-round without-replacement
#'   matching.
#' @param corstr   GEE correlation structure passed to \code{geepack::geeglm}.
#'   One of \code{"independence"} (default) or \code{"exchangeable"}.
#'
#' @return A \code{MatchOmics} object (list) with components:
#'   \describe{
#'     \item{ps_model}{Fitted PS model (MatchIt object).}
#'     \item{match_result}{Raw output of the matching step.}
#'     \item{matched_data}{Data frame of matched subjects with weights.}
#'     \item{gee_fit}{Fitted \code{geeglm} object.}
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
                       corstr   = c("independence", "exchangeable")) {

  method <- match.arg(method)
  corstr <- match.arg(corstr)

  n <- length(marker)
  if (length(outcome) != n)
    stop("'marker' and 'outcome' must have the same length.")
  if (!is.data.frame(heterogeneity) || nrow(heterogeneity) != n)
    stop("'heterogeneity' must be a data frame with one row per subject.")
  if (!all(outcome %in% c(0, 1, NA)))
    stop("'outcome' must be binary (0/1).")

  # Assign integer rownames for consistent subsetting
  outcome <- as.integer(outcome)
  rnames  <- as.character(seq_len(n))
  names(marker)  <- rnames
  names(outcome) <- rnames
  rownames(heterogeneity) <- rnames

  # Dichotomise marker at median
  marker_class <- ifelse(marker >= median(marker, na.rm = TRUE), "up", "down")
  names(marker_class) <- rnames

  # Build PS model data
  ps_data   <- data.frame(marker_class = as.factor(marker_class),
                          heterogeneity,
                          row.names = rnames)
  ps_form   <- as.formula(paste("marker_class ~",
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

  # Mean-normalise weights (guarantees sum(w) = N for GEE)
  matched_df$weight <- matched_df$weight / mean(matched_df$weight)

  # Effective sample size
  neff_val <- compute_neff(matched_df$weight)

  # Weighted GEE
  gee_fit <- geepack::geeglm(
    outcome ~ marker_class,
    data    = matched_df,
    id      = cluster,
    weights = weight,
    family  = binomial(),
    corstr  = corstr
  )

  coef_tbl <- as.data.frame(summary(gee_fit)$coefficients)

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
      corstr       = corstr
    ),
    class = "MatchOmics"
  )
}
