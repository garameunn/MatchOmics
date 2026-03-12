#' MatchOmics main analysis function
#'
#' @param marker numeric vector. Target marker to test.
#' @param outcome binary outcome vector (0/1).
#' @param omics_type one of "microbiome", "proteome", "metabolome", "transcriptome", "custom".
#' @param heterogeneity data.frame or matrix of heterogeneity covariates.
#' @param covariates optional data.frame of additional adjustment covariates.
#' @param ids optional subject IDs.
#' @param dichotomize logical; if TRUE, dichotomize marker by median.
#' @param match_ratio integer; default 1 for 1:1 matching.
#' @param replace logical; default TRUE.
#' @param caliper numeric or NULL; default 0.2.
#' @param distance distance method for MatchIt; default "glm".
#' @param corstr working correlation for GEE; default "exchangeable".
#' @param sandwich logical; use sandwich SE; default TRUE.
#' @param gee_family family for GEE; default binomial.
#' @param return_matched logical; if TRUE, return matched data.
#'
#' @return A list with propensity score model, match object, matched data,
#'   GEE fit, coefficient table, and metadata.
#' @export
MatchOmics <- function(
  marker,
  outcome,
  omics_type = c("microbiome", "proteome", "metabolome", "transcriptome", "custom"),
  heterogeneity,
  covariates = NULL,
  ids = NULL,
  dichotomize = TRUE,
  match_ratio = 1L,
  replace = TRUE,
  caliper = 0.2,
  distance = "glm",
  corstr = "exchangeable",
  sandwich = TRUE,
  gee_family = stats::binomial(),
  return_matched = TRUE
) {
  omics_type <- match.arg(omics_type)

  if (is.null(ids)) {
    ids <- seq_along(marker)
  }

  if (length(marker) != length(outcome)) {
    stop("marker and outcome must have the same length.")
  }

  if (!is.null(covariates) && nrow(as.data.frame(covariates)) != length(marker)) {
    stop("covariates must have the same number of rows as marker length.")
  }

  heterogeneity <- as.data.frame(heterogeneity)
  if (nrow(heterogeneity) != length(marker)) {
    stop("heterogeneity must have the same number of rows as marker length.")
  }

  exposure <- if (dichotomize) {
    as.integer(marker > stats::median(marker, na.rm = TRUE))
  } else {
    marker
  }

  ps_data <- tibble::tibble(
    .id = ids,
    exposure = exposure,
    outcome = outcome,
    marker = marker
  )

  if (!is.null(covariates)) {
    ps_data <- dplyr::bind_cols(ps_data, as.data.frame(covariates))
  }

  ps_data <- dplyr::bind_cols(ps_data, heterogeneity)

  ps_formula <- stats::as.formula(
    paste("exposure ~", paste(colnames(heterogeneity), collapse = " + "))
  )

  ps_fit <- stats::glm(ps_formula, data = ps_data, family = stats::binomial())
  ps_data$ps <- stats::predict(ps_fit, type = "response")

  match_formula <- stats::as.formula("exposure ~ ps")

  mobj <- MatchIt::matchit(
    formula = match_formula,
    data = ps_data,
    method = "nearest",
    distance = ps_data$ps,
    replace = replace,
    ratio = match_ratio,
    caliper = caliper
  )

  matched <- MatchIt::match.data(mobj)
  matched$.cluster <- build_matching_cluster(mobj)
  matched$.weight_norm <- normalize_match_weights(mobj, matched$.id)

  gee_formula <- stats::as.formula(
    paste("outcome ~ marker",
          if (!is.null(covariates)) paste("+", paste(colnames(covariates), collapse = " + ")) else "")
  )

  gee_fit <- geepack::geeglm(
    formula = gee_formula,
    data = matched,
    id = matched$.cluster,
    family = gee_family,
    corstr = corstr,
    weights = matched$.weight_norm,
    std.err = if (sandwich) "san.se" else "naive"
  )

  coef_tab <- summary(gee_fit)$coefficients

  out <- list(
    call = match.call(),
    omics_type = omics_type,
    ps_fit = ps_fit,
    matchit = mobj,
    matched_data = if (return_matched) matched else NULL,
    gee_fit = gee_fit,
    coefficients = coef_tab
  )

  class(out) <- "MatchOmics"
  out
}
