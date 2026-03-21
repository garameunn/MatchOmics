#' MatchOmics main analysis function
#'
#' @param marker Numeric vector. Target marker to test.
#' @param outcome Binary outcome vector coded as 0/1.
#' @param omics_type One of "microbiome", "proteome", "metabolome",
#'   "transcriptome", or "custom".
#' @param heterogeneity Data frame or matrix of heterogeneity covariates.
#' @param covariates Optional data frame of additional adjustment covariates.
#' @param ids Optional subject IDs.
#' @param dichotomize Logical; if TRUE, dichotomize marker by median.
#' @param match_ratio Integer; default 1 for 1:1 matching.
#' @param replace Logical; default TRUE.
#' @param caliper Numeric or NULL; default 0.2.
#' @param distance Distance method for MatchIt; default "glm".
#' @param corstr Working correlation for GEE; default "exchangeable".
#' @param sandwich Logical; use sandwich SE; default TRUE.
#' @param gee_family Family for GEE; default binomial.
#' @param return_matched Logical; if TRUE, return matched data.
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
  
  ## ----------------------------
  ## input validation
  ## ----------------------------
  if (missing(marker)) {
    stop("`marker` must be provided.")
  }
  if (missing(outcome)) {
    stop("`outcome` must be provided.")
  }
  if (missing(heterogeneity)) {
    stop("`heterogeneity` must be provided.")
  }
  
  if (!is.numeric(marker)) {
    stop("`marker` must be a numeric vector.")
  }
  if (is.matrix(marker) || is.data.frame(marker)) {
    stop("`marker` must be a one-dimensional numeric vector, not a matrix or data frame.")
  }
  
  if (!is.numeric(outcome) && !is.integer(outcome)) {
    stop("`outcome` must be a binary vector coded as 0/1.")
  }
  if (is.matrix(outcome) || is.data.frame(outcome)) {
    stop("`outcome` must be a one-dimensional vector, not a matrix or data frame.")
  }
  
  if (length(marker) == 0) {
    stop("`marker` must have positive length.")
  }
  if (length(outcome) == 0) {
    stop("`outcome` must have positive length.")
  }
  if (length(marker) != length(outcome)) {
    stop("`marker` and `outcome` must have the same length.")
  }
  
  if (anyNA(outcome)) {
    stop("`outcome` must not contain missing values.")
  }
  if (!all(outcome %in% c(0, 1))) {
    stop("`outcome` must be coded as 0 and 1 only.")
  }
  
  if (!is.logical(dichotomize) || length(dichotomize) != 1) {
    stop("`dichotomize` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(replace) || length(replace) != 1) {
    stop("`replace` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(sandwich) || length(sandwich) != 1) {
    stop("`sandwich` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(return_matched) || length(return_matched) != 1) {
    stop("`return_matched` must be a single logical value (TRUE/FALSE).")
  }
  
  if (!is.numeric(match_ratio) || length(match_ratio) != 1 ||
      is.na(match_ratio) || match_ratio < 1 || match_ratio %% 1 != 0) {
    stop("`match_ratio` must be a positive integer.")
  }
  
  if (!is.null(caliper)) {
    if (!is.numeric(caliper) || length(caliper) != 1 || is.na(caliper) || caliper <= 0) {
      stop("`caliper` must be NULL or a single positive numeric value.")
    }
  }
  
  if (!is.character(distance) || length(distance) != 1) {
    stop("`distance` must be a single character string.")
  }
  if (!is.character(corstr) || length(corstr) != 1) {
    stop("`corstr` must be a single character string.")
  }
  
  if (!inherits(gee_family, "family")) {
    stop("`gee_family` must be a valid family object, e.g. stats::binomial().")
  }
  
  if (!is.data.frame(heterogeneity) && !is.matrix(heterogeneity)) {
    stop("`heterogeneity` must be a data frame or matrix.")
  }
  heterogeneity <- as.data.frame(heterogeneity)
  
  if (nrow(heterogeneity) != length(marker)) {
    stop("`heterogeneity` must have the same number of rows as the length of `marker`.")
  }
  if (ncol(heterogeneity) < 1) {
    stop("`heterogeneity` must contain at least one covariate column.")
  }
  
  if (is.null(colnames(heterogeneity))) {
    colnames(heterogeneity) <- paste0("H", seq_len(ncol(heterogeneity)))
  }
  if (any(colnames(heterogeneity) == "")) {
    colnames(heterogeneity)[colnames(heterogeneity) == ""] <-
      paste0("H", which(colnames(heterogeneity) == ""))
  }
  
  if (anyNA(heterogeneity)) {
    stop("`heterogeneity` must not contain missing values.")
  }
  
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates) && !is.matrix(covariates)) {
      stop("`covariates` must be a data frame or matrix when provided.")
    }
    covariates <- as.data.frame(covariates)
    
    if (nrow(covariates) != length(marker)) {
      stop("`covariates` must have the same number of rows as the length of `marker`.")
    }
    if (is.null(colnames(covariates))) {
      colnames(covariates) <- paste0("C", seq_len(ncol(covariates)))
    }
    if (anyNA(covariates)) {
      stop("`covariates` must not contain missing values.")
    }
  }
  
  if (is.null(ids)) {
    ids <- seq_along(marker)
  } else {
    if (length(ids) != length(marker)) {
      stop("`ids` must have the same length as `marker`.")
    }
    if (anyNA(ids)) {
      stop("`ids` must not contain missing values.")
    }
    if (anyDuplicated(ids)) {
      stop("`ids` must be unique.")
    }
  }
  
  if (all(is.na(marker))) {
    stop("`marker` contains only missing values.")
  }
  
  if (dichotomize) {
    marker_median <- stats::median(marker, na.rm = TRUE)
    if (is.na(marker_median)) {
      stop("Unable to compute the median of `marker`.")
    }
    
    exposure <- as.integer(marker > marker_median)
    
    if (length(unique(stats::na.omit(exposure))) < 2) {
      stop("Median dichotomization of `marker` did not produce two groups.")
    }
  } else {
    exposure <- marker
  }
  
  ## remove marker missingness only after defining exposure rule
  keep <- !is.na(marker)
  if (!all(keep)) {
    marker <- marker[keep]
    outcome <- outcome[keep]
    heterogeneity <- heterogeneity[keep, , drop = FALSE]
    ids <- ids[keep]
    exposure <- exposure[keep]
    if (!is.null(covariates)) {
      covariates <- covariates[keep, , drop = FALSE]
    }
  }
  
  ps_data <- tibble::tibble(
    .id = ids,
    exposure = exposure,
    outcome = outcome,
    marker = marker
  )
  
  if (!is.null(covariates)) {
    ps_data <- dplyr::bind_cols(ps_data, covariates)
  }
  
  ps_data <- dplyr::bind_cols(ps_data, heterogeneity)
  
  ps_formula <- stats::as.formula(
    paste("exposure ~", paste(colnames(heterogeneity), collapse = " + "))
  )
  
  ps_fit <- stats::glm(ps_formula, data = ps_data, family = stats::binomial())
  ps_data$ps <- stats::predict(ps_fit, type = "response")
  
  mobj <- MatchIt::matchit(
    formula = ps_formula,
    data = ps_data,
    method = "nearest",
    distance = ps_data$ps,
    replace = replace,
    ratio = match_ratio,
    caliper = caliper
  )
  
  matched <- MatchIt::match.data(mobj)
  
  if (nrow(matched) == 0) {
    stop("No matched subjects were obtained. Try relaxing the caliper or checking the input data.")
  }
  
  matched$.cluster <- build_matching_cluster(mobj)
  matched$.weight_norm <- normalize_match_weights(mobj, matched$.id)
  
  covar_term <- if (!is.null(covariates)) {
    paste(colnames(covariates), collapse = " + ")
  } else {
    NULL
  }
  
  gee_rhs <- c("marker", covar_term)
  gee_rhs <- gee_rhs[!is.null(gee_rhs) & gee_rhs != ""]
  gee_formula <- stats::as.formula(
    paste("outcome ~", paste(gee_rhs, collapse = " + "))
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
    coefficients = coef_tab,
    n_input = length(keep),
    n_used = nrow(ps_data),
    n_matched = nrow(matched)
  )
  
  class(out) <- "MatchOmics"
  out
}