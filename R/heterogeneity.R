#' Compute heterogeneity indices from an omics matrix
#'
#' Derives the heterogeneity covariates used as propensity score predictors.
#' For proteomics/metabolomics: first two principal components (PCA).
#' For metagenomics: Shannon and Simpson diversity indices (requires
#' \pkg{vegan}).
#'
#' @param omics_matrix Numeric matrix, features × subjects (rows = features,
#'   columns = subjects).
#' @param omics_type Character. One of \code{"proteomics"}, \code{"metabolomics"},
#'   or \code{"metagenomics"}.
#'
#' @return Data frame with one row per subject and two heterogeneity columns.
#'   Column names are \code{het1}/\code{het2} (proteomics/metabolomics) or
#'   \code{shannon}/\code{simpson} (metagenomics).
#'
#' @examples
#' data(toy_omics)
#' het <- compute_heterogeneity(toy_omics, omics_type = "proteomics")
#' head(het)
#'
#' @export
compute_heterogeneity <- function(omics_matrix,
                                  omics_type = c("proteomics",
                                                 "metabolomics",
                                                 "metagenomics")) {
  omics_type <- match.arg(omics_type)

  if (!is.matrix(omics_matrix) && !is.data.frame(omics_matrix))
    stop("'omics_matrix' must be a matrix or data frame (features x subjects).")

  if (omics_type %in% c("proteomics", "metabolomics")) {
    pp  <- stats::prcomp(t(omics_matrix), scale. = TRUE)
    out <- data.frame(het1 = pp$x[, 1], het2 = pp$x[, 2])
    rownames(out) <- colnames(omics_matrix)
    return(out)
  }

  # metagenomics
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("Package 'vegan' is required for metagenomics heterogeneity. ",
         "Install with: install.packages('vegan')")
  shannon <- vegan::diversity(t(omics_matrix), index = "shannon")
  simpson <- vegan::diversity(t(omics_matrix), index = "simpson")
  out <- data.frame(shannon = shannon, simpson = simpson)
  rownames(out) <- colnames(omics_matrix)
  out
}
