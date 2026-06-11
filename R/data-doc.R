#' Simulated proteomics matrix
#'
#' A 50 × 200 numeric matrix of log-normalised protein abundances for 200
#' simulated subjects and 50 proteins. Rows are proteins, columns are subjects.
#'
#' @format A numeric matrix with 50 rows and 200 columns.
#' @examples
#' data(toy_omics)
#' dim(toy_omics)
"toy_omics"

#' Simulated binary outcome
#'
#' Integer vector of length 200 (0/1) representing disease status for the
#' subjects in \code{\link{toy_omics}}.
#'
#' @format An integer vector of length 200.
#' @examples
#' data(toy_outcome)
#' table(toy_outcome)
"toy_outcome"

#' Simulated heterogeneity indices
#'
#' Data frame with 200 rows and two columns (\code{het1}, \code{het2}) derived
#' from PCA of \code{\link{toy_omics}}. Used as propensity score covariates in
#' \code{\link{MatchOmics}}.
#'
#' @format A data frame with 200 rows and 2 columns:
#'   \describe{
#'     \item{het1}{First principal component score.}
#'     \item{het2}{Second principal component score.}
#'   }
#' @examples
#' data(toy_heterogeneity)
#' head(toy_heterogeneity)
"toy_heterogeneity"
