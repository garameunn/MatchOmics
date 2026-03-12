#' Normalize matching weights
#' @param mobj MatchIt object
#' @param ids subject IDs in matched data
#' @return normalized numeric vector
#' @keywords internal
normalize_match_weights <- function(mobj, ids) {
  w <- mobj$weights
  w <- w[match(ids, names(w))]
  w[is.na(w)] <- 0

  s <- sum(w)
  if (s == 0) {
    stop("All matching weights are zero.")
  }

  w / s
}
