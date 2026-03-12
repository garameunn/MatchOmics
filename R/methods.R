#' @export
print.MatchOmics <- function(x, ...) {
  cat("MatchOmics object\n")
  cat("Omics type:", x$omics_type, "\n")
  cat("Matched subjects:", if (!is.null(x$matched_data)) nrow(x$matched_data) else "not returned", "\n")
  cat("GEE coefficient rows:", nrow(x$coefficients), "\n")
  invisible(x)
}

#' @export
summary.MatchOmics <- function(object, ...) {
  list(
    omics_type = object$omics_type,
    ps_model = summary(object$ps_fit),
    gee_model = summary(object$gee_fit),
    coefficients = object$coefficients
  )
}
