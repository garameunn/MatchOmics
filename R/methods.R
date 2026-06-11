#' @export
print.MatchOmics <- function(x, ...) {
  cat("MatchOmics result\n")
  cat("  Method   :", x$method,
      if (x$method == "two_round") "(two-round rescue matching)"
      else "(single-round without replacement)", "\n")
  cat("  Caliper  :", if (is.null(x$caliper)) "none" else x$caliper, "\n")
  cat("  Corstr   :", x$corstr, "\n")
  cat("  N input  :", x$n_input, "\n")
  cat("  N matched:", x$n_matched, "\n")
  cat("  Neff     :", round(x$neff, 1), "\n")
  invisible(x)
}

#' @export
summary.MatchOmics <- function(object, ...) {
  ct <- object$coef_table
  marker_row <- ct[grep("marker_class", rownames(ct), fixed = TRUE), , drop = FALSE]

  structure(
    list(
      coef_table = ct,
      marker_row = marker_row,
      neff       = object$neff,
      n_matched  = object$n_matched,
      method     = object$method,
      caliper    = object$caliper,
      corstr     = object$corstr
    ),
    class = "summary.MatchOmics"
  )
}

#' @export
print.summary.MatchOmics <- function(x, ...) {
  cat("MatchOmics summary\n")
  cat("  Method:", x$method, "| Caliper:", if (is.null(x$caliper)) "none" else x$caliper,
      "| Corstr:", x$corstr, "\n")
  cat("  N matched:", x$n_matched, "| Neff:", round(x$neff, 1), "\n\n")
  cat("GEE coefficients:\n")
  print(x$coef_table)
  cat("\nMarker effect (marker_class):\n")
  print(x$marker_row)
  invisible(x)
}

#' @export
coef.MatchOmics <- function(object, ...) {
  stats::coef(object$gee_fit)
}
