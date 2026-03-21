#' @export
print.MatchOmics <- function(x, ...) {
  cat("MatchOmics object\n")
  cat("-----------------\n")
  cat("Omics type      :", x$omics_type, "\n")
  cat("Input samples   :", x$n_input, "\n")
  cat("Samples used    :", x$n_used, "\n")
  cat("Matched samples :", x$n_matched, "\n")
  cat("GEE fitted      :", inherits(x$gee_fit, "geeglm"), "\n")
  invisible(x)
}

#' @export
summary.MatchOmics <- function(object, ...) {
  
  coef_tab <- as.data.frame(object$coefficients)
  
  if (nrow(coef_tab) > 0) {
    coef_tab$term <- rownames(coef_tab)
    rownames(coef_tab) <- NULL
    coef_tab <- coef_tab[, c("term", setdiff(colnames(coef_tab), "term"))]
  }
  
  marker_row <- coef_tab[coef_tab$term == "marker", , drop = FALSE]
  
  out <- list(
    omics_type = object$omics_type,
    n_input = object$n_input,
    n_used = object$n_used,
    n_matched = object$n_matched,
    marker_result = marker_row,
    coefficients = coef_tab
  )
  
  class(out) <- "summary.MatchOmics"
  out
}

#' @export
print.summary.MatchOmics <- function(x, ...) {
  cat("Summary of MatchOmics\n")
  cat("---------------------\n")
  cat("Omics type      :", x$omics_type, "\n")
  cat("Input samples   :", x$n_input, "\n")
  cat("Samples used    :", x$n_used, "\n")
  cat("Matched samples :", x$n_matched, "\n\n")
  
  cat("Marker effect\n")
  cat("-------------\n")
  
  if (nrow(x$marker_result) == 0) {
    cat("Marker term was not found in the fitted GEE model.\n")
  } else {
    print(x$marker_result, row.names = FALSE)
  }
  
  invisible(x)
}