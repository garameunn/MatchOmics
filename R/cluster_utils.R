#' Build matching clusters from a MatchIt object
#' @param mobj A MatchIt object
#' @return integer vector of cluster IDs named by subject ID
#' @keywords internal
build_matching_cluster <- function(mobj) {
  mm <- mobj$match.matrix

  if (is.null(mm)) {
    stop("match.matrix is missing in MatchIt object.")
  }

  treated_ids <- rownames(mm)
  control_ids <- as.vector(mm[, 1])

  edges <- stats::na.omit(data.frame(
    from = treated_ids,
    to = control_ids,
    stringsAsFactors = FALSE
  ))

  nodes <- unique(c(edges$from, edges$to))
  parent <- stats::setNames(as.list(nodes), nodes)

  find_root <- function(x) {
    while (!identical(parent[[x]], x)) {
      x <- parent[[x]]
    }
    x
  }

  union_root <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (!identical(ra, rb)) parent[[rb]] <<- ra
  }

  for (i in seq_len(nrow(edges))) {
    union_root(edges$from[i], edges$to[i])
  }

  roots <- vapply(nodes, find_root, character(1))
  cluster <- as.integer(factor(roots))
  stats::setNames(cluster, nodes)
}
