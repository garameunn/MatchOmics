# Internal matching functions — not exported

# ── two-round matching ─────────────────────────────────────────────────────────
# Round 1 : nearest-neighbour without replacement (strict, high quality)
# Round 2 : unmatched treated  +  ALL controls re-eligible
#           rescues cases blocked by round-1 selection
two_round_match <- function(data, formula, caliper) {

  m1 <- MatchIt::matchit(formula, data = data,
                         method  = "nearest",
                         replace = FALSE,
                         caliper = caliper)

  matched1_ids    <- rownames(MatchIt::match.data(m1))
  resp_col        <- as.character(formula[[2]])
  treated_level   <- levels(data[[resp_col]])[2]
  control_level   <- levels(data[[resp_col]])[1]
  treated_all     <- rownames(data)[data[[resp_col]] == treated_level]
  unmatched_treat <- setdiff(treated_all, matched1_ids)

  pairs1 <- .extract_pairs(m1$match.matrix, round = 1L)

  result <- list(m1 = m1, m2 = NULL, pairs1 = pairs1, pairs2 = NULL)
  if (length(unmatched_treat) == 0) return(result)

  controls_all <- rownames(data)[data[[resp_col]] == control_level]
  data2        <- data[c(unmatched_treat, controls_all), ]
  n_t2 <- sum(data2[[resp_col]] == treated_level)
  n_c2 <- sum(data2[[resp_col]] == control_level)

  if (n_t2 > 0 && n_c2 > 0) {
    m2 <- tryCatch(
      MatchIt::matchit(formula, data = data2,
                       method  = "nearest",
                       replace = FALSE,
                       caliper = caliper),
      error = function(e) NULL
    )
    if (!is.null(m2)) {
      result$m2     <- m2
      result$pairs2 <- .extract_pairs(m2$match.matrix, round = 2L)
    }
  }
  result
}

# ── standard single-round without replacement ──────────────────────────────────
standard_match <- function(data, formula, caliper) {
  m <- MatchIt::matchit(formula, data = data,
                        method  = "nearest",
                        replace = FALSE,
                        caliper = caliper)
  list(m1 = m, m2 = NULL,
       pairs1 = .extract_pairs(m$match.matrix, round = 1L),
       pairs2 = NULL)
}

# ── helpers ────────────────────────────────────────────────────────────────────
.extract_pairs <- function(match_matrix, round) {
  df <- data.frame(treated = rownames(match_matrix),
                   control = match_matrix[, 1],
                   round   = round,
                   stringsAsFactors = FALSE)
  df[!is.na(df$control), ]
}

# Union-find cluster builder (merges pairs from both rounds)
build_clusters_from_pairs <- function(pairs_df) {
  if (is.null(pairs_df) || nrow(pairs_df) == 0)
    return(stats::setNames(integer(0), character(0)))

  edges  <- unique(pairs_df[, c("treated", "control")])
  nodes  <- unique(c(edges$treated, edges$control))
  parent <- stats::setNames(as.list(nodes), nodes)

  find_root <- function(x) {
    while (!identical(parent[[x]], x)) x <- parent[[x]]
    x
  }
  union_root <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (!identical(ra, rb)) parent[[rb]] <<- ra
  }
  for (i in seq_len(nrow(edges))) union_root(edges$treated[i], edges$control[i])

  roots   <- vapply(nodes, find_root, character(1))
  cluster <- as.integer(factor(roots))
  stats::setNames(cluster, nodes)
}

# Mean-normalised weights: w_j = m_j / mean(m)  →  sum(w) = N
# ensures GEE sees the correct effective sample size
build_weights <- function(all_ids, pairs_df) {
  freq <- table(c(pairs_df$treated, pairs_df$control))
  m    <- as.integer(freq[match(all_ids, names(freq))])
  m[is.na(m)] <- 0L
  active <- m[m > 0]
  w <- m / mean(active)
  stats::setNames(w, all_ids)
}
