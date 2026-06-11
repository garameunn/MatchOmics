## ============================================================
##  02_simulation_generator.R
##  Changes vs original:
##  1. Scaling bug fix: mm = m / mean(m)  [was m / sum(m)]
##  2. wrep functions: two-round matching (rescue unmatched)
##  3. round2 flag column added to matched data
##  4. Neff = (sum w)^2 / sum(w^2) computed and returned
##  5. TEST_MODE support (rep=10, dimm=5, S1/N200/rate0.2 only)
##  6. prevalence rates: 0.1, 0.2, 0.3 (= 10%, 20%, 30%)
## ============================================================

# ── TEST MODE flag ─────────────────────────────────────────
# To run a quick sanity check:
#   TEST_MODE <- TRUE; source("02_simulation_generator.R")
# Or from terminal:
#   Rscript 03_Type1Error_Calculating.R --test
if (!exists("TEST_MODE")) TEST_MODE <- FALSE

# ── helpers ────────────────────────────────────────────────
pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  myTaxa  <- allTaxa[!(allTaxa %in% badTaxa)]
  prune_taxa(myTaxa, physeq)
}

f <- function(idx) {
  set.seed(seeds[idx])
}

indc <- function(dt, omics) {
  if (omics == "Metagenomics") {
    shannon <- vegan::diversity(t(dt), index = "shannon")
    simpson <- vegan::diversity(t(dt), index = "simpson")
    d <- data.frame(shannon, simpson)
  } else if (omics %in% c("Proteomics", "Metabolomics")) {
    pp <- prcomp(t(dt), scale = TRUE)
    d  <- data.frame(pp$x[, 1:2])
  }
  return(d)
}

trnsf <- function(dt, omics, lb.size) {
  if (omics == "Metagenomics") {
    d <- cpm(dt, log = TRUE, lib.size = lb.size)
  } else {
    d <- dt
  }
  # force plain matrix (otu_table / other S4 objects break cov/unlist)
  d <- matrix(as.numeric(d),
              nrow     = nrow(d),
              ncol     = ncol(d),
              dimnames = list(rownames(d), colnames(d)))
  return(d)
}

# ── SMD helper (for diagnostics) ───────────────────────────
compute_smd <- function(df, covnames) {
  # returns named vector of SMD per covariate
  sapply(covnames, function(v) {
    x <- df[[v]]
    g <- df[["taxa_class_binary"]]   # 1=up, 0=down
    m1 <- mean(x[g == 1], na.rm = TRUE)
    m0 <- mean(x[g == 0], na.rm = TRUE)
    s1 <- var(x[g == 1],  na.rm = TRUE)
    s0 <- var(x[g == 0],  na.rm = TRUE)
    pooled_sd <- sqrt((s1 + s0) / 2)
    if (pooled_sd == 0) return(NA)
    (m1 - m0) / pooled_sd
  })
}

# ── Neff helper ────────────────────────────────────────────
compute_neff <- function(w) {
  # Kish (1965) effective sample size
  (sum(w))^2 / sum(w^2)
}

# ── two-round matching (wrep only) ─────────────────────────
# Round 1: nearest without replacement (strict, high quality)
# Round 2: unmatched treated + ALL controls re-eligible
#          → rescues cases blocked by prior selection
# Cluster: igraph connected components across both rounds
two_round_match <- function(formatch, matchform, caliper) {

  # ---------- round 1: without replacement ----------
  m1 <- MatchIt::matchit(matchform, data = formatch,
                         method   = "nearest",
                         replace  = FALSE,
                         caliper  = caliper)

  matched1_ids <- rownames(MatchIt::match.data(m1))
  # treated units that failed to match in round 1
  treated_all  <- rownames(formatch)[formatch$taxa_class == levels(formatch$taxa_class)[2]]
  unmatched_treated <- setdiff(treated_all, matched1_ids)

  result <- list(m1 = m1, m2 = NULL,
                 pairs1 = NULL, pairs2 = NULL)

  if (length(unmatched_treated) == 0) {
    # nothing to rescue
    pairs1 <- data.frame(
      treated = rownames(m1$match.matrix),
      control = m1$match.matrix[, 1],
      round   = 1L,
      stringsAsFactors = FALSE
    )
    pairs1 <- pairs1[!is.na(pairs1$control), ]
    result$pairs1 <- pairs1
    return(result)
  }

  # ---------- round 2 ----------
  # pool: unmatched treated  +  ALL controls (can be re-used)
  controls_all <- rownames(formatch)[formatch$taxa_class == levels(formatch$taxa_class)[1]]
  pool_ids     <- c(unmatched_treated, controls_all)
  formatch2    <- formatch[pool_ids, ]
  # need at least one treated and one control
  n_treated2 <- sum(formatch2$taxa_class == levels(formatch$taxa_class)[2])
  n_control2 <- sum(formatch2$taxa_class == levels(formatch$taxa_class)[1])

  if (n_treated2 > 0 && n_control2 > 0) {
    m2 <- tryCatch(
      MatchIt::matchit(matchform, data = formatch2,
                       method  = "nearest",
                       replace = FALSE,
                       caliper = caliper),
      error = function(e) NULL
    )
    result$m2 <- m2
  }

  # ---------- build pair tables ----------
  pairs1 <- data.frame(
    treated = rownames(m1$match.matrix),
    control = m1$match.matrix[, 1],
    round   = 1L,
    stringsAsFactors = FALSE
  )
  pairs1 <- pairs1[!is.na(pairs1$control), ]

  pairs2 <- NULL
  if (!is.null(result$m2)) {
    mm2 <- result$m2$match.matrix
    pairs2 <- data.frame(
      treated = rownames(mm2),
      control = mm2[, 1],
      round   = 2L,
      stringsAsFactors = FALSE
    )
    pairs2 <- pairs2[!is.na(pairs2$control), ]
  }

  result$pairs1 <- pairs1
  result$pairs2 <- pairs2
  return(result)
}

# ── cluster builder (igraph, merges both rounds) ───────────
build_clusters_from_pairs <- function(pairs_df) {
  # pairs_df: data.frame with columns treated, control, round
  if (is.null(pairs_df) || nrow(pairs_df) == 0)
    return(setNames(integer(0), character(0)))

  edges <- unique(data.frame(
    from = pairs_df$treated,
    to   = pairs_df$control,
    stringsAsFactors = FALSE
  ))

  nodes  <- unique(c(edges$from, edges$to))
  parent <- setNames(as.list(nodes), nodes)

  find_root <- function(x) {
    while (!identical(parent[[x]], x)) x <- parent[[x]]
    x
  }
  union_root <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (!identical(ra, rb)) parent[[rb]] <<- ra
  }
  for (i in seq_len(nrow(edges))) union_root(edges$from[i], edges$to[i])

  roots   <- vapply(nodes, find_root, character(1))
  cluster <- as.integer(factor(roots))
  setNames(cluster, nodes)
}

# ── weight builder ─────────────────────────────────────────
# m_j  = number of times subject j appears across all pairs
# w_j  = m_j / mean(m)   →  sum(w) = N  (GEE sees correct N)
build_weights <- function(all_ids, pairs_df) {
  freq <- table(c(pairs_df$treated, pairs_df$control))
  m    <- as.integer(freq[match(all_ids, names(freq))])
  m[is.na(m)] <- 0L
  w <- m / mean(m[m > 0])          # mean-normalised (NOT sum-to-1)
  setNames(w, all_ids)
}

# ═══════════════════════════════════════════════════════════
#  tablemaking_t1   (WITH replacement – two-round matching)
# ═══════════════════════════════════════════════════════════
tablemaking_t1 <- function(omics, trsftable, datatable, indicator,
                            rate, b, x, cl, scenario,
                            k = NULL, nonlin = "inv") {

  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))

  # ── two-round matching ──
  tmr    <- two_round_match(formatch, matchform, caliper = cl)
  pairs  <- rbind(tmr$pairs1,
                  if (!is.null(tmr$pairs2)) tmr$pairs2 else NULL)
  if (is.null(pairs) || nrow(pairs) == 0)
    return(NULL)

  all_ids  <- unique(c(pairs$treated, pairs$control))
  clusters <- build_clusters_from_pairs(pairs)
  w_vec    <- build_weights(all_ids, pairs)

  # round flag per subject (1 = round1 only, 2 = appeared in round2)
  round2_ids <- if (!is.null(tmr$pairs2))
    unique(c(tmr$pairs2$treated, tmr$pairs2$control)) else character(0)

  # PS from round-1 matchit (covers all subjects via $distance)
  ps_vec <- tmr$m1$distance

  # subclass (cluster) vector aligned to all_ids
  subclass <- clusters[all_ids]

  samdat <- data.frame(
    ps        = ps_vec[match(all_ids, names(ps_vec))],
    pair_taxa = subclass,
    variable  = taxa_class[match(all_ids, rownames(taxa_class)), ],
    value     = all_ids,
    round2    = as.integer(all_ids %in% round2_ids),
    stringsAsFactors = FALSE
  )

  # ── liability score setup ──
  trsftable_matched <- trsftable[, match(samdat[, 4], colnames(trsftable))]

  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS <- if (nonlin == "inv") 1 / samdat$ps else
        log(exp(samdat$ps) + log(samdat$ps) * ((samdat$ps)^-2 -
              ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS <- indicator[match(samdat[, 4], rownames(indicator)), 1] +
             indicator[match(samdat[, 4], rownames(indicator)), 2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x, ])), mean(PPS), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),   sd(PPS),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[4], sd = sds[4])
    v23     <- cov(unlist(trsftable_matched[x, ]), PPS)

  } else if (scenario == "S3") {
    otu_for_generate0 <- trsftable[k, match(samdat[, 4], colnames(trsftable))]
    otu_for_generate  <- matrix(ncol = ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_for_generate <- rbind(otu_for_generate,
        lm(t(otu_for_generate0)[, ii] ~ trsftable_matched[x, ])$residuals)
    otu_for_generate <- otu_for_generate[-1, ]
    means <- c(1, mean(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, mean), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, sd),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[length(means)],
                     sd = sds[length(sds)])
    vresid  <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_for_generate[y[1], ], otu_for_generate[y[2], ])))
    TT  <- sum(apply(otu_for_generate, 1, var)) + 2 * vresid
    v23 <- cov(unlist(trsftable_matched[x, ]),
               apply(otu_for_generate, 2, sum))
  }

  b1 <- 0.1; b2 <- 0

  listtable <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds[3]^2)) else
                                         sqrt(3.03 / (7  * sds[3]^2))
      L     <- b1 * x1 + b2 * unlist(trsftable_matched[x, ]) + b3 * PPS + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] + b3*means[3],
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3*sds[3])^2 + 2*b2*b3*v23 + 1))
    } else {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                         sqrt(3.03 / (7  * TT))
      L     <- b1 * x1 + b2 * unlist(trsftable_matched[x, ]) +
               b3 * apply(otu_for_generate, 2, sum) + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] +
                            b3*sum(apply(otu_for_generate, 1, mean)),
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3^2)*TT + 2*b2*b3*v23 + 1))
    }

    case    <- ifelse(L >= thres, "1", "0")
    samdat1 <- data.frame(
      ID        = samdat[, 4],
      case      = as.numeric(case[match(samdat[, 4], names(case))]),
      pair_taxa = samdat$pair_taxa,
      indOTU    = ifelse(samdat[, 3] == "up", "large", "small"),
      x1        = x1,
      PS        = samdat$ps,
      round2    = samdat$round2
    )
    samdat1 <- samdat1[!duplicated(samdat1$ID), ]
    rownames(samdat1) <- samdat1$ID

    trsftable1 <- trsftable[, match(rownames(samdat1), colnames(trsftable))]
    matchdata  <- data.frame(
      samdat1,
      hetero1 = indicator[match(rownames(samdat1), rownames(indicator)), 1],
      hetero2 = indicator[match(rownames(samdat1), rownames(indicator)), 2],
      marker  = unlist(trsftable1[x, ])
    )
    matchdata1 <- matchdata[order(matchdata$pair_taxa), ]
    matchdata1$pair_taxa <- as.factor(matchdata1$pair_taxa)

    # ── weights (mean-normalised) ──
    w_sub <- w_vec[match(rownames(matchdata1), names(w_vec))]
    w_sub[is.na(w_sub)] <- 1
    matchdata1$mm  <- w_sub / mean(w_sub)   # mean-normalised: sum = N
    matchdata1$mmm <- tmr$m1$weights[match(rownames(matchdata1),
                                           names(tmr$m1$weights))]

    # ── Neff ──
    matchdata1$.neff <- compute_neff(matchdata1$mm)

    return(matchdata1)
  })

  rm(means, sds, x1, epsilon, v23)

  # ── entire (unmatched) data for logit/combat baseline ──
  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS1 <- if (nonlin == "inv") 1 / tmr$m1$distance else
        log(exp(tmr$m1$distance) +
            log(tmr$m1$distance) *
            ((tmr$m1$distance)^-2 -
             ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS1 <- indicator[, 1] + indicator[, 2]
    }
    means1 <- c(1, mean(unlist(trsftable[x, ])), mean(PPS1), 0)
    sds1   <- c(1, sd(unlist(trsftable[x, ])),   sd(PPS1),   1)
    x11      <- rnorm(length(PPS1), mean = means1[1], sd = sds1[1])
    epsilon1 <- rnorm(length(PPS1), mean = means1[4], sd = sds1[4])
    v231     <- cov(unlist(trsftable[x, ]), PPS1)
    b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds1[3]^2)) else
           if (omics == "Proteomics")  sqrt(1.01 / (99 * sds1[3]^2)) else
                                       sqrt(3.03 / (7  * sds1[3]^2))
    L1     <- b1 * x11 + b2 * unlist(trsftable[x, ]) + b3 * PPS1 + epsilon1
    thres1 <- qnorm(1 - rate,
                    mean = b1*means1[1] + b2*means1[2] + b3*means1[3],
                    sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                (b3*sds1[3])^2 + 2*b2*b3*v231 + 1))
  } else {
    otu_fg0 <- trsftable[k, ]
    otu_fg  <- matrix(ncol = ncol(trsftable))
    for (ii in seq_along(k))
      otu_fg <- rbind(otu_fg,
        lm(t(otu_fg0)[, ii] ~ trsftable[x, ])$residuals)
    otu_fg   <- otu_fg[-1, ]
    means1   <- c(1, mean(unlist(trsftable[x, ])),
                  apply(otu_fg, 1, mean), 0)
    sds1     <- c(1, sd(unlist(trsftable[x, ])),
                  apply(otu_fg, 1, sd),   1)
    x11      <- rnorm(ncol(trsftable), mean = means1[1], sd = sds1[1])
    epsilon1 <- rnorm(ncol(trsftable), mean = means1[length(means1)],
                      sd = sds1[length(sds1)])
    vresid1  <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_fg[y[1], ], otu_fg[y[2], ])))
    TT1  <- sum(apply(otu_fg, 1, var)) + 2 * vresid1
    v231 <- cov(unlist(trsftable[x, ]), apply(otu_fg, 2, sum))
    b3   <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT1)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT1)) else
                                         sqrt(3.03 / (7  * TT1))
    L1     <- b1 * x11 + b2 * unlist(trsftable[x, ]) +
               b3 * apply(otu_fg, 2, sum) + epsilon1
    thres1 <- qnorm(1 - rate,
                    mean = b1*means1[1] + b2*means1[2] +
                           b3*sum(apply(otu_fg, 1, mean)),
                    sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                (b3^2)*TT1 + 2*b2*b3*v231 + 1))
  }

  case1      <- ifelse(L1 >= thres1, "1", "0")
  entiredata <- data.frame(
    ID      = names(tmr$m1$distance),
    case    = as.numeric(case1),
    x1      = x11,
    PS      = tmr$m1$distance,
    hetero1 = indicator[match(names(tmr$m1$distance), rownames(indicator)), 1],
    hetero2 = indicator[match(names(tmr$m1$distance), rownames(indicator)), 2],
    marker  = unlist(trsftable[x, ])
  )

  # mean reuse count across both rounds (controls used more than once)
  control_freq <- table(pairs$control)
  reused       <- control_freq[control_freq > 1]
  meanrep      <- if (length(reused) > 0) mean(reused) else 1

  append(listtable, list(entiredata, meanrep))
}


# ═══════════════════════════════════════════════════════════
#  tablemaking_t1_wo   (WITHOUT replacement – unchanged logic,
#                       only scaling bug fixed)
# ═══════════════════════════════════════════════════════════
tablemaking_t1_wo <- function(omics, trsftable, datatable, indicator,
                               rate, b, x, cl, scenario,
                               k = NULL, nonlin = "inv") {

  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))

  matmat  <- MatchIt::matchit(matchform, formatch,
                               method  = "nearest",
                               replace = FALSE,
                               caliper = cl)
  samdat0 <- data.frame(
    ps        = matmat$distance,
    pair_taxa = matmat$subclass,
    variable  = taxa_class[match(names(matmat$subclass),
                                 rownames(taxa_class)), ],
    value     = names(matmat$subclass)
  )
  samdat <- samdat0[!is.na(samdat0$pair_taxa), ]

  trsftable_matched <- trsftable[, match(samdat[, 4], colnames(trsftable))]

  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS <- if (nonlin == "inv") 1 / samdat$ps else
        log(exp(samdat$ps) + log(samdat$ps) *
            ((samdat$ps)^-2 - ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS <- indicator[match(samdat[, 4], rownames(indicator)), 1] +
             indicator[match(samdat[, 4], rownames(indicator)), 2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x, ])), mean(PPS), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),   sd(PPS),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[4], sd = sds[4])
    v23     <- cov(unlist(trsftable_matched[x, ]), PPS)
  } else {
    otu_for_generate0 <- trsftable[k, match(samdat[, 4], colnames(trsftable))]
    otu_for_generate  <- matrix(ncol = ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_for_generate <- rbind(otu_for_generate,
        lm(t(otu_for_generate0)[, ii] ~ trsftable_matched[x, ])$residuals)
    otu_for_generate <- otu_for_generate[-1, ]
    means <- c(1, mean(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, mean), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, sd),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[length(means)],
                     sd = sds[length(sds)])
    vresid <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_for_generate[y[1], ], otu_for_generate[y[2], ])))
    TT  <- sum(apply(otu_for_generate, 1, var)) + 2 * vresid
    v23 <- cov(unlist(trsftable_matched[x, ]),
               apply(otu_for_generate, 2, sum))
  }

  b1 <- 0.1; b2 <- 0

  listtable <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds[3]^2)) else
                                         sqrt(3.03 / (7  * sds[3]^2))
      L     <- b1 * x1 + b2 * unlist(trsftable_matched[x, ]) + b3 * PPS + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] + b3*means[3],
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3*sds[3])^2 + 2*b2*b3*v23 + 1))
    } else {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                         sqrt(3.03 / (7  * TT))
      L     <- b1 * x1 + b2 * unlist(trsftable_matched[x, ]) +
               b3 * apply(otu_for_generate, 2, sum) + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] +
                            b3 * sum(apply(otu_for_generate, 1, mean)),
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3^2)*TT + 2*b2*b3*v23 + 1))
    }

    case    <- ifelse(L >= thres, "1", "0")
    samdat1 <- data.frame(
      ID        = samdat[, 4],
      case      = as.numeric(case[match(samdat[, 4], names(case))]),
      pair_taxa = samdat$pair_taxa,
      indOTU    = ifelse(samdat[, 3] == "up", "large", "small"),
      x1        = x1,
      PS        = samdat$ps
    )
    samdat1 <- samdat1[!duplicated(samdat1$ID), ]
    rownames(samdat1) <- samdat1$ID

    trsftable1 <- trsftable[, match(rownames(samdat1), colnames(trsftable))]
    matchdata  <- data.frame(
      samdat1,
      hetero1 = indicator[match(rownames(samdat1), rownames(indicator)), 1],
      hetero2 = indicator[match(rownames(samdat1), rownames(indicator)), 2],
      marker  = unlist(trsftable1[x, ])
    )
    matchdata1 <- matchdata[order(matchdata$pair_taxa), ]
    matchdata1$pair_taxa <- as.factor(matchdata1$pair_taxa)

    # ── weights (mean-normalised, bug fixed) ──
    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m  <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair - 1
    matchdata1$mm  <- matchdata1$m / mean(matchdata1$m)  # FIX: was /sum
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),
                                           names(matmat$weights))]
    matchdata1$.neff <- compute_neff(matchdata1$mm)

    return(matchdata1)
  })

  # entire data
  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS1 <- if (nonlin == "inv") 1 / matmat$distance else
        log(exp(matmat$distance) +
            log(matmat$distance) *
            ((matmat$distance)^-2 -
             ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS1 <- indicator[, 1] + indicator[, 2]
    }
    means1   <- c(1, mean(unlist(trsftable[x, ])), mean(PPS1), 0)
    sds1     <- c(1, sd(unlist(trsftable[x, ])),   sd(PPS1),   1)
    x11      <- rnorm(length(PPS1), mean = means1[1], sd = sds1[1])
    epsilon1 <- rnorm(length(PPS1), mean = means1[4], sd = sds1[4])
    v231     <- cov(unlist(trsftable[x, ]), PPS1)
    b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds1[3]^2)) else
           if (omics == "Proteomics")  sqrt(1.01 / (99 * sds1[3]^2)) else
                                       sqrt(3.03 / (7  * sds1[3]^2))
    L1     <- b1 * x11 + b2 * unlist(trsftable[x, ]) + b3 * PPS1 + epsilon1
    thres1 <- qnorm(1 - rate,
                    mean = b1*means1[1] + b2*means1[2] + b3*means1[3],
                    sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                (b3*sds1[3])^2 + 2*b2*b3*v231 + 1))
  } else {
    otu_fg0  <- trsftable[k, ]
    otu_fg   <- matrix(ncol = ncol(trsftable))
    for (ii in seq_along(k))
      otu_fg <- rbind(otu_fg,
        lm(t(otu_fg0)[, ii] ~ trsftable[x, ])$residuals)
    otu_fg   <- otu_fg[-1, ]
    means1   <- c(1, mean(unlist(trsftable[x, ])),
                  apply(otu_fg, 1, mean), 0)
    sds1     <- c(1, sd(unlist(trsftable[x, ])),
                  apply(otu_fg, 1, sd),   1)
    x11      <- rnorm(ncol(trsftable), mean = means1[1], sd = sds1[1])
    epsilon1 <- rnorm(ncol(trsftable), mean = means1[length(means1)],
                      sd = sds1[length(sds1)])
    vresid1  <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_fg[y[1], ], otu_fg[y[2], ])))
    TT1  <- sum(apply(otu_fg, 1, var)) + 2 * vresid1
    v231 <- cov(unlist(trsftable[x, ]), apply(otu_fg, 2, sum))
    b3   <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT1)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT1)) else
                                         sqrt(3.03 / (7  * TT1))
    L1     <- b1 * x11 + b2 * unlist(trsftable[x, ]) +
               b3 * apply(otu_fg, 2, sum) + epsilon1
    thres1 <- qnorm(1 - rate,
                    mean = b1*means1[1] + b2*means1[2] +
                           b3 * sum(apply(otu_fg, 1, mean)),
                    sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                (b3^2)*TT1 + 2*b2*b3*v231 + 1))
  }

  case1      <- ifelse(L1 >= thres1, "1", "0")
  entiredata <- data.frame(
    ID      = names(matmat$distance),
    case    = as.numeric(case1),
    x1      = x11,
    PS      = matmat$distance,
    hetero1 = indicator[match(names(matmat$distance), rownames(indicator)), 1],
    hetero2 = indicator[match(names(matmat$distance), rownames(indicator)), 2],
    marker  = unlist(trsftable[x, ])
  )

  append(listtable, list(entiredata))
}


# ═══════════════════════════════════════════════════════════
#  tablemaking_power  (WITH replacement – two-round matching)
# ═══════════════════════════════════════════════════════════
tablemaking_power <- function(omics, trsftable, datatable, indicator,
                               rate, b, x, cl, scenario,
                               k = NULL, nonlin = "inv") {

  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))

  tmr    <- two_round_match(formatch, matchform, caliper = cl)
  pairs  <- rbind(tmr$pairs1,
                  if (!is.null(tmr$pairs2)) tmr$pairs2 else NULL)
  if (is.null(pairs) || nrow(pairs) == 0) return(NULL)

  all_ids    <- unique(c(pairs$treated, pairs$control))
  clusters   <- build_clusters_from_pairs(pairs)
  w_vec      <- build_weights(all_ids, pairs)
  round2_ids <- if (!is.null(tmr$pairs2))
    unique(c(tmr$pairs2$treated, tmr$pairs2$control)) else character(0)
  ps_vec     <- tmr$m1$distance
  subclass   <- clusters[all_ids]

  samdat <- data.frame(
    ps        = ps_vec[match(all_ids, names(ps_vec))],
    pair_taxa = subclass,
    variable  = taxa_class[match(all_ids, rownames(taxa_class)), ],
    value     = all_ids,
    round2    = as.integer(all_ids %in% round2_ids),
    stringsAsFactors = FALSE
  )

  trsftable_matched <- trsftable[, match(samdat[, 4], colnames(trsftable))]

  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS <- if (nonlin == "inv") 1 / samdat$ps else
        log(exp(samdat$ps) + log(samdat$ps) *
            ((samdat$ps)^-2 - ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS <- indicator[match(samdat[, 4], rownames(indicator)), 1] +
             indicator[match(samdat[, 4], rownames(indicator)), 2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x, ])), mean(PPS), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),   sd(PPS),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[4], sd = sds[4])
    v23     <- cov(unlist(trsftable_matched[x, ]), PPS)
  } else {
    otu_for_generate0 <- trsftable[k, match(samdat[, 4], colnames(trsftable))]
    otu_for_generate  <- matrix(ncol = ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_for_generate <- rbind(otu_for_generate,
        lm(t(otu_for_generate0)[, ii] ~ trsftable_matched[x, ])$residuals)
    otu_for_generate <- otu_for_generate[-1, ]
    means <- c(1, mean(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, mean), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, sd),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[length(means)],
                     sd = sds[length(sds)])
    vresid <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_for_generate[y[1], ], otu_for_generate[y[2], ])))
    TT  <- sum(apply(otu_for_generate, 1, var)) + 2 * vresid
    v23 <- cov(unlist(trsftable_matched[x, ]),
               apply(otu_for_generate, 2, sum))
  }

  b1 <- 0.1

  listtable <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds[3]^2)) else
                                         sqrt(3.03 / (7  * sds[3]^2))
      b2    <- Re(polyroot(c(-1.01*B + B*(b3*sds[3])^2,
                             -2*B*b3*v23,
                             (1-B)*(sds[2]^2))))[1]
      L     <- b1*x1 + b2*unlist(trsftable_matched[x, ]) + b3*PPS + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] + b3*means[3],
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3*sds[3])^2 + 2*b2*b3*v23 + 1))
    } else {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                         sqrt(3.03 / (7  * TT))
      b2    <- Re(polyroot(c(-1.01*B + B*TT*(b3^2),
                             -2*B*b3*v23,
                             (1-B)*(sds[2]^2))))[1]
      L     <- b1*x1 + b2*unlist(trsftable_matched[x, ]) +
               b3*apply(otu_for_generate, 2, sum) + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] +
                            b3*sum(apply(otu_for_generate, 1, mean)),
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3^2)*TT + 2*b2*b3*v23 + 1))
    }

    case    <- ifelse(L >= thres, "1", "0")
    samdat1 <- data.frame(
      ID        = samdat[, 4],
      case      = as.numeric(case[match(samdat[, 4], names(case))]),
      pair_taxa = samdat$pair_taxa,
      indOTU    = ifelse(samdat[, 3] == "up", "large", "small"),
      x1        = x1,
      PS        = samdat$ps,
      round2    = samdat$round2
    )
    samdat1 <- samdat1[!duplicated(samdat1$ID), ]
    rownames(samdat1) <- samdat1$ID

    trsftable1 <- trsftable[, match(rownames(samdat1), colnames(trsftable))]
    matchdata  <- data.frame(
      samdat1,
      hetero1 = indicator[match(rownames(samdat1), rownames(indicator)), 1],
      hetero2 = indicator[match(rownames(samdat1), rownames(indicator)), 2],
      marker  = unlist(trsftable1[x, ])
    )
    matchdata1 <- matchdata[order(matchdata$pair_taxa), ]
    matchdata1$pair_taxa <- as.factor(matchdata1$pair_taxa)

    w_sub <- w_vec[match(rownames(matchdata1), names(w_vec))]
    w_sub[is.na(w_sub)] <- 1
    matchdata1$mm  <- w_sub / mean(w_sub)
    matchdata1$mmm <- tmr$m1$weights[match(rownames(matchdata1),
                                           names(tmr$m1$weights))]
    matchdata1$.neff <- compute_neff(matchdata1$mm)

    return(matchdata1)
  })

  rm(means, sds, x1, epsilon, v23)

  # entire data (per effect size level)
  listtable1 <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      if (scenario == "S1") {
        PPS1 <- if (nonlin == "inv") 1 / tmr$m1$distance else
          log(exp(tmr$m1$distance) +
              log(tmr$m1$distance) *
              ((tmr$m1$distance)^-2 -
               ifelse(omics == "Metagenomics", 6, 8)) / -2)
      } else {
        PPS1 <- indicator[, 1] + indicator[, 2]
      }
      means1   <- c(1, mean(unlist(trsftable[x, ])), mean(PPS1), 0)
      sds1     <- c(1, sd(unlist(trsftable[x, ])),   sd(PPS1),   1)
      x11      <- rnorm(length(PPS1), mean = means1[1], sd = sds1[1])
      epsilon1 <- rnorm(length(PPS1), mean = means1[4], sd = sds1[4])
      v231     <- cov(unlist(trsftable[x, ]), PPS1)
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds1[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds1[3]^2)) else
                                         sqrt(3.03 / (7  * sds1[3]^2))
      b2    <- Re(polyroot(c(-1.01*B + B*(b3*sds1[3])^2,
                             -2*B*b3*v231,
                             (1-B)*(sds1[2]^2))))[1]
      L1     <- b1*x11 + b2*unlist(trsftable[x, ]) + b3*PPS1 + epsilon1
      thres1 <- qnorm(1 - rate,
                      mean = b1*means1[1] + b2*means1[2] + b3*means1[3],
                      sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                  (b3*sds1[3])^2 + 2*b2*b3*v231 + 1))
    } else {
      otu_fg0  <- trsftable[k, ]
      otu_fg   <- matrix(ncol = ncol(trsftable))
      for (ii in seq_along(k))
        otu_fg <- rbind(otu_fg,
          lm(t(otu_fg0)[, ii] ~ trsftable[x, ])$residuals)
      otu_fg   <- otu_fg[-1, ]
      means1   <- c(1, mean(unlist(trsftable[x, ])),
                    apply(otu_fg, 1, mean), 0)
      sds1     <- c(1, sd(unlist(trsftable[x, ])),
                    apply(otu_fg, 1, sd),   1)
      x11      <- rnorm(ncol(trsftable), mean = means1[1], sd = sds1[1])
      epsilon1 <- rnorm(ncol(trsftable), mean = means1[length(means1)],
                        sd = sds1[length(sds1)])
      vresid1  <- do.call(sum, lapply(
        data.frame(combn(seq_along(k), 2)),
        function(y) cov(otu_fg[y[1], ], otu_fg[y[2], ])))
      TT1  <- sum(apply(otu_fg, 1, var)) + 2 * vresid1
      v231 <- cov(unlist(trsftable[x, ]), apply(otu_fg, 2, sum))
      b3   <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT1)) else
               if (omics == "Proteomics")  sqrt(1.01 / (99 * TT1)) else
                                           sqrt(3.03 / (7  * TT1))
      b2    <- Re(polyroot(c(-1.01*B + B*TT1*(b3^2),
                             -2*B*b3*v231,
                             (1-B)*(sds1[2]^2))))[1]
      L1     <- b1*x11 + b2*unlist(trsftable[x, ]) +
                 b3*apply(otu_fg, 2, sum) + epsilon1
      thres1 <- qnorm(1 - rate,
                      mean = b1*means1[1] + b2*means1[2] +
                             b3*sum(apply(otu_fg, 1, mean)),
                      sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                  (b3^2)*TT1 + 2*b2*b3*v231 + 1))
    }
    case1 <- ifelse(L1 >= thres1, "1", "0")
    data.frame(
      ID      = names(tmr$m1$distance),
      case    = as.numeric(case1),
      x1      = x11,
      PS      = tmr$m1$distance,
      hetero1 = indicator[match(names(tmr$m1$distance), rownames(indicator)), 1],
      hetero2 = indicator[match(names(tmr$m1$distance), rownames(indicator)), 2],
      marker  = unlist(trsftable[x, ])
    )
  })

  list(listtable, listtable1)
}


# ═══════════════════════════════════════════════════════════
#  tablemaking_power_wo  (WITHOUT replacement – bug fix only)
# ═══════════════════════════════════════════════════════════
tablemaking_power_wo <- function(omics, trsftable, datatable, indicator,
                                  rate, b, x, cl, scenario,
                                  k = NULL, nonlin = "inv") {

  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))

  matmat  <- MatchIt::matchit(matchform, formatch,
                               method  = "nearest",
                               replace = FALSE,
                               caliper = cl)
  samdat0 <- data.frame(
    ps        = matmat$distance,
    pair_taxa = matmat$subclass,
    variable  = taxa_class[match(names(matmat$subclass),
                                 rownames(taxa_class)), ],
    value     = names(matmat$subclass)
  )
  samdat <- samdat0[!is.na(samdat0$pair_taxa), ]

  trsftable_matched <- trsftable[, match(samdat[, 4], colnames(trsftable))]

  if (scenario == "S1" | scenario == "S2") {
    if (scenario == "S1") {
      PPS <- if (nonlin == "inv") 1 / samdat$ps else
        log(exp(samdat$ps) + log(samdat$ps) *
            ((samdat$ps)^-2 - ifelse(omics == "Metagenomics", 6, 8)) / -2)
    } else {
      PPS <- indicator[match(samdat[, 4], rownames(indicator)), 1] +
             indicator[match(samdat[, 4], rownames(indicator)), 2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x, ])), mean(PPS), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),   sd(PPS),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[4], sd = sds[4])
    v23     <- cov(unlist(trsftable_matched[x, ]), PPS)
  } else {
    otu_for_generate0 <- trsftable[k, match(samdat[, 4], colnames(trsftable))]
    otu_for_generate  <- matrix(ncol = ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_for_generate <- rbind(otu_for_generate,
        lm(t(otu_for_generate0)[, ii] ~ trsftable_matched[x, ])$residuals)
    otu_for_generate <- otu_for_generate[-1, ]
    means <- c(1, mean(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, mean), 0)
    sds   <- c(1, sd(unlist(trsftable_matched[x, ])),
               apply(otu_for_generate, 1, sd),   1)
    x1      <- rnorm(nrow(samdat), mean = means[1], sd = sds[1])
    epsilon <- rnorm(nrow(samdat), mean = means[length(means)],
                     sd = sds[length(sds)])
    vresid <- do.call(sum, lapply(
      data.frame(combn(seq_along(k), 2)),
      function(y) cov(otu_for_generate[y[1], ], otu_for_generate[y[2], ])))
    TT  <- sum(apply(otu_for_generate, 1, var)) + 2 * vresid
    v23 <- cov(unlist(trsftable_matched[x, ]),
               apply(otu_for_generate, 2, sum))
  }

  b1 <- 0.1

  listtable <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds[3]^2)) else
                                         sqrt(3.03 / (7  * sds[3]^2))
      b2    <- Re(polyroot(c(-1.01*B + B*(b3*sds[3])^2,
                             -2*B*b3*v23,
                             (1-B)*(sds[2]^2))))[1]
      L     <- b1*x1 + b2*unlist(trsftable_matched[x, ]) + b3*PPS + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] + b3*means[3],
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3*sds[3])^2 + 2*b2*b3*v23 + 1))
    } else {
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                         sqrt(3.03 / (7  * TT))
      b2    <- Re(polyroot(c(-1.01*B + B*TT*(b3^2),
                             -2*B*b3*v23,
                             (1-B)*(sds[2]^2))))[1]
      L     <- b1*x1 + b2*unlist(trsftable_matched[x, ]) +
               b3*apply(otu_for_generate, 2, sum) + epsilon
      thres <- qnorm(1 - rate,
                     mean = b1*means[1] + b2*means[2] +
                            b3*sum(apply(otu_for_generate, 1, mean)),
                     sd   = sqrt((b1*sds[1])^2 + (b2*sds[2])^2 +
                                 (b3^2)*TT + 2*b2*b3*v23 + 1))
    }

    case    <- ifelse(L >= thres, "1", "0")
    samdat1 <- data.frame(
      ID        = samdat[, 4],
      case      = as.numeric(case[match(samdat[, 4], names(case))]),
      pair_taxa = samdat$pair_taxa,
      indOTU    = ifelse(samdat[, 3] == "up", "large", "small"),
      x1        = x1,
      PS        = samdat$ps
    )
    samdat1 <- samdat1[!duplicated(samdat1$ID), ]
    rownames(samdat1) <- samdat1$ID

    trsftable1 <- trsftable[, match(rownames(samdat1), colnames(trsftable))]
    matchdata  <- data.frame(
      samdat1,
      hetero1 = indicator[match(rownames(samdat1), rownames(indicator)), 1],
      hetero2 = indicator[match(rownames(samdat1), rownames(indicator)), 2],
      marker  = unlist(trsftable1[x, ])
    )
    matchdata1 <- matchdata[order(matchdata$pair_taxa), ]
    matchdata1$pair_taxa <- as.factor(matchdata1$pair_taxa)

    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m   <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair - 1
    matchdata1$mm  <- matchdata1$m / mean(matchdata1$m)  # FIX
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),
                                           names(matmat$weights))]
    matchdata1$.neff <- compute_neff(matchdata1$mm)

    return(matchdata1)
  })

  listtable1 <- lapply(seq_along(b), function(blen) {
    B <- b[[blen]]
    if (scenario == "S1" | scenario == "S2") {
      if (scenario == "S1") {
        PPS1 <- if (nonlin == "inv") 1 / matmat$distance else
          log(exp(matmat$distance) +
              log(matmat$distance) *
              ((matmat$distance)^-2 -
               ifelse(omics == "Metagenomics", 6, 8)) / -2)
      } else {
        PPS1 <- indicator[, 1] + indicator[, 2]
      }
      means1   <- c(1, mean(unlist(trsftable[x, ])), mean(PPS1), 0)
      sds1     <- c(1, sd(unlist(trsftable[x, ])),   sd(PPS1),   1)
      x11      <- rnorm(length(PPS1), mean = means1[1], sd = sds1[1])
      epsilon1 <- rnorm(length(PPS1), mean = means1[4], sd = sds1[4])
      v231     <- cov(unlist(trsftable[x, ]), PPS1)
      b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * sds1[3]^2)) else
             if (omics == "Proteomics")  sqrt(1.01 / (99 * sds1[3]^2)) else
                                         sqrt(3.03 / (7  * sds1[3]^2))
      b2    <- Re(polyroot(c(-1.01*B + B*(b3*sds1[3])^2,
                             -2*B*b3*v231,
                             (1-B)*(sds1[2]^2))))[1]
      L1     <- b1*x11 + b2*unlist(trsftable[x, ]) + b3*PPS1 + epsilon1
      thres1 <- qnorm(1 - rate,
                      mean = b1*means1[1] + b2*means1[2] + b3*means1[3],
                      sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                  (b3*sds1[3])^2 + 2*b2*b3*v231 + 1))
    } else {
      otu_fg0  <- trsftable[k, ]
      otu_fg   <- matrix(ncol = ncol(trsftable))
      for (ii in seq_along(k))
        otu_fg <- rbind(otu_fg,
          lm(t(otu_fg0)[, ii] ~ trsftable[x, ])$residuals)
      otu_fg   <- otu_fg[-1, ]
      means1   <- c(1, mean(unlist(trsftable[x, ])),
                    apply(otu_fg, 1, mean), 0)
      sds1     <- c(1, sd(unlist(trsftable[x, ])),
                    apply(otu_fg, 1, sd),   1)
      x11      <- rnorm(ncol(trsftable), mean = means1[1], sd = sds1[1])
      epsilon1 <- rnorm(ncol(trsftable), mean = means1[length(means1)],
                        sd = sds1[length(sds1)])
      vresid1  <- do.call(sum, lapply(
        data.frame(combn(seq_along(k), 2)),
        function(y) cov(otu_fg[y[1], ], otu_fg[y[2], ])))
      TT1  <- sum(apply(otu_fg, 1, var)) + 2 * vresid1
      v231 <- cov(unlist(trsftable[x, ]), apply(otu_fg, 2, sum))
      b3   <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT1)) else
               if (omics == "Proteomics")  sqrt(1.01 / (99 * TT1)) else
                                           sqrt(3.03 / (7  * TT1))
      b2    <- Re(polyroot(c(-1.01*B + B*TT1*(b3^2),
                             -2*B*b3*v231,
                             (1-B)*(sds1[2]^2))))[1]
      L1     <- b1*x11 + b2*unlist(trsftable[x, ]) +
                 b3*apply(otu_fg, 2, sum) + epsilon1
      thres1 <- qnorm(1 - rate,
                      mean = b1*means1[1] + b2*means1[2] +
                             b3*sum(apply(otu_fg, 1, mean)),
                      sd   = sqrt((b1*sds1[1])^2 + (b2*sds1[2])^2 +
                                  (b3^2)*TT1 + 2*b2*b3*v231 + 1))
    }
    case1 <- ifelse(L1 >= thres1, "1", "0")
    data.frame(
      ID      = names(matmat$distance),
      case    = as.numeric(case1),
      x1      = x11,
      PS      = matmat$distance,
      hetero1 = indicator[match(names(matmat$distance), rownames(indicator)), 1],
      hetero2 = indicator[match(names(matmat$distance), rownames(indicator)), 2],
      marker  = unlist(trsftable[x, ])
    )
  })

  list(listtable, listtable1)
}
