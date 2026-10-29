## ============================================================
##  02_generator_v3.R
##  Revision (BIB-26-1896) simulation engine.
##
##  Changes vs 02_simulation_generator.R / _S3fix.R:
##  1. [R3-4] Outcome is generated ONCE on the full sample per
##     (replicate, taxon). Every method (full-sample regressions,
##     weighting, matched GEE) analyses the SAME outcome vector;
##     matched analyses simply subset it. b2/b3 are never
##     recomputed inside a matched subsample.
##  2. [R3-3] h^2 -> b2 quadratic uses the correct constant term
##     -h^2 (1.01 + b3^2 var(conf)); the positive root is taken
##     explicitly (closed form) instead of polyroot()[1].
##  3. Liability threshold = empirical (1 - rate) quantile of L,
##     so each replicate has prevalence = rate (as in Methods).
##  4. Matching depends only on the marker split and heterogeneity
##     indices (not on the outcome), so it is computed once per
##     (N, taxon, caliper) and cached across replicates/scenarios.
##     Without-replacement matching = round 1 of two_round_match
##     (identical matchit call as the original *_wo functions).
##  5. conf_share: proportion of (z1 + confounder + error) liability
##     variance attributed to the confounder. 0.1 reproduces the
##     original Metagenomics setting b3 = sqrt(1.01 / (9 var(conf))).
##  6. S3 confounder = S3fix definition (dynamic low-correlation
##     taxa, mean-centred, no residualisation).
## ============================================================

if (!exists("TEST_MODE")) TEST_MODE <- FALSE

# ── helpers (unchanged from original) ──────────────────────
f <- function(idx) set.seed(seeds[idx])

trnsf <- function(dt, omics, lb.size) {
  d <- if (omics == "Metagenomics") edgeR::cpm(dt, log = TRUE, lib.size = lb.size) else dt
  matrix(as.numeric(d), nrow = nrow(d), ncol = ncol(d),
         dimnames = list(rownames(d), colnames(d)))
}

compute_neff <- function(w) (sum(w))^2 / sum(w^2)

# two-round matching: identical to 02_simulation_generator.R
two_round_match <- function(formatch, matchform, caliper) {
  m1 <- MatchIt::matchit(matchform, data = formatch, method = "nearest",
                         replace = FALSE, caliper = caliper)
  matched1_ids      <- rownames(MatchIt::match.data(m1))
  treated_all       <- rownames(formatch)[formatch$taxa_class == levels(formatch$taxa_class)[2]]
  unmatched_treated <- setdiff(treated_all, matched1_ids)

  pairs1 <- data.frame(treated = rownames(m1$match.matrix),
                       control = m1$match.matrix[, 1],
                       round = 1L, stringsAsFactors = FALSE)
  pairs1 <- pairs1[!is.na(pairs1$control), ]
  result <- list(m1 = m1, m2 = NULL, pairs1 = pairs1, pairs2 = NULL)
  if (length(unmatched_treated) == 0) return(result)

  controls_all <- rownames(formatch)[formatch$taxa_class == levels(formatch$taxa_class)[1]]
  formatch2    <- formatch[c(unmatched_treated, controls_all), ]
  if (sum(formatch2$taxa_class == levels(formatch$taxa_class)[2]) > 0 &&
      sum(formatch2$taxa_class == levels(formatch$taxa_class)[1]) > 0) {
    result$m2 <- tryCatch(
      MatchIt::matchit(matchform, data = formatch2, method = "nearest",
                       replace = FALSE, caliper = caliper),
      error = function(e) NULL)
  }
  if (!is.null(result$m2)) {
    mm2    <- result$m2$match.matrix
    pairs2 <- data.frame(treated = rownames(mm2), control = mm2[, 1],
                         round = 2L, stringsAsFactors = FALSE)
    result$pairs2 <- pairs2[!is.na(pairs2$control), ]
  }
  result
}

build_clusters_from_pairs <- function(pairs_df) {
  edges  <- unique(data.frame(from = pairs_df$treated, to = pairs_df$control,
                              stringsAsFactors = FALSE))
  nodes  <- unique(c(edges$from, edges$to))
  parent <- setNames(as.list(nodes), nodes)
  find_root <- function(x) { while (!identical(parent[[x]], x)) x <- parent[[x]]; x }
  for (i in seq_len(nrow(edges))) {
    ra <- find_root(edges$from[i]); rb <- find_root(edges$to[i])
    if (!identical(ra, rb)) parent[[rb]] <- ra
  }
  roots <- vapply(nodes, find_root, character(1))
  setNames(as.integer(factor(roots)), nodes)
}

# S3fix confounder (dynamic k, mean-centred, no residualisation)
select_k_dynamic <- function(trsftable, x, pct = 0.05) {
  x_vec <- as.numeric(trsftable[x, ])
  cors  <- sapply(seq_len(nrow(trsftable)), function(i)
    if (i == x) NA_real_ else cor(x_vec, as.numeric(trsftable[i, ]), use = "complete.obs"))
  abs_cors <- abs(cors)
  k_idx <- which(!is.na(abs_cors) & abs_cors <= quantile(abs_cors, pct, na.rm = TRUE))
  if (length(k_idx) == 0) k_idx <- which.min(abs_cors)
  k_idx
}

# ═══════════════════════════════════════════════════════════
#  1. Matching cache  (per N, taxon, caliper; outcome-free)
# ═══════════════════════════════════════════════════════════
precompute_matching <- function(datatable, indicator, x, calipers) {
  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                 paste(colnames(indicator), collapse = "+"), sep = "~"))

  cal_res <- lapply(calipers, function(cl) {
    tmr <- tryCatch(two_round_match(formatch, matchform, caliper = cl),
                    error = function(e) NULL)
    if (is.null(tmr)) return(NULL)

    ## with replacement (two-round)
    pairs <- rbind(tmr$pairs1, tmr$pairs2)
    wrep  <- NULL
    if (nrow(pairs) > 0) {
      ids  <- unique(c(pairs$treated, pairs$control))
      freq <- table(c(pairs$treated, pairs$control))
      m    <- as.numeric(freq[ids])
      cf   <- table(pairs$control)
      wrep <- list(ids     = ids,
                   cluster = build_clusters_from_pairs(pairs)[ids],
                   w       = m / mean(m),                     # sum(w) = N
                   meanrep = if (any(cf > 1)) mean(cf[cf > 1]) else 1)
    }

    ## without replacement (= round 1)
    md    <- MatchIt::match.data(tmr$m1)
    worep <- if (nrow(md) > 0)
      list(ids = rownames(md), cluster = setNames(as.integer(md$subclass), rownames(md)),
           w = rep(1, nrow(md)), meanrep = 1) else NULL

    list(wrep = wrep, worep = worep, ps = tmr$m1$distance)
  })

  ps_full <- NULL
  for (cr in cal_res) if (!is.null(cr)) { ps_full <- cr$ps; break }
  if (!is.null(ps_full) && is.null(names(ps_full))) names(ps_full) <- rownames(formatch)
  for (nm in names(cal_res)) if (!is.null(cal_res[[nm]])) cal_res[[nm]]$ps <- NULL

  list(ps         = ps_full,                                   # named by sample id
       treat      = setNames(as.integer(taxa_class == "up"), rownames(formatch)),
       cal        = cal_res)
}

# ═══════════════════════════════════════════════════════════
#  2. Confounder (full sample, aligned to colnames(trsftable))
# ═══════════════════════════════════════════════════════════
make_confounder <- function(scenario, ps, indicator, trsftable, x) {
  ids <- colnames(trsftable)
  if (scenario == "S1") {
    1 / ps[ids]
  } else if (scenario == "S2") {
    setNames(indicator[ids, 1] + indicator[ids, 2], ids)
  } else if (scenario == "S3") {
    k_idx <- select_k_dynamic(trsftable, x)
    k_mat <- trsftable[k_idx, ids, drop = FALSE]
    colSums(k_mat - rowMeans(k_mat))
  } else stop("unknown scenario: ", scenario)
}

# ═══════════════════════════════════════════════════════════
#  3. Effect sizes
# ═══════════════════════════════════════════════════════════
# b3 such that var(b3*conf) = 1.01 * s / (1 - s)
#   s = 0.1 -> b3^2 var(conf) = 1.01/9  (original Metagenomics value)
get_b3 <- function(conf, conf_share = 0.1)
  sqrt(1.01 * conf_share / (1 - conf_share) / var(conf))

# Solve h2 = b2^2 sx2 / (1.01 + b3^2 sc2 + b2^2 sx2 + 2 b2 b3 sxc + ...)
#   (1-h2) sx2 b2^2 - 2 h2 b3 sxc b2 - h2 (1.01 + b3^2 sc2) = 0
# 1.01 = b1^2 var(z1) + var(eps) with b1 = 0.1, var(z1) = var(eps) = 1.
# Positive root taken explicitly.
get_b2 <- function(h2, marker, conf, b3) {
  if (h2 == 0) return(0)
  sx2 <- var(marker); sc2 <- var(conf); sxc <- cov(marker, conf)
  qa <- (1 - h2) * sx2
  qb <- -2 * h2 * b3 * sxc
  qc <- -h2 * (1.01 + b3^2 * sc2)
  (-qb + sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)
}

# realised h2 (for audit tables)
realised_h2 <- function(b2, marker, conf, b3) {
  num <- b2^2 * var(marker)
  num / (1.01 + b3^2 * var(conf) + num + 2 * b2 * b3 * cov(marker, conf))
}

# ═══════════════════════════════════════════════════════════
#  4. Outcome (ONE draw on the full sample)
# ═══════════════════════════════════════════════════════════
generate_case <- function(marker, conf, x1, eps, b2, b3, rate, b1 = 0.1) {
  L <- b1 * x1 + b2 * marker + b3 * conf + eps
  as.integer(L > quantile(L, 1 - rate, names = FALSE))
}

# ═══════════════════════════════════════════════════════════
#  5. Fitters  — each returns c(p, beta, se[, rho])
# ═══════════════════════════════════════════════════════════
NA3 <- c(1, NA_real_, NA_real_)

glm_marker <- function(formula, data) {
  fit <- tryCatch(glm(formula, data = data, family = binomial), error = function(e) NULL)
  if (is.null(fit) || !fit$converged) return(NA3)
  cs <- coef(summary(fit))
  if (!"marker" %in% rownames(cs) && !"cmbt" %in% rownames(cs)) return(NA3)
  r  <- if ("marker" %in% rownames(cs)) "marker" else "cmbt"
  c(cs[r, 4], cs[r, 1], cs[r, 2])
}

# weighted logistic regression with HC0 sandwich SE (IPW / overlap weights)
# PS treated as known (standard, conservative for IPW).
wglm_marker <- function(data, w) {
  fit <- tryCatch(suppressWarnings(
    glm(case ~ x1 + marker, data = data, family = quasibinomial, weights = w)),
    error = function(e) NULL)
  if (is.null(fit) || !fit$converged) return(NA3)
  X  <- model.matrix(fit); mu <- fitted(fit); r <- data$case - mu
  bread <- tryCatch(solve(crossprod(X * (w * mu * (1 - mu)), X)), error = function(e) NULL)
  if (is.null(bread)) return(NA3)
  meat <- crossprod(X * (w * r))
  V    <- bread %*% meat %*% bread
  j    <- which(colnames(X) == "marker")
  b    <- coef(fit)[j]; se <- sqrt(V[j, j])
  c(2 * pnorm(-abs(b / se)), b, se)
}

# matched GEE (same model as original: case ~ x1 + PS + marker)
fit_geem <- function(tble, corstr) {
  if (is.null(tble) || nrow(tble) == 0) return(c(NA3, NA_real_))
  fit <- tryCatch(
    geeM::geem(case ~ x1 + PS + marker, data = tble, family = "binomial",
               corstr = corstr, id = tble$pair_taxa, weights = tble$mm,
               sandwich = TRUE, useP = TRUE, maxit = 25, tol = 1e-3),
    error = function(e) NULL)
  if (is.null(fit)) return(c(NA3, NA_real_))
  ok <- isTRUE(fit$converged) && abs(fit$beta[4]) < 100
  if (!ok) return(c(NA3, NA_real_))
  c(summary(fit)$p[4], fit$beta[4], sqrt(diag(fit$var))[4],
    if (corstr == "exchangeable") fit$alpha else NA_real_)
}

matched_table <- function(full, mres) {
  if (is.null(mres)) return(NULL)
  tb <- full[mres$ids, c("case", "x1", "PS", "marker")]
  tb$pair_taxa <- mres$cluster
  tb$mm        <- mres$w
  tb <- tb[order(tb$pair_taxa), ]
  tb$pair_taxa <- as.factor(tb$pair_taxa)
  tb
}

# ═══════════════════════════════════════════════════════════
#  6. One (replicate, taxon): generate once, analyse with all
# ═══════════════════════════════════════════════════════════
FULL_METHODS <- c("LR_unadj", "LR_PS", "LR_Combat", "LR_PSspline", "IPW", "OW")

analyse_one <- function(x, mcache, trsftable, indicator, cmbt_row,
                        scenario, rate, b_levels, conf_share, cal_names) {
  ids    <- colnames(trsftable)
  marker <- setNames(as.numeric(trsftable[x, ]), ids)
  ps     <- mcache$ps[ids]
  conf   <- make_confounder(scenario, mcache$ps, indicator, trsftable, x)
  b3     <- get_b3(conf, conf_share)

  # noise drawn once per (rep, taxon), shared across effect sizes
  n   <- length(ids)
  x1  <- rnorm(n, 1, 1)
  eps <- rnorm(n, 0, 1)

  trt <- mcache$treat[ids]
  w_ipw <- ifelse(trt == 1, 1 / ps, 1 / (1 - ps))
  w_ow  <- ifelse(trt == 1, 1 - ps, ps)
  lps   <- qlogis(pmin(pmax(ps, 1e-6), 1 - 1e-6))
  spl   <- splines::ns(lps, df = 4)

  lapply(b_levels, function(B) {
    b2   <- get_b2(B, marker, conf, b3)
    case <- generate_case(marker, conf, x1, eps, b2, b3, rate)
    full <- data.frame(case = case, x1 = x1, PS = ps, marker = marker,
                       hetero1 = indicator[ids, 1], hetero2 = indicator[ids, 2],
                       cmbt = cmbt_row, row.names = ids)

    res_full <- list(
      LR_unadj    = glm_marker(case ~ x1 + marker, full),
      LR_PS       = glm_marker(case ~ x1 + PS + marker, full),
      LR_Combat   = glm_marker(case ~ x1 + hetero1 + hetero2 + cmbt, full),
      LR_PSspline = glm_marker(case ~ x1 + spl + marker, full),   # spl from enclosing env
      IPW         = wglm_marker(full, w_ipw),
      OW          = wglm_marker(full, w_ow))

    res_m <- list()
    dg    <- c()
    for (cn in cal_names) {
      cr <- mcache$cal[[cn]]
      for (type in c("wrep", "worep")) {
        mres <- if (is.null(cr)) NULL else cr[[type]]
        tb   <- matched_table(full, mres)
        res_m[[paste0(type, "_", cn, "_ind")]] <- fit_geem(tb, "independence")
        res_m[[paste0(type, "_", cn, "_ex")]]  <- fit_geem(tb, "exchangeable")
        dg[paste0("nrow_",  type, "_", cn)] <- if (is.null(tb)) 0 else nrow(tb)
        dg[paste0("neff_",  type, "_", cn)] <- if (is.null(tb)) NA else compute_neff(tb$mm)
        dg[paste0("csize_", type, "_", cn)] <- if (is.null(tb)) NA else mean(table(tb$pair_taxa))
        if (type == "wrep")
          dg[paste0("meanrep_wrep_", cn)] <- if (is.null(mres)) NA else mres$meanrep
      }
    }

    all_res <- c(res_full, res_m)
    out <- c(setNames(sapply(all_res, `[`, 1), paste0("sig_",  names(all_res))),
             setNames(sapply(all_res, `[`, 2), paste0("coef_", names(all_res))),
             setNames(sapply(all_res, `[`, 3), paste0("se_",   names(all_res))),
             setNames(sapply(res_m[grep("_ex$", names(res_m))], `[`, 4),
                      paste0("rho_", grep("_ex$", names(res_m), value = TRUE))),
             dg,
             b2_true = b2, b3_true = b3,
             h2_real = realised_h2(b2, marker, conf, b3),
             prev    = mean(case))
    out
  })
}
