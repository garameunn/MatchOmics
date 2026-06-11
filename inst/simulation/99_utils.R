## ============================================================
##  99_utils.R
##  Changes vs original:
##  1. compute_neff() added   – Kish effective sample size
##  2. love_plot()    added   – Love plot (SMD pre vs post)
##  3. smd_summary()  added   – summarise SMD across taxa/reps
##  4. Bioconductor package check added to options
## ============================================================

options(digits = 3)

# ── column name reference ───────────────────────────────────
cnames <- c("LR_unadj", "LR_PS", "LR_Combat", "LR_Limma",
            "wrep_nocal_ind",   "wrep_nocal_ex",
            "wrep_cal0.3_ind",  "wrep_cal0.3_ex",
            "wrep_cal0.2_ind",  "wrep_cal0.2_ex",
            "wrep_cal0.1_ind",  "wrep_cal0.1_ex",
            "wrep_cal0.05_ind", "wrep_cal0.05_ex",
            "worep_nocal_ind",   "worep_nocal_ex",
            "worep_cal0.3_ind",  "worep_cal0.3_ex",
            "worep_cal0.2_ind",  "worep_cal0.2_ex",
            "worep_cal0.1_ind",  "worep_cal0.1_ex",
            "worep_cal0.05_ind", "worep_cal0.05_ex")

# ── OTU filtering ───────────────────────────────────────────
filtOTU <- function(x) {
  cr1 <- which(apply(x, 1, function(x) sum(x == 0)) / ncol(x) > 0.5)
  cr2 <- which(apply(x, 1, sum) / sum(apply(x, 1, sum)) < 0.001)
  if (length(c(cr1, cr2)) != 0) x[-unique(c(cr1, cr2)), ] else x
}

# ── type 1 error ────────────────────────────────────────────
type1 <- function(x, lev) {
  t001 <- sum(ifelse(is.na(unlist(x)), 1, unlist(x)) < 0.01)  / length(unlist(x))
  t005 <- sum(ifelse(is.na(unlist(x)), 1, unlist(x)) < 0.05)  / length(unlist(x))
  t01  <- sum(ifelse(is.na(unlist(x)), 1, unlist(x)) < 0.1)   / length(unlist(x))
  t02  <- sum(ifelse(is.na(unlist(x)), 1, unlist(x)) < 0.2)   / length(unlist(x))
  round(c(t001, t005, t01, t02), 4)
}

type1_1 <- function(x) {
  v <- ifelse(is.na(x), 1, x)
  round(c(
    sum(v < 0.001) / length(v),
    sum(v < 0.005) / length(v),
    sum(v < 0.01)  / length(v),
    sum(v < 0.05)  / length(v),
    sum(v < 0.1)   / length(v)
  ), 5)
}

# ── confidence interval for proportions ─────────────────────
CIfun <- function(pvals, level) {
  testp <- prop.test(x = sum(pvals < level, na.rm = TRUE),
                     n = length(pvals), correct = TRUE)
  testp$conf.int
}

CIfun1 <- function(pvals, level) {
  if (sum(is.na(pvals)) > 0) pvals[which(is.na(pvals))] <- 1
  xx    <- round((sum(pvals < level) / length(pvals)) * 5000, 0)
  testp <- prop.test(x = xx, n = 5000, correct = TRUE)
  testp$conf.int
}

# ── power helpers ───────────────────────────────────────────
powerfold <- function(p, nmodel) {
  p0001 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[1, ]))),
                  byrow = TRUE, ncol = nmodel)
  p0005 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[2, ]))),
                  byrow = TRUE, ncol = nmodel)
  p001  <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[3, ]))),
                  byrow = TRUE, ncol = nmodel)
  p002  <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[4, ]))),
                  byrow = TRUE, ncol = nmodel)
  list(p0001, p0005, p001, p002)
}

powerfold1 <- function(x) {
  sum(ifelse(is.na(x), 1, x) < 0.05) / length(x)
}

powerbias <- function(pw, nmodel) {
  pvalues_p  <- lapply(pw, function(x) lapply(x, function(xx) xx[, 1:nmodel]))
  pvalues1   <- powerfold(pvalues_p, nmodel)
  powers     <- do.call(rbind, lapply(pvalues1, function(xx) apply(xx, 2, powerfold1)))
  conf.int   <- do.call(rbind, lapply(pvalues1, function(y)
    apply(y, 2, FUN = CIfun, level = 0.05)))
  powers1    <- rbind(powers, conf.int)

  betas_p    <- lapply(pw, function(x) lapply(x, function(xx)
    xx[, ((nmodel + 1):(2 * nmodel))]))
  betas1     <- powerfold(betas_p, nmodel)
  betamean   <- matrix(unlist(lapply(betas1, function(x)
    apply(x, 2, mean, na.rm = TRUE))), ncol = nmodel, byrow = TRUE)

  rownames(powers1)[5:12]  <- paste0("CI_", c(rep(c("var.001","var.002","var.005","var.01"), each = 2)))
  rownames(powers1)[1:4]   <- paste0("Power_", c("var.001","var.002","var.005","var.01"))
  rownames(betamean)        <- paste0("Beta_",  c("var.001","var.002","var.005","var.01"))
  colnames(powers1) <- colnames(betamean) <- cnames

  list(powers1, betamean)
}

# ── relative bias ───────────────────────────────────────────
relbiasfun <- function(dt, ref, map) {
  dt1 <- dt[, -ncol(dt)]
  dtt <- c()
  for (i in seq_len(ncol(dt1))) {
    cols <- (dt1[, i] - ref[, map[i]]) / ref[, map[i]]
    dtt  <- cbind(dtt, cols)
  }
  rownames(dtt) <- paste("RelBias", rownames(dtt), sep = "_")
  colnames(dtt) <- colnames(dt1)
  dtt           <- data.frame(dtt)
  dtt$name      <- dt$name
  dtt
}

# ── nominal b2 ──────────────────────────────────────────────
nominalb2 <- function(blist, diffbig) {
  bblist <- blist[diffbig]
  nms    <- c("whole",
              "wrep_nocal",   "wrep_0.3cal",  "wrep_0.2cal",
              "wrep_0.1cal",  "wrep_0.05cal",
              "worep_nocal",  "worep_0.3cal", "worep_0.2cal",
              "worep_0.1cal", "worep_0.05cal")
  lapply(seq_along(nms), function(j)
    apply(do.call(rbind, lapply(bblist, function(y) y[, j])), 2, mean))
}

# ── misc helpers ─────────────────────────────────────────────
inormal   <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

psdiffun  <- function(otus, k) {
  taxa_class <- ifelse(otus >= median(as.numeric(unlist(otus))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))
  ps       <- MatchIt::matchit(matchform, formatch, method = "nearest")$distance
  otusup   <- mean(ps[which(otus >= median(as.numeric(unlist(otus))))])
  otusdown <- mean(ps[which(otus <  median(as.numeric(unlist(otus))))])
  otusup - otusdown
}

pickp <- function(x, otus) x[otus]

# ── power result organiser ───────────────────────────────────
nsamchar     <- c(1000, 500, 200)
mappingchar  <- c(rep(1, 4), rep(2:11, each = 2))

tablesfun <- function(norder, scene, bin = TRUE, prevalence = 0.2,
                      append = TRUE, nonlin) {
  nsamp <- nsamchar[norder]
  fpath <- if (bin)
    sprintf("../data/Power/NSAMP%d_REP1000_RATE%.1f_%s_%s.RData",
            nsamp, prevalence, nonlin, scene)
  else
    sprintf("../data/Power_conti/NSAMP%d_REP1000_%s_%s.RData",
            nsamp, nonlin, scene)
  load(fpath)
  pvalues_test <- pvalues; rm(pvalues)

  # 1. Power & Bias (all taxa)
  totres       <- data.frame(do.call(rbind, powerbias(pvalues_test, nmodel = 24)))
  totres$name  <- "Overall"

  # 2. Top 20
  psdiff       <- psdifflist[[norder]]
  pvalues_pp   <- lapply(pvalues_test, pickp,
                         otus = order(psdiff, decreasing = TRUE)[1:20])
  t20res       <- data.frame(do.call(rbind, powerbias(pvalues_pp, nmodel = 24)))
  t20res$name  <- "Top20"

  # 3. Beta relative bias
  combi       <- paste(scene, nsamp, sep = "_")
  resbeta     <- totres[13:16, ]
  refbias     <- b2list[[combi]][[1]]
  relb        <- relbiasfun(dt = resbeta,  ref = refbias,  map = mappingchar)
  resbeta20   <- t20res[13:16, ]
  refbias20   <- b2list[[combi]][[2]]
  relb20      <- relbiasfun(dt = resbeta20, ref = refbias20, map = mappingchar)

  tables       <- data.frame(rbind(totres, t20res, relb, relb20))
  globname     <- if (bin)
    paste("Bin", scene, nsamp, prevalence, nonlin, sep = "_")
  else
    paste("Con", scene, nsamp, sep = "_")
  tables$name2 <- globname
  message(globname)

  out_csv <- file.path(dirname(dirname(getwd())), "results", "powerbias.csv")
  # 서버: /home2/nekim/Matching2026/results/powerbias.csv
  # 로컬: 상대경로로 자동 설정
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.table(tables, file = out_csv, append = append,
              col.names = !append, quote = FALSE, sep = ",")
}

# ════════════════════════════════════════════════════════════
#  NEW: Neff (Kish 1965 effective sample size)
# ════════════════════════════════════════════════════════════
#' @param w numeric vector of weights (mean-normalised, sum = N)
#' @return scalar Neff
compute_neff <- function(w) {
  w <- w[!is.na(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  (sum(w))^2 / sum(w^2)
}

# ════════════════════════════════════════════════════════════
#  NEW: SMD summary across taxa and reps
# ════════════════════════════════════════════════════════════
#' Summarise pre/post SMD stored in _SMD.RData files
#' @param smd_list list[[rep]][[taxa]]$pre / $post
#' @return data.frame with mean abs SMD pre and post per covariate
smd_summary <- function(smd_list) {
  # flatten to list of pre and post vectors
  pre_mat  <- do.call(rbind, lapply(smd_list, function(h)
    do.call(rbind, lapply(h, function(x) abs(x$pre)))))
  post_mat <- do.call(rbind, lapply(smd_list, function(h)
    do.call(rbind, lapply(h, function(x) abs(x$post)))))

  data.frame(
    covariate    = colnames(pre_mat),
    mean_smd_pre  = colMeans(pre_mat,  na.rm = TRUE),
    mean_smd_post = colMeans(post_mat, na.rm = TRUE),
    row.names     = NULL
  )
}

# ════════════════════════════════════════════════════════════
#  NEW: Love plot
# ════════════════════════════════════════════════════════════
#' Draw a Love plot from smd_summary() output
#' @param smd_df   output of smd_summary()
#' @param title    plot title
#' @param threshold dashed line at this SMD value (default 0.1)
love_plot <- function(smd_df, title = "Love Plot", threshold = 0.1) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for love_plot()")

  # reshape to long
  df_long <- rbind(
    data.frame(covariate = smd_df$covariate,
               smd       = smd_df$mean_smd_pre,
               timing    = "Before matching"),
    data.frame(covariate = smd_df$covariate,
               smd       = smd_df$mean_smd_post,
               timing    = "After matching")
  )
  df_long$timing <- factor(df_long$timing,
                            levels = c("Before matching", "After matching"))

  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = smd, y = covariate,
                               colour = timing, shape = timing)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = threshold,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c("Before matching" = "#E74C3C",
                 "After matching"  = "#2ECC71")) +
    ggplot2::labs(title  = title,
                  x      = "Absolute Standardised Mean Difference",
                  y      = NULL,
                  colour = NULL, shape = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

# ════════════════════════════════════════════════════════════
#  NEW: Neff summary across scenarios (for table/figure)
# ════════════════════════════════════════════════════════════
#' Extract mean Neff from pvalues_comb
#' @param pvalues_comb matrix with neff_* columns
#' @return data.frame: model x mean_neff
# ════════════════════════════════════════════════════════════
#  CI from beta + SE  (compute later if needed in 06_)
# ════════════════════════════════════════════════════════════
#' Compute confidence intervals from beta and SE matrices
#' @param beta_mat  matrix of beta estimates (rows=obs, cols=models)
#' @param se_mat    matrix of SE estimates   (same dims as beta_mat)
#' @param level     confidence level (default 0.95)
#' @return list(lower, upper) each same dim as input
compute_ci <- function(beta_mat, se_mat, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  list(
    lower = beta_mat - z * se_mat,
    upper = beta_mat + z * se_mat
  )
}

#' Summarise mean CI width per model column
#' @param beta_mat  matrix of beta estimates
#' @param se_mat    matrix of SE estimates
#' @param level     confidence level
#' @return named vector of mean CI widths
ci_width_summary <- function(beta_mat, se_mat, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  round(colMeans(2 * z * se_mat, na.rm = TRUE), 5)
}

# ── Neff summary ─────────────────────────────────────────────
neff_summary <- function(pvalues_comb) {
  neff_cols <- grep("^neff_", colnames(pvalues_comb), value = TRUE)
  data.frame(
    model     = neff_cols,
    mean_neff = colMeans(pvalues_comb[, neff_cols, drop = FALSE], na.rm = TRUE),
    row.names = NULL
  )
}
