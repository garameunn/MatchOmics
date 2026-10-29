## ============================================================
##  11_summarize_v3b.R — summarise v3b (DR weighting) next to v3
##  Usage: Rscript 11_summarize_v3b.R            (from code/)
##  Needs: results/v3/{Power,T1} and results/v3b/{Power,T1}
##  1) Same-data check: LR_PS / LR_PSspline p-values in v3b must
##     equal those in v3 (max |diff| printed per file).
##  2) RawPw / T1 / AdjPw for IPW_DR, OW_DR, and ratios
##     wrep_cal0.3_ind / {IPW_DR, OW_DR, LR_PSspline, LR_PS}.
##  Output: ../results/v3b/summary_v3b_long.csv, summary_v3b_key_{val}.csv
##  (definitions identical to 10_summarize_v3.R)
## ============================================================
base_dir <- "/path/to/Matching2026"
v3_dir   <- file.path(base_dir, "results", "v3")
v3b_dir  <- file.path(base_dir, "results", "v3b")
load(file.path(base_dir, "data", "psdifflist.RData"))
norder <- c("1000" = 1, "500" = 2, "200" = 3)

pick_file <- function(dir, N, rate, sc) {
  f <- list.files(dir, sprintf("^NSAMP%d_REP[0-9]+_RATE%.1f_%s\\.RData$", N, rate, sc), full.names = TRUE)
  if (length(f) == 0) return(NULL)
  reps <- as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(f)))
  f[which.max(reps)]
}
load_mat <- function(f) { e <- new.env(); load(f, envir = e); e$pvalues[[length(e$pvalues)]] }

same_data <- function(a, b, f) {
  if (nrow(a) != nrow(b)) { message("  [CHECK FAIL] row count differs: ", basename(f)); return(Inf) }
  d <- max(abs(a[, c("sig_LR_PS", "sig_LR_PSspline")] - b[, c("sig_LR_PS", "sig_LR_PSspline")]),
           na.rm = TRUE)
  message(sprintf("  [CHECK] %-35s max|dp| = %.2e", basename(f), d))
  d
}

meths <- c("LR_PS", "LR_PSspline", "IPW_DR", "OW_DR", "wrep_cal0.3_ind", "wrep_nocal_ind")
rows <- list(); checks <- c()
for (N in c(1000, 500, 200)) for (sc in c("S1", "S2", "S3")) for (rate in c(0.1, 0.2, 0.3)) {
  fp  <- pick_file(file.path(v3_dir,  "Power"), N, rate, sc)
  fpb <- pick_file(file.path(v3b_dir, "Power"), N, rate, sc)
  if (is.null(fp) || is.null(fpb)) next
  ft  <- pick_file(file.path(v3_dir,  "T1"), N, rate, sc)
  ftb <- pick_file(file.path(v3b_dir, "T1"), N, rate, sc)

  pw <- load_mat(fp); pwb <- load_mat(fpb)
  checks[basename(fpb)] <- same_data(pw, pwb, fpb)
  t1 <- t1b <- NULL
  if (!is.null(ft) && !is.null(ftb)) {
    t1 <- load_mat(ft); t1b <- load_mat(ftb)
    checks[paste0("T1_", basename(ftb))] <- same_data(t1, t1b, ftb)
  }

  psd  <- psdifflist[[norder[as.character(N)]]]
  dimm <- length(psd)
  taxon_of <- function(m) (seq_len(nrow(m)) - 1) %% dimm + 1
  top20 <- order(psd, decreasing = TRUE)[1:20]

  for (meth in meths) {
    src  <- if (meth %in% c("IPW_DR", "OW_DR")) pwb else pw
    src0 <- if (meth %in% c("IPW_DR", "OW_DR")) t1b else t1
    s <- paste0("sig_", meth)
    p  <- src[, s]; p[is.na(p)] <- 1
    p0 <- if (is.null(src0)) NULL else { x <- src0[, s]; x[is.na(x)] <- 1; x }
    for (ss in c("all", "top20")) {
      keep  <- if (ss == "all") rep(TRUE, length(p)) else taxon_of(src) %in% top20
      keep0 <- if (is.null(p0)) NULL else if (ss == "all") rep(TRUE, length(p0)) else taxon_of(src0) %in% top20
      thr   <- if (is.null(p0)) NA else quantile(p0[keep0], 0.05, names = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        N = N, scenario = sc, rate = rate, subset = ss, method = meth,
        RawPw = mean(p[keep] < 0.05),
        T1    = if (is.null(p0)) NA else mean(p0[keep0] < 0.05),
        AdjPw = if (is.na(thr)) NA else mean(p[keep] < thr),
        na_rate = mean(is.na(src[keep, paste0("coef_", meth)])))
    }
  }
}
long <- do.call(rbind, rows)
write.csv(long, file.path(v3b_dir, "summary_v3b_long.csv"), row.names = FALSE)

for (val in c("RawPw", "T1", "AdjPw")) {
  if (all(is.na(long[[val]]))) { cat(sprintf("\n[%s] not available yet (no T1 files)\n", val)); next }
  d <- long[, c("subset", "scenario", "N", "rate", "method", val)]
  w <- reshape(d, idvar = c("subset", "scenario", "N", "rate"), timevar = "method", direction = "wide")
  names(w) <- sub(paste0("^", val, "\\."), "", names(w))
  w <- w[, c("subset", "scenario", "N", "rate", intersect(meths, names(w)))]
  w$ratio_vs_IPW_DR   <- w$wrep_cal0.3_ind / w$IPW_DR
  w$ratio_vs_OW_DR    <- w$wrep_cal0.3_ind / w$OW_DR
  w$ratio_vs_PSspline <- w$wrep_cal0.3_ind / w$LR_PSspline
  w$ratio_vs_LR_PS    <- w$wrep_cal0.3_ind / w$LR_PS
  w <- w[order(w$subset, w$scenario, -w$N, w$rate), ]
  write.csv(w, file.path(v3b_dir, sprintf("summary_v3b_key_%s.csv", val)), row.names = FALSE)
  cat(sprintf("\n==================== %s ====================\n", val))
  num <- sapply(w, is.numeric); w[num & !names(w) %in% c("N", "rate")] <- round(w[num & !names(w) %in% c("N", "rate")], 3)
  print(w, row.names = FALSE)
}
cat("\nsame-data check, max |dp| over all files:", signif(max(checks), 3),
    if (max(checks) < 1e-8) "(OK)" else "(MISMATCH - do not use)", "\n")
cat("max NA rate across methods:", round(max(long$na_rate), 4), "\n")
