## ============================================================
##  10_summarize_v3.R — summarise v3 power (and T1, if available)
##  Usage: Rscript 10_summarize_v3.R            (from code/)
##  Output: ../results/v3/summary_v3_long.csv   (all methods)
##          ../results/v3/summary_v3_key.csv    (key methods, wide)
##  Raw power  = mean(p < 0.05)
##  T1         = mean(p_null < 0.05)                       (needs T1 file)
##  Adj power  = mean(p < q_0.05(p_null)), method-specific (needs T1 file)
##  Subsets: all taxa / top 20 taxa by psdifflist (as in analyze_top20.R)
## ============================================================
base_dir <- "/path/to/Matching2026"
res_dir  <- file.path(base_dir, "results", "v3")
load(file.path(base_dir, "data", "psdifflist.RData"))
norder <- c("1000" = 1, "500" = 2, "200" = 3)

pick_file <- function(dir, N, rate, sc) {
  f <- list.files(dir, sprintf("^NSAMP%d_REP[0-9]+_RATE%.1f_%s\\.RData$", N, rate, sc), full.names = TRUE)
  if (length(f) == 0) return(NULL)
  reps <- as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(f)))
  f[which.max(reps)]                         # largest replicate count
}
load_mat <- function(f) { e <- new.env(); load(f, envir = e); e$pvalues[[length(e$pvalues)]] }

rows <- list()
for (N in c(1000, 500, 200)) for (sc in c("S1", "S2", "S3")) for (rate in c(0.1, 0.2, 0.3)) {
  fp <- pick_file(file.path(res_dir, "Power"), N, rate, sc)
  if (is.null(fp)) next
  pw <- load_mat(fp)
  ft <- pick_file(file.path(res_dir, "T1"), N, rate, sc)
  t1 <- if (is.null(ft)) NULL else load_mat(ft)

  psd  <- psdifflist[[norder[as.character(N)]]]
  dimm <- length(psd)
  taxon_of <- function(m) (seq_len(nrow(m)) - 1) %% dimm + 1
  top20 <- order(psd, decreasing = TRUE)[1:20]

  for (s in grep("^sig_", colnames(pw), value = TRUE)) {
    meth <- sub("^sig_", "", s)
    p  <- pw[, s]; p[is.na(p)] <- 1
    p0 <- if (!is.null(t1) && s %in% colnames(t1)) { x <- t1[, s]; x[is.na(x)] <- 1; x } else NULL
    for (ss in c("all", "top20")) {
      keep  <- if (ss == "all") rep(TRUE, length(p)) else taxon_of(pw) %in% top20
      keep0 <- if (is.null(p0)) NULL else if (ss == "all") rep(TRUE, length(p0)) else taxon_of(t1) %in% top20
      thr   <- if (is.null(p0)) NA else quantile(p0[keep0], 0.05, names = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        N = N, scenario = sc, rate = rate, subset = ss, method = meth,
        RawPw = mean(p[keep] < 0.05),
        T1    = if (is.null(p0)) NA else mean(p0[keep0] < 0.05),
        AdjPw = if (is.na(thr)) NA else mean(p[keep] < thr),
        na_rate = mean(is.na(pw[keep, sub("^sig_", "coef_", s)])),
        rep_power = as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(fp))),
        rep_t1    = if (is.null(ft)) NA else as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(ft))))
    }
  }
}
long <- do.call(rbind, rows)
write.csv(long, file.path(res_dir, "summary_v3_long.csv"), row.names = FALSE)

key <- c("LR_unadj", "LR_PS", "LR_Combat", "LR_PSspline", "IPW", "OW",
         "wrep_nocal_ind", "wrep_cal0.3_ind", "wrep_cal0.05_ind",
         "worep_nocal_ind", "worep_cal0.3_ind")
wide_of <- function(val) {
  d <- long[long$method %in% key, c("subset", "scenario", "N", "rate", "method", val)]
  w <- reshape(d, idvar = c("subset", "scenario", "N", "rate"), timevar = "method", direction = "wide")
  names(w) <- sub(paste0("^", val, "\\."), "", names(w))
  w <- w[, c("subset", "scenario", "N", "rate", intersect(key, names(w)))]
  w$ratio_vs_LR_PS     <- w$wrep_cal0.3_ind / w$LR_PS
  w$ratio_vs_PSspline  <- w$wrep_cal0.3_ind / w$LR_PSspline
  w$ratio_vs_Combat    <- w$wrep_cal0.3_ind / w$LR_Combat
  w$ratio_wrep_worep   <- w$wrep_cal0.3_ind / w$worep_cal0.3_ind
  w <- w[order(w$subset, w$scenario, -w$N, w$rate), ]
  w
}
for (val in c("RawPw", "T1", "AdjPw")) {
  if (all(is.na(long[[val]]))) { cat(sprintf("\n[%s] not available yet (no T1 files)\n", val)); next }
  w <- wide_of(val)
  write.csv(w, file.path(res_dir, sprintf("summary_v3_key_%s.csv", val)), row.names = FALSE)
  cat(sprintf("\n==================== %s ====================\n", val))
  num <- sapply(w, is.numeric); w[num & !names(w) %in% c("N", "rate")] <- round(w[num & !names(w) %in% c("N", "rate")], 3)
  print(w, row.names = FALSE)
}
cat("\nmax NA rate across methods:", round(max(long$na_rate), 4), "\n")
