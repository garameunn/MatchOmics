## ============================================================
##  15_extract_fig_data_v3.R — small CSVs for figures/tables (run on server)
##  Reads the large .RData results and writes compact CSVs that
##  16_make_figures_tables_v3.R turns into figures and tables locally.
##
##  Output (results/v3/fig_data/):
##    t1_alpha_v3.csv      T1 at alpha = 0.001, 0.005, 0.01, 0.05, 0.1
##                         (all N x scenario x rate x subset, v3 + v3b methods)
##                         -> Figure 1, Supplementary Figure 2
##    hseries_adjpw_v3.csv adjusted power at h2 = 0.001, 0.002, 0.005, 0.01
##                         (S1, N = 1000; from 13_Sim_v3_hseries.R)
##                         -> Supplementary Figure 3
##    matched_size_v3.csv  mean matched n, effective n, cluster size, reuse
##                         by N x caliper x matching type -> Supplementary Table 1
##  Same definitions as 10_summarize_v3.R (largest replicate file; NA p = 1;
##  adjusted power threshold = 5th percentile of the method's null p-values).
##  Usage (from code/): /usr/local/R-4.4/bin/Rscript 15_extract_fig_data_v3.R
## ============================================================
base_dir <- "/path/to/Matching2026"
res      <- file.path(base_dir, "results")
out_dir  <- file.path(res, "v3", "fig_data"); dir.create(out_dir, showWarnings = FALSE)
load(file.path(base_dir, "data", "psdifflist.RData"))
norder <- c("1000" = 1, "500" = 2, "200" = 3)
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

pick_file <- function(dir, N, rate, sc) {
  f <- list.files(dir, sprintf("^NSAMP%d_REP[0-9]+_RATE%.1f_%s\\.RData$", N, rate, sc), full.names = TRUE)
  if (length(f) == 0) return(NULL)
  f[which.max(as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(f))))]
}
load_all <- function(f) { e <- new.env(); load(f, envir = e); e$pvalues }
load_last <- function(f) { p <- load_all(f); p[[length(p)]] }
pclean <- function(x) { x[is.na(x)] <- 1; x }
subset_rows <- function(m, N, ss) {
  if (ss == "all") return(rep(TRUE, nrow(m)))
  psd <- psdifflist[[norder[as.character(N)]]]
  ((seq_len(nrow(m)) - 1) %% length(psd) + 1) %in% order(psd, decreasing = TRUE)[1:20]
}

## ── 1) T1 at several alpha ────────────────────────────────────
rows <- list()
for (N in c(1000, 500, 200)) for (sc in c("S1", "S2", "S3")) for (rate in c(0.1, 0.2, 0.3)) {
  srcs <- list(v3 = pick_file(file.path(res, "v3", "T1"), N, rate, sc),
               v3b = pick_file(file.path(res, "v3b", "T1"), N, rate, sc))
  for (nm in names(srcs)) {
    if (is.null(srcs[[nm]])) next
    m <- load_last(srcs[[nm]])
    sig <- grep("^sig_", colnames(m), value = TRUE)
    if (nm == "v3b") sig <- intersect(sig, c("sig_IPW_DR", "sig_OW_DR"))
    for (ss in c("all", "top20")) {
      k <- subset_rows(m, N, ss)
      for (s in sig) {
        p <- pclean(m[k, s])
        rows[[length(rows) + 1]] <- data.frame(N = N, scenario = sc, rate = rate, subset = ss,
          method = sub("^sig_", "", s), alpha = alphas, T1 = sapply(alphas, function(a) mean(p < a)))
      }
    }
  }
  message(sprintf("  T1 alpha: N=%d %s rate=%.1f", N, sc, rate))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "t1_alpha_v3.csv"), row.names = FALSE)

## ── 2) h2 series (SFig 3) ─────────────────────────────────────
rows <- list(); checks <- c()
for (rate in c(0.1, 0.2, 0.3)) {
  fh <- pick_file(file.path(res, "v3_hseries", "Power"), 1000, rate, "S1")
  ft <- pick_file(file.path(res, "v3", "T1"), 1000, rate, "S1")
  fp <- pick_file(file.path(res, "v3", "Power"), 1000, rate, "S1")
  if (is.null(fh) || is.null(ft)) { message("  [hseries] missing for rate ", rate); next }
  ph <- load_all(fh); t1 <- load_last(ft)
  # same-data check: b = 0.01 columns must equal the main v3 power file
  if (!is.null(fp) && "b0.01" %in% names(ph)) {
    pm <- load_last(fp); cc <- intersect(c("sig_LR_PS", "sig_wrep_cal0.3_ind"), colnames(ph[["b0.01"]]))
    checks[sprintf("rate%.1f", rate)] <- max(abs(ph[["b0.01"]][, cc] - pm[, cc]), na.rm = TRUE)
  }
  for (bn in names(ph)) {
    m <- ph[[bn]]; b <- as.numeric(sub("^b", "", bn))
    for (s in intersect(grep("^sig_", colnames(m), value = TRUE), colnames(t1))) {
      for (ss in c("all", "top20")) {
        thr <- quantile(pclean(t1[subset_rows(t1, 1000, ss), s]), 0.05, names = FALSE)
        p <- pclean(m[subset_rows(m, 1000, ss), s])
        rows[[length(rows) + 1]] <- data.frame(N = 1000, scenario = "S1", rate = rate, subset = ss,
          h2 = b, method = sub("^sig_", "", s), RawPw = mean(p < 0.05), AdjPw = mean(p < thr))
      }
    }
  }
}
if (length(rows)) write.csv(do.call(rbind, rows), file.path(out_dir, "hseries_adjpw_v3.csv"), row.names = FALSE)
if (length(checks)) message(sprintf("  [hseries same-data check] max|dp| = %.2e %s", max(checks),
                                    if (max(checks) < 1e-8) "(OK)" else "(MISMATCH)"))

## ── 3) matched sample size (outcome-independent; S1 rate 0.1 file) ──
rows <- list()
for (N in c(1000, 500, 200)) {
  fp <- pick_file(file.path(res, "v3", "Power"), N, 0.1, "S1"); if (is.null(fp)) next
  m <- load_last(fp)
  for (type in c("wrep", "worep")) for (cn in c("nocal", "cal0.3", "cal0.2", "cal0.1", "cal0.05")) {
    g <- function(pre) { col <- paste0(pre, "_", type, "_", cn); if (col %in% colnames(m)) mean(m[, col], na.rm = TRUE) else NA }
    rows[[length(rows) + 1]] <- data.frame(N = N, type = type, caliper = cn,
      matched_n = g("nrow"), neff = g("neff"), cluster_size = g("csize"),
      mean_reuse = if (type == "wrep") g("meanrep") else NA)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "matched_size_v3.csv"), row.names = FALSE)

message("15_extract_fig_data_v3.R complete. Files in: ", out_dir)
print(list.files(out_dir))
