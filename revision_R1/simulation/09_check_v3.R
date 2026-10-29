## ============================================================
##  09_check_v3.R — quick summary of one v3 output file
##  Usage: Rscript 09_check_v3.R <path/to/NSAMP..._S1.RData> [alpha=0.05]
##  Prints: failures, prevalence, realised h2, raw rejection rate per
##          method, mean estimated exchangeable rho, matched sizes.
## ============================================================
args  <- commandArgs(trailingOnly = TRUE)
alpha <- if (length(args) > 1) as.numeric(args[2]) else 0.05
load(args[1])   # pvalues, meta
str(meta[c("mode", "nsamp", "scenario", "rate", "conf_share", "rep_n", "dimm")])

for (bn in names(pvalues)) {
  m <- pvalues[[bn]]
  cat(sprintf("\n==== %s | rows=%d | failed rows=%d ====\n", bn, nrow(m), sum(is.na(m[, "prev"]))))
  cat(sprintf("prevalence %.3f | realised h2 %.5f | b2 %.4f | b3 %.4f\n",
              mean(m[, "prev"], na.rm = TRUE), mean(m[, "h2_real"], na.rm = TRUE),
              mean(m[, "b2_true"], na.rm = TRUE), mean(m[, "b3_true"], na.rm = TRUE)))

  sig <- grep("^sig_", colnames(m), value = TRUE)
  rr  <- sapply(sig, function(s) { p <- m[, s]; p[is.na(p)] <- 1; mean(p < alpha) })
  na  <- sapply(sig, function(s) mean(is.na(m[, sub("^sig_", "coef_", s)])))
  cat("\nrejection rate at alpha =", alpha, "(NA fits counted as non-rejection)\n")
  print(round(data.frame(reject = rr, na_rate = na, row.names = sub("^sig_", "", sig)), 4))

  rho <- grep("^rho_", colnames(m), value = TRUE)
  if (length(rho)) {
    cat("\nestimated exchangeable rho (mean / median)\n")
    print(round(t(sapply(rho, function(r) c(mean = mean(m[, r], na.rm = TRUE),
                                            median = median(m[, r], na.rm = TRUE)))), 4))
  }
  sz <- grep("^(nrow|neff|csize)_", colnames(m), value = TRUE)
  cat("\nmatched sample diagnostics (mean)\n")
  print(round(colMeans(m[, sz, drop = FALSE], na.rm = TRUE), 2))
}
