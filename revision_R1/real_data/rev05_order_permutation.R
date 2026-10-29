## rev05_order_permutation.R -- [R2-6] order dependence of greedy two-round matching
##  A) row order shuffled 20x (MatchIt default m.order = "largest" PS first;
##     row order then only matters for ties)
##  B) extra: m.order = "random" with 20 seeds (treated processed in random
##     order in both rounds) -- the stronger test of greedy-order dependence
##  cal0.3, GEE-EXW, 26 concordant markers
## Usage: Rscript rev05_order_permutation.R [n_perm=20]  (env REV_RAW, REV_BASE, REV_LIB)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
args  <- commandArgs(TRUE)
NPERM <- if (length(args)) as.integer(args[1]) else 20L
D  <- load_realdata()
t0 <- Sys.time()
cal <- calipers[["cal0.3"]]

rows <- list()
for (pl in names(D)) {
  d <- D[[pl]]
  for (mk in d$conc) {
    cat(sprintf("[%s] %s (%.1f min)\n", pl, mk, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    base_df <- prep_marker(d$df, mk)
    ref <- run_wrep(base_df, cal)
    ref_ids <- sort(rownames(ref))
    rr  <- fit_gee(ref)
    rows[[length(rows) + 1]] <- data.frame(platform = pl, marker = mk, scheme = "original", rep = 0L,
      n_matched = nrow(ref), same_matched_set = TRUE, beta = rr["beta"],
      pval = 2 * pnorm(-abs(rr["beta"] / rr["se"])), stringsAsFactors = FALSE)
    for (scheme in c("row_shuffle", "m.order_random")) {
      for (b in seq_len(NPERM)) {
        set.seed(1000 * b + nchar(mk))
        if (scheme == "row_shuffle") {
          dfp <- base_df[sample(nrow(base_df)), ]
          sub <- run_wrep(dfp, cal)
        } else {
          sub <- run_wrep(base_df, cal, m.order = "random")
        }
        r <- fit_gee(sub)
        rows[[length(rows) + 1]] <- data.frame(platform = pl, marker = mk, scheme = scheme, rep = b,
          n_matched = if (is.null(sub)) NA else nrow(sub),
          same_matched_set = !is.null(sub) && identical(sort(rownames(sub)), ref_ids),
          beta = r["beta"], pval = if (is.na(r["se"])) NA else 2 * pnorm(-abs(r["beta"] / r["se"])),
          stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(rev_dir, out_name("rev05_order_permutation_runs.csv")), row.names = FALSE)

orig <- res[res$scheme == "original", c("platform", "marker", "beta", "pval")]
names(orig)[3:4] <- c("beta_original", "p_original")
rp <- res[res$scheme != "original", ]
summ <- do.call(rbind, lapply(split(rp, list(rp$platform, rp$marker, rp$scheme), drop = TRUE), function(x)
  data.frame(platform = x$platform[1], marker = x$marker[1], scheme = x$scheme[1], n_rep = nrow(x),
             n_converged = sum(!is.na(x$pval)),
             n_matched_min = min(x$n_matched, na.rm = TRUE), n_matched_max = max(x$n_matched, na.rm = TRUE),
             pct_same_matched_set = 100 * mean(x$same_matched_set),
             beta_min = min(x$beta, na.rm = TRUE), beta_max = max(x$beta, na.rm = TRUE),
             p_min = min(x$pval, na.rm = TRUE), p_median = median(x$pval, na.rm = TRUE),
             p_max = max(x$pval, na.rm = TRUE),
             n_sig_p05 = sum(x$pval < 0.05, na.rm = TRUE),
             stringsAsFactors = FALSE)))
summ <- merge(orig, summ, by = c("platform", "marker"))
summ$sig_decision_changed <- summ$n_sig_p05 < summ$n_converged
summ$beta_sign_stable <- sign(summ$beta_min) == sign(summ$beta_max)
summ <- summ[order(summ$scheme, summ$platform != "Protein", summ$marker), ]
write.csv(summ, file.path(rev_dir, out_name("rev05_order_permutation_summary.csv")), row.names = FALSE)
print(summ[, c("scheme", "platform", "marker", "p_original", "p_min", "p_max", "n_sig_p05", "n_converged", "pct_same_matched_set")])
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
