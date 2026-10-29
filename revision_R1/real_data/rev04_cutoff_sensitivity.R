## rev04_cutoff_sensitivity.R -- [R2-7] exposure cutoff sensitivity
##  median (published) vs tertile (top 1/3 vs bottom 1/3) [vs quartile, extra]
##  PS re-estimated in the retained subjects; two-round matching cal0.3 -> GEE-EXW
##  (GEE-EXWO also reported). Marker enters GEE as continuous in all cases.
## Usage: Rscript rev04_cutoff_sensitivity.R   (env REV_RAW, REV_BASE, REV_LIB)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
D <- load_realdata()
t0 <- Sys.time()
cal <- calipers[["cal0.3"]]

rows <- list()
for (pl in names(D)) {
  d   <- D[[pl]]
  fin <- read.csv(d$final, stringsAsFactors = FALSE)$marker     # SMD-passed set
  for (i in seq_along(fin)) {
    mk <- fin[i]
    if (i %% 50 == 1) cat(sprintf("[%s] %d/%d (%.1f min)\n", pl, i, length(fin),
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    for (ex in c("median", "tertile", "quartile")) {
      df <- prep_marker(d$df, mk, exposure = ex)
      if (is.null(df)) next
      sdp <- pre_sd(df)
      for (rt in c("wrep", "worep")) {
        sub <- if (rt == "wrep") run_wrep(df, cal) else run_worep(df, cal)
        r   <- fit_gee(sub)
        s   <- if (is.null(sub)) rep(NA, 4) else smd_set(sub, sdp)
        rows[[length(rows) + 1]] <- data.frame(platform = pl, marker = mk, exposure = ex,
          method = ifelse(rt == "wrep", "GEE-EXW", "GEE-EXWO"),
          n_analysed = nrow(df), n_matched = if (is.null(sub)) NA else nrow(sub),
          beta = r["beta"], se = r["se"],
          pval = if (is.na(r["se"])) NA else 2 * pnorm(-abs(r["beta"] / r["se"])),
          smd_max_w_inv = s[2], smd_max_w_freq = s[3],
          concordant26 = mk %in% d$conc, stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(rev_dir, out_name("rev04_cutoff_per_marker.csv")), row.names = FALSE)

## wide table for the 26 concordant markers
key <- c("platform", "marker", "method")
w <- reshape(res[res$concordant26, c(key, "exposure", "n_matched", "beta", "pval")],
             idvar = key, timevar = "exposure", direction = "wide")
w$sign_agree_tertile  <- sign(w$beta.tertile)  == sign(w$beta.median)
w$sign_agree_quartile <- sign(w$beta.quartile) == sign(w$beta.median)
w$sig_tertile  <- w$pval.tertile  < 0.05
w$sig_quartile <- w$pval.quartile < 0.05
w$beta_ratio_tertile <- w$beta.tertile / w$beta.median
w <- w[order(w$method, w$platform != "Protein", w$marker), ]
write.csv(w, file.path(rev_dir, out_name("rev04_cutoff_concordant26.csv")), row.names = FALSE)

## summary
summ <- list()
for (m in c("GEE-EXW", "GEE-EXWO")) for (pl in names(D)) {
  x <- w[w$method == m & w$platform == pl, ]
  a <- res[res$method == m & res$platform == pl, ]
  aw <- reshape(a[, c("marker", "exposure", "beta", "pval")], idvar = "marker", timevar = "exposure", direction = "wide")
  for (ex in c("tertile", "quartile")) {
    bm <- aw$beta.median; be <- aw[[paste0("beta.", ex)]]
    pm <- aw$pval.median; pe <- aw[[paste0("pval.", ex)]]
    summ[[length(summ) + 1]] <- data.frame(method = m, platform = pl, cutoff = ex,
      conc26_n = nrow(x),
      conc26_sign_agree = sum(sign(x[[paste0("beta.", ex)]]) == sign(x$beta.median), na.rm = TRUE),
      conc26_p05 = sum(x[[paste0("pval.", ex)]] < 0.05, na.rm = TRUE),
      conc26_median_p = signif(median(x[[paste0("pval.", ex)]], na.rm = TRUE), 3),
      conc26_max_p = signif(max(x[[paste0("pval.", ex)]], na.rm = TRUE), 3),
      all_n = nrow(aw),
      all_sig_median = sum(pm < 0.05, na.rm = TRUE), all_sig_cutoff = sum(pe < 0.05, na.rm = TRUE),
      all_sig_both = sum(pm < 0.05 & pe < 0.05, na.rm = TRUE),
      all_sign_agree_among_sig_median = sum(pm < 0.05 & sign(be) == sign(bm), na.rm = TRUE),
      all_beta_cor = round(cor(bm, be, use = "complete.obs"), 3),
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, summ)
write.csv(summ, file.path(rev_dir, out_name("rev04_cutoff_summary.csv")), row.names = FALSE)
print(summ)
cat("\nconcordant markers with p >= 0.05 under tertile:\n")
print(w[!w$sig_tertile | is.na(w$sig_tertile), c("method", "platform", "marker", "pval.median", "pval.tertile", "pval.quartile")])
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
