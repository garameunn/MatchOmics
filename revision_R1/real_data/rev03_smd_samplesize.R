## rev03_smd_samplesize.R
##  [R3-6] post-matching SMD with matching weights; pass counts under max SMD < 0.2
##  [R1-4] matched sample size by caliper (n, Neff, % treated excluded)
##  type "wrep_pkg" (cal0.3 only): package/simulation matching for comparison
##          (round 2 without replacement, w = m/mean(m)); GEE variants -> rev07
## Usage: Rscript rev03_smd_samplesize.R     (env REV_RAW, REV_BASE, REV_LIB)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
SMD_THRESH <- 0.2
D <- load_realdata()
t0 <- Sys.time()

rows <- list()
for (pl in names(D)) {
  d <- D[[pl]]
  for (i in seq_along(d$mk)) {
    mk <- d$mk[i]
    if (i %% 50 == 1) cat(sprintf("[%s] %d/%d %s  (%.1f min)\n", pl, i, length(d$mk), mk,
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    df <- prep_marker(d$df, mk)
    if (is.null(df)) next
    sdp <- pre_sd(df)
    df$w_inv <- 1; df$w_freq <- 1
    pre <- smd_set(df, sdp)

    for (cal in names(calipers)) {
      for (rt in c("wrep", "worep", "wrep_pkg")) {
        if (rt == "wrep_pkg" && cal != "cal0.3") next
        sub <- switch(rt,
          wrep     = run_wrep(df, calipers[[cal]]),
          worep    = run_worep(df, calipers[[cal]]),
          wrep_pkg = run_wrep(df, calipers[[cal]], weight_type = "freq", round2_replace = FALSE))
        if (is.null(sub)) next
        inf <- attr(sub, "info")
        s   <- smd_set(sub, sdp)
        n_excl <- inf$n_trt_total - inf$n_trt_r1 - inf$n_trt_r2
        rows[[length(rows) + 1]] <- data.frame(
          platform = pl, marker = mk, caliper = cal, type = rt,
          n_total = nrow(df), n_trt_total = inf$n_trt_total, n_ctrl_total = inf$n_ctrl_total,
          n_matched = nrow(sub), n_trt_matched_r1 = inf$n_trt_r1, n_trt_rescued_r2 = inf$n_trt_r2,
          n_trt_excluded = n_excl, pct_trt_excluded = 100 * n_excl / inf$n_trt_total,
          n_ctrl_unique = inf$n_ctrl_unique, n_pairs = inf$n_pairs, n_clusters = inf$n_clusters,
          max_ctrl_reuse = max(sub$m[sub$taxa_class == 0]),
          neff_freq = kish(sub$m), neff_inv = kish(sub$w_inv),
          smd_pre_max = pre["max_uw"],
          smd_max_unweighted = s["max_uw"], smd_max_w_inv = s["max_w_inv"],
          smd_max_w_freq = s["max_w_freq"], smd_max_w_freq_preSD = s["max_wfix_freq"],
          stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(rev_dir, "rev03_matching_per_marker.csv"), row.names = FALSE)

## ---------- R3-6: pass counts ------------------------------------------------
smd_defs <- c(published_w_inv = "smd_max_w_inv", unweighted = "smd_max_unweighted",
              weighted_freq = "smd_max_w_freq", weighted_freq_preSD = "smd_max_w_freq_preSD")
pass <- list()
for (pl in names(D)) {
  d   <- D[[pl]]
  old <- read.csv(d$final, stringsAsFactors = FALSE)$marker
  for (tp in c("wrep", "wrep_pkg")) {
    r <- res[res$platform == pl & res$caliper == "cal0.3" & res$type == tp, ]
    for (nm in names(smd_defs)) {
      if (tp == "wrep_pkg" && nm == "published_w_inv") next
      ok <- r$marker[!is.na(r[[smd_defs[nm]]]) & r[[smd_defs[nm]]] < SMD_THRESH]
      pass[[length(pass) + 1]] <- data.frame(platform = pl, matching = tp, smd_definition = nm,
        n_markers = length(d$mk), n_pass = length(ok),
        n_pass_published = length(old), n_overlap_published = length(intersect(ok, old)),
        n_concordant_pass = sum(d$conc %in% ok), n_concordant = length(d$conc),
        concordant_failing = paste(setdiff(d$conc, ok), collapse = ";"),
        stringsAsFactors = FALSE)
    }
  }
}
pass <- do.call(rbind, pass)
write.csv(pass, file.path(rev_dir, "rev03_smd_pass_counts.csv"), row.names = FALSE)
print(pass[, 1:9])

## 26 concordant markers: SMD detail at cal0.3
conc_tab <- res[res$caliper == "cal0.3" & res$type %in% c("wrep", "worep", "wrep_pkg") &
                ((res$platform == "Protein" & res$marker %in% conc_prot) |
                 (res$platform == "Metabolite" & res$marker %in% conc_meta)),
                c("platform", "marker", "type", "n_matched", "neff_freq", "smd_pre_max",
                  "smd_max_unweighted", "smd_max_w_inv", "smd_max_w_freq", "smd_max_w_freq_preSD")]
write.csv(conc_tab, file.path(rev_dir, "rev03_smd_concordant26.csv"), row.names = FALSE)

## ---------- R1-4: sample size by caliper ---------------------------------
q <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, .25), quantile(x, .75))
q1 <- function(x) sprintf("%.1f [%.1f-%.1f]", median(x), quantile(x, .25), quantile(x, .75))
summ <- list()
for (pl in names(D)) {
  fin <- read.csv(D[[pl]]$final, stringsAsFactors = FALSE)$marker
  for (set in c("SMD_passed", "all_markers")) {
    for (tp in c("wrep", "worep")) {
      for (cal in names(calipers)) {
        r <- res[res$platform == pl & res$type == tp & res$caliper == cal, ]
        if (set == "SMD_passed") r <- r[r$marker %in% fin, ]
        summ[[length(summ) + 1]] <- data.frame(platform = pl, marker_set = set, n_markers_set = nrow(r),
          method = ifelse(tp == "wrep", "GEE-EXW (two-round)", "GEE-EXWO (one-round)"), caliper = cal,
          N_total = q(r$n_total), N_treated = q(r$n_trt_total),
          n_matched = q(r$n_matched), neff = q(r$neff_freq),
          pct_treated_excluded = q1(r$pct_trt_excluded),
          pct_treated_rescued_r2 = q1(100 * r$n_trt_rescued_r2 / r$n_trt_total),
          n_matched_mean = round(mean(r$n_matched), 1), neff_mean = round(mean(r$neff_freq), 1),
          pct_excl_mean = round(mean(r$pct_trt_excluded), 2),
          stringsAsFactors = FALSE)
      }
    }
  }
}
summ <- do.call(rbind, summ)
write.csv(summ, file.path(rev_dir, "rev06_caliper_samplesize_summary.csv"), row.names = FALSE)
print(summ[summ$marker_set == "SMD_passed", c("platform", "method", "caliper", "n_matched", "neff", "pct_treated_excluded")])

## 26 concordant markers by caliper
cz <- res[res$type %in% c("wrep", "worep") &
          ((res$platform == "Protein" & res$marker %in% conc_prot) |
           (res$platform == "Metabolite" & res$marker %in% conc_meta)),
          c("platform", "marker", "type", "caliper", "n_total", "n_trt_total", "n_matched",
            "neff_freq", "n_trt_excluded", "pct_trt_excluded", "n_trt_rescued_r2")]
write.csv(cz, file.path(rev_dir, "rev06_caliper_samplesize_concordant26.csv"), row.names = FALSE)

cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
