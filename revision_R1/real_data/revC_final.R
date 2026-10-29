## revC_final.R -- specification C (final, = MatchOmics package / manuscript Methods)
##  GEE-EXW : round 1 nearest without replacement; round 2 = unmatched treated +
##            ALL controls, nearest WITHOUT replacement; w = m / mean(m)
##            (m = number of matched pairs a subject belongs to); rows sorted by
##            cluster before geem() (geeM 0.10.1 weight-alignment fix)
##  GEE-EXWO: one round nearest without replacement; all weights 1; sorted
##  SMD     : computed on the C matched sample with w = m / mean(m) (smd_w)
##  LR_unadj / LR_PS / LR_COM do not depend on matching -> taken from the
##  validated protein_results.csv / metabolite_results.csv (rev00).
## Output: realdata/revision/C_final/C_all_markers.csv, C_matching_long.csv
## Usage : Rscript revC_final.R   (env REV_RAW, REV_BASE, REV_LIB)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
## c_dir comes from rev_common.R (env REV_CDIR, default C_final)
D  <- load_realdata()
## LR_* source: recomputed in this environment by rev08_lr_recompute.R if present
## (c_dir/LR_recomputed.csv), otherwise the validated published *_results.csv
lr_file <- file.path(c_dir, "LR_recomputed.csv")
LR <- if (file.exists(lr_file)) read.csv(lr_file, stringsAsFactors = FALSE) else NULL
cat("LR source:", if (is.null(LR)) "published *_results.csv" else lr_file, "\n")
t0 <- Sys.time()
wald <- function(r) if (is.na(r["se"])) NA_real_ else unname(2 * pnorm(-abs(r["beta"] / r["se"])))

wide <- list(); long <- list()
for (pl in names(D)) {
  d   <- D[[pl]]
  old <- read.csv(d$results, stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_along(d$mk)) {
    mk <- d$mk[i]
    if (i %% 50 == 1) cat(sprintf("[%s] %d/%d (%.1f min)\n", pl, i, length(d$mk),
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    df <- prep_marker(d$df, mk)
    if (is.null(df)) next
    sdp <- pre_sd(df); df$w_inv <- 1; df$w_freq <- 1
    pre <- smd_set(df, sdp)
    o <- if (is.null(LR)) old[old$marker == mk, ] else LR[LR$platform == pl & LR$marker == mk, ]
    row <- data.frame(platform = pl, marker = mk, stringsAsFactors = FALSE)
    for (m in c("LR_unadj", "LR_PS", "LR_COM"))
      for (s in c("pval", "beta", "se", "n")) row[[paste0(m, "_", s)]] <- o[[paste0(m, "_", s)]]
    row$smd_pre_max <- pre["max_uw"]
    for (cal in names(calipers)) {
      for (rt in c("wrep", "worep")) {
        sub <- if (rt == "wrep") run_wrep(df, calipers[[cal]], weight_type = "freq", round2_replace = FALSE)
               else run_worep(df, calipers[[cal]])
        r   <- fit_gee(sub, sort_id = TRUE)
        key <- paste0(rt, "_", cal)
        row[[paste0(key, "_pval")]] <- wald(r)
        row[[paste0(key, "_beta")]] <- unname(r["beta"])
        row[[paste0(key, "_se")]]   <- unname(r["se"])
        row[[paste0(key, "_n")]]    <- if (is.null(sub)) NA else nrow(sub)
        if (is.null(sub)) next
        s   <- smd_set(sub, sdp)
        inf <- attr(sub, "info")
        ## C SMD definition: weighted (w = m/mean(m)); for worep weights are 1 (= unweighted, pop. variance)
        row[[paste0(key, "_smd_max")]] <- unname(s["max_w_freq"])
        n_excl <- inf$n_trt_total - inf$n_trt_r1 - inf$n_trt_r2
        long[[length(long) + 1]] <- data.frame(platform = pl, marker = mk, caliper = cal,
          method = ifelse(rt == "wrep", "GEE-EXW", "GEE-EXWO"),
          n_total = nrow(df), n_trt_total = inf$n_trt_total, n_ctrl_total = inf$n_ctrl_total,
          n_matched = nrow(sub), n_trt_matched_r1 = inf$n_trt_r1, n_trt_rescued_r2 = inf$n_trt_r2,
          n_trt_excluded = n_excl, pct_trt_excluded = 100 * n_excl / inf$n_trt_total,
          n_ctrl_unique = inf$n_ctrl_unique, n_pairs = inf$n_pairs, n_clusters = inf$n_clusters,
          max_ctrl_reuse = max(sub$m[sub$taxa_class == 0]), neff = kish(sub$m),
          smd_max_weighted = unname(s["max_w_freq"]), smd_max_unweighted = unname(s["max_uw"]),
          smd_max_weighted_preSD = unname(s["max_wfix_freq"]), stringsAsFactors = FALSE)
      }
    }
    wide[[length(wide) + 1]] <- row
  }
}
## rbind with possibly differing columns (failed matchings) -> fill NA
cols <- unique(unlist(lapply(wide, names)))
wide <- do.call(rbind, lapply(wide, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
long <- do.call(rbind, long); rownames(wide) <- rownames(long) <- NULL
write.csv(wide, file.path(c_dir, "C_all_markers.csv"), row.names = FALSE)
write.csv(long, file.path(c_dir, "C_matching_long.csv"), row.names = FALSE)

## check: GEE-EXWO under C equals the published worep fit (same matching, weights 1)
for (pl in names(D)) {
  old <- read.csv(D[[pl]]$results, stringsAsFactors = FALSE)
  w <- wide[wide$platform == pl, ]
  dd <- sapply(names(calipers), function(cal) max(abs(w[[paste0("worep_", cal, "_beta")]] -
          old[match(w$marker, old$marker), paste0("worep_", cal, "_beta")]), na.rm = TRUE))
  cat(pl, "GEE-EXWO C vs published, max |beta diff| by caliper:", signif(dd, 2), "\n")
}
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
