## rev07_gee_sortfix.R -- GEE-EXW re-estimated with the geeM weight-alignment fix
##  (see rev_common.R: geem() misaligns weights unless rows are sorted by id)
##  Variants of GEE-EXW (two-round matching, median exposure):
##    A  published : realdata_analysis.R as run (unsorted, w = 1/k for reused controls)
##    B  sortfix   : same matching and weights, rows sorted by cluster  (bug fix only)
##    B2 sort_freq : same matching, sorted, w = m/mean(m)
##    C  package   : round 2 without replacement, w = m/mean(m), sorted (cal0.3 only;
##                   = MatchOmics package / simulation definition)
##  GEE-EXWO (weights all 1) is checked to be unaffected by sorting.
##  Concordance re-derived with LR_PS, LR_COM, GEE-EXWO from protein/metabolite_final.csv
## Usage: Rscript rev07_gee_sortfix.R   (env REV_RAW, REV_BASE, REV_LIB)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
D <- load_realdata()
t0 <- Sys.time()
wald <- function(r) if (is.na(r["se"])) NA_real_ else 2 * pnorm(-abs(r["beta"] / r["se"]))

rows <- list()
add <- function(pl, mk, cal, variant, sub, r)
  rows[[length(rows) + 1]] <<- data.frame(platform = pl, marker = mk, caliper = cal, variant = variant,
    n = if (is.null(sub)) NA else nrow(sub), beta = r["beta"], se = r["se"], pval = wald(r),
    stringsAsFactors = FALSE)

for (pl in names(D)) {
  d <- D[[pl]]
  for (i in seq_along(d$mk)) {
    mk <- d$mk[i]
    if (i %% 50 == 1) cat(sprintf("[%s] %d/%d (%.1f min)\n", pl, i, length(d$mk),
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    df <- prep_marker(d$df, mk)
    if (is.null(df)) next
    for (cal in names(calipers)) {
      sw <- run_wrep(df, calipers[[cal]])
      add(pl, mk, cal, "A_published", sw, fit_gee(sw, sort_id = FALSE))
      add(pl, mk, cal, "B_sortfix",   sw, fit_gee(sw, sort_id = TRUE))
      if (!is.null(sw)) sw$wt <- sw$w_freq
      add(pl, mk, cal, "B2_sort_freq", sw, fit_gee(sw, sort_id = TRUE))
      if (cal == "cal0.3") {
        sc <- run_wrep(df, calipers[[cal]], weight_type = "freq", round2_replace = FALSE)
        add(pl, mk, cal, "C_package", sc, fit_gee(sc, sort_id = TRUE))
        so <- run_worep(df, calipers[[cal]])
        add(pl, mk, cal, "EXWO_unsorted", so, fit_gee(so, sort_id = FALSE))
        add(pl, mk, cal, "EXWO_sorted",   so, fit_gee(so, sort_id = TRUE))
      }
    }
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(rev_dir, "rev07_gee_variants_per_marker.csv"), row.names = FALSE)

## ---- checks -----------------------------------------------------------------
chk <- function(v1, v2, cal = "cal0.3") {
  a <- res[res$variant == v1 & res$caliper == cal, ]; b <- res[res$variant == v2 & res$caliper == cal, ]
  m <- merge(a, b, by = c("platform", "marker")); max(abs(m$beta.x - m$beta.y), na.rm = TRUE)
}
cat(sprintf("EXWO sorted vs unsorted: max |beta diff| = %.2e\n", chk("EXWO_unsorted", "EXWO_sorted")))
for (pl in names(D)) {
  old <- read.csv(D[[pl]]$results, stringsAsFactors = FALSE)
  a <- res[res$platform == pl & res$variant == "A_published", ]
  dif <- sapply(names(calipers), function(cal) {
    x <- a[a$caliper == cal, ]; y <- old[match(x$marker, old$marker), paste0("wrep_", cal, "_beta")]
    max(abs(x$beta - y), na.rm = TRUE) })
  cat(pl, "A_published vs published wrep beta, max |diff| by caliper:", signif(dif, 2), "\n")
}

## ---- concordance on the SMD-passed sets ----------------------------------------
variants <- c("A_published", "B_sortfix", "B2_sort_freq", "C_package")
pp <- function(b, s) 2 * pnorm(-abs(b / s))
cnt <- list(); memb <- list()
for (pl in names(D)) {
  fin <- read.csv(D[[pl]]$final, stringsAsFactors = FALSE, check.names = FALSE)
  base <- data.frame(marker = fin$marker, p_PS = fin$LR_PS_pval, p_COM = fin$LR_COM_pval,
                     p_EXWO = pp(fin$`worep_cal0.3_beta`, fin$`worep_cal0.3_se`), stringsAsFactors = FALSE)
  for (v in variants) {
    e <- res[res$platform == pl & res$variant == v & res$caliper == "cal0.3", ]
    b <- base; b$p_EXW <- e$pval[match(b$marker, e$marker)]; b$beta_EXW <- e$beta[match(b$marker, e$marker)]
    q <- lapply(b[, c("p_PS", "p_COM", "p_EXWO", "p_EXW")], p.adjust, method = "BH")
    raw4 <- with(b, !is.na(p_PS + p_COM + p_EXWO + p_EXW) & p_PS < .05 & p_COM < .05 & p_EXWO < .05 & p_EXW < .05)
    fdr4 <- Reduce(`&`, lapply(q, function(z) !is.na(z) & z < .05))
    cnt[[length(cnt) + 1]] <- data.frame(platform = pl, variant = v, n_markers = nrow(b),
      EXW_sig_raw = sum(b$p_EXW < .05, na.rm = TRUE), EXW_sig_fdr = sum(q$p_EXW < .05, na.rm = TRUE),
      EXW_not_converged = sum(is.na(b$p_EXW)),
      concordant_raw = sum(raw4), concordant_fdr = sum(fdr4),
      lost_vs_26 = paste(setdiff(D[[pl]]$conc, b$marker[raw4]), collapse = ";"),
      gained_vs_26 = paste(setdiff(b$marker[raw4], D[[pl]]$conc), collapse = ";"),
      fdr_set = paste(sort(b$marker[fdr4]), collapse = ";"), stringsAsFactors = FALSE)
    memb[[length(memb) + 1]] <- data.frame(platform = pl, variant = v, marker = b$marker,
      beta_EXW = b$beta_EXW, p_EXW = b$p_EXW, fdr_EXW = q$p_EXW, concordant_raw = raw4, concordant_fdr = fdr4,
      in_published26 = b$marker %in% D[[pl]]$conc, stringsAsFactors = FALSE)
  }
}
cnt <- do.call(rbind, cnt); memb <- do.call(rbind, memb)
write.csv(cnt, file.path(rev_dir, "rev07_concordance_by_variant.csv"), row.names = FALSE)
write.csv(memb[memb$concordant_raw | memb$in_published26 | memb$concordant_fdr, ],
          file.path(rev_dir, "rev07_concordant_markers_by_variant.csv"), row.names = FALSE)
print(cnt[, 1:8]); cat("\n"); print(cnt[, c("platform", "variant", "lost_vs_26", "gained_vs_26")])

## ---- the 26 published concordant markers under each variant -------------------
c26 <- memb[memb$in_published26, c("platform", "marker", "variant", "beta_EXW", "p_EXW")]
c26w <- reshape(c26, idvar = c("platform", "marker"), timevar = "variant", direction = "wide")
write.csv(c26w, file.path(rev_dir, "rev07_concordant26_EXW_variants.csv"), row.names = FALSE)

## ---- STable 6 wrep columns with the sort fix (B) --------------------------------
fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.1e", p)))
st6 <- do.call(rbind, lapply(names(D), function(pl) {
  o <- data.frame(type = pl, marker = sort(D[[pl]]$conc), stringsAsFactors = FALSE)
  for (cal in c("nocal", "cal0.3", "cal0.05")) for (v in c("A_published", "B_sortfix", "B2_sort_freq")) {
    e <- res[res$platform == pl & res$variant == v & res$caliper == cal, ]
    e <- e[match(o$marker, e$marker), ]
    o[[paste0("wrep_", cal, "_", v)]] <- ifelse(is.na(e$beta), NA, paste0(sprintf("%.2f", e$beta), " (", fmt_p(e$pval), ")"))
  }
  o }))
write.csv(st6, file.path(rev_dir, "rev07_stable6_wrep_sortfix.csv"), row.names = FALSE)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
