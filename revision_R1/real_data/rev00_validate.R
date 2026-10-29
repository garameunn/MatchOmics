## rev00_validate.R -- check that this R / MatchIt / geeM build reproduces the
## published per-marker results (protein_results.csv, metabolite_results.csv)
## Usage: Rscript rev00_validate.R [n_markers_per_platform]  (0 = all)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
args <- commandArgs(TRUE)
nmk  <- if (length(args)) as.integer(args[1]) else 10L

cat("R", R.version$major, R.version$minor, "| MatchIt", as.character(packageVersion("MatchIt")),
    "| geeM", as.character(packageVersion("geeM")), "\n")
D <- load_realdata()
out <- list()
t0 <- Sys.time()
for (pl in names(D)) {
  d   <- D[[pl]]
  old <- read.csv(d$results, stringsAsFactors = FALSE, check.names = FALSE)
  mks <- if (nmk == 0) d$mk else unique(c(head(d$conc, ceiling(nmk/2)), head(setdiff(d$mk, d$conc), floor(nmk/2))))
  for (mk in mks) {
    df <- prep_marker(d$df, mk)
    for (cal in c("nocal", "cal0.3", "cal0.05")) {
      for (rt in c("wrep", "worep")) {
        sub <- if (rt == "wrep") run_wrep(df, calipers[[cal]]) else run_worep(df, calipers[[cal]])
        r <- fit_gee(sub)
        key <- paste0(rt, "_", cal)
        o <- old[old$marker == mk, ]
        out[[length(out) + 1]] <- data.frame(platform = pl, marker = mk, method = key,
          n_new = if (is.null(sub)) NA else nrow(sub), n_old = o[[paste0(key, "_n")]],
          beta_new = r["beta"], beta_old = o[[paste0(key, "_beta")]],
          p_new = r["pval"], p_old = o[[paste0(key, "_pval")]], stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
res$n_match <- res$n_new == res$n_old
res$beta_absdiff <- abs(res$beta_new - res$beta_old)
res$p_absdiff <- abs(res$p_new - res$p_old)
write.csv(res, if (c_dir_name == "C_final") file.path(rev_dir, "rev00_validation.csv")
               else file.path(c_dir, "rev00_validation.csv"), row.names = FALSE)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("n identical: %d / %d ; max |beta diff| = %.2e ; max |p diff| = %.2e\n",
            sum(res$n_match, na.rm = TRUE), nrow(res),
            max(res$beta_absdiff, na.rm = TRUE), max(res$p_absdiff, na.rm = TRUE)))
print(aggregate(cbind(beta_absdiff, p_absdiff) ~ method, data = res, FUN = max))
