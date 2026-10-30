## ============================================================
##  21_stable_readable_headers.R — readable column names and exact p-values
##  for Supplementary Tables 1, 10, 12, 13 and 14 (decided 2026-10-30)
##  The earlier sheets showed code variable names (e.g. n_raw_p05, wrep_cal0.3) and
##  p-values rounded to 3 decimals (p < 0.001 printed as 0). Content is unchanged;
##  sheets are rebuilt from the real-data result files.
##  Inputs : 결과파일_20260517/realdata/revision/C_final/
##             STable_order_sensitivity_C.csv, C_sig_counts_raw_vs_fdr.csv,
##             rev04_cutoff_concordant26_C.csv, STable6_C.csv, C_caliper_samplesize_summary.csv
##           manuscript/Tables_rev/STables_final.xlsx
##  Output : manuscript/Tables_rev/STables_final.xlsx (sheets above replaced;
##           previous workbook kept as STables_final_before_headers.xlsx)
##  Run from the project root:
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\21_stable_readable_headers.R
## ============================================================
suppressPackageStartupMessages(library(openxlsx))
rd   <- "결과파일_20260517/realdata/revision/C_final"
book <- "manuscript/Tables_rev/STables_final.xlsx"
file.copy(book, "manuscript/Tables_rev/STables_final_before_headers.xlsx", overwrite = TRUE)
rcsv <- function(f) read.csv(file.path(rd, f), check.names = FALSE, stringsAsFactors = FALSE)

fp  <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
fb  <- function(b) sprintf("%.2f", b)
yn  <- function(x) ifelse(as.logical(x), "Yes", "No")
meth <- c(LR_unadj = "LR", LR_PS = "LR-PS", LR_COM = "LR-COM", `GEE-EXW` = "GEE-EXW", `GEE-EXWO` = "GEE-EXWO")
cal  <- c(nocal = "None", cal0.3 = "0.3", cal0.2 = "0.2", cal0.1 = "0.1", cal0.05 = "0.05")

wb <- loadWorkbook(book)
bold <- createStyle(textDecoration = "bold")
put <- function(name, legend, df) {
  title <- read.xlsx(book, sheet = name, rows = 1, cols = 1, colNames = FALSE)[1, 1]
  removeWorksheet(wb, name); addWorksheet(wb, name)
  writeData(wb, name, title, startRow = 1); addStyle(wb, name, bold, rows = 1, cols = 1)
  writeData(wb, name, legend, startRow = 2)
  top <- sub("\\|.*", "", names(df)); bot <- sub(".*\\|", "", names(df))
  two <- any(grepl("\\|", names(df)))
  if (two) {
    top[duplicated(top)] <- ""; bot[!grepl("\\|", names(df))] <- ""
    writeData(wb, name, t(top), startRow = 4, colNames = FALSE)
    writeData(wb, name, t(bot), startRow = 5, colNames = FALSE)
    addStyle(wb, name, bold, rows = 4:5, cols = seq_along(top), gridExpand = TRUE)
    writeData(wb, name, df, startRow = 6, colNames = FALSE)
  } else {
    writeData(wb, name, df, startRow = 4, headerStyle = bold)
  }
  setColWidths(wb, name, cols = seq_along(df), widths = "auto")
}

## Supplementary Table 1: matching-order sensitivity
o <- rcsv("STable_order_sensitivity_C.csv")
stopifnot(all(o$n_not_converged == 0), all(o$n_random_orders == 20))
st1 <- data.frame(
  Platform = o$platform, Marker = o$marker, Borderline = yn(o$borderline),
  `Default order|β` = fb(o$beta_default), `Default order|p` = fp(o$p_default),
  `Default order|Matched n` = o$n_matched_default,
  `20 random orders|Significant runs (p < 0.05)` = paste0(o$n_sig_p05, "/20"),
  `20 random orders|β range` = paste0(fb(o$beta_min), " to ", fb(o$beta_max)),
  `20 random orders|p range` = paste0(fp(o$p_min), " to ", fp(o$p_max)),
  `20 random orders|Same sign in all runs` = yn(o$beta_sign_consistent),
  `20 random orders|Matched n range` = sub("-", "\u2013", o$n_matched_range), check.names = FALSE)
put("STable 1", paste("GEE-EXW with a 0.3 SD caliper for the 22 concordant markers. Default order: high-abundance subjects processed",
                      "in decreasing order of the PS (primary result). Random orders: 20 runs with a randomised processing order in both",
                      "matching rounds (all converged). Borderline: significant at p < 0.05 in fewer than all 20 runs."), st1)

## Supplementary Table 10: significant markers at raw p and FDR
s <- rcsv("C_sig_counts_raw_vs_fdr.csv")
st10 <- data.frame(
  Platform = s$platform,
  Method = ifelse(grepl("^Concordant", s$method), "All four adjusted methods (concordant)", meth[s$method]),
  `Markers tested` = s$n_markers, `Significant at raw p < 0.05` = s$n_raw_p05,
  `Significant at FDR < 0.05` = s$n_fdr05, check.names = FALSE)
put("STable 10", paste("Markers passing the balance criterion (191 proteins, 34 metabolites). FDR by the Benjamini-Hochberg procedure",
                       "within each platform. Matching methods use a 0.3 SD caliper. Concordant: significant under all four adjusted",
                       "methods (LR-PS, LR-COM, GEE-EXW and GEE-EXWO)."), st10)

## Supplementary Table 12: cutoff sensitivity
k <- rcsv("rev04_cutoff_concordant26_C.csv")
st12 <- data.frame(
  Platform = k$platform, Marker = k$marker, Method = k$method,
  `Median split (primary)|Matched n` = k$n_matched.median, `Median split (primary)|β` = fb(k$beta.median),
  `Median split (primary)|p` = fp(k$pval.median),
  `Tertiles|Matched n` = k$n_matched.tertile, `Tertiles|β` = fb(k$beta.tertile), `Tertiles|p` = fp(k$pval.tertile),
  `Tertiles|Same sign as median` = yn(k$sign_agree_tertile),
  `Quartiles|Matched n` = k$n_matched.quartile, `Quartiles|β` = fb(k$beta.quartile), `Quartiles|p` = fp(k$pval.quartile),
  `Quartiles|Same sign as median` = yn(k$sign_agree_quartile), check.names = FALSE)
put("STable 12", paste("The 22 concordant markers with the matched groups defined by the median split (primary), by tertiles (upper vs",
                       "lower third) or by quartiles (upper vs lower quarter); the PS was re-estimated in each analysis set and the",
                       "marker entered the outcome model as a continuous variable. 0.3 SD caliper."), st12)

## Supplementary Table 13: caliper sensitivity
c13 <- rcsv("STable6_C.csv")
st13 <- data.frame(Platform = c13$type, Marker = c13$marker, check.names = FALSE)
for (rp in c("wrep", "worep")) for (cl in c("nocal", "cal0.3", "cal0.05")) {
  lab <- if (rp == "wrep") "GEE-EXW" else "GEE-EXWO"
  st13[[paste0(lab, ", β (p)|Caliper ", cal[[cl]])]] <- c13[[paste0(rp, "_", cl)]]
}
for (rp in c("wrep", "worep")) for (cl in c("nocal", "cal0.3", "cal0.05")) {
  lab <- if (rp == "wrep") "GEE-EXW" else "GEE-EXWO"
  st13[[paste0(lab, ", matched n|Caliper ", cal[[cl]])]] <- c13[[paste0(rp, "_", cl, "_n")]]
}
put("STable 13", paste("Coefficient (exact p) and matched sample size of the 22 concordant markers (Table 5) under with-replacement",
                       "(two-round, GEE-EXW) and without-replacement (GEE-EXWO) matching, without a caliper and with 0.3 and 0.05 SD calipers."), st13)

## Supplementary Table 14: matched sample size by caliper
z <- rcsv("C_caliper_samplesize_summary.csv"); z <- z[z$marker_set == "C_marker_set", ]
st14 <- data.frame(
  Platform = z$platform, Markers = z$n_markers_set,
  Matching = ifelse(z$method == "GEE-EXW", "With replacement (two-round)", "Without replacement (single round)"),
  Caliper = cal[z$caliper], `Total N` = z$N_total, `High-abundance N` = z$N_treated,
  `Matched n` = z$n_matched, `Effective n` = z$neff, `High-abundance excluded (%)` = z$pct_treated_excluded,
  `High-abundance matched in round 2 (%)` = z$pct_treated_rescued_r2, `Maximum control reuse` = z$max_control_reuse,
  check.names = FALSE)
put("STable 14", paste("Median [interquartile range] over the markers passing the balance criterion. Effective n by the Kish",
                       "approximation. High-abundance excluded: high-abundance subjects without a match."), st14)

worksheetOrder(wb) <- order(as.integer(sub("^STable ", "", names(wb))))
saveWorkbook(wb, book, overwrite = TRUE)
message("Supplementary Tables 1, 10, 12, 13, 14 rewritten with readable headers")
