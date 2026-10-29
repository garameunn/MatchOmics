## ============================================================
##  19_stable11_excluded_markers.R — Supplementary Table 11 = all markers
##  excluded by the weighted post-matching balance criterion (decided 2026-10-29)
##  Input : 결과파일_20260517/realdata/revision/C_final/C_all_markers.csv
##          manuscript/Tables_rev/STables_final.xlsx (from 18_renumber_stables_citation_order.R)
##  Output: manuscript/Tables_rev/STables_final.xlsx, sheet "STable 11" replaced
##          (previous workbook kept as STables_final_before_STable11.xlsx)
##  Run from the project root:
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\19_stable11_excluded_markers.R
## ============================================================
suppressPackageStartupMessages(library(openxlsx))
src  <- "결과파일_20260517/realdata/revision/C_final/C_all_markers.csv"
book <- "manuscript/Tables_rev/STables_final.xlsx"
file.copy(book, "manuscript/Tables_rev/STables_final_before_STable11.xlsx", overwrite = TRUE)

d <- read.csv(src, check.names = FALSE)
ex <- d[d[["wrep_cal0.3_smd_max"]] >= 0.2, ]
ex <- ex[order(factor(ex$platform, levels = c("Protein", "Metabolite")), -ex[["wrep_cal0.3_smd_max"]]), ]
tab <- data.frame(
  Platform = ex$platform,
  Marker = ex$marker,
  `Max SMD before matching` = round(ex$smd_pre_max, 3),
  `Max weighted SMD after matching` = round(ex[["wrep_cal0.3_smd_max"]], 3),
  `Matched n` = ex[["wrep_cal0.3_n"]],
  check.names = FALSE)
n_tab <- table(ex$platform)
stopifnot(n_tab[["Protein"]] == 266 - 191, n_tab[["Metabolite"]] == 127 - 34)

wb <- loadWorkbook(book)
pos <- which(names(wb) == "STable 11")
removeWorksheet(wb, "STable 11")
addWorksheet(wb, "STable 11")
bold <- createStyle(textDecoration = "bold")
writeData(wb, "STable 11", "Supplementary Table 11. Markers excluded by the post-matching balance criterion in the FROM cohort.", startRow = 1)
addStyle(wb, "STable 11", bold, rows = 1, cols = 1)
writeData(wb, "STable 11", paste0(
  "Markers whose maximum weighted standardized mean difference (SMD) over age, sex, BMI, smoking status and duration of residence was >= 0.2 ",
  "after two-round matching with a 0.3 SD caliper (", n_tab[["Protein"]], " of 266 proteins and ", n_tab[["Metabolite"]], " of 127 metabolites). ",
  "Weights are the matching frequency weights used in the GEE (w = m / mean(m)); SMD before matching is unweighted in the full sample. ",
  "These markers were not included in the concordance analysis."), startRow = 2)
writeData(wb, "STable 11", tab, startRow = 4, headerStyle = bold)
setColWidths(wb, "STable 11", cols = 1:5, widths = c(12, 22, 22, 28, 12))
worksheetOrder(wb) <- order(as.integer(sub("^STable ", "", names(wb))))
saveWorkbook(wb, book, overwrite = TRUE)
message(sprintf("STable 11 replaced: %d proteins, %d metabolites excluded", n_tab[["Protein"]], n_tab[["Metabolite"]]))
