## ============================================================
##  18_renumber_stables_citation_order.R — supplementary tables numbered by
##  order of first citation in the revised manuscript (decided 2026-10-28)
##  Input : manuscript/Tables_rev/STables_rev.xlsx   (from 16_make_figures_tables_v3.R)
##  Output: manuscript/Tables_rev/STables_final.xlsx (STable 1-14 in new order)
##  The earlier output numbered 7-14 after the old 1-6 is kept as
##  STables_final_old7to14.xlsx (renamed by hand before running this script).
##  Run from the project root:
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\18_renumber_stables_citation_order.R
## ============================================================
suppressPackageStartupMessages(library(openxlsx))
src <- "manuscript/Tables_rev/STables_rev.xlsx"
out <- "manuscript/Tables_rev/STables_final.xlsx"

## source sheet in STables_rev.xlsx -> new number (old number in the 1st submission, if any)
map <- data.frame(
  sheet = c("New - matching order", "STable 1", "STable 2", "STable 3", "New - Top20 T1", "New - Top20 AdjPw",
            "STable 4", "STable 5", "New - Working corr", "New - raw vs FDR", "New - balance excl",
            "New - cutoff", "STable 6", "New - real data n"),
  new = 1:14,
  old = c(NA, 1, 2, 3, NA, NA, 4, 5, NA, NA, NA, NA, 6, NA),
  stringsAsFactors = FALSE)

wb <- loadWorkbook(src)
missing <- setdiff(map$sheet, names(wb))
if (length(missing)) stop("sheets not found: ", paste(missing, collapse = ", "))
for (i in seq_len(nrow(map))) {           # temporary names first to avoid clashes (e.g. "STable 2")
  renameWorksheet(wb, map$sheet[i], paste0("tmp_", i))
}
for (i in seq_len(nrow(map))) {
  s <- paste0("tmp_", i)
  title <- read.xlsx(wb, sheet = s, rows = 1, cols = 1, colNames = FALSE)[1, 1]
  title <- sub("^Supplementary Table (\\[No\\.\\]|[0-9]+)\\.", paste0("Supplementary Table ", map$new[i], "."), title)
  writeData(wb, s, title, startRow = 1, startCol = 1)
  renameWorksheet(wb, s, paste("STable", map$new[i]))
}
worksheetOrder(wb) <- order(as.integer(sub("^STable ", "", names(wb))))
saveWorkbook(wb, out, overwrite = TRUE)
write.csv(transform(map, old = ifelse(is.na(old), "new", paste0("STable ", old)))[, c("new", "old", "sheet")],
          "manuscript/Tables_rev/STable_numbering_map.csv", row.names = FALSE)
message("saved ", out, " and STable_numbering_map.csv")
