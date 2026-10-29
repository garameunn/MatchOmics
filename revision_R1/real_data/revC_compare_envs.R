## revC_compare_envs.R -- compare every CSV present in two C output folders
## (default: C_final [R 3.6.3] vs C_final_R43 [R 4.3.1]).
## Numeric columns: max |diff| and max relative diff; other columns: n mismatches.
## Output: <dir2>/compare_vs_<dir1>.csv
## Usage : Rscript revC_compare_envs.R [dir1] [dir2]   (env REV_BASE)
base_dir <- Sys.getenv("REV_BASE"); if (!nzchar(base_dir)) stop("Set REV_BASE")
a <- commandArgs(TRUE)
d1 <- file.path(base_dir, "realdata", "revision", if (length(a) >= 1) a[1] else "C_final")
d2 <- file.path(base_dir, "realdata", "revision", if (length(a) >= 2) a[2] else "C_final_R43")
f  <- intersect(list.files(d1, "\\.csv$"), list.files(d2, "\\.csv$"))
rows <- list()
for (fn in f) {
  x <- read.csv(file.path(d1, fn), stringsAsFactors = FALSE, check.names = FALSE)
  y <- read.csv(file.path(d2, fn), stringsAsFactors = FALSE, check.names = FALSE)
  same_shape <- identical(dim(x), dim(y)) && identical(names(x), names(y))
  if (!same_shape) {
    rows[[length(rows) + 1]] <- data.frame(file = fn, same_shape = FALSE, n_num_cols = NA, max_abs_diff = NA,
      max_rel_diff = NA, worst_column = NA, n_char_mismatch = NA,
      note = sprintf("dims %s vs %s", paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x")),
      stringsAsFactors = FALSE)
    next
  }
  num <- names(x)[sapply(x, is.numeric) & sapply(y, is.numeric)]
  oth <- setdiff(names(x), num)
  ad <- vapply(num, function(cn) { d <- abs(x[[cn]] - y[[cn]]); if (all(is.na(d))) 0 else max(d, na.rm = TRUE) }, 0)
  rd <- vapply(num, function(cn) { d <- abs(x[[cn]] - y[[cn]]) / pmax(abs(x[[cn]]), 1e-300)
                                   if (all(is.na(d))) 0 else max(d, na.rm = TRUE) }, 0)
  na_mis <- sum(vapply(num, function(cn) sum(is.na(x[[cn]]) != is.na(y[[cn]])), 0))
  ch_mis <- sum(vapply(oth, function(cn) sum(!(x[[cn]] == y[[cn]] | (is.na(x[[cn]]) & is.na(y[[cn]]))), na.rm = TRUE), 0))
  rows[[length(rows) + 1]] <- data.frame(file = fn, same_shape = TRUE, n_num_cols = length(num),
    max_abs_diff = if (length(ad)) max(ad) else 0, max_rel_diff = if (length(rd)) max(rd) else 0,
    worst_column = if (length(rd) && max(rd) > 0) names(rd)[which.max(rd)] else "",
    n_char_mismatch = ch_mis, note = if (na_mis) sprintf("%d NA pattern mismatches", na_mis) else "",
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(d2, paste0("compare_vs_", basename(d1), ".csv")), row.names = FALSE)
print(res, right = FALSE)
cat("\nfiles only in", basename(d1), ":", paste(setdiff(list.files(d1, "\\.csv$"), f), collapse = ", "), "\n")
cat("files only in", basename(d2), ":", paste(setdiff(list.files(d2, "\\.csv$"), f), collapse = ", "), "\n")
