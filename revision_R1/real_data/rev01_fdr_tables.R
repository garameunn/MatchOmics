## rev01_fdr_tables.R -- [R2-5] raw p vs BH-FDR counts, FDR concordance
##                       [R3-7] Table 5 / STable 6 regenerated with exact p
## Input : realdata/protein_final.csv, metabolite_final.csv (SMD-passed markers)
## Output: realdata/revision/rev01_*.csv
## Base R only.  Usage: Rscript rev01_fdr_tables.R   (env REV_BASE)
base_dir <- Sys.getenv("REV_BASE"); if (!nzchar(base_dir)) stop("Set REV_BASE")
res_dir  <- file.path(base_dir, "realdata")
rev_dir  <- file.path(res_dir, "revision"); dir.create(rev_dir, showWarnings = FALSE)

methods <- c(LR_unadj = "LR_unadj", LR_PS = "LR_PS", LR_COM = "LR_COM",
             `GEE-EXW` = "wrep_cal0.3", `GEE-EXWO` = "worep_cal0.3")
conc4   <- c("LR_PS", "LR_COM", "wrep_cal0.3", "worep_cal0.3")
ALPHA   <- 0.05

D <- list(Protein    = read.csv(file.path(res_dir, "protein_final.csv"),    stringsAsFactors = FALSE, check.names = FALSE),
          Metabolite = read.csv(file.path(res_dir, "metabolite_final.csv"), stringsAsFactors = FALSE, check.names = FALSE))

## geeM's summary() rounds p to 7 decimals (p = 0 means p < 5e-8). Recompute the
## Wald p from the stored beta / sandwich SE so exact p can be reported.
for (pl in names(D)) {
  d <- D[[pl]]
  for (pc in grep("^w(o)?rep_.*_pval$", names(d), value = TRUE)) {
    pre  <- sub("_pval$", "", pc)
    prec <- 2 * pnorm(-abs(d[[paste0(pre, "_beta")]] / d[[paste0(pre, "_se")]]))
    dd   <- max(abs(prec - d[[pc]]), na.rm = TRUE)
    if (dd > 1e-6) warning(sprintf("%s %s: recomputed Wald p differs (max %.2e)", pl, pc, dd))
    d[[pc]] <- prec
  }
  D[[pl]] <- d
}

## BH within platform over the SMD-passed marker set (all p columns)
for (pl in names(D)) {
  d <- D[[pl]]
  for (pc in grep("_pval$", names(d), value = TRUE)) {
    fc <- sub("_pval$", "_fdr_rev", pc)
    d[[fc]] <- p.adjust(d[[pc]], method = "BH")
    old <- sub("_pval$", "_fdr", pc)
    if (old %in% names(d)) {
      dd <- max(abs(d[[fc]] - d[[old]]), na.rm = TRUE)
      if (dd > 1e-10) warning(sprintf("%s %s: stored _fdr differs (max %.2e)", pl, pc, dd))
    }
  }
  D[[pl]] <- d
}

## ---- 1. counts by method --------------------------------------------------
cnt <- do.call(rbind, lapply(names(D), function(pl) {
  d <- D[[pl]]
  do.call(rbind, lapply(names(methods), function(nm) {
    p <- d[[paste0(methods[nm], "_pval")]]; q <- d[[paste0(methods[nm], "_fdr_rev")]]
    data.frame(platform = pl, method = nm, n_markers = nrow(d), n_tested = sum(!is.na(p)),
               n_raw_p05 = sum(p < ALPHA, na.rm = TRUE), n_fdr05 = sum(q < ALPHA, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}))
## concordance rows
for (pl in names(D)) {
  d <- D[[pl]]
  raw4 <- apply(d[, paste0(conc4, "_pval")], 1, function(r) all(!is.na(r) & r < ALPHA))
  fdr4 <- apply(d[, paste0(conc4, "_fdr_rev")], 1, function(r) all(!is.na(r) & r < ALPHA))
  cnt <- rbind(cnt, data.frame(platform = pl, method = "Concordant (LR_PS, LR_COM, EXW, EXWO)",
                               n_markers = nrow(d), n_tested = NA,
                               n_raw_p05 = sum(raw4), n_fdr05 = sum(fdr4)))
  D[[pl]]$conc_raw <- raw4; D[[pl]]$conc_fdr <- fdr4
}
write.csv(cnt, file.path(rev_dir, "rev01_sig_counts_raw_vs_fdr.csv"), row.names = FALSE)
print(cnt)

## ---- formatting -------------------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                     ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.1e", p)))
cell  <- function(b, p) ifelse(is.na(b) | is.na(p), NA_character_, paste0(sprintf("%.2f", b), " (", fmt_p(p), ")"))

## ---- 1b. FDR-concordant marker list (supplementary) ------------------------
fdr_list <- do.call(rbind, lapply(names(D), function(pl) {
  d <- D[[pl]][D[[pl]]$conc_fdr, ]
  if (!nrow(d)) return(NULL)
  o <- data.frame(platform = pl, marker = d$marker, in_raw_concordant_set = d$conc_raw,
                  stringsAsFactors = FALSE)
  for (nm in names(methods)) {
    pre <- methods[nm]
    o[[paste0(nm, "_beta")]] <- round(d[[paste0(pre, "_beta")]], 4)
    o[[paste0(nm, "_p")]]    <- signif(d[[paste0(pre, "_pval")]], 3)
    o[[paste0(nm, "_fdr")]]  <- signif(d[[paste0(pre, "_fdr_rev")]], 3)
  }
  o[order(o$marker), ]
}))
write.csv(fdr_list, file.path(rev_dir, "rev01_fdr_concordant_markers.csv"), row.names = FALSE)

## full per-marker raw/FDR table for the concordance-relevant methods
full <- do.call(rbind, lapply(names(D), function(pl) {
  d <- D[[pl]]
  o <- data.frame(platform = pl, marker = d$marker, conc_raw = d$conc_raw, conc_fdr = d$conc_fdr,
                  stringsAsFactors = FALSE)
  for (nm in names(methods)) {
    o[[paste0(nm, "_p")]]   <- d[[paste0(methods[nm], "_pval")]]
    o[[paste0(nm, "_fdr")]] <- d[[paste0(methods[nm], "_fdr_rev")]]
  }
  o
}))
write.csv(full, file.path(rev_dir, "rev01_all_markers_p_fdr.csv"), row.names = FALSE)

## ---- 2. Table 5 (raw p concordance, exact p) -------------------------------
robust <- do.call(rbind, lapply(names(D), function(pl) {
  d <- D[[pl]][D[[pl]]$conc_raw, ]; d$type <- pl; d }))
robust <- robust[order(robust$type == "Metabolite", robust$marker), ]

main <- data.frame(type = robust$type, marker = robust$marker, stringsAsFactors = FALSE)
for (nm in names(methods))
  main[[nm]] <- cell(robust[[paste0(methods[nm], "_beta")]], robust[[paste0(methods[nm], "_pval")]])
write.csv(main, file.path(rev_dir, "robust_main_table_rev.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## numeric companion (beta, exact p, FDR) for the manuscript text
num <- data.frame(type = robust$type, marker = robust$marker, stringsAsFactors = FALSE)
for (nm in names(methods)) {
  num[[paste0(nm, "_beta")]] <- robust[[paste0(methods[nm], "_beta")]]
  num[[paste0(nm, "_p")]]    <- robust[[paste0(methods[nm], "_pval")]]
  num[[paste0(nm, "_fdr")]]  <- robust[[paste0(methods[nm], "_fdr_rev")]]
}
write.csv(num, file.path(rev_dir, "robust_main_table_rev_numeric.csv"), row.names = FALSE)

## ---- 2b. STable 6 (caliper sensitivity, exact p) ---------------------------
supp_m <- c("wrep_nocal", "wrep_cal0.3", "wrep_cal0.05", "worep_nocal", "worep_cal0.3", "worep_cal0.05")
supp <- data.frame(type = robust$type, marker = robust$marker, stringsAsFactors = FALSE)
for (m in supp_m) supp[[m]] <- cell(robust[[paste0(m, "_beta")]], robust[[paste0(m, "_pval")]])
for (m in supp_m) supp[[paste0(m, "_n")]] <- robust[[paste0(m, "_n")]]
write.csv(supp, file.path(rev_dir, "robust_suppl_table_rev.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---- compare with the existing tables --------------------------------------
old_main <- read.csv(file.path(res_dir, "robust_main_table.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cat("\nexisting robust_main_table.csv markers:", nrow(old_main),
    "| new (raw p):", nrow(main),
    "| identical marker set:", setequal(paste(old_main$type, old_main$marker), paste(main$type, main$marker)), "\n")
sig_by_method <- sapply(supp_m, function(m) sum(robust[[paste0(m, "_pval")]] < ALPHA, na.rm = TRUE))
cat("STable 6: markers with p<0.05 out of", nrow(robust), "\n"); print(sig_by_method)
cat("not significant in STable 6 cells:\n")
for (m in supp_m) { k <- robust$marker[!(robust[[paste0(m, "_pval")]] < ALPHA) | is.na(robust[[paste0(m, "_pval")]])]
  if (length(k)) cat(" ", m, ":", paste(sprintf("%s(p=%s)", k, fmt_p(robust[[paste0(m, "_pval")]][robust$marker %in% k])), collapse = ", "), "\n") }
