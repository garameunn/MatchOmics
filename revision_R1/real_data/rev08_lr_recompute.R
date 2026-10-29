## rev08_lr_recompute.R -- recompute LR_unadj, LR_PS, LR_COM for all markers in
## the current R environment (code copied from realdata_analysis.R) and compare
## with protein_results.csv / metabolite_results.csv.
##  protein   LR_COM: sva::ComBat on log2 protein matrix (batch = Batch, par.prior)
##  metabolite LR_COM: pre-corrected file 20250428_Batch_Corrected_Metabolites_All.txt
## Output: <c_dir>/LR_recomputed.csv (used by revC_final.R), LR_recompute_check.csv
## Usage : Rscript rev08_lr_recompute.R   (env REV_RAW, REV_BASE, REV_CDIR)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
suppressPackageStartupMessages(library(sva))
cat("R", R.version$major, R.version$minor, "| sva", as.character(packageVersion("sva")), "\n")
D <- load_realdata()

fit_glm_safe <- function(formula_str, dat) {
  fit <- tryCatch(glm(as.formula(formula_str), data = dat, family = "binomial"), error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_, NA_real_))
  cs <- coef(summary(fit))
  if (!"mk_val" %in% rownames(cs)) return(c(NA_real_, NA_real_, NA_real_))
  c(cs["mk_val", "Pr(>|z|)"], cs["mk_val", "Estimate"], cs["mk_val", "Std. Error"])
}

## protein ComBat (identical call to realdata_analysis.R)
pd <- D$Protein$df
prot_mat <- as.matrix(pd[, D$Protein$mk])
cmbt_prot <- t(ComBat(t(prot_mat), batch = pd$Batch, par.prior = TRUE, prior.plots = FALSE))
rownames(cmbt_prot) <- pd$ID; colnames(cmbt_prot) <- D$Protein$mk

## metabolite pre-corrected file (identical to realdata_analysis.R)
mc_raw <- read.delim(file.path(raw_dir, "20250428_Batch_Corrected_Metabolites_All.txt"),
                     stringsAsFactors = FALSE, check.names = FALSE)
mc_raw$ID <- gsub("-", "_", mc_raw$ID)
md <- D$Metabolite$df
mc_mk <- intersect(D$Metabolite$mk, names(mc_raw))
cmbt_meta <- as.data.frame(mc_raw[, mc_mk, drop = FALSE]); rownames(cmbt_meta) <- mc_raw$ID
cmbt_meta <- cmbt_meta[intersect(rownames(cmbt_meta), rownames(md)), , drop = FALSE]
CM <- list(Protein = cmbt_prot, Metabolite = cmbt_meta)

rows <- list()
for (pl in names(D)) {
  d <- D[[pl]]
  for (mk in d$mk) {
    df <- prep_marker(d$df, mk)          # mk_val, taxa_class, PS (PC1 + PC2)
    r1 <- fit_glm_safe("s_area_bin ~ mk_val", df)
    r2 <- fit_glm_safe(paste("s_area_bin ~ mk_val + PS +", cov_formula), df)
    cm <- CM[[pl]]
    if (mk %in% colnames(cm)) {
      df$mk_val <- cm[rownames(df), mk]
      r3 <- fit_glm_safe(paste("s_area_bin ~ mk_val +", cov_formula, "+ PC1 + PC2"), df)
    } else r3 <- c(NA_real_, NA_real_, NA_real_)
    rows[[length(rows) + 1]] <- data.frame(platform = pl, marker = mk,
      LR_unadj_pval = r1[1], LR_unadj_beta = r1[2], LR_unadj_se = r1[3], LR_unadj_n = nrow(df),
      LR_PS_pval = r2[1], LR_PS_beta = r2[2], LR_PS_se = r2[3], LR_PS_n = nrow(df),
      LR_COM_pval = r3[1], LR_COM_beta = r3[2], LR_COM_se = r3[3], LR_COM_n = nrow(df),
      stringsAsFactors = FALSE)
  }
}
LR <- do.call(rbind, rows); rownames(LR) <- NULL
write.csv(LR, file.path(c_dir, "LR_recomputed.csv"), row.names = FALSE)

chk <- list()
for (pl in names(D)) {
  old <- read.csv(D[[pl]]$results, stringsAsFactors = FALSE)
  x <- LR[LR$platform == pl, ]; o <- old[match(x$marker, old$marker), ]
  for (m in c("LR_unadj", "LR_PS", "LR_COM")) for (s in c("beta", "pval")) {
    a <- x[[paste0(m, "_", s)]]; b <- o[[paste0(m, "_", s)]]
    chk[[length(chk) + 1]] <- data.frame(platform = pl, method = m, stat = s,
      max_abs_diff = max(abs(a - b), na.rm = TRUE), max_rel_diff = max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE),
      n_sig_changed = if (s == "pval") sum((a < .05) != (b < .05), na.rm = TRUE) else NA,
      stringsAsFactors = FALSE)
  }
}
chk <- do.call(rbind, chk)
write.csv(chk, file.path(c_dir, "LR_recompute_check.csv"), row.names = FALSE)
print(chk)
