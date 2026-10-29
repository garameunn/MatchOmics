## revC_supp.R -- supplementary tables for specification C
##  1. STable_SMDfail_C.csv : markers failing the C weighted-SMD filter
##     (a) the 7 published-concordant markers dropped by the filter
##     (b) any other filter-failing marker with raw p < 0.05 in all 4 methods
##     per-covariate weighted SMD (C matching, cal0.3), round-2 rescue,
##     5 methods (C definition) beta + exact p
##  2. STable_order_sensitivity_C.csv : default (m.order = "largest") vs 20
##     random greedy orders for the 22 C-concordant markers (from rev05 _C runs)
## Usage: Rscript revC_supp.R   (env REV_RAW, REV_BASE, REV_LIB; REV_SPEC unset)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rev_common.R"))
## c_dir comes from rev_common.R (env REV_CDIR, default C_final)
W <- read.csv(file.path(c_dir, "C_all_markers.csv"), stringsAsFactors = FALSE, check.names = FALSE)
L <- read.csv(file.path(c_dir, "C_matching_long.csv"), stringsAsFactors = FALSE)
cl <- read.csv(file.path(c_dir, "C_concordant_raw.csv"), stringsAsFactors = FALSE)

## ---- 1. SMD-fail table ---------------------------------------------------------
dropped7 <- data.frame(platform = c("Protein", rep("Metabolite", 6)),
                       marker = c("LYAM1", "AKG", "Alanine", "G6P_F6P", "Glycine", "PE_16_0", "PE_17_0"),
                       stringsAsFactors = FALSE)
W$smd_pass_C <- !is.na(W$`wrep_cal0.3_smd_max`) & W$`wrep_cal0.3_smd_max` < 0.2
sig4 <- with(W, !is.na(LR_PS_pval + LR_COM_pval + `wrep_cal0.3_pval` + `worep_cal0.3_pval`) &
               LR_PS_pval < .05 & LR_COM_pval < .05 & `wrep_cal0.3_pval` < .05 & `worep_cal0.3_pval` < .05)
key  <- paste(W$platform, W$marker)
k7   <- paste(dropped7$platform, dropped7$marker)
stopifnot(all(k7 %in% key[!W$smd_pass_C]))
extra <- W[!W$smd_pass_C & sig4 & !(key %in% k7), c("platform", "marker")]
cat("other C-filter failures with all 4 methods p<0.05:", nrow(extra), "\n")
if (nrow(extra)) print(extra)
sel <- rbind(data.frame(dropped7, category = "A-concordant, dropped by C SMD filter"),
             if (nrow(extra)) data.frame(extra, category = "not A-concordant; fails C SMD filter; 4 methods p<0.05"))
A_sets <- list(Protein = read.csv(file.path(res_dir, "protein_final.csv"))$marker,
               Metabolite = read.csv(file.path(res_dir, "metabolite_final.csv"))$marker)

D <- load_realdata()
cov_lab <- c(s_age = "age", s_gender = "sex", BMI = "BMI", smoke = "smoking", fix_dur = "fix_dur")
out <- list()
for (i in seq_len(nrow(sel))) {
  pl <- sel$platform[i]; mk <- sel$marker[i]
  df <- prep_marker(D[[pl]]$df, mk)
  sc <- run_wrep(df, 0.3, weight_type = "freq", round2_replace = FALSE)
  s  <- sapply(covariates, function(cv) smd_w(as.numeric(sc[[cv]]), sc$taxa_class, sc$w_freq))
  stopifnot(abs(max(s) - W$`wrep_cal0.3_smd_max`[key == paste(pl, mk)]) < 1e-10)
  lg <- L[L$platform == pl & L$marker == mk & L$caliper == "cal0.3" & L$method == "GEE-EXW", ]
  w  <- W[key == paste(pl, mk), ]
  r  <- data.frame(category = sel$category[i], platform = pl, marker = mk,
                   in_A_SMD_set = mk %in% A_sets[[pl]],
                   C_max_weighted_SMD = round(max(s), 3), max_SMD_covariate = cov_lab[names(which.max(s))],
                   stringsAsFactors = FALSE)
  for (cv in covariates) r[[paste0("SMD_", cov_lab[cv])]] <- round(s[cv], 3)
  r$n_treated <- lg$n_trt_total; r$n_treated_matched_round1 <- lg$n_trt_matched_r1
  r$n_treated_rescued_round2 <- lg$n_trt_rescued_r2
  r$pct_treated_rescued_round2 <- round(100 * lg$n_trt_rescued_r2 / lg$n_trt_total, 1)
  r$pct_round1_unmatched_rescued <- round(100 * lg$n_trt_rescued_r2 / (lg$n_trt_total - lg$n_trt_matched_r1), 1)
  r$n_matched_EXW <- lg$n_matched
  meth <- c(LR_unadj = "LR_unadj", LR_PS = "LR_PS", LR_COM = "LR_COM",
            `GEE-EXW` = "wrep_cal0.3", `GEE-EXWO` = "worep_cal0.3")
  for (nm in names(meth)) {
    r[[paste0(nm, "_beta")]] <- round(w[[paste0(meth[nm], "_beta")]], 3)
    r[[paste0(nm, "_p")]]    <- signif(w[[paste0(meth[nm], "_pval")]], 3)
  }
  r$all4_p05_C <- sig4[key == paste(pl, mk)]
  out[[i]] <- r
}
tab1 <- do.call(rbind, out); rownames(tab1) <- NULL
write.csv(tab1, file.path(c_dir, "STable_SMDfail_C.csv"), row.names = FALSE)
print(tab1[, c("category", "platform", "marker", "C_max_weighted_SMD", "max_SMD_covariate",
               "pct_treated_rescued_round2", "GEE-EXW_p", "all4_p05_C")])

## C-passing markers for reference: round-2 rescue among the C marker set
ref <- L[L$caliper == "cal0.3" & L$method == "GEE-EXW", ]
ref$pass <- paste(ref$platform, ref$marker) %in% key[W$smd_pass_C]
cat("\nmedian % treated rescued in round 2 (C, cal0.3):\n")
print(aggregate(100 * n_trt_rescued_r2 / n_trt_total ~ platform + pass, data = ref, FUN = median))

## ---- 2. order sensitivity table -------------------------------------------------
R <- read.csv(file.path(c_dir, "rev05_order_permutation_runs_C.csv"), stringsAsFactors = FALSE)
border <- c("HERC1", "PE_16_0_18_1", "PE_18_0_18_1")
tab2 <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
  pl <- cl$platform[i]; mk <- cl$marker[i]
  o <- R[R$platform == pl & R$marker == mk & R$scheme == "original", ]
  x <- R[R$platform == pl & R$marker == mk & R$scheme == "m.order_random", ]
  ok <- !is.na(x$pval)
  data.frame(platform = pl, marker = mk, borderline = mk %in% border,
             beta_default = round(o$beta, 3), p_default = signif(o$pval, 3),
             n_random_orders = nrow(x), n_not_converged = sum(!ok),
             n_sig_p05 = sum(x$pval[ok] < 0.05), n_converged = sum(ok),
             p_min = signif(min(x$pval[ok]), 3), p_median = signif(median(x$pval[ok]), 3),
             p_max = signif(max(x$pval[ok]), 3),
             beta_min = round(min(x$beta[ok]), 3), beta_max = round(max(x$beta[ok]), 3),
             beta_sign_consistent = all(sign(x$beta[ok]) == sign(o$beta)),
             n_matched_default = o$n_matched,
             n_matched_range = sprintf("%d-%d", min(x$n_matched, na.rm = TRUE), max(x$n_matched, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
tab2 <- tab2[order(tab2$platform != "Protein", tab2$marker), ]
write.csv(tab2, file.path(c_dir, "STable_order_sensitivity_C.csv"), row.names = FALSE)
print(tab2[, c("platform", "marker", "borderline", "p_default", "n_sig_p05", "n_converged", "n_not_converged",
               "p_min", "p_median", "p_max", "beta_sign_consistent")])
