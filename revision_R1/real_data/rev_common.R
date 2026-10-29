## ============================================================
##  rev_common.R  -- shared helpers for the BIB-26-1896 revision
##  (real data only; sourced by rev*_*.R scripts)
##
##  Matching / GEE helpers are copied from realdata_analysis.R so that
##  the default settings reproduce the published real-data results:
##    - round 1: nearest, replace = FALSE, caliper (SD of PS)
##    - round 2: unmatched treated + all controls, replace = TRUE
##    - weights: reused control gets 1/k, then mean-normalised   ("inv")
##  Extra options (sensitivity only):
##    - weight_type = "freq": w = m / mean(m), m = number of pairs a
##      subject appears in (definition used in the MatchOmics package
##      and in the simulation generator)
##    - round2_replace = FALSE: round 2 without replacement (package)
##    - m.order: passed to matchit() when not NULL
##
##  Paths: set env vars before running (ASCII-only script on purpose,
##  R 3.6 on Windows cannot parse UTF-8 Korean literals reliably)
##    REV_RAW  = <omics root>/RawData_FullSamples
##    REV_BASE = <project root>/<results_20260517 folder>
##    REV_LIB  = optional extra library path (prepended to .libPaths())
## ============================================================

if (nzchar(Sys.getenv("REV_LIB"))) .libPaths(c(Sys.getenv("REV_LIB"), .libPaths()))

suppressPackageStartupMessages({
  library(MatchIt)
  library(geeM)
  library(readxl)
})

raw_dir    <- Sys.getenv("REV_RAW")
base_dir   <- Sys.getenv("REV_BASE")
if (!nzchar(raw_dir) || !nzchar(base_dir)) stop("Set REV_RAW and REV_BASE env vars")
raw_dir_bc <- file.path(raw_dir, "BeforeBatchCorrection")
res_dir    <- file.path(base_dir, "realdata")
rev_dir    <- file.path(res_dir, "revision")
dir.create(rev_dir, recursive = TRUE, showWarnings = FALSE)

covariates  <- c("s_age", "s_gender", "BMI", "smoke", "fix_dur")
cov_formula <- paste(covariates, collapse = " + ")
calipers    <- list(nocal = NULL, cal0.3 = 0.3, cal0.2 = 0.2, cal0.1 = 0.1, cal0.05 = 0.05)
matchform   <- formula(taxa_class ~ PC1 + PC2)

## 26 concordant markers (original submission, raw p < 0.05 in LR_PS, LR_COM, EXW, EXWO)
conc_prot <- c("K2C1","CXCL7","CO8G","APOL1","RET4","LRP1","K1C10","FINC",
               "HERC1","LYAM1","C1R","HRG","IBP5")
conc_meta <- c("Pyroglutamic_acid","C8_CN","LAC","G6P_F6P","C5_CN","C12_CN",
               "Glycine","AKG","Alanine","PE_18_1_Lyso","PE_16_0_18_1",
               "PE_16_0","PE_17_0")

## ---------------- union-find (identical to realdata_analysis.R) ----------
build_cluster_ids <- function(pairs) {
  all_ids <- unique(c(pairs$treat, pairs$control))
  parent  <- setNames(all_ids, all_ids)
  find <- function(x) {
    while (parent[[x]] != x) { parent[[x]] <<- parent[[parent[[x]]]]; x <- parent[[x]] }
    x
  }
  for (i in seq_len(nrow(pairs))) {
    px <- find(pairs$treat[i]); py <- find(pairs$control[i])
    if (px != py) parent[[px]] <- py
  }
  setNames(as.integer(factor(sapply(all_ids, find))), all_ids)
}

do_matchit <- function(dat, caliper, replace, m.order) {
  args <- list(formula = matchform, data = dat, method = "nearest",
               replace = replace, caliper = caliper)
  if (!is.null(m.order)) args$m.order <- m.order
  tryCatch(suppressWarnings(do.call(matchit, args)), error = function(e) NULL)
}

## ---------------- two-round wrep -----------------------------------------
## Returns matched data with columns pair_id, wt (chosen weight), w_inv,
## w_freq, m (pair count) and attribute "info" (matching bookkeeping).
run_wrep <- function(dat, caliper, weight_type = if (SPEC_C) "freq" else "inv",
                     round2_replace = !SPEC_C, m.order = NULL) {
  m1 <- do_matchit(dat, caliper, FALSE, m.order)
  if (is.null(m1)) return(NULL)

  mm1      <- m1$match.matrix
  matched1 <- rownames(mm1)[!is.na(mm1[, 1])]
  all_trt  <- rownames(dat)[dat$taxa_class == 1]
  unmatched <- setdiff(all_trt, matched1)

  pairs1 <- data.frame(treat = matched1, control = mm1[matched1, 1],
                       stringsAsFactors = FALSE)
  pairs2 <- NULL
  if (length(unmatched) > 0) {
    ctrl_ids <- rownames(dat)[dat$taxa_class == 0]
    sub2 <- dat[c(unmatched, ctrl_ids), ]
    m2 <- do_matchit(sub2, caliper, round2_replace, m.order)
    if (!is.null(m2)) {
      mm2   <- m2$match.matrix
      done2 <- rownames(mm2)[!is.na(mm2[, 1])]
      pairs2 <- data.frame(treat = done2, control = mm2[done2, 1],
                           stringsAsFactors = FALSE)
    }
  }
  ap <- rbind(pairs1, pairs2)
  if (is.null(ap) || nrow(ap) == 0) return(NULL)

  clust   <- build_cluster_ids(ap)
  cf      <- table(ap$control)
  all_ids <- unique(c(ap$treat, ap$control))
  sub     <- dat[all_ids, ]
  sub$pair_id <- clust[rownames(sub)]

  ## "inv": published real-data weights (reused control -> 1/k)
  w_inv <- rep(1.0, nrow(sub)); names(w_inv) <- rownames(sub)
  reused <- names(cf)[cf > 1]; reused <- reused[reused %in% rownames(sub)]
  w_inv[reused] <- 1 / as.numeric(cf[reused])
  sub$w_inv <- w_inv / mean(w_inv)

  ## "freq": w = m / mean(m)  (package / simulation definition)
  m_all <- table(c(ap$treat, ap$control))
  sub$m      <- as.numeric(m_all[rownames(sub)])
  sub$w_freq <- sub$m / mean(sub$m)

  sub$wt <- if (weight_type == "freq") sub$w_freq else sub$w_inv
  attr(sub, "info") <- list(
    n_trt_total   = length(all_trt),
    n_ctrl_total  = sum(dat$taxa_class == 0),
    n_trt_r1      = length(matched1),
    n_trt_r2      = if (is.null(pairs2)) 0L else nrow(pairs2),
    n_pairs       = nrow(ap),
    n_ctrl_unique = length(unique(ap$control)),
    n_clusters    = length(unique(clust)))
  sub
}

## ---------------- worep (single round, no replacement) -------------------
run_worep <- function(dat, caliper, m.order = NULL) {
  m <- do_matchit(dat, caliper, FALSE, m.order)
  if (is.null(m)) return(NULL)
  sub <- match.data(m)
  sub$pair_id <- as.integer(sub$subclass)
  sub$wt <- 1.0; sub$w_inv <- 1.0; sub$w_freq <- 1.0; sub$m <- 1
  attr(sub, "info") <- list(
    n_trt_total   = sum(dat$taxa_class == 1),
    n_ctrl_total  = sum(dat$taxa_class == 0),
    n_trt_r1      = sum(sub$taxa_class == 1),
    n_trt_r2      = 0L,
    n_pairs       = sum(sub$taxa_class == 1),
    n_ctrl_unique = sum(sub$taxa_class == 0),
    n_clusters    = length(unique(sub$pair_id)))
  sub
}

## ---------------- GEE (identical to realdata_analysis.R) -----------------
## geeM 0.10.1 BUG: geem() sorts the model frame by id (dat[order(id), ]) but
## builds W <- Diagonal(x = weights) from the UNSORTED weight vector, so when
## rows are not already ordered by cluster the weights are attached to the
## wrong subjects. realdata_analysis.R did not sort -> GEE-EXW (non-constant
## weights) is affected; GEE-EXWO (all weights 1) is not. The simulation code
## and the MatchOmics package sort by cluster first and are not affected.
## sort_id = TRUE applies the fix; default follows env REV_GEE_SORT ("1" = fix),
## otherwise FALSE (= published behaviour).
## REV_SPEC = "C": final specification (package / Methods): run_wrep defaults to
## round 2 without replacement + w = m/mean(m), GEE sorted, marker set and
## concordant lists read from realdata/revision/C_final, outputs -> C_final/*_C.csv
SPEC_C   <- Sys.getenv("REV_SPEC") == "C"
GEE_SORT <- SPEC_C || Sys.getenv("REV_GEE_SORT") == "1"
## REV_CDIR: name of the C output folder under revision/ (default "C_final";
## e.g. "C_final_R43" for the R 4.3.1 reproduction run)
c_dir_name <- Sys.getenv("REV_CDIR", "C_final")
c_dir      <- file.path(rev_dir, c_dir_name)
dir.create(c_dir, showWarnings = FALSE)
out_name <- function(x) {
  if (SPEC_C) file.path(c_dir_name, sub("\\.csv$", "_C.csv", x))
  else if (GEE_SORT) sub("\\.csv$", "_sorted.csv", x) else x
}
fit_gee <- function(sub, corstr = "exchangeable", sort_id = GEE_SORT) {
  if (is.null(sub) || nrow(sub) < 4 || length(unique(sub$pair_id)) < 2)
    return(c(pval = NA_real_, beta = NA_real_, se = NA_real_))
  if (sort_id) sub <- sub[order(sub$pair_id), ]
  fml <- as.formula(paste("s_area_bin ~ mk_val + PS +", cov_formula))
  fit <- tryCatch(
    suppressWarnings(geem(fml, data = sub, family = binomial(),
         corstr = corstr, id = sub$pair_id,
         weights = sub$wt, sandwich = TRUE,
         useP = TRUE, maxit = 25, tol = 1e-3)),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged))
    return(c(pval = NA_real_, beta = NA_real_, se = NA_real_))
  sm <- summary(fit)
  ## unname(): with newer Matrix/R builds diag(fit$var) carries names, which would
  ## turn the element name "se" into "se.mk_val" and break r["se"] lookups
  c(pval = unname(sm$p[2]), beta = unname(fit$beta[2]), se = unname(sqrt(diag(fit$var))[2]))
}

## ---------------- SMD helpers --------------------------------------------
## weighted mean / weighted variance within each group (smd_audit.R formula)
smd_w <- function(x, grp, w) {
  g1 <- grp == 1; g0 <- grp == 0
  wm1 <- sum(w[g1]*x[g1], na.rm=TRUE)/sum(w[g1], na.rm=TRUE)
  wm0 <- sum(w[g0]*x[g0], na.rm=TRUE)/sum(w[g0], na.rm=TRUE)
  wv1 <- sum(w[g1]*(x[g1]-wm1)^2, na.rm=TRUE)/sum(w[g1], na.rm=TRUE)
  wv0 <- sum(w[g0]*(x[g0]-wm0)^2, na.rm=TRUE)/sum(w[g0], na.rm=TRUE)
  sp  <- sqrt((wv1+wv0)/2)
  if (is.na(sp) || sp == 0) return(NA_real_)
  abs(wm1-wm0)/sp
}
## unweighted (smd_audit.R smd_uw; sample variances)
smd_u <- function(x, grp) {
  g1 <- x[grp == 1]; g0 <- x[grp == 0]
  sp <- sqrt((var(g1, na.rm=TRUE) + var(g0, na.rm=TRUE)) / 2)
  if (is.na(sp) || sp == 0) return(NA_real_)
  abs(mean(g1, na.rm=TRUE) - mean(g0, na.rm=TRUE)) / sp
}
## weighted mean difference / fixed denominator (cobalt-style: pre-match pooled SD)
smd_wfix <- function(x, grp, w, sd_pre) {
  g1 <- grp == 1; g0 <- grp == 0
  wm1 <- sum(w[g1]*x[g1], na.rm=TRUE)/sum(w[g1], na.rm=TRUE)
  wm0 <- sum(w[g0]*x[g0], na.rm=TRUE)/sum(w[g0], na.rm=TRUE)
  if (is.na(sd_pre) || sd_pre == 0) return(NA_real_)
  abs(wm1-wm0)/sd_pre
}
pre_sd <- function(df) sapply(covariates, function(cv) {
  x <- as.numeric(df[[cv]]); g <- df$taxa_class
  sqrt((var(x[g==1], na.rm=TRUE) + var(x[g==0], na.rm=TRUE)) / 2)
})
## max over covariates for each weighting scheme
smd_set <- function(sub, sdp) {
  grp <- sub$taxa_class
  one <- function(f) max(sapply(covariates, function(cv) f(as.numeric(sub[[cv]]), cv)), na.rm = TRUE)
  c(max_uw      = one(function(x, cv) smd_u(x, grp)),
    max_w_inv   = one(function(x, cv) if (all(sub$w_inv == 1)) smd_u(x, grp) else smd_w(x, grp, sub$w_inv)),
    max_w_freq  = one(function(x, cv) smd_w(x, grp, sub$w_freq)),
    max_wfix_freq = one(function(x, cv) smd_wfix(x, grp, sub$w_freq, sdp[cv])))
}

kish <- function(w) sum(w)^2 / sum(w^2)

## ---------------- per-marker setup ---------------------------------------
## exposure: "median" (>= median vs < median, as published) or
##           "tertile" (top third = 1 vs bottom third = 0, middle dropped)
##           "quartile" (top quarter vs bottom quarter, middle half dropped)
## PS is re-estimated (PC1 + PC2) within the retained subjects.
prep_marker <- function(df, mk, exposure = "median") {
  df$mk_val <- df[[mk]]
  if (exposure == "median") {
    df$taxa_class <- as.integer(df$mk_val >= median(df$mk_val, na.rm = TRUE))
  } else if (exposure %in% c("tertile", "quartile")) {
    pr <- if (exposure == "tertile") c(1/3, 2/3) else c(1/4, 3/4)
    q <- quantile(df$mk_val, pr, na.rm = TRUE)
    df$taxa_class <- ifelse(df$mk_val >= q[2], 1L, ifelse(df$mk_val <= q[1], 0L, NA_integer_))
    df <- df[!is.na(df$taxa_class), ]
  }
  ps_fit <- tryCatch(glm(matchform, data = df, family = "binomial"), error = function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  df$PS <- fitted(ps_fit)
  df
}

## ---------------- data loading (identical to realdata_analysis.R) --------
load_realdata <- function() {
  ## NOTE: the local copy of this xlsx got two pivot sheets ("Sheet1", "Sheet2")
  ## inserted in front on 2026-06-09; the data sheet is "Sheet 1".
  cov_xlsx <- file.path(raw_dir, "Cleaned_Data_ADDeGFR_imputed_jm_fixdur.xlsx")
  cov_sh   <- if ("Sheet 1" %in% excel_sheets(cov_xlsx)) "Sheet 1" else 1
  cov_raw  <- as.data.frame(read_excel(cov_xlsx, sheet = cov_sh))
  names(cov_raw)[tolower(names(cov_raw)) == "id"] <- "ID"
  cov_dat <- cov_raw[, c("ID", covariates, "s_area")]
  cov_dat$ID <- gsub("-", "_", cov_dat$ID)
  cov_dat <- cov_dat[complete.cases(cov_dat), ]
  cov_dat$s_area_bin <- as.integer(cov_dat$s_area == 1)

  prot_raw <- read.delim(file.path(raw_dir_bc, "20250410_log2_norm_Proteomics_All.txt"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  prot_mk <- setdiff(names(prot_raw), c("ID", "Area", "Exposure", "Batch"))
  prot_mk <- prot_mk[!grepl("^IST00", prot_mk)]
  prot_df <- merge(prot_raw[, c("ID", "Batch", prot_mk)], cov_dat, by = "ID")
  rownames(prot_df) <- prot_df$ID
  pca <- prcomp(as.matrix(prot_df[, prot_mk]), scale. = TRUE)
  prot_df$PC1 <- pca$x[, 1]; prot_df$PC2 <- pca$x[, 2]

  meta_raw <- read.delim(file.path(raw_dir_bc, "20250428_log2_norm_Metabolite_All.txt"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  meta_mk <- setdiff(names(meta_raw), "ID")
  meta_df <- merge(meta_raw, cov_dat, by = "ID")
  rownames(meta_df) <- meta_df$ID
  pca <- prcomp(as.matrix(meta_df[, meta_mk]), scale. = TRUE)
  meta_df$PC1 <- pca$x[, 1]; meta_df$PC2 <- pca$x[, 2]

  cat(sprintf("covariates: %d (area1=%d) | protein: %d x %d | metabolite: %d x %d\n",
              nrow(cov_dat), sum(cov_dat$s_area_bin), nrow(prot_df), length(prot_mk),
              nrow(meta_df), length(meta_mk)))

  out <- list(Protein    = list(df = prot_df, mk = prot_mk, conc = conc_prot,
                         final = file.path(res_dir, "protein_final.csv"),
                         results = file.path(res_dir, "protein_results.csv")),
       Metabolite = list(df = meta_df, mk = meta_mk, conc = conc_meta,
                         final = file.path(res_dir, "metabolite_final.csv"),
                         results = file.path(res_dir, "metabolite_results.csv")))
  if (SPEC_C) {   # C marker set (C-weighted SMD filter) and C raw-p concordant set
    cd <- c_dir
    cc <- read.csv(file.path(cd, "C_concordant_raw.csv"), stringsAsFactors = FALSE)
    out$Protein$final    <- file.path(cd, "protein_final_C.csv")
    out$Metabolite$final <- file.path(cd, "metabolite_final_C.csv")
    out$Protein$conc     <- cc$marker[cc$platform == "Protein"]
    out$Metabolite$conc  <- cc$marker[cc$platform == "Metabolite"]
    cat("SPEC C: concordant", length(out$Protein$conc), "+", length(out$Metabolite$conc), "\n")
  }
  out
}
