## ============================================================
##  06_Final_Stat_Organizing.R  [SERVER VERSION]
##  Usage:
##    Rscript 06_Final_Stat_Organizing.R
##    Rscript 06_Final_Stat_Organizing.R --test
## ============================================================

args      <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args
if (TEST_MODE) message("*** TEST MODE ***")

# ── paths ───────────────────────────────────────────────────
base_dir <- "/home2/nekim/Matching2026"
code_dir <- file.path(base_dir, "code")
data_dir <- file.path(base_dir, "data")
t1_dir   <- file.path(base_dir, "results", "T1")
pw_dir   <- file.path(base_dir, "results", "Power")
out_dir  <- file.path(base_dir, "results", "Organized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── libraries & utils ───────────────────────────────────────
source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "99_utils.R"))

# ── settings ────────────────────────────────────────────────
nsamchar  <- c(1000, 500, 200)
scenarios <- c("S1", "S2", "S3")
rates     <- c(0.1, 0.2, 0.3)
b_levels  <- c(0.001, 0.002, 0.005, 0.01)
b_labels  <- paste0("h2_", b_levels)
alpha_levels <- c(0.001, 0.005, 0.01, 0.05, 0.1)

if (TEST_MODE) {
  nsamchar  <- 200
  scenarios <- "S1"
  rates     <- 0.2
  rep_t1    <- 10
  rep_pw    <- 10
} else {
  rep_t1 <- 5000
  rep_pw <- 1000
}

load(file.path(data_dir, "psdifflist.RData"))
load(file.path(data_dir, "EffectSizeList.RData"))

# ── b2_map: coef_col → b2list colname ───────────────────────
cal_names <- c("nocal","cal0.3","cal0.2","cal0.1","cal0.05")
b2_map <- c(
  setNames(rep("whole",        3),
           c("coef_LR_unadj","coef_LR_PS","coef_LR_Combat")),
  setNames(rep(paste0("wrep_",  cal_names), each=2),
           as.vector(sapply(cal_names, function(cn)
             c(paste0("coef_wrep_",cn,"_ind"), paste0("coef_wrep_",cn,"_ex"))))),
  setNames(rep(paste0("worep_", cal_names), each=2),
           as.vector(sapply(cal_names, function(cn)
             c(paste0("coef_worep_",cn,"_ind"), paste0("coef_worep_",cn,"_ex")))))
)

total_loops <- length(nsamchar) * length(scenarios) * length(rates)
loop_count  <- 0L

# ════════════════════════════════════════════════════════════
#  SECTION 1: Type 1 Error
# ════════════════════════════════════════════════════════════
message("=== Section 1: Type 1 Error ===")
t1_results      <- list()
emp_thresh_list <- list()

for (nsamp in nsamchar) {
  norder <- which(c(1000,500,200) == nsamp)
  for (scenario in scenarios) {
    for (rate in rates) {
      loop_count <- loop_count + 1L
      tag     <- sprintf("NSAMP%d_REP%d_RATE%.1f_inv_%s",
                         nsamp, rep_t1, rate, scenario)
      t1_file <- file.path(t1_dir, paste0(tag, ".RData"))
      if (!file.exists(t1_file)) {
        message(sprintf("  [%d/%d] SKIP: %s", loop_count, total_loops, tag))
        next
      }
      load(t1_file)   # pvalues_comb, pvalues_se
      sig_cols <- grep("^sig_", colnames(pvalues_comb), value=TRUE)

      t1_mat <- do.call(rbind, lapply(alpha_levels, function(a)
        round(colMeans(pvalues_comb[, sig_cols] < a, na.rm=TRUE), 5)))
      rownames(t1_mat) <- paste0("alpha_", alpha_levels)

      emp_thresh <- apply(pvalues_comb[, sig_cols], 2,
                          function(x) quantile(x, 0.05, na.rm=TRUE))

      nrow_cols <- grep("^nrow_", colnames(pvalues_comb), value=TRUE)
      neff_cols <- grep("^neff_", colnames(pvalues_comb), value=TRUE)

      key <- paste(nsamp, scenario, rate, sep="_")
      t1_results[[key]] <- list(
        t1         = t1_mat,
        nrow_mean  = round(colMeans(pvalues_comb[,nrow_cols], na.rm=TRUE), 1),
        neff_mean  = round(colMeans(pvalues_comb[,neff_cols], na.rm=TRUE), 1),
        tag        = tag
      )
      emp_thresh_list[[key]] <- emp_thresh
      message(sprintf("  [%d/%d] T1 done: %s", loop_count, total_loops, tag))
    }
  }
}

save(t1_results, emp_thresh_list,
     file=file.path(out_dir, "T1_organized.RData"))

# T1 combined Excel
t1_combined <- do.call(rbind, lapply(names(t1_results), function(key) {
  df       <- as.data.frame(t1_results[[key]]$t1)
  df$alpha <- rownames(df)
  df$key   <- key
  df
}))
openxlsx::write.xlsx(t1_combined,
                     file=file.path(out_dir, "T1_combined.xlsx"),
                     rowNames=FALSE)

# Nrow/Neff combined
nrow_neff_combined <- do.call(rbind, lapply(names(t1_results), function(key) {
  res <- t1_results[[key]]
  data.frame(key=key, stat=c("nrow_mean","neff_mean"),
             rbind(res$nrow_mean, res$neff_mean))
}))
openxlsx::write.xlsx(nrow_neff_combined,
                     file=file.path(out_dir, "Nrow_Neff_combined.xlsx"),
                     rowNames=FALSE)
message("T1 results saved.")

# ════════════════════════════════════════════════════════════
#  SECTION 2: Power & Relative Bias
# ════════════════════════════════════════════════════════════
message("=== Section 2: Power & Relative Bias ===")
pw_results <- list()
loop_count <- 0L

for (nsamp in nsamchar) {
  norder <- which(c(1000,500,200) == nsamp)
  psdiff <- psdifflist[[norder]]

  for (scenario in scenarios) {
    for (rate in rates) {
      loop_count <- loop_count + 1L
      tag     <- sprintf("NSAMP%d_REP%d_RATE%.1f_inv_%s",
                         nsamp, rep_pw, rate, scenario)
      pw_file <- file.path(pw_dir, paste0(tag, ".RData"))
      if (!file.exists(pw_file)) {
        message(sprintf("  [%d/%d] SKIP: %s", loop_count, total_loops, tag))
        next
      }
      load(pw_file)   # pvalues (list of b-level matrices), pvalues_se

      t1_key     <- paste(nsamp, scenario, rate, sep="_")
      emp_thresh <- emp_thresh_list[[t1_key]]

      b2_key_all   <- paste(scenario, nsamp, rate, "all",   sep="_")
      b2_key_top20 <- paste(scenario, nsamp, rate, "top20", sep="_")
      nom_b2_all   <- b2list[[b2_key_all]]
      nom_b2_top20 <- b2list[[b2_key_top20]]

      pw_list <- lapply(seq_along(b_levels), function(bi) {
        mat <- pvalues[[bi]]

        sig_only     <- grep("^sig_",     colnames(mat), value=TRUE)
        coef_only    <- grep("^coef_",    colnames(mat), value=TRUE)
        nrow_cols_m  <- grep("^nrow_",    colnames(mat), value=TRUE)
        neff_cols_m  <- grep("^neff_",    colnames(mat), value=TRUE)
        meanrep_cols_m <- grep("^meanrep_", colnames(mat), value=TRUE)

        pw_nom <- round(colMeans(mat[, sig_only, drop=FALSE] < 0.05,
                                 na.rm=TRUE), 4)
        pw_emp <- if (!is.null(emp_thresh))
          round(sapply(sig_only, function(cn) {
            thr <- if (cn %in% names(emp_thresh)) emp_thresh[cn] else 0.05
            mean(mat[, cn] < thr, na.rm=TRUE)
          }), 4)
        else pw_nom

        beta_mean <- if (length(coef_only) > 0)
          round(colMeans(mat[, coef_only, drop=FALSE], na.rm=TRUE), 5)
        else NULL

        # SE mean from pvalues_se
        se_mat   <- pvalues_se[[bi]]
        se_only  <- colnames(se_mat)
        se_mean  <- if (!is.null(se_mat))
          round(colMeans(se_mat, na.rm=TRUE), 5) else NULL

        nrow_mean    <- if (length(nrow_cols_m)>0)
          round(colMeans(mat[,nrow_cols_m,drop=FALSE], na.rm=TRUE), 1) else NULL
        neff_mean    <- if (length(neff_cols_m)>0)
          round(colMeans(mat[,neff_cols_m,drop=FALSE], na.rm=TRUE), 1) else NULL
        meanrep_mean <- if (length(meanrep_cols_m)>0)
          round(colMeans(mat[,meanrep_cols_m,drop=FALSE], na.rm=TRUE), 2) else NULL

        list(pw_nom=pw_nom, pw_emp=pw_emp,
             beta_mean=beta_mean, se_mean=se_mean,
             nrow_mean=nrow_mean, neff_mean=neff_mean,
             meanrep_mean=meanrep_mean)
      })
      names(pw_list) <- b_labels

      # relative bias
      calc_relb <- function(nom_b2) {
        if (is.null(nom_b2) || is.null(pw_list[[1]]$beta_mean)) return(NULL)
        do.call(rbind, lapply(seq_along(b_levels), function(bi) {
          bm       <- pw_list[[bi]]$beta_mean
          cols     <- names(bm)
          ref_cols <- b2_map[cols]
          valid    <- !is.na(ref_cols)
          ref_vals <- nom_b2[bi, ref_cols[valid]]
          result   <- rep(NA_real_, length(cols))
          result[valid] <- round((bm[valid]-ref_vals)/ref_vals, 4)
          setNames(result, cols)
        }))
      }
      relb_all   <- calc_relb(nom_b2_all)
      relb_top20 <- calc_relb(nom_b2_top20)
      if (!is.null(relb_all))   rownames(relb_all)   <- b_labels
      if (!is.null(relb_top20)) rownames(relb_top20) <- b_labels

      key <- paste(nsamp, scenario, rate, sep="_")
      pw_results[[key]] <- list(tag=tag, pw_list=pw_list,
                                relb_all=relb_all, relb_top20=relb_top20)
      message(sprintf("  [%d/%d] Power done: %s", loop_count, total_loops, tag))
    }
  }
}

save(pw_results, file=file.path(out_dir, "Power_organized.RData"))

# Power combined
pw_combined <- do.call(rbind, lapply(names(pw_results), function(key) {
  res <- pw_results[[key]]
  do.call(rbind, lapply(b_labels, function(bl) {
    pw <- res$pw_list[[bl]]
    rbind(
      data.frame(key=key, h2=bl, type="power_nominal",
                 t(pw$pw_nom), check.names=FALSE),
      data.frame(key=key, h2=bl, type="power_empirical",
                 t(pw$pw_emp), check.names=FALSE)
    )
  }))
}))
openxlsx::write.xlsx(pw_combined,
                     file=file.path(out_dir, "Power_combined.xlsx"),
                     rowNames=FALSE)

# Beta & SE combined
beta_se_combined <- do.call(rbind, lapply(names(pw_results), function(key) {
  res <- pw_results[[key]]
  do.call(rbind, lapply(b_labels, function(bl) {
    pw <- res$pw_list[[bl]]
    df_beta <- if (!is.null(pw$beta_mean))
      data.frame(key=key, h2=bl, stat="beta_mean",
                 t(pw$beta_mean), check.names=FALSE)
    df_se   <- if (!is.null(pw$se_mean))
      data.frame(key=key, h2=bl, stat="se_mean",
                 t(pw$se_mean), check.names=FALSE)
    rbind(df_beta, df_se)
  }))
}))
openxlsx::write.xlsx(beta_se_combined,
                     file=file.path(out_dir, "Beta_SE_combined.xlsx"),
                     rowNames=FALSE)

# RelBias combined
rb_combined <- do.call(rbind, lapply(names(pw_results), function(key) {
  res <- pw_results[[key]]
  rbind(
    if (!is.null(res$relb_all)) {
      rb        <- res$relb_all
      coef_only <- grep("^coef_", colnames(rb), value=TRUE)
      data.frame(key=key, group="all", h2=rownames(rb),
                 rb[, coef_only, drop=FALSE], check.names=FALSE)
    },
    if (!is.null(res$relb_top20)) {
      rb        <- res$relb_top20
      coef_only <- grep("^coef_", colnames(rb), value=TRUE)
      data.frame(key=key, group="top20", h2=rownames(rb),
                 rb[, coef_only, drop=FALSE], check.names=FALSE)
    }
  )
}))
openxlsx::write.xlsx(rb_combined,
                     file=file.path(out_dir, "RelBias_combined.xlsx"),
                     rowNames=FALSE)

message("Power & RelBias results saved.")
message("06_Final_Stat_Organizing.R complete.")
message(sprintf("Output: %s", out_dir))
