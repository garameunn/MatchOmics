## ============================================================
##  07_Diagnostics.R  [SERVER VERSION]
##  Usage:
##    Rscript 07_Diagnostics.R
##    Rscript 07_Diagnostics.R --cores 20
##    Rscript 07_Diagnostics.R --test
## ============================================================

args      <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args
mc_cores  <- if (TEST_MODE) 1L else {
  cores_idx <- which(args == "--cores")
  if (length(cores_idx) > 0 && length(args) > cores_idx)
    as.integer(args[cores_idx + 1])
  else 6L
}
message(sprintf("Cores: %d | TEST_MODE: %s", mc_cores, TEST_MODE))

# ── paths ───────────────────────────────────────────────────
base_dir <- "/home2/nekim/Matching2026"
code_dir <- file.path(base_dir, "code")
data_dir <- file.path(base_dir, "data")
t1_dir   <- file.path(base_dir, "results", "T1")
out_dir  <- file.path(base_dir, "results", "Diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── libraries ───────────────────────────────────────────────
source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "99_utils.R"))

library(ggplot2)
library(ggpubr)

# ── settings ────────────────────────────────────────────────
nsamchar  <- c(1000, 500, 200)
scenarios <- c("S1", "S2", "S3")
rates     <- c(0.1, 0.2, 0.3)

if (TEST_MODE) {
  nsamchar  <- 200
  scenarios <- "S1"
  rates     <- 0.2
  rep_t1    <- 10
  n_rep_r2  <- 5
  dimm_diag <- 5
} else {
  rep_t1    <- 5000
  n_rep_r2  <- 50
  dimm_diag <- 20
}

# ════════════════════════════════════════════════════════════
#  SECTION 1: SMD & Love Plot
# ════════════════════════════════════════════════════════════
message("=== Section 1: SMD & Love Plot ===")
smd_all     <- list()
total_loops <- length(nsamchar) * length(scenarios) * length(rates)
loop_count  <- 0L

for (nsamp in nsamchar) {
  for (scenario in scenarios) {
    for (rate in rates) {
      loop_count <- loop_count + 1L
      tag      <- sprintf("NSAMP%d_REP%d_RATE%.1f_inv_%s",
                          nsamp, rep_t1, rate, scenario)
      smd_file <- file.path(t1_dir, paste0(tag, "_SMD.RData"))
      if (!file.exists(smd_file)) {
        message(sprintf("  [%d/%d] SKIP SMD: %s", loop_count, total_loops, tag))
        next
      }
      load(smd_file)   # smd_list

      pre_vals     <- do.call(rbind, lapply(smd_list, function(h)
        do.call(rbind, lapply(h, function(x) abs(x$pre)))))
      smd_pre_mean <- round(colMeans(pre_vals, na.rm=TRUE), 4)

      method_cal_names <- names(smd_list[[1]][[1]]$post)
      smd_df_list <- lapply(method_cal_names, function(cn) {
        post_vals <- do.call(rbind, lapply(smd_list, function(h)
          do.call(rbind, lapply(h, function(x) abs(x$post[[cn]])))))
        data.frame(
          covariate     = colnames(post_vals),
          mean_smd_pre  = smd_pre_mean,
          mean_smd_post = round(colMeans(post_vals, na.rm=TRUE), 4),
          method_cal    = cn,
          method        = ifelse(grepl("^wrep", cn), "wrep", "worep"),
          caliper       = sub("^w[or]*ep_", "", cn),
          key           = tag,
          stringsAsFactors = FALSE
        )
      })
      smd_df <- do.call(rbind, smd_df_list)

      key <- paste(nsamp, scenario, rate, sep="_")
      smd_all[[key]] <- smd_df

      # Love plot per method
      for (meth in c("wrep", "worep")) {
        smd_sub <- smd_df[smd_df$method == meth, ]
        df_pre  <- data.frame(
          covariate = unique(smd_sub$covariate),
          smd       = smd_pre_mean,
          caliper   = "Before matching")
        df_post <- data.frame(
          covariate = smd_sub$covariate,
          smd       = smd_sub$mean_smd_post,
          caliper   = smd_sub$caliper)
        df_long <- rbind(df_pre, df_post)
        df_long$caliper <- factor(df_long$caliper,
                                   levels=c("Before matching",
                                            unique(smd_sub$caliper)))

        p_love <- ggplot(df_long,
                         aes(x=smd, y=covariate,
                             colour=caliper, shape=caliper)) +
          geom_point(size=3) +
          geom_vline(xintercept=0.1, linetype="dashed", colour="grey40") +
          labs(title  = sprintf("Love Plot (%s): N=%d | %s | rate=%.1f",
                                meth, nsamp, scenario, rate),
               x      = "Absolute SMD", y=NULL,
               colour = "Caliper", shape="Caliper") +
          theme_bw() + theme(legend.position="bottom")

        ggsave(
          filename = file.path(out_dir,
                               sprintf("LovePlot_%s_%s.png", key, meth)),
          plot=p_love, width=7, height=4, dpi=150)
      }
      message(sprintf("  [%d/%d] Love plots saved: %s",
                      loop_count, total_loops, key))
    }
  }
}

smd_combined <- do.call(rbind, smd_all)
openxlsx::write.xlsx(smd_combined,
                     file=file.path(out_dir, "SMD_combined.xlsx"),
                     rowNames=FALSE)
message("SMD table saved.")

# ════════════════════════════════════════════════════════════
#  SECTION 2: Neff Summary
# ════════════════════════════════════════════════════════════
message("=== Section 2: Neff Summary ===")
neff_all   <- list()
loop_count <- 0L

for (nsamp in nsamchar) {
  for (scenario in scenarios) {
    for (rate in rates) {
      loop_count <- loop_count + 1L
      tag     <- sprintf("NSAMP%d_REP%d_RATE%.1f_inv_%s",
                         nsamp, rep_t1, rate, scenario)
      t1_file <- file.path(t1_dir, paste0(tag, ".RData"))
      if (!file.exists(t1_file)) next
      load(t1_file)   # pvalues_comb

      neff_cols <- grep("^neff_", colnames(pvalues_comb), value=TRUE)
      nrow_cols <- grep("^nrow_", colnames(pvalues_comb), value=TRUE)

      neff_df <- data.frame(
        key      = paste(nsamp, scenario, rate, sep="_"),
        nsamp    = nsamp, scenario=scenario, rate=rate,
        model    = neff_cols,
        mean_neff = round(colMeans(pvalues_comb[,neff_cols], na.rm=TRUE), 1),
        mean_nrow = round(colMeans(pvalues_comb[,nrow_cols], na.rm=TRUE), 1),
        stringsAsFactors = FALSE
      )
      neff_all[[paste(nsamp, scenario, rate, sep="_")]] <- neff_df
    }
  }
}

neff_combined <- do.call(rbind, neff_all)
openxlsx::write.xlsx(neff_combined,
                     file=file.path(out_dir, "Neff_summary.xlsx"),
                     rowNames=FALSE)

# Neff barplot (S1, rate=0.1 대표)
neff_sub <- neff_combined[neff_combined$scenario=="S1" &
                          neff_combined$rate==0.1, ]
if (nrow(neff_sub) > 0) {
  neff_sub$caliper <- gsub("neff_(w[or]*ep)_(.*)", "\\2", neff_sub$model)
  neff_sub$method  <- gsub("neff_(w[or]*ep)_.*",  "\\1", neff_sub$model)
  p_neff <- ggplot(neff_sub, aes(x=caliper, y=mean_neff, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~nsamp, labeller=label_both) +
    labs(title="Mean Neff: S1 | rate=0.1",
         x="Caliper", y="Mean Neff", fill="Method") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(out_dir, "Neff_barplot.png"),
         plot=p_neff, width=8, height=4, dpi=150)
}
message("Neff summary saved.")

# ════════════════════════════════════════════════════════════
#  SECTION 3: Round2 Comparison
# ════════════════════════════════════════════════════════════
message("=== Section 3: Round2 Comparison ===")

source(file.path(code_dir, "02_simulation_generator.R"))
load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "seed.RData"))

mc_apply <- function(X, FUN, ...) {
  if (mc_cores == 1L) lapply(X, FUN, ...)
  else parallel::mclapply(X, FUN, ..., mc.cores = mc_cores)
}

nsamp_rep <- 200
rate_rep  <- 0.2
cl_rep    <- 0.2
otutable  <- otulist[[3]]
indicator <- indiclist[[3]]
lbs       <- otulist[[4]][seq_len(nsamp_rep)]
trsftable <- trnsf(dt=otutable, omics="Metagenomics", lb.size=lbs)
datatable <- otutable

round2_summary <- do.call(rbind,
  mc_apply(seq_len(n_rep_r2), function(h) {
    do.call(rbind, lapply(seq_len(dimm_diag), function(x) {
      indx <- dimm_diag*(h-1)+x; f(indx)
      tble <- tablemaking_t1(
        "Metagenomics", trsftable, datatable, indicator,
        rate_rep, 0, x, cl_rep, "S1", NULL, "inv")
      if (is.null(tble)||is.null(tble[[1]])) return(NULL)
      matched <- tble[[1]]
      n_total  <- nrow(matched)
      n_r1     <- sum(matched$round2==0)
      n_r2     <- sum(matched$round2==1)
      data.frame(
        rep=h, taxa=x,
        n_total=n_total, n_round1=n_r1, n_round2=n_r2,
        pct_round2   = round(n_r2/n_total*100, 1),
        neff         = round(matched$.neff[1], 1),
        case_rate_r1 = round(mean(matched$case[matched$round2==0], na.rm=TRUE), 3),
        case_rate_r2 = round(mean(matched$case[matched$round2==1], na.rm=TRUE), 3),
        weight_mean_r1 = round(mean(matched$mm[matched$round2==0], na.rm=TRUE), 3),
        weight_mean_r2 = round(mean(matched$mm[matched$round2==1], na.rm=TRUE), 3)
      )
    }))
  })
)
round2_summary <- round2_summary[!sapply(
  seq_len(nrow(round2_summary)), function(i)
    is.null(round2_summary[i,])), ]

cat(sprintf("\n=== Round2 Summary (N=%d, cal=%.2f) ===\n",
            nsamp_rep, cl_rep))
cat(sprintf("Mean %% round2:        %.1f%%\n",
            mean(round2_summary$pct_round2, na.rm=TRUE)))
cat(sprintf("Mean case rate r1:    %.3f\n",
            mean(round2_summary$case_rate_r1, na.rm=TRUE)))
cat(sprintf("Mean case rate r2:    %.3f\n",
            mean(round2_summary$case_rate_r2, na.rm=TRUE)))
cat(sprintf("Mean weight r1:       %.3f\n",
            mean(round2_summary$weight_mean_r1, na.rm=TRUE)))
cat(sprintf("Mean weight r2:       %.3f\n",
            mean(round2_summary$weight_mean_r2, na.rm=TRUE)))

openxlsx::write.xlsx(round2_summary,
                     file=file.path(out_dir, "Round2_comparison.xlsx"),
                     rowNames=FALSE)

p_r2 <- ggplot(round2_summary, aes(x=pct_round2)) +
  geom_histogram(bins=15, fill="#3498DB", colour="white") +
  labs(title=sprintf("Round2 sample %% (N=%d, cal=%.2f)", nsamp_rep, cl_rep),
       x="% Round2 samples", y="Count") +
  theme_bw()
ggsave(file.path(out_dir, "Round2_pct_hist.png"),
       plot=p_r2, width=6, height=4, dpi=150)

message("Round2 diagnostics saved.")
message("07_Diagnostics.R complete.")
message(sprintf("Output: %s", out_dir))
