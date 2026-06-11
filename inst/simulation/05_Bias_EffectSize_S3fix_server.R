## ============================================================
##  05_Bias_EffectSize_S3fix_server.R  [SERVER VERSION]
##  Usage:
##    Rscript 05_Bias_EffectSize_S3fix_server.R
##    Rscript 05_Bias_EffectSize_S3fix_server.R --cores 20
##    Rscript 05_Bias_EffectSize_S3fix_server.R --test
##  Purpose:
##  Recomputes EffectSizeList for S3 only, using the S3fix algorithm:
##    - Dynamic k selection (select_k_dynamic, bottom 5% corr taxa)
##    - Mean-centering instead of lm residual orthogonalization
##    - Updated two_round_match (sourced from 02_simulation_generator_S3fix.R)
##  Output: EffectSizeList_S3fix.RData + EffectSizeList_S3fix.csv
## ============================================================

# ── arguments ───────────────────────────────────────────────
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
res_dir  <- file.path(base_dir, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ── libraries & utils ───────────────────────────────────────
source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "99_utils.R"))
# S3fix generator provides: select_k_dynamic, compute_s3_confounder,
# two_round_match (updated), trnsf
source(file.path(code_dir, "02_simulation_generator_S3fix.R"))

load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "psdifflist.RData"))
load(file.path(data_dir, "seed.RData"))

# ── parallel wrapper ─────────────────────────────────────────
mc_apply <- function(X, FUN, ...) {
  if (mc_cores == 1L) lapply(X, FUN, ...)
  else parallel::mclapply(X, FUN, ..., mc.cores = mc_cores)
}

# ── settings ────────────────────────────────────────────────
nsamchar  <- c(1000, 500, 200)
scenarios <- "S3"      # S3fix changes apply to S3 only
rates     <- c(0.1, 0.2, 0.3)
b         <- c(0.001, 0.002, 0.005, 0.01)
omics     <- "Metagenomics"
rep_n     <- 1

if (TEST_MODE) {
  nsamchar  <- 200
  rates     <- 0.2
  dimm_max  <- 5
} else {
  dimm_max <- NULL
}

# ── calipers ────────────────────────────────────────────────
calipers  <- list(nocal = NULL, cal0.3 = 0.3, cal0.2 = 0.2,
                  cal0.1 = 0.1, cal0.05 = 0.05)
cal_names <- names(calipers)
b2_colnames <- c("whole",
                 paste0("wrep_",  cal_names),
                 paste0("worep_", cal_names))

# ════════════════════════════════════════════════════════════
#  S3fix EffectSize helpers
# ════════════════════════════════════════════════════════════

# Compute TT and v23 from mean-centered k-taxa matrix (S3fix structure)
compute_TT_v23_s3fix <- function(trsftable_matched, trsftable, x, sample_ids) {
  s3_res   <- compute_s3_confounder(trsftable, x, sample_ids = sample_ids)
  k_mat_mc <- s3_res$k_mat  # mean-centered k-taxa × matched-samples

  nk <- nrow(k_mat_mc)
  if (nk == 0) return(list(TT = NA_real_, v23 = NA_real_))

  if (nk == 1) {
    TT  <- var(as.numeric(k_mat_mc))
    v23 <- cov(as.numeric(trsftable_matched[x, ]),
               as.numeric(k_mat_mc))
  } else {
    var_each  <- apply(k_mat_mc, 1, var)
    cov_pairs <- do.call(sum, lapply(
      data.frame(combn(seq_len(nk), 2)),
      function(y) cov(k_mat_mc[y[1], ], k_mat_mc[y[2], ])))
    TT  <- sum(var_each) + 2 * cov_pairs
    v23 <- cov(as.numeric(trsftable_matched[x, ]),
               colSums(k_mat_mc))
  }
  list(TT = TT, v23 = v23)
}

# WITH replacement (GEE-EXW / GEE-EXWO caliper variants)
EffectSizeFun_s3fix <- function(omics, trsftable, datatable, indicator,
                                 rate, b, x, cl) {
  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))
  tmr   <- two_round_match(formatch, matchform, caliper = cl)
  pairs <- rbind(tmr$pairs1,
                 if (!is.null(tmr$pairs2)) tmr$pairs2 else NULL)
  if (is.null(pairs) || nrow(pairs) == 0) return(rep(NA_real_, length(b)))

  all_ids <- unique(c(pairs$treated, pairs$control))
  ps_vec  <- tmr$m1$distance
  samdat  <- data.frame(ps    = ps_vec[match(all_ids, names(ps_vec))],
                        value = all_ids, stringsAsFactors = FALSE)
  trsftable_matched <- trsftable[, match(samdat$value, colnames(trsftable))]

  stats <- compute_TT_v23_s3fix(trsftable_matched, trsftable, x,
                                  sample_ids = samdat$value)
  TT  <- stats$TT
  v23 <- stats$v23
  if (is.na(TT)) return(rep(NA_real_, length(b)))

  sds2 <- sd(as.numeric(trsftable_matched[x, ]))

  sapply(b, function(B) {
    b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
           if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                       sqrt(3.03 / (7  * TT))
    Re(polyroot(c(-1.01*B + B*TT*(b3^2), -2*B*b3*v23,
                  (1-B) * sds2^2)))[1]
  })
}

# WITHOUT replacement (GEE-EXWO no-replacement variants)
EffectSizeFun_wo_s3fix <- function(omics, trsftable, datatable, indicator,
                                    rate, b, x, cl) {
  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class = as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse = "+"),
                                  sep = "~"))
  matmat  <- MatchIt::matchit(matchform, formatch,
                               method = "nearest", replace = FALSE, caliper = cl)
  samdat0 <- data.frame(ps    = matmat$distance,
                        subcl = matmat$subclass,
                        value = names(matmat$subclass))
  samdat  <- samdat0[!is.na(samdat0$subcl), ]
  trsftable_matched <- trsftable[, match(samdat$value, colnames(trsftable))]

  stats <- compute_TT_v23_s3fix(trsftable_matched, trsftable, x,
                                  sample_ids = samdat$value)
  TT  <- stats$TT
  v23 <- stats$v23
  if (is.na(TT)) return(rep(NA_real_, length(b)))

  sds2 <- sd(as.numeric(trsftable_matched[x, ]))

  sapply(b, function(B) {
    b3 <- if (omics == "Metagenomics") sqrt(1.01 / (9  * TT)) else
           if (omics == "Proteomics")  sqrt(1.01 / (99 * TT)) else
                                       sqrt(3.03 / (7  * TT))
    Re(polyroot(c(-1.01*B + B*TT*(b3^2), -2*B*b3*v23,
                  (1-B) * sds2^2)))[1]
  })
}

# ════════════════════════════════════════════════════════════
#  MAIN LOOP  (S3 only)
# ════════════════════════════════════════════════════════════
b2list <- list()
total_loops <- length(nsamchar) * length(rates)
loop_count  <- 0L

for (nsamp in nsamchar) {
  targetnum <- which(c(1000, 500, 200) == nsamp)
  otutable  <- otulist[[targetnum]]
  indicator <- indiclist[[targetnum]]
  psdiff    <- psdifflist[[targetnum]]
  lbs       <- otulist[[4]][seq_len(nsamp)]
  trsftable <- trnsf(dt = otutable, omics = omics, lb.size = lbs)
  datatable <- otutable
  dimm      <- if (!is.null(dimm_max)) min(dimm_max, nrow(otutable)) else nrow(otutable)

  for (rate in rates) {
    loop_count <- loop_count + 1L
    message(sprintf("[%d/%d] N=%d | S3 | rate=%.1f",
                    loop_count, total_loops, nsamp, rate))
    set.seed(1)

    effectsize <- mc_apply(seq_len(rep_n), function(h) {
      mc_apply(seq_len(dimm), function(x) {
        indx <- dimm * (h - 1) + x
        f(indx)

        b2_whole <- tryCatch(
          EffectSizeFun_s3fix(omics, trsftable, datatable, indicator,
                              rate, b, x, NULL),
          error = function(e) {
            if (x == 1) message("  [ERR] s3fix whole x=1: ", conditionMessage(e))
            rep(NA_real_, length(b))
          })

        b2_wrep <- lapply(calipers, function(cl)
          tryCatch(
            EffectSizeFun_s3fix(omics, trsftable, datatable, indicator,
                                rate, b, x, cl),
            error = function(e) {
              if (x == 1) message("  [ERR] s3fix wrep x=1: ", conditionMessage(e))
              rep(NA_real_, length(b))
            }))

        b2_worep <- lapply(calipers, function(cl)
          tryCatch(
            EffectSizeFun_wo_s3fix(omics, trsftable, datatable, indicator,
                                   rate, b, x, cl),
            error = function(e) {
              if (x == 1) message("  [ERR] s3fix worep x=1: ", conditionMessage(e))
              rep(NA_real_, length(b))
            }))

        out <- cbind(
          whole = b2_whole,
          do.call(cbind, setNames(b2_wrep,  paste0("wrep_",  cal_names))),
          do.call(cbind, setNames(b2_worep, paste0("worep_", cal_names)))
        )
        colnames(out) <- b2_colnames
        rownames(out) <- paste0("h2_", b)
        out
      })
    })

    # ── summary ───────────────────────────────────────────
    tag_rp <- paste("S3", nsamp, rate, sep = "_")

    b2_all <- Reduce("+", effectsize[[1]]) / length(effectsize[[1]])

    top20_idx <- order(psdiff, decreasing = TRUE)[1:20]
    top20_idx <- top20_idx[top20_idx <= dimm]
    b2_top20  <- if (length(top20_idx) > 0)
      Reduce("+", effectsize[[1]][top20_idx]) / length(top20_idx)
    else b2_all

    b2list[[paste0(tag_rp, "_all")]]   <- b2_all
    b2list[[paste0(tag_rp, "_top20")]] <- b2_top20
    message(sprintf("  [DONE] %s", tag_rp))

  }  # rate
}  # nsamp

# ── save RData ──────────────────────────────────────────────
save(b2list, file = file.path(data_dir, "EffectSizeList_S3fix.RData"))
message(sprintf("Saved EffectSizeList_S3fix.RData (%d entries)", length(b2list)))

# ── save CSV ────────────────────────────────────────────────
rows <- lapply(names(b2list), function(nm) {
  parts  <- strsplit(nm, "_")[[1]]
  nsamp  <- as.integer(parts[2])
  rate   <- as.numeric(parts[3])
  subset <- parts[4]
  mat    <- b2list[[nm]]
  df <- as.data.frame(mat)
  df$b_level  <- as.numeric(sub("h2_", "", rownames(mat)))
  df$scenario <- "S3"
  df$N        <- nsamp
  df$rate     <- rate
  df$subset   <- subset
  df
})
csv_long <- do.call(rbind, rows)
write.csv(csv_long, file = file.path(res_dir, "EffectSizeList_S3fix.csv"),
          row.names = FALSE)
message(sprintf("Saved EffectSizeList_S3fix.csv (%d rows)", nrow(csv_long)))
message("05_Bias_EffectSize_S3fix_server.R complete.")
