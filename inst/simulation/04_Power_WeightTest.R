## ============================================================
##  04_Power_WeightTest.R  [검증용 - 가중치 버그 비교]
##
##  목적: mm = m/sum(m) [버그] vs mm = m/mean(m) [수정]
##        power 차이가 얼마나 나는지 빠르게 확인
##
##  Usage:
##    /usr/local/R-4.4/bin/Rscript 04_Power_WeightTest.R
##    /usr/local/R-4.4/bin/Rscript 04_Power_WeightTest.R --cores 3
## ============================================================

options(warn = -1)
args     <- commandArgs(trailingOnly = TRUE)
mc_cores <- {
  idx <- which(args == "--cores")
  if (length(idx) > 0 && length(args) > idx) as.integer(args[idx + 1]) else 2L
}
message(sprintf("Cores: %d", mc_cores))

# ── paths ───────────────────────────────────────────────────
base_dir <- "/home2/nekim/Matching2026"
code_dir <- file.path(base_dir, "code")
data_dir <- file.path(base_dir, "data")
res_dir  <- file.path(base_dir, "results", "WeightTest")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ── libraries ───────────────────────────────────────────────
source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "99_utils.R"))
source(file.path(code_dir, "02_simulation_generator.R"))

load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "seed.RData"))
bat.dat <- read.csv(file.path(data_dir, "batchinfo_full.csv"), header = TRUE)

# ── fixed settings ───────────────────────────────────────────
nsamp    <- 1000
scenario <- "S2"
rate     <- 0.2
b        <- c(0.001, 0.002, 0.005, 0.01)
omics    <- "Metagenomics"
rep_n    <- 200   # 빠른 검증용
dimm_max <- 30    # taxa 30개만

targetnum <- 1
otutable  <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
lbs       <- otulist[[4]][seq_len(nsamp)]
trsftable <- trnsf(dt = otutable, omics = omics, lb.size = lbs)
datatable <- otutable
dimm      <- min(dimm_max, nrow(otutable))

batord    <- match(colnames(otutable), bat.dat[, 1])
bat.dat1  <- bat.dat[batord, ]
rownames(bat.dat1) <- bat.dat1[, 1]
cmbatable <- sva::ComBat(trsftable,
                          batch       = as.character(bat.dat1$Batch),
                          par.prior   = TRUE,
                          prior.plots = FALSE)

calipers  <- list(nocal = NULL, cal0.2 = 0.2)
cal_names <- names(calipers)
b_labels  <- paste0("h2=", b)

# ── GEE fitter factory: weight_type = "mean"(수정) or "sum"(버그) ──
fit_geem_wt <- function(tble, corstr, weight_type = "mean") {
  if (is.null(tble) || nrow(tble) == 0) return(c(1, NA_real_))
  # 가중치 재계산
  tbpair <- table(tble$pair_taxa)
  m      <- rep(1, nrow(tble))
  m[cumsum(tbpair)] <- tbpair - 1
  mm <- if (weight_type == "sum") m / sum(m) else m / mean(m)

  fit <- tryCatch(
    geeM::geem(
      formula = case ~ x1 + PS + marker,
      data    = tble, family = "binomial",
      corstr  = corstr, id = tble$pair_taxa,
      weights = mm, sandwich = TRUE,
      useP = TRUE, maxit = 25, tol = 1e-3
    ), error = function(e) NULL)
  if (is.null(fit)) return(c(1, NA_real_))
  ok <- isTRUE(fit$converged) && abs(fit$beta[4]) < 100
  c(if (ok) summary(fit)$p[4] else 1,
    if (ok) fit$beta[4]        else NA_real_)
}

get_glm <- function(fit, idx) {
  cs <- coef(summary(fit))
  if (!fit$converged || nrow(cs) < idx) return(c(1, NA_real_))
  c(cs[idx, 4], cs[idx, 1])
}

message(sprintf("=== WeightTest: N=%d | %s | rate=%.1f | rep=%d | dimm=%d ===",
                nsamp, scenario, rate, rep_n, dimm))
set.seed(1)

pvalues_raw <- parallel::mclapply(seq_len(rep_n), function(h) {
  if (h %% 50 == 1)
    message(sprintf("  rep %d/%d", h, rep_n))

  parallel::mclapply(seq_len(dimm), function(x) {
    indx <- dimm * (h - 1) + x
    f(indx)

    tbls_w <- lapply(calipers, function(cl)
      tablemaking_power(omics, trsftable, datatable, indicator,
                        rate, b, x, cl, scenario, NULL, "inv"))

    # entiredata (b-level별)
    entire_list <- lapply(seq_along(b), function(bi) {
      ed <- tbls_w[["nocal"]][[2]][[bi]]
      ed$cmbt <- as.numeric(cmbatable[x, ])
      ed
    })

    lapply(seq_along(b), function(bi) {
      tble_entire <- entire_list[[bi]]

      # GLM (weight 무관)
      r_unadj <- get_glm(glm(case ~ x1 + marker,
                              data = tble_entire, family = "binomial"), 3)
      r_PS    <- get_glm(glm(case ~ x1 + PS + marker,
                              data = tble_entire, family = "binomial"), 4)
      r_cmbt  <- get_glm(glm(case ~ x1 + hetero1 + hetero2 + cmbt,
                              data = tble_entire, family = "binomial"), 5)

      # GEE: mean(수정) vs sum(버그) × nocal/cal0.2 × ind/ex
      gee_res <- do.call(c, lapply(cal_names, function(cn) {
        tbl <- tbls_w[[cn]][[1]][[bi]]
        c(
          mean_ind = fit_geem_wt(tbl, "independence",  "mean"),
          mean_ex  = fit_geem_wt(tbl, "exchangeable",  "mean"),
          sum_ind  = fit_geem_wt(tbl, "independence",  "sum"),
          sum_ex   = fit_geem_wt(tbl, "exchangeable",  "sum")
        )
      }))
      # gee_res: 각 cal × 4(mean_ind p/b, mean_ex p/b, sum_ind p/b, sum_ex p/b) = 16*2=16개
      c(r_unadj, r_PS, r_cmbt, gee_res)
    })
  }, mc.cores = 3)
}, mc.cores = mc_cores)

# ── 집계 ──────────────────────────────────────────────────────
# 컬럼명 구성
glm_cols <- c("p_unadj","b_unadj","p_PS","b_PS","p_cmbt","b_cmbt")
gee_cols <- as.vector(sapply(cal_names, function(cn)
  c(paste0("p_mean_ind_",cn), paste0("b_mean_ind_",cn),
    paste0("p_mean_ex_", cn), paste0("b_mean_ex_", cn),
    paste0("p_sum_ind_", cn), paste0("b_sum_ind_", cn),
    paste0("p_sum_ex_",  cn), paste0("b_sum_ex_",  cn))))
all_cols <- c(glm_cols, gee_cols)

pw_by_b <- lapply(seq_along(b), function(bi) {
  mat <- do.call(rbind, lapply(pvalues_raw, function(h_res)
    do.call(rbind, lapply(h_res, function(x_res)
      matrix(as.numeric(x_res[[bi]]), nrow = 1)))))
  colnames(mat) <- all_cols
  mat
})
names(pw_by_b) <- b_labels

# ── power 요약 출력 ───────────────────────────────────────────
p_cols <- grep("^p_", all_cols, value = TRUE)

cat("\n========================================\n")
cat(" Power @ alpha=0.05  |  N=1000 S2 rate=0.2\n")
cat("========================================\n")
pw_summary <- do.call(rbind, lapply(b_labels, function(bl) {
  mat <- pw_by_b[[bl]]
  pw  <- round(colMeans(mat[, p_cols] < 0.05, na.rm = TRUE), 3)
  data.frame(h2 = bl, t(pw), check.names = FALSE)
}))
print(pw_summary, row.names = FALSE)

cat("\n--- mean(수정) vs sum(버그) 직접 비교 ---\n")
compare <- do.call(rbind, lapply(b_labels, function(bl) {
  mat <- pw_by_b[[bl]]
  data.frame(
    h2           = bl,
    LR_PS        = round(mean(mat[,"p_PS"]    < 0.05, na.rm=TRUE), 3),
    LR_Combat    = round(mean(mat[,"p_cmbt"]  < 0.05, na.rm=TRUE), 3),
    mean_ind_nocal = round(mean(mat[,"p_mean_ind_nocal"] < 0.05, na.rm=TRUE), 3),
    mean_ex_nocal  = round(mean(mat[,"p_mean_ex_nocal"]  < 0.05, na.rm=TRUE), 3),
    sum_ind_nocal  = round(mean(mat[,"p_sum_ind_nocal"]  < 0.05, na.rm=TRUE), 3),
    sum_ex_nocal   = round(mean(mat[,"p_sum_ex_nocal"]   < 0.05, na.rm=TRUE), 3),
    mean_ind_cal02 = round(mean(mat[,"p_mean_ind_cal0.2"] < 0.05, na.rm=TRUE), 3),
    sum_ind_cal02  = round(mean(mat[,"p_sum_ind_cal0.2"]  < 0.05, na.rm=TRUE), 3)
  )
}))
print(compare, row.names = FALSE)

# ── save ──────────────────────────────────────────────────────
save(pw_by_b, pw_summary, compare,
     file = file.path(res_dir, "WeightTest_S2_N1000_rate0.2.RData"))
message(sprintf("\nSaved to %s", res_dir))
message("04_Power_WeightTest.R complete.")
