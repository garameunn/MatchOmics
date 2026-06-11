## ============================================================
##  03_Type1Error_Calculating.R  [SERVER VERSION]
##  Usage:
##    Rscript 03_Type1Error_Calculating_server.R --cores 20
##    Rscript 03_Type1Error_Calculating_server.R --cores 20 --scenario S1
##    Rscript 03_Type1Error_Calculating_server.R --cores 20 --scenario S1,S2
##    Rscript 03_Type1Error_Calculating_server.R --test
## ============================================================

# ── arguments ───────────────────────────────────────────────
options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args

mc_cores <- if (TEST_MODE) 1L else {
  cores_idx <- which(args == "--cores")
  if (length(cores_idx) > 0 && length(args) > cores_idx)
    as.integer(args[cores_idx + 1])
  else 6L
}

# --scenario S1 또는 --scenario S1,S2,S3 형식
scenario_arg <- {
  idx <- which(args == "--scenario")
  if (length(idx) > 0 && length(args) > idx)
    strsplit(args[idx + 1], ",")[[1]]
  else NULL
}

message(sprintf("Cores: %d | TEST_MODE: %s | Scenarios: %s",
                mc_cores, TEST_MODE,
                if (is.null(scenario_arg)) "all" else paste(scenario_arg, collapse=",")))

# ── paths ───────────────────────────────────────────────────
base_dir <- "/home2/nekim/Matching2026"
code_dir <- file.path(base_dir, "code")
data_dir <- file.path(base_dir, "data")
res_dir  <- file.path(base_dir, "results", "T1")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ── libraries & utils ───────────────────────────────────────
source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "99_utils.R"))
source(file.path(code_dir, "02_simulation_generator.R"))

load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "psdifflist.RData"))
load(file.path(data_dir, "seed.RData"))
bat.dat <- read.csv(file.path(data_dir, "batchinfo_full.csv"), header = TRUE)

# ── parallel wrapper ─────────────────────────────────────────
mc_apply <- function(X, FUN, ...) {
  if (mc_cores == 1L) lapply(X, FUN, ...)
  else parallel::mclapply(X, FUN, ..., mc.cores = mc_cores)
}

# ── settings ────────────────────────────────────────────────
nsamchar  <- c(1000, 500, 200)
scenarios <- c("S1", "S2", "S3")
rates     <- c(0.1, 0.2, 0.3)
k_s3      <- c(1, 8, 61, 77)
b         <- 0
omics     <- "Metagenomics"
# S1: inv only (comp excluded)
# S2, S3: inv

if (TEST_MODE) {
  nsamchar  <- 200
  scenarios <- "S1"
  rates     <- 0.2
  rep_n     <- 10
  dimm_max  <- 5
} else {
  rep_n    <- 5000
  dimm_max <- NULL
  # --scenario 인자로 필터링
  if (!is.null(scenario_arg))
    scenarios <- intersect(scenarios, scenario_arg)
}
message(sprintf("Running scenarios: %s", paste(scenarios, collapse=", ")))

# ── calipers ────────────────────────────────────────────────
calipers  <- list(nocal=NULL, cal0.3=0.3, cal0.2=0.2, cal0.1=0.1, cal0.05=0.05)
cal_names <- names(calipers)

# ── column name builder ──────────────────────────────────────
make_colnames <- function() {
  sig_base  <- c("sig_LR_unadj", "sig_LR_PS", "sig_LR_Combat")
  coef_base <- c("coef_LR_unadj", "coef_LR_PS", "coef_LR_Combat")
  se_base   <- c("se_LR_unadj",  "se_LR_PS",  "se_LR_Combat")

  wrep_sig   <- as.vector(sapply(cal_names, function(cn)
    c(paste0("sig_wrep_",  cn, "_ind"), paste0("sig_wrep_",  cn, "_ex"))))
  worep_sig  <- as.vector(sapply(cal_names, function(cn)
    c(paste0("sig_worep_", cn, "_ind"), paste0("sig_worep_", cn, "_ex"))))
  wrep_coef  <- gsub("^sig_", "coef_", wrep_sig)
  worep_coef <- gsub("^sig_", "coef_", worep_sig)
  wrep_se    <- gsub("^sig_", "se_",   wrep_sig)
  worep_se   <- gsub("^sig_", "se_",   worep_sig)

  nrow_cols    <- c(paste0("nrow_wrep_",  cal_names),
                    paste0("nrow_worep_", cal_names))
  neff_cols    <- c(paste0("neff_wrep_",  cal_names),
                    paste0("neff_worep_", cal_names))
  meanrep_cols <- paste0("meanrep_wrep_", cal_names)

  list(
    main = c(sig_base,  wrep_sig,  worep_sig,
             coef_base, wrep_coef, worep_coef,
             nrow_cols, neff_cols, meanrep_cols),
    se   = c(se_base, wrep_se, worep_se)
  )
}
COLNAMES <- make_colnames()

# ── SMD helper ───────────────────────────────────────────────
smd_between_groups <- function(df, covnames) {
  sapply(covnames, function(v) {
    x <- as.numeric(df[[v]])
    g <- as.integer(df$taxa_class == levels(df$taxa_class)[2])
    m1 <- mean(x[g==1], na.rm=TRUE); m0 <- mean(x[g==0], na.rm=TRUE)
    s1 <- var(x[g==1],  na.rm=TRUE); s0 <- var(x[g==0],  na.rm=TRUE)
    psd <- sqrt((s1+s0)/2)
    if (is.na(psd) || psd==0) NA_real_ else (m1-m0)/psd
  })
}

# ── GEE fitter: returns c(p, beta, se) ──────────────────────
fit_geem <- function(tble, corstr) {
  if (is.null(tble) || nrow(tble) == 0) return(c(1, NA_real_, NA_real_))
  fit <- tryCatch(
    geeM::geem(
      formula  = case ~ x1 + PS + marker,
      data     = tble, family = "binomial",
      corstr   = corstr, id = tble$pair_taxa,
      weights  = tble$mm, sandwich = TRUE,
      useP     = TRUE, maxit = 25, tol = 1e-3
    ), error = function(e) NULL)
  if (is.null(fit)) return(c(1, NA_real_, NA_real_))
  ok <- isTRUE(fit$converged) && abs(fit$beta[4]) < 100
  sm <- if (ok) summary(fit) else NULL
  c(if (ok) sm$p[4]                     else 1,
    if (ok) fit$beta[4]                  else NA_real_,
    if (ok) sqrt(diag(fit$var))[4]       else NA_real_)
}

# ── GLM extractor: returns c(p, beta, se) ───────────────────
get_glm <- function(fit, idx) {
  cs <- coef(summary(fit))
  if (!fit$converged || nrow(cs) < idx) return(c(1, NA_real_, NA_real_))
  c(cs[idx, 4], cs[idx, 1], cs[idx, 2])
}

# ════════════════════════════════════════════════════════════
#  MAIN LOOP
# ════════════════════════════════════════════════════════════
total_loops <- length(nsamchar) * length(scenarios) * length(rates)
loop_count  <- 0L

for (nsamp in nsamchar) {
  targetnum <- which(c(1000, 500, 200) == nsamp)
  otutable  <- otulist[[targetnum]]
  indicator <- indiclist[[targetnum]]
  lbs       <- otulist[[4]][seq_len(nsamp)]
  trsftable <- trnsf(dt=otutable, omics=omics, lb.size=lbs)
  datatable <- otutable
  dimm      <- if (!is.null(dimm_max)) min(dimm_max, nrow(otutable)) else nrow(otutable)

  batord   <- match(colnames(otutable), bat.dat[, 1])
  bat.dat1 <- bat.dat[batord, ]
  rownames(bat.dat1) <- bat.dat1[, 1]
  cmbatable <- sva::ComBat(trsftable,
                            batch       = as.character(bat.dat1$Batch),
                            par.prior   = TRUE,
                            prior.plots = FALSE)

  for (scenario in scenarios) {
    # S1: inv only
    nl <- "inv"
    k  <- if (scenario == "S3") k_s3 else NULL

    for (rate in rates) {
      loop_count <- loop_count + 1L

      # ── 스킵 로직: 계산 전 파일 존재 여부 확인 ────────────
      tag      <- sprintf("NSAMP%d_REP%d_RATE%.1f_%s_%s",
                          nsamp, rep_n, rate, nl, scenario)
      out_file <- file.path(res_dir, paste0(tag, ".RData"))
      if (file.exists(out_file)) {
        message(sprintf("  [SKIP] %s already exists", tag))
        next
      }

      message(sprintf("[%d/%d] N=%d | %s | rate=%.1f | cores=%d",
                      loop_count, total_loops, nsamp, scenario, rate, mc_cores))

      set.seed(1)

      pvalues <- mc_apply(seq_len(rep_n), function(h) {

        if (h %% 100 == 1)
          message(sprintf("  [N=%d|%s|rate=%.1f] rep %d/%d",
                          nsamp, scenario, rate, h, rep_n))

        per_taxa <- parallel::mclapply(seq_len(dimm), function(x) {
          tryCatch({

          indx <- dimm * (h-1) + x
          f(indx)

          tbls_w  <- lapply(calipers, function(cl)
            tablemaking_t1(omics, trsftable, datatable, indicator,
                           rate, b, x, cl, scenario, k, nl))
          tbls_wo <- lapply(calipers, function(cl)
            tablemaking_t1_wo(omics, trsftable, datatable, indicator,
                              rate, b, x, cl, scenario, k, nl))

          tble_entire      <- tbls_w[["nocal"]][[2]]
          tble_entire$cmbt <- as.numeric(cmbatable[x, ])

          # GEE: returns c(p, beta, se)
          wrep_res <- setNames(lapply(cal_names, function(cn) {
            tbl <- tbls_w[[cn]][[1]]
            list(ind = fit_geem(tbl, "independence"),
                 ex  = fit_geem(tbl, "exchangeable"))
          }), cal_names)

          worep_res <- setNames(lapply(cal_names, function(cn) {
            tbl <- tbls_wo[[cn]][[1]]
            list(ind = fit_geem(tbl, "independence"),
                 ex  = fit_geem(tbl, "exchangeable"))
          }), cal_names)

          # GLM: returns c(p, beta, se)
          r_unadj  <- get_glm(glm(case~x1+marker,
                                   data=tble_entire, family="binomial"), 3)
          r_PS     <- get_glm(glm(case~x1+PS+marker,
                                   data=tble_entire, family="binomial"), 4)
          r_combat <- get_glm(glm(case~x1+hetero1+hetero2+cmbt,
                                   data=tble_entire, family="binomial"), 5)

          # nrow, neff, meanrep
          nrow_w  <- sapply(cal_names, function(cn) {
            tbl <- tbls_w[[cn]][[1]]; if(is.null(tbl)) 0L else nrow(tbl) })
          nrow_wo <- sapply(cal_names, function(cn) {
            tbl <- tbls_wo[[cn]][[1]]; if(is.null(tbl)) 0L else nrow(tbl) })
          neff_w  <- sapply(cal_names, function(cn) {
            tbl <- tbls_w[[cn]][[1]]
            if(is.null(tbl)||nrow(tbl)==0) NA_real_ else tbl$.neff[1] })
          neff_wo <- sapply(cal_names, function(cn) {
            tbl <- tbls_wo[[cn]][[1]]
            if(is.null(tbl)||nrow(tbl)==0) NA_real_ else tbl$.neff[1] })
          meanrep_w <- sapply(cal_names, function(cn) {
            mr <- tbls_w[[cn]][[3]]
            if(is.null(mr)||length(mr)==0) NA_real_ else mr })

          # SMD
          cov_cols    <- colnames(indicator)
          hetero_cols <- c("hetero1", "hetero2")
          otu_x       <- as.numeric(datatable[x, ])
          pre_df      <- data.frame(
            taxa_class = as.factor(ifelse(otu_x>=median(otu_x),"up","down")),
            indicator)
          smd_pre <- smd_between_groups(pre_df, cov_cols)

          get_smd_post <- function(tbl) {
            if(is.null(tbl)||nrow(tbl)==0) return(rep(NA_real_, length(cov_cols)))
            pd <- data.frame(taxa_class=as.factor(tbl$indOTU),
                             tbl[, hetero_cols, drop=FALSE])
            colnames(pd)[-1] <- cov_cols
            smd_between_groups(pd, cov_cols)
          }
          smd_post_wrep  <- setNames(lapply(cal_names, function(cn)
            get_smd_post(tbls_w[[cn]][[1]])),  paste0("wrep_",  cal_names))
          smd_post_worep <- setNames(lapply(cal_names, function(cn)
            get_smd_post(tbls_wo[[cn]][[1]])), paste0("worep_", cal_names))

          # ── assemble main (p + beta) ──────────────────────
          sig_vec  <- c(r_unadj[1], r_PS[1], r_combat[1],
                        unlist(lapply(cal_names, function(cn)
                          c(wrep_res[[cn]]$ind[1], wrep_res[[cn]]$ex[1]))),
                        unlist(lapply(cal_names, function(cn)
                          c(worep_res[[cn]]$ind[1], worep_res[[cn]]$ex[1]))))
          coef_vec <- c(r_unadj[2], r_PS[2], r_combat[2],
                        unlist(lapply(cal_names, function(cn)
                          c(wrep_res[[cn]]$ind[2], wrep_res[[cn]]$ex[2]))),
                        unlist(lapply(cal_names, function(cn)
                          c(worep_res[[cn]]$ind[2], worep_res[[cn]]$ex[2]))))

          # ── assemble SE ───────────────────────────────────
          se_vec <- c(r_unadj[3], r_PS[3], r_combat[3],
                      unlist(lapply(cal_names, function(cn)
                        c(wrep_res[[cn]]$ind[3], wrep_res[[cn]]$ex[3]))),
                      unlist(lapply(cal_names, function(cn)
                        c(worep_res[[cn]]$ind[3], worep_res[[cn]]$ex[3]))))

          out <- c(sig_vec, coef_vec, nrow_w, nrow_wo, neff_w, neff_wo, meanrep_w)
          attr(out, "se")       <- se_vec
          attr(out, "smd_pre")  <- smd_pre
          attr(out, "smd_post") <- c(smd_post_wrep, smd_post_worep)
          out

          }, error = function(e) {
            # 에러 시 NA 벡터 반환
            out <- rep(NA_real_, length(COLNAMES$main))
            attr(out, "se")       <- rep(NA_real_, length(COLNAMES$se))
            attr(out, "smd_pre")  <- rep(NA_real_, length(colnames(indicator)))
            attr(out, "smd_post") <- NULL
            out
          })
        }, mc.cores=3)  # per_taxa
     })  # pvalues

      # ── combine ───────────────────────────────────────────
      pvalues_comb <- do.call(rbind,
        lapply(pvalues, function(h_res)
          do.call(rbind, lapply(h_res, function(x_res) {
            v <- as.numeric(x_res)
            if (length(v) != length(COLNAMES$main)) v <- rep(NA_real_, length(COLNAMES$main))
            matrix(v, nrow=1)
          }))))
      colnames(pvalues_comb) <- COLNAMES$main

      # SE matrix (separate object, same row structure)
      pvalues_se <- do.call(rbind,
        lapply(pvalues, function(h_res)
          do.call(rbind, lapply(h_res, function(x_res) {
            v <- as.numeric(attr(x_res, "se"))
            if (length(v) != length(COLNAMES$se)) v <- rep(NA_real_, length(COLNAMES$se))
            matrix(v, nrow=1)
          }))))
      colnames(pvalues_se) <- COLNAMES$se

      # SMD list
      smd_list <- lapply(pvalues, function(h_res)
        lapply(h_res, function(x_res)
          list(pre  = attr(x_res, "smd_pre"),
               post = attr(x_res, "smd_post"))))

      # ── save ──────────────────────────────────────────────
      save(pvalues_comb, pvalues_se,
           file = file.path(res_dir, paste0(tag, ".RData")))
      save(smd_list,
           file = file.path(res_dir, paste0(tag, "_SMD.RData")))
      message(sprintf("  [SAVED] %s | ncol_main=%d | ncol_se=%d",
                      tag, ncol(pvalues_comb), ncol(pvalues_se)))

    }  # rate
  }  # scenario
}  # nsamp

message("03_Type1Error_Calculating.R complete.")
