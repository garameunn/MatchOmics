## ============================================================
##  05_Bias_EffectSize_Calculating.R  [SERVER VERSION]
##  Usage:
##    Rscript 05_Bias_EffectSize_Calculating.R
##    Rscript 05_Bias_EffectSize_Calculating.R --cores 20
##    Rscript 05_Bias_EffectSize_Calculating.R --test
##  Changes vs local version:
##  1. Server paths (/home2/nekim/Matching2026/)
##  2. S1: inv only (comp excluded)
##  3. --cores argument support
##  4. Loop progress print
##  5. rep=1 always (bias is deterministic across reps)
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
source(file.path(code_dir, "02_simulation_generator.R"))

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
scenarios <- c("S1", "S2", "S3")
rates     <- c(0.1, 0.2, 0.3)
k_s3      <- c(1, 8, 61, 77)
b         <- c(0.001, 0.002, 0.005, 0.01)
omics     <- "Metagenomics"
rep_n     <- 1   # bias is deterministic across reps

if (TEST_MODE) {
  nsamchar  <- 200
  scenarios <- "S1"
  rates     <- 0.2
  dimm_max  <- 5
} else {
  dimm_max <- NULL
}

# ── calipers ────────────────────────────────────────────────
calipers  <- list(nocal=NULL, cal0.3=0.3, cal0.2=0.2, cal0.1=0.1, cal0.05=0.05)
cal_names <- names(calipers)
b2_colnames <- c("whole",
                 paste0("wrep_",  cal_names),
                 paste0("worep_", cal_names))

# ════════════════════════════════════════════════════════════
#  EffectSize functions
# ════════════════════════════════════════════════════════════
EffectSizeFun <- function(omics, trsftable, datatable, indicator,
                           rate, b, x, cl, scenario, k=NULL) {
  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse="+"), sep="~"))
  tmr   <- two_round_match(formatch, matchform, caliper=cl)
  pairs <- rbind(tmr$pairs1,
                 if (!is.null(tmr$pairs2)) tmr$pairs2 else NULL)
  if (is.null(pairs) || nrow(pairs) == 0) return(rep(NA_real_, length(b)))

  all_ids <- unique(c(pairs$treated, pairs$control))
  ps_vec  <- tmr$m1$distance
  samdat  <- data.frame(ps=ps_vec[match(all_ids, names(ps_vec))],
                        value=all_ids, stringsAsFactors=FALSE)
  trsftable_matched <- trsftable[, match(samdat$value, colnames(trsftable))]

  if (scenario=="S1"|scenario=="S2") {
    PPS <- if (scenario=="S1") 1/samdat$ps else
      indicator[match(samdat$value, rownames(indicator)),1] +
      indicator[match(samdat$value, rownames(indicator)),2]
    sds <- c(1, sd(as.numeric(trsftable_matched[x,])), sd(PPS), 1)
    v23 <- cov(as.numeric(trsftable_matched[x,]), PPS)
  } else {
    otu_fg0 <- trsftable[k, match(samdat$value, colnames(trsftable))]
    otu_fg  <- matrix(ncol=ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_fg <- rbind(otu_fg,
        lm(t(otu_fg0)[,ii] ~ trsftable_matched[x,])$residuals)
    otu_fg <- otu_fg[-1,]
    sds    <- c(1, sd(as.numeric(trsftable_matched[x,])),
                apply(otu_fg,1,sd), 1)
    vresid <- do.call(sum, lapply(data.frame(combn(seq_along(k),2)),
      function(y) cov(otu_fg[y[1],], otu_fg[y[2],])))
    TT  <- sum(apply(otu_fg,1,var)) + 2*vresid
    v23 <- cov(as.numeric(trsftable_matched[x,]), apply(otu_fg,2,sum))
  }

  sapply(b, function(B) {
    if (scenario=="S1"|scenario=="S2") {
      b3 <- if (omics=="Metagenomics") sqrt(1.01/(9*sds[3]^2)) else
             if (omics=="Proteomics")  sqrt(1.01/(99*sds[3]^2)) else
                                       sqrt(3.03/(7*sds[3]^2))
      Re(polyroot(c(-1.01*B+B*(b3*sds[3])^2, -2*B*b3*v23,
                    (1-B)*(sds[2]^2))))[1]
    } else {
      b3 <- if (omics=="Metagenomics") sqrt(1.01/(9*TT)) else
             if (omics=="Proteomics")  sqrt(1.01/(99*TT)) else
                                       sqrt(3.03/(7*TT))
      Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23,
                    (1-B)*(sds[2]^2))))[1]
    }
  })
}

EffectSizeFun_wo <- function(omics, trsftable, datatable, indicator,
                              rate, b, x, cl, scenario, k=NULL) {
  otu0       <- data.frame(t(datatable)[, x])
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  formatch   <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform  <- as.formula(paste("taxa_class",
                                  paste(colnames(indicator), collapse="+"), sep="~"))
  matmat  <- MatchIt::matchit(matchform, formatch,
                               method="nearest", replace=FALSE, caliper=cl)
  samdat0 <- data.frame(ps=matmat$distance, subcl=matmat$subclass,
                        value=names(matmat$subclass))
  samdat  <- samdat0[!is.na(samdat0$subcl),]
  trsftable_matched <- trsftable[, match(samdat$value, colnames(trsftable))]

  if (scenario=="S1"|scenario=="S2") {
    PPS <- if (scenario=="S1") 1/samdat$ps else
      indicator[match(samdat$value, rownames(indicator)),1] +
      indicator[match(samdat$value, rownames(indicator)),2]
    sds <- c(1, sd(as.numeric(trsftable_matched[x,])), sd(PPS), 1)
    v23 <- cov(as.numeric(trsftable_matched[x,]), PPS)
  } else {
    otu_fg0 <- trsftable[k, match(samdat$value, colnames(trsftable))]
    otu_fg  <- matrix(ncol=ncol(trsftable_matched))
    for (ii in seq_along(k))
      otu_fg <- rbind(otu_fg,
        lm(t(otu_fg0)[,ii] ~ trsftable_matched[x,])$residuals)
    otu_fg <- otu_fg[-1,]
    sds    <- c(1, sd(as.numeric(trsftable_matched[x,])),
                apply(otu_fg,1,sd), 1)
    vresid <- do.call(sum, lapply(data.frame(combn(seq_along(k),2)),
      function(y) cov(otu_fg[y[1],], otu_fg[y[2],])))
    TT  <- sum(apply(otu_fg,1,var)) + 2*vresid
    v23 <- cov(as.numeric(trsftable_matched[x,]), apply(otu_fg,2,sum))
  }

  sapply(b, function(B) {
    if (scenario=="S1"|scenario=="S2") {
      b3 <- if (omics=="Metagenomics") sqrt(1.01/(9*sds[3]^2)) else
             if (omics=="Proteomics")  sqrt(1.01/(99*sds[3]^2)) else
                                       sqrt(3.03/(7*sds[3]^2))
      Re(polyroot(c(-1.01*B+B*(b3*sds[3])^2, -2*B*b3*v23,
                    (1-B)*(sds[2]^2))))[1]
    } else {
      b3 <- if (omics=="Metagenomics") sqrt(1.01/(9*TT)) else
             if (omics=="Proteomics")  sqrt(1.01/(99*TT)) else
                                       sqrt(3.03/(7*TT))
      Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23,
                    (1-B)*(sds[2]^2))))[1]
    }
  })
}

# ════════════════════════════════════════════════════════════
#  MAIN LOOP
# ════════════════════════════════════════════════════════════
b2list <- list()
total_loops <- length(nsamchar) * length(scenarios) * length(rates)
loop_count  <- 0L

for (nsamp in nsamchar) {
  targetnum <- which(c(1000, 500, 200) == nsamp)
  otutable  <- otulist[[targetnum]]
  indicator <- indiclist[[targetnum]]
  psdiff    <- psdifflist[[targetnum]]
  lbs       <- otulist[[4]][seq_len(nsamp)]
  trsftable <- trnsf(dt=otutable, omics=omics, lb.size=lbs)
  datatable <- otutable
  dimm      <- if (!is.null(dimm_max)) min(dimm_max, nrow(otutable)) else nrow(otutable)

  for (scenario in scenarios) {
    nl <- "inv"   # S1: inv only
    k  <- if (scenario=="S3") k_s3 else NULL

    for (rate in rates) {
      loop_count <- loop_count + 1L
      message(sprintf("[%d/%d] N=%d | %s | rate=%.1f",
                      loop_count, total_loops, nsamp, scenario, rate))
      set.seed(1)

      effectsize <- mc_apply(seq_len(rep_n), function(h) {
        mc_apply(seq_len(dimm), function(x) {
          indx <- dimm*(h-1)+x
          f(indx)

          # whole = wrep nocal (caliper=NULL covers full N)
          b2_whole <- tryCatch(
            EffectSizeFun(omics, trsftable, datatable, indicator,
                          rate, b, x, NULL, scenario, k),
            error=function(e) rep(NA_real_, length(b)))

          b2_wrep <- lapply(calipers, function(cl)
            tryCatch(
              EffectSizeFun(omics, trsftable, datatable, indicator,
                            rate, b, x, cl, scenario, k),
              error=function(e) rep(NA_real_, length(b))))

          b2_worep <- lapply(calipers, function(cl)
            tryCatch(
              EffectSizeFun_wo(omics, trsftable, datatable, indicator,
                               rate, b, x, cl, scenario, k),
              error=function(e) rep(NA_real_, length(b))))

          out <- cbind(
            whole  = b2_whole,
            do.call(cbind, setNames(b2_wrep,  paste0("wrep_",  cal_names))),
            do.call(cbind, setNames(b2_worep, paste0("worep_", cal_names)))
          )
          colnames(out) <- b2_colnames
          rownames(out) <- paste0("h2_", b)
          out
        })
      })

      # ── summary ───────────────────────────────────────────
      tag_rp <- paste(scenario, nsamp, rate, sep="_")

      # all taxa mean
      b2_all <- Reduce("+", effectsize[[1]]) / length(effectsize[[1]])

      # top 20 taxa
      top20_idx <- order(psdiff, decreasing=TRUE)[1:20]
      top20_idx <- top20_idx[top20_idx <= dimm]
      b2_top20  <- if (length(top20_idx) > 0)
        Reduce("+", effectsize[[1]][top20_idx]) / length(top20_idx)
      else b2_all

      b2list[[paste0(tag_rp, "_all")]]   <- b2_all
      b2list[[paste0(tag_rp, "_top20")]] <- b2_top20
      message(sprintf("  [DONE] %s", tag_rp))

    }  # rate
  }  # scenario
}  # nsamp

# ── save ──────────────────────────────────────────────────────
save(b2list, file=file.path(data_dir, "EffectSizeList.RData"))
message(sprintf("Saved EffectSizeList.RData (%d entries)", length(b2list)))
message("05_Bias_EffectSize_Calculating.R complete.")
