## ============================================================
##  12_diag_S3_ex.R — why do exchangeable GEE fits return NA in S3?
##  Re-creates (rep 1, taxon x) exactly as 08_Sim_v3_server.R does
##  (same data loading, indicator alignment, matching cache, seed
##  f(dimm*(h-1)+x) with h = 1, b = 0.01) and calls geem() WITHOUT
##  tryCatch so the error message is shown.
##  Self-check: LR_PS p-value of the rebuilt data must equal row x of
##  the saved v3 power file (b = 0.01).
##  Usage (from code/):
##    /usr/local/R-4.4/bin/Rscript 12_diag_S3_ex.R [N=1000] [taxon=1]
## ============================================================
args  <- commandArgs(trailingOnly = TRUE)
nsamp <- if (length(args) >= 1) as.integer(args[1]) else 1000L
x     <- if (length(args) >= 2) as.integer(args[2]) else 1L

cat("R:", R.version.string, "| Rscript:", file.path(R.home("bin"), "Rscript"), "\n")
cat(".libPaths():", paste(.libPaths(), collapse = " ; "), "\n")
for (p in c("Matrix", "geeM", "MatchIt"))
  cat(sprintf("  %-8s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")))

base_dir <- "/path/to/Matching2026"
source(file.path(base_dir, "code", "00_library_setup.R"))
source(file.path(base_dir, "code", "02_generator_v3.R"))
load(file.path(base_dir, "data", "otulist.RData"))
load(file.path(base_dir, "data", "indiclist.RData"))
load(file.path(base_dir, "data", "seed.RData"))
stopifnot(is.function(f))                         # seed helper from 02_generator_v3.R

# ── same setup as 08_Sim_v3_server.R (lines 85-101) ──
targetnum <- which(c(1000, 500, 200) == nsamp)
otutable  <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
lbs       <- otulist[[4]][seq_len(nsamp)]
trsftable <- trnsf(dt = otutable, omics = "Metagenomics", lb.size = lbs)
ids       <- colnames(trsftable)
dimm      <- nrow(otutable)
if (is.null(rownames(indicator)) || all(rownames(indicator) == as.character(seq_len(nrow(indicator))))) {
  rownames(indicator) <- ids
} else if (!identical(rownames(indicator), ids)) {
  indicator <- indicator[ids, , drop = FALSE]
}
mcache <- readRDS(file.path(base_dir, "results", "v3", "match_cache",
                            sprintf("match_N%d_nocal-cal0.3-cal0.2-cal0.1-cal0.05.rds", nsamp)))[[x]]

# same file choice as pick_file() in 10_summarize_v3.R: largest replicate count
saved_row <- function(dir, sc) {
  fl <- list.files(file.path(base_dir, "results", "v3", dir),
                   sprintf("^NSAMP%d_REP[0-9]+_RATE0.2_%s\\.RData$", nsamp, sc), full.names = TRUE)
  if (length(fl) == 0) return(NULL)
  reps <- as.integer(sub(".*_REP([0-9]+)_.*", "\\1", basename(fl)))
  fl <- fl[which.max(reps)]
  cat("  saved file:", basename(fl), "\n")
  e <- new.env(); load(fl, envir = e)
  e$pvalues[[length(e$pvalues)]][x, ]
}

for (sc in c("S1", "S3")) {
  # ── same steps and RNG order as analyse_one() (02_generator_v3.R) ──
  f(dimm * (1 - 1) + x)
  marker <- setNames(as.numeric(trsftable[x, ]), ids)
  ps     <- mcache$ps[ids]
  conf   <- make_confounder(sc, mcache$ps, indicator, trsftable, x)
  b3     <- get_b3(conf, 0.1)
  x1     <- rnorm(length(ids), 1, 1)
  eps    <- rnorm(length(ids), 0, 1)
  b2     <- get_b2(0.01, marker, conf, b3)
  case   <- generate_case(marker, conf, x1, eps, b2, b3, 0.2)
  full   <- data.frame(case = case, x1 = x1, PS = ps, marker = marker, row.names = ids)

  cat(sprintf("\n===== %s (rep 1, taxon %d, b = 0.01) =====\n", sc, x))
  sv <- saved_row("Power", sc)
  p_rebuilt <- glm_marker(case ~ x1 + PS + marker, full)[1]
  if (is.null(sv)) {
    cat(sprintf("self-check LR_PS p: rebuilt %.6g | saved file not found (skipped)\n", p_rebuilt))
  } else {
    cat(sprintf("self-check LR_PS p: rebuilt %.6g | saved %.6g | %s\n", p_rebuilt, sv["sig_LR_PS"],
                if (isTRUE(all.equal(unname(p_rebuilt), unname(sv["sig_LR_PS"])))) "MATCH" else "DIFFERENT"))
    cat(sprintf("saved wrep_cal0.3: ind p = %.4g, ex p = %s\n", sv["sig_wrep_cal0.3_ind"],
                format(sv["sig_wrep_cal0.3_ex"])))
  }

  tb <- matched_table(full, mcache$cal$cal0.3$wrep)
  cat(sprintf("matched n = %d, clusters = %d, prev = %.3f\n", nrow(tb), nlevels(tb$pair_taxa), mean(tb$case)))
  for (cs in c("independence", "exchangeable")) {
    r <- try(geeM::geem(case ~ x1 + PS + marker, data = tb, family = "binomial",
                        corstr = cs, id = tb$pair_taxa, weights = tb$mm,
                        sandwich = TRUE, useP = TRUE, maxit = 25, tol = 1e-3), silent = TRUE)
    if (inherits(r, "try-error")) { cat(sprintf("%s: ERROR %s", cs, r)); next }
    cat(sprintf("%s: converged=%s beta_marker=%.4f p=%.4g alpha=%s\n", cs, r$converged,
                r$beta[4], summary(r)$p[4], if (is.null(r$alpha)) "NULL" else round(r$alpha, 4)))
  }
}
