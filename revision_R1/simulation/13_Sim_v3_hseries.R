## ============================================================
##  13_Sim_v3_hseries.R  — effect-size series for Supplementary Figure 3
##  Copy of 08_Sim_v3_server.R (engine 02_generator_v3.R, same seeds, same
##  matching cache) with only these differences:
##   * defaults: S1, N = 1000, all rates, h2 = 0.001,0.002,0.005,0.01,
##     calipers nocal, cal0.3, cal0.05 (the three shown in SFig 3)
##   * output:   results/v3_hseries/Power  (does not touch results/v3)
##   * cache:    reuses results/v3/match_cache full-caliper cache (subset)
##  h2 = 0.01 is included as a same-data check: its columns must equal the
##  v3 power file for the same (N, rate, scenario).
##  Usage (from code/):
##    /usr/local/R-4.4/bin/Rscript 13_Sim_v3_hseries.R --cores 30
## ============================================================

options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args

get_arg <- function(name, default = NULL, num = FALSE) {
  i <- which(args == name)
  if (length(i) == 0 || length(args) <= i) return(default)
  v <- strsplit(args[i + 1], ",")[[1]]
  if (num) as.numeric(v) else v
}

mode       <- get_arg("--mode", "power")
mc_cores   <- if (TEST_MODE) 1L else as.integer(get_arg("--cores", 6, num = TRUE))
scenarios  <- get_arg("--scenario", "S1")
nsamchar   <- get_arg("--nsamp", 1000, num = TRUE)
rates      <- get_arg("--rate", c(0.1, 0.2, 0.3), num = TRUE)
conf_share <- get_arg("--conf_share", 0.1, num = TRUE)
b_levels   <- if (mode == "t1") 0 else get_arg("--b", c(0.001, 0.002, 0.005, 0.01), num = TRUE)
rep_n      <- as.integer(get_arg("--rep", if (mode == "t1") 5000 else 1000, num = TRUE))
dimm_max   <- NULL

all_calipers <- list(nocal = NULL, cal0.3 = 0.3, cal0.2 = 0.2, cal0.1 = 0.1, cal0.05 = 0.05)
cal_names    <- get_arg("--calipers", c("nocal", "cal0.3", "cal0.05"))
calipers     <- all_calipers[cal_names]

if (TEST_MODE) {
  mode <- "power"; scenarios <- c("S1", "S2", "S3"); nsamchar <- 200; rates <- 0.2
  b_levels <- c(0, 0.01); rep_n <- 3L; dimm_max <- 5L; conf_share <- 0.1
}
stopifnot(mode %in% c("t1", "power"))

message(sprintf("v3-hseries | mode=%s | N=%s | %s | rate=%s | b=%s | cs=%s | rep=%d | cores=%d",
                mode, paste(nsamchar, collapse = ","), paste(scenarios, collapse = ","),
                paste(rates, collapse = ","), paste(b_levels, collapse = ","),
                paste(conf_share, collapse = ","), rep_n, mc_cores))

# ── paths ───────────────────────────────────────────────────
base_dir  <- "/path/to/Matching2026"
code_dir  <- file.path(base_dir, "code")
data_dir  <- file.path(base_dir, "data")
res_dir   <- file.path(base_dir, "results", "v3_hseries", if (mode == "t1") "T1" else "Power")
cache_dir <- file.path(base_dir, "results", "v3", "match_cache")
if (TEST_MODE) res_dir <- file.path(base_dir, "results", "v3_hseries", "TEST")
dir.create(res_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "02_generator_v3.R"))

load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "seed.RData"))
bat.dat <- read.csv(file.path(data_dir, "batchinfo_full.csv"), header = TRUE)
omics   <- "Metagenomics"

mc_apply <- function(X, FUN) {
  if (mc_cores == 1L) lapply(X, FUN)
  else parallel::mclapply(X, FUN, mc.cores = mc_cores, mc.preschedule = TRUE)
}

# ════════════════════════════════════════════════════════════
for (nsamp in nsamchar) {
  targetnum <- which(c(1000, 500, 200) == nsamp)
  otutable  <- otulist[[targetnum]]
  indicator <- indiclist[[targetnum]]
  lbs       <- otulist[[4]][seq_len(nsamp)]
  trsftable <- trnsf(dt = otutable, omics = omics, lb.size = lbs)
  ids       <- colnames(trsftable)
  dimm      <- if (!is.null(dimm_max)) min(dimm_max, nrow(otutable)) else nrow(otutable)

  # indicator must be aligned with sample columns (formatch is built positionally)
  if (is.null(rownames(indicator)) || all(rownames(indicator) == as.character(seq_len(nrow(indicator))))) {
    message("  [WARN] indicator has no sample-ID rownames -> assuming column order of otutable")
    rownames(indicator) <- ids
  } else if (!identical(rownames(indicator), ids)) {
    message("  [WARN] indicator rows not aligned with otutable columns -> reordering by ID")
    indicator <- indicator[ids, , drop = FALSE]
  }
  if (length(seeds) < dimm * rep_n) stop("seed vector shorter than dimm * rep_n")

  batord   <- match(ids, bat.dat[, 1])
  cmbatable <- sva::ComBat(trsftable, batch = as.character(bat.dat[batord, ]$Batch),
                           par.prior = TRUE, prior.plots = FALSE)

  # ── matching cache (outcome-free; reused by every scenario/rate/rep) ──
  cache_file <- file.path(cache_dir, sprintf("match_N%d_%s.rds", nsamp,
                                             paste(cal_names, collapse = "-")))
  full_cache <- file.path(cache_dir, sprintf("match_N%d_%s.rds", nsamp, paste(names(all_calipers), collapse = "-")))
  if (!file.exists(cache_file) && file.exists(full_cache) && !TEST_MODE) {
    mcache <- lapply(readRDS(full_cache), function(m) { if (!is.null(m)) m$cal <- m$cal[cal_names]; m })
    message(sprintf("  [CACHE] loaded %s (subset: %s)", basename(full_cache), paste(cal_names, collapse = ",")))
  } else if (file.exists(cache_file) && !TEST_MODE) {
    mcache <- readRDS(cache_file)
    message(sprintf("  [CACHE] loaded %s", basename(cache_file)))
  } else {
    message(sprintf("  [CACHE] matching N=%d, %d taxa x %d calipers ...",
                    nsamp, dimm, length(calipers)))
    mcache <- parallel::mclapply(seq_len(dimm), function(x)
      tryCatch(precompute_matching(otutable, indicator, x, calipers),
               error = function(e) NULL),
      mc.cores = if (TEST_MODE) 1L else mc_cores)
    if (!TEST_MODE) saveRDS(mcache, cache_file)
  }

  for (scenario in scenarios) for (rate in rates) for (cs in conf_share) {
    cs_tag   <- if (abs(cs - 0.1) < 1e-9) "" else sprintf("_cs%.2f", cs)
    tag      <- sprintf("NSAMP%d_REP%d_RATE%.1f_%s%s", nsamp, rep_n, rate, scenario, cs_tag)
    out_file <- file.path(res_dir, paste0(tag, ".RData"))
    if (file.exists(out_file) && !TEST_MODE) { message("  [SKIP] ", tag); next }
    message(sprintf("[RUN] %s", tag)); t0 <- Sys.time()

    res <- mc_apply(seq_len(rep_n), function(h) {
      if (h %% 50 == 1) message(sprintf("  [%s] rep %d/%d", tag, h, rep_n))
      lapply(seq_len(dimm), function(x) {
        if (is.null(mcache[[x]])) return(NULL)
        f(dimm * (h - 1) + x)
        tryCatch(analyse_one(x, mcache[[x]], trsftable, indicator,
                             cmbt_row = as.numeric(cmbatable[x, ids]),
                             scenario = scenario, rate = rate, b_levels = b_levels,
                             conf_share = cs, cal_names = cal_names),
                 error = function(e) NULL)
      })
    })

    # ── combine: list over b of matrix, row = (h-1)*dimm + x ──
    res <- lapply(res, function(h_res) if (is.list(h_res)) h_res else vector("list", dimm))
    template <- NULL
    for (h_res in res) { for (x_res in h_res) if (!is.null(x_res)) { template <- names(x_res[[1]]); break }
                         if (!is.null(template)) break }
    if (is.null(template)) { message("  [FAIL] all NULL: ", tag); next }

    pvalues <- lapply(seq_along(b_levels), function(bi) {
      mat <- do.call(rbind, lapply(res, function(h_res)
        do.call(rbind, lapply(seq_len(dimm), function(x) {
          v <- if (is.null(h_res[[x]])) NULL else h_res[[x]][[bi]]
          if (is.null(v) || length(v) != length(template)) rep(NA_real_, length(template))
          else as.numeric(v[template])
        }))))
      colnames(mat) <- template
      mat
    })
    names(pvalues) <- paste0("b", b_levels)
    meta <- list(mode = mode, nsamp = nsamp, scenario = scenario, rate = rate,
                 conf_share = cs, b_levels = b_levels, rep_n = rep_n, dimm = dimm,
                 calipers = cal_names, engine = "02_generator_v3.R", date = Sys.time())
    save(pvalues, meta, file = out_file)
    message(sprintf("  [SAVED] %s | %.1f min | n_fail=%d", tag,
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    sum(is.na(pvalues[[1]][, "prev"]))))
  }
}
message("13_Sim_v3_hseries.R complete.")
