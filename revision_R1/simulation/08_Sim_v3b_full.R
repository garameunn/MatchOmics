## ============================================================
##  08_Sim_v3b_full.R — full-sample comparators only (R2-1 add-on)
##  Same data as 08_Sim_v3_server.R: same seeds, same matching cache
##  (PS, exposure), same x1/eps draws -> same outcome per (rep, taxon).
##  Matched GEE and ComBat are NOT fitted here.
##
##  New methods (doubly robust weighted regression, "AIPW-type"):
##    IPW_DR : IPW weights + outcome model case ~ x1 + ns(logit PS, 4) + marker
##    OW_DR  : overlap weights + same outcome model
##  Check columns (must equal v3 output): LR_PS, LR_PSspline
##
##  Usage (on server, from code/):
##    Rscript 08_Sim_v3b_full.R --test
##    Rscript 08_Sim_v3b_full.R --mode power --scenario S1,S2 --b 0.01 --cores 10
##    Rscript 08_Sim_v3b_full.R --mode t1    --scenario S1,S2 --rep 2000 --cores 10
##  Output: ../results/v3b/{T1,Power}/<same tag as v3>.RData
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

mode      <- get_arg("--mode", "power")
mc_cores  <- if (TEST_MODE) 1L else as.integer(get_arg("--cores", 6, num = TRUE))
scenarios <- get_arg("--scenario", c("S1", "S2", "S3"))
nsamchar  <- get_arg("--nsamp", c(1000, 500, 200), num = TRUE)
rates     <- get_arg("--rate", c(0.1, 0.2, 0.3), num = TRUE)
b_levels  <- if (mode == "t1") 0 else get_arg("--b", 0.01, num = TRUE)
rep_n     <- as.integer(get_arg("--rep", if (mode == "t1") 2000 else 1000, num = TRUE))
dimm_max  <- NULL
cache_tag <- "nocal-cal0.3-cal0.2-cal0.1-cal0.05"     # cache written by 08_Sim_v3_server.R

if (TEST_MODE) {
  mode <- "power"; scenarios <- "S1"; nsamchar <- 200; rates <- 0.2
  b_levels <- c(0, 0.01); rep_n <- 3L; dimm_max <- 5L
}
stopifnot(mode %in% c("t1", "power"))

message(sprintf("v3b | mode=%s | N=%s | %s | rate=%s | b=%s | rep=%d | cores=%d",
                mode, paste(nsamchar, collapse = ","), paste(scenarios, collapse = ","),
                paste(rates, collapse = ","), paste(b_levels, collapse = ","), rep_n, mc_cores))

# ── paths ───────────────────────────────────────────────────
base_dir  <- "/path/to/Matching2026"
code_dir  <- file.path(base_dir, "code")
data_dir  <- file.path(base_dir, "data")
res_dir   <- file.path(base_dir, "results", "v3b", if (mode == "t1") "T1" else "Power")
cache_dir <- file.path(base_dir, "results", "v3", "match_cache")
if (TEST_MODE) res_dir <- file.path(base_dir, "results", "v3b", "TEST")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(code_dir, "00_library_setup.R"))
source(file.path(code_dir, "02_generator_v3.R"))     # f, trnsf, make_confounder, get_b2/b3, ...

load(file.path(data_dir, "otulist.RData"))
load(file.path(data_dir, "indiclist.RData"))
load(file.path(data_dir, "seed.RData"))
omics <- "Metagenomics"

mc_apply <- function(X, FUN) {
  if (mc_cores == 1L) lapply(X, FUN)
  else parallel::mclapply(X, FUN, mc.cores = mc_cores, mc.preschedule = TRUE)
}

# weighted logistic regression on a design matrix, HC0 sandwich SE
# (PS and weights treated as known, as for IPW/OW in 02_generator_v3.R)
wfit_marker <- function(case, X, w) {
  fit <- tryCatch(suppressWarnings(glm.fit(X, case, weights = w, family = quasibinomial())),
                  error = function(e) NULL)
  if (is.null(fit) || !fit$converged || anyNA(fit$coefficients)) return(NA3)
  mu <- fit$fitted.values; r <- case - mu
  bread <- tryCatch(solve(crossprod(X * (w * mu * (1 - mu)), X)), error = function(e) NULL)
  if (is.null(bread)) return(NA3)
  V <- bread %*% crossprod(X * (w * r)) %*% bread
  j <- which(colnames(X) == "marker"); b <- fit$coefficients[j]; se <- sqrt(V[j, j])
  c(2 * pnorm(-abs(b / se)), b, se)
}

# same data-generating steps (and RNG order) as analyse_one() in 02_generator_v3.R
analyse_one_full <- function(x, mcache, trsftable, indicator, scenario, rate, b_levels) {
  ids    <- colnames(trsftable)
  marker <- setNames(as.numeric(trsftable[x, ]), ids)
  ps     <- mcache$ps[ids]
  conf   <- make_confounder(scenario, mcache$ps, indicator, trsftable, x)
  b3     <- get_b3(conf, 0.1)

  n   <- length(ids)
  x1  <- rnorm(n, 1, 1)
  eps <- rnorm(n, 0, 1)

  trt   <- mcache$treat[ids]
  w_ipw <- ifelse(trt == 1, 1 / ps, 1 / (1 - ps))
  w_ow  <- ifelse(trt == 1, 1 - ps, ps)
  lps   <- qlogis(pmin(pmax(ps, 1e-6), 1 - 1e-6))
  spl   <- splines::ns(lps, df = 4)
  Xdr   <- cbind("(Intercept)" = 1, x1 = x1, spl, marker = marker)
  colnames(Xdr)[3:6] <- paste0("spl", 1:4)

  lapply(b_levels, function(B) {
    b2   <- get_b2(B, marker, conf, b3)
    case <- generate_case(marker, conf, x1, eps, b2, b3, rate)
    full <- data.frame(case = case, x1 = x1, PS = ps, marker = marker, row.names = ids)

    all_res <- list(
      LR_PS       = glm_marker(case ~ x1 + PS + marker, full),
      LR_PSspline = glm_marker(case ~ x1 + spl + marker, full),
      IPW_DR      = wfit_marker(case, Xdr, w_ipw),
      OW_DR       = wfit_marker(case, Xdr, w_ow))

    c(setNames(sapply(all_res, `[`, 1), paste0("sig_",  names(all_res))),
      setNames(sapply(all_res, `[`, 2), paste0("coef_", names(all_res))),
      setNames(sapply(all_res, `[`, 3), paste0("se_",   names(all_res))),
      max_w_ipw = max(w_ipw), prev = mean(case))
  })
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

  # indicator alignment: identical to 08_Sim_v3_server.R
  if (is.null(rownames(indicator)) || all(rownames(indicator) == as.character(seq_len(nrow(indicator))))) {
    rownames(indicator) <- ids
  } else if (!identical(rownames(indicator), ids)) {
    indicator <- indicator[ids, , drop = FALSE]
  }
  if (length(seeds) < dimm * rep_n) stop("seed vector shorter than dimm * rep_n")

  cache_file <- file.path(cache_dir, sprintf("match_N%d_%s.rds", nsamp, cache_tag))
  if (!file.exists(cache_file)) stop("matching cache not found: ", cache_file)
  mcache <- readRDS(cache_file)
  message(sprintf("  [CACHE] loaded %s", basename(cache_file)))

  for (scenario in scenarios) for (rate in rates) {
    tag      <- sprintf("NSAMP%d_REP%d_RATE%.1f_%s", nsamp, rep_n, rate, scenario)
    out_file <- file.path(res_dir, paste0(tag, ".RData"))
    if (file.exists(out_file) && !TEST_MODE) { message("  [SKIP] ", tag); next }
    message(sprintf("[RUN] %s", tag)); t0 <- Sys.time()

    res <- mc_apply(seq_len(rep_n), function(h) {
      if (h %% 100 == 1) message(sprintf("  [%s] rep %d/%d", tag, h, rep_n))
      lapply(seq_len(dimm), function(x) {
        if (is.null(mcache[[x]])) return(NULL)
        f(dimm * (h - 1) + x)
        tryCatch(analyse_one_full(x, mcache[[x]], trsftable, indicator,
                                  scenario = scenario, rate = rate, b_levels = b_levels),
                 error = function(e) NULL)
      })
    })

    # ── combine: list over b of matrix, row = (h-1)*dimm + x (same layout as v3) ──
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
                 conf_share = 0.1, b_levels = b_levels, rep_n = rep_n, dimm = dimm,
                 engine = "08_Sim_v3b_full.R + 02_generator_v3.R", date = Sys.time())
    save(pvalues, meta, file = out_file)
    message(sprintf("  [SAVED] %s | %.1f min | n_fail=%d", tag,
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    sum(is.na(pvalues[[1]][, "prev"]))))
  }
}
message("08_Sim_v3b_full.R complete.")
