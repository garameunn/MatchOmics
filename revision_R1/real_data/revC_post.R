## revC_post.R -- specification C: marker filter, FDR, concordance, tables,
##                Fig 3 inputs, comparison with the published analysis (A),
##                caliper sample sizes. Base R only.
## Input : C_final/C_all_markers.csv, C_matching_long.csv (revC_final.R)
##         realdata/protein_final.csv, metabolite_final.csv (A, published)
##         revision/rev03_matching_per_marker.csv (A SMD values)
## Usage : Rscript revC_post.R   (env REV_BASE)
base_dir <- Sys.getenv("REV_BASE"); if (!nzchar(base_dir)) stop("Set REV_BASE")
res_dir <- file.path(base_dir, "realdata"); rev_dir <- file.path(res_dir, "revision")
c_dir   <- file.path(rev_dir, Sys.getenv("REV_CDIR", "C_final"))
ALPHA <- 0.05; SMD_THRESH <- 0.2
W <- read.csv(file.path(c_dir, "C_all_markers.csv"), stringsAsFactors = FALSE, check.names = FALSE)
L <- read.csv(file.path(c_dir, "C_matching_long.csv"), stringsAsFactors = FALSE)
pls  <- c("Protein", "Metabolite")
cals <- c("nocal", "cal0.3", "cal0.2", "cal0.1", "cal0.05")
prefixes <- c("LR_unadj", "LR_PS", "LR_COM", unlist(lapply(cals, function(c) paste0(c("wrep_", "worep_"), c))))
methods <- c(LR_unadj = "LR_unadj", LR_PS = "LR_PS", LR_COM = "LR_COM",
             `GEE-EXW` = "wrep_cal0.3", `GEE-EXWO` = "worep_cal0.3")
conc4 <- c("LR_PS", "LR_COM", "wrep_cal0.3", "worep_cal0.3")
fmt_p <- function(p) ifelse(is.na(p), NA_character_, ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.1e", p)))
cell  <- function(b, p) ifelse(is.na(b) | is.na(p), NA_character_, paste0(sprintf("%.2f", b), " (", fmt_p(p), ")"))
allp  <- function(df, pre, suf) apply(df[, paste0(pre, suf), drop = FALSE], 1, function(r) all(!is.na(r) & r < ALPHA))

## ---- 1. marker filter (C weighted SMD at cal0.3, two-round matching) --------
W$smd_pass_C <- !is.na(W$`wrep_cal0.3_smd_max`) & W$`wrep_cal0.3_smd_max` < SMD_THRESH

## ---- 2. BH within platform over the C marker set; final files (Fig 3 format) --
FIN <- list()
for (pl in pls) {
  d <- W[W$platform == pl & W$smd_pass_C, ]
  out <- data.frame(marker = d$marker, stringsAsFactors = FALSE)
  for (pre in prefixes) {
    out[[paste0(pre, "_pval")]] <- d[[paste0(pre, "_pval")]]
    out[[paste0(pre, "_fdr")]]  <- p.adjust(d[[paste0(pre, "_pval")]], method = "BH")
    for (s in c("_beta", "_se", "_n")) out[[paste0(pre, s)]] <- d[[paste0(pre, s)]]
  }
  out$smd_pre_max <- d$smd_pre_max
  out$smd_max_weighted_EXW_cal0.3  <- d$`wrep_cal0.3_smd_max`
  out$smd_max_EXWO_cal0.3 <- d$`worep_cal0.3_smd_max`
  out$conc_raw <- allp(out, conc4, "_pval")
  out$conc_fdr <- allp(out, conc4, "_fdr")
  FIN[[pl]] <- out
  write.csv(out[, !names(out) %in% c("conc_raw", "conc_fdr")],
            file.path(c_dir, paste0(tolower(pl), "_final_C.csv")), row.names = FALSE)
}

## ---- counts ---------------------------------------------------------------------
cntC <- do.call(rbind, lapply(pls, function(pl) {
  d <- FIN[[pl]]
  rbind(do.call(rbind, lapply(names(methods), function(nm) data.frame(platform = pl, method = nm,
          n_markers = nrow(d), n_tested = sum(!is.na(d[[paste0(methods[nm], "_pval")]])),
          n_raw_p05 = sum(d[[paste0(methods[nm], "_pval")]] < ALPHA, na.rm = TRUE),
          n_fdr05 = sum(d[[paste0(methods[nm], "_fdr")]] < ALPHA, na.rm = TRUE), stringsAsFactors = FALSE))),
        data.frame(platform = pl, method = "Concordant (LR_PS, LR_COM, GEE-EXW, GEE-EXWO)",
          n_markers = nrow(d), n_tested = NA, n_raw_p05 = sum(d$conc_raw), n_fdr05 = sum(d$conc_fdr)))
}))
write.csv(cntC, file.path(c_dir, "C_sig_counts_raw_vs_fdr.csv"), row.names = FALSE)
print(cntC)

## concordant lists
conc_tab <- function(flag) do.call(rbind, lapply(pls, function(pl) {
  d <- FIN[[pl]][FIN[[pl]][[flag]], ]
  if (!nrow(d)) return(NULL)
  o <- data.frame(platform = pl, marker = d$marker, stringsAsFactors = FALSE)
  for (nm in names(methods)) for (s in c("_beta", "_pval", "_fdr"))
    o[[paste0(nm, s)]] <- d[[paste0(methods[nm], s)]]
  o[order(o$marker), ] }))
cr <- conc_tab("conc_raw"); cf <- conc_tab("conc_fdr")
write.csv(cr, file.path(c_dir, "C_concordant_raw.csv"), row.names = FALSE)
write.csv(cf, file.path(c_dir, "C_concordant_fdr.csv"), row.names = FALSE)

## ---- 3. Table 5 / STable 6 (C) ----------------------------------------------------
rob <- do.call(rbind, lapply(pls, function(pl) { d <- FIN[[pl]][FIN[[pl]]$conc_raw, ]; d$type <- pl; d }))
rob <- rob[order(rob$type == "Metabolite", rob$marker), ]
t5 <- data.frame(type = rob$type, marker = rob$marker, stringsAsFactors = FALSE)
for (nm in names(methods)) t5[[nm]] <- cell(rob[[paste0(methods[nm], "_beta")]], rob[[paste0(methods[nm], "_pval")]])
t5$FDR_concordant <- rob$conc_fdr
write.csv(t5, file.path(c_dir, "Table5_C.csv"), row.names = FALSE, fileEncoding = "UTF-8")
t5n <- data.frame(type = rob$type, marker = rob$marker, stringsAsFactors = FALSE)
for (nm in names(methods)) for (s in c("_beta", "_se", "_pval", "_fdr", "_n"))
  t5n[[paste0(nm, s)]] <- rob[[paste0(methods[nm], s)]]
write.csv(t5n, file.path(c_dir, "Table5_C_numeric.csv"), row.names = FALSE)

st6m <- c("wrep_nocal", "wrep_cal0.3", "wrep_cal0.05", "worep_nocal", "worep_cal0.3", "worep_cal0.05")
st6 <- data.frame(type = rob$type, marker = rob$marker, stringsAsFactors = FALSE)
for (m in st6m) st6[[m]] <- cell(rob[[paste0(m, "_beta")]], rob[[paste0(m, "_pval")]])
for (m in st6m) st6[[paste0(m, "_n")]] <- rob[[paste0(m, "_n")]]
write.csv(st6, file.path(c_dir, "STable6_C.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSTable 6 (C) cells with p >= 0.05 or NA:\n")
for (m in st6m) { p <- rob[[paste0(m, "_pval")]]; k <- which(is.na(p) | p >= ALPHA)
  if (length(k)) cat(" ", m, ":", paste(sprintf("%s(%s)", rob$marker[k], fmt_p(p[k])), collapse = ", "), "\n") }

## ---- Fig 3 inputs -------------------------------------------------------------------
venn <- do.call(rbind, lapply(pls, function(pl) {
  d <- FIN[[pl]]
  do.call(rbind, lapply(c(`LR-PS` = "LR_PS", `LR-COM` = "LR_COM", `GEE-EXW` = "wrep_cal0.3", `GEE-EXWO` = "worep_cal0.3"),
    function(pre) { k <- d$marker[!is.na(d[[paste0(pre, "_pval")]]) & d[[paste0(pre, "_pval")]] < ALPHA]
      if (length(k)) data.frame(platform = pl, method = pre, marker = k, stringsAsFactors = FALSE) else NULL }))
}))
venn$method <- c(LR_PS = "LR-PS", LR_COM = "LR-COM", wrep_cal0.3 = "GEE-EXW", worep_cal0.3 = "GEE-EXWO")[venn$method]
write.csv(venn, file.path(c_dir, "Fig3_venn_sets_C.csv"), row.names = FALSE)
reg <- do.call(rbind, lapply(pls, function(pl) {
  d <- FIN[[pl]]
  s <- sapply(c("LR_PS", "LR_COM", "wrep_cal0.3", "worep_cal0.3"), function(pre) !is.na(d[[paste0(pre, "_pval")]]) & d[[paste0(pre, "_pval")]] < ALPHA)
  pat <- apply(s, 1, function(r) paste0(ifelse(r, c("P", "C", "W", "O"), "."), collapse = ""))
  tb <- table(pat); data.frame(platform = pl, pattern_PCWO = names(tb), n = as.integer(tb), stringsAsFactors = FALSE) }))
write.csv(reg, file.path(c_dir, "Fig3_venn_region_counts_C.csv"), row.names = FALSE)
## representative marker (make_fig3.R rule: most methods significant, ties -> smallest LR_PS p)
rep_m <- do.call(rbind, lapply(pls, function(pl) {
  d <- FIN[[pl]]
  d$n_sig <- rowSums(sapply(conc4, function(pre) !is.na(d[[paste0(pre, "_pval")]]) & d[[paste0(pre, "_pval")]] < ALPHA))
  d <- d[order(-d$n_sig, d$LR_PS_pval), ]
  prev <- c(Protein = "K2C1", Metabolite = "Pyroglutamic_acid")[pl]
  data.frame(platform = pl, representative_C = d$marker[1], n_sig = d$n_sig[1],
             previous_representative = prev, previous_in_C_set = prev %in% d$marker,
             previous_concordant_raw_C = prev %in% d$marker[d$conc_raw],
             previous_concordant_fdr_C = prev %in% d$marker[d$conc_fdr], stringsAsFactors = FALSE) }))
write.csv(rep_m, file.path(c_dir, "Fig3_representative_markers_C.csv"), row.names = FALSE)
print(rep_m)

## ---- 5. comparison with the published analysis (A) -----------------------------------
A <- list(); pp <- function(b, s) 2 * pnorm(-abs(b / s))
smdA <- read.csv(file.path(rev_dir, "rev03_matching_per_marker.csv"), stringsAsFactors = FALSE)
smdA <- smdA[smdA$caliper == "cal0.3" & smdA$type == "wrep", c("platform", "marker", "smd_max_w_inv")]
for (pl in pls) {
  a <- read.csv(file.path(res_dir, paste0(tolower(pl), "_final.csv")), stringsAsFactors = FALSE, check.names = FALSE)
  for (pre in c("wrep_cal0.3", "worep_cal0.3")) a[[paste0(pre, "_pval")]] <- pp(a[[paste0(pre, "_beta")]], a[[paste0(pre, "_se")]])
  for (pre in names(methods)) a[[paste0(methods[pre], "_fdr")]] <- p.adjust(a[[paste0(methods[pre], "_pval")]], "BH")
  a$conc_raw <- allp(a, conc4, "_pval"); a$conc_fdr <- allp(a, conc4, "_fdr")
  A[[pl]] <- a
}
cmp <- list()
for (pl in pls) {
  a <- A[[pl]]; cc <- FIN[[pl]]
  add <- function(metric, va, vc) cmp[[length(cmp) + 1]] <<- data.frame(platform = pl, metric = metric,
                                     A_published = va, C_final = vc, change = vc - va, stringsAsFactors = FALSE)
  add("markers tested", sum(W$platform == pl), sum(W$platform == pl))
  add("markers passing SMD filter", nrow(a), nrow(cc))
  add("  of which in both sets", length(intersect(a$marker, cc$marker)), length(intersect(a$marker, cc$marker)))
  for (nm in names(methods)) {
    add(paste(nm, "raw p<0.05"), sum(a[[paste0(methods[nm], "_pval")]] < ALPHA, na.rm = TRUE),
        sum(cc[[paste0(methods[nm], "_pval")]] < ALPHA, na.rm = TRUE))
    add(paste(nm, "FDR<0.05"), sum(a[[paste0(methods[nm], "_fdr")]] < ALPHA, na.rm = TRUE),
        sum(cc[[paste0(methods[nm], "_fdr")]] < ALPHA, na.rm = TRUE))
  }
  add("GEE-EXW not converged", sum(is.na(a$`wrep_cal0.3_pval`)), sum(is.na(cc$`wrep_cal0.3_pval`)))
  add("concordant (raw p)", sum(a$conc_raw), sum(cc$conc_raw))
  add("concordant (FDR)", sum(a$conc_fdr), sum(cc$conc_fdr))
}
cmp <- do.call(rbind, cmp)
write.csv(cmp, file.path(c_dir, "Compare_A_vs_C_counts.csv"), row.names = FALSE)
print(cmp)

## marker-level changes with reasons
mk_chg <- list()
lab <- c(LR_PS = "LR_PS", LR_COM = "LR_COM", wrep_cal0.3 = "GEE-EXW", worep_cal0.3 = "GEE-EXWO")
for (pl in pls) for (crit in c("raw", "fdr")) {
  a <- A[[pl]]; cc <- FIN[[pl]]; w <- W[W$platform == pl, ]
  flag <- paste0("conc_", crit); suf <- if (crit == "raw") "_pval" else "_fdr"
  inA <- a$marker[a[[flag]]]; inC <- cc$marker[cc[[flag]]]
  for (mk in union(inA, inC)) {
    ra <- a[a$marker == mk, ]; rc <- cc[cc$marker == mk, ]; rw <- w[w$marker == mk, ]
    va <- if (nrow(ra)) sapply(conc4, function(p) ra[[paste0(p, suf)]]) else rep(NA, 4)
    vc <- if (nrow(rc)) sapply(conc4, function(p) rc[[paste0(p, suf)]]) else rep(NA, 4)
    sA <- smdA$smd_max_w_inv[smdA$platform == pl & smdA$marker == mk]
    status <- if (mk %in% inA && mk %in% inC) "kept" else if (mk %in% inA) "dropped" else "added"
    reason <- ""
    if (status == "dropped") {
      if (!nrow(rc)) reason <- sprintf("C SMD filter: max weighted SMD %.3f >= 0.2", rw$`wrep_cal0.3_smd_max`)
      else { bad <- which(is.na(vc) | vc >= ALPHA)
        reason <- paste(sprintf("%s %s %s -> %s", lab[conc4[bad]], if (crit == "raw") "p" else "FDR",
                                fmt_p(va[bad]), fmt_p(vc[bad])), collapse = "; ") }
    } else if (status == "added") {
      if (!nrow(ra)) reason <- sprintf("newly passes SMD filter (A max SMD %.3f, C %.3f)", sA, rw$`wrep_cal0.3_smd_max`)
      else { bad <- which(is.na(va) | va >= ALPHA)
        reason <- paste(sprintf("%s %s %s -> %s", lab[conc4[bad]], if (crit == "raw") "p" else "FDR",
                                fmt_p(va[bad]), fmt_p(vc[bad])), collapse = "; ") }
    }
    r <- data.frame(platform = pl, criterion = crit, marker = mk, status = status, reason = reason,
                    in_A_set = nrow(ra) > 0, in_C_set = nrow(rc) > 0,
                    smd_max_A = round(sA, 3), smd_max_C = round(rw$`wrep_cal0.3_smd_max`, 3), stringsAsFactors = FALSE)
    for (j in seq_along(conc4)) { r[[paste0("A_", lab[conc4[j]])]] <- signif(va[j], 3); r[[paste0("C_", lab[conc4[j]])]] <- signif(vc[j], 3) }
    mk_chg[[length(mk_chg) + 1]] <- r
  }
}
mk_chg <- do.call(rbind, mk_chg)
mk_chg <- mk_chg[order(mk_chg$criterion != "raw", mk_chg$platform != "Protein", mk_chg$status, mk_chg$marker), ]
write.csv(mk_chg, file.path(c_dir, "Compare_A_vs_C_concordant_markers.csv"), row.names = FALSE)
cat("\nconcordant set changes:\n"); print(mk_chg[mk_chg$status != "kept", c("platform", "criterion", "marker", "status", "reason")])

## SMD filter changes overall
fch <- do.call(rbind, lapply(pls, function(pl) {
  a <- A[[pl]]$marker; cm <- FIN[[pl]]$marker
  data.frame(platform = pl, marker = c(setdiff(a, cm), setdiff(cm, a)),
             change = c(rep("A only (dropped by C filter)", length(setdiff(a, cm))), rep("C only (added by C filter)", length(setdiff(cm, a)))),
             stringsAsFactors = FALSE) }))
fch$smd_max_A <- round(smdA$smd_max_w_inv[match(paste(fch$platform, fch$marker), paste(smdA$platform, smdA$marker))], 3)
fch$smd_max_C <- round(W$`wrep_cal0.3_smd_max`[match(paste(fch$platform, fch$marker), paste(W$platform, W$marker))], 3)
write.csv(fch, file.path(c_dir, "Compare_A_vs_C_marker_set.csv"), row.names = FALSE)

## ---- 6. caliper sample sizes (C) ----------------------------------------------------
q  <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, .25), quantile(x, .75))
q1 <- function(x) sprintf("%.1f [%.1f-%.1f]", median(x), quantile(x, .25), quantile(x, .75))
ss <- list()
for (pl in pls) for (set in c("C_marker_set", "all_markers")) for (m in c("GEE-EXW", "GEE-EXWO")) for (cal in cals) {
  r <- L[L$platform == pl & L$method == m & L$caliper == cal, ]
  if (set == "C_marker_set") r <- r[r$marker %in% FIN[[pl]]$marker, ]
  ss[[length(ss) + 1]] <- data.frame(platform = pl, marker_set = set, n_markers_set = nrow(r), method = m, caliper = cal,
    N_total = q(r$n_total), N_treated = q(r$n_trt_total), n_matched = q(r$n_matched), neff = q(r$neff),
    pct_treated_excluded = q1(r$pct_trt_excluded), pct_treated_rescued_r2 = q1(100 * r$n_trt_rescued_r2 / r$n_trt_total),
    max_control_reuse = q(r$max_ctrl_reuse),
    n_matched_mean = round(mean(r$n_matched), 1), neff_mean = round(mean(r$neff), 1),
    pct_excl_mean = round(mean(r$pct_trt_excluded), 2), stringsAsFactors = FALSE)
}
ss <- do.call(rbind, ss)
write.csv(ss, file.path(c_dir, "C_caliper_samplesize_summary.csv"), row.names = FALSE)
print(ss[ss$marker_set == "C_marker_set", c("platform", "method", "caliper", "n_matched", "neff", "pct_treated_excluded", "max_control_reuse")])
cz <- L[paste(L$platform, L$marker) %in% paste(cr$platform, cr$marker), ]
write.csv(cz, file.path(c_dir, "C_caliper_samplesize_concordant.csv"), row.names = FALSE)
