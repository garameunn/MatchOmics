## ============================================================
##  16_make_figures_tables_v3.R — revised figures and tables (run locally)
##  Run from the project root with R >= 4.3 (ggplot2, dplyr, tidyr, openxlsx):
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\16_make_figures_tables_v3.R
##
##  Inputs
##    결과파일_20260517/v3_summary/summary_v3_long.csv, summary_v3b_long.csv  (10/11_summarize, 10/26)
##    결과파일_20260517/v3_summary/t1_alpha_v3.csv, hseries_adjpw_v3.csv,
##                                matched_size_v3.csv                     (15_extract_fig_data_v3.R)
##    결과파일_20260517/realdata/revision/C_final/*                          (real data, option C)
##  Outputs (new folders; nothing existing is overwritten)
##    manuscript/Figures_rev/  Fig1.png Fig2.png Fig3.png SFig2.png SFig3.png
##    NOTE (2026-10-29): the final Figure 2 (Scenarios 1 and 2) and Suppl. Figure 4 come from
##    20b_fig2_matchomics_pair.R; re-run it after this script, which writes the earlier S1-only Fig2.png.
##    manuscript/Tables_rev/   Tables_main_rev.xlsx (Tables 2-5), STables_rev.xlsx
##  Each table sheet starts with a title row and a legend row (Reviewer 2, Comment 10).
##  Figures/tables whose input file is missing are skipped with a message.
## ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(openxlsx) })

CS <- "ind"   # working correlation of the GEE curves in the figures ("ind" or "ex");
              # tables always list both (INDW/EXW rows). ind and ex differ by <= 0.004.

sum_dir <- "결과파일_20260517/v3_summary"
rd_dir  <- "결과파일_20260517/realdata/revision/C_final"
fig_out <- "manuscript/Figures_rev"; tab_out <- "manuscript/Tables_rev"
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_out, showWarnings = FALSE, recursive = TRUE)

opt_csv <- function(dir, f) {
  p <- file.path(dir, f)
  if (file.exists(p)) read.csv(p, check.names = FALSE) else { message("  [skip] not found: ", p); NULL }
}
long <- opt_csv(sum_dir, "summary_v3_long.csv")
lb   <- opt_csv(sum_dir, "summary_v3b_long.csv")
allm <- bind_rows(long, if (!is.null(lb)) filter(lb, method %in% c("IPW_DR", "OW_DR")))
t1a   <- opt_csv(sum_dir, "t1_alpha_v3.csv")
hser  <- opt_csv(sum_dir, "hseries_adjpw_v3.csv")
msize <- opt_csv(sum_dir, "matched_size_v3.csv")

## ── labels ────────────────────────────────────────────────────
cal_lab <- c(nocal = "None", cal0.3 = "0.3", cal0.2 = "0.2", cal0.1 = "0.1", cal0.05 = "0.05")
fixed_lab <- c(LR_unadj = "LR", LR_PS = "LR-PS", LR_Combat = "LR-COM", LR_PSspline = "LR-PSspline",
               IPW = "IPW", OW = "OW", IPW_DR = "IPW-REG", OW_DR = "OW-REG")
lab_method <- function(m) {
  out <- unname(fixed_lab[m]); g <- is.na(out) & grepl("^(wrep|worep)_", m)
  out[g] <- vapply(strsplit(m[g], "_"), function(p)
    sprintf("GEE-%s%s (Cal=%s)", ifelse(p[3] == "ind", "IND", "EX"),
            ifelse(p[1] == "wrep", "W", "WO"), cal_lab[[p[2]]]), "")
  out
}
row_order <- c(names(fixed_lab),
               unlist(lapply(c("wrep", "worep"), function(t)
                 unlist(lapply(names(cal_lab), function(cl) paste0(t, "_", cl, c("_ind", "_ex")))))))

## ── table helpers ─────────────────────────────────────────────
wide_tab <- function(scen, val, ss = "all", Ns = c(1000, 500, 200)) {
  d <- allm %>% filter(scenario == scen, subset == ss, method %in% row_order, N %in% Ns) %>%
    mutate(col = sprintf("N = %d|%.1f", N, rate))
  w <- d %>% select(method, col, value = all_of(val)) %>% pivot_wider(names_from = col, values_from = value)
  cols <- as.vector(outer(sprintf("%.1f", c(0.1, 0.2, 0.3)), sprintf("N = %d", Ns), function(r, n) paste0(n, "|", r)))
  w <- w[match(intersect(row_order, w$method), w$method), c("method", intersect(cols, names(w)))]
  w$method <- lab_method(w$method); names(w)[1] <- "Method|Prevalence"; as.data.frame(w)
}
bold <- createStyle(textDecoration = "bold")
add_sheet <- function(wb, name, title, legend, df, digits = 3) {
  addWorksheet(wb, name)
  writeData(wb, name, title, startRow = 1); addStyle(wb, name, bold, rows = 1, cols = 1)
  writeData(wb, name, legend, startRow = 2)
  df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)
  if (any(grepl("\\|", names(df)))) {                        # two header rows: block | sub-column
    top <- sub("\\|.*", "", names(df)); bot <- sub(".*\\|", "", names(df))
    top[duplicated(top)] <- ""
    writeData(wb, name, t(top), startRow = 4, colNames = FALSE)
    writeData(wb, name, t(bot), startRow = 5, colNames = FALSE)
    addStyle(wb, name, bold, rows = 4:5, cols = seq_along(top), gridExpand = TRUE)
    writeData(wb, name, df, startRow = 6, colNames = FALSE)
  } else {
    writeData(wb, name, df, startRow = 4, headerStyle = bold)
  }
  setColWidths(wb, name, cols = 1, widths = 26)
}
abbr <- paste("LR, unadjusted logistic regression; LR-PS, logistic regression adjusted for the propensity score (PS);",
              "LR-COM, logistic regression on ComBat-corrected data; LR-PSspline, LR-PS with a natural cubic spline (4 df) of the logit PS;",
              "IPW/OW, inverse-probability/overlap weighting; IPW-REG/OW-REG, weighting combined with regression on the PS spline;",
              "GEE-INDW/EXW, with-replacement (two-round) matching analysed by weighted GEE with independence/exchangeable working correlation;",
              "GEE-INDWO/EXWO, the same for without-replacement matching; Cal, caliper in SD units of the PS.")

## ── main tables ───────────────────────────────────────────────
wb <- createWorkbook()
if (nrow(allm)) {
  add_sheet(wb, "Table 2", "Table 2. Empirical type-1 error under Scenario 1 at the 5% significance level.",
            paste("Proportion of p-values < 0.05 over 2,000 replicates, each analysing every taxon as the marker.", abbr), wide_tab("S1", "T1"))
  add_sheet(wb, "Table 3", "Table 3. Adjusted empirical power under Scenario 1 with marker effect size h2 = 0.01.",
            paste("Proportion of 1,000 replicates (all taxa) with p below the method-specific empirical type-1 error threshold at the 5% level.", abbr),
            wide_tab("S1", "AdjPw"))
  t4 <- allm %>% filter(N == 1000, rate == 0.2, subset == "all", scenario %in% c("S1", "S2"), method %in% row_order) %>%
    select(scenario, method, T1, AdjPw) %>% pivot_wider(names_from = scenario, values_from = c(T1, AdjPw)) %>%
    select(method, any_of(c("T1_S1", "AdjPw_S1", "T1_S2", "AdjPw_S2")))
  t4 <- t4[match(intersect(row_order, t4$method), t4$method), ]; t4$method <- lab_method(t4$method)
  names(t4) <- c("Method|", "Scenario 1 (PS-specified)|Type-1 error", "Scenario 1 (PS-specified)|Adjusted power",
                 "Scenario 2 (PS-misspecified)|Type-1 error", "Scenario 2 (PS-misspecified)|Adjusted power")[seq_along(t4)]
  add_sheet(wb, "Table 4", "Table 4. Type-1 error and adjusted power under PS-specified (Scenario 1) and PS-misspecified (Scenario 2) conditions.",
            paste("N = 1,000, prevalence 0.2, h2 = 0.01.", abbr), as.data.frame(t4))
}
t5 <- opt_csv(rd_dir, "Table5_C.csv")
if (!is.null(t5)) add_sheet(wb, "Table 5",
  "Table 5. Protein and metabolite markers concordantly associated with residence in a pollutant-exposed area in the FROM cohort.",
  paste("Markers significant at raw p < 0.05 under all four adjusted methods (LR-PS, LR-COM, GEE-EXW, GEE-EXWO; 0.3 SD caliper) among 191 proteins and 34 metabolites",
        "passing the weighted post-matching balance criterion (maximum SMD < 0.2). Cells: coefficient (exact p). FDR_concordant: significant at FDR < 0.05",
        "(Benjamini-Hochberg within platform) under all four methods. LR-COM coefficients are on the ComBat-corrected scale."), t5)
saveWorkbook(wb, file.path(tab_out, "Tables_main_rev.xlsx"), overwrite = TRUE)
message("  saved Tables_main_rev.xlsx (Table 1 is the literature table and is unchanged)")

## ── supplementary tables ──────────────────────────────────────
wb <- createWorkbook()
if (!is.null(msize)) {
  ms <- msize %>% mutate(row = paste0(ifelse(type == "wrep", "With replacement", "Without replacement"), " (Cal=", cal_lab[caliper], ")"),
                         col_n = sprintf("N = %d|Matched n", N), col_e = sprintf("N = %d|Effective n", N))
  st1 <- bind_rows(ms %>% transmute(row, col = col_n, v = matched_n), ms %>% transmute(row, col = col_e, v = neff)) %>%
    pivot_wider(names_from = col, values_from = v)
  st1 <- st1[, c("row", as.vector(t(outer(sprintf("N = %d", c(1000, 500, 200)), c("Matched n", "Effective n"), paste, sep = "|"))))]
  names(st1)[1] <- "Matching|"
  add_sheet(wb, "STable 1", "Supplementary Table 1. Mean matched and effective sample size by caliper (simulation).",
            "Means over taxa (outcome-independent). Effective sample size by the Kish approximation, (sum w)^2 / sum w^2, with w = m / mean(m).",
            as.data.frame(st1), digits = 1)
}
if (nrow(allm)) {
  add_sheet(wb, "STable 2", "Supplementary Table 2. Empirical type-1 error under Scenario 2 (PS misspecification) at the 5% level.", paste("Layout as in Table 2.", abbr), wide_tab("S2", "T1"))
  add_sheet(wb, "STable 3", "Supplementary Table 3. Empirical type-1 error under Scenario 3 (near-unmeasured confounding) at the 5% level.", paste("Layout as in Table 2.", abbr), wide_tab("S3", "T1"))
  add_sheet(wb, "STable 4", "Supplementary Table 4. Adjusted empirical power under Scenario 2 with h2 = 0.01.", paste("Layout as in Table 3.", abbr), wide_tab("S2", "AdjPw"))
  add_sheet(wb, "STable 5", "Supplementary Table 5. Adjusted empirical power under Scenario 3 with h2 = 0.01.", paste("Layout as in Table 3.", abbr), wide_tab("S3", "AdjPw"))
}
st6 <- opt_csv(rd_dir, "STable6_C.csv")
if (!is.null(st6)) add_sheet(wb, "STable 6", "Supplementary Table 6. Caliper sensitivity for the concordant protein and metabolite markers in the FROM cohort.",
  "Coefficient (exact p) of the 22 concordant markers (Table 5) under with- and without-replacement matching without a caliper and with 0.3 and 0.05 SD calipers; *_n, matched sample size.", st6)
if (nrow(allm)) {
  add_sheet(wb, "New - Top20 T1", "Supplementary Table [No.]. Type-1 error under Scenario 1 for the 20 most confounded taxa.",
            paste("Taxa with the largest PS difference between matching groups. Layout as in Table 2.", abbr), wide_tab("S1", "T1", "top20"))
  add_sheet(wb, "New - Top20 AdjPw", "Supplementary Table [No.]. Adjusted power under Scenario 1 for the 20 most confounded taxa.",
            paste("Layout as in Table 3.", abbr), wide_tab("S1", "AdjPw", "top20"))
  ie <- long %>% filter(grepl("_(ind|ex)$", method)) %>% mutate(base = sub("_(ind|ex)$", "", method), cs = sub(".*_", "", method)) %>%
    select(N, scenario, rate, subset, base, cs, AdjPw, T1) %>% pivot_wider(names_from = cs, values_from = c(AdjPw, T1)) %>%
    group_by(Scenario = scenario) %>% summarise(`Max |AdjPw ind - ex|` = max(abs(AdjPw_ind - AdjPw_ex), na.rm = TRUE),
                                                `Max |T1 ind - ex|` = max(abs(T1_ind - T1_ex), na.rm = TRUE), .groups = "drop")
  rho <- data.frame(`Confounder share of liability variance` = c("10% (main setting)", "20%", "40%"),
                    `Mean estimated rho` = c(0.0003, 0.0014, 0.0036),
                    `Rejection rate, independence` = c(0.543, 0.571, 0.642),
                    `Rejection rate, exchangeable` = c(0.542, 0.572, 0.644),
                    `Mean cluster size` = c(2.21, 2.21, 2.21), check.names = FALSE)
  addWorksheet(wb, "New - Working corr")
  writeData(wb, "New - Working corr", "Supplementary Table [No.]. Working correlation: estimated within-cluster correlation and independence vs exchangeable results.", startRow = 1)
  addStyle(wb, "New - Working corr", bold, rows = 1, cols = 1)
  writeData(wb, "New - Working corr", paste("Upper: diagnostic run under Scenario 1 (N = 1,000, prevalence 0.2, h2 = 0.01, 50 replicates x all taxa), with-replacement GEE with a 0.3 caliper;",
            "rejection rates at the nominal 0.05 level. Lower: maximum difference between the two working correlations over all main simulation conditions."), startRow = 2)
  writeData(wb, "New - Working corr", rho, startRow = 4, headerStyle = bold)
  writeData(wb, "New - Working corr", as.data.frame(ie), startRow = 9, headerStyle = bold)
}
fdr <- opt_csv(rd_dir, "C_sig_counts_raw_vs_fdr.csv")
if (!is.null(fdr)) add_sheet(wb, "New - raw vs FDR", "Supplementary Table [No.]. Number of significant markers by method at raw p < 0.05 and FDR < 0.05 in the FROM cohort.",
  "FDR by the Benjamini-Hochberg procedure within each platform among markers passing the balance criterion (191 proteins, 34 metabolites); 0.3 SD caliper for matching methods.", fdr)
smd <- opt_csv(rd_dir, "STable_SMDfail_C.csv")
if (!is.null(smd)) add_sheet(wb, "New - balance excl", "Supplementary Table [No.]. Previously reported markers that do not meet the weighted balance criterion.",
  paste("Markers concordant in the original submission whose maximum weighted SMD (w = m / mean(m)) after two-round matching (0.3 caliper) is >= 0.2;",
        "they are excluded from the concordance analysis because of covariate imbalance, not because of their p-values."),
  smd %>% filter(grepl("^A-concordant", category)) %>% select(-category, -in_A_SMD_set))
ord <- opt_csv(rd_dir, "STable_order_sensitivity_C.csv")
if (!is.null(ord)) add_sheet(wb, "New - matching order", "Supplementary Table [No.]. Sensitivity of GEE-EXW results to the greedy matching order.",
  paste("Primary result: default deterministic order (exposed subjects processed by decreasing PS).",
        "Compared with 20 runs in which the processing order in both matching rounds was randomised; borderline = significant in fewer than all runs."), ord)
cut <- opt_csv(rd_dir, "rev04_cutoff_concordant26_C.csv")
if (!is.null(cut)) add_sheet(wb, "New - cutoff", "Supplementary Table [No.]. Sensitivity of the concordant markers to the cutoff used to define the matched groups.",
  "Median split (primary) vs tertiles (upper vs lower third) and quartiles (upper vs lower quarter); PS re-estimated in each analysis set; marker continuous in the outcome model.", cut)
rn <- opt_csv(rd_dir, "C_caliper_samplesize_summary.csv")
if (!is.null(rn)) add_sheet(wb, "New - real data n", "Supplementary Table [No.]. Matched sample size by caliper in the FROM cohort.",
  "Median [IQR] over markers passing the balance criterion. Effective n by the Kish approximation; % treated excluded = high-abundance subjects without a match.",
  rn %>% filter(marker_set == "C_marker_set"))
saveWorkbook(wb, file.path(tab_out, "STables_rev.xlsx"), overwrite = TRUE)
message("  saved STables_rev.xlsx (new tables are labelled 'New - ...' with [No.] to be numbered by the authors)")

## ── figures ───────────────────────────────────────────────────
base_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9), legend.key.width = unit(1.8, "cm"),
        strip.background = element_rect(fill = "gray92"), strip.text = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(), axis.text = element_text(size = 8))
N_lev <- c("N = 200", "N = 500", "N = 1000")
alpha_lev <- c("0.001", "0.005", "0.010", "0.050", "0.100")
m_cs <- function(x) paste0(x, "_", CS)

## Figure 1: T1 vs alpha, S1, averaged over prevalence (as in the original Figure 1)
if (!is.null(t1a)) {
  f1 <- t1a %>% filter(scenario == "S1", subset == "all",
                       method %in% m_cs(c("wrep_nocal", "wrep_cal0.3", "wrep_cal0.05", "worep_nocal", "worep_cal0.3", "worep_cal0.05"))) %>%
    group_by(N, method, alpha) %>% summarise(T1 = mean(T1), .groups = "drop") %>%
    mutate(rep_type = factor(ifelse(grepl("^wrep", method), "With Replacement", "Without Replacement")),
           caliper = factor(case_when(grepl("nocal", method) ~ "No Caliper", grepl("cal0\\.3", method) ~ "Caliper = 0.3",
                                      TRUE ~ "Caliper = 0.05"), levels = c("Caliper = 0.05", "Caliper = 0.3", "No Caliper")),
           alpha = factor(sprintf("%.3f", alpha), levels = alpha_lev), N_lab = factor(paste0("N = ", N), levels = N_lev))
  ref <- data.frame(alpha = factor(alpha_lev, levels = alpha_lev), T1 = c(0.001, 0.005, 0.01, 0.05, 0.1))
  p1 <- ggplot(f1, aes(alpha, T1, color = N_lab, shape = N_lab, group = N_lab)) +
    geom_line(data = ref, aes(alpha, T1, group = 1), inherit.aes = FALSE, color = "gray65", linetype = "dotted", linewidth = 0.4) +
    geom_hline(yintercept = 0.05, color = "gray40", linetype = "dashed", linewidth = 0.27) +
    geom_line(linewidth = 0.67) + geom_point(size = 2.5) + facet_grid(caliper ~ rep_type) +
    scale_color_manual(values = c("N = 200" = "#56B4E9", "N = 500" = "#009E73", "N = 1000" = "#E69F00")) +
    scale_shape_manual(values = c("N = 200" = 16, "N = 500" = 17, "N = 1000" = 15)) +
    labs(x = "Significance alpha", y = "Empirical type-1 error", color = "Sample size", shape = "Sample size") +
    base_theme + theme(legend.position = "right")
  ggsave(file.path(fig_out, "Fig1.png"), p1, width = 7, height = 7, dpi = 300); message("  saved Fig1.png")
}

## Figure 2: adjusted power vs prevalence, all taxa vs top 20, by N
f2_m <- c("LR_PS", "LR_PSspline", m_cs("wrep_nocal"), m_cs("wrep_cal0.3"), m_cs("worep_cal0.3"))
f2_lab <- setNames(c("LR-PS", "LR-PSspline", "GEE-EXW (no caliper)", "GEE-EXW (cal = 0.3)", "GEE-EXWO (cal = 0.3)"), f2_m)
f2_col <- setNames(c("#222222", "#0072B2", "#E69F00", "#D32F2F", "#009E73"), f2_m)
f2_lty <- setNames(c("solid", "dotdash", "dashed", "solid", "dashed"), f2_m)
f2_lw  <- setNames(c(0.53, 0.6, 0.53, 1.07, 0.53), f2_m)
if (nrow(allm)) {
  f2 <- allm %>% filter(scenario == "S1", method %in% f2_m) %>%
    mutate(method = factor(method, levels = f2_m), N_lab = factor(paste0("N = ", N), levels = N_lev),
           subset = factor(ifelse(subset == "all", "All taxa", "Top 20"), levels = c("All taxa", "Top 20")))
  p2 <- ggplot(f2, aes(rate, AdjPw, color = method, linetype = method, linewidth = method, group = method)) +
    geom_line() + geom_point(size = 2.2, show.legend = FALSE) + facet_grid(N_lab ~ subset) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
    scale_color_manual(values = f2_col, labels = f2_lab) + scale_linetype_manual(values = f2_lty, labels = f2_lab) +
    scale_linewidth_manual(values = f2_lw, labels = f2_lab) +
    labs(x = "Prevalence rate", y = "Adjusted power", color = NULL, linetype = NULL, linewidth = NULL) +
    base_theme + guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), linewidth = guide_legend(nrow = 2))
  ggsave(file.path(fig_out, "Fig2.png"), p2, width = 7.5, height = 8, dpi = 300); message("  saved Fig2.png")
}

## Figure 3: real data (option C) — copy of the figure produced by make_fig3_C.R
if (file.exists(file.path(rd_dir, "Fig3_C.png"))) {
  file.copy(file.path(rd_dir, "Fig3_C.png"), file.path(fig_out, "Fig3.png"), overwrite = TRUE); message("  copied Fig3.png (option C)")
}

## Supplementary Figures 2 and 3: N = 1000, S1, rows = prevalence, cols = caliper
mt_col <- c("GEE-EXW" = "#D32F2F", "GEE-EXWO" = "#009E73", "LR-PS" = "gray30", "LR-COM" = "gray60", "LR-PSspline" = "#0072B2")
mt_lty <- c("GEE-EXW" = "solid", "GEE-EXWO" = "dashed", "LR-PS" = "solid", "LR-COM" = "dashed", "LR-PSspline" = "dotdash")
mt_lw  <- c("GEE-EXW" = 0.93, "GEE-EXWO" = 0.93, "LR-PS" = 0.47, "LR-COM" = 0.47, "LR-PSspline" = 0.6)
cal3 <- c("None", "0.3", "0.05")
panel_data <- function(d) {
  ref <- d %>% filter(method %in% c("LR_PS", "LR_Combat", "LR_PSspline")) %>%
    mutate(mtype = unname(fixed_lab[method])) %>% tidyr::crossing(cal_lab = cal3)
  gee <- d %>% filter(grepl(paste0("^(wrep|worep)_(nocal|cal0\\.3|cal0\\.05)_", CS, "$"), method)) %>%
    mutate(mtype = ifelse(grepl("^wrep", method), "GEE-EXW", "GEE-EXWO"),
           cal_lab = unname(cal_lab[sub("^(wrep|worep)_(.*)_(ind|ex)$", "\\2", method)]))
  bind_rows(ref, gee) %>% mutate(cal_lab = factor(cal_lab, levels = cal3), mtype = factor(mtype, levels = names(mt_col)))
}
sfig <- function(d, x, y, xlab, ylab, file) {
  p <- ggplot(d, aes(.data[[x]], .data[[y]], color = mtype, linetype = mtype, linewidth = mtype, group = mtype)) +
    geom_line() + geom_point(size = 1.8, show.legend = FALSE) +
    facet_grid(rate ~ cal_lab, labeller = labeller(cal_lab = function(v) paste0("cal=", v), rate = function(v) paste0("rate=", v))) +
    scale_color_manual(values = mt_col) + scale_linetype_manual(values = mt_lty) + scale_linewidth_manual(values = mt_lw) +
    labs(x = xlab, y = ylab, color = NULL, linetype = NULL, linewidth = NULL, caption = "N = 1000, Scenario 1") + base_theme +
    guides(color = guide_legend(nrow = 1))
  ggsave(file.path(fig_out, file), p, width = 10, height = 6, dpi = 300); message("  saved ", file)
}
if (!is.null(t1a)) {
  d <- panel_data(t1a %>% filter(N == 1000, scenario == "S1", subset == "all")) %>%
    mutate(alpha_lab = factor(sprintf("%.3f", alpha), levels = alpha_lev))
  sfig(d, "alpha_lab", "T1", "Nominal alpha", "Type-1 error", "SFig2.png")
}
if (!is.null(hser)) {
  d <- panel_data(hser %>% filter(subset == "all")) %>% mutate(h2_lab = factor(sprintf("%.3f", h2)))
  sfig(d, "h2_lab", "AdjPw", "Marker effect size (h2)", "Adjusted power", "SFig3.png")
}
message("16_make_figures_tables_v3.R complete.")
