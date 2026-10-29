## ============================================================
##  20_fig2_S1S2.R — Figure 2 with Scenario 2 added (decided 2026-10-29)
##  Figure 2      : N = 1,000, rows = Scenario 1 / 2, columns = all taxa / top 20
##  Suppl. Fig. 4 : the same for N = 200, 500 and 1,000 (former SFig 4 was deleted, R3-5)
##  Suppl. Tables 5-6 : Scenario 2 rows added below Scenario 1 (top 20 taxa)
##  Inputs : 결과파일_20260517/v3_summary/summary_v3_long.csv (as in 16_make_figures_tables_v3.R)
##           manuscript/Tables_rev/STables_final.xlsx (from 18 and 19)
##  Outputs: manuscript/Figures_rev/Fig2_S1S2.png, SFig4_S1S2.png (Fig2.png is kept)
##           manuscript/Tables_rev/STables_final.xlsx, sheets "STable 5" and "STable 6" replaced
##           (previous workbook kept as STables_final_before_S2top20.xlsx)
##  NOTE: the figures of this script were superseded by 20b_fig2_matchomics_pair.R (final Fig2.png, SFig4.png)
##        and removed; this script is kept for Supplementary Tables 5-6.
##  Run from the project root:
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\20_fig2_S1S2.R
## ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(openxlsx) })

CS <- "ind"   # working correlation of the GEE curves, as in 16_make_figures_tables_v3.R
fig_out <- "manuscript/Figures_rev"
book <- "manuscript/Tables_rev/STables_final.xlsx"
allm <- read.csv("결과파일_20260517/v3_summary/summary_v3_long.csv", check.names = FALSE)
lb   <- read.csv("결과파일_20260517/v3_summary/summary_v3b_long.csv", check.names = FALSE)
allm <- bind_rows(allm, filter(lb, method %in% c("IPW_DR", "OW_DR")))

## ── figures ───────────────────────────────────────────────────
m_cs <- function(x) paste0(x, "_", CS)
f2_m   <- c("LR_PS", "LR_PSspline", m_cs("wrep_nocal"), m_cs("wrep_cal0.3"), m_cs("worep_cal0.3"))
f2_lab <- setNames(c("LR-PS", "LR-PSspline", "GEE-EXW (no caliper)", "GEE-EXW (cal = 0.3)", "GEE-EXWO (cal = 0.3)"), f2_m)
f2_col <- setNames(c("#222222", "#0072B2", "#E69F00", "#D32F2F", "#009E73"), f2_m)
f2_lty <- setNames(c("solid", "dotdash", "dashed", "solid", "dashed"), f2_m)
f2_lw  <- setNames(c(0.53, 0.6, 0.53, 1.07, 0.53), f2_m)
scen_lab <- c(S1 = "Scenario 1 (PS-specified)", S2 = "Scenario 2 (PS-misspecified)")
base_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9), legend.key.width = unit(1.8, "cm"),
        strip.background = element_rect(fill = "gray92"), strip.text = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(), axis.text = element_text(size = 8))

fd <- allm %>% filter(scenario %in% c("S1", "S2"), method %in% f2_m) %>%
  mutate(method = factor(method, levels = f2_m),
         scen = factor(scen_lab[scenario], levels = scen_lab),
         N_lab = factor(paste0("N = ", N), levels = c("N = 200", "N = 500", "N = 1000")),
         subset = factor(ifelse(subset == "all", "All taxa", "Top 20"), levels = c("All taxa", "Top 20")))
pw_plot <- function(d) {
  ggplot(d, aes(rate, AdjPw, color = method, linetype = method, linewidth = method, group = method)) +
    geom_line() + geom_point(size = 2.2, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
    scale_color_manual(values = f2_col, labels = f2_lab) + scale_linetype_manual(values = f2_lty, labels = f2_lab) +
    scale_linewidth_manual(values = f2_lw, labels = f2_lab) +
    labs(x = "Prevalence rate", y = "Adjusted power", color = NULL, linetype = NULL, linewidth = NULL) +
    base_theme + guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), linewidth = guide_legend(nrow = 2))
}
p2 <- pw_plot(filter(fd, N == 1000)) + facet_grid(scen ~ subset)
ggsave(file.path(fig_out, "Fig2_S1S2.png"), p2, width = 7.5, height = 6.5, dpi = 300); message("  saved Fig2_S1S2.png")
ps4 <- pw_plot(fd) + facet_grid(N_lab ~ scen + subset)
ggsave(file.path(fig_out, "SFig4_S1S2.png"), ps4, width = 11, height = 8, dpi = 300); message("  saved SFig4_S1S2.png")

## ── Supplementary Tables 5-6: add Scenario 2 rows (top 20 taxa) ─
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
wide_tab <- function(scen, val, ss = "top20", Ns = c(1000, 500, 200)) {
  d <- allm %>% filter(scenario == scen, subset == ss, method %in% row_order, N %in% Ns) %>%
    mutate(col = sprintf("N = %d|%.1f", N, rate))
  w <- d %>% select(method, col, value = all_of(val)) %>% pivot_wider(names_from = col, values_from = value)
  cols <- as.vector(outer(sprintf("%.1f", c(0.1, 0.2, 0.3)), sprintf("N = %d", Ns), function(r, n) paste0(n, "|", r)))
  w <- w[match(intersect(row_order, w$method), w$method), c("method", intersect(cols, names(w)))]
  w$method <- lab_method(w$method); names(w)[1] <- "Method|Prevalence"; as.data.frame(w)
}
two_block <- function(val) {
  a <- wide_tab("S1", val); b <- wide_tab("S2", val)
  hdr <- function(lbl) { h <- a[1, ]; h[1, ] <- NA; h[1, 1] <- lbl; h }
  rbind(hdr(scen_lab[["S1"]]), a, hdr(scen_lab[["S2"]]), b)
}
abbr <- paste("LR, unadjusted logistic regression; LR-PS, logistic regression adjusted for the propensity score (PS);",
              "LR-COM, logistic regression on ComBat-corrected data; LR-PSspline, LR-PS with a natural cubic spline (4 df) of the logit PS;",
              "IPW/OW, inverse-probability/overlap weighting; IPW-REG/OW-REG, weighting combined with regression on the PS spline;",
              "GEE-INDW/EXW, with-replacement (two-round) matching analysed by weighted GEE with independence/exchangeable working correlation;",
              "GEE-INDWO/EXWO, the same for without-replacement matching; Cal, caliper in SD units of the PS.")
bold <- createStyle(textDecoration = "bold")
put_sheet <- function(wb, name, title, legend, df) {
  removeWorksheet(wb, name); addWorksheet(wb, name)
  writeData(wb, name, title, startRow = 1); addStyle(wb, name, bold, rows = 1, cols = 1)
  writeData(wb, name, legend, startRow = 2)
  df[-1] <- lapply(df[-1], function(x) round(x, 3))
  top <- sub("\\|.*", "", names(df)); bot <- sub(".*\\|", "", names(df)); top[duplicated(top)] <- ""
  writeData(wb, name, t(top), startRow = 4, colNames = FALSE)
  writeData(wb, name, t(bot), startRow = 5, colNames = FALSE)
  addStyle(wb, name, bold, rows = 4:5, cols = seq_along(top), gridExpand = TRUE)
  writeData(wb, name, df, startRow = 6, colNames = FALSE)
  addStyle(wb, name, bold, rows = 5 + which(df[[1]] %in% scen_lab), cols = 1, gridExpand = TRUE)
  setColWidths(wb, name, cols = 1, widths = 28)
}
file.copy(book, "manuscript/Tables_rev/STables_final_before_S2top20.xlsx", overwrite = TRUE)
wb <- loadWorkbook(book)
put_sheet(wb, "STable 5", "Supplementary Table 5. Type-1 error under Scenarios 1 and 2 for the 20 most confounded taxa.",
          paste("Taxa with the largest PS difference between matching groups (the same 20 taxa in both scenarios). Layout as in Table 2.", abbr),
          two_block("T1"))
put_sheet(wb, "STable 6", "Supplementary Table 6. Adjusted power under Scenarios 1 and 2 for the 20 most confounded taxa.",
          paste("Same taxa as Supplementary Table 5, h2 = 0.01. Layout as in Table 3.", abbr), two_block("AdjPw"))
worksheetOrder(wb) <- order(as.integer(sub("^STable ", "", names(wb))))
saveWorkbook(wb, book, overwrite = TRUE)
message("  STables 5-6 now include Scenario 2 (top 20 taxa)")
