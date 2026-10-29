## ============================================================
##  20b_fig2_matchomics_pair.R — Figure 2 / Suppl. Fig. 4, MatchOmics shown as a pair (decided 2026-10-29)
##  Same data and panels as 20_fig2_S1S2.R; only the styling differs:
##  the two MatchOmics settings recommended in the paper (0.3 caliper when the PS is adequate,
##  no caliper when it may be misspecified) are drawn with equal weight in red (solid / dashed),
##  the comparators LR-PS and LR-PSspline as thin lines. GEE-EXWO is left to the tables.
##  Figures only (Supplementary Tables 5-6 were updated by 20_fig2_S1S2.R).
##  Outputs: manuscript/Figures_rev/Fig2.png, SFig4.png (final; supersede Fig2.png of 16_make_figures_tables_v3.R
##           and the figures of 20_fig2_S1S2.R, which were removed)
##  Run from the project root:
##    & "C:\R\R-4.3.1\bin\Rscript.exe" code\20b_fig2_matchomics_pair.R
## ============================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

CS <- "ind"   # working correlation of the GEE curves, as in 16_make_figures_tables_v3.R
fig_out <- "manuscript/Figures_rev"
allm <- read.csv("결과파일_20260517/v3_summary/summary_v3_long.csv", check.names = FALSE)

m_cs <- function(x) paste0(x, "_", CS)
f2_m   <- c(m_cs("wrep_cal0.3"), m_cs("wrep_nocal"), "LR_PSspline", "LR_PS")
f2_lab <- setNames(c("GEE-EXW (cal = 0.3)", "GEE-EXW (no caliper)", "LR-PSspline", "LR-PS"), f2_m)
f2_col <- setNames(c("#D32F2F", "#D32F2F", "#0072B2", "#222222"), f2_m)
f2_lty <- setNames(c("solid", "22", "dotdash", "solid"), f2_m)
f2_lw  <- setNames(c(1.1, 1.1, 0.55, 0.55), f2_m)
f2_pt  <- setNames(c(2.4, 2.4, 1.6, 1.6), f2_m)
scen_lab <- c(S1 = "Scenario 1 (PS-specified)", S2 = "Scenario 2 (PS-misspecified)")
base_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9), legend.key.width = unit(1.8, "cm"),
        strip.background = element_rect(fill = "gray92"), strip.text = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(), axis.text = element_text(size = 8))

fd <- allm %>% filter(scenario %in% c("S1", "S2"), method %in% f2_m) %>%
  mutate(method = factor(method, levels = f2_m),           # thin comparator lines drawn last, so LR-PS stays visible where it overlaps no-caliper matching
         scen = factor(scen_lab[scenario], levels = scen_lab),
         N_lab = factor(paste0("N = ", N), levels = c("N = 200", "N = 500", "N = 1000")),
         subset = factor(ifelse(subset == "all", "All taxa", "Top 20"), levels = c("All taxa", "Top 20")))
pw_plot <- function(d) {
  ggplot(d, aes(rate, AdjPw, color = method, linetype = method, linewidth = method, group = method)) +
    geom_line() + geom_point(aes(size = method), show.legend = FALSE) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
    scale_color_manual(values = f2_col, labels = f2_lab, breaks = f2_m) +
    scale_linetype_manual(values = f2_lty, labels = f2_lab, breaks = f2_m) +
    scale_linewidth_manual(values = f2_lw, labels = f2_lab, breaks = f2_m) +
    scale_size_manual(values = f2_pt, guide = "none") +
    labs(x = "Prevalence rate", y = "Adjusted power", color = NULL, linetype = NULL, linewidth = NULL) +
    base_theme + guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), linewidth = guide_legend(nrow = 1))
}
p2 <- pw_plot(filter(fd, N == 1000)) + facet_grid(scen ~ subset)
ggsave(file.path(fig_out, "Fig2.png"), p2, width = 7.5, height = 6.5, dpi = 300); message("  saved Fig2.png")
ps4 <- pw_plot(fd) + facet_grid(N_lab ~ scen + subset)
ggsave(file.path(fig_out, "SFig4.png"), ps4, width = 11, height = 8, dpi = 300); message("  saved SFig4.png")
