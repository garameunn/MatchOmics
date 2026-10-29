## make_fig3_C.R -- Fig 3 under specification C (same layout/code as make_fig3.R)
##   (A) protein Venn (LR-PS, LR-COM, GEE-EXW, GEE-EXWO; raw p < 0.05, C marker set)
##   (B) protein representative marker boxplot
##   (C) metabolite Venn   (D) metabolite representative marker boxplot
## Input : realdata/revision/C_final/{protein,metabolite}_final_C.csv
## Output: realdata/revision/C_final/Fig3_C*.png   (figures/Fig3.png untouched)
## Usage : Rscript make_fig3_C.R   (env REV_RAW, REV_BASE, REV_LIB)
if (nzchar(Sys.getenv("REV_LIB"))) .libPaths(c(Sys.getenv("REV_LIB"), .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(readxl); library(ggVennDiagram)
})
raw_dir    <- Sys.getenv("REV_RAW"); base_dir <- Sys.getenv("REV_BASE")
raw_dir_bc <- file.path(raw_dir, "BeforeBatchCorrection")
c_dir      <- file.path(base_dir, "realdata", "revision", Sys.getenv("REV_CDIR", "C_final"))
ALPHA <- 0.05

base_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "gray92"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
        axis.title = element_text(size = 10), plot.title = element_text(size = 11, face = "bold"),
        plot.tag = element_text(size = 13, face = "bold"))

cov_raw <- as.data.frame(read_excel(file.path(raw_dir, "Cleaned_Data_ADDeGFR_imputed_jm_fixdur.xlsx"), sheet = "Sheet 1"))
names(cov_raw)[tolower(names(cov_raw)) == "id"] <- "ID"
cov_raw$ID <- gsub("-", "_", cov_raw$ID)
cov_raw$s_area_bin <- as.integer(cov_raw$s_area == 1)
cov_exp <- cov_raw[, c("ID", "s_area_bin")]

make_venn_panel <- function(df, tag_label, fill_high) {
  sig <- function(col) df$marker[!is.na(df[[col]]) & df[[col]] < ALPHA]
  vl <- list("LR-PS" = sig("LR_PS_pval"), "LR-COM" = sig("LR_COM_pval"),
             "GEE-EXW" = sig("wrep_cal0.3_pval"), "GEE-EXWO" = sig("worep_cal0.3_pval"))
  cat(sprintf("[%s] sig: %s (of %d)\n", tag_label, paste(names(vl), lengths(vl), collapse = " "), nrow(df)))
  ggVennDiagram(x = vl, label = "count", label_alpha = 0, edge_size = 0.6) +
    scale_fill_gradient(low = "white", high = fill_high, name = "Count") +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(tag = tag_label) + base_theme +
    theme(legend.position = "right", axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())
}

make_box_panel <- function(mk, raw_data, tag_label, y_label_suffix = "(log2 intensity)") {
  df_raw <- raw_data[, c("ID", mk), drop = FALSE]; names(df_raw)[2] <- "value"
  df_box <- merge(df_raw, cov_exp, by = "ID")
  df_box <- df_box[!is.na(df_box$s_area_bin) & !is.na(df_box$value), ]
  df_box$group <- factor(ifelse(df_box$s_area_bin == 1, "Exposed", "Control"), levels = c("Control", "Exposed"))
  grp_n <- df_box %>% group_by(group) %>% summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0(as.character(group), "\n(n=", n, ")"))
  grp_cols <- c("Control" = "#4E9CC4", "Exposed" = "#D32F2F")
  ggplot(df_box, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65, size = 0.5, color = "gray30") +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, color = "gray25") +
    scale_fill_manual(values = grp_cols, guide = "none") +
    scale_x_discrete(labels = setNames(grp_n$label, grp_n$group)) +
    labs(x = "Exposure group", y = paste0(mk, " ", y_label_suffix), tag = tag_label) +
    base_theme + theme(panel.grid.major.x = element_blank())
}

pick_best <- function(fin) {
  fin$n_sig <- with(fin, (!is.na(LR_PS_pval) & LR_PS_pval < ALPHA) + (!is.na(LR_COM_pval) & LR_COM_pval < ALPHA) +
    (!is.na(wrep_cal0.3_pval) & wrep_cal0.3_pval < ALPHA) + (!is.na(worep_cal0.3_pval) & worep_cal0.3_pval < ALPHA))
  fin$marker[order(-fin$n_sig, fin$LR_PS_pval)][1]
}

prot_fin <- read.csv(file.path(c_dir, "protein_final_C.csv"), stringsAsFactors = FALSE)
meta_fin <- read.csv(file.path(c_dir, "metabolite_final_C.csv"), stringsAsFactors = FALSE)
prot_best <- pick_best(prot_fin); meta_best <- pick_best(meta_fin)
cat("representative markers:", prot_best, meta_best, "\n")

prot_raw <- read.delim(file.path(raw_dir_bc, "20250410_log2_norm_Proteomics_All.txt"), check.names = FALSE, stringsAsFactors = FALSE)
prot_raw$ID <- gsub("-", "_", prot_raw$ID)
meta_raw <- read.delim(file.path(raw_dir_bc, "20250428_log2_norm_Metabolite_All.txt"), check.names = FALSE, stringsAsFactors = FALSE)
meta_raw$ID <- gsub("-", "_", meta_raw$ID)

figA <- make_venn_panel(prot_fin, "(A)", "#2C4770"); figB <- make_box_panel(prot_best, prot_raw, "(B)")
figC <- make_venn_panel(meta_fin, "(C)", "#5A8A5E"); figD <- make_box_panel(meta_best, meta_raw, "(D)")
fig3 <- (figA + figB) / (figC + figD) + plot_layout(widths = c(1.4, 1))
ggsave(file.path(c_dir, "Fig3_C.png"), fig3, width = 9.5, height = 10.4, dpi = 300)
ggsave(file.path(c_dir, "Fig3_C_A_ProtVenn.png"), figA, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(c_dir, "Fig3_C_B_ProtBox.png"), figB, width = 4, height = 5.2, dpi = 300)
ggsave(file.path(c_dir, "Fig3_C_C_MetaVenn.png"), figC, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(c_dir, "Fig3_C_D_MetaBox.png"), figD, width = 4, height = 5.2, dpi = 300)
cat("saved Fig3_C*.png in", c_dir, "\n")
