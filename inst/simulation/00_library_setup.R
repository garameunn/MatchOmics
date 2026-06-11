## ============================================================
##  00_library_setup.R
##  Changes vs original:
##  1. Bioconductor packages added (phyloseq, microbiome, edgeR, sva)
##  2. Auto-install for both CRAN and Bioconductor packages
##  3. Package load verification with clear error messages
## ============================================================

# library(Matrix, lib.loc = "/home2/nekim/Rpackage/") # 호환성 문제

# ── CRAN packages ───────────────────────────────────────────
cran_packages <- c(
  "MatchIt", "ggplot2", "ggpubr", "gridExtra",
  "data.table", "parallel", "tidyr", "doParallel", "foreach",
  "reshape2", "survival", "dplyr", "MASS", "readxl",
  "Matrix", "MatrixExtra", "geeM", "limma", "openxlsx",
  "vegan", "igraph"
)

# ── Bioconductor packages ───────────────────────────────────
bioc_packages <- c(
  "phyloseq",    # OTU table handling
  "microbiome",  # microbiome utilities
  "edgeR",       # cpm() function
  "sva",         # ComBat batch correction
  "limma"        # removeBatchEffect
)

# ── install BiocManager if needed ───────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# ── install & load CRAN ─────────────────────────────────────
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing CRAN package: %s", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ── install & load Bioconductor ─────────────────────────────
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing Bioconductor package: %s", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

message("All packages loaded successfully.")
