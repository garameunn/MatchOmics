# library setup

required_packages <- c("magic", "MatchIt", "microbiome","phyloseq", "ggplot2", "ggpubr", "gridExtra", "ape",
"picante","data.table","parallel","edgeR","tidyr","doParallel","foreach","reshape2","survival",
"dplyr","MASS","readxl","Matrix","MatrixExtra","geeM","limma","openxlsx")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
