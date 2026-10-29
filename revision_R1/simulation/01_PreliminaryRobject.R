## ============================================================
##  01_PreliminaryRobject.R
##  Changes vs original:
##  1. seed generation added (최초 1회 실행 후 seed.RData 저장)
##  2. paths unified to relative (../data/)
##  3. tree file removed (not needed)
##  4. Bioconductor packages handled via 00_library_setup.R
##  NOTE: Run this only once to generate RData objects.
##        On server, otulist/indiclist/psdifflist already exist.
##        Only generate_seed.R needs to be run on server.
## ============================================================

source("00_library_setup.R")
source("99_utils.R")

# ── data paths ──────────────────────────────────────────────
data_dir <- "../data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ── load raw OTU data ────────────────────────────────────────
physeq <- as.data.frame(data.table::fread(
  file.path(data_dir, "otu_table_mc2_w.txt")))

otu <- phyloseq::otu_table(
  physeq[, 2:(ncol(physeq)-1)], taxa_are_rows = TRUE)
rownames(otu) <- physeq[, 1]

tax <- phyloseq::tax_table(as.matrix(data.frame(
  tx = physeq[, ncol(physeq)]) |>
  tidyr::separate(tx,
    c("Domain","Phylum","Class","Order","Family","Genus","Species"), ";")))
rownames(tax) <- physeq[, 1]

phy <- phyloseq::phyloseq(phyloseq::otu_table(otu), phyloseq::tax_table(tax))

# ── subset by sample size ────────────────────────────────────
make_subset <- function(phy, n) {
  toysam <- c(rep(TRUE, n), rep(FALSE, phyloseq::nsamples(phy) - n))
  sub    <- phyloseq::prune_samples(toysam, phy)
  list(
    otutable  = filtOTU(phyloseq::otu_table(sub)),
    indicator = indc(phyloseq::otu_table(sub), omics = "Metagenomics")
  )
}

s1000 <- make_subset(phy, 1000)
s500  <- make_subset(phy, 500)
s200  <- make_subset(phy, 200)

lbsize <- apply(phyloseq::otu_table(
  phyloseq::prune_samples(
    c(rep(TRUE, 1000), rep(FALSE, phyloseq::nsamples(phy) - 1000)), phy)),
  2, sum)

# otulist[[1]]=1000, [[2]]=500, [[3]]=200, [[4]]=lbsize(1000 기준)
otulist  <- list(s1000$otutable, s500$otutable, s200$otutable, lbsize)
indiclist <- list(s1000$indicator, s500$indicator, s200$indicator)

save(otulist,  file = file.path(data_dir, "otulist.RData"))
save(indiclist, file = file.path(data_dir, "indiclist.RData"))
message("otulist, indiclist saved.")

# ── psdiff (PS difference between up/down groups per taxa) ───
compute_psdiff <- function(otutable, indicator) {
  unlist(lapply(data.frame(t(otutable)), psdiffun, k = indicator))
}

psdiff1000 <- compute_psdiff(s1000$otutable, s1000$indicator)
psdiff500  <- compute_psdiff(s500$otutable,  s500$indicator)
psdiff200  <- compute_psdiff(s200$otutable,  s200$indicator)

psdifflist <- list(psdiff1000, psdiff500, psdiff200)
save(psdifflist, file = file.path(data_dir, "psdifflist.RData"))
message("psdifflist saved.")

# ── seed generation (최초 1회만 실행) ────────────────────────
# rep_max × dimm_max 보다 충분히 크게 생성
# T1:    5000 rep × max(nrow(otutable)) taxa
# Power: 1000 rep × max(nrow(otutable)) taxa
# → 5000 × nrow(otulist[[1]]) 로 설정

dimm_max <- nrow(otulist[[1]])   # 1000 샘플 기준 taxa 수
rep_max  <- 5000

set.seed(42)
seeds <- sample.int(.Machine$integer.max,
                    size = rep_max * dimm_max)

save(seeds, file = file.path(data_dir, "seed.RData"))
message(sprintf("seeds saved: %d seeds (rep=%d × dimm=%d)",
                length(seeds), rep_max, dimm_max))

message("01_PreliminaryRobject.R complete.")
