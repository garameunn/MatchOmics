pkg <- "C:/Users/jhnam/OneDrive/문서/OneDrive/매칭방법론_Binary/260330_알고리즘변경/MatchOmics"

for (f in list.files(file.path(pkg, "R"), full.names = TRUE)) source(f)

load(file.path(pkg, "data/toy_omics.rda"))
load(file.path(pkg, "data/toy_outcome.rda"))
load(file.path(pkg, "data/toy_heterogeneity.rda"))

fit <- MatchOmics(toy_omics[1, ], toy_outcome, toy_heterogeneity,
                  caliper = 0.3, method = "two_round")
cat("two_round : n_matched =", fit$n_matched, "| Neff =", round(fit$neff, 1), "\n")
print(fit)

fit2 <- MatchOmics(toy_omics[1, ], toy_outcome, toy_heterogeneity,
                   caliper = 0.3, method = "standard")
cat("standard  : n_matched =", fit2$n_matched, "\n")

cat("\nAll OK\n")
