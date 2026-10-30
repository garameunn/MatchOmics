# Revision analyses (Briefings in Bioinformatics, BIB-26-1896, revision 1)

Scripts used for the simulations, real-data analyses, tables and figures of the revised manuscript
"MatchOmics: a propensity score matching framework for heterogeneity-robust multi-omics association analysis".
The MatchOmics R package itself is in the root of this repository.

Individual-level data (SNUHGC microbiome data and FROM proteomic/metabolomic data) are not included;
access is subject to the original studies' data-use agreements.

## simulation/

| Script | Purpose |
|---|---|
| `00_library_setup.R`, `99_utils.R` | packages and helper functions |
| `01_PreliminaryRobject.R` | taxa filtering (detected in >= 50% of subjects and relative abundance >= 0.1%), diversity indices, PS-difference ranking used for the top 20 taxa |
| `02_generator_v3.R` | data generation and analysis for one replicate: log-CPM markers; liability-threshold outcome generated once on the full sample; two-round matching; weighted GEE; comparators LR, LR-PS, LR-COM, LR-PSspline, IPW, OW |
| `08_Sim_v3_server.R` | main simulation (Scenarios 1-3, N = 200/500/1,000, prevalence 0.1/0.2/0.3; 2,000 replicates for type-1 error, 1,000 for power) |
| `08_Sim_v3b_full.R` | doubly robust weighting comparators (IPW-REG, OW-REG) on the same simulated data |
| `13_Sim_v3_hseries.R` | Scenario 1 power across marker effect sizes (h^2 = 0.001-0.01; Supplementary Figure 3) |
| `09_check_v3.R`, `10_summarize_v3.R`, `11_summarize_v3b.R` | checks and summaries (type-1 error, raw and adjusted power, all taxa and top 20) |
| `12_diag_S3_ex.R` | diagnostic run for the exchangeable working correlation under Scenario 3 |
| `15_extract_fig_data_v3.R` | data for Figure 1 and Supplementary Figures 2-3 |

Set `base_dir` at the top of each script to the project folder (`code/`, `data/`, `results/`).

## real_data/

| Script | Purpose |
|---|---|
| `rev_common.R` | shared functions (data loading, PS from the first two PCs, two-round matching, weights, weighted SMD, GEE) |
| `revC_final.R` | **final analysis reported in the revised manuscript**: round 2 without replacement after re-estimating the PS, weights w = m / mean(m), weighted SMD balance criterion (< 0.2) |
| `revC_post.R`, `revC_supp.R` | concordant markers, raw vs FDR counts, excluded markers, caliper and sample-size tables |
| `rev04_cutoff_sensitivity.R` | tertile/quartile cutoffs instead of the median split |
| `rev05_order_permutation.R` | sensitivity to the greedy matching order (20 random orders) |
| `rev00_validate.R`, `rev01_fdr_tables.R`, `rev03_smd_samplesize.R`, `rev07_gee_sortfix.R`, `rev08_lr_recompute.R`, `revC_compare_envs.R` | validation against the original results and supporting checks |
| `make_fig3_C.R` | Figure 3 |

The default options in `rev_common.R` reproduce the analysis of the original submission;
the final specification is selected in `revC_final.R`. Set the environment variables `REV_RAW`
(raw omics data folder) and `REV_BASE` (results folder) before running.

## figures_tables/

`16_make_figures_tables_v3.R` (Figures 1, 3, Supplementary Figures 2-3, Tables 2-5),
`20b_fig2_matchomics_pair.R` (Figure 2 and Supplementary Figure 4),
`18_renumber_stables_citation_order.R`, `19_stable11_excluded_markers.R`, `20_fig2_S1S2.R` and `21_stable_readable_headers.R`
(Supplementary Tables 1-14). Run from the project root after the summaries above.

## Software

R 4.3.1 (tables, figures, real data) and R 4.4 (simulations); MatchIt, geeM, sva (ComBat), edgeR, splines, ggplot2, openxlsx.
