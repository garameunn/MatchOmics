# Simulation scripts

These scripts reproduce the Type I error, power, and bias simulations reported
in the MatchOmics paper. They were developed for a Linux cluster (R 4.4) but
can be adapted to any environment with sufficient memory.

## ⚠️ S3 scenario: use the `*_S3fix_*` scripts only

The original S3 confounder-generation logic (`02_simulation_generator.R`,
inherited by `03_Type1Error_Calculating_server.R`,
`04_Power_Calculating_server.R`, `05_Bias_EffectSize_Calculating_server.R`)
had two bugs: a hardcoded confounder-taxa index, and a residualisation step
that made the confounder mathematically orthogonal to the matching variable —
destroying the intended confounding structure for **S3 only**.

The `two_round_match()` matching algorithm itself is identical in both
versions and was never affected — this was purely a synthetic-data-generation
bug in the S3 benchmark scenario.

- **S1 / S2 results**: use the original (non-`S3fix`) scripts — unaffected, valid as-is.
- **S3 results**: use the `*_S3fix_*` scripts exclusively (`02_simulation_generator_S3fix.R`,
  `03_Type1Error_S3fix_server.R`, `04_Power_S3fix_server.R`,
  `05_Bias_EffectSize_S3fix_server.R`). Discard any S3 output from the
  non-fix scripts.

## File overview

| File | Purpose |
|------|---------|
| `00_library_setup.R` | Install / load all required packages |
| `02_simulation_generator.R` | Core simulation functions — **valid for S1/S2 only** |
| `02_simulation_generator_S3fix.R` | Same, with corrected S3 confounder generation — **required for S3** |
| `03_Type1Error_Calculating_server.R` | Type I error simulation (5000 reps) — **valid for S1/S2 only** |
| `03_Type1Error_S3fix_server.R` | Type I error, S3 scenario only — **required for S3** |
| `04_Power_Calculating_server.R` | Power simulation (1000 reps) — **valid for S1/S2 only** |
| `04_Power_LargeB_server.R` | Power — large effect sizes (b = 0.03, 0.05) |
| `04_Power_S3fix_server.R` | Power, S3 scenario only — **required for S3** |
| `04_Power_WeightTest.R` | Weighting-scheme sensitivity check (auxiliary) |
| `05_Bias_EffectSize_Calculating_server.R` | Relative bias and effect size — **valid for S1/S2 only** |
| `05_Bias_EffectSize_Calculating_server_v2.R` | Revised bias/effect-size calculation (supersedes the above for S1/S2) |
| `05_Bias_EffectSize_S3fix_server.R` | Relative bias and effect size, S3 scenario only — **required for S3** |
| `06_Final_Stat_Organizing_server.R` | Organise results → xlsx |
| `07_Diagnostics_server.R` | SMD, Neff, round-2 diagnostics |
| `99_utils.R` | Shared utility functions |

## Directory structure expected by scripts

```
Matching2026/
  code/      ← scripts (this directory)
  data/      ← otulist, indiclist, psdifflist, seed, batchinfo (.RData)
  results/
    T1/
    Power/
    Organized/
    Diagnostics/
```

## Running

```bash
# Type I error (all scenarios, 3 cores)
Rscript 03_Type1Error_Calculating_server.R --cores 3

# Power (all scenarios)
Rscript 04_Power_Calculating_server.R --cores 3

# Post-processing (after 03 + 04 complete)
Rscript 05_Bias_EffectSize_Calculating_server.R
Rscript 06_Final_Stat_Organizing_server.R
Rscript 07_Diagnostics_server.R
```

## Simulation design

| Parameter | Values |
|-----------|--------|
| Scenarios | S1 (PS-based confounding) / S2 (direct confounding) / S3 (multi-taxon PS) |
| Sample sizes | 200 / 500 / 1000 |
| Prevalence | 0.1 / 0.2 / 0.3 |
| Caliper | NULL / 0.3 / 0.2 / 0.1 / 0.05 |
| GEE corstr | independence / exchangeable |
| T1 reps | 5000 |
| Power reps | 1000 |
