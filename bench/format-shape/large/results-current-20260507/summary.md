# Current FORMAT Planner Large Corpus Results (2026-05-07)

Current-state `test/test_view -b -l 0` benchmark for product commit `bd643182c8fa722abbc0cb89860263a90bb97020` (`feature/vcf-parsing-speedup`). Each row is the mean of five sequential baseline/planned repetitions. Baseline uses `HTS_VCF_FORMAT_PLAN=0`; planned uses `HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1`.

All 50 baseline-vs-plan BCF comparisons in `checks.tsv` are `ok`.

| Input | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup | Plan hits/fallback |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ccdg_10k` | 2.928 s | 2.682 s | 1.09x | 2.640 s | 2.356 s | 1.12x | 9861 / 139 |
| `1000g_chr22_full_genotypes` | 28.092 s | 11.980 s | 2.34x | 26.174 s | 10.750 s | 2.43x | 1103547 / 0 |
| `large_ccdg_likelihood_2048s` | 4.576 s | 4.148 s | 1.10x | 4.208 s | 3.824 s | 1.10x | 20000 / 0 |
| `large_reordered_likelihood_2048s` | 3.252 s | 2.744 s | 1.19x | 3.022 s | 2.532 s | 1.19x | 20000 / 0 |
| `large_multiallelic_likelihood_2048s` | 3.492 s | 3.116 s | 1.12x | 3.238 s | 2.860 s | 1.13x | 16000 / 0 |
| `large_float_string_2048s` | 3.386 s | 3.342 s | 1.01x | 3.026 s | 3.022 s | 1.00x | 0 / 16000 |
| `large_phase_width_variation_2048s` | 2.982 s | 2.834 s | 1.05x | 2.680 s | 2.536 s | 1.06x | 12000 / 0 |
| `large_mixed_likelihood_2048s` | 2.424 s | 2.116 s | 1.15x | 2.238 s | 1.952 s | 1.15x | 12000 / 0 |
| `large_gt_first_reordered_2048s` | 1.938 s | 1.612 s | 1.20x | 1.778 s | 1.488 s | 1.19x | 12000 / 0 |
| `large_two_string_float_2048s` | 2.608 s | 2.608 s | 1.00x | 2.344 s | 2.340 s | 1.00x | 0 / 12000 |

Standard errors for real/user time are recorded in `summary.tsv`; raw per-repetition rows are in `timings.tsv`.
