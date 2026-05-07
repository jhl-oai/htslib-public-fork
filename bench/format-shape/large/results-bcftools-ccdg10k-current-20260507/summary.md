# Current bcftools FORMAT Planner Results

Input: `bench/format-shape/public/ccdg_chr22_10k.vcf.gz`

Binary: `Using htslib 1.23.1-58-gbd643182`

Each row is the mean of five streamed baseline/planned repetitions. Outputs were piped through `cksum`; all baseline-vs-plan checksum comparisons in `checks.tsv` are `ok`.

| Command | Threads | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 0 | 4.382 s | 4.032 s | 1.09x | 2.716 s | 2.412 s | 1.13x |
| `view_bcf` | 2 | 3.160 s | 2.810 s | 1.12x | 2.806 s | 2.480 s | 1.13x |
| `view_bcf` | 4 | 3.056 s | 2.852 s | 1.07x | 2.776 s | 2.480 s | 1.12x |
| `view_keep2_bcf` | 0 | 2.682 s | 2.414 s | 1.11x | 2.588 s | 2.296 s | 1.13x |
| `view_keep2_bcf` | 2 | 2.362 s | 2.060 s | 1.15x | 2.670 s | 2.376 s | 1.12x |
| `view_keep2_bcf` | 4 | 2.340 s | 2.046 s | 1.14x | 2.662 s | 2.376 s | 1.12x |
| `view_sites_bcf` | 0 | 2.718 s | 2.402 s | 1.13x | 2.600 s | 2.298 s | 1.13x |
| `view_sites_bcf` | 2 | 2.312 s | 2.008 s | 1.15x | 2.644 s | 2.340 s | 1.13x |
| `view_sites_bcf` | 4 | 2.310 s | 2.020 s | 1.14x | 2.634 s | 2.338 s | 1.13x |
| `query_sites` | 0 | 0.430 s | 0.430 s | 1.00x | 0.396 s | 0.392 s | 1.01x |
| `query_gt_keep2` | 0 | 0.716 s | 0.462 s | 1.55x | 0.668 s | 0.422 s | 1.58x |
| `query_gt_all` | 0 | 3.056 s | 2.748 s | 1.11x | 2.942 s | 2.634 s | 1.12x |
| `stats` | 0 | 0.470 s | 0.486 s | 0.97x | 0.432 s | 0.434 s | 1.00x |
| `filter_gt_keep2` | 0 | 2.994 s | 2.668 s | 1.12x | 2.872 s | 2.572 s | 1.12x |
| `filter_gt_keep2` | 2 | 2.690 s | 2.330 s | 1.15x | 2.988 s | 2.654 s | 1.13x |
| `filter_gt_keep2` | 4 | 2.632 s | 2.330 s | 1.13x | 2.952 s | 2.662 s | 1.11x |
