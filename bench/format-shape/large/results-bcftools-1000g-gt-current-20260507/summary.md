# Current bcftools FORMAT Planner Results

Input: `bench/format-shape/large/public/1000g_chr22_full_genotypes.vcf.gz`

Binary: `Using htslib 1.23.1-58-gbd643182`

Each row is the mean of five streamed baseline/planned repetitions. Outputs were piped through `cksum`; all baseline-vs-plan checksum comparisons in `checks.tsv` are `ok`.

| Command | Threads | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 0 | 33.056 s | 15.922 s | 2.08x | 26.946 s | 10.316 s | 2.61x |
| `view_bcf` | 2 | 24.924 s | 11.700 s | 2.13x | 27.088 s | 10.662 s | 2.54x |
| `view_bcf` | 4 | 24.674 s | 11.708 s | 2.11x | 26.886 s | 10.846 s | 2.48x |
| `view_keep2_bcf` | 0 | 27.076 s | 10.522 s | 2.57x | 26.680 s | 10.260 s | 2.60x |
| `view_keep2_bcf` | 2 | 25.132 s | 8.850 s | 2.84x | 27.016 s | 10.768 s | 2.51x |
| `view_keep2_bcf` | 4 | 25.130 s | 8.734 s | 2.88x | 27.034 s | 10.662 s | 2.54x |
| `view_sites_bcf` | 0 | 26.262 s | 9.984 s | 2.63x | 25.894 s | 9.704 s | 2.67x |
| `view_sites_bcf` | 2 | 24.432 s | 8.310 s | 2.94x | 26.346 s | 10.244 s | 2.57x |
| `view_sites_bcf` | 4 | 24.418 s | 8.242 s | 2.96x | 26.410 s | 10.150 s | 2.60x |
| `query_sites` | 0 | 2.542 s | 2.524 s | 1.01x | 2.466 s | 2.456 s | 1.00x |
| `query_gt_keep2` | 0 | 5.836 s | 2.998 s | 1.95x | 5.728 s | 2.888 s | 1.98x |
| `query_gt_all` | 0 | 58.820 s | 42.396 s | 1.39x | 57.230 s | 40.938 s | 1.40x |
| `stats` | 0 | 2.970 s | 2.970 s | 1.00x | 2.876 s | 2.872 s | 1.00x |
| `filter_gt_keep2` | 0 | 49.376 s | 33.010 s | 1.50x | 48.782 s | 32.420 s | 1.50x |
| `filter_gt_keep2` | 2 | 48.132 s | 31.410 s | 1.53x | 49.516 s | 33.116 s | 1.50x |
| `filter_gt_keep2` | 4 | 49.984 s | 33.958 s | 1.47x | 51.412 s | 35.188 s | 1.46x |
