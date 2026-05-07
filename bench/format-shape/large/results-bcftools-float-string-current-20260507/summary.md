# Current bcftools FORMAT Planner Negative-Control Results

Input: `bench/format-shape/large/synthetic/large_float_string_2048s.vcf.gz`

Binary: `Using htslib 1.23.1-58-gbd643182`

Each row is the mean of five streamed baseline/planned repetitions. Outputs were piped through `cksum`; all baseline-vs-plan checksum comparisons in `checks.tsv` are `ok`.

| Command | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 4.718 s | 4.584 s | 1.03x | 3.188 s | 3.146 s | 1.01x |
| `view_keep2_bcf` | 3.162 s | 3.190 s | 0.99x | 3.038 s | 3.066 s | 0.99x |
| `view_sites_bcf` | 3.322 s | 3.170 s | 1.05x | 3.024 s | 3.010 s | 1.00x |
| `query_sites` | 0.566 s | 0.564 s | 1.00x | 0.510 s | 0.510 s | 1.00x |
| `query_gt_keep2` | 1.064 s | 1.042 s | 1.02x | 0.984 s | 0.972 s | 1.01x |
| `query_gt_all` | 3.594 s | 4.088 s | 0.88x | 3.390 s | 3.394 s | 1.00x |
| `stats` | 0.560 s | 0.562 s | 1.00x | 0.504 s | 0.508 s | 0.99x |
| `filter_gt_keep2` | 3.684 s | 3.446 s | 1.07x | 3.372 s | 3.336 s | 1.01x |
