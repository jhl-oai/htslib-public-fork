# Current bcftools Command Smoke Results

Single-pass representative command benchmark with `bcftools 1.23.1-63-gd9375776` using `htslib 1.23.1-58-gbd643182`. Outputs were compared with `cmp`; `checks.tsv` has only `ok` and expected `skipped_no_samples` rows.

| Input | Command | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---|---:|---:|---:|---:|---:|---:|
| `ccdg_10k` | `view_bcf` | 2.870 s | 2.560 s | 1.12x | 2.590 s | 2.300 s | 1.13x |
| `ccdg_10k` | `view_sites` | 2.720 s | 2.490 s | 1.09x | 2.580 s | 2.290 s | 1.13x |
| `ccdg_10k` | `query_sites` | 0.430 s | 0.500 s | 0.86x | 0.400 s | 0.390 s | 1.03x |
| `ccdg_10k` | `query_format` | 0.700 s | 0.450 s | 1.56x | 0.670 s | 0.420 s | 1.60x |
| `ccdg_10k` | `stats` | 0.470 s | 0.460 s | 1.02x | 0.430 s | 0.430 s | 1.00x |
| `ccdg_10k` | `filter_gt` | 2.980 s | 2.650 s | 1.12x | 2.860 s | 2.560 s | 1.12x |
| `1000g_chr22_full_genotypes` | `view_bcf` | 28.430 s | 11.260 s | 2.52x | 26.440 s | 9.980 s | 2.65x |
| `1000g_chr22_full_genotypes` | `view_sites` | 27.100 s | 10.000 s | 2.71x | 26.230 s | 9.690 s | 2.71x |
| `1000g_chr22_full_genotypes` | `query_sites` | 2.570 s | 2.580 s | 1.00x | 2.470 s | 2.470 s | 1.00x |
| `1000g_chr22_full_genotypes` | `query_format` | 5.950 s | 3.070 s | 1.94x | 5.740 s | 2.920 s | 1.97x |
| `1000g_chr22_full_genotypes` | `stats` | 3.050 s | 3.030 s | 1.01x | 2.930 s | 2.930 s | 1.00x |
| `1000g_chr22_full_genotypes` | `filter_gt` | 50.180 s | 33.420 s | 1.50x | 49.170 s | 32.670 s | 1.51x |
| `large_reordered_likelihood_2048s` | `view_bcf` | 3.260 s | 2.720 s | 1.20x | 2.980 s | 2.500 s | 1.19x |
| `large_reordered_likelihood_2048s` | `view_sites` | 3.120 s | 2.610 s | 1.20x | 2.980 s | 2.480 s | 1.20x |
| `large_reordered_likelihood_2048s` | `query_sites` | 0.830 s | 0.820 s | 1.01x | 0.760 s | 0.760 s | 1.00x |
| `large_reordered_likelihood_2048s` | `query_format` | 1.200 s | 0.860 s | 1.40x | 1.110 s | 0.770 s | 1.44x |
| `large_reordered_likelihood_2048s` | `stats` | 0.820 s | 0.830 s | 0.99x | 0.760 s | 0.760 s | 1.00x |
| `large_reordered_likelihood_2048s` | `filter_gt` | 3.550 s | 3.120 s | 1.14x | 3.420 s | 3.000 s | 1.14x |
| `large_float_string_2048s` | `view_bcf` | 3.380 s | 3.540 s | 0.95x | 3.000 s | 3.070 s | 0.98x |
| `large_float_string_2048s` | `view_sites` | 3.120 s | 3.060 s | 1.02x | 2.960 s | 2.940 s | 1.01x |
| `large_float_string_2048s` | `query_sites` | 0.550 s | 0.540 s | 1.02x | 0.500 s | 0.500 s | 1.00x |
| `large_float_string_2048s` | `query_format` | 1.020 s | 1.210 s | 0.84x | 0.960 s | 1.020 s | 0.94x |
| `large_float_string_2048s` | `stats` | 0.610 s | 0.620 s | 0.98x | 0.510 s | 0.520 s | 0.98x |
| `large_float_string_2048s` | `filter_gt` | 3.600 s | 3.480 s | 1.03x | 3.360 s | 3.360 s | 1.00x |
| `gnomad_sites_chr22` | `view_bcf` | 0.620 s | 0.580 s | 1.07x | 0.440 s | 0.450 s | 0.98x |
| `gnomad_sites_chr22` | `view_sites` | 0.520 s | 0.600 s | 0.87x | 0.430 s | 0.450 s | 0.96x |
| `gnomad_sites_chr22` | `query_sites` | 0.150 s | 0.180 s | 0.83x | 0.140 s | 0.140 s | 1.00x |
| `gnomad_sites_chr22` | `stats` | 0.460 s | 0.490 s | 0.94x | 0.440 s | 0.460 s | 0.96x |
