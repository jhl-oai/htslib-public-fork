# bcftools 1000G GT-Only Threaded Ad Hoc Benchmark

This run measures realistic bcftools operations on the full 1000 Genomes Phase
3 chromosome 22 genotype VCF.  The input is a BGZF-compressed VCF text file with
2,504 samples and 205,612,353 compressed bytes, dominated by `GT` FORMAT rows.

The benchmark compares the generic parser (`HTS_VCF_FORMAT_PLAN=0`) against the
planned FORMAT parser (`HTS_VCF_FORMAT_PLAN=1`) using five repetitions.  Output
streams were piped through `cksum`; all 80 baseline-vs-plan checksum comparisons
were `ok`.

`--threads` was applied only to `bcftools view` command shapes.  `bcftools query`
and `bcftools stats` rows are unthreaded and therefore use `threads=0`.  This
run was resumed once with `RESUME=1`; the completed rows were preserved and the
missing tail was appended.

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
INPUT=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/htslib-vcf-avx-sanity/bench/format-shape/large/public/1000g_chr22_full_genotypes.vcf.gz \
NAME=1000g_chr22_full_gt \
OUTDIR=bench/format-shape/large/results-bcftools-1000g-gt-ad-hoc-threaded \
REPS=5 \
THREADS_LIST="0 2 4" \
  bash bench/format-shape/scripts/run_bcftools_ad_hoc_bench.sh
```

Resume command used after interruption:

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
INPUT=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/htslib-vcf-avx-sanity/bench/format-shape/large/public/1000g_chr22_full_genotypes.vcf.gz \
NAME=1000g_chr22_full_gt \
OUTDIR=bench/format-shape/large/results-bcftools-1000g-gt-ad-hoc-threaded \
REPS=5 \
THREADS_LIST="0 2 4" \
RESUME=1 \
  bash bench/format-shape/scripts/run_bcftools_ad_hoc_bench.sh
```

The bcftools binary reported `bcftools 1.23.1-63-gd9375776` using
`htslib 1.23.1-56-g8c809baf`.  Times are mean +/- standard error across five
sequential runs.

| Command | Threads | Baseline user | Planned user | User speedup | Baseline real | Planned real | Real speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 0 | 26.384 +/- 0.113 s | 10.310 +/- 0.036 s | 2.56x | 32.446 +/- 0.225 s | 15.942 +/- 0.040 s | 2.04x |
| `view_bcf` | 2 | 28.112 +/- 0.148 s | 11.952 +/- 0.091 s | 2.35x | 25.744 +/- 0.133 s | 12.344 +/- 0.028 s | 2.09x |
| `view_bcf` | 4 | 28.202 +/- 0.043 s | 11.878 +/- 0.058 s | 2.37x | 25.782 +/- 0.055 s | 12.328 +/- 0.024 s | 2.09x |
| `view_keep2_bcf` | 0 | 26.514 +/- 0.096 s | 10.454 +/- 0.106 s | 2.54x | 27.170 +/- 0.141 s | 11.236 +/- 0.492 s | 2.42x |
| `view_keep2_bcf` | 2 | 27.956 +/- 0.091 s | 11.956 +/- 0.037 s | 2.34x | 25.938 +/- 0.138 s | 9.688 +/- 0.022 s | 2.68x |
| `view_keep2_bcf` | 4 | 28.030 +/- 0.087 s | 11.812 +/- 0.077 s | 2.37x | 25.918 +/- 0.090 s | 9.522 +/- 0.070 s | 2.72x |
| `view_sites_bcf` | 0 | 25.734 +/- 0.185 s | 9.780 +/- 0.064 s | 2.63x | 26.246 +/- 0.216 s | 10.108 +/- 0.081 s | 2.60x |
| `view_sites_bcf` | 2 | 27.238 +/- 0.124 s | 11.338 +/- 0.075 s | 2.40x | 25.192 +/- 0.111 s | 9.076 +/- 0.063 s | 2.78x |
| `view_sites_bcf` | 4 | 27.260 +/- 0.101 s | 11.276 +/- 0.018 s | 2.42x | 25.284 +/- 0.143 s | 9.008 +/- 0.020 s | 2.81x |
| `query_sites` | 0 | 2.464 +/- 0.007 s | 2.460 +/- 0.003 s | 1.00x | 2.554 +/- 0.007 s | 2.544 +/- 0.004 s | 1.00x |
| `query_gt_keep2` | 0 | 5.708 +/- 0.053 s | 2.856 +/- 0.014 s | 2.00x | 5.836 +/- 0.061 s | 2.946 +/- 0.027 s | 1.98x |
| `query_gt_all` | 0 | 57.454 +/- 0.196 s | 41.258 +/- 0.114 s | 1.39x | 58.538 +/- 0.298 s | 42.248 +/- 0.267 s | 1.39x |
| `stats` | 0 | 2.770 +/- 0.011 s | 2.754 +/- 0.004 s | 1.01x | 2.836 +/- 0.016 s | 2.808 +/- 0.005 s | 1.01x |
| `filter_gt_keep2` | 0 | 47.540 +/- 0.285 s | 31.272 +/- 0.120 s | 1.52x | 48.006 +/- 0.319 s | 31.612 +/- 0.135 s | 1.52x |
| `filter_gt_keep2` | 2 | 49.904 +/- 0.172 s | 33.490 +/- 0.148 s | 1.49x | 47.874 +/- 0.206 s | 31.302 +/- 0.176 s | 1.53x |
| `filter_gt_keep2` | 4 | 49.658 +/- 0.145 s | 33.718 +/- 0.074 s | 1.47x | 47.638 +/- 0.170 s | 31.642 +/- 0.114 s | 1.51x |

Raw files in this directory:

- `timings.tsv`: parsed `/usr/bin/time -p` output for every repetition.
- `checks.tsv`: baseline-vs-plan checksum comparison status.
- `checksums.tsv`: streamed output checksums and byte counts.
- `commands.tsv`: command descriptions.
- `metadata.tsv`: input, binary, sample, repeat, thread, and resume metadata.
