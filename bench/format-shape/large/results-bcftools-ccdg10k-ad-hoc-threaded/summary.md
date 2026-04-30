# bcftools CCDG 10k Threaded Ad Hoc Benchmark

This run measures realistic bcftools operations on the CCDG chromosome 22
10,000-record slice.  The input is a BGZF-compressed VCF text file with 3,202
samples and 123,431,420 compressed bytes.

The benchmark compares the generic parser (`HTS_VCF_FORMAT_PLAN=0`) against the
planned FORMAT parser (`HTS_VCF_FORMAT_PLAN=1`) using five repetitions.  Output
streams were piped through `cksum`; all 80 baseline-vs-plan checksum comparisons
were `ok`.

`--threads` was applied only to `bcftools view` command shapes.  `bcftools query`
and `bcftools stats` rows are unthreaded and therefore use `threads=0`.

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
INPUT=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/htslib-vcf-avx-sanity/bench/format-shape/public/ccdg_chr22_10k.vcf.gz \
NAME=ccdg_10k \
OUTDIR=bench/format-shape/large/results-bcftools-ccdg10k-ad-hoc-threaded \
REPS=5 \
THREADS_LIST="0 2 4" \
  bash bench/format-shape/scripts/run_bcftools_ad_hoc_bench.sh
```

The bcftools binary reported `bcftools 1.23.1-63-gd9375776` using
`htslib 1.23.1-56-g8c809baf`.  Times are mean +/- standard error across five
sequential runs.

| Command | Threads | Baseline user | Planned user | User speedup | Baseline real | Planned real | Real speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 0 | 2.398 +/- 0.007 s | 2.200 +/- 0.005 s | 1.09x | 3.914 +/- 0.011 s | 3.718 +/- 0.009 s | 1.05x |
| `view_bcf` | 2 | 2.908 +/- 0.004 s | 2.712 +/- 0.002 s | 1.07x | 3.116 +/- 0.005 s | 2.932 +/- 0.005 s | 1.06x |
| `view_bcf` | 4 | 2.918 +/- 0.002 s | 2.720 +/- 0.000 s | 1.07x | 3.110 +/- 0.000 s | 2.932 +/- 0.002 s | 1.06x |
| `view_keep2_bcf` | 0 | 2.516 +/- 0.004 s | 2.314 +/- 0.002 s | 1.09x | 2.574 +/- 0.004 s | 2.366 +/- 0.002 s | 1.09x |
| `view_keep2_bcf` | 2 | 2.672 +/- 0.006 s | 2.474 +/- 0.002 s | 1.08x | 2.274 +/- 0.005 s | 2.074 +/- 0.002 s | 1.10x |
| `view_keep2_bcf` | 4 | 2.684 +/- 0.004 s | 2.476 +/- 0.004 s | 1.08x | 2.284 +/- 0.004 s | 2.074 +/- 0.002 s | 1.10x |
| `view_sites_bcf` | 0 | 2.512 +/- 0.007 s | 2.304 +/- 0.009 s | 1.09x | 2.564 +/- 0.012 s | 2.350 +/- 0.010 s | 1.09x |
| `view_sites_bcf` | 2 | 2.648 +/- 0.004 s | 2.448 +/- 0.002 s | 1.08x | 2.244 +/- 0.007 s | 2.044 +/- 0.002 s | 1.10x |
| `view_sites_bcf` | 4 | 2.650 +/- 0.006 s | 2.450 +/- 0.003 s | 1.08x | 2.244 +/- 0.005 s | 2.042 +/- 0.004 s | 1.10x |
| `query_sites` | 0 | 0.380 +/- 0.000 s | 0.382 +/- 0.002 s | 0.99x | 0.404 +/- 0.002 s | 0.410 +/- 0.000 s | 0.99x |
| `query_gt_keep2` | 0 | 0.640 +/- 0.000 s | 0.410 +/- 0.000 s | 1.56x | 0.670 +/- 0.000 s | 0.438 +/- 0.002 s | 1.53x |
| `query_gt_all` | 0 | 2.880 +/- 0.005 s | 2.664 +/- 0.004 s | 1.08x | 2.954 +/- 0.009 s | 2.728 +/- 0.006 s | 1.08x |
| `stats` | 0 | 0.420 +/- 0.000 s | 0.420 +/- 0.000 s | 1.00x | 0.450 +/- 0.000 s | 0.452 +/- 0.002 s | 1.00x |
| `filter_gt_keep2` | 0 | 2.790 +/- 0.003 s | 2.580 +/- 0.000 s | 1.08x | 2.846 +/- 0.004 s | 2.638 +/- 0.002 s | 1.08x |
| `filter_gt_keep2` | 2 | 2.940 +/- 0.003 s | 2.732 +/- 0.002 s | 1.08x | 2.544 +/- 0.006 s | 2.340 +/- 0.004 s | 1.09x |
| `filter_gt_keep2` | 4 | 2.944 +/- 0.006 s | 2.740 +/- 0.003 s | 1.07x | 2.546 +/- 0.009 s | 2.350 +/- 0.018 s | 1.08x |

Raw files in this directory:

- `timings.tsv`: parsed `/usr/bin/time -p` output for every repetition.
- `checks.tsv`: baseline-vs-plan checksum comparison status.
- `checksums.tsv`: streamed output checksums and byte counts.
- `commands.tsv`: command descriptions.
- `metadata.tsv`: input, binary, sample, repeat, and thread metadata.
