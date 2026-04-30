# bcftools CCDG 10k Ad Hoc Benchmark

This run measures realistic bcftools operations on the CCDG chromosome 22
10,000-record slice.  The input is a BGZF-compressed VCF text file with 3,202
samples, 123,431,420 compressed bytes, derived from the 1000G/CCDG high-coverage
chromosome 22 VCF.

The benchmark compares the generic parser (`HTS_VCF_FORMAT_PLAN=0`) against the
planned FORMAT parser (`HTS_VCF_FORMAT_PLAN=1`) using:

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
INPUT=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/htslib-vcf-avx-sanity/bench/format-shape/public/ccdg_chr22_10k.vcf.gz \
OUTDIR=bench/format-shape/large/results-bcftools-ccdg10k-ad-hoc \
REPS=5 \
  bash bench/format-shape/scripts/run_bcftools_ad_hoc_bench.sh
```

The bcftools binary reported `bcftools 1.23.1-63-gd9375776` using
`htslib 1.23.1-56-g8c809baf`.  Output streams were piped through `cksum`;
all 40 baseline-vs-plan checksum comparisons were `ok`.

Times are mean +/- standard error across five sequential runs.

| Command | Baseline user | Planned user | User speedup | Baseline real | Planned real | Real speedup | Notes |
|---|---:|---:|---:|---:|---:|---:|---|
| `view_bcf` | 2.386 +/- 0.005 s | 2.194 +/- 0.007 s | 1.09x | 3.894 +/- 0.010 s | 3.702 +/- 0.010 s | 1.05x | Full `bcftools view -Ob -l 0` conversion for all samples. |
| `view_keep2_bcf` | 2.440 +/- 0.003 s | 2.248 +/- 0.002 s | 1.09x | 2.484 +/- 0.004 s | 2.288 +/- 0.005 s | 1.09x | `bcftools view` retaining the first two samples. |
| `view_sites_bcf` | 2.434 +/- 0.014 s | 2.234 +/- 0.005 s | 1.09x | 2.480 +/- 0.018 s | 2.272 +/- 0.005 s | 1.09x | `bcftools view -G`, dropping genotype/sample columns. |
| `query_sites` | 0.370 +/- 0.000 s | 0.370 +/- 0.000 s | 1.00x | 0.398 +/- 0.002 s | 0.390 +/- 0.000 s | 1.02x | Fixed site-field query, mostly non-FORMAT parsing. |
| `query_gt_keep2` | 0.628 +/- 0.004 s | 0.400 +/- 0.000 s | 1.57x | 0.654 +/- 0.004 s | 0.424 +/- 0.002 s | 1.54x | `%GT` query for the first two samples. |
| `query_gt_all` | 2.832 +/- 0.007 s | 2.632 +/- 0.005 s | 1.08x | 2.888 +/- 0.010 s | 2.686 +/- 0.007 s | 1.08x | `%GT` query across all 3,202 samples. |
| `stats` | 0.410 +/- 0.000 s | 0.414 +/- 0.002 s | 0.99x | 0.440 +/- 0.000 s | 0.440 +/- 0.000 s | 1.00x | `bcftools stats`, effectively parity on this slice. |
| `filter_gt_keep2` | 2.762 +/- 0.007 s | 2.554 +/- 0.004 s | 1.08x | 2.824 +/- 0.013 s | 2.608 +/- 0.005 s | 1.08x | `bcftools view -s first-two -i 'GT="alt"' -Ob -l 0`. |

Raw files in this directory:

- `timings.tsv`: parsed `/usr/bin/time -p` output for every repetition.
- `checks.tsv`: baseline-vs-plan checksum comparison status.
- `checksums.tsv`: streamed output checksums and byte counts.
- `commands.tsv`: command descriptions.
- `metadata.tsv`: input, binary, sample, and repeat metadata.
