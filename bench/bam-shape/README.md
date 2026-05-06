# BAM Shape Benchmark Corpus

This directory holds the corpus-side benchmark scaffolding for HTSlib-level BAM
throughput work inspired by Sambamba. It is intentionally separate from the
product branch so benchmark scripts and result artifacts do not have to travel
with the upstreamable implementation.

The running docs are:

- `docs/BAM_THROUGHPUT_OVERVIEW.md`
- `docs/BAM_THROUGHPUT_CURRENT.md`
- `docs/BAM_THROUGHPUT_EXPERIMENT_LOG.md`

## Layout

```text
bench/bam-shape/
  inputs.tsv                  manifest of local benchmark BAM slices
  scripts/run_bam_bench.sh    repeat-capable Samtools/Sambamba timing harness
  results/                    generated timing/check summaries, ignored by git
```

The BAM payloads are not stored in this branch. The initial local inputs live in
the umbrella workspace's `data/` directory and are documented there.

## Running

From this worktree root:

```sh
SAMTOOLS=./test/test_view \
SAMBAMBA=../../sambamba/bin/sambamba-1.0.1 \
REPEATS=5 \
THREADS_LIST="0 2 4 8" \
sh bench/bam-shape/scripts/run_bam_bench.sh
```

Useful variants:

```sh
# Benchmark the sibling Samtools checkout.
SAMTOOLS=../../samtools/samtools sh bench/bam-shape/scripts/run_bam_bench.sh

# Benchmark a product-worktree Samtools build, if one is linked against the
# product HTSlib worktree.
SAMTOOLS=/path/to/product-linked/samtools sh bench/bam-shape/scripts/run_bam_bench.sh
```

The script writes:

```text
bench/bam-shape/results/metadata.tsv
bench/bam-shape/results/timings.tsv
bench/bam-shape/results/checks.tsv
```

`timings.tsv` records one row per tool/input/thread/repeat. `checks.tsv`
records read counts so timing comparisons are not silently measuring different
work. Use medians from repeated runs for product decisions.

## Current Inputs

See `../../data/README.md` from the umbrella workspace for exact source URLs,
remote regions, retrieval commands, sizes, and read counts.
