# BAM Throughput Overview

This branch is the corpus and evidence side of the HTSlib-level BAM throughput
workstream. The goal is to identify which Sambamba-inspired optimizations are
worth porting into HTSlib/Samtools without coupling the benchmark corpus to the
product implementation branch.

## Goal

Improve BAM read/write throughput in HTSlib so the gains apply broadly across
Samtools commands and other HTSlib consumers. The product branch should preserve
public API and ABI unless a deliberate exception is made later.

The first candidate areas are internal, low-friction changes:

- BGZF queue-depth and scheduling behavior;
- allocation/reallocation churn in BAM record reads;
- decode/encode batching where public APIs can remain unchanged;
- iterator and indexed-read overheads;
- compression/decompression backend usage, especially libdeflate behavior.

## Branch Model

The workstream intentionally uses two independent branches from `develop`:

- `feature/bam-throughput-product`: product implementation and small tests.
- `feature/bam-throughput-corpus`: benchmark corpus, scripts, result summaries,
  and running analysis docs.

The corpus branch must not depend on product-branch commits. It should compare
arbitrary builds by explicit binary paths, for example stock Samtools, a
product-worktree-linked Samtools, and Sambamba.

## Benchmark Corpus

The initial local corpus lives in the umbrella workspace's `data/` directory and
is documented in `data/README.md`. It includes:

- low-coverage short-read WGS;
- high-coverage short-read WGS;
- exome/sparse indexed access;
- Oxford Nanopore ultra-long reads.

The corpus branch tracks manifests and scripts under `bench/bam-shape/`, while
large BAM payloads stay out of git.

## Evidence Standard

Before an optimization graduates into the product branch, collect:

- repeated timings with medians, not one-off runs;
- read-count or checksum checks proving equal work;
- thread-scaling data where threading is relevant;
- build/version metadata for Samtools, HTSlib, Sambamba, compression libraries,
  and compiler flags;
- notes on neutral or negative cases, especially long-read and sparse indexed
  workloads.

The product branch should only carry changes that survive this evidence pass.
