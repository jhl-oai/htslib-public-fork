# BAM Throughput: Current State

This document records the current implementation and benchmark setup for the
BAM throughput workstream. Update it whenever the product branch changes or the
benchmark interpretation changes.

## Branches

```text
feature/bam-throughput-product   product implementation branch
feature/bam-throughput-corpus    benchmark/evidence branch
```

Both branches currently start at HTSlib `develop` commit `2b3ddfd0`.

## Current Product Implementation

The product worktree currently has an HTSlib-internal BAM batch/record-view
reader plus a Samtools `view` consumer used to validate wholesale migration
away from the BAM `sam_read1()` loop. No public HTSlib structs or exported
function signatures are changed.

Current product worktree:

```text
worktrees/htslib-bam-product/
```

Core implementation files:

```text
worktrees/htslib-bam-product/sam.c
worktrees/htslib-bam-product/sam_internal.h
worktrees/htslib-bam-product/test/sam.c
samtools/bedidx.c
samtools/bedidx.h
samtools/region_plan.c
samtools/region_plan.h
samtools/bam2depth.c
samtools/sam_view.c
samtools/test/test.pl
```

Implementation summary:

- adds internal `sam_bam_read_batch()` and `sam_bam_itr_next_batch()` entry
  points for BAM record-view batches;
- adds internal `sam_bam_read_batch_count()` and `sam_bam_read_core_batch()`
  for exact count-only scans and MAPQ/flag count filters without full
  `bam1_t` materialization;
- reuses the internal core-batch buffer across reads to avoid per-batch
  `bam1_core_t` allocation churn;
- adds internal `sam_bam_batch_record_write1()` for unchanged BAM output
  records, so no-transform BAM output can write the raw input frame from the
  batch record when output indexing and mutations are absent;
- adds internal `sam_bam_batch_write1()` and segmented batch storage so a
  batch can describe multiple decoded BGZF-backed byte ranges while preserving
  record order and batch lifetime;
- adds internal record-view accessors for core fields, qname, CIGAR, seq/qual,
  aux lookup, aux length, virtual offsets, CG-aware query length, and end
  position;
- adds `sam_bam_prepare_batch_reader()` so Samtools can request the batch
  reader before applying `-@`, preserving thread-pool intent without routing
  input through the ordinary BGZF threaded reader;
- supports streaming BAM, indexed regions, multi-region/BED, `--fetch-pairs`,
  `-U`, count/filter paths, SAM/BAM/CRAM output materialization, tag edits,
  sanitize/unmap, filter expressions, library filtering, and remove-B;
- materializes to `bam1_t` only where existing Samtools semantics require
  mutation, complex filter evaluation, or output formatting;
- keeps SAM and CRAM input on their existing reader paths.
- keeps an ordered fused decode+parse prototype behind
  `HTS_BAM_BATCH_FUSED=1`; it is correctness-clean in targeted checks but
  serialized at fused-job boundaries and therefore not a default candidate.
- adds private HTSlib region helpers that coalesce copied `hts_reglist_t`
  query intervals and apply a caller-supplied exact overlap predicate during
  serial or BAM batch iterator reads;
- migrates the messy-BED `view -L` exact-filter path to those HTSlib internal
  wrappers while keeping BED parsing and command eligibility in Samtools.

Current constraints:

- preserve public HTSlib headers and exported symbols;
- keep changes at the HTSlib BAM/BGZF/SAM layer when possible;
- add small deterministic tests only when they directly cover changed behavior;
- keep benchmark corpus files and exploratory notes out of the product branch.

## Current Region-Based Speedup Status

The messy-BED speedup is now expressed as HTSlib-private iterator plumbing:

- Samtools still owns BED parsing, `bed_overlap()`, special `*`/`.` reference
  checks, local-index eligibility, and fallback policy.
- HTSlib internally owns the generic pieces: copied interval coalescing and
  exact-overlap filtered iterator reads.
- `samtools view -b/-u -L` uses the shared HTSlib wrappers for the local
  indexed coordinate-sorted BAM fast path.
- `samtools merge -L` remains on the normalized region-list helper because its
  semantics are already region-list based.
- Empty BEDs and special-reference BEDs (`*`/`.`) stay on the existing
  streaming/fallback path rather than entering the auto-index shortcut.
- Count-only `view -c -L` now uses the same exact filtered iterator when
  auto-indexing is eligible.
- Default `samtools depth -b` now uses the same shared region-plan machinery
  to build one merged indexed query per input when the BED is non-empty,
  non-special, `-r` is absent, `-a/-aa` is absent, and every BAM/CRAM input has
  an index.  The existing per-position `bed_overlap()` predicate still controls
  output rows.
- `depth -b` falls back to streaming when a BED names references absent from any
  input header, preserving legacy quiet-stderr behavior.

Fresh single warm-cache messy-BED timings after the HTSlib-internal migration
are directional only.  They have not yet been promoted into the repeated median
corpus harness output.

| Input | Command shape | Stock | Sambamba | Product |
| --- | --- | ---: | ---: | ---: |
| dense highcov | `view -@8 -b -L ... -o /dev/null` | 1.90s | 2.04s | 0.87s |
| dense highcov | `view -@8 -u -L ... -o /dev/null` | 1.24s | 0.27s | 0.26s |
| exome | `view -@8 -b -L ... -o /dev/null` | 0.16s | 0.10s | 0.06s |
| exome | `view -@8 -u -L ... -o /dev/null` | 0.14s | 0.03s | 0.03s |
| ONT | `view -@8 -b -L ... -o /dev/null` | 1.36s | 1.65s | 0.68s |
| ONT | `view -@8 -u -L ... -o /dev/null` | 0.81s | 0.22s | 0.18s |

Validation:

- HTSlib `./test/sam` passed, including coalesced-region and filtered-iterator
  checks.
- HTSlib `make -j8` passed.
- Samtools `make -j8` passed.
- Samtools `env PATH=.../worktrees/htslib-bam-product:$PATH ./test/test.pl view merge`
  passed: 1141 total, 1091 passed, 50 expected failures, 0 unexpected failures.
- Dense, exome, and ONT messy-BED decoded SAM outputs matched installed
  production Samtools exactly.
- After the `depth -b` migration, full Samtools
  `env PATH=.../worktrees/htslib-bam-product:$PATH ./test/test.pl` passed:
  1153 total, 1103 passed, 50 expected failures, 0 unexpected failures.
- Dense, exome, and ONT `depth -@8 -b messy` outputs matched installed
  production Samtools byte-for-byte.

Repeated-median `depth -b` timings are now decision-quality enough to promote
from directional notes.  Command shape was
`depth -@ THREADS -b messy.bed -o /dev/null input.bam`, five repeats per row.

| Threads | Input | Stock | Product | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 0 | dense highcov | 7.77s | 0.99s | 7.85x |
| 0 | exome | 0.81s | 0.11s | 7.36x |
| 0 | ONT | 6.97s | 0.92s | 7.58x |
| 8 | dense highcov | 1.36s | 0.32s | 4.25x |
| 8 | exome | 0.28s | 0.07s | 4.00x |
| 8 | ONT | 1.88s | 0.37s | 5.08x |

The detailed command log and Sambamba context are in
`data/benchmarks/depth-bed-region-plan-20260505.md`.

Next candidates to investigate: BED-driven `bedcov` and `coverage`, where the
same coarse merged query plus exact per-position/per-base predicate should be
possible if each command's zero/summary semantics are preserved.

## Current Corpus Implementation

The corpus branch has scaffold and HTSlib-specific product-vs-base timing
support:

```text
bench/bam-shape/
  inputs.tsv
  README.md
  scripts/run_bam_bench.sh
  scripts/run_htslib_test_view_bench.pl
```

The Samtools/Sambamba harness supports:

- Samtools-compatible binaries through `SAMTOOLS=/path/to/samtools`;
- Sambamba through `SAMBAMBA=/path/to/sambamba`;
- repeated runs through `REPEATS=N`;
- thread sweeps through `THREADS_LIST="0 2 4 8"`;
- result output under `bench/bam-shape/results/`.

The generated result directory is ignored by git. Copy only curated summaries
into tracked docs or tracked result-summary directories when they are worth
keeping.

The HTSlib `test_view` harness supports:

- explicit base and product binaries through `BASE_TEST_VIEW` and
  `PRODUCT_TEST_VIEW`;
- decode-only mode with `test_view -B`;
- decode+recompress mode with `test_view -b -p /dev/null`;
- repeated in-row loops through `LOOPS=N` for second-scale timings;
- metadata capture for paths, git commits/status, command lines, libdeflate
  version, checks, and timing status.

The BGZF qsize matrix harness additionally records macOS maximum resident set
size via `/usr/bin/time -l`:

```text
bench/bam-shape/scripts/run_bgzf_qsize_matrix.pl
```

## Current Inputs

The active local inputs are listed in `bench/bam-shape/inputs.tsv` and point to
the umbrella workspace's `data/` directory.

| Input | Shape |
| --- | --- |
| `HG00096.lowcov.chr20_1mb.bam` | short-read low-coverage smoke slice |
| `HG00096.lowcov.chr20_10-20Mb.bam` | short-read low-coverage moderate slice |
| `HG00096.highcov.chr20_10-11Mb.bam` | short-read high-coverage dense slice |
| `HG00096.exome.chr20.bam` | exome/sparse indexed-access slice |
| `HG002.ont_ul.chr20_10-10.2Mb.bam` | Oxford Nanopore ultra-long-read slice |

Source URLs, regions, sizes, and retrieval commands are recorded in
`data/README.md` at the umbrella workspace root.

## Current Smoke Result

The initial scaffold smoke run verified that Samtools and Sambamba count the
same records for all five local inputs with `THREADS_LIST=0`.

The smoke harness wrote ignored local files:

```text
bench/bam-shape/results/metadata.tsv
bench/bam-shape/results/checks.tsv
bench/bam-shape/results/timings.tsv
```

This smoke run is not a decision-quality benchmark. It was used to validate
paths, count checks, and Sambamba banner handling.

For the first product prototype, base and product HTSlib `test_view` produced
byte-identical BAM output on the high-coverage chr20 1 Mb slice at `-@ 4`, and
both outputs passed `samtools quickcheck` with 216,942 records.

The product worktree also passes `make test`: 358 HTSlib tests passed, 0
failed, and the libcurl retry tests passed.

## Current Benchmark Summary

### Wholesale `samtools view` BAM Batch Reader

Correctness status:

- HTSlib `make test` passed outside the sandbox: 361 passed, 0 failed.
- Samtools `make test` passed when linked against
  `worktrees/htslib-bam-product/libhts.a` and compiled with
  `-DHTS_BAM_BATCH_READER_CONSUMER`.
- Manual parity checks covered dense streaming counts, filtered counts,
  indexed regions, multi-region/BED, `--fetch-pairs`, `-U`, and SAM output
  byte comparisons.
- Post-review correctness fixes cover sanitize-before-filter parity,
  missing-NUL QNAME repair for QNAME-dependent filters, and iterator
  CIGAR/query-length validation before out-of-region records are skipped.

Current large-slice timings, stock installed Samtools vs migrated local build:

| Input | Command shape | Stock | Migrated | Speedup |
| --- | --- | ---: | ---: | ---: |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -c` | 1.02s | 0.92s | 1.11x |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -q 20 -c` | 1.05s | 0.92s | 1.14x |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -O SAM -o /dev/null` | 2.60s | 2.42s | 1.07x |
| `HG002.ont_ul.chr20_10-40Mb.bam` | `view -@8 -c` | 0.78s | 0.70s | 1.11x |

Interpretation: the feature coverage is now broad enough for `samtools view`
BAM migration testing, but the current design is not yet a Sambamba-scale
speedup under `-@`. The next high-ROI work is a correct ordered parse-stage
scheduler that parallelizes record-view construction/validation after BGZF
inflate. A first attempt at this was rejected because it dropped about one
queue-depth of records from streaming input.

Latest continuation status:

- exact count and MAPQ/flag count have internal HTSlib count/core helpers;
- raw unchanged BAM output has an internal frame writer and passed a decoded
  SAM byte-compare on the 1 Mb highcov slice;
- HTSlib `make test` passed after the latest changes: 361 passed, 0 failed;
- Samtools `make test` passed against product HTSlib: 1,068 total, 0 failed,
  46 expected failures;
- CIGAR-length caching and multi-block aggregation attempts were removed after
  they regressed local timings;
- fresh large-slice timings were taken under heavy host CPU load and are not
  decision-quality, but they still show the target is unmet.

Representative continuation timings:

| Input | Command shape | Stock | Product | Interpretation |
| --- | --- | ---: | ---: | --- |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -c` | 1.46-1.69s | 1.54-1.62s | Neutral within noise before host saturation. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -q 20 -c` | 1.44-1.50s | 1.41-1.50s | Neutral within noise before host saturation. |
| `HG002.ont_ul.chr20_10-40Mb.bam` | `view -@8 -c` | 1.18-1.48s | 0.94-1.01s | Modest ONT win in this run. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -b -o /dev/null` | 15.25s | 16.84s | Raw frame writer did not beat stock with output threads. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@1 -b -o /dev/null` | 72.46s | 66.06s | Slight no-thread output win, not enough. |

Current interpretation: useful internal plumbing is in place, but the design is
still not sufficiently general or fast under threaded `-@` workloads.  The next
endpoint should be a fused ordered HTSlib pipeline that shares one scheduler
for BGZF inflate, frame scan, optional core/view construction, and raw BAM
frame output.

Latest fused/segmented checkpoint:

- `HTS_BAM_BATCH_FUSED=1` now uses segmented batches for multi-block decoded
  storage and keeps cross-job carry correct by allowing only one fused job in
  flight.
- Targeted checks pass for dense streaming, split/large records, indexed
  regions, and fused BAM output decoded back to SAM.
- Final full suites pass after this checkpoint: HTSlib `make test` reports
  361 passed, 0 failed; Samtools `make test` reports 1,068 total, 0 failed,
  46 expected failures, with regression suites passing.
- Product default, not fused, carries the useful current win:

| Input | Command shape | Stock | Product default | Fused segmented |
| --- | --- | ---: | ---: | ---: |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -b -o /dev/null` | 14.79s | 9.63s | 12.83s |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -O SAM -o /dev/null` | 4.57s | 3.28s | 11.44s |
| `HG002.ont_ul.chr20_10-40Mb.bam` | `view -@8 -b -o /dev/null` | 11.13-11.45s | 10.93-11.50s | not kept |

Current decision: do not default the fused path.  The next product endpoint is
safe multi-in-flight fused scheduling: preserve ordered emission, but either
prove a physical job boundary ends on a record boundary or replay prefetched
jobs when the prior job returns carry.

### Threaded BGZF Libdeflate Object Reuse

Focused result files:

```text
data/benchmarks/libdeflate-cache-focused-readwrite/
data/benchmarks/libdeflate-cache-focused-decode/
```

Setup:

- input: `exome_chr20`, region `20`, 1,933,817 records;
- base/product HTSlib commit: `2b3ddfd0`;
- product status: `M bgzf.c`;
- libdeflate: `1.25`;
- read/write: `REPEATS=5`, `LOOPS=3`, `THREADS_LIST="2 4"`;
- decode-only: `REPEATS=5`, `LOOPS=12`, `THREADS_LIST="2 4"`.

Median seconds per `test_view` invocation:

| Mode | Threads | Base | Product | Speedup |
| --- | ---: | ---: | ---: | ---: |
| decode | 2 | 0.323798 | 0.322885 | 0.28% |
| decode | 4 | 0.175482 | 0.173781 | 0.97% |
| readwrite | 2 | 2.479780 | 2.469370 | 0.42% |
| readwrite | 4 | 1.292620 | 1.280490 | 0.94% |

Interpretation: correct and low-risk, but only marginally positive on the
largest current local workload. Treat this as a possible cleanup-level
optimization, not a high-ROI throughput change unless larger or more
compressor-heavy workloads show a stronger effect.

### Opt-In Queue Size

Rejected/reverted experiment: changing the default BGZF queue depth from
`2 * threads` to `4 * threads`. It showed a promising 5.46% decode-only result
at `-@ 8` on `exome_chr20`, but reviewers flagged the global memory/behavior
cost and the benchmark did not isolate queue depth from the libdeflate-cache
prototype. The product branch is back to the stock `qsize=0` default.

Current opt-in queue experiment: `HTS_BGZF_QUEUE_SIZE` can override the BGZF
queue size only when callers pass `qsize=0`. The documented default remains
`2 * threads` when the variable is absent or invalid. In an isolated same-binary
run at `-@ 8` on `exome_chr20`, decode medians improved from 0.096737 seconds
to 0.092395 seconds at q32 and 0.090937 seconds at q64. Read/write remained
neutral.

The broader RSS matrix is in:

```text
data/benchmarks/bgzf-qsize-matrix-20260501/
```

Key readout: q32/q64 help the largest 8-thread exome decode case by 5.03% and
6.87%, respectively, at about +3.1 MiB and +9.3 MiB max RSS. Other decode
workloads are mixed, and read/write remains mostly neutral except for smaller
tail-latency-sensitive runs. Keep this as an advanced opt-in knob, not a
default.

### Threaded Writer Flush Condition Variable

Focused result files:

```text
data/benchmarks/bgzf-cond-flush-focused-readwrite-final/
```

Setup:

- inputs: all five local corpus entries;
- mode: `test_view -b -p /dev/null`;
- threads: 4 and 8;
- `REPEATS=5`, `LOOPS=3`;
- product stack includes libdeflate cache, direct threaded BGZF read, qsize env
  support, and the condition-variable flush change.

Median seconds per `test_view` invocation:

| Threads | Input | Base | Product | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 4 | `exome_chr20` | 1.321340 | 1.253463 | 5.14% |
| 4 | `highcov_chr20_10_11mb` | 0.512260 | 0.480912 | 6.12% |
| 4 | `lowcov_chr20_10_20mb` | 0.379690 | 0.351527 | 7.42% |
| 4 | `lowcov_chr20_1mb` | 0.071738 | 0.045310 | 36.84% |
| 4 | `ont_ul_chr20_10_10_2mb` | 0.111300 | 0.088792 | 20.22% |
| 8 | `exome_chr20` | 0.681446 | 0.656175 | 3.71% |
| 8 | `highcov_chr20_10_11mb` | 0.278821 | 0.252144 | 9.57% |
| 8 | `lowcov_chr20_10_20mb` | 0.209740 | 0.182391 | 13.04% |
| 8 | `lowcov_chr20_1mb` | 0.052884 | 0.026395 | 50.09% |
| 8 | `ont_ul_chr20_10_10_2mb` | 0.072370 | 0.049106 | 32.15% |

Interpretation: this is the first clear high-ROI HTSlib-level win in the BAM
workstream. The largest absolute workloads improve by 3.7-13.0%, while smaller
outputs see larger relative gains because the old `hts_usleep(10000)` polling
added visible flush/close tail latency.

### BAM Write Packing Fast Path

Focused result files:

```text
data/benchmarks/bam-write-pack-focused-readwrite-final/
```

This is an incremental benchmark against the saved product binary before the
BAM-local fast path, so it isolates the `sam.c` change from the larger BGZF
threading work.

Median seconds per `test_view` invocation:

| Threads | Input | Before | After | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 4 | `exome_chr20` | 1.256242 | 1.254703 | 0.12% |
| 4 | `highcov_chr20_10_11mb` | 0.472840 | 0.471703 | 0.24% |
| 4 | `lowcov_chr20_10_20mb` | 0.341788 | 0.343143 | -0.40% |
| 4 | `lowcov_chr20_1mb` | 0.047099 | 0.045361 | 3.69% |
| 4 | `ont_ul_chr20_10_10_2mb` | 0.089909 | 0.088485 | 1.58% |
| 8 | `exome_chr20` | 0.651759 | 0.654103 | -0.36% |
| 8 | `highcov_chr20_10_11mb` | 0.250566 | 0.249691 | 0.35% |
| 8 | `lowcov_chr20_10_20mb` | 0.184770 | 0.182363 | 1.30% |
| 8 | `lowcov_chr20_1mb` | 0.027072 | 0.026187 | 3.27% |
| 8 | `ont_ul_chr20_10_10_2mb` | 0.051217 | 0.049555 | 3.25% |

Interpretation: correct and low blast radius, but marginal on larger current
workloads. Keep it as an optional small cleanup only if we are comfortable with
the extra test surface; it is not a headline Sambamba-derived win.

After the qsize docs and direction-specific cache allocation cleanup, targeted
`make check TEST_OPTS='test_view sam'` passed: 361 tests passed, 0 failed, and
the libcurl retry tests passed.

After the direct-read and condition-variable changes:

- `make lib-static test/test_view htsfile test/test_bgzf` passed;
- `./test/test_bgzf test/bgziptest.txt` passed;
- `git diff --check` passed;
- full `make test` passed: 361 tests passed, 0 failed, and the libcurl retry
  tests passed;
- highcov `test_view -@ 8 -b` output remained byte-identical to the base
  output and passed `samtools quickcheck`;
- the same highcov write smoke with `HTS_BGZF_QUEUE_SIZE=1` also remained
  byte-identical and passed `samtools quickcheck`.

## Current BAM view migration status, 2026-05-03 follow-up

The product and Samtools worktrees now represent a no-fallback BAM `view`
migration for this branch:

- `samtools/sam_view.c` uses `sam_bam_prepare_batch_reader()`,
  `sam_bam_read_batch()`, and `sam_bam_itr_next_batch()` unconditionally for
  BAM input.  The previous `HTS_BAM_BATCH_READER_CONSUMER` compile-time guard
  has been removed.
- `samtools/Makefile` tracks `$(HTSDIR)/sam_internal.h` unconditionally, so
  clean no-macro builds have the right dependency graph.
- SAM and CRAM input paths are unchanged.
- BAM batch-reader unavailability is a hard error, not a feature fallback to
  the old BAM `sam_read1()` loop.
- Plain streaming batch reads do not compute per-record virtual offsets; region
  and multi-region iterator reads still do.
- An ordered async record-view parser exists behind `HTS_BAM_BATCH_PARSE=1`,
  but remains experimental because it is slower on the current dense count
  benchmark despite passing targeted parity checks.

Validation:

- product HTSlib: `make libhts.a test/sam` and `make test` passed; full suite
  result was 361 passed, 0 failed;
- Samtools linked against product HTSlib without any batch-consumer CPP flag:
  full `make test` passed with 1,068 tests, 0 failed, 46 expected failures, and
  the regression suites reported `PASS`;
- manual parity covered highcov count (`7,579,904`), highcov `-q 20 -c`
  (`7,396,382`), ONT count (`156,695`), highcov indexed region count
  (`216,942`), async parser count/region parity with `HTS_BAM_BATCH_PARSE=1`,
  and a `--no-PG -O SAM` byte comparison on a highcov region.

Latest single-run timing remains noisy on this machine.  After the no-fallback
cleanup, Makefile dependency fix, and configured no-macro rebuild, a fresh
paired run had installed stock `view -@8 -c` on the highcov 40 Mb slice at
1.18s and product batch at 1.16s; ONT 30 Mb was 0.87s stock vs 0.79s product.
Earlier paired runs showed up to 1.29s stock vs 1.11s product on the same
highcov command, so current evidence should be read as modest/noisy rather
than a headline speedup.  The async record-view parser was slower at
2.40-2.58s and is not the default.

## Next Product Candidates

Candidate work should be evaluated in roughly this order:

1. Establish stable baseline timings across Samtools and Sambamba on the local
   corpus.
2. Identify which workloads are dominated by BAM decode, BGZF inflate/deflate,
   allocation churn, iterator overhead, or output writing.
3. Prototype the smallest HTSlib-internal change in the product worktree.
4. Re-run the same corpus against stock and product-linked builds.
5. Promote only changes with consistent wins and no correctness drift.

## Current Follow-up: Selective Compressed Filtered BAM Output

The current product branch additionally improves selective compressed
`samtools view -b FILTERS -o out.bam` output through small internal changes:

- core-field-only filters avoid unnecessary batch record-view construction;
- `sam_bam_batch_write1_ranges()` resolves BGZF output state once per range
  write and can pack sparse selected raw ranges into a temporary batch-local
  buffer before writing;
- public HTSlib/Samtools APIs remain unchanged.  This experimental Samtools
  branch still intentionally builds against the paired HTSlib source tree and
  private `sam_internal.h`.

Correctness is clean on the current tree: product HTSlib `./test/sam` passed;
Samtools `make check` passed with 1155 total tests, 1105 passed, 50 expected
failures, and all regression scripts PASS; decoded SAM equality was checked for
dense, exome, and ONT real-output BAM cases.

Performance status is split by output selectivity.  Dense high-pass compressed
BAM output is still BGZF-compression-bound: dense `-@8 -q20 -F4` was 8.54s
stock vs 8.06s product, and ONT `-@8 -q20 -F4` was 4.67s vs 4.52s.  Selective
compressed output is the useful win: dense `-@8 -e 'mapq < 20'` was 1.49s
stock vs 1.01s product, and exome `-@0 -e 'mapq < 20'` was 0.79s vs 0.65s.

Candidate changes tried and removed:

- making `HTS_BAM_BATCH_PARSE=1` the default for threaded batch views was not a
  stable win, so async record-view parsing remains opt-in;
- reusing stock BGZF MT from the batch reader was correct in small-slice parity
  tests but slower in both copy-based and retained-block prototypes;
- enabling raw ranges adaptively for compressed pure flag-only filters regressed
  dense `-f 4`, so the pure-flag compressed path stays conservative.

Current conclusion: the batch migration is broad and correctness-clean.
Selective filtered BAM output is meaningfully faster, but dense high-pass
compressed BAM rewriting should be described as neutral to slightly faster, not
as a headline compressed-output speedup.
