# BAM Throughput Experiment Log

This log records experiments, results, and decisions for the BAM throughput
workstream. Keep it chronological. Prefer short entries with enough detail to
reproduce or discard the result later.

## 2026-05-01: Workstream Split

Created two independent HTSlib worktrees from `develop`:

```text
worktrees/htslib-bam-product   feature/bam-throughput-product
worktrees/htslib-bam-corpus    feature/bam-throughput-corpus
```

Decision: keep the corpus branch independent from the product branch. The
corpus harness benchmarks explicit tool paths instead of inheriting product
commits.

Reasoning: this mirrors the FORMAT planner packaging model while making it
easier to compare stock HTSlib/Samtools, product prototypes, and Sambamba
without rebasing benchmark infrastructure through the implementation branch.

## 2026-05-01: Initial BAM Corpus

Downloaded bounded public BAM slices into the umbrella `data/` directory:

- HG00096 low-coverage WGS, chr20 1 Mb smoke interval;
- HG00096 low-coverage WGS, chr20 10 Mb interval;
- HG00096 high-coverage WGS, chr20 1 Mb interval;
- HG00096 exome, chromosome 20;
- HG002 Oxford Nanopore ultra-long reads, chr20 200 kb interval.

Decision: keep payloads out of git and track provenance in `data/README.md`.

Reasoning: the local files are large enough to produce useful timing signal but
small enough for repeated local runs. The long-read slice exercises a separate
record-shape regime from short-read WGS/exome BAMs.

## 2026-05-01: Initial Harness Scaffold

Added `bench/bam-shape/scripts/run_bam_bench.sh`.

Current behavior:

- reads `bench/bam-shape/inputs.tsv`;
- runs count checks for Samtools and Sambamba;
- runs repeated BAM-to-`/dev/null` streaming timings;
- supports `REPEATS`, `THREADS_LIST`, `SAMTOOLS`, `SAMBAMBA`, and `OUTDIR`;
- writes metadata, checks, and timing TSVs under an ignored `results/`
  directory.

Smoke result: Samtools and Sambamba record counts matched across all five local
inputs with `THREADS_LIST=0`.

Notes:

- Sambamba is available locally as `../../sambamba/bin/sambamba-1.0.1`.
- Sambamba prints version/banner text on normal runs, so the harness filters
  count output and redirects timing-run banner output.
- The smoke timings from the first scaffold run are not decision-quality
  benchmark evidence.

## Open Questions

- Which local command mix best isolates HTSlib BAM decode from Samtools command
  overhead?
- Should the first product prototype target allocation churn in `bam_read1()`,
  BGZF queue depth, or batch decode?
- What is the right product-linked Samtools build layout for comparing stock and
  product HTSlib without changing the corpus branch?
- Do the long-read and exome slices reveal different bottlenecks from dense
  short-read WGS?

## 2026-05-03: Wholesale BAM `samtools view` Batch-Reader Migration

Implemented the first full BAM `view` consumer on top of HTSlib internal
record-view batches.  The intent is no feature-based fallback to the old BAM
`sam_read1()` loop for `samtools view`: BAM streaming and indexed-region input
go through the HTSlib batch reader, while SAM and CRAM keep their existing
paths.

Implementation summary:

- HTSlib product branch adds internal batch APIs in `sam_internal.h` and
  `sam.c`: streaming `sam_bam_read_batch()`, indexed
  `sam_bam_itr_next_batch()`, record-view accessors, materialization helpers,
  CG-aware query/end-position helpers, virtual-offset tracking, seek/reset
  support, and `sam_bam_prepare_batch_reader()` so `-@` config is routed to the
  batch reader instead of ordinary BGZF threading.
- Samtools `sam_view.c` now compiles a BAM batch-reader consumer under
  `HTS_BAM_BATCH_READER_CONSUMER`.  It handles count/filter, indexed regions,
  multi-region/BED, `--fetch-pairs`, `-U`, SAM/BAM/CRAM output materialization,
  tag operations, sanitize/unmap, and expression/library/remove-B cases by
  materializing only where existing semantics require `bam1_t`.
- A count-only batch loop avoids selected/unselected output machinery when
  `view -c` has no write-side effects, while preserving raw filters and
  CIGAR/query validation.

Correctness checks:

- HTSlib `test/sam` targeted batch-reader tests pass, including invalid-TID
  terminal error handling for both unthreaded and threaded inputs.
- Post-review fixes preserve sanitize-before-filter behavior, materialize
  malformed missing-NUL QNAME records before QNAME-dependent filters, and make
  indexed batch end-position checks validate CIGAR/query-length parity before
  records can be skipped as out-of-region.
- Full HTSlib `make test` passed outside the sandbox: 361 passed, 0 failed.
  The sandboxed run failed only because `ref-cache` could not bind local ports.
- Full Samtools `make test` passed with `samtools` linked to
  `worktrees/htslib-bam-product/libhts.a` and compiled with
  `-DHTS_BAM_BATCH_READER_CONSUMER`.
- Manual parity checks covered dense streaming counts, filtered counts,
  indexed regions, multi-region/BED, `--fetch-pairs`, `-U`, and SAM output
  byte comparisons with `--no-PG -O SAM`.

Large-slice single-run timings with `/usr/bin/time -p`, stock installed
Samtools vs the migrated local build:

| Input | Command shape | Stock | Migrated | Speedup |
| --- | --- | ---: | ---: | ---: |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -c` | 1.02s | 0.92s | 1.11x |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -q 20 -c` | 1.05s | 0.92s | 1.14x |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -O SAM -o /dev/null` | 2.60s | 2.42s | 1.07x |
| `HG002.ont_ul.chr20_10-40Mb.bam` | `view -@8 -c` | 0.78s | 0.70s | 1.11x |

Interpretation: the migration is now broad enough to exercise `samtools view`
as a real consumer, and it reduces sys time substantially by avoiding normal
BGZF threaded input handoff.  It is not yet the Sambamba-scale win.  Under
`-@8`, decompression and CIGAR/query validation dominate enough that the
current ordered batch reader is only about 1.1x faster on these large local
slices.

Rejected experiment: enabling a parse-stage worker queue for batch record-view
construction dropped about one queue-depth of streaming records
(`216,942` stock vs `213,743` migrated on the highcov 1 Mb slice).  It was
backed out.  The next version needs explicit tests for post-header block
ownership, decoded-block queue draining, and EOF ordering before it can replace
the synchronous record-view builder.

## 2026-05-03: Count/Core Helpers and Raw BAM Frame Writer

Continued the no-fallback BAM `view` migration with HTSlib-internal helpers
intended to reduce avoidable per-record work without changing public API/ABI:

- `sam_bam_read_batch_count()` scans BAM frames for exact count-only consumers
  without building batch record views or materializing `bam1_t`;
- exact count header `tid`/`mtid` validation is now folded into the first frame
  scan instead of re-walking every frame in a second count-only pass;
- `sam_bam_read_core_batch()` exposes compact `bam1_core_t` batches for
  MAPQ/flag count filters and now reuses its caller-owned core buffer across
  reads to avoid per-batch core-array allocation churn;
- `sam_bam_batch_record_write1()` writes an unchanged BAM input frame directly
  to BAM output when there are no output mutations and no output index.  If the
  output handle is unsuitable, the same batch path materializes and uses the
  existing `sam_write1()` writer.

Rejected or bounded sub-ideas in this pass:

- CIGAR length caching was removed because it made the current workloads
  slower.
- Multi-BGZF-block aggregation was removed because the extra copy/lifetime
  bookkeeping cost more than it saved.
- Count/core-only fast paths intentionally do not cover BED, RG/tag, read-name,
  subsampling, expression/library, or `-m` cases.  Those still use the fuller
  batch record-view path.

Correctness checks performed in this pass:

- Product and stock BAM output on `HG00096.highcov.chr20_10-11Mb.bam` were
  decoded to SAM and `cmp` passed.
- HTSlib `make test` passed after this pass: 361 passed, 0 failed, and the
  libcurl retry tests passed.
- Samtools `make test` passed after this pass when linked against product
  HTSlib: 1,068 total, 1,022 passed, 0 failed, 46 expected failures; regression
  suites `mpileup.reg`, `depth.reg`, `cram-size.reg`, `consensus.reg`,
  threaded `consensus.reg`, and `cram_size.reg` all passed.

Benchmark caveat: by the time the raw-writer and reusable-core changes were
timed, the host was under heavy unrelated CPU load (`CyberhavenSystemMonitor`,
multiple Python processes, VS Code helpers, and other Codex sessions).  Treat
these numbers as guardrails, not final benchmark evidence.

Representative timings:

| Input | Command shape | Stock | Product | Interpretation |
| --- | --- | ---: | ---: | --- |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -c` | 1.46-1.69s | 1.54-1.62s | Neutral before host saturation. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -q 20 -c` | 1.44-1.50s | 1.41-1.50s | Neutral before host saturation. |
| `HG002.ont_ul.chr20_10-40Mb.bam` | `view -@8 -c` | 1.18-1.48s | 0.94-1.01s | Modest ONT win persisted. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@8 -b -o /dev/null` | 15.25s | 16.84s | Raw writer did not beat stock with output threads. |
| `HG00096.highcov.chr20_10-50Mb.bam` | `view -@1 -b -o /dev/null` | 72.46s | 66.06s | Raw writer helped slightly without output threads. |

Decision: keep the raw writer and reusable core-buffer pieces as useful
internal building blocks only if the full suites stay green.  They do not meet
the performance goal.  The next serious attempt should move below these
Samtools-visible consumers and fuse ordered BGZF inflate, frame scan, optional
core/view construction, and unchanged BAM output so the pipeline composes with
`-@` instead of competing with it.

## 2026-05-03: Fused Reader Correctness Scaffold and Segmented Batches

Implemented an opt-in ordered fused reader behind `HTS_BAM_BATCH_FUSED=1`.
The current version:

- defers `-@` into the batch reader when the fused env is enabled;
- keeps count/core/batch reads on the same fused queue when fused is active;
- serializes fused jobs (`BAM_STREAM_FUSED_MAX_IN_FLIGHT=1`) so cross-block
  BAM-record carry is known before the next physical block range is parsed;
- adds an internal segmented `bam_batch_t` so multi-block fused batches can
  reference decoded BGZF block storage without copying every record frame into
  one contiguous buffer;
- updates the raw batch writer and `test/sam` to validate segmented record
  ownership and live-batch lifetime.

Correctness checks passed after targeted fixes:

- `./test/sam`;
- full HTSlib `make test`: 361 passed, 0 failed, plus libcurl retry tests
  passed;
- full Samtools `make test` against product HTSlib: 1,068 total, 1,022
  passed, 0 failed, 46 expected failures; regression suites passed;
- fused dense SAM byte-compare on `HG00096.highcov.chr20_10-11Mb.bam`;
- fused split/large-record SAM byte-compare using `ce#large_seq.bam`;
- fused indexed-region SAM byte-compare for `20:10020000-10070000`;
- fused BAM output decoded back to SAM and matched default output.

Performance result: segmented fused is still not a candidate default because
the serialized boundary loses input parallelism.  On
`HG00096.highcov.chr20_10-50Mb.bam`, single runs were:

| Command shape | Installed stock | Product default | Fused segmented |
| --- | ---: | ---: | ---: |
| `view -@8 -b -o /dev/null` | 14.79s | 9.63s | 12.83s |
| `view -@8 -O SAM -o /dev/null` | 4.57s | 3.28s | 11.44s |

The useful win in this checkpoint is the default batch/raw-writer path, not the
fused experiment.  The fused code is useful mainly because it identifies the
next required design: safe parallel fused scheduling with cross-job carry
handling, likely via boundary proof or replay of prefetched jobs when carry is
reported.

## 2026-05-01: Prototype 1, Threaded BGZF Libdeflate Object Reuse

Implemented the first product prototype in `worktrees/htslib-bam-product`:
cache one libdeflate compressor/decompressor per HTSlib thread-pool worker for
threaded BGZF compression/decompression.

Implementation details:

- changed only `bgzf.c`;
- kept `bgzf_compress()` public API unchanged;
- kept all new helpers `static`;
- stored cache pointers on internal `mtaux_t`;
- used `hts_tpool_worker_id()` for per-worker lookup;
- tracked compression level per worker to avoid stale compressor reuse;
- freed cached objects in `mt_destroy()`;
- fixed setup-failure cleanup after the extra process-queue reference is added.

Correctness checks:

- `make lib-static test/test_view htsfile` passed in the product worktree;
- `git diff --check` passed in product and corpus worktrees;
- base/product `test_view -@ 4 -b` produced byte-identical BAM output for
  `HG00096.highcov.chr20_10-11Mb.bam` region `20:10000000-11000000`;
- both smoke outputs passed `samtools quickcheck`;
- both smoke outputs had 216,942 records.
- `make test` passed in the product worktree: 358 HTSlib tests passed, 0
  failed, and the libcurl retry tests passed.

Review notes:

- code review found no API/ABI or thread-safety blocker in the cache design;
- reviewer flagged a setup-failure leak in `bgzf_thread_pool()` after
  `hts_tpool_process_ref_incr()`;
- fixed by dropping the extra process-queue reference before destroying the
  queue on pre-thread-start errors.

Initial ad hoc benchmark:

- files: `data/benchmarks/libdeflate-cache-initial.tsv` and
  `data/benchmarks/libdeflate-cache-initial-hires.tsv`;
- conclusion: non-hires timing was too coarse; high-resolution data was
  neutral-to-marginal and needed a captured harness.

Added `bench/bam-shape/scripts/run_htslib_test_view_bench.pl` to capture
reproducible product-vs-base `test_view` timings with metadata, checks, command
strings, repeated loops, and timing status. Also corrected the existing
Samtools/Sambamba harness default so `SAMTOOLS` points at the local Samtools
binary rather than `test_view`.

Focused benchmark:

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

Decision: keep the patch available for now, but do not count it as a high-ROI
optimization. The evidence says this is correct and low-friction, but the
measured win on the current corpus is under 1%. Revisit only if a larger or
more compressor-heavy workload shows a clearer effect.

## 2026-05-01: Prototype 2, BGZF Default Queue Depth

Tested changing the default threaded BGZF queue depth from `2 * threads` to
`4 * threads` when callers pass `qsize=0`.

Correctness checks:

- product rebuild passed with `make lib-static test/test_view htsfile`;
- `git diff --check` passed;
- base/product `test_view -@ 4 -b` produced byte-identical BAM output for
  `HG00096.lowcov.chr20_10-20Mb.bam` region `20:10000000-20000000`;
- both outputs passed `samtools quickcheck`;
- the product output had 486,569 records.

Focused benchmark:

```text
data/benchmarks/bgzf-qdepth-focused-readwrite/
data/benchmarks/bgzf-qdepth-focused-decode/
```

Setup:

- input: `exome_chr20`, region `20`, 1,933,817 records;
- base/product HTSlib commit: `2b3ddfd0`;
- product status during benchmark: `M bgzf.c; M htslib/bgzf.h`;
- libdeflate: `1.25`;
- read/write: `REPEATS=5`, `LOOPS=2`, `THREADS_LIST="2 4 8"`;
- decode-only: `REPEATS=5`, `LOOPS=10`, `THREADS_LIST="2 4 8"`.

Median seconds per `test_view` invocation:

| Mode | Threads | Base | Product | Speedup |
| --- | ---: | ---: | ---: | ---: |
| decode | 2 | 0.323232 | 0.322005 | 0.38% |
| decode | 4 | 0.175695 | 0.173606 | 1.19% |
| decode | 8 | 0.092243 | 0.087210 | 5.46% |
| readwrite | 2 | 2.480360 | 2.466350 | 0.56% |
| readwrite | 4 | 1.292460 | 1.281930 | 0.81% |
| readwrite | 8 | 0.682189 | 0.674737 | 1.09% |

Review outcome:

- code review found no ABI break, but rated the global behavior/memory default
  change as too broad for the current evidence;
- each `bgzf_job` contains compressed and uncompressed BGZF buffers, so a
  higher default can materially increase memory when many files are threaded;
- `qsize=0` semantics would also diverge between BGZF and SAM/text threading;
- benchmark review noted that the product branch also contained prototype 1, so
  queue-depth attribution was not isolated.

Decision: revert the default queue-depth change. Do not keep this as an
upstreamable default without an isolated qsize A/B, 8+ thread coverage across
more workloads, and a memory/RSS table. Future queue-depth work should use an
explicit qsize path or an experiment-only harness, not altered `qsize=0`
semantics.

## 2026-05-01: Prototype 2b, Opt-In BGZF Queue Size

Reworked the queue-depth idea into an opt-in `HTS_BGZF_QUEUE_SIZE` override
used only when callers pass `qsize=0`. The default remains `2 * threads` when
the environment variable is absent or invalid.

Implementation details:

- changed `bgzf.c`, `htslib/bgzf.h`, `htslib/faidx.h`, and `htslib/hts.h`;
- public default behavior is unchanged unless the environment variable is set;
- public docs now mention the environment override wherever `qsize=0` behavior
  is documented;
- invalid environment values fall back to the default and log only at debug
  level;
- libdeflate cache arrays are now direction-specific, so read handles allocate
  decompressor slots only and write handles allocate compressor/level slots
  only.

Correctness checks:

- product rebuild passed with `make lib-static test/test_view htsfile`;
- `git diff --check` passed;
- targeted `make check TEST_OPTS='test_view sam'` passed: 361 tests passed,
  0 failed, and the libcurl retry tests passed;
- product `test_view -@ 8 -b` with default queue size and with
  `HTS_BGZF_QUEUE_SIZE=32` produced byte-identical BAM output for
  `HG00096.highcov.chr20_10-11Mb.bam` region `20:10000000-11000000`;
- both outputs passed `samtools quickcheck`;
- the q32 output had 216,942 records.

Isolated same-binary benchmark:

```text
data/benchmarks/bgzf-qsize-env-isolated/
```

Setup:

- input: `HG00096.exome.chr20.bam`, region `20`;
- product binary only, so libdeflate-cache code is held constant;
- `threads=8`;
- queue cases interleaved: default, `HTS_BGZF_QUEUE_SIZE=32`,
  `HTS_BGZF_QUEUE_SIZE=64`;
- decode-only: `REPEATS=7`, `LOOPS=10`;
- read/write: `REPEATS=7`, `LOOPS=2`.

Median seconds per `test_view` invocation:

| Mode | Queue | Median |
| --- | ---: | ---: |
| decode | default | 0.096737 |
| decode | 32 | 0.092395 |
| decode | 64 | 0.090937 |
| readwrite | default | 0.678635 |
| readwrite | 32 | 0.676319 |
| readwrite | 64 | 0.677719 |

Decision: keep this as an opt-in experiment/advanced tuning knob, not a default
change. The qsize signal is real for 8-thread decode on the exome workload, but
read/write is neutral and RSS/many-file costs still need measurement before any
stronger recommendation.

## 2026-05-01: Rejected Experiment, BGZF Worker-ID TLS Cache

Tested caching `hts_tpool_worker_id()` in `bgzf.c` using C11 `_Thread_local`
state to avoid a per-block linear scan of thread-pool worker IDs.

Review outcome:

- code review found no direct thread-safety or ABI issue;
- reviewer flagged C11 TLS portability/configuration risk and noted that the
  helper was only needed in `HAVE_LIBDEFLATE` builds;
- benchmark review noted that the first base-vs-product timing did not isolate
  the worker-id cache from the rest of the product patch.

Isolated same-binary benchmark:

```text
data/benchmarks/bgzf-worker-id-cache-isolated/
```

Setup:

- input: `HG00096.exome.chr20.bam`, region `20`;
- product binary only, with a temporary `HTS_BGZF_WORKER_ID_CACHE=0` toggle;
- `threads=8`;
- decode-only: `REPEATS=9`, `LOOPS=12`;
- read/write: `REPEATS=9`, `LOOPS=3`.

Median seconds per `test_view` invocation:

| Mode | Worker-ID cache | Median |
| --- | --- | ---: |
| decode | off | 0.096948 |
| decode | on | 0.097073 |
| readwrite | off | 0.688695 |
| readwrite | on | 0.688795 |

Decision: remove the worker-id TLS cache. It is not measurably beneficial on
the current workload and adds portability noise.

## 2026-05-01: Queue Size RSS Matrix

Added `bench/bam-shape/scripts/run_bgzf_qsize_matrix.pl`, which runs the
product `test_view` binary with rotated qsize order, warmups, record-count
checks, and `/usr/bin/time -l` maximum resident set size.

Current result:

```text
data/benchmarks/bgzf-qsize-matrix-20260501/
```

Setup:

- inputs: all five local BAM corpus entries;
- modes: decode and read/write;
- threads: 4 and 8;
- qsizes: default, 32, 64;
- `REPEATS=5`, `LOOPS=3`, `WARMUPS=1`.

Key results:

- `exome_chr20`, decode, `-@ 8`: q32 was 5.03% faster than default with
  +3.12 MiB max RSS; q64 was 6.87% faster with +9.28 MiB max RSS.
- other 8-thread decode workloads were mixed: highcov improved about 1.6-1.9%,
  lowcov moderate improved about 0.8-1.5%, tiny lowcov regressed, and the ONT
  slice was neutral.
- read/write remained mostly neutral on larger workloads; q64 sometimes helped
  smaller tail-latency-sensitive runs but at larger RSS cost.

Decision: keep `HTS_BGZF_QUEUE_SIZE` as an advanced opt-in tuning knob. The
matrix supports the knob for some decode-heavy 8-thread cases, but not a
default change.

## 2026-05-01: Prototype 3, Direct Threaded BGZF Block Read

Implemented a smaller version of the threaded BGZF read simplification:

- kept the non-consuming `hpeek()` header probe so non-BGZF gzip fallback
  remains possible for existing control flow;
- after validating a BGZF header, read the full compressed block directly into
  `j->comp_data` with one `hread()`;
- removed the threaded `load_block_from_cache()` shortcut because it mutates
  `fp` state rather than populating the `bgzf_job`.

Review notes:

- first review rejected a more aggressive no-peek version because it consumed
  bytes before the `BGZF_ERR_MT` fallback path;
- revised review found no blocking correctness issue;
- reviewer requested raw gzip and cache-before-MT coverage.

Tests added:

- `test_write_read(..., "wg", ..., nthreads=2)` covers requesting MT on a raw
  gzip stream, which should remain single-threaded;
- `test_cache_before_mt()` populates the BGZF block cache before enabling MT,
  seeks back to the beginning, and reads with MT enabled.

Focused decode benchmark:

```text
data/benchmarks/bgzf-direct-read-focused-decode/
```

Median decode seconds per `test_view` invocation:

| Threads | Input | Base | Product | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 4 | `exome_chr20` | 0.173316 | 0.173869 | -0.32% |
| 4 | `highcov_chr20_10_11mb` | 0.067686 | 0.067415 | 0.40% |
| 4 | `lowcov_chr20_10_20mb` | 0.053344 | 0.052928 | 0.78% |
| 8 | `exome_chr20` | 0.095374 | 0.094776 | 0.63% |
| 8 | `highcov_chr20_10_11mb` | 0.039651 | 0.038595 | 2.66% |
| 8 | `ont_ul_chr20_10_10_2mb` | 0.014619 | 0.014290 | 2.25% |

Tiny inputs were noisy and sometimes slower. Decision: keep for now as a
small correctness-cleanup/perf-neutral BGZF internal improvement, not as a
headline optimization.

## 2026-05-01: Prototype 4, Threaded Writer Flush Condition Variable

Replaced `mt_flush_queue()` polling on `jobs_pending` with an internal
`pthread_cond_t` signalled by `bgzf_mt_writer()` when queued writer jobs
complete. The wait remains timed, but normal completion no longer pays the
fixed `hts_usleep(10000)` granularity.

Implementation details:

- added `job_pool_c` to internal `mtaux_t`;
- signal after normal `jobs_pending--`;
- signal after `mt_queue()` dispatch failure decrements `jobs_pending`;
- broadcast after writer error shuts down the process queue, so a flush waiter
  can observe shutdown;
- destroy the condition only after `pthread_join()` in `mt_destroy()`.

Review outcome:

- first concurrency review found no normal-path deadlock but noted the error
  path should wake the condition too;
- follow-up review approved moving the error-path broadcast after
  `hts_tpool_process_destroy()` and removing an unnecessary `mt_destroy`
  broadcast.

Correctness checks:

- `make lib-static test/test_view htsfile test/test_bgzf` passed;
- `./test/test_bgzf test/bgziptest.txt` passed;
- `git diff --check` passed;
- full `make test` passed: 361 tests passed, 0 failed, and the libcurl retry
  tests passed;
- highcov `test_view -@ 8 -b` output matched base byte-for-byte and passed
  `samtools quickcheck`;
- highcov `HTS_BGZF_QUEUE_SIZE=1 test_view -@ 8 -b` output also matched base
  byte-for-byte and passed `samtools quickcheck`.

Focused benchmark:

```text
data/benchmarks/bgzf-cond-flush-focused-readwrite-final/
```

Setup:

- inputs: all five local BAM corpus entries;
- mode: read/write (`test_view -b -p /dev/null`);
- threads: 4 and 8;
- `REPEATS=5`, `LOOPS=3`.

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

Decision: keep. This is the strongest HTSlib-level win so far. It is
ABI-neutral, applies to all threaded BGZF writers, and especially improves
tail latency for shorter BAM outputs while still helping larger workloads.

## 2026-05-01: Prototype 5, BAM Write Packing Fast Path

Implemented a narrow BAM-local fast path in `bam_write1()` for the common case:

- little-endian host;
- BGZF-compressed, non-gzip output;
- `n_cigar <= 0xffff`;
- the complete BAM record fits in the current BGZF uncompressed block with
  spare room.

The fast path copies block length, BAM core, qname, and payload tail directly
into `fp->uncompressed_block` and advances `block_offset` once. Exact-fill
cases, big-endian, gzip, uncompressed output, and long-CIGAR records keep using
the old sequence of `bgzf_write_small()` calls.

Review outcome:

- initial review found no byte-order/API issue, but requested an in-tree parity
  test and noted a long-CIGAR indexed-write preflush size mismatch;
- added a `test/sam.c` payload parity test that writes identical BAM records as
  normal BGZF BAM and gzip BAM, then compares the decompressed BAM byte stream;
- fixed `bam_write_idx1()` to include the `+16` long-CIGAR expansion before
  `bgzf_flush_try()` and index amendment.

Correctness checks:

- `make lib-static test/test_view htsfile test/sam test/test_bgzf` passed;
- `./test/sam` passed;
- `git diff --check` passed;
- highcov `test_view -@ 8 -b` output matched the pre-fast-path product
  byte-for-byte, including with `HTS_BGZF_QUEUE_SIZE=1`.

Focused benchmark:

```text
data/benchmarks/bam-write-pack-focused-readwrite-final/
```

Setup:

- baseline: saved product `test_view` binary from before the `sam.c` fast path;
- product: current product `test_view`;
- inputs: all five local BAM corpus entries;
- mode: read/write (`test_view -b -p /dev/null`);
- threads: 4 and 8;
- `REPEATS=5`, `LOOPS=3`.

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

Decision: marginal. This is low blast radius after the parity test, but it is
not a high-ROI optimization on the current corpus. Keep it only as a small
cleanup candidate; the flush condition variable remains the stronger upstream
story.

## 2026-05-03: No-fallback Samtools view Migration and Async Parser v2

Samtools `view` now treats BAM batch reading as the normal BAM path in this
product branch rather than a compile-time experiment:

- removed the `HTS_BAM_BATCH_READER_CONSUMER` guards from `sam_view.c`;
- include `sam_internal.h` unconditionally in `sam_view.c`;
- track `$(HTSDIR)/sam_internal.h` unconditionally in `samtools/Makefile`;
- request `sam_bam_prepare_batch_reader()` after reading a BAM header and
  before `-@` is applied;
- route BAM streaming, indexed region, multi-region, BED, `--fetch-pairs`,
  `-U`, mutation, and output cases through `sam_bam_read_batch()` /
  `sam_bam_itr_next_batch()`;
- keep SAM/CRAM input on their existing reader paths;
- error if the BAM batch reader cannot be used, rather than falling back to the
  old BAM `sam_read1()` loop.

HTSlib follow-up:

- added an ordered async record-view parser path behind
  `HTS_BAM_BATCH_PARSE=1`;
- fixed EOF draining, queued decoded-block ownership, post-header serial block
  copying, valid-prefix error ordering, and worker-side core-validation prefix
  handling;
- kept it off by default because it is slower on the dense count workload;
- added a stream-vs-iterator distinction so plain streaming batches skip
  per-record virtual-offset calculation, while iterator reads still request
  offsets.

Validation:

- `worktrees/htslib-bam-product`: `make libhts.a test/sam` passed;
- `worktrees/htslib-bam-product`: full `make test` passed with 361 passed,
  0 failed;
- `samtools`: full `make test` passed linked against product HTSlib without
  any batch-consumer CPP flag: 1,068 tests, 0 failed, 46 expected failures, and
  regression suites `PASS`;
- manual parity:
  - highcov streaming count: `7,579,904`;
  - highcov `-q 20 -c`: `7,396,382`;
  - ONT count: `156,695`;
  - highcov indexed 1 Mb region: `216,942`;
  - `HTS_BAM_BATCH_PARSE=1` matched the same highcov and ONT counts/region;
  - `--no-PG -O SAM` highcov region output matched stock byte-for-byte.

Single-run timing is noisy on this machine.  The latest paired run after the
no-fallback cleanup, Makefile dependency fix, and configured no-macro rebuild
had installed stock highcov `view -@8 -c` at 1.18s and product batch at 1.16s;
ONT 30 Mb was 0.87s stock vs 0.79s product.  Earlier paired runs showed up to
1.29s stock vs 1.11s product on the same highcov command.  The async parser
mode was slower at 2.40-2.58s for the same command, so the next performance
step should be a lighter internal count/scan mode or a fused decode+view
worker rather than a second queue that reparses decoded blocks.

## 2026-05-03: Whole-batch Raw BAM Writer and Pruned Threading Experiments

Implemented a narrower all-pass streaming writer on top of the migrated batch
reader:

- HTSlib internal API: `sam_bam_batch_write1()` writes a full contiguous
  `bam_batch_t` as one raw BAM frame span.
- Samtools consumer: `view` uses this only when every record in a streaming
  batch should pass and output is unchanged BAM with no output index.
- Guardrails: the path rejects non-BAM output, output indexing, filters,
  sanitization, tag edits, `-U`, `-p`, flag mutation, and count mode.  It also
  validates QNAME NUL state and CIGAR/query-length consistency before writing.
  QNAME repair falls back to the materialized batch path for that batch, not
  to old BAM `sam_read1()`.

Correctness:

- highcov all-pass product BAM quickchecked;
- decoded product SAM matched stock decoded SAM on
  `HG00096.highcov.chr20_10-11Mb.bam`;
- indexed region output matched stock for `20:10020000-10070000`;
- product HTSlib `./test/sam` passed;
- `git diff --check` passed in product HTSlib, Samtools, and corpus worktrees.

Directional dense write timings before heavy host saturation:

| Mode | Dense `view -@8 -b -o /dev/null` real time |
| --- | ---: |
| product whole-batch raw writer | 12.54s |
| stock installed Samtools | 13.10s |
| product after raw validation guard | 16.03s |
| stock paired with validation-guard run | 17.32s |

Interpretation: the whole-batch writer is correct and low friction, but the
win is small in the threaded dense regime and still far below the 1.5x target.

Pruned experiments:

- Auto-defaulting the async record-view parser was reverted.  It remained
  parity-clean, but a same-run dense write A/B was not favorable:
  `HTS_BAM_BATCH_PARSE=0` 11.21s, `HTS_BAM_BATCH_PARSE=1` 11.52s, stock
  11.42s.  The parser remains opt-in only.
- Reusing stock BGZF MT through the batch reader was prototyped behind an env
  flag.  The safe copy-based version took 12.49s.  A private retained decoded
  block handle avoided the copy and passed small-slice parity, but took 30.58s
  on the dense write benchmark.  Both variants were removed from product code.

Final benchmark attempt in this session was discarded as non-decision-quality:
the host had multiple unrelated Python/AWK/Codex processes active, and a paired
dense write run measured 49.47s product vs 37.46s stock.  Treat earlier
pre-saturation pairs as more informative, but still directional.

## 2026-05-05: HTSlib-Internal Region Iterator Generalization

Goal: move the messy-BED `view -L` speedup out of a Samtools-only helper and
into HTSlib internals without changing public HTSlib API/ABI.

Implementation:

- added private `sam_internal.h` helpers for copied `hts_reglist_t`
  coalescing and exact-overlap filtered serial/BAM-batch iterator reads;
- implemented the helpers in `sam.c`;
- changed `samtools/region_plan.c` to build the ordinary BED `hts_reglist_t`
  and ask HTSlib for a coalesced query copy;
- changed `samtools/sam_view.c` to pass the exact BED predicate into the
  HTSlib filtered iterator wrappers for the auto-indexed messy-BED path;
- kept BED parsing, special-reference handling, auto-index eligibility, and
  fallback policy in Samtools.
- tightened the review guardrails so empty BEDs and special-reference BEDs
  stay on the streaming path, and count-only `view -c -L` uses the same exact
  filtered iterator when auto-indexing is eligible.

Validation:

- HTSlib `./test/sam` passed, including a new coalesced-reglist plus filtered
  serial/batch iterator test.
- HTSlib `make -j8` passed.
- Samtools `make -j8` passed.
- Samtools `env PATH=.../worktrees/htslib-bam-product:$PATH ./test/test.pl view merge`
  passed: 1141 total, 1091 passed, 50 expected failures, 0 unexpected failures.
- `git diff --check` passed in product HTSlib and Samtools.
- Dense WGS, exome, and ONT messy-BED decoded SAM outputs matched installed
  production Samtools exactly.

Fresh single warm-cache timings.  These are directional local timings, not yet
the repeated median corpus harness output:

| Input | Command | Stock | Sambamba | Product |
| --- | --- | ---: | ---: | ---: |
| dense highcov | `view -@8 -b -L messy -o /dev/null` | 1.90s | 2.04s | 0.87s |
| dense highcov | `view -@8 -u -L messy -o /dev/null` | 1.24s | 0.27s | 0.26s |
| exome | `view -@8 -b -L messy -o /dev/null` | 0.16s | 0.10s | 0.06s |
| exome | `view -@8 -u -L messy -o /dev/null` | 0.14s | 0.03s | 0.03s |
| ONT | `view -@8 -b -L messy -o /dev/null` | 1.36s | 1.65s | 0.68s |
| ONT | `view -@8 -u -L messy -o /dev/null` | 0.81s | 0.22s | 0.18s |

Interpretation:

- The migration preserves the previous region speedup while making the fast
  shape reusable by other Samtools commands.
- The retained abstraction is deliberately narrow: HTSlib owns generic
  coalescing and exact-filtered ordered iterator emission; Samtools owns BED
  semantics and command policy.
- `merge -L` is already the safest additional consumer.  The next practical
  not-yet-migrated consumer to investigate is probably `depth -b`, but only as
  an indexed query shortcut that keeps per-position BED filtering and no-index
  streaming.

## 2026-05-05: `samtools depth -b` Region-Plan Migration

Implemented the next BED-driven region consumer on top of the shared Samtools
region-plan machinery.

Implementation:

- `bam2depth.c` includes `region_plan.h` and opportunistically builds
  `sam_itr_regions()` iterators from merged BED reglists for default
  `depth -b`.
- The optimization is gated to `-b` without `-r`, without `-a/-aa`, with a
  non-empty and non-special BED, and with an index available for every BAM/CRAM
  input.
- It also requires all BED reference names to exist in each input header; if not,
  the command keeps the legacy streaming path and quiet stderr.
- If the fast path is ineligible or an index is missing, `depth` silently keeps
  the existing streaming path.  The existing required-index behavior for `-r`
  is unchanged.
- Multi-region iterators are prevented from initializing the single-region
  `beg/end/tid` clipping state.
- Exact output semantics remain in the old pileup path: `zero_region()` and
  `add_depth()` still call `bed_overlap()` per output position.

Tests:

- Added `test_large_positions` coverage comparing indexed `depth -b` output
  against forced no-index streaming output for messy/overlapping/duplicate BEDs,
  filters, `-a`, `-aa`, empty BEDs, special-reference BED mixed with a real
  interval, and multiple inputs.
- Added explicit empty-BED `-H` header coverage.
- Added `-r + -b` and CRAM `-X` coverage.
- Added missing-BED-reference coverage that compares fallback output and checks
  stderr stays empty.

Validation:

- `git diff --check` passed in Samtools.
- `make -j8` passed in product HTSlib and Samtools.
- Full Samtools `env PATH=.../worktrees/htslib-bam-product:$PATH ./test/test.pl`
  passed: 1153 total, 1103 passed, 50 expected failures, 0 unexpected
  failures.
- Production/product output matched byte-for-byte for dense highcov, exome,
  and ONT messy-BED depth outputs.

Repeated-median timings, five warm-cache repeats:

| Threads | Input | Stock | Product | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 0 | dense highcov | 7.77s | 0.99s | 7.85x |
| 0 | exome | 0.81s | 0.11s | 7.36x |
| 0 | ONT | 6.97s | 0.92s | 7.58x |
| 8 | dense highcov | 1.36s | 0.32s | 4.25x |
| 8 | exome | 0.28s | 0.07s | 4.00x |
| 8 | ONT | 1.88s | 0.37s | 5.08x |

Sambamba `depth base -t8 -L` was also timed as context with a filter closer to
Samtools depth defaults, but its output schema and defaults differ.  Three-run
medians were 4.45s dense, 1.06s exome, and 4.23s ONT, all slower than the
product path.

Decision: this is the clearest common-operation win so far.  It does not make
the pileup algorithm faster; it avoids reading most irrelevant records for the
common target-panel or messy-BED depth workflow.  Keep `-a/-aa` on the streaming
path until zero-fill semantics are explicitly redesigned and proven.

## 2026-05-05: Selective Compressed Filtered BAM Output

Goal: improve the common `samtools view -b FILTERS -o out.bam` path without
changing public HTSlib/Samtools APIs or weakening fallback semantics.

Implementation kept:

- `sam_view.c` delays batch record-view construction until a filter needs
  QNAME or aux data; core-field-only filters avoid the per-record view setup.
- The internal HTSlib raw-range writer resolves the output `BGZF *` once per
  range-write call.
- Sparse selected raw ranges can be packed into one temporary buffer per batch
  before writing to BGZF.  This is private to `sam_internal.h`; no public header
  or export map changed.

Implementation tried and removed:

- An adaptive raw path for compressed pure flag-only filters.  It regressed the
  dense `-f 4` case, so pure flag-only compressed filters keep the existing
  materialized write path unless another selective planner/filter path applies.

Validation:

- Product HTSlib `./test/sam` passed.
- Samtools `make check` passed: 1155 total, 1105 passed, 50 expected failures,
  0 unexpected failures, and all regression scripts PASS.
- Added Samtools tests for threaded compressed BAM output with a sparse
  expression filter and a pure flag filter.
- Dense, exome, and ONT real-output BAMs quickchecked; decoded SAM equality was
  explicitly checked for dense `-e 'mapq < 20'`, exome `-e 'mapq < 20'`, and
  ONT `-q 20 -F 4`.
- `git diff --check` passed in product HTSlib and Samtools.

Directional paired timings, real BAM output under `/private/tmp`:

| Input | Threads | Filter | Stock | Product | Result |
| --- | ---: | --- | ---: | ---: | --- |
| dense | 0 | `-q 20 -F 4` | 60.86s | 60.31s | tied |
| dense | 8 | `-q 20 -F 4` | 8.54s | 8.06s | 1.06x |
| dense | 0 | `-e 'mapq < 20'` | 8.31s | 7.84s | 1.06x |
| dense | 8 | `-e 'mapq < 20'` | 1.49s | 1.01s | 1.48x |
| exome | 0 | `-q 20 -F 4` | 5.02s | 4.99s | tied |
| exome | 8 | `-q 20 -F 4` | 0.65s | 0.65s | tied/noisy |
| exome | 0 | `-e 'mapq < 20'` | 0.79s | 0.65s | 1.22x |
| exome | 8 | `-e 'mapq < 20'` | 0.23s | 0.10s | positive but short/noisy |
| ONT | 0 | `-q 20 -F 4` | 34.66s | 34.14s | tied |
| ONT | 8 | `-q 20 -F 4` | 4.67s | 4.52s | 1.03x |
| ONT | 0 | `-e 'mapq < 20'` | 18.52s | 18.44s | tied |
| ONT | 8 | `-e 'mapq < 20'` | 2.49s | 2.41s | small win |

Conclusion: selective filtered BAM output is now meaningfully faster when the
filter is sparse and can use the batch/raw record path.  Dense high-pass
compressed BAM output is still BGZF-compression-bound and should be described
as neutral to slightly faster, not as a broad compressed-BAM rewrite win.
