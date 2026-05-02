# Ordered BAM Reader Benchmark Notes

## Scope

This note records the current state of the experimental HTSlib-level ordered
BAM reader.  The reader is internal and opt-in via:

- `HTS_BAM_ORDERED_READER=1`
- `-@ N`, `hts_set_threads()`, or `HTS_OPT_THREAD_POOL` for the thread budget
- `HTS_BAM_ORDERED_READER_THREADS=N` as a fallback when the caller did not
  configure threads
- `HTS_BAM_ORDERED_READER_REQUIRE=1` for strict activation during tests

The implementation preserves the public API.  It uses local BAI/CSI virtual
offsets to split a BAM into non-overlapping file spans, parses those spans on
worker file handles, and emits records through the existing `sam_read1()` stream
in original file order.

## Chunked Drain Experiment

A chunked worker-output variant was tested to see whether the reader could start
draining large jobs before a full span finished.  Correctness checks passed, but
benchmarking showed that the per-chunk queue and synchronization overhead
dominated while callers still consumed one `bam1_t` at a time through
`sam_read1()`.

Reviewer findings also identified two design problems in the chunked variant:

- The current drain job could enqueue unbounded chunks, so the configured global
  chunk bound was not a true memory bound.
- Worker-level failures used a global error flag and could overtake already
  queued earlier records.

Decision: prune the chunked queue layer from the default implementation.  The
remaining reader uses the lower-overhead full-job result slots.

## Benchmarks

Original command shape:

```sh
./test/test_view -@ THREADS -B -p /dev/null BAM
```

Ordered mode additionally sets:

```sh
HTS_BAM_ORDERED_READER=1 \
HTS_BAM_ORDERED_READER_REQUIRE=1 \
HTS_BAM_ORDERED_READER_THREADS=THREADS
```

Each row is the median of five runs.  Raw timings:

- `/tmp/htslib_bam_ordered_reader_chunk65536_timings.tsv`
- `/tmp/htslib_bam_ordered_reader_pruned_default_timings.tsv`
- `/tmp/htslib_bam_ordered_reader_fulljob_final_timings.tsv`

Final full-job reader results:

| input | threads | off | ordered | delta |
| --- | ---: | ---: | ---: | ---: |
| exome_chr20 | 1 | 0.601 | 0.730 | -21.46% |
| exome_chr20 | 2 | 0.321 | 0.416 | -29.60% |
| exome_chr20 | 4 | 0.174 | 0.311 | -78.74% |
| exome_chr20 | 8 | 0.179 | 0.366 | -104.47% |
| highcov_1mb | 1 | 0.224 | 0.239 | -6.70% |
| highcov_1mb | 2 | 0.126 | 0.132 | -4.76% |
| highcov_1mb | 4 | 0.069 | 0.085 | -23.19% |
| highcov_1mb | 8 | 0.039 | 0.080 | -105.13% |
| lowcov_10_20mb | 1 | 0.170 | 0.189 | -11.18% |
| lowcov_10_20mb | 2 | 0.097 | 0.112 | -15.46% |
| lowcov_10_20mb | 4 | 0.055 | 0.080 | -45.45% |
| lowcov_10_20mb | 8 | 0.054 | 0.103 | -90.74% |
| lowcov_1mb | 1 | 0.029 | 0.029 | +0.00% |
| lowcov_1mb | 2 | 0.019 | 0.019 | +0.00% |
| lowcov_1mb | 4 | 0.012 | 0.019 | -58.33% |
| lowcov_1mb | 8 | 0.015 | 0.020 | -33.33% |
| ont_ul_200kb | 1 | 0.053 | 0.052 | +1.89% |
| ont_ul_200kb | 2 | 0.031 | 0.033 | -6.45% |
| ont_ul_200kb | 4 | 0.022 | 0.022 | +0.00% |
| ont_ul_200kb | 8 | 0.015 | 0.018 | -20.00% |

Those historical results compared ordered parsing against an already
BGZF-threaded baseline before `-@` was wired into the ordered reader.  In that
mode the ordered reader used `HTS_BAM_ORDERED_READER_THREADS`, so it did not
consume the caller's `-@` thread budget.

## Ordered Reader Without Input BGZF Threads

After reviewer feedback, `bam_ordered_process_job()` was changed to stop taking
the shared reader mutex once per record just to check for shutdown.  It now
checks for close requests every 1024 records.  The ordinary non-`-@` read path
then shows a broad win on the local benchmark slices.

Command shape:

```sh
./test/test_view -B -p /tmp/out BAM
HTS_BAM_ORDERED_READER=1 \
HTS_BAM_ORDERED_READER_REQUIRE=1 \
HTS_BAM_ORDERED_READER_THREADS=4 \
./test/test_view -B -p /tmp/out BAM
```

Each row is the median of seven runs.  Raw timings:

- `/tmp/htslib_bam_ordered_mutex_before.tsv`
- `/tmp/htslib_bam_ordered_mutex_after_pruned.tsv`
- `/tmp/htslib_bam_ordered_mutex_final.tsv`

| input | ordinary | ordered | speedup |
| --- | ---: | ---: | ---: |
| highcov chr20 10-11Mb | 0.20 | 0.08 | 2.5x |
| exome chr20 | 0.64 | 0.19 | 3.4x |
| ONT ultra-long chr20 10-10.2Mb | 0.04 | 0.02 | 2.0x |

## Interpretation

The ordered reader now clears the 2x target for default single-stream
`sam_read1()` read-discard workloads on the local BAM corpus.  Separately,
`test/sam` and a decoded SAM `cmp` on the high-coverage slice verified that the
reader preserves record order and contents through the existing API.  It should
still remain opt-in for now: it depends on a local index, has not shown
write-pipeline wins, and still pays per-record ownership handoff costs.  The
next ordered-reader work should focus on worker buffer reuse and output/write
workloads.

## `-@` Integration

`hts_set_threads()` and `HTS_OPT_THREAD_POOL` now route BAM through the
SAM-layer thread setter instead of attaching BGZF threading immediately.  For
BAM input with `HTS_BAM_ORDERED_READER=1`, HTSlib stores a small deferred thread
configuration.  On the first `sam_read1()` call it tries to open the ordered
reader with that `-@` thread count.  If the ordered reader cannot activate and
strict mode is not set, HTSlib converts the deferred configuration into the
ordinary BGZF threaded path and continues reading normally.

Current caveat: when the caller supplies `HTS_OPT_THREAD_POOL`, the ordered
reader uses the pool size as its worker count but still creates its own worker
threads and worker file handles internally.  It does not dispatch jobs onto the
caller's shared `hts_tpool`, so this is a thread-budget integration rather than
a full shared-pool integration.

Correctness checks:

```sh
./test/sam
./test/test_view -@ 4 -p /tmp/base.sam BAM
HTS_BAM_ORDERED_READER=1 HTS_BAM_ORDERED_READER_REQUIRE=1 \
  ./test/test_view -@ 4 -p /tmp/ordered.sam BAM
cmp /tmp/base.sam /tmp/ordered.sam
```

Local median-of-five smoke timings for read-discard with `-@ 4` after this
integration:

| input | BGZF `-@ 4` | ordered `-@ 4` | result |
| --- | ---: | ---: | ---: |
| highcov chr20 10-11Mb | 0.06 | 0.09 | slower |
| exome chr20 | 0.17 | 0.27 | slower |
| ONT ultra-long chr20 10-10.2Mb | 0.02 | 0.03 | noisy/slower |

These results mean `-@` now controls the ordered reader, but the existing BGZF
threaded path remains faster for already-threaded read-discard workloads when
the ordered reader is disabled or cannot activate.

## BGZF Threaded Libdeflate Cache

The BGZF worker path now reuses one libdeflate compressor/decompressor per
worker thread instead of allocating one per BGZF block.  This is an internal
optimization for ordinary threaded BGZF read/write workloads, not a headline
Sambamba-style win.

Benchmark command shapes:

```sh
./test/test_view -@ THREADS -B -p /tmp/out BAM
./test/test_view -@ THREADS -b -p /tmp/out.bam BAM
```

A temporary opt-out gate was used only for measurement:
`HTS_BGZF_LIBDEFLATE_CACHE=0`.  It was removed after benchmarking.  Each row is
the median of five runs.  Raw timings:

- `/tmp/htslib_bgzf_libdeflate_cache_timings.tsv`

| input | workload | threads | cache off | cache on | speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| highcov chr20 10-11Mb | read discard | 4 | 0.082 | 0.082 | 1.00x |
| highcov chr20 10-11Mb | read discard | 8 | 0.055 | 0.055 | 1.01x |
| highcov chr20 10-11Mb | BAM write | 4 | 0.488 | 0.483 | 1.01x |
| highcov chr20 10-11Mb | BAM write | 8 | 0.266 | 0.263 | 1.01x |
| exome chr20 | read discard | 4 | 0.189 | 0.184 | 1.03x |
| exome chr20 | read discard | 8 | 0.205 | 0.194 | 1.06x |
| exome chr20 | BAM write | 4 | 1.269 | 1.262 | 1.01x |
| exome chr20 | BAM write | 8 | 0.667 | 0.659 | 1.01x |

Decision: keep as a neutral-to-small local improvement if review remains clean,
but do not position it as a major standalone result.

## Pruned Write Fast Path

`bam_write1_block_copy()` was tested as a BAM writer micro-optimization that
collapses several `bgzf_write_small()` calls into direct writes into the current
BGZF uncompressed block.  It did not produce a reliable win in normal
`test/test_view -b` BAM-copy workloads, so it was pruned from the product patch.

Raw timings:

- `/tmp/htslib_bam_write_block_copy_gate_timings.tsv`

Median results were neutral or slightly slower in most cells, with only a tiny
highcov 8-thread improvement.  The result is not worth carrying next to the much
larger raw BGZF pass-through primitive.

## Raw BGZF BAM Pass-Through Experiment

After the ordered-reader and per-record raw-copy experiments failed to produce a
large win, an internal BAM-to-BAM pass-through primitive was added:

- `sam_bam_raw_copy_blocks(htsFile *in, htsFile *out)` in `sam.c`
- declaration only in `sam_internal.h`
- no public header or symbol-map exposure

The primitive is intended for whole-file BAM copy cases where records are not
filtered, decoded, mutated, reheadered record-by-record, indexed on the fly, or
written with a caller-selected compression mode.  It reads and writes the header
normally, writes any residual uncompressed bytes from the current input BGZF
block, flushes the output BGZF stream, then copies the remaining compressed BGZF
blocks directly with `bgzf_raw_read()`/`bgzf_raw_write()`.

The block-copy loop is BGZF-block-aware rather than a blind byte-stream copy. It
checks the BGZF header shape, copies complete non-empty compressed blocks, stops
at the exact canonical 28-byte BGZF EOF marker, rejects trailing bytes after
that marker, and lets `bgzf_close()` write the output EOF marker.  It also
rejects active BGZF thread state, non-compressed/gzip streams, active output
indexing, and non-default output compression level.  This pass-through path is
stricter than ordinary BAM reading: an input missing the terminal BGZF EOF marker
is rejected instead of tolerated with a warning.

The `test/test_view.c` hook is intentionally narrow and opt-in via
`HTS_BAM_RAW_BGZF_COPY=1`.  It currently requires:

- BAM input and BAM output
- compressed BAM output via `-b`
- whole-file copy, no regions
- no `-B` benchmark-discard mode
- no `-N` record limit
- no on-the-fly index
- no explicit `-l` compression level
- no `-@` thread pool

Correctness checks performed so far:

- decoded SAM parity for high-coverage chr20 output
- `test/sam` unit coverage via `test_bam_raw_block_copy()`
- negative tests for uncompressed BAM output, missing BGZF EOF, and trailing
  bytes after BGZF EOF
- negative test for malformed terminal empty BGZF block framing
- build with `make -j4 libhts.a test/test_view test/sam`
- full `make test`: 361 passed, 0 failed, libcurl tests passed

Benchmark command shape:

```sh
./test/test_view -b -p /tmp/out.bam BAM
HTS_BAM_RAW_BGZF_COPY=1 ./test/test_view -b -p /tmp/out.bam BAM
```

Each row is the median of seven runs.  Raw timings:

- `/tmp/htslib_bam_raw_bgzf_copy_safe_timings.tsv`
- `/tmp/htslib_bam_raw_bgzf_copy_eofsafe_timings.tsv` for a shorter rerun
  after exact EOF-marker validation

| input | normal | raw BGZF | speedup |
| --- | ---: | ---: | ---: |
| highcov chr20 10-11Mb | 1.674 | 0.046 | 36.7x |
| exome chr20 | 4.713 | 0.083 | 56.5x |
| ONT ultra-long chr20 10-10.2Mb | 0.311 | 0.024 | 12.9x |

This is the first htslib-level experiment in this branch that comfortably clears
the 2x target.  Its scope is narrower than a general ordered parser: it is a
pass-through primitive for exact BAM copy workloads.  The main remaining design
question is how a real samtools command should prove that no record-level
transformation is active before calling it.  A fully semantic variant would need
to decompress and validate copied BGZF blocks, which would reduce the measured
win; the current version validates BGZF block framing and EOF/trailing-data
structure, but intentionally does not parse each copied BAM record.
