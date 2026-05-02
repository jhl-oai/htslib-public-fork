# BAM Batch Reader Plan

## Goal

Build the next BAM throughput path as an internal/tool-facing batch reader,
not another transparent `sam_read1()` replacement.  The earlier stream-reader
prototype showed that ordered BGZF decode and framing are feasible, but the
public `sam_read1()` contract still requires one independent `bam1_t` per
record, which prevents a broad win over existing BGZF `-@`.

## Batch Shape

The first batch surface is internal-only and keeps public API/ABI unchanged.

```text
tool/internal consumer
  -> sam_bam_read_batch()
     -> bam_batch_t
        -> contiguous BAM frame bytes
        -> record count
        -> ownership handle for decoded BGZF block or carry buffer
```

Each frame is the on-disk BAM record payload:

```text
int32_t block_len
32-byte core
variable BAM data
```

For the common case, a batch points directly into one ordered decompressed BGZF
block and owns the corresponding pool result until the consumer releases the
batch.  For the uncommon split-record case, the batch owns a carry buffer with
one complete assembled frame.

## Ownership Rules

- A batch is valid until `sam_bam_batch_destroy()` or the next explicit reuse
  by the caller.
- Record frame pointers are views; consumers must not retain them after batch
  release.
- Record accessors derive raw qname, CIGAR, sequence, quality, and aux offsets
  from the descriptor without constructing a `bam1_t`.
- The batch reader owns decoded BGZF block results while a batch references
  them.
- Split records are copied into a batch-owned carry buffer.
- `sam_bam_batch_record_to_bam1()` materializes one descriptor into a normal
  `bam1_t`, so future consumers can batch-scan all records and copy only the
  records they need to emit.
- `sam_read1()` remains the compatibility path and materializes normal
  independent `bam1_t` records.

## Threading

`HTS_BAM_BATCH_READER=1` uses the same deferred BAM thread budget as the stream
reader:

- `hts_set_threads()` creates one reader-owned pool.
- `HTS_OPT_THREAD_POOL` borrows the caller's pool.
- The batch reader does not enable BGZF `fp->mt`.
- Compressed block reads stay on the caller thread; BGZF inflate workers fill
  ordered decoded block results.

## First Implementation Slice

1. Add an internal `bam_batch_t` and `sam_bam_read_batch()` declaration in
   `sam_internal.h`.
2. Reuse the stream-reader activation and ordered BGZF decode path.
3. Reuse the existing frame builder to return one complete batch at a time
   without worker-side `bam1_t` allocation.
4. Add a `test/test_view.c` env-gated benchmark consumer that counts/hashes
   batch frames without materializing every record.
5. Add focused tests for batch count/hash parity, split-record handling,
   thread-budget activation, borrowed thread-pool activation, and fallback.

## Current Status

- Added internal `bam_batch_t`, `sam_bam_read_batch()`, and
  `sam_bam_batch_destroy()` in `sam_internal.h`.
- Added `HTS_BAM_BATCH_READER=1` as an experimental internal gate.
- Batch mode reuses the existing stream-reader activation path and therefore
  consumes the normal `-@` / `HTS_OPT_THREAD_POOL` budget without enabling BGZF
  `fp->mt`.
- `sam_bam_read_batch()` returns one complete ordered raw-frame batch at a
  time.  Batch data is backed either by an owned decoded BGZF block result or
  by an owned carry buffer for split records.
- Each batch includes lightweight `bam_batch_record_t` views with frame/body
  pointers, raw payload length, and decoded core fields.  This lets internal
  consumers inspect per-record metadata without allocating `bam1_t` objects.
- `sam_internal.h` also provides record-view accessors for qname, CIGAR, seq,
  qual, and aux offsets, allowing simple consumers to evaluate filters beyond
  core flag/MAPQ fields while still avoiding full materialization.
- Descriptor construction is fused into the frame scan for batch reads, so the
  batch path no longer scans each raw frame once to count and again to build
  record views.
- Batch frame scanning validates BAM core layout before exposing record views,
  and serial batches detach their backing data so live batches survive
  `htsFile` close just like threaded block-backed batches.
- `test/test_view.c` has a benchmark-only batch consumer for BAM input under
  `HTS_BAM_BATCH_READER=1` and `-B`; it walks record views without
  materializing every record as a `bam1_t`.
- `test/sam.c` covers batch parity by rewriting returned raw frames into a
  temporary BAM and comparing ordinary `sam_read1()` hash/count results.  It
  covers serial batch reads, `hts_set_threads()`, borrowed
  `HTS_OPT_THREAD_POOL`, split-record payload handling, malformed core
  rejection, and live serial/threaded batch data surviving `htsFile` close
  until `sam_bam_batch_destroy()`.

Verification so far:

```sh
make -j4 libhts.a test/test_view test/test_bgzf test/sam
./test/sam
./test/test_bgzf test/bgziptest.txt
make test
git diff --check
```

Repeated median-of-five benchmark notes from `/tmp/bam_batch_bench.tsv`:

```text
test_view -B, lower is better

Input                         BGZF -@4  Batch -@4  BGZF -@8  Batch -@8
HG00096.exome.chr20.bam          0.28s      0.20s      0.28s      0.16s
HG00096.lowcov.chr20_10-20Mb     0.09s      0.06s      0.08s      0.06s
HG00096.highcov.chr20_10-11Mb    0.07s      0.07s      0.06s      0.06s
HG002.ont_ul.chr20_10-10.2Mb     0.03s      0.03s      0.02s      0.02s
```

The current slice is faster on the medium short-read workloads and neutral on
the tiny highcov/ONT slices.  It clears the local "not slower than BGZF `-@`"
benchmark sweep, but only the larger exome slice reaches roughly 1.5x.  Fusing
header TID validation into descriptor construction removes the remaining
validation pass over views, but did not materially move this local benchmark.

## Samtools Consumer Experiment

A sibling `samtools` branch, `feature/bam-batch-reader-consumer`, adds an
opt-in `samtools view -c` batch consumer compiled with
`HTS_BAM_BATCH_READER_CONSUMER` and run under `HTS_BAM_BATCH_READER=1`.
It only handles streaming BAM count mode with simple flag/MAPQ filters and
minimum query-length filtering from raw CIGAR views; complex filters and output
modes fall back to the existing `sam_read1()` path.

Median-of-five local timings from `/tmp/samtools_batch_view_bench.tsv`:

```text
samtools view -c, lower is better

Input                         BGZF -@1  Batch -@1  BGZF -@4  Batch -@4  BGZF -@8  Batch -@8
HG00096.exome.chr20.bam          0.60s      0.62s      0.28s      0.21s      0.28s      0.18s
HG00096.lowcov.chr20_10-20Mb     0.16s      0.17s      0.07s      0.06s      0.08s      0.05s
HG00096.highcov.chr20_10-11Mb    0.22s      0.22s      0.06s      0.06s      0.05s      0.05s
HG002.ont_ul.chr20_10-10.2Mb     0.05s      0.05s      0.02s      0.02s      0.01s      0.02s
```

The end-to-end count consumer is neutral at one thread, faster on the larger
short-read slices with `-@`, and still dominated by timing noise on the tiny
ONT slice.  The plain count path now validates long-CIGAR `CG` candidates and
mapped raw-CIGAR query lengths to match the normal decode error contract.

Median-of-five local timings for the raw-CIGAR `-m 75` count fast path from
`/tmp/samtools_batch_minqlen_bench.tsv`:

```text
samtools view -c -m 75, lower is better

Input                         BGZF -@4  Batch -@4  BGZF -@8  Batch -@8
HG00096.exome.chr20.bam          0.28s      0.20s      0.28s      0.19s
HG00096.lowcov.chr20_10-20Mb     0.07s      0.05s      0.07s      0.06s
HG00096.highcov.chr20_10-11Mb    0.06s      0.06s      0.05s      0.05s
HG002.ont_ul.chr20_10-10.2Mb     0.01s      0.02s      0.01s      0.02s
```

This extends the fast path beyond core-only predicates while keeping the same
short-read benefit.  The ONT slice is too small for stable timing and remains a
prototype-only datapoint.  `-B` remains excluded from this fast path because
normal `samtools view` applies `bam_remove_B()` before query-length and flag
filtering.  Long-CIGAR `CG` candidates are selectively materialized for
correct query-length semantics and malformed-input parity; ordinary records
validate raw CIGAR query length before filtering.

## Benchmark Gate

Compare the batch consumer against serial `sam_read1()`, existing BGZF `-@`,
and the previous stream/parse prototype on the local BAM corpus at `1/4/8`
threads.  The batch path should not be slower than BGZF `-@` at `4/8` threads
and should target at least 1.5x on medium/large read-throughput workloads.

If it cannot meet that gate, keep it prototype-only and document why.
