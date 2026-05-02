# BAM Parallel Stream Reader Plan

## Goal

Build a general BAM input pipeline that uses the caller's single `-@` thread
budget and still emits the normal ordered `sam_read1()` stream.  This replaces
the product direction of the indexed ordered-reader prototype, which only wins
against unthreaded `sam_read1()` and loses to existing BGZF `-@`.

## Target Dataflow

```text
sam_read1()
  -> ordered parsed-record queue
     <- parallel BAM parse workers
        <- ordered record-framing layer
           <- parallel BGZF block decompress
              <- compressed BGZF block reader
```

The stream reader splits by physical BAM stream order, not genomic regions, so
it does not require BAI/CSI.  BGZF blocks can be decompressed independently, but
BAM records must be framed from the ordered uncompressed byte stream because a
record can cross a BGZF block boundary.

## Compatibility Rules

- No public API or ABI changes.
- `sam_read1(fp, hdr, b)` remains the caller-facing API.
- `hts_set_threads()` and `HTS_OPT_THREAD_POOL` provide the single thread
  budget.
- Do not enable ordinary BGZF `fp->mt` on the same input while the stream reader
  is active.
- Use bounded queues for compressed blocks, decompressed blocks, parse jobs, and
  parsed result batches.
- Preserve record order, EOF behavior, and error ordering.
- Fall back to existing BGZF/serial paths for unsupported cases.

## Internal State

The eventual implementation should use one tagged BAM input state owned from
`htsFile->state`.

```text
bam_stream_reader_state
  magic/state tag
  borrowed or owned hts_tpool
  thread budget and queue limits
  compressed block window
  decompressed block table
  record-framing carry buffer
  parse job queue
  ordered result table
  current drained batch
  EOF/error state
```

Pools created by `hts_set_threads()` are owned by the BAM stream reader and
destroyed on close.  Pools supplied through `HTS_OPT_THREAD_POOL` are borrowed
and must not be destroyed.

## Staged Implementation

1. **Serial framer.**  Add an internal BAM record framer that reads
   decompressed BGZF bytes in order and emits complete record frames.  It must
   handle split `block_len`, split core, and split payload.
2. **Opt-in stream-reader state.**  Add a BAM-specific state behind an env gate,
   separate from the existing indexed ordered reader.  Activation must happen
   only before records are read.
3. **Parallel BGZF decompress.**  Use the `-@` budget to decompress BGZF blocks
   in parallel, then feed decompressed blocks to the ordered framer.
4. **Parallel parse batches.**  Emit complete-record jobs from the framer and
   parse batches on workers.
5. **Ordered drain.**  Drain parsed batches by sequence number through
   `sam_read1_bam()`, one `bam1_t` at a time.
6. **Fallback and hardening.**  Add non-strict fallback tests, strict activation
   tests, early-close tests, malformed/truncated-input tests, and benchmark
   gates.

## Correctness Cases

- Records crossing BGZF block boundaries.
- Partial `block_len`, partial 32-byte BAM core, and partial variable payload.
- Empty BGZF EOF block, missing EOF marker, truncated BGZF block, and complete
  BGZF block with incomplete BAM record.
- Later worker success must not overtake an earlier worker error.
- Early close must stop workers and free queued batches.
- `tid` and `mtid` validation must use the caller's active header.
- Fallback must leave existing BGZF/serial behavior intact.

## Benchmark Gates

Benchmarks compare serial, current BGZF `-@`, indexed ordered-reader prototype,
and the new stream reader across the local BAM corpus:

- high-coverage short-read slice
- low-coverage short-read slice
- exome chr20 slice
- ONT ultra-long-read slice

Acceptance:

- 1 thread: within 3% of serial.
- `-@ 4` and `-@ 8`: not materially slower than current BGZF `-@`.
- Medium/large read-throughput workloads: target 1.5x or better over BGZF
  `-@ 4`.
- Decoded SAM output must match baseline.
- Memory must remain bounded by documented queue limits.

If the stream reader cannot beat BGZF `-@` broadly, document why and keep it as
an experiment rather than a product feature.

## Current Status

- Added `HTS_BAM_STREAM_READER=1` as a separate opt-in from the indexed ordered
  reader.
- Added a serial stream framer that reads arbitrary decompressed byte chunks,
  assembles complete BAM records, decodes them into the caller's `bam1_t`, and
  emits through `sam_read1_bam()`.
- Added `HTS_BAM_STREAM_READER_CHUNK` as an internal stress knob; tests use a
  7-byte chunk to force split `block_len`, core, and payload handling.
- Thread settings are deferred when stream-reader mode is enabled, so this path
  does not stack on top of BGZF `fp->mt`.
- Added malformed-frame parity tests for partial `block_len`, partial core, and
  partial body, both as the first record and after a valid record.
- Added physical BGZF split-record tests by flushing between pieces of one BAM
  record.
- Added non-strict fallback coverage for the case where stream-reader activation
  declines because BGZF threading is already active.

Verification so far:

```sh
make -j4 libhts.a test/sam
./test/sam
./test/test_view -p /tmp/htslib_stream_base.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
HTS_BAM_STREAM_READER=1 HTS_BAM_STREAM_READER_REQUIRE=1 \
  HTS_BAM_STREAM_READER_CHUNK=7 \
  ./test/test_view -p /tmp/htslib_stream_on.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
cmp /tmp/htslib_stream_base.sam /tmp/htslib_stream_on.sam
git diff --check
```
