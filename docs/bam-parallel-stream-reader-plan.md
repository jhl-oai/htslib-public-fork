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
- Added a serial stream framer that consumes explicit decompressed BGZF block
  bytes, assembles complete BAM records, decodes them into the caller's
  `bam1_t`, and emits through `sam_read1_bam()`.
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
- Added ordered parallel BGZF block decompression behind the same
  `HTS_BAM_STREAM_READER=1` opt-in.  `hts_set_threads()` creates one
  reader-owned pool; `HTS_OPT_THREAD_POOL` borrows the caller's pool; neither
  path creates BGZF `fp->mt`.
- Compressed BGZF reads remain on the caller thread.  Workers only inflate
  independent BGZF blocks into per-job buffers; the ordered drain is the only
  code that updates visible `BGZF` error/index/position state.
- Added `HTS_BAM_STREAM_PARSE=1` as a separate experimental gate for parse
  batches.  The current experimental parser avoids copying raw frames that are
  fully contained within one decompressed BGZF block; only records spanning
  BGZF block boundaries use a carry buffer.  It still allocates decoded
  `bam1_t` records in workers before moving them into the caller's `bam1_t`.

Verification so far:

```sh
make -j4 libhts.a test/sam
./test/sam
./test/test_bgzf test/bgziptest.txt
./test/test_view -p /tmp/htslib_stream_base.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
HTS_BAM_STREAM_READER=1 HTS_BAM_STREAM_READER_REQUIRE=1 \
  HTS_BAM_STREAM_READER_CHUNK=7 \
  ./test/test_view -p /tmp/htslib_stream_on.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
cmp /tmp/htslib_stream_base.sam /tmp/htslib_stream_on.sam
HTS_BAM_STREAM_READER=1 HTS_BAM_STREAM_READER_REQUIRE=1 \
  ./test/test_view -@4 -p /tmp/htslib_stream_threaded.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
cmp /tmp/htslib_stream_base.sam /tmp/htslib_stream_threaded.sam
HTS_BAM_STREAM_READER=1 HTS_BAM_STREAM_READER_REQUIRE=1 \
  HTS_BAM_STREAM_PARSE=1 \
  ./test/test_view -@4 -p /tmp/htslib_stream_parse.sam DATA/HG00096.highcov.chr20_10-11Mb.bam
cmp /tmp/htslib_stream_base.sam /tmp/htslib_stream_parse.sam
git diff --check
```

Smoke benchmark notes on `HG00096.exome.chr20.bam` with `test_view -B`
(single-run, noisy):

- After zero-copy block references: BGZF `-@4` 0.27s, stream `-@4` 0.33s,
  parse `-@4` 0.27s; BGZF `-@8` 0.28s, stream `-@8` 0.31s, parse `-@8`
  0.31s.
- Parallel BGZF inflate alone is roughly at parity in this slice, but not a
  product-level win yet.
- The first raw-copy parse-batch implementation was slower and was pruned.
- A descriptor/materialize experiment that avoided worker-side `bam1_t`
  allocation was also slower in smoke tests and was pruned.
- A later local experiment added recycled `bam1_t::data` buffers, grouped
  multiple BGZF-block segments into one parse job, and cached libdeflate
  decompressors for the custom decode workers.  It preserved parity but still
  lost to current BGZF `-@` on repeated corpus sweeps, so it was pruned rather
  than committed as product code.
- The remaining experimental parse path is correctness-useful but still not a
  product-level win because worker-side `bam1_t` allocation and queue overhead
  remain significant.

Next likely useful parse design:

- Replace per-record worker allocations with a batch arena or reusable
  `bam1_t` pool only if the design preserves normal `bam1_t` lifetime after
  file close.  Reader-owned arenas are unsafe for the existing `sam_read1()`
  contract unless records are copied before returning or arena ownership is
  transferred to the caller.
- Keep frame references into ordered decompressed block buffers plus a small
  carry buffer for split records.
- Keep `HTS_BAM_STREAM_PARSE=1` as a correctness scaffold until it beats BGZF
  `-@` on repeated median benchmarks.
