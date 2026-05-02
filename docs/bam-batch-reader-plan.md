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
- The batch reader owns decoded BGZF block results while a batch references
  them.
- Split records are copied into a batch-owned carry buffer.
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
- `test/test_view.c` has a benchmark-only batch consumer for BAM input under
  `HTS_BAM_BATCH_READER=1` and `-B`; it counts batches without materializing
  every record as a `bam1_t`.
- `test/sam.c` covers batch parity by rewriting returned raw frames into a
  temporary BAM and comparing ordinary `sam_read1()` hash/count results.  It
  covers serial batch reads, `hts_set_threads()`, borrowed
  `HTS_OPT_THREAD_POOL`, split-record payload handling, and live batch data
  surviving `htsFile` close until `sam_bam_batch_destroy()`.

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
HG00096.exome.chr20.bam          0.28s      0.20s      0.28s      0.17s
HG00096.lowcov.chr20_10-20Mb     0.09s      0.07s      0.08s      0.06s
HG00096.highcov.chr20_10-11Mb    0.07s      0.08s      0.07s      0.06s
HG002.ont_ul.chr20_10-10.2Mb     0.02s      0.03s      0.02s      0.03s
```

The first slice is directionally useful on medium short-read workloads, but it
does not yet clear the product gate broadly.  The tiny ONT slice is too small
to amortize batch setup and header validation, and highcov remains effectively
neutral.  The next step is either a real tool-facing batch consumer or a richer
record-view surface that lets tools do useful work without materializing every
record.

## Benchmark Gate

Compare the batch consumer against serial `sam_read1()`, existing BGZF `-@`,
and the previous stream/parse prototype on the local BAM corpus at `1/4/8`
threads.  The batch path should not be slower than BGZF `-@` at `4/8` threads
and should target at least 1.5x on medium/large read-throughput workloads.

If it cannot meet that gate, keep it prototype-only and document why.
