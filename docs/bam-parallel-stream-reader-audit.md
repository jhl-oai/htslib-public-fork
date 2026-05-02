# BAM Parallel Stream Reader Completion Audit

## Objective

Build an HTSlib-level BAM input pipeline that uses one caller-provided thread
budget, preserves the ordered `sam_read1(fp, hdr, b)` stream, avoids public
API/ABI changes, and is not slower than existing BGZF `-@` at `-@4`/`-@8`
while targeting a 1.5x+ read-throughput win.

## Checklist

| Requirement | Evidence | Status |
|---|---|---|
| No public API/ABI break | Changes are in `sam.c`, `bgzf.c`, `bgzf_internal.h`, and tests/docs. No public `htslib/sam.h` API change or export-map change. | Pass |
| `sam_read1(fp, hdr, b)` unchanged | Stream reader activates internally from `sam_read1_bam()` and returns through the existing `bam1_t *b`. | Pass |
| `hts_set_threads()` / `HTS_OPT_THREAD_POOL` provide one budget | BAM threading is deferred when `HTS_BAM_STREAM_READER=1`; owned pools are created for `hts_set_threads()`, borrowed pools are used for `HTS_OPT_THREAD_POOL`. | Pass |
| Do not stack on BGZF `fp->mt` | Stream open rejects existing `fp->fp.bgzf->mt`; deferred config is consumed by the stream reader before ordinary BGZF threading is enabled. | Pass |
| No BAI/CSI required | Stream reader consumes physical BGZF block order and does not use an index. | Pass |
| Preserve record order | `hts_tpool_process` ordered result drain is used for BGZF decode and parse batches; `sam_read1_bam()` drains one record at a time. | Pass |
| Preserve error/EOF behavior | Tests cover partial `block_len`, partial core, partial body, valid-then-partial records, strict activation, and fallback. | Mostly pass |
| Handle records spanning BGZF blocks | Tests write records split across physical BGZF flush boundaries at block_len/core/body positions. | Pass |
| Bounded queues | Decode and parse queues are bounded by `qsize`; parse batches are capped by record/byte limits. | Pass |
| Clean fallback | Non-strict fallback activates existing BGZF threading if stream open declines before installation. | Pass |
| Full test gates | `make -j4 libhts.a test/test_view test/test_bgzf test/sam`, `./test/sam`, `./test/test_bgzf test/bgziptest.txt`, and `make test` have passed after the implementation commits. | Pass |
| Decoded SAM parity | Baseline vs stream and parse outputs compare cleanly on highcov and exome slices. | Pass |
| Not slower than BGZF `-@` at `-@4`/`-@8` | Repeated `test_view -B` sweeps show the stream/parse paths still lose or tie noisily; highcov and small ONT slices are worse. | Fail |
| Target 1.5x+ on medium/large read-throughput workloads | No benchmark reaches this. | Fail |

## Current Conclusion

The branch has a correct opt-in scaffold for ordered stream reading and
parallel BGZF/parse experiments, but it does not meet the product performance
gate.  It should remain experimental.

The main blocker is the existing `sam_read1()` contract: callers receive a
normal `bam1_t` whose `data` remains valid independently of the file reader.
Reader-owned arenas or block-backed records would reduce copies and allocation
cost, but they would change lifetime semantics unless data is copied before
returning or ownership is transferred to the caller.  The attempted low-risk
variants did not beat current BGZF `-@` broadly.

## Pruned Experiments

- Raw-copy parse batches: correct but slower due frame copies and per-record
  allocations.
- Descriptor/materialize parse batches: avoided worker-side `bam1_t`
  allocation but was slower in smoke tests.
- Local uncommitted experiment with recycled `bam1_t::data`, grouped parse
  segments, and libdeflate decode-cache reuse: parity preserved, but repeated
  corpus sweeps still lost to current BGZF `-@`.

## Product Path Implication

A broad Sambamba-style win likely needs an internal batch-oriented read surface
or tool-layer batch consumption so parsed records can be retained in arenas or
block-backed storage without materializing one independent `bam1_t` per
`sam_read1()` call.  That can be added without breaking the public API if
`sam_read1()` remains a compatibility wrapper, but it is a larger design than
the transparent reader replacement attempted here.
