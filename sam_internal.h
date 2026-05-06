/*  sam_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2019-2020, 2023-2024 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_SAM_INTERNAL_H
#define HTSLIB_SAM_INTERNAL_H

#include <errno.h>
#include <stdint.h>

#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

// Used internally in the SAM format multi-threading.
int sam_state_destroy(samFile *fp);
int sam_bam_state_destroy(samFile *fp);
int sam_set_thread_pool(htsFile *fp, htsThreadPool *p);
int sam_set_threads(htsFile *fp, int nthreads);
int sam_bam_raw_copy_blocks(htsFile *in, htsFile *out);

enum {
    // Internal batch-owned frame/body storage.  Callers must still treat public
    // record pointers as read-only; this flag lets HTSlib internals know the
    // backing buffer may be temporarily patched and restored before return.
    BAM_BATCH_RECORD_F_OWNED = 1u,
    BAM_BATCH_RECORD_F_RAW_WRITE_SAFE = 2u,
    BAM_BATCH_RECORD_F_NEEDS_MATERIALIZE = 4u,
    BAM_BATCH_RECORD_F_ENDPOS_VALID = 8u,
    BAM_BATCH_RECORD_F_RAW_LAYOUT_VALID = 16u,
    BAM_BATCH_SEGMENT_F_OWNED = 1u
};

typedef struct bam_batch_record_t {
    const uint8_t *frame;
    size_t frame_len;
    const uint8_t *body;
    uint32_t raw_l_data;
    uint64_t voff_beg;
    uint64_t voff_end;
    hts_pos_t endpos;
    uint32_t flags;
    // Core values are decoded from the raw BAM frame.  l_qname remains the
    // on-wire length, so body + core.l_qname points at raw CIGAR data.
    bam1_core_t core;
} bam_batch_record_t;

typedef struct bam_batch_segment_t {
    const uint8_t *data;
    size_t len;
    uint32_t flags;
} bam_batch_segment_t;

typedef struct bam_batch_t {
    const uint8_t *data;
    size_t len;
    int n_records;
    bam_batch_record_t *records;
    int n_segments;
    const bam_batch_segment_t *segments;
    void *impl;
} bam_batch_t;

typedef struct sam_bam_voff_span_t {
    uint64_t beg;
    uint64_t end;
} sam_bam_voff_span_t;

typedef struct sam_bam_record_view_t {
    const bam1_core_t *core;
    const uint8_t *body;
    uint32_t l_data;
    const bam_batch_record_t *batch_record;
    const bam1_t *bam_record;
} sam_bam_record_view_t;

typedef int (*sam_bam_aux_filter_f)(const char tag[2],
                                    const uint8_t *value,
                                    void *data);

enum {
    SAM_BAM_FILTER_PLAN_MAX_PREDICATES = 16
};

typedef enum sam_bam_filter_plan_class_t {
    SAM_BAM_FILTER_PLAN_UNUSABLE = 0,
    SAM_BAM_FILTER_PLAN_RAW_VIEW_SAFE = 1,
    SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE = 2,
    SAM_BAM_FILTER_PLAN_MATERIALIZE_REQUIRED = 3
} sam_bam_filter_plan_class_t;

typedef struct sam_bam_filter_plan_pred_t {
    int kind;
    int field;
    int cmp;
    int negate;
    char tag[2];
    double number;
    char *string;
} sam_bam_filter_plan_pred_t;

typedef struct sam_bam_filter_plan_t {
    sam_bam_filter_plan_class_t class_;
    int n_predicates;
    sam_bam_filter_plan_pred_t predicates[SAM_BAM_FILTER_PLAN_MAX_PREDICATES];
} sam_bam_filter_plan_t;

static inline const uint8_t *sam_bam_batch_record_qname(
        const bam_batch_record_t *record)
{
    return record->body;
}

static inline const uint8_t *sam_bam_batch_record_cigar(
        const bam_batch_record_t *record)
{
    return record->body + record->core.l_qname;
}

static inline const uint8_t *sam_bam_batch_record_seq(
        const bam_batch_record_t *record)
{
    return sam_bam_batch_record_cigar(record) + ((size_t)record->core.n_cigar << 2);
}

static inline const uint8_t *sam_bam_batch_record_qual(
        const bam_batch_record_t *record)
{
    return sam_bam_batch_record_seq(record) + (((size_t)record->core.l_qseq + 1) >> 1);
}

static inline const uint8_t *sam_bam_batch_record_aux(
        const bam_batch_record_t *record)
{
    return sam_bam_batch_record_qual(record) + record->core.l_qseq;
}

static inline size_t sam_bam_batch_record_aux_len(
        const bam_batch_record_t *record)
{
    const uint8_t *aux = sam_bam_batch_record_aux(record);
    return record->body + record->raw_l_data >= aux
           ? (size_t)(record->body + record->raw_l_data - aux) : 0;
}

static inline void sam_bam_record_view_from_batch(
        sam_bam_record_view_t *view, const bam_batch_record_t *record)
{
    view->core = &record->core;
    view->body = record->body;
    view->l_data = record->raw_l_data;
    view->batch_record = record;
    view->bam_record = NULL;
}

static inline void sam_bam_record_view_from_bam(
        sam_bam_record_view_t *view, const bam1_t *bam)
{
    view->core = &bam->core;
    view->body = bam->data;
    view->l_data = bam->l_data;
    view->batch_record = NULL;
    view->bam_record = bam;
}

static inline const uint8_t *sam_bam_record_view_qname(
        const sam_bam_record_view_t *view)
{
    return view->body;
}

static inline const uint8_t *sam_bam_record_view_cigar(
        const sam_bam_record_view_t *view)
{
    return view->body + view->core->l_qname;
}

static inline const uint8_t *sam_bam_record_view_seq(
        const sam_bam_record_view_t *view)
{
    return sam_bam_record_view_cigar(view) + ((size_t)view->core->n_cigar << 2);
}

static inline const uint8_t *sam_bam_record_view_qual(
        const sam_bam_record_view_t *view)
{
    return sam_bam_record_view_seq(view) + (((size_t)view->core->l_qseq + 1) >> 1);
}

static inline const uint8_t *sam_bam_record_view_aux(
        const sam_bam_record_view_t *view)
{
    return sam_bam_record_view_qual(view) + view->core->l_qseq;
}

static inline size_t sam_bam_record_view_aux_len(
        const sam_bam_record_view_t *view)
{
    const uint8_t *aux = sam_bam_record_view_aux(view);
    return view->body + view->l_data >= aux
           ? (size_t)(view->body + view->l_data - aux) : 0;
}

static inline int sam_bam_record_view_aux_type_size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline const uint8_t *sam_bam_record_view_aux_skip(
        const uint8_t *s, const uint8_t *end)
{
    int size;
    uint32_t n;

    if (s >= end)
        return end;
    size = sam_bam_record_view_aux_type_size(*s++);
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s)
            s++;
        return s < end ? s + 1 : NULL;
    case 'B':
        if (end - s < 5)
            return NULL;
        size = sam_bam_record_view_aux_type_size(*s++);
        n = le_to_u32(s);
        s += 4;
        if (size == 0 ||
            (size_t)n > (size_t)(end - s) / (size_t)size)
            return NULL;
        return s + (size_t)size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size)
            return NULL;
        return s + size;
    }
}

static inline const uint8_t *sam_bam_record_view_aux_get_status(
        const sam_bam_record_view_t *view, const char tag[2], int *status)
{
    const uint8_t *s = sam_bam_record_view_aux(view);
    const uint8_t *end = view->body + view->l_data;

    if (status)
        *status = ENOENT;
    while (end - s > 2) {
        const uint8_t *value = s + 2;
        const uint8_t *next;

        next = sam_bam_record_view_aux_skip(value, end);
        if (!next) {
            if (status)
                *status = EINVAL;
            return NULL;
        }
        if (s[0] == tag[0] && s[1] == tag[1]) {
            if ((*value == 'Z' || *value == 'H') && next[-1] != '\0') {
                if (status)
                    *status = EINVAL;
                return NULL;
            }
            if (status)
                *status = 0;
            return value;
        }
        s = next;
    }
    return NULL;
}

static inline const uint8_t *sam_bam_record_view_aux_get(
        const sam_bam_record_view_t *view, const char tag[2])
{
    int status = 0;
    const uint8_t *s = sam_bam_record_view_aux_get_status(view, tag, &status);

    if (!s)
        errno = status;
    return s;
}

int sam_bam_read_batch(htsFile *fp, sam_hdr_t *h, bam_batch_t *batch);
int sam_bam_read_batch_voff(htsFile *fp, sam_hdr_t *h, bam_batch_t *batch);
int sam_bam_read_batch_count(htsFile *fp, sam_hdr_t *h, int *n_records);
int sam_bam_itr_next_batch(htsFile *fp, hts_itr_t *iter, sam_hdr_t *h,
                           bam_batch_t *batch);
int sam_bam_batch_seek(htsFile *fp, uint64_t voff);
void sam_bam_batch_destroy(bam_batch_t *batch);
typedef int (*sam_region_overlap_f)(void *data, sam_hdr_t *h, int tid,
                                    hts_pos_t beg, hts_pos_t end);
hts_reglist_t *sam_reglist_dup_merged(const hts_reglist_t *reglist,
                                      int count, int *out_count);
int sam_itr_next_filtered(htsFile *fp, hts_itr_t *iter, sam_hdr_t *h,
                          bam1_t *record, sam_region_overlap_f filter,
                          void *filter_data);
// Returns a record-list batch only.  Filtering can break full-batch raw storage
// contiguity, so callers should use records/range writers or materialize.
int sam_bam_itr_next_batch_filtered(htsFile *fp, hts_itr_t *iter,
                                    sam_hdr_t *h, bam_batch_t *batch,
                                    sam_region_overlap_f filter,
                                    void *filter_data, bam1_t *scratch);
typedef struct sam_bam_voff_span_stats_t {
    uint64_t selected_uncomp;
    uint64_t full_uncomp;
    uint64_t partial_uncomp;
    uint64_t full_blocks;
    uint64_t partial_blocks;
} sam_bam_voff_span_stats_t;
int sam_bam_raw_copy_voff_spans_stats(htsFile *in,
                                      const sam_bam_voff_span_t *spans,
                                      int n_spans,
                                      sam_bam_voff_span_stats_t *stats);
// Raw byte primitive for BAM record spans collected from sam_bam_read_batch_voff().
// The caller must pass sorted, non-overlapping spans that begin/end on raw-write
// safe BAM record boundaries.  Writing requires default-compressed BAM output;
// stats-only preflight does not inspect output state.
int sam_bam_raw_copy_voff_spans(htsFile *in, htsFile *out,
                                const sam_bam_voff_span_t *spans,
                                int n_spans);
// Materialize a record view into bam1_t.  Mutation-heavy callers should use
// this and then apply the existing bam1_t mutation APIs, including aux append
// and aux update helpers, so public bam1_t semantics stay centralized.
int sam_bam_batch_record_to_bam1(const bam_batch_record_t *record, bam1_t *bam);
int sam_bam_batch_record_cg_candidate(const bam_batch_record_t *record);
// Validates the decode-sensitive checks and fixups performed by bam_read1().
// If materialized is set to 1, scratch holds the bam_read1()-equivalent record
// and raw frame output would bypass a decode-time fixup.
int sam_bam_batch_record_validate_decode(bam_batch_record_t *record,
                                         bam1_t *scratch, int *materialized);
// Returns 0 when the raw BAM frame can be written unchanged, 1 when the
// record must be materialized to preserve bam_read1() decode fixups, and -1
// when the raw record layout or decode-sensitive state is invalid.
int sam_bam_batch_record_raw_write_status(bam_batch_record_t *record);
// Raw writers preserve the input frame only when decode-sensitive fixups are
// unnecessary.  They return -2 for records that need materialization, leaving
// callers responsible for materializing and writing through the ordinary path.
int sam_bam_batch_record_write1(htsFile *fp, const sam_hdr_t *h,
                                bam_batch_record_t *record);
int sam_bam_batch_record_write1_aux_filtered(htsFile *fp, const sam_hdr_t *h,
                                             bam_batch_record_t *record,
                                             sam_bam_aux_filter_f filter,
                                             void *filter_data);
int sam_bam_batch_write1(htsFile *fp, const sam_hdr_t *h,
                         const bam_batch_t *batch);
int sam_bam_batch_write1_range(htsFile *fp, const sam_hdr_t *h,
                               const bam_batch_t *batch,
                               int beg, int end);
int sam_bam_batch_write1_ranges(htsFile *fp, const sam_hdr_t *h,
                                const bam_batch_t *batch,
                                const int *ranges, int n_ranges);
// Internal SAM text writer for unchanged BAM batch records.  It formats from
// batch record views and returns -2 when callers must materialize and fall back
// to sam_write1().
int sam_bam_batch_write_sam_ranges(htsFile *fp, const sam_hdr_t *h,
                                   const bam_batch_t *batch,
                                   const int *ranges, int n_ranges);
int sam_bam_batch_write_sam_ranges_steal(htsFile *fp, const sam_hdr_t *h,
                                         bam_batch_t *batch,
                                         const int *ranges, int n_ranges);
int sam_bam_batch_write1_range_flags(htsFile *fp, const sam_hdr_t *h,
                                     const bam_batch_t *batch,
                                     int beg, int end,
                                     uint16_t set_flags,
                                     uint16_t clear_flags);
int sam_bam_batch_record_endpos(const bam_batch_record_t *record,
                                bam1_t *scratch, hts_pos_t *endpos);
// Computes CIGAR query length from a batch-lifetime record view.  If query_len
// is NULL, this only validates decode-sensitive CIGAR state.  Long-CIGAR CG
// candidates are materialized into scratch, which must be non-NULL for those
// records and may be overwritten.
int sam_bam_batch_record_query_len(const bam_batch_record_t *record,
                                   bam1_t *scratch, int include_hard_clip,
                                   hts_pos_t *query_len);
int sam_bam_batch_record_passes_filter(const sam_hdr_t *h,
                                       bam_batch_record_t *record,
                                       bam1_t *scratch,
                                       struct hts_filter_t *filt,
                                       int *materialized);
int sam_bam_filter_plan_init(sam_bam_filter_plan_t *plan, const char *expr);
void sam_bam_filter_plan_destroy(sam_bam_filter_plan_t *plan);
int sam_bam_filter_plan_is_usable(const sam_bam_filter_plan_t *plan);
sam_bam_filter_plan_class_t sam_bam_filter_plan_class(
        const sam_bam_filter_plan_t *plan);
int sam_bam_batch_record_passes_filter_plan(const sam_hdr_t *h,
                                            bam_batch_record_t *record,
                                            bam1_t *scratch,
                                            const sam_bam_filter_plan_t *plan,
                                            int *materialized);
int sam_bam_prepare_batch_reader(htsFile *fp);

// Fastq state
int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...);
void fastq_state_destroy(samFile *fp);

// bam1_t data (re)allocation
int sam_realloc_bam_data(bam1_t *b, size_t desired);

static inline int realloc_bam_data(bam1_t *b, size_t desired)
{
    if (desired <= b->m_data) return 0;
    return sam_realloc_bam_data(b, desired);
}

static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    size_t new_len = (size_t) b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < bytes) { // Too big or overflow
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return sam_realloc_bam_data(b, new_len);
}

/*
 * Convert a nibble encoded BAM sequence to a string of bases.
 *
 * We do this 2 bp at a time for speed. Equiv to:
 *
 * for (i = 0; i < len; i++)
 *    seq[i] = seq_nt16_str[bam_seqi(nib, i)];
 */
static inline void nibble2base_default(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len/2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i*2], &code2base[(size_t)nib[i]*2], 2);

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

#if defined HAVE_ATTRIBUTE_CONSTRUCTOR && \
    ((defined __x86_64__ && \
      defined HAVE_X86INTRIN_H && HAVE_X86INTRIN_H && \
      defined HAVE_ATTRIBUTE_TARGET_SSSE3 && \
      defined HAVE_BUILTIN_CPU_SUPPORT_SSSE3) || \
     (defined __ARM_NEON))
#define BUILDING_SIMD_NIBBLE2BASE
#endif

static inline void nibble2base(uint8_t *nib, char *seq, int len) {
#ifdef BUILDING_SIMD_NIBBLE2BASE
    extern void (*htslib_nibble2base)(uint8_t *nib, char *seq, int len);
    htslib_nibble2base(nib, seq, len);
#else
    nibble2base_default(nib, seq, len);
#endif
}

#ifdef __cplusplus
}
#endif

#endif
