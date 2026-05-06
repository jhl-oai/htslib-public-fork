/*  sam.c -- SAM and BAM file I/O and manipulation.

    Copyright (C) 2008-2010, 2012-2025 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include <signal.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <unistd.h>
#include <regex.h>
#include <ctype.h>

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
#include "fuzz_settings.h"
#endif

// Suppress deprecation message for cigar_tab, which we initialise
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hts_internal.h"
#include "bgzf_internal.h"
#include "sam_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/hts_expr.h"
#include "header.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_SET_INIT_INT(tag)

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif

/**********************
 *** BAM header I/O ***
 **********************/

HTSLIB_EXPORT
const int8_t bam_cigar_table[256] = {
    // 0 .. 47
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 48 .. 63  (including =)
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, BAM_CEQUAL, -1, -1,

    // 64 .. 79  (including MIDNHB)
    -1, -1, BAM_CBACK, -1,  BAM_CDEL, -1, -1, -1,
        BAM_CHARD_CLIP, BAM_CINS, -1, -1,  -1, BAM_CMATCH, BAM_CREF_SKIP, -1,

    // 80 .. 95  (including SPX)
    BAM_CPAD, -1, -1, BAM_CSOFT_CLIP,  -1, -1, -1, -1,
        BAM_CDIFF, -1, -1, -1,  -1, -1, -1, -1,

    // 96 .. 127
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 128 .. 255
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
};

sam_hdr_t *sam_hdr_init(void)
{
    sam_hdr_t *bh = (sam_hdr_t*)calloc(1, sizeof(sam_hdr_t));
    if (bh == NULL) return NULL;

    bh->cigar_tab = bam_cigar_table;
    return bh;
}

void sam_hdr_destroy(sam_hdr_t *bh)
{
    int32_t i;

    if (bh == NULL) return;

    if (bh->ref_count > 0) {
        --bh->ref_count;
        return;
    }

    if (bh->target_name) {
        for (i = 0; i < bh->n_targets; ++i)
            free(bh->target_name[i]);
        free(bh->target_name);
        free(bh->target_len);
    }
    free(bh->text);
    if (bh->hrecs)
        sam_hrecs_free(bh->hrecs);
    if (bh->sdict)
        kh_destroy(s2i, (khash_t(s2i) *) bh->sdict);
    free(bh);
}

// Copy the sam_hdr_t::sdict hash, used to store the real lengths of long
// references before sam_hdr_t::hrecs is populated
int sam_hdr_dup_sdict(const sam_hdr_t *h0, sam_hdr_t *h)
{
    const khash_t(s2i) *src_long_refs = (khash_t(s2i) *) h0->sdict;
    khash_t(s2i) *dest_long_refs = kh_init(s2i);
    int i;
    if (!dest_long_refs) return -1;

    for (i = 0; i < h->n_targets; i++) {
        int ret;
        khiter_t ksrc, kdest;
        if (h->target_len[i] < UINT32_MAX) continue;
        ksrc = kh_get(s2i, src_long_refs, h->target_name[i]);
        if (ksrc == kh_end(src_long_refs)) continue;
        kdest = kh_put(s2i, dest_long_refs, h->target_name[i], &ret);
        if (ret < 0) {
            kh_destroy(s2i, dest_long_refs);
            return -1;
        }
        kh_val(dest_long_refs, kdest) = kh_val(src_long_refs, ksrc);
    }

    h->sdict = dest_long_refs;
    return 0;
}

sam_hdr_t *sam_hdr_dup(const sam_hdr_t *h0)
{
    if (h0 == NULL) return NULL;
    sam_hdr_t *h;
    if ((h = sam_hdr_init()) == NULL) return NULL;
    // copy the simple data
    h->n_targets = 0;
    h->ignore_sam_err = h0->ignore_sam_err;
    h->l_text = 0;

    // Then the pointery stuff

    if (!h0->hrecs) {
        h->target_len = (uint32_t*)calloc(h0->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto fail;
        h->target_name = (char**)calloc(h0->n_targets, sizeof(char*));
        if (!h->target_name) goto fail;

        int i;
        for (i = 0; i < h0->n_targets; ++i) {
            h->target_len[i] = h0->target_len[i];
            h->target_name[i] = strdup(h0->target_name[i]);
            if (!h->target_name[i]) break;
        }
        h->n_targets = i;
        if (i < h0->n_targets) goto fail;

        if (h0->sdict) {
            if (sam_hdr_dup_sdict(h0, h) < 0) goto fail;
        }
    }

    if (h0->hrecs) {
        kstring_t tmp = { 0, 0, NULL };
        if (sam_hrecs_rebuild_text(h0->hrecs, &tmp) != 0) {
            free(ks_release(&tmp));
            goto fail;
        }

        h->l_text = tmp.l;
        h->text   = ks_release(&tmp);

        if (sam_hdr_update_target_arrays(h, h0->hrecs, 0) != 0)
            goto fail;
    } else {
        h->l_text = h0->text ? h0->l_text : 0;
        h->text = malloc(h->l_text + 1);
        if (!h->text) goto fail;
        if (h0->text)
            memcpy(h->text, h0->text, h->l_text);
        h->text[h->l_text] = '\0';
    }

    return h;

 fail:
    sam_hdr_destroy(h);
    return NULL;
}

sam_hdr_t *bam_hdr_read(BGZF *fp)
{
    sam_hdr_t *h;
    uint8_t buf[4];
    int magic_len, has_EOF;
    int32_t i, name_len, num_names = 0;
    size_t bufsize;
    ssize_t bytes;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) {
        perror("[W::bam_hdr_read] bgzf_check_EOF");
    } else if (has_EOF == 0) {
        hts_log_warning("EOF marker is absent. The input is probably truncated");
    }
    // read "BAM1"
    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || memcmp(buf, "BAM\1", 4)) {
        hts_log_error("Invalid BAM binary header");
        return 0;
    }
    h = sam_hdr_init();
    if (!h) goto nomem;

    // read plain text and the number of reference sequences
    bytes = bgzf_read(fp, buf, 4);
    if (bytes != 4) goto read_err;
    h->l_text = le_to_u32(buf);

    bufsize = h->l_text + 1;
    if (bufsize < h->l_text) goto nomem; // so large that adding 1 overflowed
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (bufsize > FUZZ_ALLOC_LIMIT) goto nomem;
#endif
    h->text = (char*)malloc(bufsize);
    if (!h->text) goto nomem;
    h->text[h->l_text] = 0; // make sure it is NULL terminated
    bytes = bgzf_read(fp, h->text, h->l_text);
    if (bytes != h->l_text) goto read_err;

    bytes = bgzf_read(fp, &h->n_targets, 4);
    if (bytes != 4) goto read_err;
    if (fp->is_be) ed_swap_4p(&h->n_targets);

    if (h->n_targets < 0) goto invalid;

    // read reference sequence names and lengths
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (h->n_targets > (FUZZ_ALLOC_LIMIT - bufsize)/(sizeof(char*)+sizeof(uint32_t)))
        goto nomem;
#endif
    if (h->n_targets > 0) {
        h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
        if (!h->target_name) goto nomem;
        h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto nomem;
    }
    else {
        h->target_name = NULL;
        h->target_len = NULL;
    }

    for (i = 0; i != h->n_targets; ++i) {
        bytes = bgzf_read(fp, &name_len, 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&name_len);
        if (name_len <= 0) goto invalid;

        h->target_name[i] = (char*)malloc(name_len);
        if (!h->target_name[i]) goto nomem;
        num_names++;

        bytes = bgzf_read(fp, h->target_name[i], name_len);
        if (bytes != name_len) goto read_err;

        if (h->target_name[i][name_len - 1] != '\0') {
            /* Fix missing NUL-termination.  Is this being too nice?
               We could alternatively bail out with an error. */
            char *new_name;
            if (name_len == INT32_MAX) goto invalid;
            new_name = realloc(h->target_name[i], name_len + 1);
            if (new_name == NULL) goto nomem;
            h->target_name[i] = new_name;
            h->target_name[i][name_len] = '\0';
        }

        bytes = bgzf_read(fp, &h->target_len[i], 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&h->target_len[i]);
    }
    return h;

 nomem:
    hts_log_error("Out of memory");
    goto clean;

 read_err:
    if (bytes < 0) {
        hts_log_error("Error reading BGZF stream");
    } else {
        hts_log_error("Truncated BAM header");
    }
    goto clean;

 invalid:
    hts_log_error("Invalid BAM binary header");

 clean:
    if (h != NULL) {
        h->n_targets = num_names; // ensure we free only allocated target_names
        sam_hdr_destroy(h);
    }
    return NULL;
}

int bam_hdr_write(BGZF *fp, const sam_hdr_t *h)
{
    int32_t i, name_len, x;
    kstring_t hdr_ks = { 0, 0, NULL };
    char *text;
    uint32_t l_text;

    if (!h) return -1;

    if (h->hrecs) {
        if (sam_hrecs_rebuild_text(h->hrecs, &hdr_ks) != 0) return -1;
        if (hdr_ks.l > UINT32_MAX) {
            hts_log_error("Header too long for BAM format");
            free(hdr_ks.s);
            return -1;
        } else if (hdr_ks.l > INT32_MAX) {
            hts_log_warning("Header too long for BAM specification (>2GB)");
            hts_log_warning("Output file may not be portable");
        }
        text = hdr_ks.s;
        l_text = hdr_ks.l;
    } else {
        if (h->l_text > UINT32_MAX) {
            hts_log_error("Header too long for BAM format");
            return -1;
        } else if (h->l_text > INT32_MAX) {
            hts_log_warning("Header too long for BAM specification (>2GB)");
            hts_log_warning("Output file may not be portable");
        }
        text = h->text;
        l_text = h->l_text;
    }
    // write "BAM1"
    if (bgzf_write(fp, "BAM\1", 4) < 0) { free(hdr_ks.s); return -1; }
    // write plain text and the number of reference sequences
    if (fp->is_be) {
        x = ed_swap_4(l_text);
        if (bgzf_write(fp, &x, 4) < 0) { free(hdr_ks.s); return -1; }
        if (l_text) {
            if (bgzf_write(fp, text, l_text) < 0) { free(hdr_ks.s); return -1; }
        }
        x = ed_swap_4(h->n_targets);
        if (bgzf_write(fp, &x, 4) < 0) { free(hdr_ks.s); return -1; }
    } else {
        if (bgzf_write(fp, &l_text, 4) < 0) { free(hdr_ks.s); return -1; }
        if (l_text) {
            if (bgzf_write(fp, text, l_text) < 0) { free(hdr_ks.s); return -1; }
        }
        if (bgzf_write(fp, &h->n_targets, 4) < 0) { free(hdr_ks.s); return -1; }
    }
    free(hdr_ks.s);
    // write sequence names and lengths
    for (i = 0; i != h->n_targets; ++i) {
        char *p = h->target_name[i];
        name_len = strlen(p) + 1;
        if (fp->is_be) {
            x = ed_swap_4(name_len);
            if (bgzf_write(fp, &x, 4) < 0) return -1;
        } else {
            if (bgzf_write(fp, &name_len, 4) < 0) return -1;
        }
        if (bgzf_write(fp, p, name_len) < 0) return -1;
        if (fp->is_be) {
            x = ed_swap_4(h->target_len[i]);
            if (bgzf_write(fp, &x, 4) < 0) return -1;
        } else {
            if (bgzf_write(fp, &h->target_len[i], 4) < 0) return -1;
        }
    }
    if (bgzf_flush(fp) < 0) return -1;
    return 0;
}

// Wrap around bam_name2id() to get the right signature for hts_name2id_f
static int bam_name2id_wrapper(void *vhdr, const char *ref) {
    return bam_name2id((sam_hdr_t *) vhdr, ref);
}

const char *sam_parse_region(sam_hdr_t *h, const char *s, int *tid,
                             hts_pos_t *beg, hts_pos_t *end, int flags) {
    return hts_parse_region(s, tid, beg, end, bam_name2id_wrapper, h, flags);
}

/*************************
 *** BAM alignment I/O ***
 *************************/

bam1_t *bam_init1(void)
{
    return (bam1_t*)calloc(1, sizeof(bam1_t));
}

int sam_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data); // next power of 2
    new_m_data += 32; // reduces malloc arena migrations?
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (new_m_data > FUZZ_ALLOC_LIMIT) {
        errno = ENOMEM;
        return -1;
    }
#endif
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

void bam_destroy1(bam1_t *b)
{
    if (b == 0) return;
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        free(b->data);
        if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) != 0) {
            // In case of reuse
            b->data = NULL;
            b->m_data = 0;
            b->l_data = 0;
        }
    }

    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) == 0)
        free(b);
}

bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
    if (realloc_bam_data(bdst, bsrc->l_data) < 0) return NULL;
    memcpy(bdst->data, bsrc->data, bsrc->l_data); // copy var-len data
    memcpy(&bdst->core, &bsrc->core, sizeof(bsrc->core)); // copy the rest
    bdst->l_data = bsrc->l_data;
    bdst->id = bsrc->id;
    return bdst;
}

bam1_t *bam_dup1(const bam1_t *bsrc)
{
    if (bsrc == NULL) return NULL;
    bam1_t *bdst = bam_init1();
    if (bdst == NULL) return NULL;
    if (bam_copy1(bdst, bsrc) == NULL) {
        bam_destroy1(bdst);
        return NULL;
    }
    return bdst;
}

static void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                             hts_pos_t *rlen, hts_pos_t *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}

static int subtract_check_underflow(size_t length, size_t *limit)
{
    if (length <= *limit) {
        *limit -= length;
        return 0;
    }

    return -1;
}

int bam_set1(bam1_t *bam,
             size_t l_qname, const char *qname,
             uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
             size_t n_cigar, const uint32_t *cigar,
             int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
             size_t l_seq, const char *seq, const char *qual,
             size_t l_aux)
{
    // use a default qname "*" if none is provided
    if (l_qname == 0) {
        l_qname = 1;
        qname = "*";
    }

    // note: the qname is stored nul terminated and padded as described in the
    // documentation for the bam1_t struct.
    size_t qname_nuls = 4 - l_qname % 4;

    // the aligment length, needed for bam_reg2bin(), is calculated as in bam_endpos().
    // can't use bam_endpos() directly as some fields not yet set up.
    hts_pos_t rlen = 0, qlen = 0;
    if (!(flag & BAM_FUNMAP)) {
        bam_cigar2rqlens((int)n_cigar, cigar, &rlen, &qlen);
    }
    if (rlen == 0) {
        rlen = 1;
    }

    // validate parameters
    if (l_qname > 254) {
        hts_log_error("Query name too long");
        errno = EINVAL;
        return -1;
    }
    if (HTS_POS_MAX - rlen <= pos) {
        hts_log_error("Read ends beyond highest supported position");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && n_cigar == 0) {
        hts_log_error("Mapped query must have a CIGAR");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && l_seq != qlen) {
        hts_log_error("CIGAR and query sequence are of different length");
        errno = EINVAL;
        return -1;
    }

    size_t limit = INT32_MAX;
    int u = subtract_check_underflow(l_qname + qname_nuls, &limit);
    u    += subtract_check_underflow(n_cigar * 4, &limit);
    u    += subtract_check_underflow((l_seq + 1) / 2, &limit);
    u    += subtract_check_underflow(l_seq, &limit);
    u    += subtract_check_underflow(l_aux, &limit);
    if (u != 0) {
        hts_log_error("Size overflow");
        errno = EINVAL;
        return -1;
    }

    // re-allocate the data buffer as needed.
    size_t data_len = l_qname + qname_nuls + n_cigar * 4 + (l_seq + 1) / 2 + l_seq;
    if (realloc_bam_data(bam, data_len + l_aux) < 0) {
        return -1;
    }

    bam->l_data = (int)data_len;
    bam->core.pos = pos;
    bam->core.tid = tid;
    bam->core.bin = bam_reg2bin(pos, pos + rlen);
    bam->core.qual = mapq;
    bam->core.l_extranul = (uint8_t)(qname_nuls - 1);
    bam->core.flag = flag;
    bam->core.l_qname = (uint16_t)(l_qname + qname_nuls);
    bam->core.n_cigar = (uint32_t)n_cigar;
    bam->core.l_qseq = (int32_t)l_seq;
    bam->core.mtid = mtid;
    bam->core.mpos = mpos;
    bam->core.isize = isize;

    uint8_t *cp = bam->data;
    strncpy((char *)cp, qname, l_qname);
    int i;
    for (i = 0; i < qname_nuls; i++) {
        cp[l_qname + i] = '\0';
    }
    cp += l_qname + qname_nuls;

    if (n_cigar > 0) {
        memcpy(cp, cigar, n_cigar * 4);
    }
    cp += n_cigar * 4;

#define NN 16
    const uint8_t *useq = (uint8_t *)seq;
    for (i = 0; i + NN < l_seq; i += NN) {
        int j;
        const uint8_t *u2 = useq+i;
        for (j = 0; j < NN/2; j++)
            cp[j] = (seq_nt16_table[u2[j*2]]<<4) | seq_nt16_table[u2[j*2+1]];
        cp += NN/2;
    }
    for (; i + 1 < l_seq; i += 2) {
        *cp++ = (seq_nt16_table[useq[i]] << 4) | seq_nt16_table[useq[i + 1]];
    }

    for (; i < l_seq; i++) {
        *cp++ = seq_nt16_table[(unsigned char)seq[i]] << 4;
    }

    if (qual) {
        memcpy(cp, qual, l_seq);
    }
    else {
        memset(cp, '\xff', l_seq);
    }

    return (int)data_len;
}

hts_pos_t bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_endpos(const bam1_t *b)
{
    hts_pos_t rlen = (b->core.flag & BAM_FUNMAP)? 0 : bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    if (rlen == 0) rlen = 1;
    return b->core.pos + rlen;
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;

    // Bail out as fast as possible for the easy case
    uint32_t test_CG = BAM_CSOFT_CLIP | (c->l_qseq << BAM_CIGAR_SHIFT);
    if (c->n_cigar == 0 || test_CG != *bam_get_cigar(b))
        return 0;

    // The above isn't fool proof - we may have old CIGAR tags that aren't used,
    // but this is much less likely so do as a secondary check.
    if (c->tid < 0 || c->pos < 0)
        return 0;

    // Do we have a CG tag?
    uint8_t *CG = bam_aux_get(b, "CG");
    int saved_errno = errno;
    if (!CG) {
        if (errno != ENOENT) return -1;  // Bad aux data
        errno = saved_errno; // restore errno on expected no-CG-tag case
        return 0;
    }

    // Now we start with the serious work migrating CG to CIGAR
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data,
        *cigar0, CG_len, fake_bytes;
    cigar0 = bam_get_cigar(b);
    fake_bytes = c->n_cigar * 4;
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    // don't move if the real CIGAR length is shorter than the fake cigar length
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0;

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    b->l_data = b->l_data - fake_bytes + n_cigar4;
    // insert c->n_cigar-fake_bytes empty space to make room
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes));
    // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4);
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_warning("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

static inline int aux_type2size(uint8_t type)
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

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

// Fix bad records where qname is not terminated correctly.
static int fixup_missing_qname_nul(bam1_t *b) {
    bam1_core_t *c = &b->core;

    // Note this is called before c->l_extranul is added to c->l_qname
    if (c->l_extranul > 0) {
        b->data[c->l_qname++] = '\0';
        c->l_extranul--;
    } else {
        if (b->l_data > INT_MAX - 4) return -1;
        if (realloc_bam_data(b, b->l_data + 4) < 0) return -1;
        b->l_data += 4;
        b->data[c->l_qname++] = '\0';
        c->l_extranul = 3;
    }
    return 0;
}

static int bam_decode1_body(BGZF *fp, bam1_t *b, int32_t block_len,
                            const uint8_t *body)
{
    bam1_core_t *c = &b->core;
    const uint8_t *x = body;
    int raw_l_qname, rest_len, i;
    uint32_t new_l_data;

    b->l_data = 0;

    if (block_len < 32)
        return -4;

    c->tid        = le_to_u32(x);
    c->pos        = le_to_i32(x+4);
    uint32_t x2   = le_to_u32(x+8);
    c->bin        = x2>>16;
    c->qual       = x2>>8&0xff;
    c->l_qname    = x2&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    uint32_t x3   = le_to_u32(x+12);
    c->flag       = x3>>16;
    c->n_cigar    = x3&0xffff;
    c->l_qseq     = le_to_u32(x+16);
    c->mtid       = le_to_u32(x+20);
    c->mpos       = le_to_i32(x+24);
    c->isize      = le_to_i32(x+28);

    raw_l_qname = c->l_qname;
    new_l_data = block_len - 32 + c->l_extranul;
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1)
        return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    if (realloc_bam_data(b, new_l_data) < 0)
        return -4;
    b->l_data = new_l_data;

    memcpy(b->data, body + 32, raw_l_qname);
    if (b->data[raw_l_qname - 1] != '\0') {
        if (fixup_missing_qname_nul(b) < 0)
            return -4;
    }
    for (i = 0; i < c->l_extranul; ++i)
        b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;

    rest_len = block_len - 32 - raw_l_qname;
    if (rest_len < 0 || b->l_data < c->l_qname ||
        b->l_data - c->l_qname != rest_len)
        return -4;
    memcpy(b->data + c->l_qname, body + 32 + raw_l_qname,
           (size_t)rest_len);

    if (fp && fp->is_be)
        swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0)
        return -4;

    if (c->n_cigar > 0) {
        hts_pos_t rlen, qlen;
        bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP) || rlen == 0)
            rlen = 1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) &&
            qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                          bam_get_qname(b));
            return -4;
        }
    }

    return 4 + block_len;
}

static int bam_validate1_body_core(int32_t block_len, const uint8_t *body)
{
    const uint8_t *x = body;
    int l_qname, l_extranul, n_cigar;
    int32_t l_qseq;
    uint32_t new_l_data;
    uint32_t x2, x3;

    if (block_len < 32)
        return -4;

    x2 = le_to_u32(x+8);
    l_qname = x2 & 0xff;
    l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
    x3 = le_to_u32(x+12);
    n_cigar = x3 & 0xffff;
    l_qseq = le_to_u32(x+16);

    new_l_data = block_len - 32 + l_extranul;
    if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1)
        return -4;
    if (((uint64_t)n_cigar << 2) + l_qname + l_extranul
        + (((uint64_t)l_qseq + 1) >> 1) + l_qseq > (uint64_t)new_l_data)
        return -4;
    return 0;
}

/*
 * Note a second interface that returns a bam pointer instead would avoid bam_copy1
 * in multi-threaded handling.  This may be worth considering for htslib2.
 */
int bam_read1(BGZF *fp, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t new_l_data;
    uint8_t tmp[32], *x;

    b->l_data = 0;

    if ((ret = bgzf_read_small(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
    if (fp->is_be)
        ed_swap_4p(&block_len);
    if (block_len < 32) return -4;  // block_len includes core data
    if (fp->block_length - fp->block_offset > 32) {
        // Avoid bgzf_read and a temporary copy to a local buffer
        x = (uint8_t *)fp->uncompressed_block + fp->block_offset;
        fp->block_offset += 32;
    } else {
        x = tmp;
        if (bgzf_read(fp, x, 32) != 32) return -3;
    }

    c->tid        = le_to_u32(x);
    c->pos        = le_to_i32(x+4);
    uint32_t x2   = le_to_u32(x+8);
    c->bin        = x2>>16;
    c->qual       = x2>>8&0xff;
    c->l_qname    = x2&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    uint32_t x3   = le_to_u32(x+12);
    c->flag       = x3>>16;
    c->n_cigar    = x3&0xffff;
    c->l_qseq     = le_to_u32(x+16);
    c->mtid       = le_to_u32(x+20);
    c->mpos       = le_to_i32(x+24);
    c->isize      = le_to_i32(x+28);

    new_l_data = block_len - 32 + c->l_extranul;
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    if (realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (bgzf_read_small(fp, b->data, c->l_qname) != c->l_qname) return -4;
    if (b->data[c->l_qname - 1] != '\0') { // try to fix missing nul termination
        if (fixup_missing_qname_nul(b) < 0) return -4;
    }
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;
    if (b->l_data < c->l_qname ||
        bgzf_read_small(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
        return -4;
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0)
        return -4;

    if (c->n_cigar > 0) {
        // recompute "bin" and check CIGAR-qlen consistency
        hts_pos_t rlen, qlen;
        bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP) || rlen == 0) rlen = 1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                    bam_get_qname(b));
            return -4;
        }
    }

    return 4 + block_len;
}

int bam_write1(BGZF *fp, const bam1_t *b)
{
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (bgzf_flush_try(fp, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (bgzf_write_small(fp, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) ok = (bgzf_write_small(fp, &block_len, 4) >= 0);
    }
    if (ok) ok = (bgzf_write_small(fp, x, 32) >= 0);
    if (ok) ok = (bgzf_write_small(fp, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (bgzf_write_small(fp, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1<<28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"PRIhts_pos
                          " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                          bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) ok = (bgzf_write_small(fp, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        if (ok) ok = (bgzf_write_small(fp, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
        if (ok) ok = (bgzf_write_small(fp, "CGBI", 4) >= 0); // write CG:B,I
        u32_to_le(c->n_cigar, buf);
        if (ok) ok = (bgzf_write_small(fp, buf, 4) >= 0); // write the true CIGAR length
        if (ok) ok = (bgzf_write_small(fp, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
}

static int bgzf_raw_read_exact(BGZF *fp, void *data, size_t length)
{
    uint8_t *dst = (uint8_t *) data;

    while (length > 0) {
        ssize_t n = bgzf_raw_read(fp, dst, length);
        if (n <= 0) {
            if (n == 0)
                errno = EIO;
            return -1;
        }
        dst += n;
        length -= n;
    }
    return 0;
}

static int bam_bgzf_header_is_valid(const uint8_t header[18])
{
    return header[0] == 31 && header[1] == 139 && header[2] == 8
        && (header[3] & 4) != 0
        && le_to_u16(header + 10) == 6
        && header[12] == 'B' && header[13] == 'C'
        && le_to_u16(header + 14) == 2;
}

static int bam_bgzf_is_eof_marker(const uint8_t *block, uint16_t block_len)
{
    static const uint8_t eof_marker[28] = {
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    };

    return block_len == sizeof(eof_marker) &&
        memcmp(block, eof_marker, sizeof(eof_marker)) == 0;
}

int sam_bam_raw_copy_blocks(htsFile *in, htsFile *out)
{
    const size_t buf_size = 1024 * 1024;
    BGZF *ib = NULL, *ob = NULL;
    uint8_t *buf = NULL;
    int ret = -1;

    if (!in || !out || in->is_write || !out->is_write ||
        in->format.format != bam || out->format.format != bam ||
        !in->is_bgzf || !out->is_bgzf || !in->fp.bgzf || !out->fp.bgzf) {
        errno = EINVAL;
        return -1;
    }

    ib = in->fp.bgzf;
    ob = out->fp.bgzf;
    if (ib->mt || ob->mt || ib->is_gzip || ob->is_gzip ||
        !ib->is_compressed || !ob->is_compressed ||
        ob->compress_level != Z_DEFAULT_COMPRESSION ||
        out->idx || ob->idx || ob->idx_build_otf) {
        errno = EINVAL;
        return -1;
    }

    if (ib->block_length == 0 && ib->block_offset > 0) {
        if (bgzf_read_block(ib) < 0)
            return -1;
        if (ib->block_offset > ib->block_length) {
            errno = EINVAL;
            return -1;
        }
    }

    if (ib->block_offset < ib->block_length) {
        int available = ib->block_length - ib->block_offset;
        uint8_t *src = (uint8_t *) ib->uncompressed_block + ib->block_offset;

        if (bgzf_write(ob, src, available) != available)
            return -1;
        ib->block_offset = ib->block_length;
    }

    if (bgzf_flush(ob) < 0)
        return -1;

    buf = malloc(BGZF_MAX_BLOCK_SIZE > buf_size ? BGZF_MAX_BLOCK_SIZE : buf_size);
    if (!buf)
        return -1;

    for (;;) {
        uint8_t header[18];
        uint16_t block_len;
        uint32_t isize;
        ssize_t n = bgzf_raw_read(ib, header, sizeof(header));

        if (n < 0) {
            break;
        } else if (n == 0) {
            hts_log_warning("BGZF EOF marker is absent in raw BAM block copy");
            ret = 0;
            break;
        } else if (n != sizeof(header)) {
            hts_log_error("Truncated BGZF header in raw BAM block copy");
            errno = EIO;
            break;
        }
        if (!bam_bgzf_header_is_valid(header)) {
            hts_log_error("Invalid BGZF header in raw BAM block copy");
            errno = EFTYPE;
            break;
        }

        block_len = le_to_u16(header + 16) + 1;
        if (block_len < 26) {
            hts_log_error("Invalid BGZF block length in raw BAM block copy");
            errno = EFTYPE;
            break;
        }

        memcpy(buf, header, sizeof(header));
        if (bgzf_raw_read_exact(ib, buf + sizeof(header),
                                block_len - sizeof(header)) < 0) {
            hts_log_error("Truncated BGZF block in raw BAM block copy");
            break;
        }

        isize = le_to_u32(buf + block_len - 4);
        if (isize == 0) {
            uint8_t trailing;

            if (!bam_bgzf_is_eof_marker(buf, block_len)) {
                hts_log_error("Invalid BGZF EOF marker in raw BAM block copy");
                errno = EFTYPE;
                break;
            }
            n = bgzf_raw_read(ib, &trailing, 1);
            if (n == 0) {
                ret = 0;
            } else if (n > 0) {
                hts_log_error("Trailing data after BGZF EOF marker in raw BAM block copy");
                errno = EFTYPE;
            }
            break;
        }

        if (bgzf_raw_write_full_block(ob, buf, block_len) < 0)
            break;
    }

    free(buf);
    return ret;
}

static uint64_t bam_voff_block_beg(const bgzf_block_data_t *block)
{
    return ((uint64_t)block->block_address) << 16;
}

static uint64_t bam_voff_block_end(const bgzf_block_data_t *block)
{
    if (block->uncomp_len >= 0x10000)
        return ((uint64_t)(block->block_address + block->comp_len)) << 16;
    return (((uint64_t)block->block_address) << 16) |
           (uint64_t)block->uncomp_len;
}

static int bam_voff_spans_valid(const sam_bam_voff_span_t *spans,
                                int n_spans)
{
    int i;

    if (n_spans < 0 || (n_spans > 0 && !spans)) {
        errno = EINVAL;
        return 0;
    }
    for (i = 0; i < n_spans; i++) {
        if (spans[i].beg >= spans[i].end ||
            (i > 0 && spans[i - 1].end > spans[i].beg)) {
            errno = EINVAL;
            return 0;
        }
    }
    return 1;
}

static int bam_voff_block_span_bounds(const bgzf_block_data_t *block,
                                      const sam_bam_voff_span_t *span,
                                      int *beg, int *end)
{
    uint64_t block_beg = bam_voff_block_beg(block);
    uint64_t block_end = bam_voff_block_end(block);
    uint64_t s_beg = span->beg > block_beg ? span->beg : block_beg;
    uint64_t s_end = span->end < block_end ? span->end : block_end;

    if (s_beg >= s_end)
        return 0;
    if ((s_beg >> 16) != (block_beg >> 16) ||
        (s_end != block_end && (s_end >> 16) != (block_beg >> 16))) {
        errno = EINVAL;
        return -1;
    }
    *beg = (int)(s_beg & 0xffff);
    *end = s_end == block_end ? block->uncomp_len : (int)(s_end & 0xffff);
    if (*beg < 0 || *end < *beg || *end > block->uncomp_len) {
        errno = EINVAL;
        return -1;
    }
    return 1;
}

static int bam_voff_block_fully_covered(const bgzf_block_data_t *block,
                                        const sam_bam_voff_span_t *spans,
                                        int n_spans, int span_i)
{
    uint64_t block_beg = bam_voff_block_beg(block);
    uint64_t block_end = bam_voff_block_end(block);
    int covered = 0;

    while (span_i < n_spans && spans[span_i].beg < block_end) {
        int beg, end, ret;

        if (spans[span_i].end <= block_beg) {
            span_i++;
            continue;
        }
        ret = bam_voff_block_span_bounds(block, &spans[span_i], &beg, &end);
        if (ret < 0)
            return 0;
        if (ret > 0) {
            if (beg != covered)
                return 0;
            covered = end;
            if (covered == block->uncomp_len)
                return 1;
        }
        if (spans[span_i].end <= block_end)
            span_i++;
        else
            break;
    }
    return covered == block->uncomp_len;
}

static int bam_voff_write_block_spans(BGZF *ob, bgzf_block_data_t *block,
                                      const sam_bam_voff_span_t *spans,
                                      int n_spans, int span_i)
{
    uint64_t block_beg = bam_voff_block_beg(block);
    uint64_t block_end = bam_voff_block_end(block);

    while (span_i < n_spans && spans[span_i].beg < block_end) {
        int beg, end, ret;

        if (spans[span_i].end <= block_beg) {
            span_i++;
            continue;
        }
        ret = bam_voff_block_span_bounds(block, &spans[span_i], &beg, &end);
        if (ret < 0)
            return -1;
        if (ret > 0 &&
            bgzf_write(ob, block->uncomp_data + beg,
                       (size_t)(end - beg)) != end - beg)
            return -1;
        if (spans[span_i].end <= block_end)
            span_i++;
        else
            break;
    }
    return 0;
}

static int bam_voff_block_selected_uncomp(const bgzf_block_data_t *block,
                                          const sam_bam_voff_span_t *spans,
                                          int n_spans, int span_i,
                                          uint64_t *selected)
{
    uint64_t block_beg = bam_voff_block_beg(block);
    uint64_t block_end = bam_voff_block_end(block);

    *selected = 0;
    while (span_i < n_spans && spans[span_i].beg < block_end) {
        int beg, end, ret;

        if (spans[span_i].end <= block_beg) {
            span_i++;
            continue;
        }
        ret = bam_voff_block_span_bounds(block, &spans[span_i], &beg, &end);
        if (ret < 0)
            return -1;
        if (ret > 0)
            *selected += (uint64_t)(end - beg);
        if (spans[span_i].end <= block_end)
            span_i++;
        else
            break;
    }
    return 0;
}

static int sam_bam_raw_copy_voff_spans_core(
        htsFile *in, htsFile *out, const sam_bam_voff_span_t *spans,
        int n_spans, sam_bam_voff_span_stats_t *stats, int do_write)
{
    BGZF *ib, *ob;
    bgzf_block_data_t *block = NULL;
    int span_i = 0, ret = -1;

    if (stats)
        memset(stats, 0, sizeof(*stats));
    if (!bam_voff_spans_valid(spans, n_spans))
        return -1;
    if (n_spans == 0)
        return 0;
    if (!in || in->is_write || in->format.format != bam ||
        !in->is_bgzf || !in->fp.bgzf) {
        errno = EINVAL;
        return -1;
    }

    ib = in->fp.bgzf;
    if (ib->mt || ib->is_gzip || !ib->is_compressed) {
        errno = EINVAL;
        return -1;
    }
    if (do_write) {
        if (!out || !out->is_write || out->format.format != bam ||
            !out->is_bgzf || !out->fp.bgzf) {
            errno = EINVAL;
            return -1;
        }
        ob = out->fp.bgzf;
        if (ob->is_gzip || !ob->is_compressed ||
            ob->compress_level != -1 ||
            out->idx || ob->idx || ob->idx_build_otf) {
            errno = EINVAL;
            return -1;
        }
    } else {
        ob = NULL;
    }

    block = malloc(sizeof(*block));
    if (!block)
        return -1;

    if (bgzf_seek(ib, (int64_t)(spans[0].beg & ~UINT64_C(0xffff)),
                  SEEK_SET) < 0)
        goto cleanup;

    while (span_i < n_spans) {
        uint64_t block_beg, block_end;

        if (bgzf_read_block_data(ib, block) < 0 || block->hit_eof) {
            errno = EIO;
            goto cleanup;
        }
        block_beg = bam_voff_block_beg(block);
        block_end = bam_voff_block_end(block);

        while (span_i < n_spans && spans[span_i].end <= block_beg)
            span_i++;
        if (span_i >= n_spans)
            break;
        if (spans[span_i].beg >= block_end)
            continue;

        if (bam_voff_block_fully_covered(block, spans, n_spans, span_i)) {
            if (stats) {
                stats->selected_uncomp += (uint64_t)block->uncomp_len;
                stats->full_uncomp += (uint64_t)block->uncomp_len;
                stats->full_blocks++;
            }
            if (do_write &&
                bgzf_raw_write_full_block(ob, block->comp_data,
                                          (size_t)block->comp_len) < 0)
                goto cleanup;
        } else {
            uint64_t selected = 0;

            if (bam_voff_block_selected_uncomp(block, spans, n_spans,
                                               span_i, &selected) < 0)
                goto cleanup;
            if (stats) {
                stats->selected_uncomp += selected;
                stats->partial_uncomp += selected;
                if (selected > 0)
                    stats->partial_blocks++;
            }
            if (do_write &&
                bam_voff_write_block_spans(ob, block, spans, n_spans,
                                           span_i) < 0)
                goto cleanup;
        }

        while (span_i < n_spans && spans[span_i].end <= block_end)
            span_i++;
    }

    ret = 0;

cleanup:
    free(block);
    return ret;
}

int sam_bam_raw_copy_voff_spans_stats(htsFile *in,
                                      const sam_bam_voff_span_t *spans,
                                      int n_spans,
                                      sam_bam_voff_span_stats_t *stats)
{
    if (!stats) {
        errno = EINVAL;
        return -1;
    }
    return sam_bam_raw_copy_voff_spans_core(in, NULL, spans, n_spans,
                                            stats, 0);
}

int sam_bam_raw_copy_voff_spans(htsFile *in, htsFile *out,
                                const sam_bam_voff_span_t *spans,
                                int n_spans)
{
    return sam_bam_raw_copy_voff_spans_core(in, out, spans, n_spans,
                                            NULL, 1);
}

/*
 * Write a BAM file and append to the in-memory index simultaneously.
 */
static int bam_write_idx1(htsFile *fp, const sam_hdr_t *h, const bam1_t *b) {
    BGZF *bfp = fp->fp.bgzf;

    if (!fp->idx)
        return bam_write1(bfp, b);

    uint32_t block_len = b->l_data - b->core.l_extranul + 32;
    if (b->core.n_cigar > 0xffff) block_len += 16;
    if (bgzf_flush_try(bfp, 4 + block_len) < 0)
        return -1;
    if (!bfp->mt)
        hts_idx_amend_last(fp->idx, bgzf_tell(bfp));

    int ret = bam_write1(bfp, b);
    if (ret < 0)
        return -1;

    if (bgzf_idx_push(bfp, fp->idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(bfp), !(b->core.flag&BAM_FUNMAP)) < 0) {
        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
        ret = -1;
    }

    return ret;
}

/*
 * Set the qname in a BAM record
 */
int bam_set_qname(bam1_t *rec, const char *qname)
{
    if (!rec) return -1;
    if (!qname || !*qname) return -1;

    size_t old_len = rec->core.l_qname;
    size_t new_len = strlen(qname) + 1;
    if (new_len < 1 || new_len > 255) return -1;

    int extranul = (new_len%4 != 0) ? (4 - new_len%4) : 0;

    size_t new_data_len = rec->l_data - old_len + new_len + extranul;
    if (realloc_bam_data(rec, new_data_len) < 0) return -1;

    // Make room
    if (new_len + extranul != rec->core.l_qname)
        memmove(rec->data + new_len + extranul, rec->data + rec->core.l_qname, rec->l_data - rec->core.l_qname);
    // Copy in new name and pad if needed
    memcpy(rec->data, qname, new_len);
    int n;
    for (n = 0; n < extranul; n++) rec->data[new_len + n] = '\0';

    rec->l_data = new_data_len;
    rec->core.l_qname = new_len + extranul;
    rec->core.l_extranul = extranul;

    return 0;
}

/********************
 *** BAM indexing ***
 ********************/

static hts_idx_t *sam_index(htsFile *fp, int min_shift)
{
    int n_lvls, i, fmt, ret;
    bam1_t *b;
    hts_idx_t *idx;
    sam_hdr_t *h;
    h = sam_hdr_read(fp);
    if (h == NULL) return NULL;
    if (min_shift > 0) {
        hts_pos_t max_len = 0;
        for (i = 0; i < h->n_targets; ++i) {
            hts_pos_t len = sam_hdr_tid2len(h, i);
            if (max_len < len) max_len = len;
        }
        n_lvls = 0;
        hts_adjust_csi_settings(max_len, &min_shift, &n_lvls);
        fmt = HTS_FMT_CSI;
    } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;
    idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
    b = bam_init1();
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        ret = hts_idx_push(idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP));
        if (ret < 0) { // unsorted or doesn't fit
            hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed", bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
            goto err;
        }
    }
    if (ret < -1) goto err; // corrupted BAM file

    hts_idx_finish(idx, bgzf_tell(fp->fp.bgzf));
    sam_hdr_destroy(h);
    bam_destroy1(b);
    return idx;

err:
    bam_destroy1(b);
    hts_idx_destroy(idx);
    return NULL;
}

int sam_index_build3(const char *fn, const char *fnidx, int min_shift, int nthreads)
{
    hts_idx_t *idx;
    htsFile *fp;
    int ret = 0;

    if ((fp = hts_open(fn, "r")) == 0) return -2;
    if (nthreads)
        hts_set_threads(fp, nthreads);

    switch (fp->format.format) {
    case cram:

        ret = cram_index_build(fp->fp.cram, fn, fnidx);
        break;

    case bam:
    case sam:
        if (fp->format.compression != bgzf) {
            hts_log_error("%s file \"%s\" not BGZF compressed",
                          fp->format.format == bam ? "BAM" : "SAM", fn);
            ret = -1;
            break;
        }
        idx = sam_index(fp, min_shift);
        if (idx) {
            ret = hts_idx_save_as(idx, fn, fnidx, (min_shift > 0)? HTS_FMT_CSI : HTS_FMT_BAI);
            if (ret < 0) ret = -4;
            hts_idx_destroy(idx);
        }
        else ret = -1;
        break;

    default:
        ret = -3;
        break;
    }
    hts_close(fp);

    return ret;
}

int sam_index_build2(const char *fn, const char *fnidx, int min_shift)
{
    return sam_index_build3(fn, fnidx, min_shift, 0);
}

int sam_index_build(const char *fn, int min_shift)
{
    return sam_index_build3(fn, NULL, min_shift, 0);
}

// Provide bam_index_build() symbol for binary compatibility with earlier HTSlib
#undef bam_index_build
int bam_index_build(const char *fn, int min_shift)
{
    return sam_index_build2(fn, NULL, min_shift);
}

// Initialise fp->idx for the current format type.
// This must be called after the header has been written but no other data.
int sam_idx_init(htsFile *fp, sam_hdr_t *h, int min_shift, const char *fnidx) {
    fp->fnidx = fnidx;
    if (fp->format.format == bam || fp->format.format == bcf ||
        (fp->format.format == sam && fp->format.compression == bgzf)) {
        int n_lvls, fmt = HTS_FMT_CSI;
        if (min_shift > 0) {
            int64_t max_len = 0;
            int i;
            for (i = 0; i < h->n_targets; ++i)
                if (max_len < h->target_len[i]) max_len = h->target_len[i];
            n_lvls = 0;
            hts_adjust_csi_settings(max_len, &min_shift, &n_lvls);
        } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;

        fp->idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
        return fp->idx ? 0 : -1;
    }

    if (fp->format.format == cram) {
        fp->fp.cram->idxfp = bgzf_open(fnidx, "wg");
        return fp->fp.cram->idxfp ? 0 : -1;
    }

    return -1;
}

// Finishes an index. Call after the last record has been written.
// Returns 0 on success, <0 on failure.
int sam_idx_save(htsFile *fp) {
    if (fp->format.format == bam || fp->format.format == bcf ||
        fp->format.format == vcf || fp->format.format == sam) {
        int ret;
        if ((ret = sam_state_destroy(fp)) < 0) {
            errno = -ret;
            return -1;
        }
        if (!fp->is_bgzf || bgzf_flush(fp->fp.bgzf) < 0)
            return -1;
        hts_idx_amend_last(fp->idx, bgzf_tell(fp->fp.bgzf));

        if (hts_idx_finish(fp->idx, bgzf_tell(fp->fp.bgzf)) < 0)
            return -1;

        return hts_idx_save_but_not_close(fp->idx, fp->fnidx, hts_idx_fmt(fp->idx));

    } else if (fp->format.format == cram) {
        // flushed and closed by cram_close
    }

    return 0;
}

static int sam_readrec(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = (htsFile *)fpv;
    bam1_t *b = bv;
    fp->line.l = 0;
    int ret = sam_read1(fp, fp->bam_header, b);
    if (ret >= 0) {
        *tid = b->core.tid;
        *beg = b->core.pos;
        *end = bam_endpos(b);
    }
    return ret;
}

// This is used only with read_rest=1 iterators, so need not set tid/beg/end.
static int sam_readrec_rest(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = (htsFile *)fpv;
    bam1_t *b = bv;
    fp->line.l = 0;
    int ret = sam_read1(fp, fp->bam_header, b);
    return ret;
}

// Internal (for now) func used by bam_sym_lookup.  This is copied from
// samtools/bam.c.
static const char *bam_get_library(const bam_hdr_t *h, const bam1_t *b)
{
    const char *rg;
    kstring_t lib = { 0, 0, NULL };
    rg = (char *)bam_aux_get(b, "RG");

    if (!rg)
        return NULL;
    else
        rg++;

    if (sam_hdr_find_tag_id((bam_hdr_t *)h, "RG", "ID", rg, "LB", &lib)  < 0)
        return NULL;

    static char LB_text[1024];
    int len = lib.l < sizeof(LB_text) - 1 ? lib.l : sizeof(LB_text) - 1;

    memcpy(LB_text, lib.s, len);
    LB_text[len] = 0;

    free(lib.s);

    return LB_text;
}


// Bam record pointer and SAM header combined
typedef struct {
    const sam_hdr_t *h;
    const bam1_t *b;
} hb_pair;

// Looks up variable names in str and replaces them with their value.
// Also supports aux tags.
//
// Note the expression parser deliberately overallocates str size so it
// is safe to use memcmp over strcmp.
static int bam_sym_lookup(void *data, char *str, char **end,
                          hts_expr_val_t *res) {
    hb_pair *hb = (hb_pair *)data;
    const bam1_t *b = hb->b;

    res->is_str = 0;
    switch(*str) {
    case 'c':
        if (memcmp(str, "cigar", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            ks_clear(&res->s);
            uint32_t *cigar = bam_get_cigar(b);
            int i, n = b->core.n_cigar, r = 0;
            if (n) {
                for (i = 0; i < n; i++) {
                    r |= kputw (bam_cigar_oplen(cigar[i]), &res->s) < 0;
                    r |= kputc_(bam_cigar_opchr(cigar[i]), &res->s) < 0;
                }
                r |= kputs("", &res->s) < 0;
            } else {
                r |= kputs("*", &res->s) < 0;
            }
            return r ? -1 : 0;
        }
        break;

    case 'e':
        if (memcmp(str, "endpos", 6) == 0) {
            *end = str+6;
            res->d = bam_endpos(b);
            return 0;
        }
        break;

    case 'f':
        if (memcmp(str, "flag", 4) == 0) {
            str = *end = str+4;
            if (*str != '.') {
                res->d = b->core.flag;
                return 0;
            } else {
                str++;
                if (!memcmp(str, "paired", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FPAIRED;
                    return 0;
                } else if (!memcmp(str, "proper_pair", 11)) {
                    *end = str+11;
                    res->d = b->core.flag & BAM_FPROPER_PAIR;
                    return 0;
                } else if (!memcmp(str, "unmap", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FUNMAP;
                    return 0;
                } else if (!memcmp(str, "munmap", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FMUNMAP;
                    return 0;
                } else if (!memcmp(str, "reverse", 7)) {
                    *end = str+7;
                    res->d = b->core.flag & BAM_FREVERSE;
                    return 0;
                } else if (!memcmp(str, "mreverse", 8)) {
                    *end = str+8;
                    res->d = b->core.flag & BAM_FMREVERSE;
                    return 0;
                } else if (!memcmp(str, "read1", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FREAD1;
                    return 0;
                } else if (!memcmp(str, "read2", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FREAD2;
                    return 0;
                } else if (!memcmp(str, "secondary", 9)) {
                    *end = str+9;
                    res->d = b->core.flag & BAM_FSECONDARY;
                    return 0;
                } else if (!memcmp(str, "qcfail", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FQCFAIL;
                    return 0;
                } else if (!memcmp(str, "dup", 3)) {
                    *end = str+3;
                    res->d = b->core.flag & BAM_FDUP;
                    return 0;
                } else if (!memcmp(str, "supplementary", 13)) {
                    *end = str+13;
                    res->d = b->core.flag & BAM_FSUPPLEMENTARY;
                    return 0;
                } else {
                    hts_log_error("Unrecognised flag string");
                    return -1;
                }
            }
        }
        break;

    case 'h':
        if (memcmp(str, "hclen", 5) == 0) {
            int hclen = 0;
            uint32_t *cigar = bam_get_cigar(b);
            uint32_t ncigar = b->core.n_cigar;

            // left
            if (ncigar > 0 && bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)
                hclen = bam_cigar_oplen(cigar[0]);

            // right
            if (ncigar > 1 && bam_cigar_op(cigar[ncigar-1]) == BAM_CHARD_CLIP)
                hclen += bam_cigar_oplen(cigar[ncigar-1]);

            *end = str+5;
            res->d = hclen;
            return 0;
        }
        break;

    case 'l':
        if (memcmp(str, "library", 7) == 0) {
            *end = str+7;
            res->is_str = 1;
            const char *lib = bam_get_library(hb->h, b);
            kputs(lib ? lib : "", ks_clear(&res->s));
            return 0;
        }
        break;

    case 'm':
        if (memcmp(str, "mapq", 4) == 0) {
            *end = str+4;
            res->d = b->core.qual;
            return 0;
        } else if (memcmp(str, "mpos", 4) == 0) {
            *end = str+4;
            res->d = b->core.mpos+1;
            return 0;
        } else if (memcmp(str, "mrname", 6) == 0) {
            *end = str+6;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.mtid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "mrefid", 6) == 0) {
            *end = str+6;
            res->d = b->core.mtid;
            return 0;
        }
        break;

    case 'n':
        if (memcmp(str, "ncigar", 6) == 0) {
            *end = str+6;
            res->d = b->core.n_cigar;
            return 0;
        }
        break;

    case 'p':
        if (memcmp(str, "pos", 3) == 0) {
            *end = str+3;
            res->d = b->core.pos+1;
            return 0;
        } else if (memcmp(str, "pnext", 5) == 0) {
            *end = str+5;
            res->d = b->core.mpos+1;
            return 0;
        }
        break;

    case 'q':
        if (memcmp(str, "qlen", 4) == 0) {
            *end = str+4;
            res->d = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
            return 0;
        } else if (memcmp(str, "qname", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            kputs(bam_get_qname(b), ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "qual", 4) == 0) {
            *end = str+4;
            ks_clear(&res->s);
            if (ks_resize(&res->s, b->core.l_qseq+1) < 0)
                return -1;
            memcpy(res->s.s, bam_get_qual(b), b->core.l_qseq);
            res->s.l = b->core.l_qseq;
            res->is_str = 1;
            return 0;
        }
        break;

    case 'r':
        if (memcmp(str, "rlen", 4) == 0) {
            *end = str+4;
            res->d = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
            return 0;
        } else if (memcmp(str, "rname", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.tid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "rnext", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.mtid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "refid", 5) == 0) {
            *end = str+5;
            res->d = b->core.tid;
            return 0;
        }
        break;

    case 's':
        if (memcmp(str, "seq", 3) == 0) {
            *end = str+3;
            ks_clear(&res->s);
            if (ks_resize(&res->s, b->core.l_qseq+1) < 0)
                return -1;
            nibble2base(bam_get_seq(b), res->s.s, b->core.l_qseq);
            res->s.s[b->core.l_qseq] = 0;
            res->s.l = b->core.l_qseq;
            res->is_str = 1;
            return 0;
        } else if (memcmp(str, "sclen", 5) == 0) {
            int sclen = 0;
            uint32_t *cigar = bam_get_cigar(b);
            int ncigar = b->core.n_cigar;
            int left = 0;

            // left
            if (ncigar > 0
                && bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
                left = 0, sclen += bam_cigar_oplen(cigar[0]);
            else if (ncigar > 1
                     && bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP
                     && bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP)
                left = 1, sclen += bam_cigar_oplen(cigar[1]);

            // right
            if (ncigar-1 > left
                && bam_cigar_op(cigar[ncigar-1]) == BAM_CSOFT_CLIP)
                sclen += bam_cigar_oplen(cigar[ncigar-1]);
            else if (ncigar-2 > left
                     && bam_cigar_op(cigar[ncigar-1]) == BAM_CHARD_CLIP
                     && bam_cigar_op(cigar[ncigar-2]) == BAM_CSOFT_CLIP)
                sclen += bam_cigar_oplen(cigar[ncigar-2]);

            *end = str+5;
            res->d = sclen;
            return 0;
        }
        break;

    case 't':
        if (memcmp(str, "tlen", 4) == 0) {
            *end = str+4;
            res->d = b->core.isize;
            return 0;
        }
        break;

    case '[':
        if (*str == '[' && str[1] && str[2] && str[3] == ']') {
            /* aux tags */
            *end = str+4;

            uint8_t *aux = bam_aux_get(b, str+1);
            if (aux) {
                // we define the truth of a tag to be its presence, even if 0.
                res->is_true = 1;
                switch (*aux) {
                case 'Z':
                case 'H':
                    res->is_str = 1;
                    kputs((char *)aux+1, ks_clear(&res->s));
                    break;

                case 'A':
                    res->is_str = 1;
                    kputsn((char *)aux+1, 1, ks_clear(&res->s));
                    break;

                case 'i': case 'I':
                case 's': case 'S':
                case 'c': case 'C':
                    res->is_str = 0;
                    res->d = bam_aux2i(aux);
                    break;

                case 'f':
                case 'd':
                    res->is_str = 0;
                    res->d = bam_aux2f(aux);
                    break;

                default:
                    hts_log_error("Aux type '%c not yet supported by filters",
                                  *aux);
                    return -1;
                }
                return 0;

            } else {
                // hence absent tags are always false (and strings)
                res->is_str = 1;
                res->s.l = 0;
                res->d = 0;
                res->is_true = 0;
                return 0;
            }
        }
        break;
    }

    // All successful matches in switch should return 0.
    // So if we didn't match, it's a parse error.
    return -1;
}

// Returns 1 when accepted by the filter, 0 if not, -1 on error.
int sam_passes_filter(const sam_hdr_t *h, const bam1_t *b, hts_filter_t *filt)
{
    hb_pair hb = {h, b};
    hts_expr_val_t res = HTS_EXPR_VAL_INIT;
    if (hts_filter_eval2(filt, &hb, bam_sym_lookup, &res)) {
        hts_log_error("Couldn't process filter expression");
        hts_expr_val_free(&res);
        return -1;
    }

    int t = res.is_true;
    hts_expr_val_free(&res);

    return t;
}

typedef struct {
    const sam_hdr_t *h;
    bam_batch_record_t *record;
    bam1_t *scratch;
} hb_view_pair;

static inline uint32_t bam_view_cigar_op_at(const bam_batch_record_t *record,
                                            int idx)
{
    return le_to_u32(sam_bam_batch_record_cigar(record) + ((size_t)idx << 2));
}

static int bam_view_cigar_to_str(const bam_batch_record_t *record,
                                 kstring_t *str)
{
    int i, ret = 0;

    ks_clear(str);
    if (record->core.n_cigar == 0)
        return kputs("*", str) < 0 ? -1 : 0;
    for (i = 0; i < record->core.n_cigar; i++) {
        uint32_t c = bam_view_cigar_op_at(record, i);
        ret |= kputw(bam_cigar_oplen(c), str) < 0;
        ret |= kputc_(bam_cigar_opchr(c), str) < 0;
    }
    ret |= kputs("", str) < 0;
    return ret ? -1 : 0;
}

static hts_pos_t bam_view_cigar_rlen(const bam_batch_record_t *record)
{
    hts_pos_t rlen = 0;
    int i;

    for (i = 0; i < record->core.n_cigar; i++) {
        uint32_t c = bam_view_cigar_op_at(record, i);
        if (bam_cigar_type(bam_cigar_op(c)) & 2)
            rlen += bam_cigar_oplen(c);
    }
    return rlen;
}

static int bam_view_cigar_hclen(const bam_batch_record_t *record)
{
    int hclen = 0, n = record->core.n_cigar;

    if (n > 0) {
        uint32_t c = bam_view_cigar_op_at(record, 0);
        if (bam_cigar_op(c) == BAM_CHARD_CLIP)
            hclen += bam_cigar_oplen(c);
    }
    if (n > 1) {
        uint32_t c = bam_view_cigar_op_at(record, n - 1);
        if (bam_cigar_op(c) == BAM_CHARD_CLIP)
            hclen += bam_cigar_oplen(c);
    }
    return hclen;
}

static int bam_view_cigar_sclen(const bam_batch_record_t *record)
{
    int left = 0, sclen = 0, n = record->core.n_cigar;

    if (n > 0) {
        uint32_t c = bam_view_cigar_op_at(record, 0);
        if (bam_cigar_op(c) == BAM_CSOFT_CLIP) {
            sclen += bam_cigar_oplen(c);
        } else if (n > 1 &&
                   bam_cigar_op(c) == BAM_CHARD_CLIP) {
            uint32_t c1 = bam_view_cigar_op_at(record, 1);
            if (bam_cigar_op(c1) == BAM_CSOFT_CLIP) {
                left = 1;
                sclen += bam_cigar_oplen(c1);
            }
        }
    }

    if (n - 1 > left) {
        uint32_t c = bam_view_cigar_op_at(record, n - 1);
        if (bam_cigar_op(c) == BAM_CSOFT_CLIP) {
            sclen += bam_cigar_oplen(c);
        } else if (n - 2 > left &&
                   bam_cigar_op(c) == BAM_CHARD_CLIP) {
            uint32_t c1 = bam_view_cigar_op_at(record, n - 2);
            if (bam_cigar_op(c1) == BAM_CSOFT_CLIP)
                sclen += bam_cigar_oplen(c1);
        }
    }
    return sclen;
}

static int bam_view_get_library(const sam_hdr_t *h,
                                const bam_batch_record_t *record,
                                kstring_t *lib)
{
    sam_bam_record_view_t view;
    const uint8_t *rg;
    int found;

    sam_bam_record_view_from_batch(&view, record);
    rg = sam_bam_record_view_aux_get(&view, "RG");
    if (!rg)
        return 0;
    if (*rg != 'Z' && *rg != 'H')
        return 0;
    found = sam_hdr_find_tag_id((sam_hdr_t *)h, "RG", "ID",
                                (const char *)rg + 1, "LB", lib) < 0 ? 0 : 1;
    if (found && lib->l > 1023) {
        lib->l = 1023;
        lib->s[lib->l] = '\0';
    }
    return found;
}

static int bam_view_sym_lookup(void *data, char *str, char **end,
                               hts_expr_val_t *res) {
    hb_view_pair *hb = (hb_view_pair *)data;
    const bam_batch_record_t *record = hb->record;
    const bam1_core_t *core = &record->core;

    res->is_str = 0;
    switch(*str) {
    case 'c':
        if (memcmp(str, "cigar", 5) == 0) {
            *end = str + 5;
            res->is_str = 1;
            return bam_view_cigar_to_str(record, ks_clear(&res->s));
        }
        break;

    case 'e':
        if (memcmp(str, "endpos", 6) == 0) {
            hts_pos_t endpos;

            *end = str + 6;
            if (sam_bam_batch_record_endpos(record, hb->scratch,
                                            &endpos) < 0)
                return -1;
            res->d = endpos;
            return 0;
        }
        break;

    case 'f':
        if (memcmp(str, "flag", 4) == 0) {
            str = *end = str + 4;
            if (*str != '.') {
                res->d = core->flag;
                return 0;
            } else {
                str++;
                if (!memcmp(str, "paired", 6)) {
                    *end = str + 6;
                    res->d = core->flag & BAM_FPAIRED;
                    return 0;
                } else if (!memcmp(str, "proper_pair", 11)) {
                    *end = str + 11;
                    res->d = core->flag & BAM_FPROPER_PAIR;
                    return 0;
                } else if (!memcmp(str, "unmap", 5)) {
                    *end = str + 5;
                    res->d = core->flag & BAM_FUNMAP;
                    return 0;
                } else if (!memcmp(str, "munmap", 6)) {
                    *end = str + 6;
                    res->d = core->flag & BAM_FMUNMAP;
                    return 0;
                } else if (!memcmp(str, "reverse", 7)) {
                    *end = str + 7;
                    res->d = core->flag & BAM_FREVERSE;
                    return 0;
                } else if (!memcmp(str, "mreverse", 8)) {
                    *end = str + 8;
                    res->d = core->flag & BAM_FMREVERSE;
                    return 0;
                } else if (!memcmp(str, "read1", 5)) {
                    *end = str + 5;
                    res->d = core->flag & BAM_FREAD1;
                    return 0;
                } else if (!memcmp(str, "read2", 5)) {
                    *end = str + 5;
                    res->d = core->flag & BAM_FREAD2;
                    return 0;
                } else if (!memcmp(str, "secondary", 9)) {
                    *end = str + 9;
                    res->d = core->flag & BAM_FSECONDARY;
                    return 0;
                } else if (!memcmp(str, "qcfail", 6)) {
                    *end = str + 6;
                    res->d = core->flag & BAM_FQCFAIL;
                    return 0;
                } else if (!memcmp(str, "dup", 3)) {
                    *end = str + 3;
                    res->d = core->flag & BAM_FDUP;
                    return 0;
                } else if (!memcmp(str, "supplementary", 13)) {
                    *end = str + 13;
                    res->d = core->flag & BAM_FSUPPLEMENTARY;
                    return 0;
                } else {
                    hts_log_error("Unrecognised flag string");
                    return -1;
                }
            }
        }
        break;

    case 'h':
        if (memcmp(str, "hclen", 5) == 0) {
            *end = str + 5;
            res->d = bam_view_cigar_hclen(record);
            return 0;
        }
        break;

    case 'l':
        if (memcmp(str, "library", 7) == 0) {
            kstring_t lib = { 0, 0, NULL };
            int found;

            *end = str + 7;
            res->is_str = 1;
            found = bam_view_get_library(hb->h, record, &lib);
            if (found < 0) {
                free(lib.s);
                return -1;
            }
            if (kputs(found ? lib.s : "", ks_clear(&res->s)) < 0) {
                free(lib.s);
                return -1;
            }
            free(lib.s);
            return 0;
        }
        break;

    case 'm':
        if (memcmp(str, "mapq", 4) == 0) {
            *end = str + 4;
            res->d = core->qual;
            return 0;
        } else if (memcmp(str, "mpos", 4) == 0) {
            *end = str + 4;
            res->d = core->mpos + 1;
            return 0;
        } else if (memcmp(str, "mrname", 6) == 0) {
            const char *rn = sam_hdr_tid2name(hb->h, core->mtid);

            *end = str + 6;
            res->is_str = 1;
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "mrefid", 6) == 0) {
            *end = str + 6;
            res->d = core->mtid;
            return 0;
        }
        break;

    case 'n':
        if (memcmp(str, "ncigar", 6) == 0) {
            *end = str + 6;
            res->d = core->n_cigar;
            return 0;
        }
        break;

    case 'p':
        if (memcmp(str, "pos", 3) == 0) {
            *end = str + 3;
            res->d = core->pos + 1;
            return 0;
        } else if (memcmp(str, "pnext", 5) == 0) {
            *end = str + 5;
            res->d = core->mpos + 1;
            return 0;
        }
        break;

    case 'q':
        if (memcmp(str, "qlen", 4) == 0) {
            hts_pos_t qlen;

            *end = str + 4;
            if (sam_bam_batch_record_query_len(record, hb->scratch, 0,
                                               &qlen) < 0)
                return -1;
            res->d = qlen;
            return 0;
        } else if (memcmp(str, "qname", 5) == 0) {
            *end = str + 5;
            res->is_str = 1;
            kputs((const char *)sam_bam_batch_record_qname(record),
                  ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "qual", 4) == 0) {
            *end = str + 4;
            ks_clear(&res->s);
            if (ks_resize(&res->s, core->l_qseq + 1) < 0)
                return -1;
            memcpy(res->s.s, sam_bam_batch_record_qual(record),
                   core->l_qseq);
            res->s.l = core->l_qseq;
            res->is_str = 1;
            return 0;
        }
        break;

    case 'r':
        if (memcmp(str, "rlen", 4) == 0) {
            *end = str + 4;
            res->d = bam_view_cigar_rlen(record);
            return 0;
        } else if (memcmp(str, "rname", 5) == 0) {
            const char *rn = sam_hdr_tid2name(hb->h, core->tid);

            *end = str + 5;
            res->is_str = 1;
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "rnext", 5) == 0) {
            const char *rn = sam_hdr_tid2name(hb->h, core->mtid);

            *end = str + 5;
            res->is_str = 1;
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "refid", 5) == 0) {
            *end = str + 5;
            res->d = core->tid;
            return 0;
        }
        break;

    case 's':
        if (memcmp(str, "seq", 3) == 0) {
            *end = str + 3;
            ks_clear(&res->s);
            if (ks_resize(&res->s, core->l_qseq + 1) < 0)
                return -1;
            nibble2base((uint8_t *)sam_bam_batch_record_seq(record),
                        res->s.s, core->l_qseq);
            res->s.s[core->l_qseq] = 0;
            res->s.l = core->l_qseq;
            res->is_str = 1;
            return 0;
        } else if (memcmp(str, "sclen", 5) == 0) {
            *end = str + 5;
            res->d = bam_view_cigar_sclen(record);
            return 0;
        }
        break;

    case 't':
        if (memcmp(str, "tlen", 4) == 0) {
            *end = str + 4;
            res->d = core->isize;
            return 0;
        }
        break;

    case '[':
        if (*str == '[' && str[1] && str[2] && str[3] == ']') {
            sam_bam_record_view_t view;
            const uint8_t *aux;

            *end = str + 4;
            sam_bam_record_view_from_batch(&view, record);
            aux = sam_bam_record_view_aux_get(&view, str + 1);
            if (aux) {
                res->is_true = 1;
                switch (*aux) {
                case 'Z':
                case 'H':
                    res->is_str = 1;
                    kputs((const char *)aux + 1, ks_clear(&res->s));
                    break;

                case 'A':
                    res->is_str = 1;
                    kputsn((const char *)aux + 1, 1, ks_clear(&res->s));
                    break;

                case 'i': case 'I':
                case 's': case 'S':
                case 'c': case 'C':
                    res->is_str = 0;
                    res->d = bam_aux2i((uint8_t *)aux);
                    break;

                case 'f':
                case 'd':
                    res->is_str = 0;
                    res->d = bam_aux2f((uint8_t *)aux);
                    break;

                default:
                    hts_log_error("Aux type '%c not yet supported by filters",
                                  *aux);
                    return -1;
                }
                return 0;
            } else {
                res->is_str = 1;
                ks_clear(&res->s);
                res->d = 0;
                res->is_true = 0;
                return 0;
            }
        }
        break;
    }

    return -1;
}

enum {
    SAM_BAM_FILTER_PRED_NUM_CMP = 1,
    SAM_BAM_FILTER_PRED_FLAG_BITS,
    SAM_BAM_FILTER_PRED_AUX_EXISTS,
    SAM_BAM_FILTER_PRED_AUX_NUM_CMP,
    SAM_BAM_FILTER_PRED_AUX_STR_CMP,
    SAM_BAM_FILTER_PRED_LIBRARY_STR_CMP
};

enum {
    SAM_BAM_FILTER_FIELD_MAPQ = 1,
    SAM_BAM_FILTER_FIELD_FLAG,
    SAM_BAM_FILTER_FIELD_REFID,
    SAM_BAM_FILTER_FIELD_POS,
    SAM_BAM_FILTER_FIELD_MREFID,
    SAM_BAM_FILTER_FIELD_MPOS,
    SAM_BAM_FILTER_FIELD_TLEN,
    SAM_BAM_FILTER_FIELD_NCIGAR,
    SAM_BAM_FILTER_FIELD_ENDPOS,
    SAM_BAM_FILTER_FIELD_QLEN,
    SAM_BAM_FILTER_FIELD_RLEN,
    SAM_BAM_FILTER_FIELD_SCLEN,
    SAM_BAM_FILTER_FIELD_HCLEN
};

enum {
    SAM_BAM_FILTER_CMP_EQ = 1,
    SAM_BAM_FILTER_CMP_NE,
    SAM_BAM_FILTER_CMP_LT,
    SAM_BAM_FILTER_CMP_LE,
    SAM_BAM_FILTER_CMP_GT,
    SAM_BAM_FILTER_CMP_GE
};

static char *sam_filter_trimdup(const char *beg, const char *end)
{
    char *out;
    size_t len;

    while (beg < end && isspace((unsigned char)*beg))
        beg++;
    while (end > beg && isspace((unsigned char)end[-1]))
        end--;
    len = (size_t)(end - beg);
    out = malloc(len + 1);
    if (!out)
        return NULL;
    memcpy(out, beg, len);
    out[len] = '\0';
    return out;
}

static const char *sam_filter_find_top_and(const char *s)
{
    int depth = 0, quote = 0;

    for (; *s; s++) {
        if (quote) {
            if (*s == quote)
                quote = 0;
            else if (*s == '\\' && s[1])
                s++;
            continue;
        }
        if (*s == '"' || *s == '\'') {
            quote = *s;
        } else if (*s == '(') {
            depth++;
        } else if (*s == ')' && depth > 0) {
            depth--;
        } else if (depth == 0 && s[0] == '&' && s[1] == '&') {
            return s;
        }
    }
    return s;
}

static int sam_filter_matching_outer_parens(const char *s)
{
    size_t len = strlen(s), i;
    int depth = 0, quote = 0;

    if (len < 2 || s[0] != '(' || s[len - 1] != ')')
        return 0;
    for (i = 0; i < len; i++) {
        if (quote) {
            if (s[i] == quote)
                quote = 0;
            else if (s[i] == '\\' && i + 1 < len)
                i++;
            continue;
        }
        if (s[i] == '"' || s[i] == '\'') {
            quote = s[i];
        } else if (s[i] == '(') {
            depth++;
        } else if (s[i] == ')') {
            depth--;
            if (depth == 0 && i != len - 1)
                return 0;
        }
    }
    return depth == 0;
}

static char *sam_filter_strip_outer(char *term)
{
    char *s = term;

    while (*s && isspace((unsigned char)*s))
        s++;
    while (sam_filter_matching_outer_parens(s)) {
        size_t len = strlen(s);
        s[len - 1] = '\0';
        s++;
        while (*s && isspace((unsigned char)*s))
            s++;
    }
    return s;
}

static int sam_filter_parse_cmp(char **sp)
{
    char *s = *sp;

    while (*s && isspace((unsigned char)*s))
        s++;
    if (s[0] == '=' && s[1] == '=') {
        *sp = s + 2;
        return SAM_BAM_FILTER_CMP_EQ;
    }
    if (s[0] == '!' && s[1] == '=') {
        *sp = s + 2;
        return SAM_BAM_FILTER_CMP_NE;
    }
    if (s[0] == '<' && s[1] == '=') {
        *sp = s + 2;
        return SAM_BAM_FILTER_CMP_LE;
    }
    if (s[0] == '>' && s[1] == '=') {
        *sp = s + 2;
        return SAM_BAM_FILTER_CMP_GE;
    }
    if (*s == '<') {
        *sp = s + 1;
        return SAM_BAM_FILTER_CMP_LT;
    }
    if (*s == '>') {
        *sp = s + 1;
        return SAM_BAM_FILTER_CMP_GT;
    }
    return 0;
}

static int sam_filter_parse_number(char **sp, double *value)
{
    char *s = *sp, *end = NULL;

    while (*s && isspace((unsigned char)*s))
        s++;
    if ((s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) ||
        (s[0] == '-' && s[1] == '0' && (s[2] == 'x' || s[2] == 'X'))) {
        long long v = strtoll(s, &end, 0);
        if (end == s)
            return -1;
        *value = (double)v;
    } else {
        *value = strtod(s, &end);
        if (end == s)
            return -1;
    }
    while (*end && isspace((unsigned char)*end))
        end++;
    *sp = end;
    return 0;
}

static char *sam_filter_parse_string(char **sp)
{
    char *s = *sp, *out, *dst;
    size_t len = 0;
    int quote;

    while (*s && isspace((unsigned char)*s))
        s++;
    if (*s != '"' && *s != '\'')
        return NULL;
    quote = *s++;
    while (s[len] && s[len] != quote) {
        if (s[len] == '\\' && s[len + 1])
            len++;
        len++;
    }
    if (s[len] != quote)
        return NULL;
    out = malloc(len + 1);
    if (!out)
        return NULL;
    dst = out;
    while (*s && *s != quote) {
        if (*s == '\\' && s[1])
            s++;
        *dst++ = *s++;
    }
    *dst = '\0';
    if (*s == quote)
        s++;
    while (*s && isspace((unsigned char)*s))
        s++;
    *sp = s;
    return out;
}

static int sam_filter_parse_field(const char *name, int len, int *field,
                                  int *lazy)
{
    *lazy = 0;
    if (len == 4 && memcmp(name, "mapq", 4) == 0)
        *field = SAM_BAM_FILTER_FIELD_MAPQ;
    else if (len == 4 && memcmp(name, "flag", 4) == 0)
        *field = SAM_BAM_FILTER_FIELD_FLAG;
    else if (len == 5 && memcmp(name, "refid", 5) == 0)
        *field = SAM_BAM_FILTER_FIELD_REFID;
    else if (len == 3 && memcmp(name, "pos", 3) == 0)
        *field = SAM_BAM_FILTER_FIELD_POS;
    else if (len == 6 && memcmp(name, "mrefid", 6) == 0)
        *field = SAM_BAM_FILTER_FIELD_MREFID;
    else if ((len == 4 && memcmp(name, "mpos", 4) == 0) ||
             (len == 5 && memcmp(name, "pnext", 5) == 0))
        *field = SAM_BAM_FILTER_FIELD_MPOS;
    else if (len == 4 && memcmp(name, "tlen", 4) == 0)
        *field = SAM_BAM_FILTER_FIELD_TLEN;
    else if (len == 6 && memcmp(name, "ncigar", 6) == 0)
        *field = SAM_BAM_FILTER_FIELD_NCIGAR;
    else if (len == 6 && memcmp(name, "endpos", 6) == 0)
        *field = SAM_BAM_FILTER_FIELD_ENDPOS, *lazy = 1;
    else if (len == 4 && memcmp(name, "qlen", 4) == 0)
        *field = SAM_BAM_FILTER_FIELD_QLEN, *lazy = 1;
    else if (len == 4 && memcmp(name, "rlen", 4) == 0)
        *field = SAM_BAM_FILTER_FIELD_RLEN, *lazy = 1;
    else if (len == 5 && memcmp(name, "sclen", 5) == 0)
        *field = SAM_BAM_FILTER_FIELD_SCLEN, *lazy = 1;
    else if (len == 5 && memcmp(name, "hclen", 5) == 0)
        *field = SAM_BAM_FILTER_FIELD_HCLEN, *lazy = 1;
    else
        return 0;
    return 1;
}

static int sam_filter_trailing_empty(const char *s)
{
    while (*s && isspace((unsigned char)*s))
        s++;
    return *s == '\0';
}

static void sam_bam_filter_plan_raise_class(sam_bam_filter_plan_t *plan,
                                            sam_bam_filter_plan_class_t cls)
{
    if (plan->class_ < cls)
        plan->class_ = cls;
}

static int sam_filter_parse_term(sam_bam_filter_plan_t *plan, char *term)
{
    sam_bam_filter_plan_pred_t *pred;
    char *s = sam_filter_strip_outer(term), *name, *after;
    int negate = 0, cmp, field, lazy, len;
    double number;

    if (plan->n_predicates >= SAM_BAM_FILTER_PLAN_MAX_PREDICATES)
        return 0;
    if (*s == '!') {
        negate = 1;
        s++;
        s = sam_filter_strip_outer(s);
    }

    pred = &plan->predicates[plan->n_predicates];
    memset(pred, 0, sizeof(*pred));
    pred->negate = negate;

    if (*s == '[' && s[1] && s[2] && s[3] == ']') {
        pred->tag[0] = s[1];
        pred->tag[1] = s[2];
        s += 4;
        while (*s && isspace((unsigned char)*s))
            s++;
        if (*s == '\0') {
            pred->kind = SAM_BAM_FILTER_PRED_AUX_EXISTS;
            plan->n_predicates++;
            sam_bam_filter_plan_raise_class(plan,
                                            SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE);
            return 1;
        }
        cmp = sam_filter_parse_cmp(&s);
        if (!cmp)
            return 0;
        pred->cmp = cmp;
        while (*s && isspace((unsigned char)*s))
            s++;
        if (*s == '"' || *s == '\'') {
            char *parsed = sam_filter_parse_string(&s);
            if (!parsed || !sam_filter_trailing_empty(s)) {
                free(parsed);
                return 0;
            }
            pred->string = parsed;
            pred->kind = SAM_BAM_FILTER_PRED_AUX_STR_CMP;
        } else {
            if (sam_filter_parse_number(&s, &number) < 0 ||
                !sam_filter_trailing_empty(s))
                return 0;
            pred->number = number;
            pred->kind = SAM_BAM_FILTER_PRED_AUX_NUM_CMP;
        }
        plan->n_predicates++;
        sam_bam_filter_plan_raise_class(plan,
                                        SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE);
        return 1;
    }

    name = s;
    while (*s && (isalnum((unsigned char)*s) || *s == '_' || *s == '.'))
        s++;
    len = (int)(s - name);
    if (len <= 0)
        return 0;

    if (len == 7 && memcmp(name, "library", 7) == 0) {
        char *parsed;

        cmp = sam_filter_parse_cmp(&s);
        if (cmp != SAM_BAM_FILTER_CMP_EQ && cmp != SAM_BAM_FILTER_CMP_NE)
            return 0;
        parsed = sam_filter_parse_string(&s);
        if (!parsed || !sam_filter_trailing_empty(s)) {
            free(parsed);
            return 0;
        }
        pred->string = parsed;
        pred->kind = SAM_BAM_FILTER_PRED_LIBRARY_STR_CMP;
        pred->cmp = cmp;
        plan->n_predicates++;
        sam_bam_filter_plan_raise_class(plan,
                                        SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE);
        return 1;
    }

    if (!sam_filter_parse_field(name, len, &field, &lazy))
        return 0;
    pred->field = field;
    while (*s && isspace((unsigned char)*s))
        s++;
    if (field == SAM_BAM_FILTER_FIELD_FLAG && *s == '&') {
        s++;
        if (sam_filter_parse_number(&s, &number) < 0 ||
            !sam_filter_trailing_empty(s))
            return 0;
        pred->kind = SAM_BAM_FILTER_PRED_FLAG_BITS;
        pred->number = number;
        plan->n_predicates++;
        sam_bam_filter_plan_raise_class(plan,
                                        SAM_BAM_FILTER_PLAN_RAW_VIEW_SAFE);
        return 1;
    }

    after = s;
    cmp = sam_filter_parse_cmp(&after);
    if (!cmp)
        return 0;
    s = after;
    if (sam_filter_parse_number(&s, &number) < 0 ||
        !sam_filter_trailing_empty(s))
        return 0;
    pred->kind = SAM_BAM_FILTER_PRED_NUM_CMP;
    pred->cmp = cmp;
    pred->number = number;
    plan->n_predicates++;
    sam_bam_filter_plan_raise_class(plan, lazy
                                    ? SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE
                                    : SAM_BAM_FILTER_PLAN_RAW_VIEW_SAFE);
    return 1;
}

int sam_bam_filter_plan_init(sam_bam_filter_plan_t *plan, const char *expr)
{
    char *owned_expr = NULL, *stripped_expr = NULL;
    const char *s;

    if (!plan) {
        errno = EINVAL;
        return -1;
    }
    memset(plan, 0, sizeof(*plan));
    if (!expr)
        return 0;
    owned_expr = sam_filter_trimdup(expr, expr + strlen(expr));
    if (!owned_expr)
        return -1;
    stripped_expr = sam_filter_strip_outer(owned_expr);
    s = stripped_expr;
    while (*s) {
        const char *end = sam_filter_find_top_and(s);
        char *term = sam_filter_trimdup(s, end);
        int ok;

        if (!term) {
            sam_bam_filter_plan_destroy(plan);
            free(owned_expr);
            return -1;
        }
        ok = sam_filter_parse_term(plan, term);
        free(term);
        if (!ok) {
            sam_bam_filter_plan_destroy(plan);
            plan->class_ = SAM_BAM_FILTER_PLAN_MATERIALIZE_REQUIRED;
            free(owned_expr);
            return 0;
        }
        s = end;
        if (s[0] == '&' && s[1] == '&')
            s += 2;
        while (*s && isspace((unsigned char)*s))
            s++;
    }
    if (plan->n_predicates == 0)
        plan->class_ = SAM_BAM_FILTER_PLAN_MATERIALIZE_REQUIRED;
    free(owned_expr);
    return 0;
}

void sam_bam_filter_plan_destroy(sam_bam_filter_plan_t *plan)
{
    int i;

    if (!plan)
        return;
    for (i = 0; i < plan->n_predicates; i++)
        free(plan->predicates[i].string);
    memset(plan, 0, sizeof(*plan));
}

int sam_bam_filter_plan_is_usable(const sam_bam_filter_plan_t *plan)
{
    return plan && plan->n_predicates > 0 &&
           (plan->class_ == SAM_BAM_FILTER_PLAN_RAW_VIEW_SAFE ||
            plan->class_ == SAM_BAM_FILTER_PLAN_LAZY_RAW_VIEW_SAFE);
}

sam_bam_filter_plan_class_t sam_bam_filter_plan_class(
        const sam_bam_filter_plan_t *plan)
{
    return plan ? plan->class_ : SAM_BAM_FILTER_PLAN_UNUSABLE;
}

static int sam_filter_cmp_number(double lhs, int cmp, double rhs)
{
    switch (cmp) {
    case SAM_BAM_FILTER_CMP_EQ: return lhs == rhs;
    case SAM_BAM_FILTER_CMP_NE: return lhs != rhs;
    case SAM_BAM_FILTER_CMP_LT: return lhs < rhs;
    case SAM_BAM_FILTER_CMP_LE: return lhs <= rhs;
    case SAM_BAM_FILTER_CMP_GT: return lhs > rhs;
    case SAM_BAM_FILTER_CMP_GE: return lhs >= rhs;
    default: return 0;
    }
}

static int sam_filter_cmp_string(const char *lhs, int cmp, const char *rhs)
{
    int c;

    if (!lhs || !rhs)
        return 0;
    c = strcmp(lhs, rhs);
    switch (cmp) {
    case SAM_BAM_FILTER_CMP_EQ: return c == 0;
    case SAM_BAM_FILTER_CMP_NE: return c != 0;
    case SAM_BAM_FILTER_CMP_LT: return c < 0;
    case SAM_BAM_FILTER_CMP_LE: return c <= 0;
    case SAM_BAM_FILTER_CMP_GT: return c > 0;
    case SAM_BAM_FILTER_CMP_GE: return c >= 0;
    default: return 0;
    }
}

static int sam_filter_record_number(const sam_hdr_t *h,
                                    const bam_batch_record_t *record,
                                    bam1_t *scratch, int field,
                                    double *value)
{
    const bam1_core_t *core = &record->core;
    hts_pos_t len;

    (void)h;
    switch (field) {
    case SAM_BAM_FILTER_FIELD_MAPQ:
        *value = core->qual; return 0;
    case SAM_BAM_FILTER_FIELD_FLAG:
        *value = core->flag; return 0;
    case SAM_BAM_FILTER_FIELD_REFID:
        *value = core->tid; return 0;
    case SAM_BAM_FILTER_FIELD_POS:
        *value = core->pos + 1; return 0;
    case SAM_BAM_FILTER_FIELD_MREFID:
        *value = core->mtid; return 0;
    case SAM_BAM_FILTER_FIELD_MPOS:
        *value = core->mpos + 1; return 0;
    case SAM_BAM_FILTER_FIELD_TLEN:
        *value = core->isize; return 0;
    case SAM_BAM_FILTER_FIELD_NCIGAR:
        *value = core->n_cigar; return 0;
    case SAM_BAM_FILTER_FIELD_ENDPOS:
        if (sam_bam_batch_record_endpos(record, scratch, &len) < 0)
            return -1;
        *value = len; return 0;
    case SAM_BAM_FILTER_FIELD_QLEN:
        if (sam_bam_batch_record_query_len(record, scratch, 0, &len) < 0)
            return -1;
        *value = len; return 0;
    case SAM_BAM_FILTER_FIELD_RLEN:
        *value = bam_view_cigar_rlen(record); return 0;
    case SAM_BAM_FILTER_FIELD_SCLEN:
        *value = bam_view_cigar_sclen(record); return 0;
    case SAM_BAM_FILTER_FIELD_HCLEN:
        *value = bam_view_cigar_hclen(record); return 0;
    default:
        errno = EINVAL;
        return -1;
    }
}

static int sam_filter_aux_number(const uint8_t *aux, double *value)
{
    if (!aux)
        return 0;
    switch (*aux) {
    case 'i': case 'I': case 's': case 'S': case 'c': case 'C':
        *value = bam_aux2i((uint8_t *)aux);
        return 1;
    case 'f': case 'd':
        *value = bam_aux2f((uint8_t *)aux);
        return 1;
    default:
        return 0;
    }
}

static const char *sam_filter_aux_string(const uint8_t *aux, char tmp[2])
{
    if (!aux)
        return NULL;
    switch (*aux) {
    case 'Z':
    case 'H':
        return (const char *)aux + 1;
    case 'A':
        tmp[0] = (char)aux[1];
        tmp[1] = '\0';
        return tmp;
    default:
        return NULL;
    }
}

static int sam_bam_filter_plan_eval_pred(const sam_hdr_t *h,
                                         bam_batch_record_t *record,
                                         bam1_t *scratch,
                                         const sam_bam_filter_plan_pred_t *pred,
                                         int *out)
{
    sam_bam_record_view_t view;
    const uint8_t *aux;
    double value;
    char tmp[2];
    int pass = 0;

    sam_bam_record_view_from_batch(&view, record);
    switch (pred->kind) {
    case SAM_BAM_FILTER_PRED_NUM_CMP:
        if (sam_filter_record_number(h, record, scratch, pred->field,
                                     &value) < 0)
            return -1;
        pass = sam_filter_cmp_number(value, pred->cmp, pred->number);
        break;
    case SAM_BAM_FILTER_PRED_FLAG_BITS:
        pass = ((uint32_t)record->core.flag & (uint32_t)pred->number) != 0;
        break;
    case SAM_BAM_FILTER_PRED_AUX_EXISTS:
        aux = sam_bam_record_view_aux_get(&view, pred->tag);
        if (!aux && errno == EINVAL)
            return -1;
        pass = aux != NULL;
        break;
    case SAM_BAM_FILTER_PRED_AUX_NUM_CMP:
        aux = sam_bam_record_view_aux_get(&view, pred->tag);
        if (!aux) {
            if (errno == EINVAL)
                return -1;
            pass = 0;
            break;
        }
        pass = sam_filter_aux_number(aux, &value)
               ? sam_filter_cmp_number(value, pred->cmp, pred->number) : 0;
        break;
    case SAM_BAM_FILTER_PRED_AUX_STR_CMP:
        aux = sam_bam_record_view_aux_get(&view, pred->tag);
        if (!aux) {
            if (errno == EINVAL)
                return -1;
            pass = 0;
            break;
        }
        pass = sam_filter_cmp_string(sam_filter_aux_string(aux, tmp),
                                     pred->cmp, pred->string);
        break;
    case SAM_BAM_FILTER_PRED_LIBRARY_STR_CMP: {
        kstring_t lib = { 0, 0, NULL };
        int found = bam_view_get_library(h, record, &lib);

        if (found < 0) {
            free(lib.s);
            return -1;
        }
        pass = found ? sam_filter_cmp_string(lib.s, pred->cmp, pred->string)
                     : 0;
        free(lib.s);
        break;
    }
    default:
        errno = EINVAL;
        return -1;
    }
    if (pred->negate)
        pass = !pass;
    *out = pass;
    return 0;
}

int sam_bam_batch_record_passes_filter_plan(const sam_hdr_t *h,
                                            bam_batch_record_t *record,
                                            bam1_t *scratch,
                                            const sam_bam_filter_plan_t *plan,
                                            int *materialized)
{
    int local_materialized = 0, *mat;
    int i, pass;

    if (!sam_bam_filter_plan_is_usable(plan))
        return -2;
    mat = materialized ? materialized : &local_materialized;
    if (*mat)
        return -2;
    if (sam_bam_batch_record_validate_decode(record, scratch, mat) < 0)
        return -1;
    if (*mat)
        return -2;
    for (i = 0; i < plan->n_predicates; i++) {
        if (sam_bam_filter_plan_eval_pred(h, record, scratch,
                                          &plan->predicates[i], &pass) < 0)
            return -1;
        if (!pass)
            return 0;
    }
    return 1;
}

int sam_bam_batch_record_passes_filter(const sam_hdr_t *h,
                                       bam_batch_record_t *record,
                                       bam1_t *scratch,
                                       hts_filter_t *filt,
                                       int *materialized)
{
    hb_view_pair hb = {h, record, scratch};
    hts_expr_val_t res = HTS_EXPR_VAL_INIT;
    int local_materialized = 0;
    int *mat = materialized ? materialized : &local_materialized;
    int t;

    if (!filt)
        return 1;
    if (*mat) {
        if (!scratch) {
            errno = EINVAL;
            return -1;
        }
        return sam_passes_filter(h, scratch, filt);
    }
    if (sam_bam_batch_record_validate_decode(record, scratch, mat) < 0)
        return -1;
    if (*mat)
        return sam_passes_filter(h, scratch, filt);

    if (hts_filter_eval2(filt, &hb, bam_view_sym_lookup, &res)) {
        hts_log_error("Couldn't process filter expression");
        hts_expr_val_free(&res);
        return -1;
    }

    t = res.is_true;
    hts_expr_val_free(&res);
    return t;
}

static int cram_readrec(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = fpv;
    bam1_t *b = bv;
    int pass_filter, ret;

    do {
        ret = cram_get_bam_seq(fp->fp.cram, &b);
        if (ret < 0)
            return cram_eof(fp->fp.cram) ? -1 : -2;

        if (bam_tag2cigar(b, 1, 1) < 0)
            return -2;

        *tid = b->core.tid;
        *beg = b->core.pos;
        *end = bam_endpos(b);

        if (fp->filter) {
            pass_filter = sam_passes_filter(fp->bam_header, b, fp->filter);
            if (pass_filter < 0)
                return -2;
        } else {
            pass_filter = 1;
        }
    } while (pass_filter == 0);

    return ret;
}

static int cram_pseek(void *fp, int64_t offset, int whence)
{
    cram_fd *fd =  (cram_fd *)fp;

    if ((0 != cram_seek(fd, offset, SEEK_SET))
     && (0 != cram_seek(fd, offset - fd->first_container, SEEK_CUR)))
        return -1;

    fd->curr_position = offset;

    if (fd->ctr) {
        cram_free_container(fd->ctr);
        if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
            cram_free_container(fd->ctr_mt);

        fd->ctr = NULL;
        fd->ctr_mt = NULL;
        fd->ooc = 0;
    }

    return 0;
}

/*
 * cram_ptell is a pseudo-tell function, because it matches the position of the disk cursor only
 *   after a fresh seek call. Otherwise it indicates that the read takes place inside the buffered
 *   container previously fetched. It was designed like this to integrate with the functionality
 *   of the iterator stepping logic.
 */

static int64_t cram_ptell(void *fp)
{
    cram_fd *fd = (cram_fd *)fp;
    cram_container *c;
    cram_slice *s;
    int64_t ret = -1L;

    if (fd) {
        if ((c = fd->ctr) != NULL) {
            if ((s = c->slice) != NULL && s->max_rec) {
                if ((c->curr_slice + s->curr_rec/s->max_rec) >= (c->max_slice + 1))
                    fd->curr_position += c->offset + c->length;
            }
        }
        ret = fd->curr_position;
    }

    return ret;
}

static int bam_pseek(void *fp, int64_t offset, int whence)
{
    BGZF *fd = (BGZF *)fp;

    return bgzf_seek(fd, offset, whence);
}

static int64_t bam_ptell(void *fp)
{
    BGZF *fd = (BGZF *)fp;
    if (!fd)
        return -1L;

    return bgzf_tell(fd);
}



static hts_idx_t *index_load(htsFile *fp, const char *fn, const char *fnidx, int flags)
{
    switch (fp->format.format) {
    case bam:
    case sam:
        return hts_idx_load3(fn, fnidx, HTS_FMT_BAI, flags);

    case cram: {
        if (cram_index_load(fp->fp.cram, fn, fnidx) < 0) return NULL;

        // Cons up a fake "index" just pointing at the associated cram_fd:
        hts_cram_idx_t *idx = malloc(sizeof (hts_cram_idx_t));
        if (idx == NULL) return NULL;
        idx->fmt = HTS_FMT_CRAI;
        idx->cram = fp->fp.cram;
        return (hts_idx_t *) idx;
        }

    default:
        return NULL; // TODO Would use tbx_index_load if it returned hts_idx_t
    }
}

hts_idx_t *sam_index_load3(htsFile *fp, const char *fn, const char *fnidx, int flags)
{
    return index_load(fp, fn, fnidx, flags);
}

hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx) {
    return index_load(fp, fn, fnidx, HTS_IDX_SAVE_REMOTE);
}

hts_idx_t *sam_index_load(htsFile *fp, const char *fn)
{
    return index_load(fp, fn, NULL, HTS_IDX_SAVE_REMOTE);
}

static hts_itr_t *cram_itr_query(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end, hts_readrec_func *readrec)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    hts_itr_t *iter = (hts_itr_t *) calloc(1, sizeof(hts_itr_t));
    if (iter == NULL) return NULL;

    // Cons up a dummy iterator for which hts_itr_next() will simply invoke
    // the readrec function:
    iter->is_cram = 1;
    iter->read_rest = 1;
    iter->off = NULL;
    iter->bins.a = NULL;
    iter->readrec = readrec;

    if (tid >= 0 || tid == HTS_IDX_NOCOOR || tid == HTS_IDX_START) {
        cram_range r = { tid, beg+1, end };
        int ret = cram_set_option(cidx->cram, CRAM_OPT_RANGE, &r);

        iter->curr_off = 0;
        // The following fields are not required by hts_itr_next(), but are
        // filled in in case user code wants to look at them.
        iter->tid = tid;
        iter->beg = beg;
        iter->end = end;

        switch (ret) {
        case 0:
            break;

        case -2:
            // No data vs this ref, so mark iterator as completed.
            // Same as HTS_IDX_NONE.
            iter->finished = 1;
            break;

        default:
            free(iter);
            return NULL;
        }
    }
    else switch (tid) {
    case HTS_IDX_REST:
        iter->curr_off = 0;
        break;
    case HTS_IDX_NONE:
        iter->curr_off = 0;
        iter->finished = 1;
        break;
    default:
        hts_log_error("Query with tid=%d not implemented for CRAM files", tid);
        abort();
        break;
    }

    return iter;
}

hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    if (idx == NULL)
        return hts_itr_query(NULL, tid, beg, end, sam_readrec_rest);
    else if (cidx->fmt == HTS_FMT_CRAI)
        return cram_itr_query(idx, tid, beg, end, sam_readrec);
    else
        return hts_itr_query(idx, tid, beg, end, sam_readrec);
}

static int cram_name2id(void *fdv, const char *ref)
{
    cram_fd *fd = (cram_fd *) fdv;
    return sam_hdr_name2tid(fd->header, ref);
}

hts_itr_t *sam_itr_querys(const hts_idx_t *idx, sam_hdr_t *hdr, const char *region)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    return hts_itr_querys(idx, region, bam_name2id_wrapper, hdr,
                          cidx->fmt == HTS_FMT_CRAI ? cram_itr_query : hts_itr_query,
                          sam_readrec);
}

hts_itr_t *sam_itr_regarray(const hts_idx_t *idx, sam_hdr_t *hdr, char **regarray, unsigned int regcount)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    hts_reglist_t *r_list = NULL;
    int r_count = 0;

    if (!cidx || !hdr)
        return NULL;

    hts_itr_t *itr = NULL;
    if (cidx->fmt == HTS_FMT_CRAI) {
        r_list = hts_reglist_create(regarray, regcount, &r_count, cidx->cram, cram_name2id);
        if (!r_list)
            return NULL;
        itr = hts_itr_regions(idx, r_list, r_count, cram_name2id, cidx->cram,
                   hts_itr_multi_cram, cram_readrec, cram_pseek, cram_ptell);
    } else {
        r_list = hts_reglist_create(regarray, regcount, &r_count, hdr, bam_name2id_wrapper);
        if (!r_list)
            return NULL;
        itr = hts_itr_regions(idx, r_list, r_count, bam_name2id_wrapper, hdr,
                   hts_itr_multi_bam, sam_readrec, bam_pseek, bam_ptell);
    }

    if (!itr)
        hts_reglist_free(r_list, r_count);

    return itr;
}

hts_itr_t *sam_itr_regions(const hts_idx_t *idx, sam_hdr_t *hdr, hts_reglist_t *reglist, unsigned int regcount)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;

    if(!cidx || !hdr || !reglist)
        return NULL;

    if (cidx->fmt == HTS_FMT_CRAI)
        return hts_itr_regions(idx, reglist, regcount, cram_name2id, cidx->cram,
                   hts_itr_multi_cram, cram_readrec, cram_pseek, cram_ptell);
    else
        return hts_itr_regions(idx, reglist, regcount, bam_name2id_wrapper, hdr,
                   hts_itr_multi_bam, sam_readrec, bam_pseek, bam_ptell);
}

static int sam_reglist_interval_cmp(const void *av, const void *bv)
{
    const hts_pair_pos_t *a = (const hts_pair_pos_t *)av;
    const hts_pair_pos_t *b = (const hts_pair_pos_t *)bv;

    if (a->beg < b->beg) return -1;
    if (a->beg > b->beg) return 1;
    if (a->end < b->end) return -1;
    if (a->end > b->end) return 1;
    return 0;
}

hts_reglist_t *sam_reglist_dup_merged(const hts_reglist_t *reglist,
                                      int count, int *out_count)
{
    hts_reglist_t *out = NULL;
    int i, used = 0;

    if (out_count)
        *out_count = 0;
    if (!reglist || count <= 0 || !out_count)
        return NULL;

    out = (hts_reglist_t *)calloc((size_t)count, sizeof(*out));
    if (!out)
        return NULL;

    for (i = 0; i < count; i++) {
        const hts_reglist_t *src = &reglist[i];
        hts_reglist_t *dst;
        uint32_t j, n = 0;

        if (src->count == 0)
            continue;
        if (!src->intervals)
            goto fail;

        dst = &out[used];
        dst->reg = src->reg;
        dst->tid = src->tid;
        dst->intervals = (hts_pair_pos_t *)malloc((size_t)src->count *
                                                  sizeof(*dst->intervals));
        if (!dst->intervals)
            goto fail;
        memcpy(dst->intervals, src->intervals,
               (size_t)src->count * sizeof(*dst->intervals));
        qsort(dst->intervals, src->count, sizeof(*dst->intervals),
              sam_reglist_interval_cmp);

        for (j = 0; j < src->count; j++) {
            hts_pair_pos_t *cur = &dst->intervals[j];

            if (n > 0 && dst->intervals[n - 1].end >= cur->beg) {
                if (dst->intervals[n - 1].end < cur->end)
                    dst->intervals[n - 1].end = cur->end;
                continue;
            }
            if (n != j)
                dst->intervals[n] = *cur;
            n++;
        }

        dst->count = n;
        dst->min_beg = dst->intervals[0].beg;
        dst->max_end = dst->intervals[n - 1].end;
        used++;
    }

    if (used == 0)
        goto fail;
    *out_count = used;
    return out;

fail:
    hts_reglist_free(out, count);
    return NULL;
}

int sam_itr_next_filtered(htsFile *fp, hts_itr_t *iter, sam_hdr_t *h,
                          bam1_t *record, sam_region_overlap_f filter,
                          void *filter_data)
{
    int ret;

    if (!filter)
        return sam_itr_next(fp, iter, record);

    while ((ret = sam_itr_next(fp, iter, record)) >= 0) {
        hts_pos_t end;
        int keep;

        if (record->core.tid < 0)
            continue;
        end = bam_endpos(record);
        keep = filter(filter_data, h, record->core.tid, record->core.pos,
                      end);
        if (keep < 0)
            return -2;
        if (keep)
            return ret;
    }
    return ret;
}

int sam_bam_itr_next_batch_filtered(htsFile *fp, hts_itr_t *iter,
                                    sam_hdr_t *h, bam_batch_t *batch,
                                    sam_region_overlap_f filter,
                                    void *filter_data, bam1_t *scratch)
{
    int ret;

    if (!filter)
        return sam_bam_itr_next_batch(fp, iter, h, batch);

    while ((ret = sam_bam_itr_next_batch(fp, iter, h, batch)) >= 0) {
        int i, keep = 0;

        for (i = 0; i < batch->n_records; i++) {
            bam_batch_record_t *record = &batch->records[i];
            hts_pos_t end;
            int pass;

            if (record->core.tid < 0)
                continue;
            if (sam_bam_batch_record_endpos(record, scratch, &end) < 0) {
                sam_bam_batch_destroy(batch);
                return -2;
            }
            pass = filter(filter_data, h, record->core.tid,
                          record->core.pos, end);
            if (pass < 0) {
                sam_bam_batch_destroy(batch);
                return -2;
            }
            if (pass) {
                if (keep != i)
                    batch->records[keep] = *record;
                keep++;
            }
        }

        batch->n_records = keep;
        batch->data = NULL;
        batch->len = 0;
        batch->n_segments = 0;
        batch->segments = NULL;
        if (keep > 0)
            return keep;
        sam_bam_batch_destroy(batch);
    }
    return ret;
}

/**********************
 *** SAM header I/O ***
 **********************/

#include "htslib/kseq.h"
#include "htslib/kstring.h"

sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text)
{
    sam_hdr_t *bh = sam_hdr_init();
    if (!bh) return NULL;

    if (sam_hdr_add_lines(bh, text, l_text) != 0) {
        sam_hdr_destroy(bh);
        return NULL;
    }

    return bh;
}

// Minimal sanitisation of a header to ensure.
// - null terminated string.
// - all lines start with @ (also implies no blank lines).
//
// Much more could be done, but currently is not, including:
// - checking header types are known (HD, SQ, etc).
// - syntax (eg checking tab separated fields).
// - validating n_targets matches @SQ records.
// - validating target lengths against @SQ records.
static sam_hdr_t *sam_hdr_sanitise(sam_hdr_t *h) {
    if (!h)
        return NULL;

    // Special case for empty headers.
    if (h->l_text == 0)
        return h;

    size_t i;
    unsigned int lnum = 0;
    char *cp = h->text, last = '\n';
    for (i = 0; i < h->l_text; i++) {
        // NB: l_text excludes terminating nul.  This finds early ones.
        if (cp[i] == 0)
            break;

        // Error on \n[^@], including duplicate newlines
        if (last == '\n') {
            lnum++;
            if (cp[i] != '@') {
                hts_log_error("Malformed SAM header at line %u", lnum);
                sam_hdr_destroy(h);
                return NULL;
            }
        }

        last = cp[i];
    }

    if (i < h->l_text) { // Early nul found.  Complain if not just padding.
        size_t j = i;
        while (j < h->l_text && cp[j] == '\0') j++;
        if (j < h->l_text)
            hts_log_warning("Unexpected NUL character in header. Possibly truncated");
    }

    // Add trailing newline and/or trailing nul if required.
    if (last != '\n') {
        hts_log_warning("Missing trailing newline on SAM header. Possibly truncated");

        if (h->l_text < 2 || i >= h->l_text - 2) {
            if (h->l_text >= SIZE_MAX - 2) {
                hts_log_error("No room for extra newline");
                sam_hdr_destroy(h);
                return NULL;
            }

            cp = realloc(h->text, (size_t) h->l_text+2);
            if (!cp) {
                sam_hdr_destroy(h);
                return NULL;
            }
            h->text = cp;
        }
        cp[i++] = '\n';

        // l_text may be larger already due to multiple nul padding
        if (h->l_text < i)
            h->l_text = i;
        cp[h->l_text] = '\0';
    }

    return h;
}

static sam_hdr_t *sam_hdr_create(htsFile* fp) {
    sam_hdr_t* h = sam_hdr_init();
    if (!h)
        return NULL;

    if (sam_hdr_build_from_sam_file(h, fp) != 0) {
        sam_hdr_destroy(h);
        return NULL;
    }

    if (fp->bam_header)
        sam_hdr_destroy(fp->bam_header);
    fp->bam_header = sam_hdr_sanitise(h);
    fp->bam_header->ref_count = 1;

    return fp->bam_header;
}

sam_hdr_t *sam_hdr_read(htsFile *fp)
{
    sam_hdr_t *h = NULL;
    if (!fp) {
        errno = EINVAL;
        return NULL;
    }

    switch (fp->format.format) {
    case bam:
        h = sam_hdr_sanitise(bam_hdr_read(fp->fp.bgzf));
        break;

    case cram:
        h = sam_hdr_sanitise(sam_hdr_dup(fp->fp.cram->header));
        break;

    case sam:
        h = sam_hdr_create(fp);
        break;

    case fastq_format:
    case fasta_format:
        return sam_hdr_init();

    case empty_format:
        errno = EPIPE;
        return NULL;

    default:
        errno = EFTYPE;
        return NULL;
    }
    //only sam,bam and cram reaches here
    if (h && !fp->bam_header) { //set except for sam which already has it
        //for cram, it is the o/p header as for rest and not the internal header
        fp->bam_header = h;
        sam_hdr_incr_ref(fp->bam_header);
    }
    return h;
}

int sam_hdr_write(htsFile *fp, const sam_hdr_t *h)
{
    if (!fp || !h) {
        errno = EINVAL;
        return -1;
    }

    switch (fp->format.format) {
    case binary_format:
        fp->format.category = sequence_data;
        fp->format.format = bam;
        /* fall-through */
    case bam:
        if (bam_hdr_write(fp->fp.bgzf, h) < 0) return -1;
        break;

    case cram: {
        cram_fd *fd = fp->fp.cram;
        if (cram_set_header2(fd, h) < 0) return -1;
        if (fp->fn_aux)
            cram_load_reference(fd, fp->fn_aux);
        if (cram_write_SAM_hdr(fd, fd->header) < 0) return -1;
        }
        break;

    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam: {
        if (!h->hrecs && !h->text)
            return 0;
        char *text;
        kstring_t hdr_ks = { 0, 0, NULL };
        size_t l_text;
        ssize_t bytes;
        int r = 0, no_sq = 0;

        if (h->hrecs) {
            if (sam_hrecs_rebuild_text(h->hrecs, &hdr_ks) != 0)
                return -1;
            text = hdr_ks.s;
            l_text = hdr_ks.l;
        } else {
            const char *p = NULL;
            do {
                const char *q = p == NULL ? h->text : p + 4;
                p = strstr(q, "@SQ\t");
            } while (!(p == NULL || p == h->text || *(p - 1) == '\n'));
            no_sq = p == NULL;
            text = h->text;
            l_text = h->l_text;
        }

        if (fp->is_bgzf) {
            bytes = bgzf_write(fp->fp.bgzf, text, l_text);
        } else {
            bytes = hwrite(fp->fp.hfile, text, l_text);
        }
        free(hdr_ks.s);
        if (bytes != l_text)
            return -1;

        if (no_sq) {
            int i;
            for (i = 0; i < h->n_targets; ++i) {
                fp->line.l = 0;
                r |= kputsn("@SQ\tSN:", 7, &fp->line) < 0;
                r |= kputs(h->target_name[i], &fp->line) < 0;
                r |= kputsn("\tLN:", 4, &fp->line) < 0;
                r |= kputw(h->target_len[i], &fp->line) < 0;
                r |= kputc('\n', &fp->line) < 0;
                if (r != 0)
                    return -1;

                if (fp->is_bgzf) {
                    bytes = bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l);
                } else {
                    bytes = hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
                }
                if (bytes != fp->line.l)
                    return -1;
            }
        }
        if (fp->is_bgzf) {
            if (bgzf_flush(fp->fp.bgzf) != 0) return -1;
        } else {
            if (hflush(fp->fp.hfile) != 0) return -1;
        }
        }
        break;

    case fastq_format:
    case fasta_format:
        // Nothing to output; FASTQ has no file headers.
        return 0;
        break;

    default:
        errno = EBADF;
        return -1;
    }
    //only sam,bam and cram reaches here
    if (h) {    //the new header
        sam_hdr_t *tmp = fp->bam_header;
        fp->bam_header = sam_hdr_dup(h);
        sam_hdr_destroy(tmp);
        if (!fp->bam_header && h)
            return -1;  //failed to duplicate
    }
    return 0;
}

static int old_sam_hdr_change_HD(sam_hdr_t *h, const char *key, const char *val)
{
    char *p, *q, *beg = NULL, *end = NULL, *newtext;
    size_t new_l_text;
    if (!h || !key)
        return -1;

    if (h->l_text > 3) {
        if (strncmp(h->text, "@HD", 3) == 0) { //@HD line exists
            if ((p = strchr(h->text, '\n')) == 0) return -1;
            *p = '\0'; // for strstr call

            char tmp[5] = { '\t', key[0], key[0] ? key[1] : '\0', ':', '\0' };

            if ((q = strstr(h->text, tmp)) != 0) { // key exists
                *p = '\n'; // change back

                // mark the key:val
                beg = q;
                for (q += 4; *q != '\n' && *q != '\t'; ++q);
                end = q;

                if (val && (strncmp(beg + 4, val, end - beg - 4) == 0)
                    && strlen(val) == end - beg - 4)
                     return 0; // val is the same, no need to change

            } else {
                beg = end = p;
                *p = '\n';
            }
        }
    }
    if (beg == NULL) { // no @HD
        new_l_text = h->l_text;
        if (new_l_text > SIZE_MAX - strlen(SAM_FORMAT_VERSION) - 9)
            return -1;
        new_l_text += strlen(SAM_FORMAT_VERSION) + 8;
        if (val) {
            if (new_l_text > SIZE_MAX - strlen(val) - 5)
                return -1;
            new_l_text += strlen(val) + 4;
        }
        newtext = (char*)malloc(new_l_text + 1);
        if (!newtext) return -1;

        if (val)
            snprintf(newtext, new_l_text + 1,
                    "@HD\tVN:%s\t%s:%s\n%s", SAM_FORMAT_VERSION, key, val, h->text);
        else
            snprintf(newtext, new_l_text + 1,
                    "@HD\tVN:%s\n%s", SAM_FORMAT_VERSION, h->text);
    } else { // has @HD but different or no key
        new_l_text = (beg - h->text) + (h->text + h->l_text - end);
        if (val) {
            if (new_l_text > SIZE_MAX - strlen(val) - 5)
                return -1;
            new_l_text += strlen(val) + 4;
        }
        newtext = (char*)malloc(new_l_text + 1);
        if (!newtext) return -1;

        if (val) {
            snprintf(newtext, new_l_text + 1, "%.*s\t%s:%s%s",
                    (int) (beg - h->text), h->text, key, val, end);
        } else { //delete key
            snprintf(newtext, new_l_text + 1, "%.*s%s",
                    (int) (beg - h->text), h->text, end);
        }
    }
    free(h->text);
    h->text = newtext;
    h->l_text = new_l_text;
    return 0;
}


int sam_hdr_change_HD(sam_hdr_t *h, const char *key, const char *val)
{
    if (!h || !key)
        return -1;

    if (!h->hrecs)
        return old_sam_hdr_change_HD(h, key, val);

    if (val) {
        if (sam_hdr_update_line(h, "HD", NULL, NULL, key, val, NULL) != 0)
            return -1;
    } else {
        if (sam_hdr_remove_tag_id(h, "HD", NULL, NULL, key) != 0)
            return -1;
    }
    return sam_hdr_rebuild(h);
}

/* releases existing header and sets new one; increments ref count if not
duplicating */
int sam_hdr_set(samFile *fp, sam_hdr_t *h, int duplicate)
{
    if (!fp)
        return -1;

    if (duplicate) {
        sam_hdr_t *tmp = fp->bam_header;
        fp->bam_header = sam_hdr_dup(h);
        sam_hdr_destroy(tmp);
        if (!fp->bam_header && h)
            return -1;  //duplicate failed
    } else {
        if (fp->bam_header != h) {  //if not the same
            sam_hdr_destroy(fp->bam_header);
            fp->bam_header = h;
            sam_hdr_incr_ref(fp->bam_header);
        }
    }

    return 0;
}

//return the bam_header, user has to use sam_hdr_incr_ref where ever required
sam_hdr_t* sam_hdr_get(samFile* fp)
{
    if (!fp)
        return NULL;
    return fp->bam_header;
}

/**********************
 *** SAM record I/O ***
 **********************/

// The speed of this code can vary considerably depending on minor code
// changes elsewhere as some of the tight loops are particularly prone to
// speed changes when the instruction blocks are split over a 32-byte
// boundary.  To protect against this, we explicitly specify an alignment
// for this function.  If this is insufficient, we may also wish to
// consider alignment of blocks within this function via
// __attribute__((optimize("align-loops=5"))) (gcc) or clang equivalents.
// However it's not very portable.
// Instead we break into separate functions so we can explicitly specify
// use __attribute__((aligned(32))) instead and force consistent loop
// alignment.
static inline int64_t grow_B_array(bam1_t *b, uint32_t *n, size_t size) {
    // Avoid overflow on 32-bit platforms, but it breaks BAM anyway
    if (*n > INT32_MAX*0.666) {
        errno = ENOMEM;
        return -1;
    }

    size_t bytes = (size_t)size * (size_t)(*n>>1);
    if (possibly_expand_bam_data(b, bytes) < 0) {
        hts_log_error("Out of memory");
        return -1;
    }

    (*n)+=*n>>1;
    return 0;
}


// This ensures that q always ends up at the next comma after
// reading a number even if it's followed by junk.  It
// prevents the possibility of trying to read more than n items.
#define skip_to_comma_(q) do { while (*(q) > '\t' && *(q) != ',') (q)++; } while (0)

HTS_ALIGN32
static char *sam_parse_Bc_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        *(b->data + b->l_data) = hts_str2int(q + 1, &q, 8, overflow);
        b->l_data++;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BC_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            *(b->data + b->l_data) = hts_str2uint(q + 1, &q, 8, overflow);
            b->l_data++;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bs_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        i16_to_le(hts_str2int(q + 1, &q, 16, overflow),
                  b->data + b->l_data);
        b->l_data += 2;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BS_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u16_to_le(hts_str2uint(q + 1, &q, 16, overflow),
                      b->data + b->l_data);
            b->l_data += 2;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bi_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        i32_to_le(hts_str2int(q + 1, &q, 32, overflow),
                  b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BI_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u32_to_le(hts_str2uint(q + 1, &q, 32, overflow),
                      b->data + b->l_data);
            b->l_data += 4;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bf_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        float_to_le(strtod(q + 1, &q), b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
}

HTS_ALIGN32
static int sam_parse_B_vals_r(char type, uint32_t nalloc, char *in,
                              char **end, bam1_t *b,
                              int *ctr) {
    // Protect against infinite recursion when dealing with invalid input.
    // An example string is "XX:B:C,-".  The lack of a number means min=0,
    // but it overflowed due to "-" and so we repeat ad-infinitum.
    //
    // Loop detection is the safest solution incase there are other
    // strange corner cases with malformed inputs.
    if (++(*ctr) > 2) {
        hts_log_error("Malformed data in B:%c array", type);
        return -1;
    }

    int orig_l = b->l_data;
    char *q = in;
    int32_t size;
    size_t bytes;
    int overflow = 0;

    size = aux_type2size(type);
    if (size <= 0 || size > 4) {
        hts_log_error("Unrecognized type B:%c", type);
        return -1;
    }

    // Ensure space for type + values.
    // The first pass through here we don't know the number of entries and
    // nalloc == 0.  We start with a small working set and then parse the
    // data, growing as needed.
    //
    // If we have a second pass through we do know the number of entries
    // and nalloc is already known.  We have no need to expand the bam data.
    if (!nalloc)
         nalloc=7;

    // Ensure allocated memory is big enough (for current nalloc estimate)
    bytes = (size_t) nalloc * (size_t) size;
    if (bytes / size != nalloc
        || possibly_expand_bam_data(b, bytes + 2 + sizeof(uint32_t))) {
        hts_log_error("Out of memory");
        return -1;
    }

    uint32_t nused = 0;

    b->data[b->l_data++] = 'B';
    b->data[b->l_data++] = type;
    // 32-bit B-array length is inserted later once we know it.
    int b_len_idx = b->l_data;
    b->l_data += sizeof(uint32_t);

    if (type == 'c') {
        if (!(q = sam_parse_Bc_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'C') {
        if (!(q = sam_parse_BC_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 's') {
        if (!(q = sam_parse_Bs_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'S') {
        if (!(q = sam_parse_BS_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'i') {
        if (!(q = sam_parse_Bi_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'I') {
        if (!(q = sam_parse_BI_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'f') {
        if (!(q = sam_parse_Bf_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    }
    if (*q != '\t' && *q != '\0') {
        // Unknown B array type or junk in the numbers
        hts_log_error("Malformed B:%c", type);
        return -1;
    }
    i32_to_le(nused, b->data + b_len_idx);

    if (!overflow) {
        *end = q;
        return 0;
    } else {
        int64_t max = 0, min = 0, val;
        // Given type was incorrect.  Try to rescue the situation.
        char *r = q;
        q = in;
        overflow = 0;
        b->l_data = orig_l;
        // Find out what range of values is present
        while (q < r) {
            val = hts_str2int(q + 1, &q, 64, &overflow);
            if (max < val) max = val;
            if (min > val) min = val;
            skip_to_comma_(q);
        }
        // Retry with appropriate type
        if (!overflow) {
            if (min < 0) {
                if (min >= INT8_MIN && max <= INT8_MAX) {
                    return sam_parse_B_vals_r('c', nalloc, in, end, b, ctr);
                } else if (min >= INT16_MIN && max <= INT16_MAX) {
                    return sam_parse_B_vals_r('s', nalloc, in, end, b, ctr);
                } else if (min >= INT32_MIN && max <= INT32_MAX) {
                    return sam_parse_B_vals_r('i', nalloc, in, end, b, ctr);
                }
            } else {
                if (max < UINT8_MAX) {
                    return sam_parse_B_vals_r('C', nalloc, in, end, b, ctr);
                } else if (max <= UINT16_MAX) {
                    return sam_parse_B_vals_r('S', nalloc, in, end, b, ctr);
                } else if (max <= UINT32_MAX) {
                    return sam_parse_B_vals_r('I', nalloc, in, end, b, ctr);
                }
            }
        }
        // If here then at least one of the values is too big to store
        hts_log_error("Numeric value in B array out of allowed range");
        return -1;
    }
#undef skip_to_comma_
}

HTS_ALIGN32
static int sam_parse_B_vals(char type, char *in, char **end, bam1_t *b)
{
    int ctr = 0;
    uint32_t nalloc = 0;
    return sam_parse_B_vals_r(type, nalloc, in, end, b, &ctr);
}

static inline unsigned int parse_sam_flag(char *v, char **rv, int *overflow) {
    if (*v >= '1' && *v <= '9') {
        return hts_str2uint(v, rv, 16, overflow);
    }
    else if (*v == '0') {
        // handle single-digit "0" directly; otherwise it's hex or octal
        if (v[1] == '\t') { *rv = v+1; return 0; }
        else {
            unsigned long val = strtoul(v, rv, 0);
            if (val > 65535) { *overflow = 1; return 65535; }
            return val;
        }
    }
    else {
        // TODO implement symbolic flag letters
        *rv = v;
        return 0;
    }
}

// Parse tag line and append to bam object b.
// Shared by both SAM and FASTQ parsers.
//
// The difference between the two is how lenient we are to recognising
// non-compliant strings.  The FASTQ parser glosses over arbitrary
// non-SAM looking strings.
static inline int aux_parse(char *start, char *end, bam1_t *b, int lenient,
                            khash_t(tag) *tag_whitelist) {
    int overflow = 0;
    int checkpoint;
    char logbuf[40];
    char *q = start, *p = end;

#define _parse_err(cond, ...)                   \
    do {                                        \
        if (cond) {                             \
            if (lenient) {                      \
                while (q < p && !isspace_c(*q))   \
                    q++;                        \
                while (q < p && isspace_c(*q))    \
                    q++;                        \
                b->l_data = checkpoint;         \
                goto loop;                      \
            } else {                            \
                hts_log_error(__VA_ARGS__);     \
                goto err_ret;                   \
            }                                   \
        }                                       \
    } while (0)

    while (q < p) loop: {
        char type;
        checkpoint = b->l_data;
        if (p - q < 5) {
            if (lenient) {
                break;
            } else {
                hts_log_error("Incomplete aux field");
                goto err_ret;
            }
        }
        _parse_err(q[0] < '!' || q[1] < '!', "invalid aux tag id");

        if (lenient && (q[2] | q[4]) != ':') {
            while (q < p && !isspace_c(*q))
                q++;
            while (q < p && isspace_c(*q))
                q++;
            continue;
        }

        if (tag_whitelist) {
            int tt = q[0]*256 + q[1];
            if (kh_get(tag, tag_whitelist, tt) == kh_end(tag_whitelist)) {
                while (q < p && *q != '\t')
                    q++;
                continue;
            }
        }

        // Copy over id
        if (possibly_expand_bam_data(b, 2) < 0) goto err_ret;
        memcpy(b->data + b->l_data, q, 2); b->l_data += 2;
        q += 3; type = *q++; ++q; // q points to value
        if (type != 'Z' && type != 'H') // the only zero length acceptable fields
            _parse_err(*q <= '\t', "incomplete aux field");

        // Ensure enough space for a double + type allocated.
        if (possibly_expand_bam_data(b, 16) < 0) goto err_ret;

        if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
            b->data[b->l_data++] = 'A';
            b->data[b->l_data++] = *q++;
        } else if (type == 'i' || type == 'I') {
            if (*q == '-') {
                int32_t x = hts_str2int(q, &q, 32, &overflow);
                if (x >= INT8_MIN) {
                    b->data[b->l_data++] = 'c';
                    b->data[b->l_data++] = x;
                } else if (x >= INT16_MIN) {
                    b->data[b->l_data++] = 's';
                    i16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'i';
                    i32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            } else {
                uint32_t x = hts_str2uint(q, &q, 32, &overflow);
                if (x <= UINT8_MAX) {
                    b->data[b->l_data++] = 'C';
                    b->data[b->l_data++] = x;
                } else if (x <= UINT16_MAX) {
                    b->data[b->l_data++] = 'S';
                    u16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'I';
                    u32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            }
        } else if (type == 'f') {
            b->data[b->l_data++] = 'f';
            float_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(float);
        } else if (type == 'd') {
            b->data[b->l_data++] = 'd';
            double_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(double);
        } else if (type == 'Z' || type == 'H') {
            char *end = strchr(q, '\t');
            if (!end) end = q + strlen(q);
            _parse_err(type == 'H' && ((end-q)&1) != 0,
                       "hex field does not have an even number of digits");
            b->data[b->l_data++] = type;
            if (possibly_expand_bam_data(b, end - q + 1) < 0) goto err_ret;
            memcpy(b->data + b->l_data, q, end - q);
            b->l_data += end - q;
            b->data[b->l_data++] = '\0';
            q = end;
        } else if (type == 'B') {
            type = *q++; // q points to the first ',' following the typing byte
            _parse_err(*q && *q != ',' && *q != '\t',
                       "B aux field type not followed by ','");

            if (sam_parse_B_vals(type, q, &q, b) < 0)
                goto err_ret;
        } else _parse_err(1, "unrecognized type %s", hts_strprint(logbuf, sizeof logbuf, '\'', &type, 1));

        while (*q > '\t') { q++; } // Skip any junk to next tab
        q++;
    }

    _parse_err(!lenient && overflow != 0, "numeric value out of allowed range");
#undef _parse_err

    return 0;

err_ret:
    return -2;
}

int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret; *tab = '\0'; (_p) = tab + 1; } while (0)

#if HTS_ALLOW_UNALIGNED != 0 && ULONG_MAX == 0xffffffffffffffff

// Macro that operates on 64-bits at a time.
#define COPY_MINUS_N(to,from,n,l,failed)                        \
    do {                                                        \
        uint64_u *from8 = (uint64_u *)(from);                   \
        uint64_u *to8 = (uint64_u *)(to);                       \
        uint64_t uflow = 0;                                     \
        size_t l8 = (l)>>3, i;                                  \
        for (i = 0; i < l8; i++) {                              \
            to8[i] = from8[i] - (n)*0x0101010101010101UL;       \
            uflow |= to8[i];                                    \
        }                                                       \
        for (i<<=3; i < (l); ++i) {                             \
            to[i] = from[i] - (n);                              \
            uflow |= to[i];                                     \
        }                                                       \
        failed = (uflow & 0x8080808080808080UL) > 0;            \
    } while (0)

#else

// Basic version which operates a byte at a time
#define COPY_MINUS_N(to,from,n,l,failed) do {                \
        uint8_t uflow = 0;                                   \
        for (i = 0; i < (l); ++i) {                          \
            (to)[i] = (from)[i] - (n);                       \
            uflow |= (uint8_t) (to)[i];                      \
        }                                                    \
        failed = (uflow & 0x80) > 0;                         \
    } while (0)

#endif

#define _get_mem(type_t, x, b, l) if (possibly_expand_bam_data((b), (l)) < 0) goto err_ret; *(x) = (type_t*)((b)->data + (b)->l_data); (b)->l_data += (l)
#define _parse_err(cond, ...) do { if (cond) { hts_log_error(__VA_ARGS__); goto err_ret; } } while (0)
#define _parse_warn(cond, ...) do { if (cond) { hts_log_warning(__VA_ARGS__); } } while (0)

    uint8_t *t;

    char *p = s->s, *q;
    int i, overflow = 0;
    char logbuf[40];
    hts_pos_t cigreflen;
    bam1_core_t *c = &b->core;

    b->l_data = 0;
    memset(c, 0, 32);

    // qname
    q = _read_token(p);

    _parse_warn(p - q <= 1, "empty query name");
    _parse_err(p - q > 255, "query name too long");
    // resize large enough for name + extranul
    if (possibly_expand_bam_data(b, (p - q) + 4) < 0) goto err_ret;
    memcpy(b->data + b->l_data, q, p-q); b->l_data += p-q;

    c->l_extranul = (4 - (b->l_data & 3)) & 3;
    memcpy(b->data + b->l_data, "\0\0\0\0", c->l_extranul);
    b->l_data += c->l_extranul;

    c->l_qname = p - q + c->l_extranul;

    // flag
    c->flag = parse_sam_flag(p, &p, &overflow);
    if (*p++ != '\t') goto err_ret; // malformated flag

    // chr
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(h->n_targets == 0, "no SQ lines present in the header");
        c->tid = bam_name2id(h, q);
        _parse_err(c->tid < -1, "failed to parse header");
        _parse_warn(c->tid < 0, "unrecognized reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    } else c->tid = -1;

    // pos
    c->pos = hts_str2uint(p, &p, 62, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->pos < 0 && c->tid >= 0) {
        _parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
        c->tid = -1;
    }
    if (c->tid < 0) c->flag |= BAM_FUNMAP;

    // mapq
    c->qual = hts_str2uint(p, &p, 8, &overflow);
    if (*p++ != '\t') goto err_ret;
    // cigar
    if (*p != '*') {
        uint32_t *cigar = NULL;
        int old_l_data = b->l_data;
        int n_cigar = bam_parse_cigar(p, &p, b);
        if (n_cigar < 1 || *p++ != '\t') goto err_ret;
        cigar = (uint32_t *)(b->data + old_l_data);

        // can't use bam_endpos() directly as some fields not yet set up
        cigreflen = (!(c->flag&BAM_FUNMAP))? bam_cigar2rlen(c->n_cigar, cigar) : 1;
        if (cigreflen == 0) cigreflen = 1;
    } else {
        _parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
        c->flag |= BAM_FUNMAP;
        q = _read_token(p);
        cigreflen = 1;
    }
    _parse_err(HTS_POS_MAX - cigreflen <= c->pos,
               "read ends beyond highest supported position");
    c->bin = hts_reg2bin(c->pos, c->pos + cigreflen, 14, 5);
    // mate chr
    q = _read_token(p);
    if (strcmp(q, "=") == 0) {
        c->mtid = c->tid;
    } else if (strcmp(q, "*") == 0) {
        c->mtid = -1;
    } else {
        c->mtid = bam_name2id(h, q);
        _parse_err(c->mtid < -1, "failed to parse header");
        _parse_warn(c->mtid < 0, "unrecognized mate reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    }
    // mpos
    c->mpos = hts_str2uint(p, &p, 62, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->mpos < 0 && c->mtid >= 0) {
        _parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
        c->mtid = -1;
    }
    // tlen
    c->isize = hts_str2int(p, &p, 63, &overflow);
    if (*p++ != '\t') goto err_ret;
    _parse_err(overflow, "number outside allowed range");
    // seq
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(p - q - 1 > INT32_MAX, "read sequence is too long");
        c->l_qseq = p - q - 1;
        hts_pos_t ql = bam_cigar2qlen(c->n_cigar, (uint32_t*)(b->data + c->l_qname));
        _parse_err(c->n_cigar && ql != c->l_qseq, "CIGAR and query sequence are of different length");
        i = (c->l_qseq + 1) >> 1;
        _get_mem(uint8_t, &t, b, i);

        unsigned int lqs2 = c->l_qseq&~1, i;
        for (i = 0; i < lqs2; i+=2)
            t[i>>1] = (seq_nt16_table[(unsigned char)q[i]] << 4) | seq_nt16_table[(unsigned char)q[i+1]];
        for (; i < c->l_qseq; ++i)
            t[i>>1] = seq_nt16_table[(unsigned char)q[i]] << ((~i&1)<<2);
    } else c->l_qseq = 0;
    // qual
    _get_mem(uint8_t, &t, b, c->l_qseq);
    if (p[0] == '*' && (p[1] == '\t' || p[1] == '\0')) {
        memset(t, 0xff, c->l_qseq);
        p += 2;
    } else {
        int failed = 0;
        _parse_err(s->l - (p - s->s) < c->l_qseq
                   || (p[c->l_qseq] != '\t' && p[c->l_qseq] != '\0'),
                   "SEQ and QUAL are of different length");
        COPY_MINUS_N(t, p, 33, c->l_qseq, failed);
        _parse_err(failed, "invalid QUAL character");
        p += c->l_qseq + 1;
    }

    // aux
    if (aux_parse(p, s->s + s->l, b, 0, NULL) < 0)
        goto err_ret;

    if (bam_tag2cigar(b, 1, 1) < 0)
        return -2;
    return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token
err_ret:
    return -2;
}

static uint32_t read_ncigar(const char *q) {
    uint32_t n_cigar = 0;
    for (; *q && *q != '\t'; ++q)
        if (!isdigit_c(*q)) ++n_cigar;
    if (!n_cigar) {
        hts_log_error("No CIGAR operations");
        return 0;
    }
    if (n_cigar >= 2147483647) {
        hts_log_error("Too many CIGAR operations");
        return 0;
    }

    return n_cigar;
}

/*! @function
 @abstract  Parse a CIGAR string into preallocated a uint32_t array
 @param  in      [in]  pointer to the source string
 @param  a_cigar [out]  address of the destination uint32_t buffer
 @return         number of processed input characters; 0 on error
 */
static int parse_cigar(const char *in, uint32_t *a_cigar, uint32_t n_cigar) {
    int i, overflow = 0;
    const char *p = in;
    for (i = 0; i < n_cigar; i++) {
        uint32_t len;
        int op;
        char *q;
        len = hts_str2uint(p, &q, 28, &overflow)<<BAM_CIGAR_SHIFT;
        if (q == p) {
            hts_log_error("CIGAR length invalid at position %d (%s)", (int)(i+1), p);
            return 0;
        }
        if (overflow) {
            hts_log_error("CIGAR length too long at position %d (%.*s)", (int)(i+1), (int)(q-p+1), p);
            return 0;
        }
        p = q;
        op = bam_cigar_table[(unsigned char)*p++];
        if (op < 0) {
            hts_log_error("Unrecognized CIGAR operator");
            return 0;
        }
        a_cigar[i] = len;
        a_cigar[i] |= op;
    }

    return p-in;
}

ssize_t sam_parse_cigar(const char *in, char **end, uint32_t **a_cigar, size_t *a_mem) {
    size_t n_cigar = 0;
    int diff;

    if (!in || !a_cigar || !a_mem) {
        hts_log_error("NULL pointer arguments");
        return -1;
    }
    if (end) *end = (char *)in;

    if (*in == '*') {
        if (end) (*end)++;
        return 0;
    }
    n_cigar = read_ncigar(in);
    if (!n_cigar) return 0;
    if (n_cigar > *a_mem) {
        uint32_t *a_tmp = realloc(*a_cigar, n_cigar*sizeof(**a_cigar));
        if (a_tmp) {
            *a_cigar = a_tmp;
            *a_mem = n_cigar;
        } else {
            hts_log_error("Memory allocation error");
            return -1;
        }
    }

    if (!(diff = parse_cigar(in, *a_cigar, n_cigar))) return -1;
    if (end) *end = (char *)in+diff;

    return n_cigar;
}

ssize_t bam_parse_cigar(const char *in, char **end, bam1_t *b) {
    size_t n_cigar = 0;
    int diff;

    if (!in || !b) {
        hts_log_error("NULL pointer arguments");
        return -1;
    }
    if (end) *end = (char *)in;

    n_cigar = (*in == '*') ? 0 : read_ncigar(in);
    if (!n_cigar && b->core.n_cigar == 0) {
        if (end) *end = (char *)in+1;
        return 0;
    }

    ssize_t cig_diff = n_cigar - b->core.n_cigar;
    if (cig_diff > 0 &&
        possibly_expand_bam_data(b, cig_diff * sizeof(uint32_t)) < 0) {
        hts_log_error("Memory allocation error");
        return -1;
    }

    uint32_t *cig = bam_get_cigar(b);
    if ((uint8_t *)cig != b->data + b->l_data) {
        // Modifying an BAM existing BAM record
        uint8_t  *seq = bam_get_seq(b);
        memmove(cig + n_cigar, seq, (b->data + b->l_data) - seq);
    }

    if (n_cigar) {
        if (!(diff = parse_cigar(in, cig, n_cigar)))
            return -1;
    } else {
        diff = 1; // handle "*"
    }

    b->l_data += cig_diff * sizeof(uint32_t);
    b->core.n_cigar = n_cigar;
    if (end) *end = (char *)in + diff;

    return n_cigar;
}

/*
 * -----------------------------------------------------------------------------
 * SAM threading
 */
// Size of SAM text block (reading)
#define SAM_NBYTES 240000

// Number of BAM records (writing, up to NB_mem in size)
#define SAM_NBAM 1000

struct SAM_state;

// Output job - a block of BAM records
typedef struct sp_bams {
    struct sp_bams *next;
    int serial;

    bam1_t *bams;
    int nbams, abams; // used and alloc for bams[] array
    size_t bam_mem;   // very approximate total size

    struct SAM_state *fd;
} sp_bams;

// Output job - a block of batch record views copied from raw BAM input.
typedef struct sp_bam_views {
    int serial;

    bam_batch_record_t *records;
    int n_records, a_records;
    uint8_t *data;
    size_t data_len, data_alloc;
    bam_batch_t batch;
    int *ranges, n_ranges;
    int owns_batch;

    struct SAM_state *fd;
} sp_bam_views;

// Input job - a block of SAM text
typedef struct sp_lines {
    struct sp_lines *next;
    int serial;

    char *data;
    int data_size;
    int alloc;

    struct SAM_state *fd;
    sp_bams *bams;
} sp_lines;

enum sam_cmd {
    SAM_NONE = 0,
    SAM_CLOSE,
    SAM_CLOSE_DONE,
    SAM_AT_EOF,
};

typedef struct SAM_state {
    sam_hdr_t *h;

    hts_tpool *p;
    int own_pool;
    pthread_mutex_t lines_m;
    hts_tpool_process *q;
    pthread_t dispatcher;
    int dispatcher_set;

    sp_lines *lines;
    sp_bams *bams;

    sp_bams *curr_bam;
    sp_bam_views *curr_views;
    int curr_idx;
    int serial;

    // Be warned: moving these mutexes around in this struct can reduce
    // threading performance by up to 70%!
    pthread_mutex_t command_m;
    pthread_cond_t command_c;
    enum sam_cmd command;

    // One of the E* errno codes
    int errcode;

    htsFile *fp;
} SAM_state;

// Returns a SAM_state struct from a generic hFILE.
//
// Returns NULL on failure.
static SAM_state *sam_state_create(htsFile *fp) {
    // Ideally sam_open wouldn't be a #define to hts_open but instead would
    // be a redirect call with an additional 'S' mode.  This in turn would
    // correctly set the designed format to sam instead of a generic
    // text_format.
    if (fp->format.format != sam && fp->format.format != text_format)
        return NULL;

    SAM_state *fd = calloc(1, sizeof(*fd));
    if (!fd)
        return NULL;

    fp->state = fd;
    fd->fp = fp;

    return fd;
}

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);
static int sam_format_batch_record_append(const bam_hdr_t *h,
                                          const bam_batch_record_t *record,
                                          kstring_t *str);
static void *sam_format_worker(void *arg);
static void *sam_format_batch_worker(void *arg);
static void *sam_dispatcher_write(void *vp);

static void sam_state_err(SAM_state *fd, int errcode) {
    pthread_mutex_lock(&fd->command_m);
    if (!fd->errcode)
        fd->errcode = errcode;
    pthread_mutex_unlock(&fd->command_m);
}

static void sam_free_sp_bams(sp_bams *b) {
    if (!b)
        return;

    if (b->bams) {
        int i;
        for (i = 0; i < b->abams; i++) {
            if (b->bams[i].data)
                free(b->bams[i].data);
        }
        free(b->bams);
    }
    free(b);
}

static void sam_free_sp_bam_views(sp_bam_views *b) {
    if (!b)
        return;
    if (b->owns_batch)
        sam_bam_batch_destroy(&b->batch);
    free(b->ranges);
    free(b->records);
    free(b->data);
    free(b);
}

static int sam_start_threaded_output(htsFile *fp, const sam_hdr_t *h) {
    SAM_state *fd = (SAM_state *)fp->state;

    if (!fd->h) {
        // NB: discard const.  We don't actually modify sam_hdr_t here,
        // just data pointed to by it, but our cached pointer must be
        // non-const as sam_hdr_destroy takes non-const.
        fd->h = (sam_hdr_t *)h;
        fd->h->ref_count++;

        if (pthread_create(&fd->dispatcher, NULL, sam_dispatcher_write,
                           fp) != 0)
            return -2;
        fd->dispatcher_set = 1;
    }

    if (fd->h != h) {
        hts_log_error("SAM multi-threaded writing does not support changing header");
        return -2;
    }
    return 0;
}

// Destroys the state produce by sam_state_create.
int sam_state_destroy(htsFile *fp) {
    int ret = 0;

    if (!fp->state)
        return 0;

    SAM_state *fd = fp->state;
    if (fd->p) {
        if (fd->h) {
            // Notify sam_dispatcher we're closing
            pthread_mutex_lock(&fd->command_m);
            if (fd->command != SAM_CLOSE_DONE)
                fd->command = SAM_CLOSE;
            pthread_cond_signal(&fd->command_c);
            ret = -fd->errcode;
            if (fd->q)
                hts_tpool_wake_dispatch(fd->q); // unstick the reader

            if (!fp->is_write && fd->q && fd->dispatcher_set) {
                for (;;) {
                    // Avoid deadlocks with dispatcher
                    if (fd->command == SAM_CLOSE_DONE)
                        break;
                    hts_tpool_wake_dispatch(fd->q);
                    pthread_mutex_unlock(&fd->command_m);
                    hts_usleep(10000);
                    pthread_mutex_lock(&fd->command_m);
                }
            }
            pthread_mutex_unlock(&fd->command_m);

            if (fp->is_write) {
                // Dispatch the last partial block.
                sp_bams *gb = fd->curr_bam;
                if (!ret && gb && gb->nbams > 0 && fd->q)
                    ret = hts_tpool_dispatch(fd->p, fd->q, sam_format_worker, gb);
                fd->curr_bam = NULL;
                sp_bam_views *gv = fd->curr_views;
                if (!ret && gv && gv->n_records > 0 && fd->q)
                    ret = hts_tpool_dispatch(fd->p, fd->q,
                                             sam_format_batch_worker, gv);
                else
                    sam_free_sp_bam_views(gv);
                fd->curr_views = NULL;

                // Flush and drain output
                if (fd->q)
                    hts_tpool_process_flush(fd->q);
                pthread_mutex_lock(&fd->command_m);
                if (!ret) ret = -fd->errcode;
                pthread_mutex_unlock(&fd->command_m);

                while (!ret && fd->q && !hts_tpool_process_empty(fd->q)) {
                    hts_usleep(10000);
                    pthread_mutex_lock(&fd->command_m);
                    ret = -fd->errcode;
                    // not empty but shutdown implies error
                    if (hts_tpool_process_is_shutdown(fd->q) && !ret)
                        ret = EIO;
                    pthread_mutex_unlock(&fd->command_m);
                }
                if (fd->q)
                    hts_tpool_process_shutdown(fd->q);
            }

            // Wait for it to acknowledge
            if (fd->dispatcher_set)
                pthread_join(fd->dispatcher, NULL);
            if (!ret) ret = -fd->errcode;
        }

        // Tidy up memory
        if (fd->q)
            hts_tpool_process_destroy(fd->q);

        if (fd->own_pool && fp->format.compression == no_compression) {
            hts_tpool_destroy(fd->p);
            fd->p = NULL;
        }
        pthread_mutex_destroy(&fd->lines_m);
        pthread_mutex_destroy(&fd->command_m);
        pthread_cond_destroy(&fd->command_c);

        sp_lines *l = fd->lines;
        while (l) {
            sp_lines *n = l->next;
            free(l->data);
            free(l);
            l = n;
        }

        sp_bams *b = fd->bams;
        while (b) {
            if (fd->curr_bam == b)
                fd->curr_bam = NULL;
            sp_bams *n = b->next;
            sam_free_sp_bams(b);
            b = n;
        }

        if (fd->curr_bam)
            sam_free_sp_bams(fd->curr_bam);
        if (fd->curr_views)
            sam_free_sp_bam_views(fd->curr_views);

        // Decrement counter by one, maybe destroying too.
        // This is to permit the caller using bam_hdr_destroy
        // before sam_close without triggering decode errors
        // in the background threads.
        bam_hdr_destroy(fd->h);
    }

    free(fp->state);
    fp->state = NULL;
    return ret;
}

// Cleanup function - job for sam_parse_worker; result for sam_format_worker
static void cleanup_sp_lines(void *arg) {
    sp_lines *gl = (sp_lines *)arg;
    if (!gl) return;

    // Should always be true for lines passed to / from thread workers.
    assert(gl->next == NULL);

    free(gl->data);
    sam_free_sp_bams(gl->bams);
    free(gl);
}

// Run from one of the worker threads.
// Convert a passed in array of lines to array of BAMs, returning
// the result back to the thread queue.
static void *sam_parse_worker(void *arg) {
    sp_lines *gl = (sp_lines *)arg;
    sp_bams *gb = NULL;
    char *lines = gl->data;
    int i;
    bam1_t *b;
    SAM_state *fd = gl->fd;

    // Use a block of BAM structs we had earlier if available.
    pthread_mutex_lock(&fd->lines_m);
    if (fd->bams) {
        gb = fd->bams;
        fd->bams = gb->next;
    }
    pthread_mutex_unlock(&fd->lines_m);

    if (gb == NULL) {
        gb = calloc(1, sizeof(*gb));
        if (!gb) {
            return NULL;
        }
        gb->abams = 100;
        gb->bams = b = calloc(gb->abams, sizeof(*b));
        if (!gb->bams) {
            sam_state_err(fd, ENOMEM);
            goto err;
        }
        gb->nbams = 0;
        gb->bam_mem = 0;
    }
    gb->serial = gl->serial;
    gb->next = NULL;

    b = (bam1_t *)gb->bams;
    if (!b) {
        sam_state_err(fd, ENOMEM);
        goto err;
    }

    i = 0;
    char *cp = lines, *cp_end = lines + gl->data_size;
    while (cp < cp_end) {
        if (i >= gb->abams) {
            int old_abams = gb->abams;
            gb->abams *= 2;
            b = (bam1_t *)realloc(gb->bams, gb->abams*sizeof(bam1_t));
            if (!b) {
                gb->abams /= 2;
                sam_state_err(fd, ENOMEM);
                goto err;
            }
            memset(&b[old_abams], 0, (gb->abams - old_abams)*sizeof(*b));
            gb->bams = b;
        }

        // Ideally we'd get sam_parse1 to return the number of
        // bytes decoded and to be able to stop on newline as
        // well as \0.
        //
        // We can then avoid the additional strchr loop.
        // It's around 6% of our CPU cost, albeit threadable.
        //
        // However this is an API change so for now we copy.

        char *nl = strchr(cp, '\n');
        char *line_end;
        if (nl) {
            line_end = nl;
            if (line_end > cp && *(line_end - 1) == '\r')
                line_end--;
            nl++;
        } else {
            nl = line_end = cp_end;
        }
        *line_end = '\0';
        kstring_t ks = { line_end - cp, gl->alloc, cp };
        if (sam_parse1(&ks, fd->h, &b[i]) < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            cleanup_sp_lines(gl);
            goto err;
        }

        cp = nl;
        i++;
    }
    gb->nbams = i;

    pthread_mutex_lock(&fd->lines_m);
    gl->next = fd->lines;
    fd->lines = gl;
    pthread_mutex_unlock(&fd->lines_m);
    return gb;

 err:
    sam_free_sp_bams(gb);
    return NULL;
}

static void *sam_parse_eof(void *arg) {
    return NULL;
}

// Cleanup function - result for sam_parse_worker; job for sam_format_worker
static void cleanup_sp_bams(void *arg) {
    sam_free_sp_bams((sp_bams *) arg);
}

static void cleanup_sp_bam_views(void *arg) {
    sam_free_sp_bam_views((sp_bam_views *) arg);
}

// Runs in its own thread.
// Reads a block of text (SAM) and sends a new job to the thread queue to
// translate this to BAM.
static void *sam_dispatcher_read(void *vp) {
    htsFile *fp = vp;
    kstring_t line = {0};
    int line_frag = 0;
    SAM_state *fd = fp->state;
    sp_lines *l = NULL;

    // Pre-allocate buffer for left-over bits of line (exact size doesn't
    // matter as it will grow if necessary).
    if (ks_resize(&line, 1000) < 0)
        goto err;

    for (;;) {
        // Check for command
        pthread_mutex_lock(&fd->command_m);
        switch (fd->command) {

        case SAM_CLOSE:
            pthread_cond_signal(&fd->command_c);
            pthread_mutex_unlock(&fd->command_m);
            hts_tpool_process_shutdown(fd->q);
            goto tidyup;

        default:
            break;
        }
        pthread_mutex_unlock(&fd->command_m);

        pthread_mutex_lock(&fd->lines_m);
        if (fd->lines) {
            // reuse existing line buffer
            l = fd->lines;
            fd->lines = l->next;
        }
        pthread_mutex_unlock(&fd->lines_m);

        if (l == NULL) {
            // none to reuse, to create a new one
            l = calloc(1, sizeof(*l));
            if (!l)
                goto err;
            l->alloc = SAM_NBYTES;
            l->data = malloc(l->alloc+8); // +8 for optimisation in sam_parse1
            if (!l->data) {
                free(l);
                l = NULL;
                goto err;
            }
            l->fd = fd;
        }
        l->next = NULL;

        if (l->alloc < line_frag+SAM_NBYTES/2) {
            char *rp = realloc(l->data, line_frag+SAM_NBYTES/2 +8);
            if (!rp)
                goto err;
            l->alloc = line_frag+SAM_NBYTES/2;
            l->data = rp;
        }
        memcpy(l->data, line.s, line_frag);

        l->data_size = line_frag;
        ssize_t nbytes;
    longer_line:
        if (fp->is_bgzf)
            nbytes = bgzf_read(fp->fp.bgzf, l->data + line_frag, l->alloc - line_frag);
        else
            nbytes = hread(fp->fp.hfile, l->data + line_frag, l->alloc - line_frag);
        if (nbytes < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            goto err;
        } else if (nbytes == 0)
            break; // EOF
        l->data_size += nbytes;

        // trim to last \n. Maybe \r\n, but that's still fine
        if (nbytes == l->alloc - line_frag) {
            char *cp_end = l->data + l->data_size;
            char *cp = cp_end-1;

            while (cp > (char *)l->data && *cp != '\n')
                cp--;

            // entire buffer is part of a single line
            if (cp == l->data) {
                line_frag = l->data_size;
                char *rp = realloc(l->data, l->alloc * 2 + 8);
                if (!rp)
                    goto err;
                l->alloc *= 2;
                l->data = rp;
                assert(l->alloc >= l->data_size);
                assert(l->alloc >= line_frag);
                assert(l->alloc >= l->alloc - line_frag);
                goto longer_line;
            }
            cp++;

            // line holds the remainder of our line.
            if (ks_resize(&line, cp_end - cp) < 0)
                goto err;
            memcpy(line.s, cp, cp_end - cp);
            line_frag = cp_end - cp;
            l->data_size = l->alloc - line_frag;
        } else {
            // out of buffer
            line_frag = 0;
        }

        l->serial = fd->serial++;
        //fprintf(stderr, "Dispatching %p, %d bytes, serial %d\n", l, l->data_size, l->serial);
        if (hts_tpool_dispatch3(fd->p, fd->q, sam_parse_worker, l,
                                cleanup_sp_lines, cleanup_sp_bams, 0) < 0)
            goto err;
        pthread_mutex_lock(&fd->command_m);
        if (fd->command == SAM_CLOSE) {
            pthread_mutex_unlock(&fd->command_m);
            l = NULL;
            goto tidyup;
        }
        l = NULL;  // Now "owned" by sam_parse_worker()
        pthread_mutex_unlock(&fd->command_m);
    }

    // Submit a NULL sp_bams entry to act as an EOF marker
    if (hts_tpool_dispatch(fd->p, fd->q, sam_parse_eof, NULL) < 0)
        goto err;

    // At EOF, wait for close request.
    // (In future if we add support for seek, this is where we need to catch it.)
    for (;;) {
        pthread_mutex_lock(&fd->command_m);
        if (fd->command == SAM_NONE)
            pthread_cond_wait(&fd->command_c, &fd->command_m);
        switch (fd->command) {
        case SAM_CLOSE:
            pthread_cond_signal(&fd->command_c);
            pthread_mutex_unlock(&fd->command_m);
            hts_tpool_process_shutdown(fd->q);
            goto tidyup;

        default:
            pthread_mutex_unlock(&fd->command_m);
            break;
        }
    }

 tidyup:
    pthread_mutex_lock(&fd->command_m);
    fd->command = SAM_CLOSE_DONE;
    pthread_cond_signal(&fd->command_c);
    pthread_mutex_unlock(&fd->command_m);

    if (l) {
        pthread_mutex_lock(&fd->lines_m);
        l->next = fd->lines;
        fd->lines = l;
        pthread_mutex_unlock(&fd->lines_m);
    }
    free(line.s);

    return NULL;

 err:
    sam_state_err(fd, errno ? errno : ENOMEM);
    hts_tpool_process_shutdown(fd->q);
    goto tidyup;
}

// Runs in its own thread.
// Takes encoded blocks of SAM off the thread results queue and writes them
// to our output stream.
static void *sam_dispatcher_write(void *vp) {
    htsFile *fp = vp;
    SAM_state *fd = fp->state;
    hts_tpool_result *r;

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(fd->q))) {
        sp_lines *gl = (sp_lines *)hts_tpool_result_data(r);
        if (!gl) {
            sam_state_err(fd, ENOMEM);
            goto err;
        }

        if (fp->idx) {
            sp_bams *gb = gl->bams;
            int i = 0, count = 0;
            while (i < gl->data_size) {
                int j = i;
                while (i < gl->data_size && gl->data[i] != '\n')
                    i++;
                if (i < gl->data_size)
                    i++;

                if (fp->is_bgzf) {
                    if (bgzf_flush_try(fp->fp.bgzf, i-j) < 0)
                        goto err;
                    if (bgzf_write(fp->fp.bgzf, &gl->data[j], i-j) != i-j)
                        goto err;
                } else {
                    if (hwrite(fp->fp.hfile, &gl->data[j], i-j) != i-j)
                        goto err;
                }

                bam1_t *b = &gb->bams[count++];
                if (fp->format.compression == bgzf) {
                    if (bgzf_idx_push(fp->fp.bgzf, fp->idx,
                                      b->core.tid, b->core.pos, bam_endpos(b),
                                      bgzf_tell(fp->fp.bgzf),
                                      !(b->core.flag&BAM_FUNMAP)) < 0) {
                        sam_state_err(fd, errno ? errno : ENOMEM);
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(fd->h, b->core.tid), sam_hdr_tid2len(fd->h, b->core.tid), b->core.flag, b->core.pos+1);
                        goto err;
                    }
                } else {
                    if (hts_idx_push(fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                     bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        sam_state_err(fd, errno ? errno : ENOMEM);
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(fd->h, b->core.tid), sam_hdr_tid2len(fd->h, b->core.tid), b->core.flag, b->core.pos+1);
                        goto err;
                    }
                }
            }

            assert(count == gb->nbams);

            // Add bam array to free-list
            pthread_mutex_lock(&fd->lines_m);
            gb->next = fd->bams;
            fd->bams = gl->bams;
            gl->bams = NULL;
            pthread_mutex_unlock(&fd->lines_m);
        } else {
            if (fp->is_bgzf) {
                // We keep track of how much in the current block we have
                // remaining => R.  We look for the last newline in input
                // [i] to [i+R], backwards => position N.
                //
                // If we find a newline, we write out bytes i to N.
                // We know we cannot fit the next record in this bgzf block,
                // so we flush what we have and copy input N to i+R into
                // the start of a new block, and recompute a new R for that.
                //
                // If we don't find a newline (i==N) then we cannot extend
                // the current block at all, so flush whatever is in it now
                // if it ends on a newline.
                // We still copy i(==N) to i+R to the next block and
                // continue as before with a new R.
                //
                // The only exception on the flush is when we run out of
                // data in the input.  In that case we skip it as we don't
                // yet know if the next record will fit.
                //
                // Both conditions share the same code here:
                // - Look for newline (pos N)
                // - Write i to N (which maybe 0)
                // - Flush if block ends on newline and not end of input
                // - write N to i+R

                int i = 0;
                BGZF *fb = fp->fp.bgzf;
                while (i < gl->data_size) {
                    // remaining space in block
                    int R = BGZF_BLOCK_SIZE - fb->block_offset;
                    int eod = 0;
                    if (R > gl->data_size-i)
                        R = gl->data_size-i, eod = 1;

                    // Find last newline in input data
                    int N = i + R;
                    while (--N > i) {
                        if (gl->data[N] == '\n')
                            break;
                    }

                    if (N != i) {
                        // Found a newline
                        N++;
                        if (bgzf_write(fb, &gl->data[i], N-i) != N-i)
                            goto err;
                    }

                    // Flush bgzf block
                    int b_off = fb->block_offset;
                    if (!eod && b_off &&
                        ((char *)fb->uncompressed_block)[b_off-1] == '\n')
                        if (bgzf_flush_try(fb, BGZF_BLOCK_SIZE) < 0)
                            goto err;

                    // Copy from N onwards into next block
                    if (i+R > N)
                        if (bgzf_write(fb, &gl->data[N], i+R - N)
                            != i+R - N)
                            goto err;

                    i = i+R;
                }
            } else {
                if (hwrite(fp->fp.hfile, gl->data, gl->data_size) != gl->data_size)
                    goto err;
            }
        }

        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&fd->lines_m);
        gl->next = fd->lines;
        fd->lines = gl;
        pthread_mutex_unlock(&fd->lines_m);
    }

    sam_state_err(fd, 0); // success
    hts_tpool_process_shutdown(fd->q);
    return NULL;

 err:
    sam_state_err(fd, errno ? errno : EIO);
    return (void *)-1;
}

// Run from one of the worker threads.
// Convert a passed in array of BAMs (sp_bams) and converts to a block
// of text SAM records (sp_lines).
static void *sam_format_worker(void *arg) {
    sp_bams *gb = (sp_bams *)arg;
    sp_lines *gl = NULL;
    int i;
    SAM_state *fd = gb->fd;
    htsFile *fp = fd->fp;

    // Use a block of SAM strings we had earlier if available.
    pthread_mutex_lock(&fd->lines_m);
    if (fd->lines) {
        gl = fd->lines;
        fd->lines = gl->next;
    }
    pthread_mutex_unlock(&fd->lines_m);

    if (gl == NULL) {
        gl = calloc(1, sizeof(*gl));
        if (!gl) {
            sam_state_err(fd, ENOMEM);
            return NULL;
        }
        gl->alloc = gl->data_size = 0;
        gl->data = NULL;
    }
    gl->serial = gb->serial;
    gl->next = NULL;

    kstring_t ks = {0, gl->alloc, gl->data};

    for (i = 0; i < gb->nbams; i++) {
        if (sam_format1_append(fd->h, &gb->bams[i], &ks) < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            goto err;
        }
        kputc('\n', &ks);
    }

    pthread_mutex_lock(&fd->lines_m);
    gl->data_size = ks.l;
    gl->alloc = ks.m;
    gl->data = ks.s;

    if (fp->idx) {
        // Keep hold of the bam array a little longer as
        // sam_dispatcher_write needs to use them for building the index.
        gl->bams = gb;
    } else {
        // Add bam array to free-list
        gb->next = fd->bams;
        fd->bams = gb;
    }
    pthread_mutex_unlock(&fd->lines_m);

    return gl;

 err:
    // Possible race between this and fd->curr_bam.
    // Easier to not free and leave it on the input list so it
    // gets freed there instead?
    // sam_free_sp_bams(gb);
    if (gl) {
        free(gl->data);
        free(gl);
    }
    return NULL;
}

// Run from one of the worker threads.
// Convert copied raw BAM record views to a block of text SAM records.
static void *sam_format_batch_worker(void *arg) {
    sp_bam_views *gv = (sp_bam_views *)arg;
    sp_lines *gl = NULL;
    int i;
    SAM_state *fd = gv->fd;

    pthread_mutex_lock(&fd->lines_m);
    if (fd->lines) {
        gl = fd->lines;
        fd->lines = gl->next;
    }
    pthread_mutex_unlock(&fd->lines_m);

    if (gl == NULL) {
        gl = calloc(1, sizeof(*gl));
        if (!gl) {
            sam_state_err(fd, ENOMEM);
            sam_free_sp_bam_views(gv);
            return NULL;
        }
        gl->alloc = gl->data_size = 0;
        gl->data = NULL;
    }
    gl->serial = gv->serial;
    gl->next = NULL;
    gl->bams = NULL;

    kstring_t ks = {0, gl->alloc, gl->data};

    if (gv->owns_batch) {
        int r;

        for (r = 0; r < gv->n_ranges; r++) {
            int beg = gv->ranges[r * 2];
            int end_i = gv->ranges[r * 2 + 1];

            for (i = beg; i < end_i; i++) {
                if (sam_format_batch_record_append(fd->h,
                                                   &gv->batch.records[i],
                                                   &ks) < 0) {
                    sam_state_err(fd, errno ? errno : EIO);
                    goto err;
                }
                kputc('\n', &ks);
            }
        }
    } else {
        for (i = 0; i < gv->n_records; i++) {
            if (sam_format_batch_record_append(fd->h, &gv->records[i], &ks) < 0) {
                sam_state_err(fd, errno ? errno : EIO);
                goto err;
            }
            kputc('\n', &ks);
        }
    }

    gl->data_size = ks.l;
    gl->alloc = ks.m;
    gl->data = ks.s;
    sam_free_sp_bam_views(gv);
    return gl;

 err:
    if (gl) {
        free(ks.s);
        free(gl);
    }
    sam_free_sp_bam_views(gv);
    return NULL;
}

#define BAM_DEFERRED_THREADS_MAGIC 0x62746872u
#define BAM_STREAM_READER_MAGIC 0x62737472u
#define BAM_BATCH_REQUEST_MAGIC 0x62627172u
#define BAM_STREAM_READER_DEFAULT_CHUNK 32768
#define BAM_STREAM_PARSE_BATCH_RECORDS 16384
#define BAM_STREAM_PARSE_BATCH_BYTES (16u << 20)
#define BAM_STREAM_DEFAULT_QSIZE 64
#define BAM_BATCH_VIEW_INIT_RECORDS 512
// Fused jobs must drain in order before dispatching the next physical block
// range: a BAM record can span BGZF blocks, and the next job cannot know its
// carry-in until the previous job has been parsed.
#define BAM_STREAM_FUSED_MAX_IN_FLIGHT 1
#define BAM_STREAM_FUSED_BLOCKS_PER_JOB 64

typedef struct bam_batch_request_t {
    uint32_t magic;
} bam_batch_request_t;

typedef struct bam_deferred_threads_t {
    uint32_t magic;
    int n_threads;
    int qsize;
    hts_tpool *pool;
    int use_pool;
} bam_deferred_threads_t;

typedef struct bam_stream_reader_t {
    uint32_t magic;
    BGZF *bgzf;
    bgzf_block_data_t *serial_block;
    bgzf_block_data_t *block;
    hts_tpool *pool;
    hts_tpool_process *decode_q;
    hts_tpool_process *parse_q;
    hts_tpool_process *fused_q;
    hts_tpool_result *block_result;
    hts_tpool_result *parse_result;
    struct bam_stream_parse_job_t *parse_batch;
    int own_pool;
    int qsize;
    int in_flight;
    int parse_in_flight;
    int fused_in_flight;
    int input_eof;
    int input_error;
    int input_paused_after_empty;
    int parse_input_eof;
    int parse_input_error;
    int parse_hit_limit;
    int pending_frame_error;
    int parse_views;
    size_t block_off;
    size_t seek_block_off;
    size_t parse_i;
    uint8_t *carry;
    size_t carry_len;
    size_t carry_cap;
    uint64_t carry_voff_beg;
    uint8_t *buf;
    size_t off;
    size_t len;
    size_t cap;
    size_t chunk_size;
    int eof;
} bam_stream_reader_t;

typedef struct bam_stream_decode_job_t {
    bgzf_block_data_t block;
} bam_stream_decode_job_t;

typedef struct bam_stream_parse_job_t {
    uint8_t *data;
    const uint8_t *ref_data;
    size_t len;
    size_t cap;
    int n_records;
    int n_parsed;
    int ret;
    int is_be;
    int hit_limit;
    int need_voff;
    uint64_t voff_beg;
    uint64_t voff_end;
    uint64_t voff_next;
    bgzf_block_data_t *fused_block;
    bgzf_block_data_t *fused_blocks;
    int n_fused_blocks;
    bam_batch_segment_t *segments;
    int n_segments;
    int m_segments;
    bam1_t *records;
    bam_batch_record_t *views;
    sam_hdr_t *view_header;
    hts_tpool_result *owned_block_result;
    uint8_t *carry_out;
    size_t carry_out_len;
    size_t carry_out_cap;
    uint64_t carry_voff_beg;
    size_t block_off_end;
    int errcode;
} bam_stream_parse_job_t;

typedef struct bam_stream_fused_job_t {
    bgzf_block_data_t block;
    bgzf_block_data_t *extra_blocks;
    int n_extra_blocks;
    uint8_t *carry;
    size_t carry_len;
    size_t carry_cap;
    uint64_t carry_voff_beg;
    size_t start_off;
    uint64_t limit_voff;
    sam_hdr_t *h;
    int is_be;
    int need_voff;
    int already_decoded;
} bam_stream_fused_job_t;

static int bam_ordered_env_enabled(void);
static int bam_batch_env_enabled(void);
static int bam_stream_env_enabled(void);
static int bam_stream_parse_env_enabled(void);
static int bam_batch_fused_env_enabled(void);

static int bam_stream_env_strict(void)
{
    const char *env = getenv("HTS_BAM_STREAM_READER_REQUIRE");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_stream_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_STREAM_READER");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_batch_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_BATCH_READER");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_batch_env_strict(void)
{
    const char *env = getenv("HTS_BAM_BATCH_READER_REQUIRE");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_stream_parse_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_STREAM_PARSE");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_batch_parse_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_BATCH_PARSE");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_batch_fused_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_BATCH_FUSED");
    return env && *env && strcmp(env, "0") != 0;
}

static size_t bam_stream_env_chunk_size(void)
{
    const char *env = getenv("HTS_BAM_STREAM_READER_CHUNK");
    char *end = NULL;
    long n;

    if (!env || !*env)
        return BAM_STREAM_READER_DEFAULT_CHUNK;

    errno = 0;
    n = strtol(env, &end, 10);
    if (errno || end == env || *end || n < 1 ||
        n > BGZF_MAX_BLOCK_SIZE)
        return BAM_STREAM_READER_DEFAULT_CHUNK;
    return (size_t)n;
}

static int bam_deferred_threads_set(htsFile *fp, int n_threads,
                                    htsThreadPool *p)
{
    bam_deferred_threads_t *cfg;

    if (fp->state)
        return -2;
    if (n_threads <= 0 && (!p || !p->pool))
        return -1;

    cfg = calloc(1, sizeof(*cfg));
    if (!cfg)
        return -1;

    cfg->magic = BAM_DEFERRED_THREADS_MAGIC;
    if (p && p->pool) {
        cfg->pool = p->pool;
        cfg->qsize = p->qsize;
        cfg->n_threads = hts_tpool_size(p->pool);
        cfg->use_pool = 1;
    } else {
        cfg->n_threads = n_threads;
    }

    if (cfg->n_threads <= 0) {
        free(cfg);
        return -1;
    }

    fp->state = cfg;
    return 0;
}

static void bam_deferred_threads_destroy(bam_deferred_threads_t *cfg)
{
    free(cfg);
}

static int bam_batch_state_is_request(const htsFile *fp)
{
    return fp && fp->state &&
           *(uint32_t *)fp->state == BAM_BATCH_REQUEST_MAGIC;
}

static void bam_batch_request_destroy(htsFile *fp)
{
    if (bam_batch_state_is_request(fp)) {
        free(fp->state);
        fp->state = NULL;
    }
}

static int bam_deferred_threads_enable_bgzf(htsFile *fp,
                                            bam_deferred_threads_t *cfg)
{
    if (!fp || !cfg || !fp->is_bgzf || !fp->fp.bgzf)
        return -1;

    if (cfg->use_pool)
        return bgzf_thread_pool(fp->fp.bgzf, cfg->pool, cfg->qsize);
    return bgzf_mt(fp->fp.bgzf, cfg->n_threads, 256);
}

static void bam_stream_decode_job_free(void *arg)
{
    free(arg);
}

static void bam_stream_parse_job_free(void *arg)
{
    bam_stream_parse_job_t *job = (bam_stream_parse_job_t *)arg;
    int i;

    if (!job)
        return;
    if (job->records) {
        for (i = 0; i < job->n_parsed; i++)
            free(job->records[i].data);
        free(job->records);
    }
    free(job->fused_block);
    free(job->fused_blocks);
    free(job->segments);
    free(job->views);
    if (job->owned_block_result)
        hts_tpool_delete_result(job->owned_block_result, 1);
    free(job->carry_out);
    free(job->data);
    free(job);
}

static void bam_stream_fused_job_free(void *arg)
{
    bam_stream_fused_job_t *job = (bam_stream_fused_job_t *)arg;

    if (!job)
        return;
    free(job->extra_blocks);
    free(job->carry);
    free(job);
}

static void *bam_stream_decode_worker(void *arg)
{
    bam_stream_decode_job_t *job = (bam_stream_decode_job_t *)arg;

    if (bgzf_decode_block_data(NULL, &job->block) < 0 && !job->block.errcode)
        job->block.errcode = BGZF_ERR_ZLIB;
    return job;
}

static int bam_stream_parse_job_build_views(bam_stream_parse_job_t *job,
                                            sam_hdr_t *h);

static void *bam_stream_parse_worker(void *arg)
{
    bam_stream_parse_job_t *job = (bam_stream_parse_job_t *)arg;
    BGZF fake_bgzf = {0};
    const uint8_t *data = job->data ? job->data : job->ref_data;
    size_t off = 0;

    fake_bgzf.is_be = job->is_be;
    job->records = calloc((size_t)job->n_records, sizeof(*job->records));
    if (!job->records) {
        job->ret = -4;
        return job;
    }

    while (job->n_parsed < job->n_records) {
        int32_t block_len;
        int ret;

        if (job->len - off < 4) {
            job->ret = -2;
            return job;
        }
        block_len = le_to_i32(data + off);
        if (block_len < 32 || job->len - off < 4 + (size_t)block_len) {
            job->ret = -4;
            return job;
        }
        ret = bam_decode1_body(job->is_be ? &fake_bgzf : NULL,
                               &job->records[job->n_parsed], block_len,
                               data + off + 4);
        if (ret < 0) {
            job->ret = ret;
            return job;
        }
        off += 4 + (size_t)block_len;
        job->n_parsed++;
    }

    job->ret = 0;
    return job;
}

static void *bam_stream_parse_view_worker(void *arg)
{
    bam_stream_parse_job_t *job = (bam_stream_parse_job_t *)arg;

    job->ret = bam_stream_parse_job_build_views(job, job->view_header);
    return job;
}

static uint64_t bam_stream_block_voff(const bgzf_block_data_t *block,
                                      size_t offset)
{
    if (offset >= 0x10000)
        return ((uint64_t)(block->block_address + block->comp_len)) << 16;
    return ((uint64_t)block->block_address << 16) | offset;
}

static int bam_stream_reader_init_threads(bam_stream_reader_t *reader,
                                          bam_deferred_threads_t *cfg,
                                          int parse_views)
{
    int qsize;
    int enable_fused_views = parse_views && bam_batch_fused_env_enabled();
    int enable_parse_views = parse_views && bam_batch_parse_env_enabled();
    int enable_parse_records = !parse_views && bam_stream_parse_env_enabled();

    if (!cfg || cfg->n_threads <= 1)
        return 0;

    if (cfg->use_pool) {
        reader->pool = cfg->pool;
        reader->own_pool = 0;
    } else {
        reader->pool = hts_tpool_init(cfg->n_threads);
        if (!reader->pool)
            return -1;
        reader->own_pool = 1;
    }

    qsize = cfg->qsize;
    if (qsize <= 0)
        qsize = BAM_STREAM_DEFAULT_QSIZE;
    if (qsize <= 0)
        qsize = 2;
    reader->qsize = qsize;
    reader->decode_q = hts_tpool_process_init(reader->pool, qsize, 0);
    if (!reader->decode_q) {
        if (reader->own_pool) {
            hts_tpool_destroy(reader->pool);
            reader->pool = NULL;
            reader->own_pool = 0;
        }
        return -1;
    }
    if (enable_fused_views) {
        int fqsize = qsize < BAM_STREAM_FUSED_MAX_IN_FLIGHT
                     ? qsize : BAM_STREAM_FUSED_MAX_IN_FLIGHT;

        if (fqsize <= 0)
            fqsize = 1;
        reader->fused_q = hts_tpool_process_init(reader->pool, fqsize, 0);
        if (!reader->fused_q) {
            hts_tpool_process_destroy(reader->decode_q);
            reader->decode_q = NULL;
            if (reader->own_pool) {
                hts_tpool_destroy(reader->pool);
                reader->pool = NULL;
                reader->own_pool = 0;
            }
            return -1;
        }
        return 0;
    }
    if (enable_parse_views || enable_parse_records) {
        reader->parse_q = hts_tpool_process_init(reader->pool, qsize, 0);
        if (!reader->parse_q) {
            hts_tpool_process_destroy(reader->decode_q);
            reader->decode_q = NULL;
            if (reader->own_pool) {
                hts_tpool_destroy(reader->pool);
                reader->pool = NULL;
                reader->own_pool = 0;
            }
            return -1;
        }
        reader->parse_views = enable_parse_views;
    }

    return 0;
}

static bam_stream_reader_t *bam_stream_reader_open(htsFile *fp,
                                                   bam_deferred_threads_t *cfg,
                                                   int parse_views)
{
    bam_stream_reader_t *reader;

    if (!fp || !fp->is_bgzf || !fp->fp.bgzf || fp->fp.bgzf->is_gzip ||
        fp->fp.bgzf->mt)
        return NULL;

    reader = calloc(1, sizeof(*reader));
    if (!reader)
        return NULL;

    reader->magic = BAM_STREAM_READER_MAGIC;
    reader->bgzf = fp->fp.bgzf;
    reader->chunk_size = bam_stream_env_chunk_size();
    reader->serial_block = calloc(1, sizeof(*reader->serial_block));
    if (!reader->serial_block) {
        free(reader);
        return NULL;
    }
    if (bam_stream_reader_init_threads(reader, cfg, parse_views) < 0) {
        free(reader->serial_block);
        free(reader);
        return NULL;
    }
    if (fp->fp.bgzf->block_length > fp->fp.bgzf->block_offset) {
        reader->block = reader->serial_block;
        reader->block->block_address = fp->fp.bgzf->block_address;
        reader->block->comp_len = fp->fp.bgzf->block_clength;
        reader->block->uncomp_len = fp->fp.bgzf->block_length;
        reader->block_off = (size_t)fp->fp.bgzf->block_offset;
        memcpy(reader->block->uncomp_data, fp->fp.bgzf->uncompressed_block,
               (size_t)fp->fp.bgzf->block_length);
    } else {
        reader->seek_block_off = (size_t)fp->fp.bgzf->block_offset;
    }
    return reader;
}

static void bam_stream_reader_destroy(bam_stream_reader_t *reader)
{
    if (!reader)
        return;
    if (reader->block_result)
        hts_tpool_delete_result(reader->block_result, 1);
    if (reader->parse_result)
        hts_tpool_delete_result(reader->parse_result, 1);
    if (reader->parse_q)
        hts_tpool_process_destroy(reader->parse_q);
    if (reader->fused_q)
        hts_tpool_process_destroy(reader->fused_q);
    if (reader->decode_q)
        hts_tpool_process_destroy(reader->decode_q);
    if (reader->own_pool && reader->pool)
        hts_tpool_destroy(reader->pool);
    free(reader->serial_block);
    free(reader->carry);
    free(reader->buf);
    free(reader);
}

static int bam_stream_reader_reserve(bam_stream_reader_t *reader,
                                     size_t extra)
{
    uint8_t *new_buf;
    size_t avail = reader->len - reader->off;
    size_t need = avail + extra;
    size_t new_cap;

    if (need <= reader->cap - reader->off)
        return 0;

    if (reader->off > 0 && avail > 0)
        memmove(reader->buf, reader->buf + reader->off, avail);
    reader->off = 0;
    reader->len = avail;
    if (need <= reader->cap)
        return 0;

    new_cap = reader->cap ? reader->cap : reader->chunk_size;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {
            errno = ENOMEM;
            return -1;
        }
        new_cap *= 2;
    }

    new_buf = realloc(reader->buf, new_cap);
    if (!new_buf) {
        errno = ENOMEM;
        return -1;
    }
    reader->buf = new_buf;
    reader->cap = new_cap;
    return 0;
}

static void bam_stream_reader_release_block(bam_stream_reader_t *reader)
{
    if (reader->block_result) {
        hts_tpool_delete_result(reader->block_result, 1);
        reader->block_result = NULL;
    }
    reader->block = NULL;
    reader->block_off = 0;
}

static int bam_stream_reader_dispatch_block(bam_stream_reader_t *reader)
{
    bam_stream_decode_job_t *job;

    if (reader->input_eof || reader->input_error)
        return reader->input_error ? -1 : 0;

    job = malloc(sizeof(*job));
    if (!job)
        return -1;
    job->block.block_address = 0;
    job->block.comp_len = 0;
    job->block.uncomp_len = 0;
    job->block.hit_eof = 0;
    job->block.errcode = 0;

    if (bgzf_read_block_compressed(reader->bgzf, &job->block) < 0) {
        free(job);
        reader->input_error = 1;
        return -1;
    }
    if (job->block.hit_eof) {
        free(job);
        reader->input_eof = 1;
        return 0;
    }

    int maybe_empty = (le_to_u32(job->block.comp_data +
                                 job->block.comp_len - 4) == 0);

    if (hts_tpool_dispatch3(reader->pool, reader->decode_q,
                            bam_stream_decode_worker, job,
                            bam_stream_decode_job_free,
                            bam_stream_decode_job_free, 0) < 0) {
        free(job);
        reader->input_error = 1;
        return -1;
    }
    reader->in_flight++;
    if (maybe_empty)
        reader->input_paused_after_empty = 1;
    return 0;
}

static int bam_stream_reader_fill_decode(bam_stream_reader_t *reader)
{
    while (!reader->input_eof && !reader->input_error &&
           !reader->input_paused_after_empty &&
           reader->in_flight < reader->qsize) {
        if (bam_stream_reader_dispatch_block(reader) < 0)
            return -1;
    }
    return reader->input_error ? -1 : 0;
}

static int bam_stream_reader_next_parallel_block(bam_stream_reader_t *reader)
{
    BGZF *bgzf = reader->bgzf;

    bam_stream_reader_release_block(reader);

    for (;;) {
        hts_tpool_result *result;
        bgzf_block_data_t *block;

        if (bam_stream_reader_fill_decode(reader) < 0 &&
            reader->in_flight == 0)
            return -1;
        if (reader->in_flight == 0)
            return reader->input_eof ? 0 : -1;

        result = hts_tpool_next_result_wait(reader->decode_q);
        if (!result) {
            reader->input_error = 1;
            return -1;
        }
        reader->in_flight--;
        block = &((bam_stream_decode_job_t *)hts_tpool_result_data(result))->block;
        reader->block_result = result;
        reader->block = block;
        reader->block_off = 0;

        if (block->errcode) {
            bgzf->errcode |= block->errcode;
            return -1;
        }

        if (bgzf_block_data_update_index(bgzf, block) < 0)
            return -1;
        if (block->uncomp_len == 0) {
            bam_stream_reader_release_block(reader);
            reader->input_paused_after_empty = 0;
            continue;
        }
        return 1;
    }
}

static int bam_stream_reader_next_serial_block(bam_stream_reader_t *reader)
{
    BGZF *bgzf = reader->bgzf;

    bam_stream_reader_release_block(reader);
    reader->block = reader->serial_block;
    if (bgzf_read_block_data(bgzf, reader->block) < 0) {
        reader->block = NULL;
        return -1;
    }
    reader->block_off = 0;
    if (reader->block->hit_eof)
        return 0;

    bgzf->block_address = reader->block->block_address;
    bgzf->block_clength = reader->block->comp_len;
    bgzf->block_length = reader->block->uncomp_len;
    bgzf->block_offset = 0;
    return 1;
}

static int bam_stream_reader_next_block(bam_stream_reader_t *reader)
{
    int ret;

    if (reader->decode_q)
        ret = bam_stream_reader_next_parallel_block(reader);
    else
        ret = bam_stream_reader_next_serial_block(reader);

    if (ret > 0 && reader->block) {
        BGZF *bgzf = reader->bgzf;
        if (reader->seek_block_off) {
            reader->block_off =
                reader->seek_block_off <= (size_t)reader->block->uncomp_len
                ? reader->seek_block_off : (size_t)reader->block->uncomp_len;
            reader->seek_block_off = 0;
        }
        bgzf->block_address = reader->block->block_address;
        bgzf->block_clength = reader->block->comp_len;
        bgzf->block_length = reader->block->uncomp_len;
        bgzf->block_offset = (int)reader->block_off;
    }
    return ret;
}

static ssize_t bam_stream_reader_read_block_bytes(bam_stream_reader_t *reader,
                                                  uint8_t *dst, size_t len)
{
    BGZF *bgzf = reader->bgzf;
    size_t copied = 0;

    while (copied < len) {
        int n;
        size_t avail;

        if (!reader->block ||
            reader->block_off >= (size_t)reader->block->uncomp_len) {
            int ret = bam_stream_reader_next_block(reader);
            if (ret < 0)
                return copied ? (ssize_t)copied : -1;
            if (ret == 0)
                break;
        }

        avail = (size_t)reader->block->uncomp_len - reader->block_off;
        if (avail == 0)
            continue;
        n = (int)avail;
        if ((size_t)n > len - copied)
            n = (int)(len - copied);

        memcpy(dst + copied, reader->block->uncomp_data + reader->block_off,
               (size_t)n);
        reader->block_off += (size_t)n;
        bgzf->block_offset = (int)reader->block_off;
        copied += (size_t)n;
    }

    return (ssize_t)copied;
}

static int bam_stream_reader_need(bam_stream_reader_t *reader, size_t need)
{
    while (reader->len - reader->off < need && !reader->eof) {
        size_t want = need - (reader->len - reader->off);
        ssize_t n;

        if (want > reader->chunk_size)
            want = reader->chunk_size;
        if (want > SIZE_MAX - reader->len) {
            errno = ENOMEM;
            return -2;
        }
        if (bam_stream_reader_reserve(reader, want) < 0)
            return -2;

        n = bam_stream_reader_read_block_bytes(reader, reader->buf + reader->len,
                                               want);
        if (n < 0)
            return -2;
        if (n == 0) {
            reader->eof = 1;
            break;
        }
        reader->len += (size_t)n;
    }

    if (reader->len - reader->off >= need)
        return 0;
    return reader->len == reader->off ? -1 : -2;
}

static int bam_stream_parse_job_reserve(bam_stream_parse_job_t *job,
                                        size_t extra)
{
    uint8_t *new_data;
    size_t need = job->len + extra;
    size_t new_cap;

    if (need <= job->cap)
        return 0;

    new_cap = job->cap ? job->cap : BAM_STREAM_PARSE_BATCH_BYTES;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {
            errno = ENOMEM;
            return -1;
        }
        new_cap *= 2;
    }

    new_data = realloc(job->data, new_cap);
    if (!new_data) {
        errno = ENOMEM;
        return -1;
    }
    job->data = new_data;
    job->cap = new_cap;
    return 0;
}

static int bam_stream_reader_carry_reserve(bam_stream_reader_t *reader,
                                           size_t need)
{
    uint8_t *new_carry;
    size_t new_cap;

    if (need <= reader->carry_cap)
        return 0;

    new_cap = reader->carry_cap ? reader->carry_cap : 256;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {
            errno = ENOMEM;
            return -1;
        }
        new_cap *= 2;
    }

    new_carry = realloc(reader->carry, new_cap);
    if (!new_carry) {
        errno = ENOMEM;
        return -1;
    }
    reader->carry = new_carry;
    reader->carry_cap = new_cap;
    return 0;
}

static int bam_stream_reader_carry_append(bam_stream_reader_t *reader,
                                          const uint8_t *src, size_t len)
{
    if (bam_stream_reader_carry_reserve(reader, reader->carry_len + len) < 0)
        return -2;
    memcpy(reader->carry + reader->carry_len, src, len);
    reader->carry_len += len;
    return 0;
}

static int bam_stream_reader_carry_reserve_frame(bam_stream_reader_t *reader,
                                                 const uint8_t *data,
                                                 size_t pos, size_t end)
{
    int32_t block_len;
    size_t frame_len;

    if (end - pos < 4)
        return 0;
    block_len = le_to_i32(data + pos);
    if (block_len < 32)
        return 0;
    frame_len = 4 + (size_t)block_len;
    if (frame_len > BAM_STREAM_PARSE_BATCH_BYTES)
        return 0;
    if (frame_len <= end - pos)
        return 0;
    if (reader->carry_len == 0 && reader->block)
        reader->carry_voff_beg = bam_stream_block_voff(reader->block, pos);
    return bam_stream_reader_carry_reserve(reader, frame_len) < 0 ? -2 : 0;
}

static int bam_stream_reader_carry_error(bam_stream_reader_t *reader)
{
    if (reader->carry_len == 0)
        return -1;
    if (reader->carry_len < 4)
        return -2;
    if (reader->carry_len < 4 + 32)
        return -3;
    return -4;
}

static int bam_stream_reader_next_frame(bam_stream_reader_t *reader,
                                        bam_stream_parse_job_t *job)
{
    int32_t block_len;
    size_t frame_len, avail;
    int ret;

    ret = bam_stream_reader_need(reader, 4);
    if (ret < 0) {
        if (ret == -2)
            reader->off = reader->len;
        return ret;
    }

    block_len = le_to_i32(reader->buf + reader->off);
    if (block_len < 32) {
        reader->off += 4;
        return -4;
    }
    frame_len = 4 + (size_t)block_len;

    ret = bam_stream_reader_need(reader, 4 + 32);
    if (ret < 0) {
        reader->off = reader->len;
        return -3;
    }

    if (bam_validate1_body_core(block_len, reader->buf + reader->off + 4) < 0) {
        reader->off += 4 + 32;
        return -4;
    }

    ret = bam_stream_reader_need(reader, frame_len);
    if (ret < 0) {
        avail = reader->len - reader->off;
        reader->off = reader->len;
        return avail < 4 + 32 ? -3 : -4;
    }

    if (job) {
        if (bam_stream_parse_job_reserve(job, frame_len) < 0)
            return -2;
        memcpy(job->data + job->len, reader->buf + reader->off, frame_len);
        job->len += frame_len;
        job->n_records++;
        reader->off += frame_len;
        return 1;
    }

    return 1;
}

static int bam_stream_reader_next_serial(bam_stream_reader_t *reader, bam1_t *b)
{
    int32_t block_len;
    size_t frame_len;
    int ret = bam_stream_reader_next_frame(reader, NULL);

    if (ret < 0)
        return ret;

    block_len = le_to_i32(reader->buf + reader->off);
    frame_len = 4 + (size_t)block_len;
    ret = bam_decode1_body(reader->bgzf, b, block_len,
                           reader->buf + reader->off + 4);
    reader->off += frame_len;
    return ret;
}

static void bam_batch_record_view_set(bam_batch_record_t *view,
                                      const uint8_t *frame, size_t frame_len,
                                      uint64_t voff_beg, uint64_t voff_end)
{
    const uint8_t *body = frame + 4;
    uint32_t x2 = le_to_u32(body + 8);
    uint32_t x3 = le_to_u32(body + 12);

    view->frame = frame;
    view->frame_len = frame_len;
    view->body = body + 32;
    view->raw_l_data = (uint32_t)(frame_len - 4 - 32);
    view->voff_beg = voff_beg;
    view->voff_end = voff_end;
    view->endpos = 0;
    view->flags = BAM_BATCH_RECORD_F_RAW_LAYOUT_VALID;
    view->core.tid = le_to_i32(body);
    view->core.pos = le_to_i32(body + 4);
    view->core.bin = x2 >> 16;
    view->core.qual = (x2 >> 8) & 0xff;
    view->core.l_qname = x2 & 0xff;
    view->core.l_extranul = (-view->core.l_qname) & 3;
    view->core.flag = x3 >> 16;
    view->core.n_cigar = x3 & 0xffff;
    view->core.l_qseq = le_to_i32(body + 16);
    view->core.mtid = le_to_i32(body + 20);
    view->core.mpos = le_to_i32(body + 24);
    view->core.isize = le_to_i32(body + 28);
}

static int bam_body_header_valid(const uint8_t *body, sam_hdr_t *h)
{
    int32_t tid, mtid;

    if (!h)
        return 1;
    tid = le_to_i32(body);
    mtid = le_to_i32(body + 20);
    return tid >= -1 && tid < h->n_targets &&
           mtid >= -1 && mtid < h->n_targets;
}

static int bam_batch_record_view_tid_valid(const bam_batch_record_t *view,
                                           sam_hdr_t *h)
{
    if (!h)
        return 1;
    return view->core.tid >= -1 && view->core.tid < h->n_targets &&
           view->core.mtid >= -1 && view->core.mtid < h->n_targets;
}

static int sam_bam_batch_record_decode_status(
        bam_batch_record_t *record, int *needs_materialize);
static int sam_bam_batch_record_materialize_validated(
        const bam_batch_record_t *record, bam1_t *scratch, int *materialized);

int sam_bam_batch_record_cg_candidate(const bam_batch_record_t *record)
{
    const uint8_t *cigar;
    uint32_t first;

    if (!record || record->core.n_cigar == 0 || record->core.tid < 0 ||
        record->core.pos < 0)
        return 0;
    cigar = sam_bam_batch_record_cigar(record);
    first = le_to_u32(cigar);
    return first == (((uint32_t)record->core.l_qseq << BAM_CIGAR_SHIFT) |
                     BAM_CSOFT_CLIP);
}

static int bam_batch_record_view_append(bam_batch_record_t **views,
                                        int *n_views, int *m_views,
                                        const uint8_t *frame,
                                        size_t frame_len,
                                        uint64_t voff_beg,
                                        uint64_t voff_end,
                                        sam_hdr_t *h)
{
    bam_batch_record_t *new_views;
    bam_batch_record_t view;
    int new_m;

    bam_batch_record_view_set(&view, frame, frame_len, voff_beg, voff_end);
    if (!bam_batch_record_view_tid_valid(&view, h)) {
        errno = ERANGE;
        return -3;
    }
    {
        int needs_materialize = 0;
        if (sam_bam_batch_record_decode_status(&view, &needs_materialize) < 0)
            return -4;
    }

    if (*n_views == *m_views) {
        new_m = *m_views ? (*m_views > INT_MAX / 2
                            ? INT_MAX : *m_views * 2)
                         : BAM_BATCH_VIEW_INIT_RECORDS;
        if (new_m == *m_views ||
            (size_t)new_m > SIZE_MAX / sizeof(*new_views)) {
            errno = ENOMEM;
            return -1;
        }
        new_views = realloc(*views, (size_t)new_m * sizeof(*new_views));
        if (!new_views) {
            errno = ENOMEM;
            return -1;
        }
        *views = new_views;
        *m_views = new_m;
    }

    (*views)[(*n_views)++] = view;
    return 0;
}

static int bam_stream_parse_job_reserve_view_data(bam_stream_parse_job_t *job,
                                                  size_t extra)
{
    const uint8_t *old_data = job->data;
    int i;

    if (bam_stream_parse_job_reserve(job, extra) < 0)
        return -2;
    if (old_data && old_data != job->data) {
        for (i = 0; i < job->n_records; i++) {
            if (!(job->views[i].flags & BAM_BATCH_RECORD_F_OWNED))
                continue;
            size_t frame_off = (size_t)(job->views[i].frame - old_data);
            size_t body_off = (size_t)(job->views[i].body - old_data);

            job->views[i].frame = job->data + frame_off;
            job->views[i].body = job->data + body_off;
        }
        for (i = 0; i < job->n_segments; i++) {
            if (job->segments[i].flags & BAM_BATCH_SEGMENT_F_OWNED) {
                size_t data_off = (size_t)(job->segments[i].data - old_data);

                job->segments[i].data = job->data + data_off;
            }
        }
    }
    return 0;
}

static int bam_stream_parse_job_append_segment(bam_stream_parse_job_t *job,
                                               const uint8_t *data,
                                               size_t len, uint32_t flags)
{
    bam_batch_segment_t *segments;
    int new_m;

    if (len == 0)
        return 0;
    if (!data) {
        errno = EINVAL;
        return -2;
    }
    if (job->n_segments > 0) {
        bam_batch_segment_t *last = &job->segments[job->n_segments - 1];

        if (last->flags == flags && last->data + last->len == data) {
            last->len += len;
            return 0;
        }
    }
    if (job->n_segments == job->m_segments) {
        new_m = job->m_segments ? (job->m_segments > INT_MAX / 2
                                   ? INT_MAX : job->m_segments * 2)
                                : 16;
        if (new_m == job->m_segments ||
            (size_t)new_m > SIZE_MAX / sizeof(*segments)) {
            errno = ENOMEM;
            return -2;
        }
        segments = realloc(job->segments, (size_t)new_m * sizeof(*segments));
        if (!segments) {
            errno = ENOMEM;
            return -2;
        }
        job->segments = segments;
        job->m_segments = new_m;
    }
    job->segments[job->n_segments].data = data;
    job->segments[job->n_segments].len = len;
    job->segments[job->n_segments].flags = flags;
    job->n_segments++;
    return 0;
}

static int bam_stream_fused_carry_reserve(uint8_t **carry, size_t *cap,
                                          size_t need)
{
    uint8_t *new_carry;
    size_t new_cap;

    if (need <= *cap)
        return 0;
    new_cap = *cap ? *cap : 256;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {
            errno = ENOMEM;
            return -2;
        }
        new_cap *= 2;
    }
    new_carry = realloc(*carry, new_cap);
    if (!new_carry) {
        errno = ENOMEM;
        return -2;
    }
    *carry = new_carry;
    *cap = new_cap;
    return 0;
}

static int bam_stream_fused_copy_carry_out(bam_stream_parse_job_t *out,
                                           const uint8_t *src, size_t len,
                                           uint64_t voff_beg)
{
    if (len == 0)
        return 0;
    out->carry_out = malloc(len);
    if (!out->carry_out) {
        errno = ENOMEM;
        return -2;
    }
    memcpy(out->carry_out, src, len);
    out->carry_out_len = len;
    out->carry_out_cap = len;
    out->carry_voff_beg = voff_beg;
    return 0;
}

static int bam_stream_fused_append_view(bam_stream_parse_job_t *out,
                                        int *m_views, const uint8_t *frame,
                                        size_t frame_len, uint64_t voff_beg,
                                        uint64_t voff_end, sam_hdr_t *h)
{
    int vret = bam_batch_record_view_append(&out->views, &out->n_records,
                                            m_views, frame, frame_len,
                                            voff_beg, voff_end, h);
    if (vret == -3 || vret == -4)
        return vret;
    return vret < 0 ? -2 : 0;
}

static int bam_stream_fused_append_copied_view(bam_stream_parse_job_t *out,
                                               int *m_views,
                                               const uint8_t *frame,
                                               size_t frame_len,
                                               uint64_t voff_beg,
                                               uint64_t voff_end,
                                               sam_hdr_t *h)
{
    size_t old_len = out->len;
    int ret;

    if (bam_stream_parse_job_reserve_view_data(out, frame_len) < 0)
        return -2;
    memcpy(out->data + old_len, frame, frame_len);
    ret = bam_stream_fused_append_view(out, m_views, out->data + old_len,
                                       frame_len, voff_beg, voff_end, h);
    if (ret < 0)
        return ret;
    out->views[out->n_records - 1].flags |= BAM_BATCH_RECORD_F_OWNED;
    out->len += frame_len;
    return bam_stream_parse_job_append_segment(out, out->data + old_len,
                                               frame_len,
                                               BAM_BATCH_SEGMENT_F_OWNED);
}

static bgzf_block_data_t *bam_stream_fused_job_block(
        bam_stream_fused_job_t *job, int i)
{
    return i == 0 ? &job->block : &job->extra_blocks[i - 1];
}

static void *bam_stream_fused_worker_multi(bam_stream_fused_job_t *job)
{
    bam_stream_parse_job_t *out;
    uint8_t *carry = NULL;
    size_t carry_len = 0, carry_cap = 0;
    uint64_t carry_voff_beg = job->carry_voff_beg;
    int n_blocks = 1 + job->n_extra_blocks;
    int m_views = 0;
    int bi;

    out = calloc(1, sizeof(*out));
    if (!out) {
        bam_stream_fused_job_free(job);
        return NULL;
    }
    out->is_be = job->is_be;
    out->need_voff = job->need_voff;
    out->fused_blocks = calloc((size_t)n_blocks, sizeof(*out->fused_blocks));
    if (!out->fused_blocks) {
        out->ret = -2;
        bam_stream_fused_job_free(job);
        return out;
    }

    carry = job->carry;
    carry_len = job->carry_len;
    carry_cap = job->carry_cap;
    job->carry = NULL;
    job->carry_len = job->carry_cap = 0;

    for (bi = 0; bi < n_blocks && out->ret == 0 && !out->hit_limit; bi++) {
        bgzf_block_data_t *job_block = bam_stream_fused_job_block(job, bi);
        bgzf_block_data_t *block;
        uint8_t *data;
        size_t pos, end;

        if (!(bi == 0 && job->already_decoded) &&
            bgzf_decode_block_data(NULL, job_block) < 0) {
            out->ret = -2;
            out->errcode = job_block->errcode ? job_block->errcode
                                               : BGZF_ERR_ZLIB;
            break;
        }

        out->fused_blocks[out->n_fused_blocks] = *job_block;
        block = &out->fused_blocks[out->n_fused_blocks++];
        end = block->uncomp_len > 0 ? (size_t)block->uncomp_len : 0;
        pos = bi == 0 && job->start_off <= end ? job->start_off : 0;
        out->block_off_end = pos;
        if (end == 0)
            continue;
        data = block->uncomp_data;

        if (carry_len) {
            int32_t block_len;
            size_t frame_len, need, avail;
            uint64_t voff_end;
            int ret;

            if (carry_voff_beg >= job->limit_voff) {
                out->hit_limit = 1;
                break;
            }
            while (carry_len < 4 && pos < end) {
                if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                                   carry_len + 1) < 0) {
                    out->ret = -2;
                    break;
                }
                carry[carry_len++] = data[pos++];
            }
            if (out->ret < 0)
                break;
            if (carry_len < 4) {
                out->block_off_end = pos;
                continue;
            }

            block_len = le_to_i32(carry);
            if (block_len < 32) {
                out->ret = -4;
                break;
            }
            frame_len = 4 + (size_t)block_len;
            if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                               frame_len) < 0) {
                out->ret = -2;
                break;
            }
            need = frame_len - carry_len;
            avail = end - pos;
            if (need > avail)
                need = avail;
            memcpy(carry + carry_len, data + pos, need);
            carry_len += need;
            pos += need;
            if (carry_len < frame_len) {
                out->block_off_end = pos;
                continue;
            }

            if (bam_validate1_body_core(block_len, carry + 4) < 0) {
                out->ret = -4;
                break;
            }
            if (!bam_body_header_valid(carry + 4, job->h)) {
                out->ret = -3;
                break;
            }
            voff_end = bam_stream_block_voff(block, pos);
            ret = bam_stream_fused_append_copied_view(
                    out, &m_views, carry, frame_len, carry_voff_beg,
                    voff_end, job->h);
            if (ret < 0) {
                out->ret = ret;
                break;
            }
            carry_len = 0;
        }

        while (pos + 4 <= end) {
            int32_t block_len = le_to_i32(data + pos);
            size_t frame_len;
            uint64_t voff_beg =
                job->need_voff || job->limit_voff != UINT64_MAX
                ? bam_stream_block_voff(block, pos) : 0;
            uint64_t voff_end;
            int ret;

            if (voff_beg >= job->limit_voff) {
                out->hit_limit = 1;
                break;
            }
            if (block_len < 32) {
                out->ret = -4;
                break;
            }
            frame_len = 4 + (size_t)block_len;
            if (frame_len > end - pos) {
                if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                                   end - pos) < 0) {
                    out->ret = -2;
                    break;
                }
                memcpy(carry, data + pos, end - pos);
                carry_len = end - pos;
                carry_voff_beg = bam_stream_block_voff(block, pos);
                pos = end;
                break;
            }
            if (bam_validate1_body_core(block_len, data + pos + 4) < 0) {
                out->ret = -4;
                break;
            }
            if (!bam_body_header_valid(data + pos + 4, job->h)) {
                out->ret = -3;
                break;
            }

            voff_end = job->need_voff
                       ? bam_stream_block_voff(block, pos + frame_len) : 0;
            ret = bam_stream_fused_append_view(
                    out, &m_views, data + pos, frame_len, voff_beg,
                    voff_end, job->h);
            if (ret == 0)
                ret = bam_stream_parse_job_append_segment(
                        out, data + pos, frame_len, 0);
            if (ret < 0) {
                out->ret = ret;
                break;
            }
            pos += frame_len;
        }

        if (out->ret == 0 && !out->hit_limit && pos < end &&
            pos + 4 > end) {
            if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                               end - pos) < 0) {
                out->ret = -2;
                break;
            }
            memcpy(carry, data + pos, end - pos);
            carry_len = end - pos;
            carry_voff_beg = bam_stream_block_voff(block, pos);
            pos = end;
        }
        out->block_off_end = pos;
    }

    if (out->ret == 0 && carry_len) {
        out->carry_out = carry;
        out->carry_out_len = carry_len;
        out->carry_out_cap = carry_cap;
        out->carry_voff_beg = carry_voff_beg;
        carry = NULL;
    }

    free(carry);
    bam_stream_fused_job_free(job);
    return out;
}

static void *bam_stream_fused_worker(void *arg)
{
    bam_stream_fused_job_t *job = (bam_stream_fused_job_t *)arg;
    bam_stream_parse_job_t *out;
    uint8_t *carry = NULL;
    size_t carry_len = 0, carry_cap = 0;
    uint8_t *data;
    size_t start, pos, end, batch_end;
    int m_views = 0;
    int copy_mode = 0;

    if (job->n_extra_blocks > 0)
        return bam_stream_fused_worker_multi(job);

    out = calloc(1, sizeof(*out));
    if (!out) {
        bam_stream_fused_job_free(job);
        return NULL;
    }
    out->is_be = job->is_be;
    out->need_voff = job->need_voff;

    if (!job->already_decoded &&
        bgzf_decode_block_data(NULL, &job->block) < 0) {
        out->ret = -2;
        out->errcode = job->block.errcode ? job->block.errcode
                                          : BGZF_ERR_ZLIB;
        bam_stream_fused_job_free(job);
        return out;
    }

    out->fused_block = malloc(sizeof(*out->fused_block));
    if (!out->fused_block) {
        out->ret = -2;
        bam_stream_fused_job_free(job);
        return out;
    }
    *out->fused_block = job->block;

    end = job->block.uncomp_len > 0 ? (size_t)job->block.uncomp_len : 0;
    start = job->start_off <= end ? job->start_off : end;
    pos = start;
    batch_end = start;
    out->block_off_end = pos;
    if (end == 0) {
        out->ret = 0;
        bam_stream_fused_job_free(job);
        return out;
    }

    carry = job->carry;
    carry_len = job->carry_len;
    carry_cap = job->carry_cap;
    job->carry = NULL;
    job->carry_len = job->carry_cap = 0;

    data = out->fused_block->uncomp_data;
    if (carry_len) {
        int32_t block_len;
        size_t frame_len, need, avail;
        uint64_t voff_end;
        int ret;

        while (carry_len < 4 && pos < end) {
            if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                               carry_len + 1) < 0) {
                out->ret = -2;
                goto done;
            }
            carry[carry_len++] = data[pos++];
        }
        if (carry_len < 4) {
            out->carry_out = carry;
            out->carry_out_len = carry_len;
            out->carry_out_cap = carry_cap;
            out->carry_voff_beg = job->carry_voff_beg;
            carry = NULL;
            out->block_off_end = pos;
            out->ret = 0;
            goto done;
        }

        block_len = le_to_i32(carry);
        if (block_len < 32) {
            out->ret = -4;
            goto done;
        }
        frame_len = 4 + (size_t)block_len;
        if (bam_stream_fused_carry_reserve(&carry, &carry_cap,
                                           frame_len) < 0) {
            out->ret = -2;
            goto done;
        }
        need = frame_len - carry_len;
        avail = end - pos;
        if (need > avail)
            need = avail;
        memcpy(carry + carry_len, data + pos, need);
        carry_len += need;
        pos += need;
        if (carry_len < frame_len) {
            out->carry_out = carry;
            out->carry_out_len = carry_len;
            out->carry_out_cap = carry_cap;
            out->carry_voff_beg = job->carry_voff_beg;
            carry = NULL;
            out->block_off_end = pos;
            out->ret = 0;
            goto done;
        }

        if (bam_validate1_body_core(block_len, carry + 4) < 0) {
            out->ret = -4;
            goto done;
        }
        if (!bam_body_header_valid(carry + 4, job->h)) {
            out->ret = -3;
            goto done;
        }
        out->data = carry;
        out->len = frame_len;
        out->cap = carry_cap;
        carry = NULL;
        copy_mode = 1;
        if (pos < end &&
            bam_stream_parse_job_reserve_view_data(out, end - pos) < 0) {
            out->ret = -2;
            goto done;
        }
        voff_end = bam_stream_block_voff(out->fused_block, pos);
        ret = bam_stream_fused_append_view(out, &m_views, out->data,
                                           frame_len, job->carry_voff_beg,
                                           voff_end, job->h);
        if (ret < 0) {
            out->ret = ret;
            goto done;
        }
        out->views[out->n_records - 1].flags |= BAM_BATCH_RECORD_F_OWNED;
        if (bam_stream_parse_job_append_segment(
                    out, out->data, frame_len,
                    BAM_BATCH_SEGMENT_F_OWNED) < 0) {
            out->ret = -2;
            goto done;
        }
        batch_end = pos;
    }

    while (pos + 4 <= end) {
        int32_t block_len = le_to_i32(data + pos);
        size_t frame_len;
        uint64_t voff_beg = job->need_voff || job->limit_voff != UINT64_MAX
                             ? bam_stream_block_voff(out->fused_block, pos) : 0;
        uint64_t voff_end;
        int ret;

        if (voff_beg >= job->limit_voff) {
            out->hit_limit = 1;
            break;
        }
        if (block_len < 32) {
            out->ret = -4;
            break;
        }
        frame_len = 4 + (size_t)block_len;
        if (frame_len > end - pos) {
            if (bam_stream_fused_copy_carry_out(
                        out, data + pos, end - pos,
                        bam_stream_block_voff(out->fused_block, pos)) < 0)
                out->ret = -2;
            pos = end;
            break;
        }
        if (bam_validate1_body_core(block_len, data + pos + 4) < 0) {
            out->ret = -4;
            break;
        }
        if (!bam_body_header_valid(data + pos + 4, job->h)) {
            out->ret = -3;
            break;
        }

        voff_end = job->need_voff
                   ? bam_stream_block_voff(out->fused_block, pos + frame_len)
                   : 0;
        if (copy_mode) {
            ret = bam_stream_fused_append_copied_view(
                    out, &m_views, data + pos, frame_len, voff_beg,
                    voff_end, job->h);
        } else {
            ret = bam_stream_fused_append_view(
                    out, &m_views, data + pos, frame_len, voff_beg,
                    voff_end, job->h);
        }
        if (ret < 0) {
            out->ret = ret;
            break;
        }
        pos += frame_len;
        batch_end = pos;
    }

    if (out->ret == 0 && pos < end && pos + 4 > end) {
        if (bam_stream_fused_copy_carry_out(
                    out, data + pos, end - pos,
                    bam_stream_block_voff(out->fused_block, pos)) < 0)
            out->ret = -2;
        pos = end;
    }
    if (!copy_mode) {
        out->ref_data = data + start;
        out->len = batch_end - start;
    }
    out->block_off_end = pos;
    out->ret = out->ret ? out->ret : 0;

done:
    free(carry);
    bam_stream_fused_job_free(job);
    return out;
}

static int bam_stream_parse_job_append_owned_frame(bam_stream_parse_job_t *job,
                                                  const uint8_t *frame,
                                                  size_t frame_len)
{
    if (bam_stream_parse_job_reserve(job, frame_len) < 0)
        return -2;
    memcpy(job->data + job->len, frame, frame_len);
    job->len += frame_len;
    job->n_records++;
    return 0;
}

static uint64_t bam_stream_parse_job_voff(const bam_stream_parse_job_t *job,
                                          size_t offset);

static int bam_stream_parse_job_build_views(bam_stream_parse_job_t *job,
                                            sam_hdr_t *h)
{
    bam_batch_record_t *views = NULL;
    const uint8_t *data = job->data ? job->data : job->ref_data;
    size_t pos = 0;
    int expected_records = job->n_records;
    int n_records = 0, m_views = 0;

    if (!data && job->len) {
        job->n_records = 0;
        return -4;
    }
    if (expected_records > 0) {
        views = malloc((size_t)expected_records * sizeof(*views));
        if (!views) {
            job->n_records = 0;
            return -2;
        }
        m_views = expected_records;
    }

    while (pos + 4 <= job->len) {
        int32_t block_len = le_to_i32(data + pos);
        size_t frame_len;
        int vret;

        if (block_len < 32) {
            free(views);
            job->n_records = 0;
            return -4;
        }
        frame_len = 4 + (size_t)block_len;
        if (frame_len > job->len - pos ||
            bam_validate1_body_core(block_len, data + pos + 4) < 0) {
            if (n_records > 0) {
                job->views = views;
                job->n_records = n_records;
                return -4;
            }
            free(views);
            job->n_records = 0;
            return -4;
        }
        vret = bam_batch_record_view_append(&views, &n_records, &m_views,
                                            data + pos, frame_len,
                                            job->need_voff
                                            ? bam_stream_parse_job_voff(job, pos) : 0,
                                            job->need_voff
                                            ? bam_stream_parse_job_voff(job, pos + frame_len) : 0,
                                            h);
        if (vret < 0) {
            if ((vret == -3 || vret == -4) && n_records > 0) {
                job->views = views;
                job->n_records = n_records;
                return vret;
            }
            free(views);
            job->n_records = 0;
            return (vret == -3 || vret == -4) ? vret : -2;
        }
        if (job->data)
            views[n_records - 1].flags |= BAM_BATCH_RECORD_F_OWNED;
        pos += frame_len;
    }

    if (pos != job->len || n_records != expected_records) {
        free(views);
        job->n_records = 0;
        return -4;
    }
    job->views = views;
    job->n_records = n_records;
    return 0;
}

static uint64_t bam_stream_parse_job_voff(const bam_stream_parse_job_t *job,
                                          size_t offset)
{
    size_t block_off;

    if (!job->voff_next) {
        if (offset == 0)
            return job->voff_beg;
        if (offset == job->len)
            return job->voff_end;
        return 0;
    }

    block_off = (size_t)(job->voff_beg & 0xffff) + offset;
    if (block_off >= 0x10000)
        return job->voff_next;
    return (job->voff_beg & ~UINT64_C(0xffff)) | block_off;
}

static int bam_stream_parse_job_detach_ref_data(bam_stream_parse_job_t *job)
{
    uint8_t *data;

    if (job->data || job->owned_block_result)
        return 0;
    if (!job->ref_data || !job->len)
        return -2;

    data = malloc(job->len);
    if (!data) {
        errno = ENOMEM;
        return -2;
    }
    memcpy(data, job->ref_data, job->len);
    job->data = data;
    job->cap = job->len;
    job->ref_data = NULL;
    return 0;
}

static int bam_stream_reader_extend_owned_job_from_block(
        bam_stream_reader_t *reader, bam_stream_parse_job_t *job,
        int build_views, sam_hdr_t *h)
{
    BGZF *bgzf = reader->bgzf;
    uint8_t *data = reader->block->uncomp_data;
    size_t pos = reader->block_off;
    size_t end = (size_t)reader->block->uncomp_len;

    while (pos + 4 <= end &&
           job->n_records < BAM_STREAM_PARSE_BATCH_RECORDS &&
           job->len < BAM_STREAM_PARSE_BATCH_BYTES) {
        int32_t block_len = le_to_i32(data + pos);
        size_t frame_len;

        if (block_len < 32) {
            reader->pending_frame_error = -4;
            break;
        }
        frame_len = 4 + (size_t)block_len;
        if (job->len + frame_len > BAM_STREAM_PARSE_BATCH_BYTES)
            break;
        if (frame_len > end - pos)
            break;
        if (bam_validate1_body_core(block_len, data + pos + 4) < 0) {
            reader->pending_frame_error = -4;
            break;
        }
        if (!build_views && !bam_body_header_valid(data + pos + 4, h)) {
            reader->pending_frame_error = -3;
            errno = ERANGE;
            break;
        }
        if (build_views) {
            bam_batch_record_t view;

            bam_batch_record_view_set(
                    &view, data + pos, frame_len,
                    bam_stream_block_voff(reader->block, pos),
                    bam_stream_block_voff(reader->block, pos + frame_len));
            if (!bam_batch_record_view_tid_valid(&view, h)) {
                reader->pending_frame_error = -3;
                errno = ERANGE;
                break;
            }
        }
        if (bam_stream_parse_job_append_owned_frame(job, data + pos,
                                                   frame_len) < 0)
            return -2;
        pos += frame_len;
    }

    reader->block_off = pos;
    bgzf->block_offset = (int)reader->block_off;

    if (pos < end && !reader->pending_frame_error &&
        job->n_records < BAM_STREAM_PARSE_BATCH_RECORDS &&
        job->len < BAM_STREAM_PARSE_BATCH_BYTES) {
        if (bam_stream_reader_carry_append(reader, data + pos, end - pos) < 0)
            return -2;
        reader->block_off = end;
        bgzf->block_offset = (int)reader->block_off;
        bam_stream_reader_release_block(reader);
    } else if (pos == end || reader->pending_frame_error) {
        bam_stream_reader_release_block(reader);
    }

    return 0;
}

static int bam_stream_reader_build_parse_job(bam_stream_reader_t *reader,
                                             bam_stream_parse_job_t **job_out,
                                             int build_views, sam_hdr_t *h,
                                             uint64_t limit_voff,
                                             int *hit_limit,
                                             int need_voff,
                                             int keep_view_voffs)
{
    BGZF *bgzf = reader->bgzf;
    bam_stream_parse_job_t *job;
    uint8_t *data;
    bam_batch_record_t *views = NULL;
    size_t start, pos, end;
    int n_records = 0;
    int m_views = 0;

    if (hit_limit)
        *hit_limit = 0;

    for (;;) {
        if (!reader->block ||
            reader->block_off >= (size_t)reader->block->uncomp_len) {
            int ret = bam_stream_reader_next_block(reader);
            if (ret < 0)
                return -2;
            if (ret == 0) {
                int err;
                reader->parse_input_eof = 1;
                err = bam_stream_reader_carry_error(reader);
                return err == -1 ? 0 : err;
            }
        }

        data = reader->block->uncomp_data;
        end = (size_t)reader->block->uncomp_len;

        if (reader->carry_len) {
            int32_t block_len;
            size_t frame_len, need, avail;

            if (reader->carry_voff_beg >= limit_voff) {
                if (hit_limit)
                    *hit_limit = 1;
                return 0;
            }
            if (reader->carry_len < 4) {
                need = 4 - reader->carry_len;
                avail = end - reader->block_off;
                if (need > avail)
                    need = avail;
                if (bam_stream_reader_carry_append(reader,
                                                   data + reader->block_off,
                                                   need) < 0)
                    return -2;
                reader->block_off += need;
                bgzf->block_offset = (int)reader->block_off;
                if (reader->carry_len < 4)
                    continue;
            }

            block_len = le_to_i32(reader->carry);
            if (block_len < 32) {
                reader->carry_len = 0;
                return -4;
            }
            frame_len = 4 + (size_t)block_len;
            if (frame_len <= BAM_STREAM_PARSE_BATCH_BYTES &&
                bam_stream_reader_carry_reserve(reader, frame_len) < 0)
                return -2;
            need = frame_len - reader->carry_len;
            avail = end - reader->block_off;
            if (need > avail)
                need = avail;
            if (bam_stream_reader_carry_append(reader, data + reader->block_off,
                                               need) < 0)
                return -2;
            reader->block_off += need;
            bgzf->block_offset = (int)reader->block_off;
            if (reader->carry_len < frame_len)
                continue;
            if (!keep_view_voffs &&
                bam_validate1_body_core(block_len, reader->carry + 4) < 0) {
                reader->carry_len = 0;
                return -4;
            }
            if (!build_views &&
                !bam_body_header_valid(reader->carry + 4, h)) {
                reader->carry_len = 0;
                errno = ERANGE;
                return -3;
            }

            job = calloc(1, sizeof(*job));
            if (!job)
                return -2;
            job->data = reader->carry;
            job->len = frame_len;
            job->cap = reader->carry_cap;
            job->n_records = 1;
            job->is_be = bgzf->is_be;
            job->need_voff = need_voff;
            job->voff_beg = reader->carry_voff_beg;
            job->voff_end = bgzf_tell(bgzf);
            reader->carry = NULL;
            reader->carry_len = reader->carry_cap = 0;
            if (!build_views && !keep_view_voffs && reader->block_off < end) {
                int eret = bam_stream_reader_extend_owned_job_from_block(
                        reader, job, build_views, h);
                if (eret < 0) {
                    bam_stream_parse_job_free(job);
                    return eret;
                }
            }
            if (build_views) {
                int vret = bam_stream_parse_job_build_views(job, h);

                if (vret < 0) {
                    bam_stream_parse_job_free(job);
                    if (vret == -3)
                        errno = ERANGE;
                    return vret;
                }
            }
            if (hit_limit)
                job->hit_limit = *hit_limit;
            *job_out = job;
            return 1;
        }

        start = reader->block_off;
        pos = start;
        while (pos + 4 <= end) {
            int32_t block_len = le_to_i32(data + pos);
            size_t frame_len;
            uint64_t voff_beg = need_voff || limit_voff != UINT64_MAX
                                 ? bam_stream_block_voff(reader->block, pos)
                                 : 0;

            if (voff_beg >= limit_voff) {
                if (hit_limit)
                    *hit_limit = 1;
                break;
            }
            if (block_len < 32) {
                if (pos == start) {
                    reader->block_off += 4;
                    bgzf->block_offset = (int)reader->block_off;
                    return -4;
                }
                reader->pending_frame_error = -4;
                break;
            }
            frame_len = 4 + (size_t)block_len;
            if (pos + frame_len > end)
                break;
            if (!keep_view_voffs &&
                bam_validate1_body_core(block_len, data + pos + 4) < 0) {
                if (pos == start) {
                    reader->block_off += 4 + 32;
                    bgzf->block_offset = (int)reader->block_off;
                    return -4;
                }
                reader->pending_frame_error = -4;
                break;
            }
            if (!build_views &&
                !bam_body_header_valid(data + pos + 4, h)) {
                if (pos == start) {
                    reader->block_off += 4 + 32;
                    bgzf->block_offset = (int)reader->block_off;
                    errno = ERANGE;
                    return -3;
                }
                reader->pending_frame_error = -3;
                errno = ERANGE;
                break;
            }
            if (build_views) {
                int vret = bam_batch_record_view_append(
                        &views, &n_records, &m_views, data + pos, frame_len,
                        voff_beg,
                        need_voff ? bam_stream_block_voff(reader->block,
                                                          pos + frame_len) : 0,
                        h);
                if (vret < 0) {
                    if ((vret == -3 || vret == -4) && n_records > 0) {
                        reader->pending_frame_error = vret;
                        break;
                    }
                    free(views);
                    return (vret == -3 || vret == -4) ? vret : -2;
                }
            }
            pos += frame_len;
            if (!build_views)
                n_records++;
        }

        if (n_records == 0 && hit_limit && *hit_limit) {
            free(views);
            return 0;
        }

        if (n_records > 0) {
            job = calloc(1, sizeof(*job));
            if (!job) {
                free(views);
                return -2;
            }
            job->ref_data = data + start;
            job->len = pos - start;
            job->n_records = n_records;
            job->is_be = bgzf->is_be;
            job->hit_limit = hit_limit ? *hit_limit : 0;
            job->need_voff = need_voff;
            if (need_voff) {
                job->voff_beg = bam_stream_block_voff(reader->block, start);
                job->voff_end = bam_stream_block_voff(reader->block, pos);
                job->voff_next =
                    ((uint64_t)(reader->block->block_address +
                                reader->block->comp_len)) << 16;
            }
            job->views = views;
            views = NULL;
            reader->block_off = pos;
            bgzf->block_offset = (int)reader->block_off;

            if (pos < end && !reader->pending_frame_error &&
                !(hit_limit && *hit_limit)) {
                if (bam_stream_reader_carry_reserve_frame(reader, data, pos,
                                                          end) < 0) {
                    bam_stream_parse_job_free(job);
                    return -2;
                }
                if (bam_stream_reader_carry_append(reader, data + pos,
                                                   end - pos) < 0) {
                    bam_stream_parse_job_free(job);
                    return -2;
                }
                reader->block_off = end;
                bgzf->block_offset = (int)reader->block_off;
            }

            if (!(hit_limit && *hit_limit)) {
                job->owned_block_result = reader->block_result;
                reader->block_result = NULL;
                reader->block = NULL;
            }
            *job_out = job;
            return 1;
        }

        if (pos < end) {
            if (bam_stream_reader_carry_reserve_frame(reader, data, pos,
                                                      end) < 0)
                return -2;
            if (bam_stream_reader_carry_append(reader, data + pos,
                                               end - pos) < 0)
                return -2;
            reader->block_off = end;
            bgzf->block_offset = (int)reader->block_off;
            bam_stream_reader_release_block(reader);
            continue;
        }

        bam_stream_reader_release_block(reader);
    }
}

static int bam_stream_reader_dispatch_parse(bam_stream_reader_t *reader,
                                            bam_stream_parse_job_t *job,
                                            int build_views)
{
    if (bam_stream_parse_job_detach_ref_data(job) < 0) {
        bam_stream_parse_job_free(job);
        return -1;
    }
    if (hts_tpool_dispatch3(reader->pool, reader->parse_q,
                            build_views ? bam_stream_parse_view_worker
                                        : bam_stream_parse_worker,
                            job,
                            bam_stream_parse_job_free,
                            bam_stream_parse_job_free, 0) < 0) {
        bam_stream_parse_job_free(job);
        return -1;
    }
    reader->parse_in_flight++;
    return 0;
}

static int bam_stream_reader_fill_parse(bam_stream_reader_t *reader)
{
    while (!reader->parse_input_eof && !reader->parse_input_error &&
           !reader->pending_frame_error &&
           reader->parse_in_flight < reader->qsize) {
        bam_stream_parse_job_t *job = NULL;
        int ret = bam_stream_reader_build_parse_job(reader, &job, 0, NULL,
                                                    UINT64_MAX, NULL, 0, 0);

        if (ret > 0) {
            if (bam_stream_reader_dispatch_parse(reader, job, 0) < 0) {
                reader->parse_input_error = -2;
                return -1;
            }
        } else if (ret == 0) {
            break;
        } else {
            reader->parse_input_error = ret;
            return -1;
        }
    }

    return reader->parse_input_error ? -1 : 0;
}

static void bam_stream_reader_release_parse_batch(bam_stream_reader_t *reader)
{
    if (reader->parse_result) {
        hts_tpool_delete_result(reader->parse_result, 1);
        reader->parse_result = NULL;
    }
    reader->parse_batch = NULL;
    reader->parse_i = 0;
}

static int bam_stream_reader_next_parsed(bam_stream_reader_t *reader, bam1_t *b)
{
    for (;;) {
        if (reader->parse_batch) {
            bam_stream_parse_job_t *job = reader->parse_batch;

            if (reader->parse_i < (size_t)job->n_parsed) {
                bam1_t *src = &job->records[reader->parse_i++];
                if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) != 0) {
                    if (bam_copy1(b, src) < 0)
                        return -4;
                } else {
                    free(b->data);
                    *b = *src;
                    src->data = NULL;
                    src->l_data = src->m_data = 0;
                }
                return 36 + b->l_data - b->core.l_extranul;
            }

            if (job->ret < 0) {
                int ret = job->ret;
                bam_stream_reader_release_parse_batch(reader);
                return ret;
            }
            bam_stream_reader_release_parse_batch(reader);
        }

        if (bam_stream_reader_fill_parse(reader) < 0 &&
            reader->parse_in_flight == 0)
            return reader->parse_input_error;

        if (reader->parse_in_flight == 0) {
            if (reader->pending_frame_error) {
                int ret = reader->pending_frame_error;
                reader->pending_frame_error = 0;
                if (ret == -3)
                    errno = ERANGE;
                return ret;
            }
            if (reader->parse_input_error)
                return reader->parse_input_error;
            if (reader->parse_input_eof)
                return -1;
        }

        if (reader->parse_in_flight > 0) {
            hts_tpool_result *result =
                hts_tpool_next_result_wait(reader->parse_q);
            if (!result) {
                reader->parse_input_error = -2;
                return -2;
            }
            reader->parse_in_flight--;
            reader->parse_result = result;
            reader->parse_batch =
                (bam_stream_parse_job_t *)hts_tpool_result_data(result);
            reader->parse_i = 0;
        }
    }
}

static int bam_stream_reader_next(bam_stream_reader_t *reader, bam1_t *b)
{
    if (reader->parse_q)
        return bam_stream_reader_next_parsed(reader, b);
    return bam_stream_reader_next_serial(reader, b);
}

static void bam_stream_reader_reset_parse_queue(bam_stream_reader_t *reader)
{
    if (!reader->parse_q)
        return;
    if (reader->parse_result) {
        hts_tpool_delete_result(reader->parse_result, 1);
        reader->parse_result = NULL;
    }
    reader->parse_batch = NULL;
    reader->parse_i = 0;
    hts_tpool_process_reset(reader->parse_q, 1);
    reader->parse_in_flight = 0;
    reader->parse_hit_limit = 0;
}

static int bam_stream_reader_fill_view_parse(bam_stream_reader_t *reader,
                                             sam_hdr_t *h,
                                             uint64_t limit_voff,
                                             int need_voff)
{
    while (!reader->parse_input_eof && !reader->parse_input_error &&
           !reader->pending_frame_error && !reader->parse_hit_limit &&
           reader->parse_in_flight < reader->qsize) {
        bam_stream_parse_job_t *job = NULL;
        int hit_limit = 0;
        int ret = bam_stream_reader_build_parse_job(reader, &job, 0, h,
                                                    limit_voff, &hit_limit,
                                                    need_voff, 1);

        if (ret > 0) {
            job->view_header = h;
            job->hit_limit = hit_limit;
            job->need_voff = need_voff;
            if (bam_stream_reader_dispatch_parse(reader, job, 1) < 0) {
                reader->parse_input_error = -2;
                reader->parse_input_eof = 1;
                return -1;
            }
            if (hit_limit) {
                reader->parse_hit_limit = 1;
                break;
            }
        } else if (ret == 0) {
            if (hit_limit)
                reader->parse_hit_limit = 1;
            break;
        } else {
            reader->parse_input_error = ret;
            reader->parse_input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
            return -1;
        }
    }

    return reader->parse_input_error ? -1 : 0;
}

static int bam_stream_reader_next_view_job(bam_stream_reader_t *reader,
                                           sam_hdr_t *h,
                                           bam_stream_parse_job_t **job_out,
                                           uint64_t limit_voff,
                                           int *hit_limit,
                                           int need_voff)
{
    *job_out = NULL;
    if (hit_limit)
        *hit_limit = 0;

    for (;;) {
        if (!reader->parse_input_eof && !reader->parse_input_error &&
            !reader->pending_frame_error && !reader->parse_hit_limit &&
            bam_stream_reader_fill_view_parse(reader, h, limit_voff,
                                              need_voff) < 0 &&
            reader->parse_in_flight == 0)
            return reader->parse_input_error;

        if (reader->parse_in_flight > 0) {
            hts_tpool_result *result =
                hts_tpool_next_result_wait(reader->parse_q);
            bam_stream_parse_job_t *job;

            if (!result) {
                reader->parse_input_error = -2;
                reader->parse_input_eof = 1;
                return -2;
            }
            reader->parse_in_flight--;
            job = (bam_stream_parse_job_t *)hts_tpool_result_data(result);
            hts_tpool_delete_result(result, 0);

            if (job->ret < 0) {
                int ret = job->ret;

                if (job->views && job->n_records > 0) {
                    reader->pending_frame_error = ret;
                    reader->parse_input_eof = 1;
                    bam_stream_reader_reset_parse_queue(reader);
                    *job_out = job;
                    return job->n_records;
                }

                bam_stream_parse_job_free(job);
                reader->parse_input_eof = 1;
                bam_stream_reader_reset_parse_queue(reader);
                if (ret == -3)
                    errno = ERANGE;
                return ret;
            }

            if (!job->views) {
                bam_stream_parse_job_free(job);
                reader->parse_input_error = -2;
                reader->parse_input_eof = 1;
                bam_stream_reader_reset_parse_queue(reader);
                return -2;
            }

            if (hit_limit)
                *hit_limit = job->hit_limit;
            if (job->hit_limit)
                reader->parse_hit_limit = 0;
            *job_out = job;
            return job->n_records;
        }

        if (reader->pending_frame_error) {
            int ret = reader->pending_frame_error;

            reader->pending_frame_error = 0;
            reader->parse_input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
            return ret;
        }
        if (reader->parse_input_error)
            return reader->parse_input_error;
        if (reader->parse_hit_limit) {
            reader->parse_hit_limit = 0;
            if (hit_limit)
                *hit_limit = 1;
            return 0;
        }
        if (reader->parse_input_eof)
            return 0;
    }
}

static int bam_stream_reader_dispatch_fused(bam_stream_reader_t *reader,
                                            sam_hdr_t *h,
                                            uint64_t limit_voff,
                                            int need_voff)
{
    bam_stream_fused_job_t *job;
    int stop_after_first = 0;

    if (reader->input_eof || reader->input_error)
        return reader->input_error ? -1 : 0;

    if (reader->block &&
        reader->block_off >= (size_t)reader->block->uncomp_len)
        bam_stream_reader_release_block(reader);

    job = calloc(1, sizeof(*job));
    if (!job)
        return -1;
    job->h = h;
    job->is_be = reader->bgzf->is_be;
    job->limit_voff = limit_voff;
    job->need_voff = need_voff;

    if (reader->block) {
        job->block = *reader->block;
        job->start_off = reader->block_off;
        job->already_decoded = 1;
        if (job->block.uncomp_len == 0) {
            stop_after_first = 1;
            reader->input_paused_after_empty = 1;
        }
        bam_stream_reader_release_block(reader);
    } else {
        job->start_off = reader->seek_block_off;
        reader->seek_block_off = 0;
        if (bgzf_read_block_compressed(reader->bgzf, &job->block) < 0) {
            free(job);
            reader->input_error = 1;
            return -1;
        }
        if (job->block.hit_eof) {
            free(job);
            reader->input_eof = 1;
            if (reader->carry_len) {
                reader->pending_frame_error =
                    bam_stream_reader_carry_error(reader);
                free(reader->carry);
                reader->carry = NULL;
                reader->carry_len = reader->carry_cap = 0;
            }
            return 0;
        }
        if (le_to_u32(job->block.comp_data + job->block.comp_len - 4) == 0) {
            stop_after_first = 1;
            reader->input_paused_after_empty = 1;
        }
    }

    while (!stop_after_first &&
           1 + job->n_extra_blocks < BAM_STREAM_FUSED_BLOCKS_PER_JOB) {
        bgzf_block_data_t *new_blocks, *block;

        new_blocks = realloc(job->extra_blocks,
                             (size_t)(job->n_extra_blocks + 1) *
                             sizeof(*job->extra_blocks));
        if (!new_blocks)
            break;
        job->extra_blocks = new_blocks;
        block = &job->extra_blocks[job->n_extra_blocks];
        memset(block, 0, sizeof(*block));

        if (bgzf_read_block_compressed(reader->bgzf, block) < 0) {
            reader->input_error = 1;
            break;
        }
        if (block->hit_eof) {
            reader->input_eof = 1;
            break;
        }
        job->n_extra_blocks++;
        if (le_to_u32(block->comp_data + block->comp_len - 4) == 0) {
            reader->input_paused_after_empty = 1;
            break;
        }
    }

    job->carry = reader->carry;
    job->carry_len = reader->carry_len;
    job->carry_cap = reader->carry_cap;
    job->carry_voff_beg = reader->carry_voff_beg;
    reader->carry = NULL;
    reader->carry_len = reader->carry_cap = 0;

    if (hts_tpool_dispatch3(reader->pool, reader->fused_q,
                            bam_stream_fused_worker, job,
                            bam_stream_fused_job_free,
                            bam_stream_parse_job_free, 0) < 0) {
        bam_stream_fused_job_free(job);
        reader->input_error = 1;
        return -1;
    }
    reader->fused_in_flight++;
    return 1;
}

static int bam_stream_reader_fill_fused(bam_stream_reader_t *reader,
                                        sam_hdr_t *h,
                                        uint64_t limit_voff,
                                        int need_voff)
{
    while (!reader->input_eof && !reader->input_error &&
           !reader->input_paused_after_empty &&
           !reader->pending_frame_error &&
           !reader->parse_hit_limit &&
           reader->fused_in_flight < BAM_STREAM_FUSED_MAX_IN_FLIGHT) {
        int ret = bam_stream_reader_dispatch_fused(reader, h, limit_voff,
                                                   need_voff);

        if (ret < 0) {
            reader->parse_input_error = -2;
            return -1;
        }
        if (ret == 0)
            break;
    }
    return reader->parse_input_error ? -1 : 0;
}

static void bam_stream_reader_accept_fused_carry(
        bam_stream_reader_t *reader, bam_stream_parse_job_t *job)
{
    if (!job->carry_out)
        return;
    free(reader->carry);
    reader->carry = job->carry_out;
    reader->carry_len = job->carry_out_len;
    reader->carry_cap = job->carry_out_cap;
    reader->carry_voff_beg = job->carry_voff_beg;
    job->carry_out = NULL;
    job->carry_out_len = job->carry_out_cap = 0;
    if (reader->input_eof && reader->carry_len) {
        reader->pending_frame_error = bam_stream_reader_carry_error(reader);
        free(reader->carry);
        reader->carry = NULL;
        reader->carry_len = reader->carry_cap = 0;
    }
}

static int bam_stream_reader_update_fused_block(bam_stream_reader_t *reader,
                                                bam_stream_parse_job_t *job)
{
    BGZF *bgzf = reader->bgzf;
    bgzf_block_data_t *block = job->fused_blocks && job->n_fused_blocks > 0
                               ? &job->fused_blocks[job->n_fused_blocks - 1]
                               : job->fused_block;
    size_t block_off;
    int i;

    if (job->errcode && job->n_records <= 0) {
        bgzf->errcode |= job->errcode;
        return -2;
    }
    if (!block)
        return -2;
    if (job->fused_blocks) {
        for (i = 0; i < job->n_fused_blocks; i++) {
            if (bgzf_block_data_update_index(bgzf, &job->fused_blocks[i]) < 0)
                return -2;
        }
    } else if (bgzf_block_data_update_index(bgzf, block) < 0) {
        return -2;
    }

    block_off = job->block_off_end <= (size_t)block->uncomp_len
                ? job->block_off_end : (size_t)block->uncomp_len;
    bgzf->block_address = block->block_address;
    bgzf->block_clength = block->comp_len;
    bgzf->block_length = block->uncomp_len;
    bgzf->block_offset = (int)block_off;
    if (job->errcode)
        bgzf->errcode |= job->errcode;
    if (block->uncomp_len == 0)
        reader->input_paused_after_empty = 0;
    return 0;
}

static int bam_stream_reader_next_fused_job(bam_stream_reader_t *reader,
                                            sam_hdr_t *h,
                                            bam_stream_parse_job_t **job_out,
                                            uint64_t limit_voff,
                                            int *hit_limit,
                                            int need_voff)
{
    *job_out = NULL;
    if (hit_limit)
        *hit_limit = 0;

    for (;;) {
        if (reader->fused_in_flight == 0 &&
            !reader->input_eof && !reader->input_error &&
            !reader->pending_frame_error && !reader->parse_hit_limit &&
            bam_stream_reader_fill_fused(reader, h, limit_voff,
                                         need_voff) < 0 &&
            reader->fused_in_flight == 0)
            return reader->parse_input_error;

        if (reader->fused_in_flight > 0) {
            hts_tpool_result *result =
                hts_tpool_next_result_wait(reader->fused_q);
            bam_stream_parse_job_t *job;

            if (!result) {
                reader->parse_input_error = -2;
                reader->input_eof = 1;
                return -2;
            }
            reader->fused_in_flight--;
            job = (bam_stream_parse_job_t *)hts_tpool_result_data(result);
            hts_tpool_delete_result(result, 0);
            if (!job) {
                reader->parse_input_error = -2;
                reader->input_eof = 1;
                return -2;
            }

            if (bam_stream_reader_update_fused_block(reader, job) < 0) {
                int ret = job->errcode ? -2 : -2;
                bam_stream_parse_job_free(job);
                reader->parse_input_error = ret;
                reader->input_eof = 1;
                return ret;
            }
            bam_stream_reader_accept_fused_carry(reader, job);

            if (job->ret < 0) {
                int ret = job->ret;

                if (job->views && job->n_records > 0) {
                    reader->pending_frame_error = ret;
                    reader->input_eof = 1;
                    *job_out = job;
                    return job->n_records;
                }

                bam_stream_parse_job_free(job);
                reader->input_eof = 1;
                if (ret == -3)
                    errno = ERANGE;
                return ret;
            }

            if (hit_limit)
                *hit_limit = job->hit_limit;
            if (job->hit_limit)
                reader->parse_hit_limit = 1;

            if (job->views && job->n_records > 0) {
                if (!job->hit_limit && !reader->pending_frame_error)
                    (void)bam_stream_reader_fill_fused(reader, h,
                                                       limit_voff, need_voff);
                *job_out = job;
                return job->n_records;
            }

            bam_stream_parse_job_free(job);
            if (reader->parse_hit_limit) {
                reader->parse_hit_limit = 0;
                if (hit_limit)
                    *hit_limit = 1;
                return 0;
            }
            continue;
        }

        if (reader->pending_frame_error) {
            int ret = reader->pending_frame_error;

            reader->pending_frame_error = 0;
            reader->input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
            return ret;
        }
        if (reader->parse_input_error)
            return reader->parse_input_error;
        if (reader->input_error)
            return -2;
        if (reader->input_eof)
            return 0;
    }
}

void sam_bam_batch_destroy(bam_batch_t *batch)
{
    if (!batch)
        return;
    if (batch->impl)
        bam_stream_parse_job_free(batch->impl);
    memset(batch, 0, sizeof(*batch));
}

int sam_bam_batch_record_to_bam1(const bam_batch_record_t *record, bam1_t *bam)
{
    BGZF fake_bgzf = {0};
    int32_t block_len;

    if (!record || !bam || !record->frame || record->frame_len < 36 ||
        record->frame_len - 4 > INT32_MAX)
        return -2;

    block_len = (int32_t)(record->frame_len - 4);
    fake_bgzf.is_be = ed_is_big();
    return bam_decode1_body(fake_bgzf.is_be ? &fake_bgzf : NULL, bam,
                            block_len, record->frame + 4);
}

static int sam_bam_batch_record_materialize_validated(
        const bam_batch_record_t *record, bam1_t *scratch, int *materialized)
{
    if (!scratch) {
        errno = EINVAL;
        return -1;
    }
    if (sam_bam_batch_record_to_bam1(record, scratch) < 0) {
        errno = EINVAL;
        return -1;
    }
    if (materialized)
        *materialized = 1;
    return 0;
}

static int sam_bam_batch_record_cg_fixup_needed(
        const bam_batch_record_t *record, int *needed)
{
    sam_bam_record_view_t view;
    const uint8_t *cg;

    *needed = 0;
    if (!sam_bam_batch_record_cg_candidate(record))
        return 0;

    sam_bam_record_view_from_batch(&view, record);
    cg = sam_bam_record_view_aux_get(&view, "CG");
    if (!cg) {
        if (errno == EINVAL)
            return -1;
        return 0;
    }
    if (*cg == 'B' && (cg[1] == 'I' || cg[1] == 'i')) {
        uint32_t cg_len = le_to_u32(cg + 2);

        if (cg_len >= record->core.n_cigar && cg_len < 1U << 29)
            *needed = 1;
    }
    return 0;
}

static int sam_bam_batch_record_decode_status(
        bam_batch_record_t *record, int *needs_materialize)
{
    const uint8_t *cigar;
    hts_pos_t rlen = 0, qlen = 0;
    int cg_fixup = 0, k;

    *needs_materialize = 0;
    if (!record || !record->body || record->core.l_qname < 1) {
        errno = EINVAL;
        return -1;
    }
    if (record->flags & BAM_BATCH_RECORD_F_RAW_WRITE_SAFE)
        return 0;
    if (record->flags & BAM_BATCH_RECORD_F_NEEDS_MATERIALIZE) {
        *needs_materialize = 1;
        return 0;
    }

    if (sam_bam_batch_record_qname(record)[record->core.l_qname - 1] != '\0') {
        *needs_materialize = 1;
        record->flags |= BAM_BATCH_RECORD_F_NEEDS_MATERIALIZE;
        return 0;
    }

    if (sam_bam_batch_record_cg_fixup_needed(record, &cg_fixup) < 0) {
        errno = EINVAL;
        return -1;
    }
    if (cg_fixup) {
        *needs_materialize = 1;
        record->flags |= BAM_BATCH_RECORD_F_NEEDS_MATERIALIZE;
        return 0;
    }

    if (record->core.n_cigar == 0) {
        record->endpos = record->core.pos + 1;
        record->flags |= BAM_BATCH_RECORD_F_RAW_WRITE_SAFE |
                         BAM_BATCH_RECORD_F_ENDPOS_VALID;
        return 0;
    }

    cigar = sam_bam_batch_record_cigar(record);
    for (k = 0; k < record->core.n_cigar; k++) {
        uint32_t c = le_to_u32(cigar + ((size_t)k << 2));
        int op = bam_cigar_op(c);

        if (bam_cigar_type(op) & 2)
            rlen += bam_cigar_oplen(c);
        if (bam_cigar_type(op) & 1)
            qlen += bam_cigar_oplen(c);
    }
    if ((record->core.flag & BAM_FUNMAP) || rlen == 0)
        rlen = 1;
    if (record->core.bin != hts_reg2bin(record->core.pos,
                                        record->core.pos + rlen, 14, 5)) {
        *needs_materialize = 1;
        record->flags |= BAM_BATCH_RECORD_F_NEEDS_MATERIALIZE;
        return 0;
    }
    if (record->core.l_qseq > 0 && !(record->core.flag & BAM_FUNMAP) &&
        qlen != record->core.l_qseq) {
        int qname_len = record->core.l_qname > 0 ? record->core.l_qname - 1 : 0;
        hts_log_error("CIGAR and query sequence lengths differ for %.*s",
                      qname_len, sam_bam_batch_record_qname(record));
        errno = EINVAL;
        return -1;
    }
    record->endpos = record->core.pos + rlen;
    record->flags |= BAM_BATCH_RECORD_F_RAW_WRITE_SAFE |
                     BAM_BATCH_RECORD_F_ENDPOS_VALID;
    return 0;
}

int sam_bam_batch_record_validate_decode(bam_batch_record_t *record,
                                         bam1_t *scratch, int *materialized)
{
    int needs_materialize = 0;

    if (materialized)
        *materialized = 0;
    if (sam_bam_batch_record_decode_status(record, &needs_materialize) < 0)
        return -1;
    if (needs_materialize)
        return sam_bam_batch_record_materialize_validated(record, scratch,
                                                         materialized);
    return 0;
}

static int sam_bam_batch_record_raw_layout_valid(bam_batch_record_t *record)
{
    size_t raw_l_data, need, seq_len;

    if (record && (record->flags & BAM_BATCH_RECORD_F_RAW_LAYOUT_VALID))
        return 0;

    if (!record || !record->frame || !record->body ||
        record->frame_len < 36 || record->frame_len > INT_MAX ||
        record->frame + 36 != record->body ||
        record->raw_l_data != record->frame_len - 36 ||
        le_to_i32(record->frame) != (int32_t)(record->frame_len - 4) ||
        record->core.l_qname < 1 || record->core.l_qseq < 0) {
        errno = EINVAL;
        return -1;
    }

    raw_l_data = record->raw_l_data;
    need = (size_t)record->core.l_qname;
    if (need > raw_l_data ||
        (size_t)record->core.n_cigar > (SIZE_MAX - need) >> 2) {
        errno = EINVAL;
        return -1;
    }
    need += (size_t)record->core.n_cigar << 2;
    seq_len = ((size_t)record->core.l_qseq + 1) >> 1;
    if (seq_len > SIZE_MAX - need) {
        errno = EINVAL;
        return -1;
    }
    need += seq_len;
    if ((size_t)record->core.l_qseq > SIZE_MAX - need) {
        errno = EINVAL;
        return -1;
    }
    need += (size_t)record->core.l_qseq;
    if (need > raw_l_data) {
        errno = EINVAL;
        return -1;
    }
    record->flags |= BAM_BATCH_RECORD_F_RAW_LAYOUT_VALID;
    return 0;
}

static int sam_bam_batch_record_raw_write_safe(
        bam_batch_record_t *record)
{
    int needs_materialize = 0;

    if (sam_bam_batch_record_raw_layout_valid(record) < 0)
        return -1;
    if (sam_bam_batch_record_decode_status(record, &needs_materialize) < 0)
        return -1;
    return needs_materialize ? 1 : 0;
}

int sam_bam_batch_record_raw_write_status(bam_batch_record_t *record)
{
    return sam_bam_batch_record_raw_write_safe(record);
}

static int sam_bam_raw_write_get_bgzf_common(htsFile *fp, BGZF **bfp,
                                             int allow_hts_index)
{
    if (!fp || !fp->is_write || !fp->is_bgzf || !fp->fp.bgzf) {
        errno = EINVAL;
        return -2;
    }
    if (fp->format.format != bam)
        return -2;
    *bfp = fp->fp.bgzf;
    if ((!allow_hts_index && fp->idx) || (*bfp)->idx ||
        (*bfp)->idx_build_otf) {
        errno = EINVAL;
        return -2;
    }
    return 0;
}

static int sam_bam_raw_write_get_bgzf(htsFile *fp, BGZF **bfp)
{
    return sam_bam_raw_write_get_bgzf_common(fp, bfp, 0);
}

static int sam_bam_raw_write_get_bgzf_indexed(htsFile *fp, BGZF **bfp)
{
    return sam_bam_raw_write_get_bgzf_common(fp, bfp, 1);
}

static int sam_bam_raw_write_frames_bgzf(BGZF *bfp, const void *data, size_t len)
{
    if (!bfp || !data || len == 0) {
        errno = EINVAL;
        return -2;
    }
#ifdef SSIZE_MAX
    if (len > (size_t)SSIZE_MAX) {
        errno = EOVERFLOW;
        return -2;
    }
#endif

    if (bfp->is_compressed && bgzf_flush_try(bfp, len) < 0)
        return -1;
    return bgzf_write_small(bfp, data, len) == (ssize_t)len ? 0 : -1;
}

static int sam_bam_raw_write_frames(htsFile *fp, const void *data, size_t len)
{
    BGZF *bfp;

    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0)
        return -2;
    return sam_bam_raw_write_frames_bgzf(bfp, data, len);
}

#define BAM_BATCH_PACKED_WRITE_MAX (1024 * 1024)

static int sam_bam_batch_write1_ranges_packed(BGZF *bfp,
                                              const bam_batch_t *batch,
                                              const int *ranges,
                                              int n_ranges,
                                              size_t total_len)
{
    uint8_t *buf, *dst;
    int ri;

    buf = malloc(total_len);
    if (!buf) {
        errno = ENOMEM;
        return -2;
    }
    dst = buf;
    for (ri = 0; ri < n_ranges; ri++) {
        int i, beg = ranges[ri * 2], end_i = ranges[ri * 2 + 1];

        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            memcpy(dst, rec->frame, rec->frame_len);
            dst += rec->frame_len;
        }
    }

    if (sam_bam_raw_write_frames_bgzf(bfp, buf, total_len) < 0) {
        free(buf);
        return -1;
    }
    free(buf);
    return 0;
}

static int sam_bam_raw_write_index_push(htsFile *fp, const sam_hdr_t *h,
                                        BGZF *bfp, bam_batch_record_t *record)
{
    hts_pos_t endpos;

    if (!fp->idx)
        return 0;
    if (sam_bam_batch_record_endpos(record, NULL, &endpos) < 0)
        return -2;
    if (bgzf_idx_push(bfp, fp->idx, record->core.tid, record->core.pos,
                      endpos, bgzf_tell(bfp),
                      !(record->core.flag & BAM_FUNMAP)) < 0) {
        int qname_len = record->core.l_qname > 0 ? record->core.l_qname - 1 : 0;
        const char *rname = h && record->core.tid >= 0
                            ? sam_hdr_tid2name(h, record->core.tid) : "*";
        hts_pos_t rlen = h && record->core.tid >= 0
                         ? sam_hdr_tid2len(h, record->core.tid) : 0;

        hts_log_error("Read '%.*s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                      qname_len, sam_bam_batch_record_qname(record),
                      rname ? rname : "*", rlen, record->core.flag,
                      record->core.pos + 1);
        return -1;
    }
    return 0;
}

static int sam_bam_raw_write_record_indexed(htsFile *fp, const sam_hdr_t *h,
                                            bam_batch_record_t *record)
{
    BGZF *bfp;
    hts_pos_t endpos;
    int ret;

    if (sam_bam_raw_write_get_bgzf_indexed(fp, &bfp) < 0)
        return -2;
    if (fp->idx && sam_bam_batch_record_endpos(record, NULL, &endpos) < 0)
        return -2;
    if (bfp->is_compressed && bgzf_flush_try(bfp, record->frame_len) < 0)
        return -1;
    if (fp->idx && !bfp->mt)
        hts_idx_amend_last(fp->idx, bgzf_tell(bfp));
    if (bgzf_write_small(bfp, record->frame, record->frame_len) !=
        (ssize_t)record->frame_len)
        return -1;
    ret = sam_bam_raw_write_index_push(fp, h, bfp, record);
    return ret < 0 ? ret : 0;
}

static int sam_bam_batch_write1_range_indexed(htsFile *fp, const sam_hdr_t *h,
                                              const bam_batch_t *batch,
                                              int beg, int end_i)
{
    int i;

    if (!batch || beg < 0 || end_i < beg || end_i > batch->n_records ||
        (end_i > beg && !batch->records)) {
        errno = EINVAL;
        return -2;
    }
    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];
        hts_pos_t endpos;

        if (sam_bam_batch_record_raw_write_safe(rec) != 0)
            return -2;
        if (fp->idx &&
            sam_bam_batch_record_endpos(rec, NULL, &endpos) < 0)
            return -2;
    }
    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];
        int raw_safe;

        raw_safe = sam_bam_raw_write_record_indexed(fp, h, rec);
        if (raw_safe < 0)
            return raw_safe;
    }
    return 0;
}

#define BAM_RAW_FLAG_OFFSET (4 + 14)

static void sam_bam_raw_store_frame_flag(uint8_t *frame, uint16_t flag)
{
    u16_to_le(flag, frame + BAM_RAW_FLAG_OFFSET);
}

static int sam_bam_batch_write1_range_flags_copy(
        htsFile *fp, const bam_batch_t *batch, int beg, int end_i,
        uint16_t set_flags, uint16_t clear_flags)
{
    uint8_t *buf = NULL, *dst;
    BGZF *bfp;
    size_t len = 0;
    int i, ret;

    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];

        if (sam_bam_batch_record_raw_layout_valid(rec) < 0)
            return -2;
        if (rec->frame_len > SIZE_MAX - len) {
            errno = EOVERFLOW;
            return -2;
        }
        len += rec->frame_len;
    }
#ifdef SSIZE_MAX
    if (len > (size_t)SSIZE_MAX) {
        errno = EOVERFLOW;
        return -2;
    }
#endif
    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0)
        return -2;

    buf = malloc(len);
    if (!buf) {
        errno = ENOMEM;
        return -2;
    }
    dst = buf;
    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];
        uint16_t flag = (uint16_t)((rec->core.flag | set_flags) & ~clear_flags);

        memcpy(dst, rec->frame, rec->frame_len);
        sam_bam_raw_store_frame_flag(dst, flag);
        dst += rec->frame_len;
    }

    if (bfp->is_compressed && bgzf_flush_try(bfp, len) < 0) {
        free(buf);
        return -1;
    }
    ret = bgzf_write_small(bfp, buf, len) == (ssize_t)len ? 0 : -1;
    free(buf);
    return ret;
}

int sam_bam_batch_record_write1(htsFile *fp, const sam_hdr_t *h,
                                bam_batch_record_t *record)
{
    int raw_safe, ret;

    (void)h;
    raw_safe = sam_bam_batch_record_raw_write_safe(record);
    if (raw_safe != 0)
        return -2;
    if (fp->idx) {
        ret = sam_bam_raw_write_record_indexed(fp, h, record);
        return ret < 0 ? ret : (int)record->frame_len;
    }

    ret = sam_bam_raw_write_frames(fp, record->frame, record->frame_len);
    return ret < 0 ? ret : (int)record->frame_len;
}

int sam_bam_batch_record_write1_aux_filtered(htsFile *fp, const sam_hdr_t *h,
                                             bam_batch_record_t *record,
                                             sam_bam_aux_filter_f filter,
                                             void *filter_data)
{
    const uint8_t *aux, *end, *s, *prefix;
    size_t prefix_len, new_aux_len = 0, block_len, frame_len;
    uint8_t block_len_le[4];
    uint8_t *keep_map = NULL;
    BGZF *bfp;
    int dropped = 0, raw_safe, ret, n_aux = 0, m_aux = 0, idx;

    (void)h;
    if (!filter)
        return sam_bam_batch_record_write1(fp, h, record);
    raw_safe = sam_bam_batch_record_raw_write_safe(record);
    if (raw_safe != 0)
        return -2;

    aux = sam_bam_batch_record_aux(record);
    end = record->body + record->raw_l_data;
    prefix = record->frame + 4;
    if (aux < record->body || aux > end) {
        errno = EINVAL;
        return -2;
    }
    prefix_len = (size_t)(aux - prefix);

    s = aux;
    while (s < end) {
        const uint8_t *value, *next;
        char tag[2];
        int keep;

        if (end - s < 3) {
            errno = EINVAL;
            free(keep_map);
            return -2;
        }
        tag[0] = s[0];
        tag[1] = s[1];
        value = s + 2;
        next = sam_bam_record_view_aux_skip(value, end);
        if (!next) {
            errno = EINVAL;
            free(keep_map);
            return -2;
        }
        keep = filter(tag, value, filter_data);
        if (keep < 0) {
            errno = EINVAL;
            free(keep_map);
            return -2;
        }
        if (n_aux == m_aux) {
            int new_m = m_aux ? (m_aux > INT_MAX / 2 ? INT_MAX : m_aux * 2) : 8;
            uint8_t *new_keep;

            if (new_m == m_aux ||
                (size_t)new_m > SIZE_MAX / sizeof(*keep_map)) {
                errno = ENOMEM;
                free(keep_map);
                return -2;
            }
            new_keep = realloc(keep_map, (size_t)new_m * sizeof(*keep_map));
            if (!new_keep) {
                errno = ENOMEM;
                free(keep_map);
                return -2;
            }
            keep_map = new_keep;
            m_aux = new_m;
        }
        keep_map[n_aux++] = keep ? 1 : 0;
        if (keep) {
            new_aux_len += (size_t)(next - s);
        } else {
            dropped = 1;
        }
        s = next;
    }

    if (!dropped) {
        free(keep_map);
        return sam_bam_batch_record_write1(fp, h, record);
    }
    if (fp->idx) {
        free(keep_map);
        return -2;
    }

    block_len = prefix_len + new_aux_len;
    frame_len = 4 + block_len;
    if (block_len > INT32_MAX || frame_len > INT_MAX) {
        errno = EOVERFLOW;
        free(keep_map);
        return -2;
    }
    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0) {
        free(keep_map);
        return -2;
    }
    if (bfp->is_compressed && bgzf_flush_try(bfp, frame_len) < 0) {
        free(keep_map);
        return -1;
    }

    u32_to_le((uint32_t)block_len, block_len_le);
    ret = bgzf_write_small(bfp, block_len_le, sizeof(block_len_le)) == 4 &&
          bgzf_write_small(bfp, prefix, prefix_len) == (ssize_t)prefix_len;
    for (s = aux, idx = 0; s < end; idx++) {
        const uint8_t *value = s + 2;
        const uint8_t *next = sam_bam_record_view_aux_skip(value, end);

        if (idx >= n_aux || !next) {
            ret = 0;
            break;
        }
        if (keep_map[idx]) {
            size_t len = (size_t)(next - s);
            ret = ret &&
                  bgzf_write_small(bfp, s, len) == (ssize_t)len;
        }
        s = next;
    }

    free(keep_map);
    return ret ? (int)frame_len : -1;
}

static int sam_bam_batch_record_sam_write_safe(bam_batch_record_t *record)
{
    int cg_fixup = 0;

    if (sam_bam_batch_record_raw_layout_valid(record) < 0)
        return -1;
    if (sam_bam_batch_record_qname(record)[record->core.l_qname - 1] != '\0')
        return 1;
    if (sam_bam_batch_record_cg_fixup_needed(record, &cg_fixup) < 0)
        return -1;
    return cg_fixup ? 1 : 0;
}

static int sam_bam_batch_record_validate_sam_cigar(
        const bam_batch_record_t *record)
{
    const uint8_t *cigar;
    hts_pos_t qlen = 0;
    int k;

    if (!record)
        return -1;
    cigar = sam_bam_batch_record_cigar(record);
    for (k = 0; k < record->core.n_cigar; k++) {
        uint32_t c = le_to_u32(cigar + ((size_t)k << 2));
        int op = bam_cigar_op(c);

        if (bam_cigar_type(op) & 1)
            qlen += bam_cigar_oplen(c);
    }
    if (record->core.n_cigar > 0 && record->core.l_qseq > 0 &&
        !(record->core.flag & BAM_FUNMAP) && qlen != record->core.l_qseq) {
        int qname_len = record->core.l_qname > 0 ? record->core.l_qname - 1 : 0;

        hts_log_error("CIGAR and query sequence lengths differ for %.*s",
                      qname_len, sam_bam_batch_record_qname(record));
        errno = EINVAL;
        return -1;
    }
    return 0;
}

static int sam_bam_batch_sam_validate_ranges(const bam_batch_t *batch,
                                             const int *ranges, int n_ranges)
{
    int r;

    if (!batch || n_ranges < 0 || (n_ranges > 0 && !ranges)) {
        errno = EINVAL;
        return -2;
    }
    for (r = 0; r < n_ranges; r++) {
        int beg = ranges[r * 2];
        int end_i = ranges[r * 2 + 1];
        int i;

        if (beg < 0 || end_i < beg || end_i > batch->n_records ||
            (end_i > beg && !batch->records)) {
            errno = EINVAL;
            return -2;
        }
        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            if (sam_bam_batch_record_sam_write_safe(rec) != 0)
                return -2;
        }
    }
    return 0;
}

static sp_bam_views *sam_bam_view_job_alloc(SAM_state *fd)
{
    sp_bam_views *gv = calloc(1, sizeof(*gv));

    if (!gv)
        return NULL;
    gv->a_records = SAM_NBAM;
    gv->records = calloc(gv->a_records, sizeof(*gv->records));
    if (!gv->records) {
        free(gv);
        return NULL;
    }
    gv->fd = fd;
    return gv;
}

static int sam_bam_view_job_append(sp_bam_views *gv,
                                   const bam_batch_record_t *rec)
{
    bam_batch_record_t *dst;
    uint8_t *new_data;
    size_t data_off, new_len, new_alloc;

    if (gv->n_records == gv->a_records) {
        int new_a = gv->a_records > INT_MAX / 2
                    ? INT_MAX : gv->a_records * 2;
        bam_batch_record_t *new_records;

        if (new_a == gv->a_records) {
            errno = ENOMEM;
            return -1;
        }
        new_records = realloc(gv->records,
                              (size_t)new_a * sizeof(*new_records));
        if (!new_records) {
            errno = ENOMEM;
            return -1;
        }
        gv->records = new_records;
        gv->a_records = new_a;
    }
    if ((size_t)rec->raw_l_data > SIZE_MAX - gv->data_len) {
        errno = EOVERFLOW;
        return -1;
    }
    new_len = gv->data_len + (size_t)rec->raw_l_data;
    if (new_len > gv->data_alloc) {
        new_alloc = gv->data_alloc ? gv->data_alloc : SAM_NBYTES;
        while (new_alloc < new_len) {
            if (new_alloc > SIZE_MAX / 2) {
                new_alloc = new_len;
                break;
            }
            new_alloc *= 2;
        }
        new_data = realloc(gv->data, new_alloc);
        if (!new_data) {
            errno = ENOMEM;
            return -1;
        }
        gv->data = new_data;
        gv->data_alloc = new_alloc;
    }

    data_off = gv->data_len;
    memcpy(gv->data + data_off, rec->body, rec->raw_l_data);
    gv->data_len = new_len;

    dst = &gv->records[gv->n_records++];
    *dst = *rec;
    dst->frame = NULL;
    dst->frame_len = 0;
    dst->body = gv->data + data_off;
    dst->voff_beg = dst->voff_end = 0;
    return 0;
}

static int sam_bam_view_job_dispatch(SAM_state *fd, sp_bam_views **job)
{
    sp_bam_views *gv = *job;
    int ret = 0;

    if (!gv || (!gv->owns_batch && gv->n_records == 0) ||
        (gv->owns_batch && gv->n_ranges == 0)) {
        sam_free_sp_bam_views(gv);
        *job = NULL;
        return 0;
    }
    gv->serial = fd->serial++;
    pthread_mutex_lock(&fd->command_m);
    if (fd->errcode != 0) {
        ret = -fd->errcode;
    } else if (hts_tpool_dispatch3(fd->p, fd->q, sam_format_batch_worker, gv,
                                   cleanup_sp_bam_views,
                                   cleanup_sp_lines, 0) < 0) {
        ret = -1;
    }
    pthread_mutex_unlock(&fd->command_m);
    if (ret < 0)
        return ret;
    *job = NULL;
    return 0;
}

static int sam_flush_pending_bams(SAM_state *fd)
{
    sp_bams *gb = fd->curr_bam;
    int ret = 0;

    if (!gb || gb->nbams == 0)
        return 0;
    gb->serial = fd->serial++;
    pthread_mutex_lock(&fd->command_m);
    if (fd->errcode != 0) {
        ret = -fd->errcode;
    } else if (hts_tpool_dispatch3(fd->p, fd->q, sam_format_worker, gb,
                                   cleanup_sp_bams,
                                   cleanup_sp_lines, 0) < 0) {
        ret = -1;
    }
    pthread_mutex_unlock(&fd->command_m);
    if (ret < 0)
        return ret;
    fd->curr_bam = NULL;
    return 0;
}

static int sam_bam_batch_write_sam_ranges_threaded(htsFile *fp,
                                                   const sam_hdr_t *h,
                                                   const bam_batch_t *batch,
                                                   const int *ranges,
                                                   int n_ranges)
{
    SAM_state *fd = (SAM_state *)fp->state;
    sp_bam_views *gv = NULL;
    int r, ret;

    if ((ret = sam_start_threaded_output(fp, h)) < 0)
        return ret;
    if (sam_flush_pending_bams(fd) < 0)
        return -1;

    for (r = 0; r < n_ranges; r++) {
        int beg = ranges[r * 2];
        int end_i = ranges[r * 2 + 1];
        int i;

        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            if (!gv && !(gv = sam_bam_view_job_alloc(fd)))
                return -1;
            if (sam_bam_view_job_append(gv, rec) < 0)
                goto err;
            if (gv->n_records == SAM_NBAM ||
                gv->data_len > SAM_NBYTES * 0.8) {
                if (sam_bam_view_job_dispatch(fd, &gv) < 0)
                    goto err;
            }
        }
    }
    if (sam_bam_view_job_dispatch(fd, &gv) < 0)
        goto err;
    return 1;

 err:
    sam_free_sp_bam_views(gv);
    return -1;
}

static int sam_bam_batch_write_sam_ranges_threaded_steal(htsFile *fp,
                                                         const sam_hdr_t *h,
                                                         bam_batch_t *batch,
                                                         const int *ranges,
                                                         int n_ranges)
{
    SAM_state *fd = (SAM_state *)fp->state;
    sp_bam_views *gv = NULL;
    int ret;

    if ((ret = sam_start_threaded_output(fp, h)) < 0)
        return ret;
    if (sam_flush_pending_bams(fd) < 0)
        return -1;

    gv = calloc(1, sizeof(*gv));
    if (!gv) {
        errno = ENOMEM;
        return -1;
    }
    gv->ranges = malloc((size_t)n_ranges * 2 * sizeof(*gv->ranges));
    if (!gv->ranges) {
        errno = ENOMEM;
        sam_free_sp_bam_views(gv);
        return -1;
    }
    memcpy(gv->ranges, ranges, (size_t)n_ranges * 2 * sizeof(*gv->ranges));
    gv->n_ranges = n_ranges;
    gv->batch = *batch;
    gv->owns_batch = 1;
    gv->fd = fd;
    memset(batch, 0, sizeof(*batch));

    if (sam_bam_view_job_dispatch(fd, &gv) < 0)
        goto err;
    return 1;

 err:
    sam_free_sp_bam_views(gv);
    return -1;
}

static int sam_bam_batch_write_sam_ranges_unthreaded(htsFile *fp,
                                                     const sam_hdr_t *h,
                                                     const bam_batch_t *batch,
                                                     const int *ranges,
                                                     int n_ranges)
{
    int r, wrote = 0;

    for (r = 0; r < n_ranges; r++) {
        int beg = ranges[r * 2];
        int end_i = ranges[r * 2 + 1];
        int i;

        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            fp->line.l = 0;
            if (sam_format_batch_record_append(h, rec, &fp->line) < 0)
                return -1;
            kputc('\n', &fp->line);
            if (fp->is_bgzf) {
                if (bgzf_flush_try(fp->fp.bgzf, fp->line.l) < 0)
                    return -1;
                if (bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l) !=
                    fp->line.l)
                    return -1;
            } else {
                if (hwrite(fp->fp.hfile, fp->line.s, fp->line.l) !=
                    fp->line.l)
                    return -1;
            }
            wrote = 1;
        }
    }
    return wrote;
}

int sam_bam_batch_write_sam_ranges(htsFile *fp, const sam_hdr_t *h,
                                   const bam_batch_t *batch,
                                   const int *ranges, int n_ranges)
{
    if (!fp || !fp->is_write || !h)
        return -2;
    if (n_ranges == 0)
        return 0;
    if (fp->idx)
        return -2;

    switch (fp->format.format) {
    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam:
        break;
    default:
        return -2;
    }

    if (sam_bam_batch_sam_validate_ranges(batch, ranges, n_ranges) < 0)
        return -2;
    if (fp->state)
        return sam_bam_batch_write_sam_ranges_threaded(fp, h, batch,
                                                       ranges, n_ranges);
    return sam_bam_batch_write_sam_ranges_unthreaded(fp, h, batch,
                                                     ranges, n_ranges);
}

int sam_bam_batch_write_sam_ranges_steal(htsFile *fp, const sam_hdr_t *h,
                                         bam_batch_t *batch,
                                         const int *ranges, int n_ranges)
{
    if (!fp || !fp->is_write || !h || !batch)
        return -2;
    if (n_ranges == 0)
        return 0;
    if (fp->idx)
        return -2;

    switch (fp->format.format) {
    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam:
        break;
    default:
        return -2;
    }

    if (sam_bam_batch_sam_validate_ranges(batch, ranges, n_ranges) < 0)
        return -2;
    if (fp->state)
        return sam_bam_batch_write_sam_ranges_threaded_steal(fp, h, batch,
                                                             ranges,
                                                             n_ranges);
    return sam_bam_batch_write_sam_ranges_unthreaded(fp, h, batch,
                                                     ranges, n_ranges);
}

int sam_bam_batch_write1(htsFile *fp, const sam_hdr_t *h,
                         const bam_batch_t *batch)
{
    const bam_batch_record_t *first, *last;
    const uint8_t *beg, *end, *expected;
    int i, j;

    (void)h;
    if (!batch)
        return -2;
    if (batch->n_records == 0)
        return 0;
    if (fp->idx)
        return sam_bam_batch_write1_range_indexed(fp, h, batch, 0,
                                                  batch->n_records);
    if (batch->n_segments > 0) {
        int rec_i = 0;

        if (!batch->segments || !batch->records) {
            errno = EINVAL;
            return -2;
        }
        for (i = 0; i < batch->n_segments; i++) {
            const bam_batch_segment_t *seg = &batch->segments[i];
            size_t off = 0;

            if (!seg->data || seg->len == 0) {
                errno = EINVAL;
                return -2;
            }
            while (off < seg->len) {
                bam_batch_record_t *rec;

                if (rec_i >= batch->n_records) {
                    errno = EINVAL;
                    return -2;
                }
                rec = &batch->records[rec_i];
                if (sam_bam_batch_record_raw_write_safe(rec) != 0 ||
                    rec->frame != seg->data + off ||
                    rec->frame_len > seg->len - off) {
                    errno = EINVAL;
                    return -2;
                }
                off += rec->frame_len;
                rec_i++;
            }
        }
        if (rec_i != batch->n_records) {
            errno = EINVAL;
            return -2;
        }
        for (i = 0; i < batch->n_segments; i++)
            if (sam_bam_raw_write_frames(fp, batch->segments[i].data,
                                         batch->segments[i].len) < 0)
                return -1;
        return 0;
    }
    if (!batch->data || !batch->records || batch->len == 0) {
        errno = EINVAL;
        return -2;
    }

    first = &batch->records[0];
    last = &batch->records[batch->n_records - 1];
    beg = first->frame;
    end = last->frame + last->frame_len;
    if (beg != batch->data || end != batch->data + batch->len || end < beg) {
        errno = EINVAL;
        return -2;
    }
    expected = batch->data;
    for (j = 0; j < batch->n_records; j++) {
        if (batch->records[j].frame != expected) {
            errno = EINVAL;
            return -2;
        }
        if (sam_bam_batch_record_raw_write_safe(&batch->records[j]) != 0)
            return -2;
        expected += batch->records[j].frame_len;
    }
    if (expected != batch->data + batch->len) {
        errno = EINVAL;
        return -2;
    }

    return sam_bam_raw_write_frames(fp, batch->data, batch->len);
}

int sam_bam_batch_write1_range(htsFile *fp, const sam_hdr_t *h,
                               const bam_batch_t *batch, int beg, int end_i)
{
    const uint8_t *run_beg = NULL, *expected = NULL;
    BGZF *bfp;
    size_t run_len = 0;
    int i;

    (void)h;
    if (!batch || beg < 0 || end_i < beg || end_i > batch->n_records ||
        (end_i > beg && !batch->records)) {
        errno = EINVAL;
        return -2;
    }
    if (beg == end_i)
        return 0;
    if (fp->idx)
        return sam_bam_batch_write1_range_indexed(fp, h, batch, beg, end_i);

    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];

        if (sam_bam_batch_record_raw_write_safe(rec) != 0)
            return -2;
        if (rec->frame_len > SIZE_MAX - run_len) {
            errno = EOVERFLOW;
            return -2;
        }
        run_len += rec->frame_len;
    }
    run_len = 0;
    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0)
        return -2;

    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];

        if (!run_beg) {
            run_beg = rec->frame;
            run_len = rec->frame_len;
        } else if (rec->frame == expected) {
            if (rec->frame_len > SIZE_MAX - run_len) {
                errno = EOVERFLOW;
                return -2;
            }
            run_len += rec->frame_len;
        } else {
            if (sam_bam_raw_write_frames_bgzf(bfp, run_beg, run_len) < 0)
                return -1;
            run_beg = rec->frame;
            run_len = rec->frame_len;
        }
        expected = rec->frame + rec->frame_len;
    }

    if (run_beg && sam_bam_raw_write_frames_bgzf(bfp, run_beg, run_len) < 0)
        return -1;
    return 0;
}

int sam_bam_batch_write1_ranges(htsFile *fp, const sam_hdr_t *h,
                                const bam_batch_t *batch,
                                const int *ranges, int n_ranges)
{
    const uint8_t *run_beg = NULL, *expected = NULL;
    BGZF *bfp;
    size_t run_len = 0, total_len = 0;
    int ri, last_end = 0;

    if (!batch || n_ranges < 0 || (n_ranges > 0 && !ranges) ||
        n_ranges > INT_MAX / 2 ||
        (batch->n_records > 0 && !batch->records)) {
        errno = EINVAL;
        return -2;
    }
    if (n_ranges == 0)
        return 0;
    if (fp->idx) {
        for (ri = 0; ri < n_ranges; ri++) {
            int beg = ranges[ri * 2], end_i = ranges[ri * 2 + 1];
            int i;

            if (beg < last_end || end_i <= beg || end_i > batch->n_records) {
                errno = EINVAL;
                return -2;
            }
            last_end = end_i;
            for (i = beg; i < end_i; i++) {
                bam_batch_record_t *rec = &batch->records[i];
                hts_pos_t endpos;

                if (sam_bam_batch_record_raw_write_safe(rec) != 0 ||
                    sam_bam_batch_record_endpos(rec, NULL, &endpos) < 0)
                    return -2;
            }
        }
        for (ri = 0; ri < n_ranges; ri++) {
            int beg = ranges[ri * 2], end_i = ranges[ri * 2 + 1];
            int wr = sam_bam_batch_write1_range_indexed(fp, h, batch,
                                                        beg, end_i);
            if (wr < 0)
                return wr;
        }
        return 0;
    }

    for (ri = 0; ri < n_ranges; ri++) {
        int beg = ranges[ri * 2], end_i = ranges[ri * 2 + 1];
        int i;

        if (beg < last_end || end_i <= beg || end_i > batch->n_records) {
            errno = EINVAL;
            return -2;
        }
        last_end = end_i;
        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            if (sam_bam_batch_record_raw_write_safe(rec) != 0)
                return -2;
            if (rec->frame_len > SIZE_MAX - total_len) {
                errno = EOVERFLOW;
                return -2;
            }
            total_len += rec->frame_len;
        }
    }
    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0)
        return -2;
    if (n_ranges > 1 && total_len > 0 &&
        total_len <= BAM_BATCH_PACKED_WRITE_MAX &&
        (!batch->len || total_len * 4 <= batch->len))
        return sam_bam_batch_write1_ranges_packed(bfp, batch, ranges,
                                                  n_ranges, total_len);

    for (ri = 0; ri < n_ranges; ri++) {
        int beg = ranges[ri * 2], end_i = ranges[ri * 2 + 1];
        int i;

        for (i = beg; i < end_i; i++) {
            bam_batch_record_t *rec = &batch->records[i];

            if (!run_beg) {
                run_beg = rec->frame;
                run_len = rec->frame_len;
            } else if (rec->frame == expected) {
                if (rec->frame_len > SIZE_MAX - run_len) {
                errno = EOVERFLOW;
                return -2;
            }
            run_len += rec->frame_len;
        } else {
            if (sam_bam_raw_write_frames_bgzf(bfp, run_beg, run_len) < 0)
                return -1;
            run_beg = rec->frame;
            run_len = rec->frame_len;
            }
            expected = rec->frame + rec->frame_len;
        }
    }

    if (run_beg && sam_bam_raw_write_frames_bgzf(bfp, run_beg, run_len) < 0)
        return -1;
    return 0;
}

int sam_bam_batch_write1_range_flags(htsFile *fp, const sam_hdr_t *h,
                                     const bam_batch_t *batch, int beg,
                                     int end_i, uint16_t set_flags,
                                     uint16_t clear_flags)
{
    BGZF *bfp;
    const uint8_t *run_beg = NULL, *expected = NULL;
    size_t run_len = 0;
    int i, ret = 0, all_mutable = 1;

    (void)h;
    if (!set_flags && !clear_flags)
        return sam_bam_batch_write1_range(fp, h, batch, beg, end_i);
    if (fp->idx)
        return -2;
    if (!batch || beg < 0 || end_i < beg || end_i > batch->n_records ||
        (end_i > beg && !batch->records)) {
        errno = EINVAL;
        return -2;
    }
    if (beg == end_i)
        return 0;

    for (i = beg; i < end_i; i++) {
        bam_batch_record_t *rec = &batch->records[i];

        if (sam_bam_batch_record_raw_write_safe(rec) != 0)
            return -2;
        if (rec->frame_len > SIZE_MAX - run_len) {
            errno = EOVERFLOW;
            return -2;
        }
        run_len += rec->frame_len;
        if (!(rec->flags & BAM_BATCH_RECORD_F_OWNED))
            all_mutable = 0;
    }
    run_len = 0;
    if (!all_mutable)
        return sam_bam_batch_write1_range_flags_copy(
            fp, batch, beg, end_i, set_flags, clear_flags);
    if (sam_bam_raw_write_get_bgzf(fp, &bfp) < 0)
        return -2;

    for (i = beg; i < end_i; i++) {
        const bam_batch_record_t *rec = &batch->records[i];
        uint16_t flag = (uint16_t)((rec->core.flag | set_flags) & ~clear_flags);

        sam_bam_raw_store_frame_flag((uint8_t *)rec->frame, flag);
    }

    for (i = beg; i < end_i; i++) {
        const bam_batch_record_t *rec = &batch->records[i];

        if (!run_beg) {
            run_beg = rec->frame;
            run_len = rec->frame_len;
        } else if (rec->frame == expected) {
            if (rec->frame_len > SIZE_MAX - run_len) {
                errno = EOVERFLOW;
                ret = -2;
                break;
            }
            run_len += rec->frame_len;
        } else {
            if (sam_bam_raw_write_frames(fp, run_beg, run_len) < 0) {
                ret = -1;
                break;
            }
            run_beg = rec->frame;
            run_len = rec->frame_len;
        }
        expected = rec->frame + rec->frame_len;
    }

    if (ret == 0 && run_beg &&
        sam_bam_raw_write_frames(fp, run_beg, run_len) < 0)
        ret = -1;

    for (i = beg; i < end_i; i++)
        sam_bam_raw_store_frame_flag((uint8_t *)batch->records[i].frame,
                                     (uint16_t)batch->records[i].core.flag);
    return ret;
}

int sam_bam_batch_record_query_len(const bam_batch_record_t *record,
                                   bam1_t *scratch, int include_hard_clip,
                                   hts_pos_t *query_len)
{
    const uint8_t *cigar;
    hts_pos_t cigar_qlen = 0, filter_qlen = 0;
    int is_cg, k, n_cigar;

    if (!record)
        return -1;

    is_cg = sam_bam_batch_record_cg_candidate(record);
    if (is_cg) {
        if (!scratch || sam_bam_batch_record_to_bam1(record, scratch) < 0)
            return -1;
        cigar = (const uint8_t *)bam_get_cigar(scratch);
        n_cigar = scratch->core.n_cigar;
    } else {
        cigar = sam_bam_batch_record_cigar(record);
        n_cigar = record->core.n_cigar;
    }

    for (k = 0; k < n_cigar; k++) {
        uint32_t c = is_cg ? bam_get_cigar(scratch)[k]
                           : le_to_u32(cigar + ((size_t)k << 2));
        int op = bam_cigar_op(c);

        if (bam_cigar_type(op) & 1)
            cigar_qlen += bam_cigar_oplen(c);
        if (query_len && ((bam_cigar_type(op) & 1) ||
                          (include_hard_clip && op == BAM_CHARD_CLIP)))
            filter_qlen += bam_cigar_oplen(c);
    }

    if (record->core.n_cigar > 0 && record->core.l_qseq > 0 &&
        !(record->core.flag & BAM_FUNMAP) &&
        cigar_qlen != record->core.l_qseq) {
        int qname_len = record->core.l_qname > 0 ? record->core.l_qname - 1 : 0;
        hts_log_error("CIGAR and query sequence lengths differ for %.*s",
                      qname_len, sam_bam_batch_record_qname(record));
        return -1;
    }
    if (query_len)
        *query_len = filter_qlen;
    return 0;
}

int sam_bam_batch_record_endpos(const bam_batch_record_t *record,
                                bam1_t *scratch, hts_pos_t *endpos)
{
    hts_pos_t rlen = 0;
    hts_pos_t qlen = 0;
    int k;

    if (!record || !endpos)
        return -1;
    if (record->flags & BAM_BATCH_RECORD_F_ENDPOS_VALID) {
        *endpos = record->endpos;
        return 0;
    }

    if (record->core.flag & BAM_FUNMAP) {
        *endpos = record->core.pos + 1;
        return 0;
    }

    if (sam_bam_batch_record_cg_candidate(record) &&
        !(record->flags & BAM_BATCH_RECORD_F_RAW_WRITE_SAFE)) {
        if (!scratch || sam_bam_batch_record_to_bam1(record, scratch) < 0)
            return -1;
        rlen = bam_cigar2rlen(scratch->core.n_cigar, bam_get_cigar(scratch));
    } else {
        const uint8_t *cigar = sam_bam_batch_record_cigar(record);

        for (k = 0; k < record->core.n_cigar; k++) {
            uint32_t c = le_to_u32(cigar + ((size_t)k << 2));
            int op = bam_cigar_op(c);
            if (bam_cigar_type(op) & 2)
                rlen += bam_cigar_oplen(c);
            if (bam_cigar_type(op) & 1)
                qlen += bam_cigar_oplen(c);
        }
        if (record->core.n_cigar > 0 && record->core.l_qseq > 0 &&
            !(record->core.flag & BAM_FUNMAP) &&
            qlen != record->core.l_qseq) {
            int qname_len = record->core.l_qname > 0 ? record->core.l_qname - 1 : 0;
            hts_log_error("CIGAR and query sequence lengths differ for %.*s",
                          qname_len, sam_bam_batch_record_qname(record));
            return -1;
        }
    }
    if (rlen == 0)
        rlen = 1;
    *endpos = record->core.pos + rlen;
    return 0;
}

static bam_stream_reader_t *bam_batch_reader_get(htsFile *fp, int direct)
{
    if (!fp || !fp->is_bgzf || !fp->fp.bgzf)
        return NULL;

    if (fp->state) {
        uint32_t magic = *(uint32_t *)fp->state;

        if (magic == BAM_STREAM_READER_MAGIC)
            return (bam_stream_reader_t *)fp->state;
        if (magic == BAM_BATCH_REQUEST_MAGIC) {
            bam_batch_request_t *req = (bam_batch_request_t *)fp->state;
            bam_stream_reader_t *reader;

            fp->state = NULL;
            free(req);
            reader = bam_stream_reader_open(fp, NULL, 1);
            if (reader) {
                fp->state = reader;
                return reader;
            }
            return NULL;
        }
        if (magic == BAM_DEFERRED_THREADS_MAGIC) {
            bam_deferred_threads_t *cfg =
                (bam_deferred_threads_t *)fp->state;
            bam_stream_reader_t *reader = bam_stream_reader_open(fp, cfg, 1);

            if (reader) {
                fp->state = NULL;
                bam_deferred_threads_destroy(cfg);
                fp->state = reader;
                return reader;
            }
            if (!bam_batch_env_strict())
                (void)bam_deferred_threads_enable_bgzf(fp, cfg);
            fp->state = NULL;
            bam_deferred_threads_destroy(cfg);
            return NULL;
        }
        return NULL;
    }

    if (direct || bam_batch_env_enabled()) {
        bam_stream_reader_t *reader = bam_stream_reader_open(fp, NULL, 1);

        if (reader) {
            fp->state = reader;
            return reader;
        }
    }
    return NULL;
}

static int bam_batch_detach_ref_data(bam_stream_parse_job_t *job)
{
    const uint8_t *old_data = job->ref_data;
    uint8_t *data;
    int i;

    if (job->data || job->owned_block_result)
        return 0;
    if (!old_data || !job->len)
        return -2;

    data = malloc(job->len);
    if (!data) {
        errno = ENOMEM;
        return -2;
    }
    memcpy(data, old_data, job->len);

    for (i = 0; i < job->n_records; i++) {
        size_t frame_off = (size_t)(job->views[i].frame - old_data);
        size_t body_off = (size_t)(job->views[i].body - old_data);

        job->views[i].frame = data + frame_off;
        job->views[i].body = data + body_off;
        job->views[i].flags |= BAM_BATCH_RECORD_F_OWNED;
    }

    job->data = data;
    job->cap = job->len;
    job->ref_data = NULL;
    return 0;
}

int sam_bam_prepare_batch_reader(htsFile *fp)
{
    bam_batch_request_t *req;

    if (!fp || fp->format.format != bam || fp->is_write)
        return -2;
    if (fp->state) {
        uint32_t magic = *(uint32_t *)fp->state;

        return magic == BAM_BATCH_REQUEST_MAGIC ||
               magic == BAM_DEFERRED_THREADS_MAGIC ||
               magic == BAM_STREAM_READER_MAGIC ? 0 : -2;
    }

    req = calloc(1, sizeof(*req));
    if (!req)
        return -1;
    req->magic = BAM_BATCH_REQUEST_MAGIC;
    fp->state = req;
    return 0;
}

static void sam_bam_batch_attach_job(bam_batch_t *batch,
                                     bam_stream_parse_job_t *job)
{
    int i;

    if (job->data || job->owned_block_result)
        for (i = 0; i < job->n_records; i++)
            job->views[i].flags |= BAM_BATCH_RECORD_F_OWNED;
    batch->data = job->data ? job->data : job->ref_data;
    batch->len = job->len;
    batch->n_records = job->n_records;
    batch->records = job->views;
    batch->n_segments = job->n_segments;
    batch->segments = job->segments;
    batch->impl = job;
}

static int sam_bam_read_batch_upto(htsFile *fp, sam_hdr_t *h,
                                   bam_batch_t *batch, uint64_t limit_voff,
                                   int *hit_limit, int need_voff)
{
    bam_stream_reader_t *reader;
    bam_stream_parse_job_t *job = NULL;
    int ret;

    if (!batch)
        return -2;
    sam_bam_batch_destroy(batch);

    if (!fp || fp->format.format != bam || fp->is_write)
        return -2;

    reader = bam_batch_reader_get(fp, 1);
    if (!reader)
        return -2;

    if (reader->fused_q) {
        ret = bam_stream_reader_next_fused_job(reader, h, &job, limit_voff,
                                               hit_limit, need_voff);
        if (ret <= 0)
            return ret == 0 ? -1 : ret;

        sam_bam_batch_attach_job(batch, job);
        return job->n_records;
    }

    if (reader->parse_q && reader->parse_views) {
        ret = bam_stream_reader_next_view_job(reader, h, &job, limit_voff,
                                              hit_limit, need_voff);
        if (ret <= 0)
            return ret == 0 ? -1 : ret;

        sam_bam_batch_attach_job(batch, job);
        return job->n_records;
    }

    if (reader->pending_frame_error) {
        ret = reader->pending_frame_error;
        reader->pending_frame_error = 0;
        reader->parse_input_eof = 1;
        if (ret == -3)
            errno = ERANGE;
        return ret;
    }
    if (reader->parse_input_error)
        return reader->parse_input_error;
    if (reader->parse_input_eof)
        return -1;

    ret = bam_stream_reader_build_parse_job(reader, &job, 1, h, limit_voff,
                                            hit_limit, need_voff, 0);
    if (ret <= 0) {
        if (ret < 0) {
            reader->parse_input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
        }
        return ret == 0 ? -1 : ret;
    }

    if (!job->views) {
        bam_stream_parse_job_free(job);
        return -2;
    }
    if (bam_batch_detach_ref_data(job) < 0) {
        bam_stream_parse_job_free(job);
        return -2;
    }

    sam_bam_batch_attach_job(batch, job);
    return job->n_records;
}

int sam_bam_read_batch(htsFile *fp, sam_hdr_t *h, bam_batch_t *batch)
{
    return sam_bam_read_batch_upto(fp, h, batch, UINT64_MAX, NULL, 0);
}

int sam_bam_read_batch_voff(htsFile *fp, sam_hdr_t *h, bam_batch_t *batch)
{
    return sam_bam_read_batch_upto(fp, h, batch, UINT64_MAX, NULL, 1);
}

static int bam_count_job_validate_decode(bam_stream_parse_job_t *job,
                                         sam_hdr_t *h, int *n_valid)
{
    const uint8_t *data;
    bam1_t *scratch = NULL;
    size_t pos = 0;
    int n = 0;
    int ret = 0;

    *n_valid = 0;
    if (!job)
        return -4;

    if (job->views) {
        for (n = 0; n < job->n_records; n++) {
            int needs_materialize = 0;

            if (!bam_batch_record_view_tid_valid(&job->views[n], h)) {
                errno = ERANGE;
                ret = -3;
                goto cleanup;
            }
            if (sam_bam_batch_record_decode_status(&job->views[n],
                                                   &needs_materialize) < 0) {
                ret = -4;
                goto cleanup;
            }
            if (needs_materialize) {
                if (!scratch) {
                    scratch = bam_init1();
                    if (!scratch)
                        return -2;
                }
                if (sam_bam_batch_record_materialize_validated(&job->views[n],
                                                               scratch,
                                                               NULL) < 0) {
                    ret = -4;
                    goto cleanup;
                }
            }
            *n_valid = n + 1;
        }
        goto cleanup;
    }

    data = job->data ? job->data : job->ref_data;
    if (!data && job->len)
        return -4;

    while (pos + 4 <= job->len) {
        int32_t block_len = le_to_i32(data + pos);
        size_t frame_len;
        bam_batch_record_t view;
        int needs_materialize = 0;

        if (block_len < 32) {
            ret = -4;
            goto cleanup;
        }
        frame_len = 4 + (size_t)block_len;
        if (frame_len > job->len - pos) {
            ret = -4;
            goto cleanup;
        }
        if (n >= job->n_records) {
            ret = -4;
            goto cleanup;
        }

        bam_batch_record_view_set(&view, data + pos, frame_len, 0, 0);
        if (!bam_batch_record_view_tid_valid(&view, h)) {
            errno = ERANGE;
            ret = -3;
            goto cleanup;
        }
        if (sam_bam_batch_record_decode_status(&view, &needs_materialize) < 0) {
            ret = -4;
            goto cleanup;
        }
        if (needs_materialize) {
            if (!scratch) {
                scratch = bam_init1();
                if (!scratch)
                    return -2;
            }
            if (sam_bam_batch_record_materialize_validated(&view, scratch,
                                                           NULL) < 0) {
                ret = -4;
                goto cleanup;
            }
        }

        n++;
        *n_valid = n;
        pos += frame_len;
    }

    if (pos != job->len || n != job->n_records)
        ret = -4;

cleanup:
    bam_destroy1(scratch);
    return ret;
}

int sam_bam_read_batch_count(htsFile *fp, sam_hdr_t *h, int *n_records)
{
    bam_stream_reader_t *reader;
    bam_stream_parse_job_t *job = NULL;
    int n_valid = 0, ret;

    if (!n_records)
        return -2;
    *n_records = 0;

    if (!fp || fp->format.format != bam || fp->is_write)
        return -2;

    reader = bam_batch_reader_get(fp, 1);
    if (!reader)
        return -2;

    if (reader->fused_q) {
        ret = bam_stream_reader_next_fused_job(reader, h, &job, UINT64_MAX,
                                               NULL, 0);
        if (ret <= 0)
            return ret == 0 ? -1 : ret;
        ret = bam_count_job_validate_decode(job, h, &n_valid);
        if (ret < 0) {
            bam_stream_parse_job_free(job);
            if (n_valid > 0) {
                reader->pending_frame_error = ret;
                reader->input_eof = 1;
                *n_records = n_valid;
                return n_valid;
            }
            reader->input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
            return ret;
        }
        *n_records = job->n_records;
        ret = job->n_records;
        bam_stream_parse_job_free(job);
        return ret;
    }

    if (reader->pending_frame_error) {
        ret = reader->pending_frame_error;
        reader->pending_frame_error = 0;
        reader->parse_input_eof = 1;
        if (ret == -3)
            errno = ERANGE;
        return ret;
    }
    if (reader->parse_input_error)
        return reader->parse_input_error;
    if (reader->parse_input_eof)
        return -1;

    ret = bam_stream_reader_build_parse_job(reader, &job, 0, h, UINT64_MAX,
                                            NULL, 0, 0);
    if (ret <= 0) {
        if (ret < 0) {
            reader->parse_input_eof = 1;
            if (ret == -3)
                errno = ERANGE;
        }
        return ret == 0 ? -1 : ret;
    }

    ret = bam_count_job_validate_decode(job, h, &n_valid);
    if (ret < 0) {
        bam_stream_parse_job_free(job);
        if (n_valid > 0) {
            reader->pending_frame_error = ret;
            *n_records = n_valid;
            return n_valid;
        }
        reader->parse_input_eof = 1;
        if (ret == -3)
            errno = ERANGE;
        return ret;
    }

    *n_records = job->n_records;
    ret = job->n_records;
    bam_stream_parse_job_free(job);
    return ret;
}

static int sam_bam_restore_deferred_threads(htsFile *fp,
                                            const bam_deferred_threads_t *src)
{
    bam_deferred_threads_t *cfg;

    if (!src || src->n_threads <= 1)
        return 0;

    cfg = calloc(1, sizeof(*cfg));
    if (!cfg)
        return -2;
    *cfg = *src;
    cfg->magic = BAM_DEFERRED_THREADS_MAGIC;
    fp->state = cfg;
    return 0;
}

int sam_bam_batch_seek(htsFile *fp, uint64_t voff)
{
    bam_deferred_threads_t cfg = {0};
    int have_cfg = 0;

    if (!fp || !fp->is_bgzf || !fp->fp.bgzf)
        return -2;

    if (fp->state) {
        uint32_t magic = *(uint32_t *)fp->state;

        if (magic == BAM_STREAM_READER_MAGIC) {
            bam_stream_reader_t *reader = (bam_stream_reader_t *)fp->state;

            if (reader->pool) {
                cfg.magic = BAM_DEFERRED_THREADS_MAGIC;
                cfg.n_threads = hts_tpool_size(reader->pool);
                cfg.qsize = reader->qsize;
                cfg.pool = reader->own_pool ? NULL : reader->pool;
                cfg.use_pool = reader->own_pool ? 0 : 1;
                have_cfg = cfg.n_threads > 1;
            }
            fp->state = NULL;
            bam_stream_reader_destroy(reader);
        } else if (magic == BAM_DEFERRED_THREADS_MAGIC) {
            bam_deferred_threads_t *old =
                (bam_deferred_threads_t *)fp->state;

            cfg = *old;
            have_cfg = cfg.n_threads > 1;
            fp->state = NULL;
            bam_deferred_threads_destroy(old);
        } else if (magic == BAM_BATCH_REQUEST_MAGIC) {
            bam_batch_request_destroy(fp);
        } else {
            return -2;
        }
    }

    if (bgzf_seek(fp->fp.bgzf, (int64_t)voff, SEEK_SET) < 0)
        return -2;
    if (have_cfg && sam_bam_restore_deferred_threads(fp, &cfg) < 0)
        return -2;
    return 0;
}

static int bam_batch_compare_regions_by_tid(const void *av, const void *bv)
{
    const hts_reglist_t *a = (const hts_reglist_t *)av;
    const hts_reglist_t *b = (const hts_reglist_t *)bv;

    if (a->tid < 0 && b->tid >= 0)
        return 1;
    if (a->tid >= 0 && b->tid < 0)
        return -1;
    return (a->tid > b->tid) - (a->tid < b->tid);
}

static int sam_bam_itr_read_rest_batch(htsFile *fp, hts_itr_t *iter,
                                       sam_hdr_t *h, bam_batch_t *batch)
{
    int ret;

    if (iter->curr_off) {
        if (sam_bam_batch_seek(fp, iter->curr_off) < 0)
            return -2;
        iter->curr_off = 0;
    }

    ret = sam_bam_read_batch_upto(fp, h, batch, UINT64_MAX, NULL, 1);
    if (ret < 0) {
        iter->finished = 1;
        return ret;
    }
    if (batch->n_records > 0) {
        const bam_batch_record_t *rec = &batch->records[batch->n_records - 1];
        hts_pos_t endpos = 0;
        bam1_t *scratch = bam_init1();

        if (!scratch ||
            sam_bam_batch_record_endpos(rec, scratch, &endpos) < 0) {
            bam_destroy1(scratch);
            sam_bam_batch_destroy(batch);
            return -2;
        }
        bam_destroy1(scratch);
        iter->curr_tid = rec->core.tid;
        iter->curr_beg = rec->core.pos;
        iter->curr_end = endpos;
        iter->curr_off = rec->voff_end;
    }
    return ret;
}

static int sam_bam_itr_next_single_batch(htsFile *fp, hts_itr_t *iter,
                                         sam_hdr_t *h, bam_batch_t *batch)
{
    bam1_t *scratch = NULL;
    int ret = -1;

    for (;;) {
        bam_batch_t raw = {0};
        bam_batch_record_t *views;
        int hit_limit = 0, keep = 0, i;

        if (iter->curr_off == 0 || iter->i < 0 ||
            iter->curr_off >= iter->off[iter->i].v) {
            if (iter->i == iter->n_off - 1) {
                ret = -1;
                break;
            }
            if (iter->i < 0 ||
                iter->off[iter->i].v != iter->off[iter->i + 1].u) {
                if (sam_bam_batch_seek(fp, iter->off[iter->i + 1].u) < 0) {
                    bam_destroy1(scratch);
                    return -2;
                }
                iter->curr_off = iter->off[iter->i + 1].u;
            }
            ++iter->i;
        }

        ret = sam_bam_read_batch_upto(fp, h, &raw, iter->off[iter->i].v,
                                      &hit_limit, 1);
        if (ret < 0) {
            if (hit_limit) {
                iter->curr_off = iter->off[iter->i].v;
                continue;
            }
            break;
        }

        views = (bam_batch_record_t *)raw.records;
        for (i = 0; i < raw.n_records; i++) {
            const bam_batch_record_t *rec = &raw.records[i];
            hts_pos_t endpos;

            if (!scratch) {
                scratch = bam_init1();
                if (!scratch) {
                    sam_bam_batch_destroy(&raw);
                    return -2;
                }
            }
            if (sam_bam_batch_record_endpos(rec, scratch, &endpos) < 0) {
                sam_bam_batch_destroy(&raw);
                bam_destroy1(scratch);
                return -2;
            }

            iter->curr_off = rec->voff_end;
            if (rec->core.tid != iter->tid || rec->core.pos >= iter->end) {
                iter->finished = 1;
                break;
            }
            if (endpos > iter->beg && iter->end > rec->core.pos) {
                iter->curr_tid = rec->core.tid;
                iter->curr_beg = rec->core.pos;
                iter->curr_end = endpos;
                views[keep++] = *rec;
            }
        }

        if (hit_limit)
            iter->curr_off = iter->off[iter->i].v;
        if (keep > 0) {
            raw.n_records = keep;
            *batch = raw;
            bam_destroy1(scratch);
            return keep;
        }
        sam_bam_batch_destroy(&raw);
        if (iter->finished) {
            ret = -1;
            break;
        }
    }

    iter->finished = 1;
    bam_destroy1(scratch);
    return ret;
}

static int sam_bam_itr_multi_no_coord_batch(htsFile *fp, hts_itr_t *iter,
                                            sam_hdr_t *h, bam_batch_t *batch)
{
    for (;;) {
        bam_batch_t raw = {0};
        bam_batch_record_t *views;
        int i, keep = 0, ret;

        ret = sam_bam_read_batch_upto(fp, h, &raw, UINT64_MAX, NULL, 1);
        if (ret < 0) {
            iter->finished = 1;
            return ret;
        }

        views = (bam_batch_record_t *)raw.records;
        for (i = 0; i < raw.n_records; i++) {
            const bam_batch_record_t *rec = &raw.records[i];

            if (rec->core.tid >= 0)
                continue;
            views[keep++] = *rec;
        }

        if (keep > 0) {
            const bam_batch_record_t *rec = &views[keep - 1];

            raw.n_records = keep;
            iter->read_rest = 1;
            iter->curr_off = 0;
            iter->curr_tid = rec->core.tid;
            iter->curr_beg = rec->core.pos;
            iter->curr_end = rec->core.pos + 1;
            *batch = raw;
            return keep;
        }
        sam_bam_batch_destroy(&raw);
    }
}

static int sam_bam_itr_next_multi_batch(htsFile *fp, hts_itr_t *iter,
                                        sam_hdr_t *h, bam_batch_t *batch)
{
    bam1_t *scratch = NULL;
    int ret = -1, next_range = 0;

    for (;;) {
        bam_batch_t raw = {0};
        bam_batch_record_t *views;
        int hit_limit = 0, keep = 0, i;

        if (next_range || iter->curr_off == 0 || iter->i < 0 ||
            iter->i >= iter->n_off ||
            iter->curr_off >= iter->off[iter->i].v ||
            (iter->off[iter->i].max >> 32 == (uint64_t)iter->curr_tid &&
             (iter->off[iter->i].max & 0xffffffff) <
             (uint32_t)iter->curr_intv)) {
            do {
                iter->i++;
            } while (iter->i < iter->n_off &&
                     (iter->curr_off >= iter->off[iter->i].v ||
                      (iter->off[iter->i].max >> 32 ==
                       (uint64_t)iter->curr_tid &&
                       (iter->off[iter->i].max & 0xffffffff) <
                       (uint32_t)iter->curr_intv)));

            if (iter->i >= iter->n_off) {
                if (iter->nocoor) {
                    if (sam_bam_batch_seek(fp, iter->nocoor_off) < 0) {
                        bam_destroy1(scratch);
                        return -2;
                    }
                    iter->nocoor = 0;
                    ret = sam_bam_itr_multi_no_coord_batch(fp, iter, h, batch);
                    bam_destroy1(scratch);
                    return ret;
                }
                ret = -1;
                break;
            }

            if (iter->curr_off < iter->off[iter->i].u || next_range) {
                iter->curr_off = iter->off[iter->i].u;
                if (sam_bam_batch_seek(fp, iter->curr_off) < 0) {
                    bam_destroy1(scratch);
                    return -2;
                }
                next_range = 0;
            }
        }

        ret = sam_bam_read_batch_upto(fp, h, &raw, iter->off[iter->i].v,
                                      &hit_limit, 1);
        if (ret < 0) {
            if (hit_limit) {
                iter->curr_off = iter->off[iter->i].v;
                continue;
            }
            break;
        }

        views = (bam_batch_record_t *)raw.records;
        for (i = 0; i < raw.n_records; i++) {
            const bam_batch_record_t *rec = &raw.records[i];
            hts_reglist_t key, *found_reg;
            hts_pos_t endpos;
            int cr, ci, j;

            if (!scratch) {
                scratch = bam_init1();
                if (!scratch) {
                    sam_bam_batch_destroy(&raw);
                    return -2;
                }
            }
            if (sam_bam_batch_record_endpos(rec, scratch, &endpos) < 0) {
                sam_bam_batch_destroy(&raw);
                bam_destroy1(scratch);
                return -2;
            }

            iter->curr_off = rec->voff_end;
            if (rec->core.tid != iter->curr_tid) {
                key.tid = rec->core.tid;
                found_reg = (hts_reglist_t *)bsearch(
                        &key, iter->reg_list, iter->n_reg,
                        sizeof(hts_reglist_t),
                        bam_batch_compare_regions_by_tid);
                if (!found_reg)
                    continue;

                iter->curr_reg = (int)(found_reg - iter->reg_list);
                iter->curr_tid = rec->core.tid;
                iter->curr_intv = 0;
            }

            cr = iter->curr_reg;
            ci = iter->curr_intv;
            for (j = ci; j < iter->reg_list[cr].count; j++) {
                if (endpos > iter->reg_list[cr].intervals[j].beg &&
                    iter->reg_list[cr].intervals[j].end > rec->core.pos) {
                    iter->curr_beg = rec->core.pos;
                    iter->curr_end = endpos;
                    iter->curr_intv = j;
                    views[keep++] = *rec;
                    break;
                }

                if (rec->core.pos > iter->reg_list[cr].intervals[j].end)
                    iter->curr_intv = j + 1;
                if (endpos < iter->reg_list[cr].intervals[j].beg)
                    break;
            }
        }

        if (hit_limit)
            iter->curr_off = iter->off[iter->i].v;
        if (keep > 0) {
            raw.n_records = keep;
            *batch = raw;
            bam_destroy1(scratch);
            return keep;
        }
        sam_bam_batch_destroy(&raw);
    }

    iter->finished = 1;
    bam_destroy1(scratch);
    return ret;
}

int sam_bam_itr_next_batch(htsFile *fp, hts_itr_t *iter, sam_hdr_t *h,
                           bam_batch_t *batch)
{
    if (!batch)
        return -2;
    sam_bam_batch_destroy(batch);

    if (!fp || !iter || iter->finished)
        return -1;
    if (fp->format.format != bam || !fp->is_bgzf || !fp->fp.bgzf)
        return -2;

    if (iter->read_rest)
        return sam_bam_itr_read_rest_batch(fp, iter, h, batch);

    if (iter->multi)
        return sam_bam_itr_next_multi_batch(fp, iter, h, batch);

    if (!iter->off)
        return -1;
    return sam_bam_itr_next_single_batch(fp, iter, h, batch);
}

int sam_set_thread_pool(htsFile *fp, htsThreadPool *p) {
    if (fp->format.format == bam) {
        if (!fp->is_write &&
            (bam_batch_state_is_request(fp) || bam_batch_env_enabled() ||
             bam_batch_fused_env_enabled() || bam_stream_env_enabled() ||
             bam_ordered_env_enabled()))
        {
            bam_batch_request_destroy(fp);
            return bam_deferred_threads_set(fp, 0, p);
        }
        if (fp->format.compression == bgzf)
            return bgzf_thread_pool(fp->fp.bgzf, p->pool, p->qsize);
        return 0;
    }

    if (fp->state)
        return -2;   //already exists!

    if (!(fp->state = sam_state_create(fp)))
        return -1;
    SAM_state *fd = (SAM_state *)fp->state;

    pthread_mutex_init(&fd->lines_m, NULL);
    pthread_mutex_init(&fd->command_m, NULL);
    pthread_cond_init(&fd->command_c, NULL);
    fd->p = p->pool;
    int qsize = p->qsize;
    if (!qsize)
        qsize = 2*hts_tpool_size(fd->p);
    fd->q = hts_tpool_process_init(fd->p, qsize, 0);
    if (!fd->q) {
        sam_state_destroy(fp);
        return -1;
    }

    if (fp->format.compression == bgzf)
        return bgzf_thread_pool(fp->fp.bgzf, p->pool, p->qsize);

    return 0;
}

int sam_set_threads(htsFile *fp, int nthreads) {
    if (nthreads <= 0)
        return 0;

    if (fp->format.format == bam) {
        if (!fp->is_write &&
            (bam_batch_state_is_request(fp) || bam_batch_env_enabled() ||
             bam_batch_fused_env_enabled() || bam_stream_env_enabled() ||
             bam_ordered_env_enabled()))
        {
            bam_batch_request_destroy(fp);
            return bam_deferred_threads_set(fp, nthreads, NULL);
        }
        if (fp->format.compression == bgzf)
            return bgzf_mt(fp->fp.bgzf, nthreads, 256);
        return 0;
    }

    htsThreadPool p;
    p.pool = hts_tpool_init(nthreads);
    p.qsize = nthreads*2;

    int ret = sam_set_thread_pool(fp, &p);
    if (ret < 0) {
        if (p.pool)
            hts_tpool_destroy(p.pool);
        return ret;
    }

    SAM_state *fd = (SAM_state *)fp->state;
    fd->own_pool = 1;

    return 0;
}

/*
 * Experimental ordered BAM partition reader.
 *
 * This is deliberately internal and opt-in.  It uses index-derived virtual
 * offsets to split a coordinate-sorted BAM into non-overlapping file spans,
 * parses those spans on worker file handles, and drains completed spans in
 * file-offset order so sam_read1() callers still see one ordered stream.
 */
#define BAM_ORDERED_READER_MAGIC 0x626f7264u
#define BAM_ORDERED_READER_DEFAULT_THREADS 4
#define BAM_ORDERED_READER_TARGET_JOBS_PER_THREAD 16
#define BAM_ORDERED_READER_MAX_AHEAD_PER_THREAD 2

typedef struct bam_ordered_job_t {
    uint64_t beg;
    uint64_t end;
} bam_ordered_job_t;

typedef struct bam_ordered_batch_t {
    bam1_t *records;
    int n_records;
    int m_records;
} bam_ordered_batch_t;

typedef struct bam_ordered_result_t {
    int ready;
    int ret;
    bam_ordered_batch_t batch;
} bam_ordered_result_t;

typedef struct bam_ordered_reader_t {
    uint32_t magic;
    char *fn;
    bam_ordered_job_t *jobs;
    int n_jobs;
    int next_job;
    int drain_job;
    int max_inflight;

    pthread_t *threads;
    int n_threads;
    int n_threads_started;
    int active_workers;
    int closing;
    int error;

    bam_ordered_result_t *results;
    bam_ordered_batch_t curr;
    int curr_job;
    int curr_idx;

    pthread_mutex_t mutex;
    pthread_cond_t result_c;
    pthread_cond_t space_c;
} bam_ordered_reader_t;

static int bam_ordered_env_enabled(void)
{
    const char *env = getenv("HTS_BAM_ORDERED_READER");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_ordered_env_strict(void)
{
    const char *env = getenv("HTS_BAM_ORDERED_READER_REQUIRE");
    return env && *env && strcmp(env, "0") != 0;
}

static int bam_ordered_env_threads(void)
{
    const char *env = getenv("HTS_BAM_ORDERED_READER_THREADS");
    char *end = NULL;
    long n;

    if (!env || !*env)
        return BAM_ORDERED_READER_DEFAULT_THREADS;

    errno = 0;
    n = strtol(env, &end, 10);
    if (errno || end == env || *end || n < 1 || n > INT_MAX / 2)
        return BAM_ORDERED_READER_DEFAULT_THREADS;
    return (int)n;
}

static void bam_ordered_batch_init_record_slots(bam1_t *records, int beg,
                                                int end)
{
    int i;

    for (i = beg; i < end; i++)
        bam_set_mempolicy(&records[i], BAM_USER_OWNS_STRUCT);
}

static void bam_ordered_batch_destroy(bam_ordered_batch_t *batch)
{
    int i;

    if (!batch)
        return;
    if (batch->records) {
        for (i = 0; i < batch->m_records; i++)
            bam_destroy1(&batch->records[i]);
    }
    free(batch->records);
    memset(batch, 0, sizeof(*batch));
}

static int bam_ordered_move1(bam1_t *dst, bam1_t *src)
{
    uint32_t dst_policy = bam_get_mempolicy(dst);

    if ((dst_policy & BAM_USER_OWNS_DATA) == 0)
        free(dst->data);

    dst->core = src->core;
    dst->id = src->id;
    dst->data = src->data;
    dst->l_data = src->l_data;
    dst->m_data = src->m_data;
    bam_set_mempolicy(dst, dst_policy & BAM_USER_OWNS_STRUCT);

    src->data = NULL;
    src->l_data = 0;
    src->m_data = 0;
    return 0;
}

static bam1_t *bam_ordered_batch_next_record(bam_ordered_batch_t *batch)
{
    if (batch->n_records == batch->m_records) {
        int new_m;
        bam1_t *new_records;

        if (batch->m_records > INT_MAX / 2) {
            errno = ENOMEM;
            return NULL;
        }
        new_m = batch->m_records ? batch->m_records * 2 : 256;
        if ((size_t)new_m > SIZE_MAX / sizeof(*new_records)) {
            errno = ENOMEM;
            return NULL;
        }

        new_records = realloc(batch->records,
                              (size_t)new_m * sizeof(*new_records));
        if (!new_records) {
            errno = ENOMEM;
            return NULL;
        }
        memset(new_records + batch->m_records, 0,
               (size_t)(new_m - batch->m_records) * sizeof(*new_records));
        bam_ordered_batch_init_record_slots(new_records, batch->m_records,
                                            new_m);
        batch->records = new_records;
        batch->m_records = new_m;
    }

    return &batch->records[batch->n_records];
}

static int bam_ordered_uint64_cmp(const void *a, const void *b)
{
    uint64_t va = *(const uint64_t *)a;
    uint64_t vb = *(const uint64_t *)b;

    return va > vb ? 1 : va < vb ? -1 : 0;
}

static int bam_ordered_append_offset(uint64_t **offsets, int *n_offsets,
                                     int *m_offsets, uint64_t offset)
{
    uint64_t *new_offsets;
    int new_m;

    if (*n_offsets == *m_offsets) {
        new_m = *m_offsets ? (*m_offsets > INT_MAX / 2
                              ? INT_MAX : *m_offsets * 2) : 1024;
        if (new_m == *m_offsets ||
            (size_t)new_m > SIZE_MAX / sizeof(*new_offsets)) {
            errno = ENOMEM;
            return -1;
        }
        new_offsets = realloc(*offsets,
                              (size_t)new_m * sizeof(*new_offsets));
        if (!new_offsets) {
            errno = ENOMEM;
            return -1;
        }
        *offsets = new_offsets;
        *m_offsets = new_m;
    }

    (*offsets)[(*n_offsets)++] = offset;
    return 0;
}

static uint64_t bam_ordered_voff_block_span(uint64_t beg, uint64_t end)
{
    uint64_t beg_block = beg >> 16;
    uint64_t end_block = end >> 16;

    return end_block > beg_block ? end_block - beg_block : 1;
}

static int bam_ordered_append_job(bam_ordered_job_t **jobs, int *n_jobs,
                                  int *m_jobs, uint64_t beg, uint64_t end)
{
    bam_ordered_job_t *new_jobs;
    int new_m;

    if (end <= beg)
        return 0;

    if (*n_jobs == *m_jobs) {
        new_m = *m_jobs ? (*m_jobs > INT_MAX / 2 ? INT_MAX : *m_jobs * 2)
                        : 1024;
        if (new_m == *m_jobs ||
            (size_t)new_m > SIZE_MAX / sizeof(*new_jobs)) {
            errno = ENOMEM;
            return -1;
        }
        new_jobs = realloc(*jobs, (size_t)new_m * sizeof(*new_jobs));
        if (!new_jobs) {
            errno = ENOMEM;
            return -1;
        }
        *jobs = new_jobs;
        *m_jobs = new_m;
    }

    (*jobs)[*n_jobs].beg = beg;
    (*jobs)[*n_jobs].end = end;
    (*n_jobs)++;
    return 0;
}

static int bam_ordered_build_offsets(sam_hdr_t *hdr, const hts_idx_t *idx,
                                     uint64_t **offsets, int *n_offsets)
{
    uint64_t *out = NULL;
    int n = 0, m = 0;
    int tid, i, j;

    *offsets = NULL;
    *n_offsets = 0;

    for (tid = 0; tid < sam_hdr_nref(hdr); tid++) {
        uint64_t *ref_offsets = NULL;
        int n_ref_offsets = 0;

        if (!hts_idx_ref_has_data(idx, tid))
            continue;
        if (hts_idx_get_ref_file_offsets(idx, tid, &ref_offsets,
                                         &n_ref_offsets) < 0)
            goto fail;
        for (i = 0; i < n_ref_offsets; i++) {
            if (bam_ordered_append_offset(&out, &n, &m,
                                          ref_offsets[i]) < 0) {
                free(ref_offsets);
                goto fail;
            }
        }
        free(ref_offsets);
    }

    if (n < 2) {
        free(out);
        return 1;
    }

    qsort(out, (size_t)n, sizeof(*out), bam_ordered_uint64_cmp);
    for (i = 1, j = 0; i < n; i++) {
        if (out[i] != out[j])
            out[++j] = out[i];
    }
    n = j + 1;
    if (n < 2) {
        free(out);
        return 1;
    }

    *offsets = out;
    *n_offsets = n;
    return 0;

fail:
    free(out);
    *offsets = NULL;
    *n_offsets = 0;
    return -1;
}

static int bam_ordered_build_jobs(sam_hdr_t *hdr, const hts_idx_t *idx,
                                  int n_threads,
                                  bam_ordered_job_t **jobs, int *n_jobs)
{
    uint64_t *offsets = NULL;
    uint64_t total_span = 0, span_per_job, span;
    int n_offsets = 0, m_jobs = 0;
    int64_t target_jobs;
    int i, ret;

    *jobs = NULL;
    *n_jobs = 0;

    ret = bam_ordered_build_offsets(hdr, idx, &offsets, &n_offsets);
    if (ret != 0)
        return ret;

    for (i = 1; i < n_offsets; i++) {
        span = bam_ordered_voff_block_span(offsets[i - 1], offsets[i]);
        if (UINT64_MAX - total_span < span)
            total_span = UINT64_MAX;
        else
            total_span += span;
    }

    target_jobs = (int64_t)n_threads * BAM_ORDERED_READER_TARGET_JOBS_PER_THREAD;
    if (target_jobs < 1)
        target_jobs = 1;
    span_per_job = total_span / (uint64_t)target_jobs +
        (total_span % (uint64_t)target_jobs != 0);
    if (span_per_job < 1)
        span_per_job = 1;

    for (i = 0; i + 1 < n_offsets;) {
        int end_i = i + 1;
        uint64_t start = offsets[i];

        while (end_i + 1 < n_offsets &&
               bam_ordered_voff_block_span(start, offsets[end_i]) <
                   span_per_job)
            end_i++;

        if (bam_ordered_append_job(jobs, n_jobs, &m_jobs, start,
                                   offsets[end_i]) < 0)
            goto fail;
        i = end_i;
    }

    if (bam_ordered_append_job(jobs, n_jobs, &m_jobs,
                               offsets[n_offsets - 1], UINT64_MAX) < 0)
        goto fail;

    free(offsets);
    return *n_jobs > 0 ? 0 : 1;

fail:
    free(offsets);
    free(*jobs);
    *jobs = NULL;
    *n_jobs = 0;
    return -1;
}

static int bam_ordered_claim_job(bam_ordered_reader_t *reader)
{
    int job = -1;

    pthread_mutex_lock(&reader->mutex);
    while (!reader->closing && !reader->error &&
           reader->next_job < reader->n_jobs &&
           reader->next_job >= reader->drain_job + reader->max_inflight)
        pthread_cond_wait(&reader->space_c, &reader->mutex);

    if (!reader->closing && !reader->error && reader->next_job < reader->n_jobs)
        job = reader->next_job++;
    pthread_mutex_unlock(&reader->mutex);

    return job;
}

static void bam_ordered_set_error(bam_ordered_reader_t *reader)
{
    pthread_mutex_lock(&reader->mutex);
    reader->error = 1;
    pthread_cond_broadcast(&reader->result_c);
    pthread_cond_broadcast(&reader->space_c);
    pthread_mutex_unlock(&reader->mutex);
}

static int bam_ordered_is_closing(bam_ordered_reader_t *reader)
{
    int closing;

    pthread_mutex_lock(&reader->mutex);
    closing = reader->closing;
    pthread_mutex_unlock(&reader->mutex);
    return closing;
}

static int bam_ordered_complete_job(bam_ordered_reader_t *reader, int job_id,
                                    bam_ordered_batch_t *batch, int ret)
{
    pthread_mutex_lock(&reader->mutex);
    if (reader->closing) {
        pthread_mutex_unlock(&reader->mutex);
        bam_ordered_batch_destroy(batch);
        return -1;
    }

    reader->results[job_id].batch = *batch;
    reader->results[job_id].ret = ret;
    reader->results[job_id].ready = 1;
    memset(batch, 0, sizeof(*batch));
    pthread_cond_broadcast(&reader->result_c);
    pthread_mutex_unlock(&reader->mutex);
    return 0;
}

static void bam_ordered_worker_done(bam_ordered_reader_t *reader)
{
    pthread_mutex_lock(&reader->mutex);
    reader->active_workers--;
    pthread_cond_broadcast(&reader->result_c);
    pthread_mutex_unlock(&reader->mutex);
}

static int bam_ordered_process_job(bam_ordered_reader_t *reader, htsFile *fp,
                                   sam_hdr_t *hdr,
                                   const bam_ordered_job_t *job,
                                   bam_ordered_batch_t *batch)
{
    int ret, records_since_close_check = 0;

    if (!fp->is_bgzf || !fp->fp.bgzf) {
        errno = EINVAL;
        return -2;
    }
    if (bgzf_seek(fp->fp.bgzf, (int64_t)job->beg, SEEK_SET) < 0)
        return -2;

    while (job->end == UINT64_MAX ||
           (uint64_t)bgzf_tell(fp->fp.bgzf) < job->end) {
        bam1_t *rec;

        if (records_since_close_check >= 1024) {
            records_since_close_check = 0;
            if (bam_ordered_is_closing(reader))
                break;
        }

        rec = bam_ordered_batch_next_record(batch);
        if (!rec)
            return -2;

        ret = bam_read1(fp->fp.bgzf, rec);
        if (ret < 0)
            return ret == -1 ? 0 : -2;

        if (hdr && (rec->core.tid >= hdr->n_targets || rec->core.tid < -1 ||
                    rec->core.mtid >= hdr->n_targets || rec->core.mtid < -1)) {
            errno = ERANGE;
            return -3;
        }
        batch->n_records++;
        records_since_close_check++;
    }

    return 0;
}

static void *bam_ordered_worker_main(void *arg)
{
    bam_ordered_reader_t *reader = (bam_ordered_reader_t *)arg;
    htsFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    int job_id;

    fp = sam_open(reader->fn, "rb");
    if (!fp)
        goto fail;
    hdr = sam_hdr_read(fp);
    if (!hdr)
        goto fail;

    while ((job_id = bam_ordered_claim_job(reader)) >= 0) {
        bam_ordered_batch_t batch = {0};
        int ret = bam_ordered_process_job(reader, fp, hdr,
                                          &reader->jobs[job_id], &batch);
        if (ret < 0) {
            bam_ordered_batch_destroy(&batch);
            goto fail;
        }
        if (bam_ordered_complete_job(reader, job_id, &batch, ret) < 0)
            goto done;
    }

done:
    sam_hdr_destroy(hdr);
    if (fp)
        sam_close(fp);
    bam_ordered_worker_done(reader);
    return NULL;

fail:
    sam_hdr_destroy(hdr);
    if (fp)
        sam_close(fp);
    bam_ordered_set_error(reader);
    bam_ordered_worker_done(reader);
    return NULL;
}

static void bam_ordered_reader_destroy(bam_ordered_reader_t *reader)
{
    int i;

    if (!reader)
        return;

    pthread_mutex_lock(&reader->mutex);
    reader->closing = 1;
    pthread_cond_broadcast(&reader->result_c);
    pthread_cond_broadcast(&reader->space_c);
    pthread_mutex_unlock(&reader->mutex);

    for (i = 0; i < reader->n_threads_started; i++)
        pthread_join(reader->threads[i], NULL);

    for (i = 0; i < reader->n_jobs; i++)
        bam_ordered_batch_destroy(&reader->results[i].batch);
    bam_ordered_batch_destroy(&reader->curr);

    pthread_cond_destroy(&reader->space_c);
    pthread_cond_destroy(&reader->result_c);
    pthread_mutex_destroy(&reader->mutex);
    free(reader->results);
    free(reader->threads);
    free(reader->jobs);
    free(reader->fn);
    free(reader);
}

int sam_bam_state_destroy(htsFile *fp)
{
    uint32_t magic;

    if (!fp || !fp->state)
        return 0;

    magic = *(uint32_t *)fp->state;
    if (magic == BAM_ORDERED_READER_MAGIC) {
        bam_ordered_reader_t *reader = (bam_ordered_reader_t *)fp->state;

        fp->state = NULL;
        bam_ordered_reader_destroy(reader);
        return 0;
    }
    if (magic == BAM_STREAM_READER_MAGIC) {
        bam_stream_reader_t *reader = (bam_stream_reader_t *)fp->state;

        fp->state = NULL;
        bam_stream_reader_destroy(reader);
        return 0;
    }
    if (magic == BAM_BATCH_REQUEST_MAGIC) {
        bam_batch_request_t *req = (bam_batch_request_t *)fp->state;

        fp->state = NULL;
        free(req);
        return 0;
    }
    if (magic == BAM_DEFERRED_THREADS_MAGIC) {
        bam_deferred_threads_t *cfg = (bam_deferred_threads_t *)fp->state;

        fp->state = NULL;
        bam_deferred_threads_destroy(cfg);
        return 0;
    }
    return 0;
}

static bam_ordered_reader_t *bam_ordered_reader_open(htsFile *fp,
                                                     sam_hdr_t *hdr,
                                                     int n_threads)
{
    bam_ordered_reader_t *reader = NULL;
    hts_idx_t *idx = NULL;
    char *fnidx = NULL;
    int64_t cur;
    int i;
    int mutex_init = 0, result_c_init = 0, space_c_init = 0;

    if (!fp || !hdr || !fp->fn || hisremote(fp->fn) || !fp->is_bgzf ||
        !fp->fp.bgzf || fp->fp.bgzf->is_gzip)
        return NULL;

    if (n_threads <= 0 || n_threads > INT_MAX / 2)
        n_threads = bam_ordered_env_threads();
    if (!hts_idx_check_local(fp->fn, HTS_FMT_BAI, &fnidx))
        return NULL;
    idx = sam_index_load2(fp, fp->fn, fnidx);
    free(fnidx);
    fnidx = NULL;
    if (!idx)
        return NULL;

    reader = calloc(1, sizeof(*reader));
    if (!reader)
        goto fail;
    reader->magic = BAM_ORDERED_READER_MAGIC;
    reader->curr_job = -1;
    reader->n_threads = n_threads;
    reader->max_inflight = n_threads * BAM_ORDERED_READER_MAX_AHEAD_PER_THREAD;
    if (reader->max_inflight < 1)
        reader->max_inflight = 1;
    reader->fn = strdup(fp->fn);
    if (!reader->fn)
        goto fail;

    if (bam_ordered_build_jobs(hdr, idx, n_threads, &reader->jobs,
                               &reader->n_jobs) != 0)
        goto fail;
    hts_idx_destroy(idx);
    idx = NULL;

    if (reader->n_jobs < n_threads) {
        errno = 0;
        goto fail;
    }

    cur = bgzf_tell(fp->fp.bgzf);
    if (cur < 0 || (uint64_t)cur > reader->jobs[0].beg) {
        errno = 0;
        goto fail;
    }

    reader->results = calloc((size_t)reader->n_jobs,
                             sizeof(*reader->results));
    reader->threads = calloc((size_t)n_threads, sizeof(*reader->threads));
    if (!reader->results || !reader->threads)
        goto fail;

    if (pthread_mutex_init(&reader->mutex, NULL) != 0)
        goto fail;
    mutex_init = 1;
    if (pthread_cond_init(&reader->result_c, NULL) != 0)
        goto fail;
    result_c_init = 1;
    if (pthread_cond_init(&reader->space_c, NULL) != 0)
        goto fail;
    space_c_init = 1;

    for (i = 0; i < n_threads; i++) {
        pthread_mutex_lock(&reader->mutex);
        reader->active_workers++;
        pthread_mutex_unlock(&reader->mutex);
        if (pthread_create(&reader->threads[i], NULL,
                           bam_ordered_worker_main, reader) != 0) {
            pthread_mutex_lock(&reader->mutex);
            reader->active_workers--;
            pthread_mutex_unlock(&reader->mutex);
            goto fail;
        }
        reader->n_threads_started++;
    }

    return reader;

fail:
    free(fnidx);
    hts_idx_destroy(idx);
    if (reader) {
        if (mutex_init) {
            pthread_mutex_lock(&reader->mutex);
            reader->closing = 1;
            reader->error = 1;
            if (result_c_init)
                pthread_cond_broadcast(&reader->result_c);
            if (space_c_init)
                pthread_cond_broadcast(&reader->space_c);
            pthread_mutex_unlock(&reader->mutex);
        }
        for (i = 0; i < reader->n_threads_started; i++)
            pthread_join(reader->threads[i], NULL);
        if (space_c_init)
            pthread_cond_destroy(&reader->space_c);
        if (result_c_init)
            pthread_cond_destroy(&reader->result_c);
        if (mutex_init)
            pthread_mutex_destroy(&reader->mutex);
        if (reader->results) {
            int j;
            for (j = 0; j < reader->n_jobs; j++)
                bam_ordered_batch_destroy(&reader->results[j].batch);
        }
        free(reader->results);
        free(reader->threads);
        free(reader->jobs);
        free(reader->fn);
        free(reader);
    }
    return NULL;
}

static int bam_ordered_reader_next(bam_ordered_reader_t *reader, bam1_t *b)
{
    for (;;) {
        bam_ordered_result_t *res;

        if (reader->curr_job >= 0 &&
            reader->curr_idx < reader->curr.n_records) {
            bam_ordered_move1(b, &reader->curr.records[reader->curr_idx++]);
            return 0;
        }

        if (reader->curr_job >= 0) {
            bam_ordered_batch_destroy(&reader->curr);
            reader->curr_job = -1;
            reader->curr_idx = 0;
            pthread_mutex_lock(&reader->mutex);
            reader->drain_job++;
            pthread_cond_broadcast(&reader->space_c);
            pthread_mutex_unlock(&reader->mutex);
        }

        pthread_mutex_lock(&reader->mutex);
        while (!reader->error && reader->drain_job < reader->n_jobs &&
               !reader->results[reader->drain_job].ready)
            pthread_cond_wait(&reader->result_c, &reader->mutex);

        if (reader->error) {
            pthread_mutex_unlock(&reader->mutex);
            return -2;
        }
        if (reader->drain_job >= reader->n_jobs) {
            pthread_mutex_unlock(&reader->mutex);
            return -1;
        }

        res = &reader->results[reader->drain_job];
        reader->curr = res->batch;
        reader->curr_job = reader->drain_job;
        reader->curr_idx = 0;
        memset(&res->batch, 0, sizeof(res->batch));
        res->ready = 0;
        pthread_mutex_unlock(&reader->mutex);

        if (res->ret < 0)
            return res->ret;
    }
}

#define UMI_TAGS 5
typedef struct {
    kstring_t name;
    kstring_t comment; // NB: pointer into name, do not free
    kstring_t seq;
    kstring_t qual;
    int casava;
    int aux;
    int rnum;
    char BC[3];         // aux tag ID for barcode
    char UMI[UMI_TAGS][3]; // aux tag list for UMIs.
    khash_t(tag) *tags; // which aux tags to use (if empty, use all).
    char nprefix;
    int sra_names;
    regex_t regex;
} fastq_state;

// Initialise fastq state.
// Name char of '@' or '>' distinguishes fastq vs fasta variant
static fastq_state *fastq_state_init(int name_char) {
    fastq_state *x = (fastq_state *)calloc(1, sizeof(*x));
    if (!x)
        return NULL;
    strcpy(x->BC, "BC");
    x->nprefix = name_char;
    // Default Illumina naming convention
    char *re = "^[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:([^:#/]+)";
    if (regcomp(&x->regex, re, REG_EXTENDED) != 0) {
        free(x);
        return NULL;
    }

    return x;
}

void fastq_state_destroy(htsFile *fp) {
    if (fp->state) {
        fastq_state *x = (fastq_state *)fp->state;
        if (x->tags)
            kh_destroy(tag, x->tags);
        ks_free(&x->name);
        ks_free(&x->seq);
        ks_free(&x->qual);
        regfree(&x->regex);
        free(fp->state);
    }
}

int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...) {
    va_list args;

    if (!fp)
        return -1;
    if (!fp->state)
        if (!(fp->state = fastq_state_init(fp->format.format == fastq_format
                                           ? '@' : '>')))
            return -1;

    fastq_state *x = (fastq_state *)fp->state;

    switch (opt) {
    case FASTQ_OPT_CASAVA:
        x->casava = 1;
        break;

    case FASTQ_OPT_NAME2:
        x->sra_names = 1;
        break;

    case FASTQ_OPT_AUX: {
        va_start(args, opt);
        x->aux = 1;
        char *tag = va_arg(args, char *);
        va_end(args);
        if (tag && strcmp(tag, "1") != 0) {
            if (!x->tags)
                if (!(x->tags = kh_init(tag)))
                    return -1;

            size_t i, tlen = strlen(tag);
            for (i = 0; i+3 <= tlen+1; i += 3) {
                if (tag[i+0] == ',' || tag[i+1] == ',' ||
                    !(tag[i+2] == ',' || tag[i+2] == '\0')) {
                    hts_log_warning("Bad tag format '%.3s'; skipping option", tag+i);
                    break;
                }
                int ret, tcode = tag[i+0]*256 + tag[i+1];
                kh_put(tag, x->tags, tcode, &ret);
                if (ret < 0)
                    return -1;
            }
        }
        break;
    }

    case FASTQ_OPT_BARCODE: {
        va_start(args, opt);
        char *bc = va_arg(args, char *);
        va_end(args);
        strncpy(x->BC, bc, 2);
        x->BC[2] = 0;
        break;
    }

    case FASTQ_OPT_UMI: {
        // UMI tag: an empty string disables UMI by setting x->UMI[0] to \0\0\0
        va_start(args, opt);
        char *bc = va_arg(args, char *), *bc_orig = bc;
        va_end(args);
        if (!bc || strcmp(bc, "1") == 0)
            bc = "RX";
        int ntags = 0, err = 0;
        for (ntags = 0; *bc && ntags < UMI_TAGS; ntags++) {
            if (!isalpha(bc[0]) || !isalnum_c(bc[1])) {
                err = 1;
                break;
            }

            strncpy(x->UMI[ntags], bc, 3);
            bc += 2;
            if (*bc && *bc != ',') {
                err = 1;
                break;
            }
            bc+=(*bc==',');
            x->UMI[ntags][2] = 0;
        }
        for (; ntags < UMI_TAGS; ntags++)
            x->UMI[ntags][0] = x->UMI[ntags][1] = x->UMI[ntags][2] = 0;


        if (err)
            hts_log_warning("Bad UMI tag list '%s'", bc_orig);

        break;
    }

    case FASTQ_OPT_UMI_REGEX: {
        va_start(args, opt);
        char *re = va_arg(args, char *);
        va_end(args);

        regfree(&x->regex);
        if (regcomp(&x->regex, re, REG_EXTENDED) != 0) {
            hts_log_error("Regular expression '%s' is not supported", re);
            return -1;
        }
        break;
    }

    case FASTQ_OPT_RNUM:
        x->rnum = 1;
        break;

    default:
        break;
    }
    return 0;
}

static int fastq_parse1(htsFile *fp, bam1_t *b) {
    fastq_state *x = (fastq_state *)fp->state;
    size_t i, l;
    int ret = 0;

    if (fp->format.format == fasta_format && fp->line.s) {
        // For FASTA we've already read the >name line; steal it
        // Not the most efficient, but we don't optimise for fasta reading.
        if (fp->line.l == 0)
            return -1; // EOF

        free(x->name.s);
        x->name = fp->line;
        fp->line.l = fp->line.m = 0;
        fp->line.s = NULL;
    } else {
        // Read a FASTQ format entry.
        ret = hts_getline(fp, KS_SEP_LINE, &x->name);
        if (ret == -1)
            return -1;  // EOF
        else if (ret < -1)
            return ret; // ERR
    }

    // Name
    if (*x->name.s != x->nprefix)
        return -2;

    // Reverse the SRA strangeness of putting the run_name.number before
    // the read name.
    i = 0;
    char *name = x->name.s+1;
    if (x->sra_names) {
        char *cp = strpbrk(x->name.s, " \t");
        if (cp) {
            while (*cp == ' ' || *cp == '\t')
                cp++;
            *--cp = '@';
            i = cp - x->name.s;
            name = cp+1;
        }
    }

    l = x->name.l;
    char *s = x->name.s;
    while (i < l && !isspace_c(s[i]))
        i++;
    if (i < l) {
        s[i] = 0;
        x->name.l = i++;
    }

    // Comment; a kstring struct, but pointer into name line.  (Do not free)
    while (i < l && isspace_c(s[i]))
        i++;
    x->comment.s = s+i;
    x->comment.l = l - i;

    // Seq
    x->seq.l = 0;
    for (;;) {
        if ((ret = hts_getline(fp, KS_SEP_LINE, &fp->line)) < 0)
            if (fp->format.format == fastq_format || ret < -1)
                return -2;
        if (ret == -1 ||
            *fp->line.s == (fp->format.format == fastq_format ? '+' : '>'))
            break;
        if (kputsn(fp->line.s, fp->line.l, &x->seq) < 0)
            return -2;
    }

    // Qual
    if (fp->format.format == fastq_format) {
        size_t remainder = x->seq.l;
        x->qual.l = 0;
        do {
            if (hts_getline(fp, KS_SEP_LINE, &fp->line) < 0)
                return -2;
            if (fp->line.l > remainder)
                return -2;
            if (kputsn(fp->line.s, fp->line.l, &x->qual) < 0)
                return -2;
            remainder -= fp->line.l;
        } while (remainder > 0);

        // Decr qual
        for (i = 0; i < x->qual.l; i++)
            x->qual.s[i] -= '!';
    }

    int flag = BAM_FUNMAP; int pflag = BAM_FMUNMAP | BAM_FPAIRED;
    if (x->name.l > 2 &&
        x->name.s[x->name.l-2] == '/' &&
        isdigit_c(x->name.s[x->name.l-1])) {
        switch(x->name.s[x->name.l-1]) {
        case '1': flag |= BAM_FREAD1 | pflag; break;
        case '2': flag |= BAM_FREAD2 | pflag; break;
        default : flag |= BAM_FREAD1 | BAM_FREAD2 | pflag; break;
        }
        x->name.s[x->name.l-=2] = 0;
    }

    // Strip Illumina formatted UMI off read-name
    char UMI_seq[256]; // maximum length in spec
    size_t UMI_len = 0;
    if (x->UMI[0][0]) {
        regmatch_t match[3];
        if (regexec(&x->regex, x->name.s, 2, match, 0) == 0
            && match[0].rm_so >= 0     // whole regex
            && match[1].rm_so >= 0) {  // bracketted UMI component
            UMI_len = match[1].rm_eo - match[1].rm_so;
            if (UMI_len > 255) {
                hts_log_error("SAM read name is too long");
                return -2;
            }

            // The SAMTags spec recommends (but not requires) separating
            // barcodes with hyphen ('-').
            size_t i;
            for (i = 0; i < UMI_len; i++)
                UMI_seq[i] = isalpha_c(x->name.s[i+match[1].rm_so])
                    ? x->name.s[i+match[1].rm_so]
                    : '-';

            // Move any trailing #num earlier in the name
            if (UMI_len) {
                UMI_seq[UMI_len++] = 0;

                x->name.l = match[1].rm_so;
                if (x->name.l > 0 && x->name.s[x->name.l-1] == ':')
                    x->name.l--; // remove colon too
                char *cp = x->name.s + match[1].rm_eo;
                while (*cp)
                    x->name.s[x->name.l++] = *cp++;
                x->name.s[x->name.l] = 0;
            }
        }
    }

    // Convert to BAM
    ret = bam_set1(b,
                   x->name.s + x->name.l - name, name,
                   flag,
                   -1, -1, 0, // ref '*', pos, mapq,
                   0, NULL,     // no cigar,
                   -1, -1, 0,    // mate
                   x->seq.l, x->seq.s, x->qual.s,
                   0);
    if (ret < 0) return -2;

    // Add UMI tag if removed from read-name above
    if (UMI_len) {
        if (bam_aux_append(b, x->UMI[0], 'Z', UMI_len, (uint8_t *)UMI_seq) < 0)
            ret = -2;
    }

    // Identify Illumina CASAVA strings.
    // <read>:<is_filtered>:<control_bits>:<barcode_sequence>
    char *barcode = NULL;
    int barcode_len = 0;
    kstring_t *kc = &x->comment;
    char *endptr;
    if (x->casava &&
        // \d:[YN]:\d+:[ACGTN]+
        kc->l > 6 && (kc->s[1] | kc->s[3]) == ':' && isdigit_c(kc->s[0]) &&
        strtol(kc->s+4, &endptr, 10) >= 0 && endptr != kc->s+4
        && *endptr == ':') {

        // read num
        switch(kc->s[0]) {
        case '1': b->core.flag |= BAM_FREAD1 | pflag; break;
        case '2': b->core.flag |= BAM_FREAD2 | pflag; break;
        default : b->core.flag |= BAM_FREAD1 | BAM_FREAD2 | pflag; break;
        }

        if (kc->s[2] == 'Y')
            b->core.flag |= BAM_FQCFAIL;

        // Barcode, maybe numeric in which case we skip it
        if (!isdigit_c(endptr[1])) {
            barcode = endptr+1;
            for (i = barcode - kc->s; i < kc->l; i++)
                if (isspace_c(kc->s[i]))
                    break;

            kc->s[i] = 0;
            barcode_len = i+1-(barcode - kc->s);
        }
    }

    if (ret >= 0 && barcode_len)
        if (bam_aux_append(b, x->BC, 'Z', barcode_len, (uint8_t *)barcode) < 0)
            ret = -2;

    if (!x->aux)
        return ret;

    // Identify any SAM style aux tags in comments too.
    if (aux_parse(&kc->s[barcode_len], kc->s + kc->l, b, 1, x->tags) < 0)
        ret = -2;

    return ret;
}

static inline int sam_read1_bam_check_header(sam_hdr_t *h, bam1_t *b, int ret)
{
    if (h && ret >= 0) {
        if (b->core.tid  >= h->n_targets || b->core.tid  < -1 ||
            b->core.mtid >= h->n_targets || b->core.mtid < -1) {
            errno = ERANGE;
            return -3;
        }
    }
    return ret;
}

// Internal component of sam_read1 below
static inline int sam_read1_bam(htsFile *fp, sam_hdr_t *h, bam1_t *b) {
    if (fp->state) {
        uint32_t magic = *(uint32_t *)fp->state;
        if (magic == BAM_ORDERED_READER_MAGIC) {
            bam_ordered_reader_t *reader = (bam_ordered_reader_t *)fp->state;
            return sam_read1_bam_check_header(h, b,
                                              bam_ordered_reader_next(reader, b));
        }
        if (magic == BAM_STREAM_READER_MAGIC) {
            bam_stream_reader_t *reader = (bam_stream_reader_t *)fp->state;
            return sam_read1_bam_check_header(h, b,
                                              bam_stream_reader_next(reader, b));
        }
        if (magic == BAM_DEFERRED_THREADS_MAGIC) {
            bam_deferred_threads_t *cfg =
                (bam_deferred_threads_t *)fp->state;

            if (bam_stream_env_enabled()) {
                bam_stream_reader_t *reader = bam_stream_reader_open(fp, cfg, 0);
                if (reader) {
                    fp->state = NULL;
                    bam_deferred_threads_destroy(cfg);
                    fp->state = reader;
                    return sam_read1_bam_check_header(h, b,
                                                      bam_stream_reader_next(reader, b));
                }
                if (bam_stream_env_strict())
                    return -2;
            }

            if (bam_ordered_env_enabled()) {
                bam_ordered_reader_t *reader =
                    bam_ordered_reader_open(fp, h, cfg->n_threads);
                if (reader) {
                    fp->state = NULL;
                    bam_deferred_threads_destroy(cfg);
                    fp->state = reader;
                    return sam_read1_bam_check_header(h, b,
                                                      bam_ordered_reader_next(reader, b));
                }
                if (bam_ordered_env_strict())
                    return -2;
            }

            int ret = bam_deferred_threads_enable_bgzf(fp, cfg);
            fp->state = NULL;
            bam_deferred_threads_destroy(cfg);
            if (ret < 0)
                return -2;
        }
    } else if (bam_stream_env_enabled()) {
        bam_stream_reader_t *reader = bam_stream_reader_open(fp, NULL, 0);
        if (reader) {
            fp->state = reader;
            return sam_read1_bam_check_header(h, b,
                                              bam_stream_reader_next(reader, b));
        }
        if (bam_stream_env_strict())
            return -2;
    } else if (bam_ordered_env_enabled()) {
        bam_ordered_reader_t *reader = bam_ordered_reader_open(fp, h, 0);
        if (reader) {
            fp->state = reader;
            return sam_read1_bam_check_header(h, b,
                                              bam_ordered_reader_next(reader, b));
        }
        if (bam_ordered_env_strict())
            return -2;
    }

    return sam_read1_bam_check_header(h, b, bam_read1(fp->fp.bgzf, b));
}

// Internal component of sam_read1 below
static inline int sam_read1_cram(htsFile *fp, sam_hdr_t *h, bam1_t **b) {
    int ret = cram_get_bam_seq(fp->fp.cram, b);
    if (ret < 0)
        return cram_eof(fp->fp.cram) ? -1 : -2;

    if (bam_tag2cigar(*b, 1, 1) < 0)
        return -2;

    return ret;
}

// Internal component of sam_read1 below
static inline int sam_read1_sam(htsFile *fp, sam_hdr_t *h, bam1_t *b) {
    int ret;

    // Consume 1st line after header parsing as it wasn't using peek
    if (fp->line.l != 0) {
        ret = sam_parse1(&fp->line, h, b);
        fp->line.l = 0;
        return ret;
    }

    if (fp->state) {
        SAM_state *fd = (SAM_state *)fp->state;

        if (fp->format.compression == bgzf && fp->fp.bgzf->seeked) {
            // We don't support multi-threaded SAM parsing with seeks yet.
            int ret;
            if ((ret = sam_state_destroy(fp)) < 0) {
                errno = -ret;
                return -2;
            }
            if (bgzf_seek(fp->fp.bgzf, fp->fp.bgzf->seeked, SEEK_SET) < 0)
                return -2;
            fp->fp.bgzf->seeked = 0;
            goto err_recover;
        }

        if (!fd->h) {
            fd->h = h;
            fd->h->ref_count++;
            // Ensure hrecs is initialised now as we don't want multiple
            // threads trying to do this simultaneously.
            if (!fd->h->hrecs && sam_hdr_fill_hrecs(fd->h) < 0)
                return -2;

            // We can only do this once we've got a header
            if (pthread_create(&fd->dispatcher, NULL, sam_dispatcher_read,
                               fp) != 0)
                return -2;
            fd->dispatcher_set = 1;
        }

        if (fd->h != h) {
            hts_log_error("SAM multi-threaded decoding does not support changing header");
            return -2;
        }

        sp_bams *gb = fd->curr_bam;
        if (!gb) {
            if (fd->errcode) {
                // In case reader failed
                errno = fd->errcode;
                return -2;
            }

            pthread_mutex_lock(&fd->command_m);
            int cmd = fd->command;
            pthread_mutex_unlock(&fd->command_m);
            if (cmd == SAM_AT_EOF)
                return -1;

            hts_tpool_result *r = hts_tpool_next_result_wait(fd->q);
            if (!r)
                return -2;
            fd->curr_bam = gb = (sp_bams *)hts_tpool_result_data(r);
            hts_tpool_delete_result(r, 0);
        }
        if (!gb) {
            pthread_mutex_lock(&fd->command_m);
            fd->command = SAM_AT_EOF;
            pthread_mutex_unlock(&fd->command_m);
            return fd->errcode ? -2 : -1;
        }
        bam1_t *b_array = (bam1_t *)gb->bams;
        if (fd->curr_idx < gb->nbams)
            if (!bam_copy1(b, &b_array[fd->curr_idx++]))
                return -2;
        if (fd->curr_idx == gb->nbams) {
            pthread_mutex_lock(&fd->lines_m);
            gb->next = fd->bams;
            fd->bams = gb;
            pthread_mutex_unlock(&fd->lines_m);

            fd->curr_bam = NULL;
            fd->curr_idx = 0;
        // Consider prefetching next record?  I.e.
        // } else {
        //     __builtin_prefetch(&b_array[fd->curr_idx], 0, 3);
        }

        ret = 0;

    } else  {
    err_recover:
        ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
        if (ret < 0) return ret;

        ret = sam_parse1(&fp->line, h, b);
        fp->line.l = 0;
        if (ret < 0) {
            hts_log_warning("Parse error at line %lld", (long long)fp->lineno);
            if (h && h->ignore_sam_err) goto err_recover;
        }
    }

    return ret;
}

// Returns 0 on success,
//        -1 on EOF,
//       <-1 on error
int sam_read1(htsFile *fp, sam_hdr_t *h, bam1_t *b)
{
    int ret, pass_filter;

    do {
        switch (fp->format.format) {
        case bam:
            ret = sam_read1_bam(fp, h, b);
            break;

        case cram:
            ret = sam_read1_cram(fp, h, &b);
            break;

        case sam:
            ret = sam_read1_sam(fp, h, b);
            break;

        case fasta_format:
        case fastq_format: {
            fastq_state *x = (fastq_state *)fp->state;
            if (!x) {
                if (!(fp->state = fastq_state_init(fp->format.format
                                                   == fastq_format ? '@' : '>')))
                    return -2;
            }

            return fastq_parse1(fp, b);
        }

        case empty_format:
            errno = EPIPE;
            return -3;

        default:
            errno = EFTYPE;
            return -3;
        }

        pass_filter = (ret >= 0 && fp->filter)
            ? sam_passes_filter(h, b, fp->filter)
            : 1;
    } while (pass_filter == 0);

    return pass_filter < 0 ? -2 : ret;
}

// With gcc, -O3 or -ftree-loop-vectorize is really key here as otherwise
// this code isn't vectorised and runs far slower than is necessary (even
// with the restrict keyword being used).
static inline void HTS_OPT3
add33(uint8_t *a, const uint8_t * b, int32_t len) {
    uint32_t i;
    for (i = 0; i < len; i++)
        a[i] = b[i]+33;
}

static const uint8_t *sam_format_aux1_append(const uint8_t *key,
                                             const uint8_t type,
                                             const uint8_t *tag,
                                             const uint8_t *end,
                                             kstring_t *ks)
{
    if (type == 'Z' || type == 'H') {
        const uint8_t *nul = memchr(tag, '\0', end - tag);
        size_t len;
        int r = 0;

        if (!nul)
            return NULL;
        len = (size_t)(nul - tag);
        if (ks_resize(ks, ks->l + 5 + len + 1) < 0)
            return NULL;
        r |= kputsn_((const char *)key, 2, ks);
        r |= kputc_(':', ks);
        r |= kputc_(type, ks);
        r |= kputc_(':', ks);
        r |= kputsn_((const char *)tag, len, ks);
        r |= kputsn("", 0, ks);
        return r < 0 ? NULL : nul + 1;
    }
    return sam_format_aux1(key, type, tag, end, ks);
}

static inline uint32_t sam_record_view_cigar_word(
        const sam_bam_record_view_t *view, int idx)
{
    const uint8_t *cigar = sam_bam_record_view_cigar(view);

    if (view->batch_record)
        return le_to_u32(cigar + ((size_t)idx << 2));
    return ((const uint32_t *)cigar)[idx];
}

static int sam_format_record_view_append(const bam_hdr_t *h,
                                         const sam_bam_record_view_t *view,
                                         kstring_t *str)
{
    int i, r = 0;
    const uint8_t *s, *end;
    const bam1_core_t *c = view->core;
    const uint8_t *qname = sam_bam_record_view_qname(view);
    int qname_len;

    if (c->l_qname == 0)
        return -1;
    qname_len = c->l_qname - 1 - (view->batch_record ? 0 : c->l_extranul);
    if (qname_len < 0)
        return -1;
    r |= kputsn_((const char *)qname, qname_len, str);
    r |= kputc_('\t', str); // query name
    r |= kputw(c->flag, str); r |= kputc_('\t', str); // flag
    if (c->tid >= 0) { // chr
        r |= kputs(h->target_name[c->tid] , str);
        r |= kputc_('\t', str);
    } else r |= kputsn_("*\t", 2, str);
    r |= kputll(c->pos + 1, str); r |= kputc_('\t', str); // pos
    r |= kputw(c->qual, str); r |= kputc_('\t', str); // qual
    if (c->n_cigar) { // cigar
        for (i = 0; i < c->n_cigar; ++i) {
            uint32_t cigar = sam_record_view_cigar_word(view, i);

            r |= kputw(bam_cigar_oplen(cigar), str);
            r |= kputc_(bam_cigar_opchr(cigar), str);
        }
    } else r |= kputc_('*', str);
    r |= kputc_('\t', str);
    if (c->mtid < 0) r |= kputsn_("*\t", 2, str); // mate chr
    else if (c->mtid == c->tid) r |= kputsn_("=\t", 2, str);
    else {
        r |= kputs(h->target_name[c->mtid], str);
        r |= kputc_('\t', str);
    }
    r |= kputll(c->mpos + 1, str); r |= kputc_('\t', str); // mate pos
    r |= kputll(c->isize, str); r |= kputc_('\t', str); // template len
    if (c->l_qseq) { // seq and qual
        const uint8_t *seq = sam_bam_record_view_seq(view);
        const uint8_t *qual;
        if (ks_resize(str, str->l+2+2*c->l_qseq) < 0) goto mem_err;
        char *cp = str->s + str->l;

        // Sequence, 2 bases at a time
        nibble2base((uint8_t *)seq, cp, c->l_qseq);
        cp[c->l_qseq] = '\t';
        cp += c->l_qseq+1;

        // Quality
        qual = sam_bam_record_view_qual(view);
        i = 0;
        if (qual[0] == 0xff) {
            cp[i++] = '*';
        } else {
            add33((uint8_t *)cp, qual, c->l_qseq); // cp[i] = s[i]+33;
            i = c->l_qseq;
        }
        cp[i] = 0;
        cp += i;
        str->l = cp - str->s;
    } else r |= kputsn_("*\t*", 3, str);

    s = sam_bam_record_view_aux(view); // aux
    end = view->body + view->l_data;

    while (end - s >= 4) {
        r |= kputc_('\t', str);
        if ((s = sam_format_aux1_append(s, s[2], s+3, end, str)) == NULL)
            goto bad_aux;
    }
    r |= kputsn("", 0, str); // nul terminate
    if (r < 0) goto mem_err;

    return str->l;

 bad_aux:
    hts_log_error("Corrupted aux data for read %.*s flag %d",
                  c->l_qname, qname, c->flag);
    errno = EINVAL;
    return -1;

 mem_err:
    hts_log_error("Out of memory");
    errno = ENOMEM;
    return -1;
}

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    sam_bam_record_view_t view;

    sam_bam_record_view_from_bam(&view, b);
    return sam_format_record_view_append(h, &view, str);
}

static int sam_format_batch_record_append(const bam_hdr_t *h,
                                          const bam_batch_record_t *record,
                                          kstring_t *str)
{
    sam_bam_record_view_t view;

    if (sam_bam_batch_record_validate_sam_cigar(record) < 0)
        return -1;
    sam_bam_record_view_from_batch(&view, record);
    return sam_format_record_view_append(h, &view, str);
}

int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    str->l = 0;
    return sam_format1_append(h, b, str);
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end);
int fastq_format1(fastq_state *x, const bam1_t *b, kstring_t *str)
{
    unsigned flag = b->core.flag;
    int i, e = 0, len = b->core.l_qseq;
    uint8_t *seq, *qual;

    str->l = 0;

    // Name
    if (kputc(x->nprefix, str) == EOF || kputs(bam_get_qname(b), str) == EOF)
        return -1;

    // UMI tag
    if (x && *x->UMI[0]) {
        // Temporary copy of '#num' if present
        char plex[256];
        size_t len = str->l;
        while (len && str->s[len] != ':' && str->s[len] != '#')
            len--;

        if (str->s[len] == '#' && str->l - len < 255) {
            memcpy(plex, &str->s[len], str->l - len);
            plex[str->l - len] = 0;
            str->l = len;
        } else {
            *plex = 0;
        }

        uint8_t *bc = NULL;
        int n;
        for (n = 0; !bc && n < UMI_TAGS; n++)
            bc = bam_aux_get(b, x->UMI[n]);
        if (bc && *bc == 'Z') {
            int err = kputc(':', str) < 0;
            // Replace any non-alpha with '+'
            while (*++bc)
                err |= kputc(isalpha_c(*bc) ? toupper_c(*bc) : '+', str) < 0;
            if (err)
                return -1;
        }

        if (*plex && kputs(plex, str) < 0)
            return -1;
    }

    // /1 or /2 suffix
    if (x && x->rnum && (flag & BAM_FPAIRED)) {
        int r12 = flag & (BAM_FREAD1 | BAM_FREAD2);
        if (r12 == BAM_FREAD1) {
            if (kputs("/1", str) == EOF)
                return -1;
        } else if (r12 == BAM_FREAD2) {
            if (kputs("/2", str) == EOF)
                return -1;
        }
    }

    // Illumina CASAVA tag.
    // This is <rnum>:<Y/N qcfail>:<control-bits>:<barcode-or-zero>
    if (x && x->casava) {
        int rnum = (flag & BAM_FREAD1)? 1 : (flag & BAM_FREAD2)? 2 : 0;
        char filtered = (flag & BAM_FQCFAIL)? 'Y' : 'N';
        uint8_t *bc = bam_aux_get(b, x->BC);
        if (ksprintf(str, " %d:%c:0:%s", rnum, filtered,
                     bc ? (char *)bc+1 : "0") < 0)
            return -1;

        if (bc && (*bc != 'Z' || (!isupper_c(bc[1]) && !islower_c(bc[1])))) {
            hts_log_warning("BC tag starts with non-sequence base; using '0'");
            str->l -= strlen((char *)bc)-2; // limit to 1 char
            str->s[str->l-1] = '0';
            str->s[str->l] = 0;
            bc = NULL;
        }

        // Replace any non-alpha with '+'.  Ie seq-seq to seq+seq
        if (bc) {
            int l = strlen((char *)bc+1);
            char *c = (char *)str->s + str->l - l;
            for (i = 0; i < l; i++) {
                if (!isalpha_c(c[i]))
                    c[i] = '+';
                else if (islower_c(c[i]))
                    c[i] = toupper_c(c[i]);
            }
        }
    }

    // Aux tags
    if (x && x->aux) {
        uint8_t *s = bam_get_aux(b), *end = b->data + b->l_data;
        while (s && end - s >= 4) {
            int tt = s[0]*256 + s[1];
            if (x->tags == NULL ||
                kh_get(tag, x->tags, tt) != kh_end(x->tags)) {
                e |= kputc_('\t', str) < 0;
                if (!(s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)))
                    return -1;
            } else {
                s = skip_aux(s+2, end);
            }
        }
        e |= kputsn("", 0, str) < 0; // nul terminate
    }

    if (ks_resize(str, str->l + 1 + len+1 + 2 + len+1 + 1) < 0) return -1;
    e |= kputc_('\n', str) < 0;

    // Seq line
    seq = bam_get_seq(b);
    if (flag & BAM_FREVERSE)
        for (i = len-1; i >= 0; i--)
            e |= kputc_("!TGKCYSBAWRDMHVN"[bam_seqi(seq, i)], str) < 0;
    else
        for (i = 0; i < len; i++)
            e |= kputc_(seq_nt16_str[bam_seqi(seq, i)], str) < 0;


    // Qual line
    if (x->nprefix == '@') {
        kputsn("\n+\n", 3, str);
        qual = bam_get_qual(b);
        if (qual[0] == 0xff)
            for (i = 0; i < len; i++)
                e |= kputc_('B', str) < 0;
        else if (flag & BAM_FREVERSE)
            for (i = len-1; i >= 0; i--)
                e |= kputc_(33 + qual[i], str) < 0;
        else
            for (i = 0; i < len; i++)
                e |= kputc_(33 + qual[i], str) < 0;

    }
    e |= kputc('\n', str) < 0;

    return e ? -1 : str->l;
}

// Sadly we need to be able to modify the bam_hdr here so we can
// reference count the structure.
int sam_write1(htsFile *fp, const sam_hdr_t *h, const bam1_t *b)
{
    switch (fp->format.format) {
    case binary_format:
        fp->format.category = sequence_data;
        fp->format.format = bam;
        /* fall-through */
    case bam:
        return bam_write_idx1(fp, h, b);

    case cram:
        return cram_put_bam_seq(fp->fp.cram, (bam1_t *)b);

    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam:
        if (fp->state) {
            SAM_state *fd = (SAM_state *)fp->state;

            // Threaded output
            if (sam_start_threaded_output(fp, h) < 0)
                return -2;

            // Find a suitable BAM array to copy to
            sp_bams *gb = fd->curr_bam;
            if (!gb) {
                pthread_mutex_lock(&fd->lines_m);
                if (fd->bams) {
                    fd->curr_bam = gb = fd->bams;
                    fd->bams = gb->next;
                    gb->next = NULL;
                    gb->nbams = 0;
                    gb->bam_mem = 0;
                    pthread_mutex_unlock(&fd->lines_m);
                } else {
                    pthread_mutex_unlock(&fd->lines_m);
                    if (!(gb = calloc(1, sizeof(*gb)))) return -1;
                    if (!(gb->bams = calloc(SAM_NBAM, sizeof(*gb->bams)))) {
                        free(gb);
                        return -1;
                    }
                    gb->nbams = 0;
                    gb->abams = SAM_NBAM;
                    gb->bam_mem = 0;
                    gb->fd = fd;
                    fd->curr_idx = 0;
                    fd->curr_bam = gb;
                }
            }

            if (!bam_copy1(&gb->bams[gb->nbams++], b))
                return -2;
            gb->bam_mem += b->l_data + sizeof(*b);

            // Dispatch if full
            if (gb->nbams == SAM_NBAM || gb->bam_mem > SAM_NBYTES*0.8) {
                gb->serial = fd->serial++;
                pthread_mutex_lock(&fd->command_m);
                if (fd->errcode != 0) {
                    pthread_mutex_unlock(&fd->command_m);
                    return -fd->errcode;
                }
                if (hts_tpool_dispatch3(fd->p, fd->q, sam_format_worker, gb,
                                        cleanup_sp_bams,
                                        cleanup_sp_lines, 0) < 0) {
                    pthread_mutex_unlock(&fd->command_m);
                    return -1;
                }
                pthread_mutex_unlock(&fd->command_m);
                fd->curr_bam = NULL;
            }

            // Dummy value as we don't know how long it really is.
            // We could track file sizes via a SAM_state field, but I don't think
            // it is necessary.
            return 1;
        } else {
            if (sam_format1(h, b, &fp->line) < 0) return -1;
            kputc('\n', &fp->line);
            if (fp->is_bgzf) {
                if (bgzf_flush_try(fp->fp.bgzf, fp->line.l) < 0)
                    return -1;
                if ( bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l) != fp->line.l ) return -1;
            } else {
                if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;
            }

            if (fp->idx) {
                if (fp->format.compression == bgzf) {
                    if (bgzf_idx_push(fp->fp.bgzf, fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                      bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
                        return -1;
                    }
                } else {
                    if (hts_idx_push(fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                     bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
                        return -1;
                    }
                }
            }

            return fp->line.l;
        }


    case fasta_format:
    case fastq_format: {
        fastq_state *x = (fastq_state *)fp->state;
        if (!x) {
            if (!(fp->state = fastq_state_init(fp->format.format
                                               == fastq_format ? '@' : '>')))
                return -2;
        }

        if (fastq_format1(fp->state, b, &fp->line) < 0)
            return -1;
        if (fp->is_bgzf) {
            if (bgzf_flush_try(fp->fp.bgzf, fp->line.l) < 0)
                return -1;
            if (bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l) != fp->line.l)
                return -1;
        } else {
            if (hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l)
                return -1;
        }
        return fp->line.l;
    }

    default:
        errno = EBADF;
        return -1;
    }
}

/************************
 *** Auxiliary fields ***
 ************************/
#ifndef HTS_LITTLE_ENDIAN
static int aux_to_le(char type, uint8_t *out, const uint8_t *in, size_t len) {
    int tsz = aux_type2size(type);

    if (tsz >= 2 && tsz <= 8 && (len & (tsz - 1)) != 0) return -1;

    switch (tsz) {
        case 'H': case 'Z': case 1:  // Trivial
            memcpy(out, in, len);
            break;

#define aux_val_to_le(type_t, store_le) do {                            \
        type_t v;                                                       \
        size_t i;                                                       \
        for (i = 0; i < len; i += sizeof(type_t), out += sizeof(type_t)) { \
            memcpy(&v, in + i, sizeof(type_t));                         \
            store_le(v, out);                                           \
        }                                                               \
    } while (0)

        case 2: aux_val_to_le(uint16_t, u16_to_le); break;
        case 4: aux_val_to_le(uint32_t, u32_to_le); break;
        case 8: aux_val_to_le(uint64_t, u64_to_le); break;

#undef aux_val_to_le

        case 'B': { // Recurse!
            uint32_t n;
            if (len < 5) return -1;
            memcpy(&n, in + 1, 4);
            out[0] = in[0];
            u32_to_le(n, out + 1);
            return aux_to_le(in[0], out + 5, in + 5, len - 5);
        }

        default: // Unknown type code
            return -1;
    }



    return 0;
}
#endif

int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data)
{
    uint32_t new_len;

    assert(b->l_data >= 0);
    new_len = b->l_data + 3 + len;
    if (new_len > INT32_MAX || new_len < b->l_data) goto nomem;

    if (realloc_bam_data(b, new_len) < 0) return -1;

    b->data[b->l_data] = tag[0];
    b->data[b->l_data + 1] = tag[1];
    b->data[b->l_data + 2] = type;

#ifdef HTS_LITTLE_ENDIAN
    memcpy(b->data + b->l_data + 3, data, len);
#else
    if (aux_to_le(type, b->data + b->l_data + 3, data, len) != 0) {
        errno = EINVAL;
        return -1;
    }
#endif

    b->l_data = new_len;

    return 0;

 nomem:
    errno = ENOMEM;
    return -1;
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        s = memchr(s, 0, end-s);
        return s ? s+1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

uint8_t *bam_aux_first(const bam1_t *b)
{
    uint8_t *s = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    if (end - s <= 2) { errno = ENOENT; return NULL; }
    return s+2;
}

uint8_t *bam_aux_next(const bam1_t *b, const uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = s? skip_aux((uint8_t *) s, end) : end;
    if (next == NULL) goto bad_aux;
    if (end - next <= 2) { errno = ENOENT; return NULL; }
    return next+2;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s flag %d",
                  bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return NULL;
}

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    for (s = bam_aux_first(b); s; s = bam_aux_next(b, s))
        if (s[-2] == tag[0] && s[-1] == tag[1]) {
            // Check the tag value is valid and complete
            uint8_t *e = skip_aux(s, b->data + b->l_data);
            if (e == NULL) goto bad_aux;
            if ((*s == 'Z' || *s == 'H') && *(e - 1) != '\0') goto bad_aux;

            return s;
        }

    // errno now as set by bam_aux_first()/bam_aux_next()
    return NULL;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s flag %d",
                  bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return NULL;
}

int bam_aux_del(bam1_t *b, uint8_t *s)
{
    s = bam_aux_remove(b, s);
    return (s || errno == ENOENT)? 0 : -1;
}

uint8_t *bam_aux_remove(bam1_t *b, uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = skip_aux(s, end);
    if (next == NULL) goto bad_aux;

    b->l_data -= next - (s-2);
    if (next >= end) { errno = ENOENT; return NULL; }

    memmove(s-2, next, end - next);
    return s;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s flag %d",
                  bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return NULL;
}

int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data)
{
    // FIXME: This is not at all efficient!
    size_t ln = len >= 0 ? len : strlen(data) + 1;
    size_t old_ln = 0;
    int need_nul = ln == 0 || data[ln - 1] != '\0';
    int save_errno = errno;
    int new_tag = 0;
    uint8_t *s = bam_aux_get(b,tag), *e;

    if (s) {  // Replacing existing tag
        char type = *s;
        if (type != 'Z') {
            hts_log_error("Called bam_aux_update_str for type '%c' instead of 'Z'", type);
            errno = EINVAL;
            return -1;
        }
        s++;
        e = memchr(s, '\0', b->data + b->l_data - s);
        old_ln = (e ? e - s : b->data + b->l_data - s) + 1;
        s -= 3;
    } else {
        if (errno != ENOENT) { // Invalid aux data, give up
            return -1;
        } else { // Tag doesn't exist - put it on the end
            errno = save_errno;
            s = b->data + b->l_data;
            new_tag = 3;
        }
    }

    if (old_ln < ln + need_nul + new_tag) {
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, ln + need_nul + new_tag - old_ln) < 0)
            return -1;
        s = b->data + s_offset;
    }
    if (!new_tag) {
        memmove(s + 3 + ln + need_nul,
                s + 3 + old_ln,
                b->l_data - (s + 3 - b->data) - old_ln);
    }
    b->l_data += new_tag + ln + need_nul - old_ln;

    s[0] = tag[0];
    s[1] = tag[1];
    s[2] = 'Z';
    memmove(s+3,data,ln);
    if (need_nul) s[3 + ln] = '\0';
    return 0;
}

int bam_aux_update_int(bam1_t *b, const char tag[2], int64_t val)
{
    uint32_t sz, old_sz = 0, new = 0;
    uint8_t *s, type;

    if (val < INT32_MIN || val > UINT32_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    if (val < INT16_MIN)       { type = 'i'; sz = 4; }
    else if (val < INT8_MIN)   { type = 's'; sz = 2; }
    else if (val < 0)          { type = 'c'; sz = 1; }
    else if (val < UINT8_MAX)  { type = 'C'; sz = 1; }
    else if (val < UINT16_MAX) { type = 'S'; sz = 2; }
    else                       { type = 'I'; sz = 4; }

    s = bam_aux_get(b, tag);
    if (s) {  // Tag present - how big was the old one?
        switch (*s) {
            case 'c': case 'C': old_sz = 1; break;
            case 's': case 'S': old_sz = 2; break;
            case 'i': case 'I': old_sz = 4; break;
            default: errno = EINVAL; return -1;  // Not an integer
        }
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            s = b->data + b->l_data;
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    if (new || old_sz < sz) {
        // Make room for new tag
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, (new ? 3 : 0) + sz - old_sz) < 0)
            return -1;
        s =  b->data + s_offset;
        if (new) { // Add tag id
            *s++ = tag[0];
            *s++ = tag[1];
        } else {   // Shift following data so we have space
            memmove(s + sz, s + old_sz, b->l_data - s_offset - old_sz);
        }
    } else {
        // Reuse old space.  Data value may be bigger than necessary but
        // we avoid having to move everything else
        sz = old_sz;
        type = (val < 0 ? "\0cs\0i" : "\0CS\0I")[old_sz];
        assert(type > 0);
    }
    *s++ = type;
#ifdef HTS_LITTLE_ENDIAN
    memcpy(s, &val, sz);
#else
    switch (sz) {
        case 4:  u32_to_le(val, s); break;
        case 2:  u16_to_le(val, s); break;
        default: *s = val; break;
    }
#endif
    b->l_data += (new ? 3 : 0) + sz - old_sz;
    return 0;
}

int bam_aux_update_float(bam1_t *b, const char tag[2], float val)
{
    uint8_t *s = bam_aux_get(b, tag);
    int shrink = 0, new = 0;

    if (s) { // Tag present - what was it?
        switch (*s) {
            case 'f': break;
            case 'd': shrink = 1; break;
            default: errno = EINVAL; return -1;  // Not a float
        }
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    if (new) { // Ensure there's room
        if (possibly_expand_bam_data(b, 3 + 4) < 0)
            return -1;
        s = b->data + b->l_data;
        *s++ = tag[0];
        *s++ = tag[1];
    } else if (shrink) { // Convert non-standard double tag to float
        memmove(s + 5, s + 9, b->l_data - ((s + 9) - b->data));
        b->l_data -= 4;
    }
    *s++ = 'f';
    float_to_le(val, s);
    if (new) b->l_data += 7;

    return 0;
}

int bam_aux_update_array(bam1_t *b, const char tag[2],
                         uint8_t type, uint32_t items, void *data)
{
    uint8_t *s = bam_aux_get(b, tag);
    size_t old_sz = 0, new_sz;
    int new = 0;

    if (s) { // Tag present
        if (*s != 'B') { errno = EINVAL; return -1; }
        old_sz = aux_type2size(s[1]);
        if (old_sz < 1 || old_sz > 4) { errno = EINVAL; return -1; }
        old_sz *= le_to_u32(s + 2);
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            s = b->data + b->l_data;
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    new_sz = aux_type2size(type);
    if (new_sz < 1 || new_sz > 4) { errno = EINVAL; return -1; }
    if (items > INT32_MAX / new_sz) { errno = ENOMEM; return -1; }
    new_sz *= items;

    if (new || old_sz < new_sz) {
        // Make room for new tag
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, (new ? 8 : 0) + new_sz - old_sz) < 0)
            return -1;
        s =  b->data + s_offset;
    }
    if (new) { // Add tag id and type
        *s++ = tag[0];
        *s++ = tag[1];
        *s = 'B';
        b->l_data += 8 + new_sz;
    } else if (old_sz != new_sz) { // shift following data if necessary
        memmove(s + 6 + new_sz, s + 6 + old_sz,
                b->l_data - ((s + 6 + old_sz) - b->data));
        b->l_data -= old_sz;
        b->l_data += new_sz;
    }

    s[1] = type;
    u32_to_le(items, s + 2);
    if (new_sz > 0) {
#ifdef HTS_LITTLE_ENDIAN
        memcpy(s + 6, data, new_sz);
#else
        return aux_to_le(type, s + 6, data, new_sz);
#endif
    }
    return 0;
}

static inline int64_t get_int_aux_val(uint8_t type, const uint8_t *s,
                                      uint32_t idx)
{
    switch (type) {
        case 'c': return le_to_i8(s + idx);
        case 'C': return s[idx];
        case 's': return le_to_i16(s + 2 * idx);
        case 'S': return le_to_u16(s + 2 * idx);
        case 'i': return le_to_i32(s + 4 * idx);
        case 'I': return le_to_u32(s + 4 * idx);
        default:
            errno = EINVAL;
            return 0;
    }
}

int64_t bam_aux2i(const uint8_t *s)
{
    int type;
    type = *s++;
    return get_int_aux_val(type, s, 0);
}

double bam_aux2f(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'd') return le_to_double(s);
    else if (type == 'f') return le_to_float(s);
    else return get_int_aux_val(type, s, 0);
}

char bam_aux2A(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'A') return *(char*)s;
    errno = EINVAL;
    return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'Z' || type == 'H') return (char*)s;
    errno = EINVAL;
    return 0;
}

uint32_t bam_auxB_len(const uint8_t *s)
{
    if (s[0] != 'B') {
        errno = EINVAL;
        return 0;
    }
    return le_to_u32(s + 2);
}

int64_t bam_auxB2i(const uint8_t *s, uint32_t idx)
{
    uint32_t len = bam_auxB_len(s);
    if (idx >= len) {
        errno = ERANGE;
        return 0;
    }
    return get_int_aux_val(s[1], s + 6, idx);
}

double bam_auxB2f(const uint8_t *s, uint32_t idx)
{
    uint32_t len = bam_auxB_len(s);
    if (idx >= len) {
        errno = ERANGE;
        return 0.0;
    }
    if (s[1] == 'f') return le_to_float(s + 6 + 4 * idx);
    else return get_int_aux_val(s[1], s + 6, idx);
}

int sam_open_mode(char *mode, const char *fn, const char *format)
{
    // TODO Parse "bam5" etc for compression level
    if (format == NULL) {
        // Try to pick a format based on the filename extension
        char extension[HTS_MAX_EXT_LEN];
        if (find_file_extension(fn, extension) < 0) return -1;
        return sam_open_mode(mode, fn, extension);
    }
    else if (strcasecmp(format, "bam") == 0) strcpy(mode, "b");
    else if (strcasecmp(format, "cram") == 0) strcpy(mode, "c");
    else if (strcasecmp(format, "sam") == 0) strcpy(mode, "");
    else if (strcasecmp(format, "sam.gz") == 0) strcpy(mode, "z");
    else if (strcasecmp(format, "fastq") == 0 ||
             strcasecmp(format, "fq") == 0) strcpy(mode, "f");
    else if (strcasecmp(format, "fastq.gz") == 0 ||
             strcasecmp(format, "fq.gz") == 0) strcpy(mode, "fz");
    else if (strcasecmp(format, "fasta") == 0 ||
             strcasecmp(format, "fa") == 0) strcpy(mode, "F");
    else if (strcasecmp(format, "fasta.gz") == 0 ||
             strcasecmp(format, "fa.gz") == 0) strcpy(mode, "Fz");
    else return -1;

    return 0;
}

// A version of sam_open_mode that can handle ,key=value options.
// The format string is allocated and returned, to be freed by the caller.
// Prefix should be "r" or "w",
char *sam_open_mode_opts(const char *fn,
                         const char *mode,
                         const char *format)
{
    char *mode_opts = malloc((format ? strlen(format) : 1) +
                             (mode   ? strlen(mode)   : 1) + 12);
    char *opts, *cp;
    int format_len;

    if (!mode_opts)
        return NULL;

    strcpy(mode_opts, mode ? mode : "r");
    cp = mode_opts + strlen(mode_opts);

    if (format == NULL) {
        // Try to pick a format based on the filename extension
        char extension[HTS_MAX_EXT_LEN];
        if (find_file_extension(fn, extension) < 0) {
            free(mode_opts);
            return NULL;
        }
        if (sam_open_mode(cp, fn, extension) == 0) {
            return mode_opts;
        } else {
            free(mode_opts);
            return NULL;
        }
    }

    if ((opts = strchr(format, ','))) {
        format_len = opts-format;
    } else {
        opts="";
        format_len = strlen(format);
    }

    if (strncmp(format, "bam", format_len) == 0) {
        *cp++ = 'b';
    } else if (strncmp(format, "cram", format_len) == 0) {
        *cp++ = 'c';
    } else if (strncmp(format, "cram2", format_len) == 0) {
        *cp++ = 'c';
        strcpy(cp, ",VERSION=2.1");
        cp += 12;
    } else if (strncmp(format, "cram3", format_len) == 0) {
        *cp++ = 'c';
        strcpy(cp, ",VERSION=3.0");
        cp += 12;
    } else if (strncmp(format, "sam", format_len) == 0) {
        ; // format mode=""
    } else if (strncmp(format, "sam.gz", format_len) == 0) {
        *cp++ = 'z';
    } else if (strncmp(format, "fastq", format_len) == 0 ||
               strncmp(format, "fq", format_len) == 0) {
        *cp++ = 'f';
    } else if (strncmp(format, "fastq.gz", format_len) == 0 ||
               strncmp(format, "fq.gz", format_len) == 0) {
        *cp++ = 'f';
        *cp++ = 'z';
    } else if (strncmp(format, "fasta", format_len) == 0 ||
               strncmp(format, "fa", format_len) == 0) {
        *cp++ = 'F';
    } else if (strncmp(format, "fasta.gz", format_len) == 0 ||
               strncmp(format, "fa", format_len) == 0) {
        *cp++ = 'F';
        *cp++ = 'z';
    } else {
        free(mode_opts);
        return NULL;
    }

    strcpy(cp, opts);

    return mode_opts;
}

#define STRNCMP(a,b,n) (strncasecmp((a),(b),(n)) || strlen(a)!=(n))
int bam_str2flag(const char *str)
{
    char *end, *beg = (char*) str;
    long int flag = strtol(str, &end, 0);
    if ( end!=str ) return flag;    // the conversion was successful
    flag = 0;
    while ( *str )
    {
        end = beg;
        while ( *end && *end!=',' ) end++;
        if ( !STRNCMP("PAIRED",beg,end-beg) ) flag |= BAM_FPAIRED;
        else if ( !STRNCMP("PROPER_PAIR",beg,end-beg) ) flag |= BAM_FPROPER_PAIR;
        else if ( !STRNCMP("UNMAP",beg,end-beg) ) flag |= BAM_FUNMAP;
        else if ( !STRNCMP("MUNMAP",beg,end-beg) ) flag |= BAM_FMUNMAP;
        else if ( !STRNCMP("REVERSE",beg,end-beg) ) flag |= BAM_FREVERSE;
        else if ( !STRNCMP("MREVERSE",beg,end-beg) ) flag |= BAM_FMREVERSE;
        else if ( !STRNCMP("READ1",beg,end-beg) ) flag |= BAM_FREAD1;
        else if ( !STRNCMP("READ2",beg,end-beg) ) flag |= BAM_FREAD2;
        else if ( !STRNCMP("SECONDARY",beg,end-beg) ) flag |= BAM_FSECONDARY;
        else if ( !STRNCMP("QCFAIL",beg,end-beg) ) flag |= BAM_FQCFAIL;
        else if ( !STRNCMP("DUP",beg,end-beg) ) flag |= BAM_FDUP;
        else if ( !STRNCMP("SUPPLEMENTARY",beg,end-beg) ) flag |= BAM_FSUPPLEMENTARY;
        else return -1;
        if ( !*end ) break;
        beg = end + 1;
    }
    return flag;
}

char *bam_flag2str(int flag)
{
    kstring_t str = {0,0,0};
    if ( flag&BAM_FPAIRED ) ksprintf(&str,"%s%s", str.l?",":"","PAIRED");
    if ( flag&BAM_FPROPER_PAIR ) ksprintf(&str,"%s%s", str.l?",":"","PROPER_PAIR");
    if ( flag&BAM_FUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","UNMAP");
    if ( flag&BAM_FMUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","MUNMAP");
    if ( flag&BAM_FREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","REVERSE");
    if ( flag&BAM_FMREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","MREVERSE");
    if ( flag&BAM_FREAD1 ) ksprintf(&str,"%s%s", str.l?",":"","READ1");
    if ( flag&BAM_FREAD2 ) ksprintf(&str,"%s%s", str.l?",":"","READ2");
    if ( flag&BAM_FSECONDARY ) ksprintf(&str,"%s%s", str.l?",":"","SECONDARY");
    if ( flag&BAM_FQCFAIL ) ksprintf(&str,"%s%s", str.l?",":"","QCFAIL");
    if ( flag&BAM_FDUP ) ksprintf(&str,"%s%s", str.l?",":"","DUP");
    if ( flag&BAM_FSUPPLEMENTARY ) ksprintf(&str,"%s%s", str.l?",":"","SUPPLEMENTARY");
    if ( str.l == 0 ) kputsn("", 0, &str);
    return str.s;
}


/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

#include <assert.h>

/*******************
 *** Memory pool ***
 *******************/

typedef struct {
    int k, y;
    hts_pos_t x, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
    bam1_t b;
    hts_pos_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
    bam_pileup_cd cd;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init(void)
{
    mempool_t *mp;
    mp = (mempool_t*)calloc(1, sizeof(mempool_t));
    return mp;
}
static void mp_destroy(mempool_t *mp)
{
    int k;
    for (k = 0; k < mp->n; ++k) {
        free(mp->buf[k]->b.data);
        free(mp->buf[k]);
    }
    free(mp->buf);
    free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
    else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
    --mp->cnt; p->next = 0; // clear lbnode_t::next here
    if (mp->n == mp->max) {
        mp->max = mp->max? mp->max<<1 : 256;
        mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
    }
    mp->buf[mp->n++] = p;
}

/**********************
 *** CIGAR resolver ***
 **********************/

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordinate of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, hts_pos_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

    bam1_t *b = p->b;
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    int k;
    // determine the current CIGAR operation
    //fprintf(stderr, "%s\tpos=%ld\tend=%ld\t(%d,%ld,%d)\n", bam_get_qname(b), pos, s->end, s->k, s->x, s->y);
    if (s->k == -1) { // never processed
        p->qpos = 0;
        if (c->n_cigar == 1) { // just one operation, save a loop
          if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
        } else { // find the first match or deletion
            for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
                int op = _cop(cigar[k]);
                int l = _cln(cigar[k]);
                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                    op == BAM_CEQUAL || op == BAM_CDIFF) break;
                else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
            }
            assert(k < c->n_cigar);
            s->k = k;
        }
    } else { // the read has been processed before
        int op, l = _cln(cigar[s->k]);
        if (pos - s->x >= l) { // jump to the next operation
            assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
            op = _cop(cigar[s->k+1]);
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                ++s->k;
            } else { // find the next M/D/N/=/X
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                for (k = s->k + 1; k < c->n_cigar; ++k) {
                    op = _cop(cigar[k]), l = _cln(cigar[k]);
                    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
                    else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
                }
                s->k = k;
            }
            assert(s->k < c->n_cigar); // otherwise a bug
        } // else, do nothing
    }
    { // collect pileup information
        int op, l;
        op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
        p->is_del = p->indel = p->is_refskip = 0;
        if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
            int op2 = _cop(cigar[s->k+1]);
            int l2 = _cln(cigar[s->k+1]);
            if (op2 == BAM_CDEL && op != BAM_CDEL) {
                // At start of a new deletion, merge e.g. 1D2D to 3D.
                // Within a deletion (the 2D in 1D2D) we keep p->indel=0
                // and rely on is_del=1 as we would for 3D.
                p->indel = -(int)l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CDEL) p->indel -= l2;
                    else break;
                }
            } else if (op2 == BAM_CINS) {
                p->indel = l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) p->indel += l2;
                    else if (op2 != BAM_CPAD) break;
                }
            } else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) {
                int l3 = 0;
                for (k = s->k + 2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) l3 += l2;
                    else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
                }
                if (l3 > 0) p->indel = l3;
            }
        }
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            p->qpos = s->y + (pos - s->x);
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
            p->is_refskip = (op == BAM_CREF_SKIP);
        } // cannot be other operations; otherwise a bug
        p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
    }
    p->cigar_ind = s->k;
    return 1;
}

/*******************************
 *** Expansion of insertions ***
 *******************************/

/*
 * Fills out the kstring with the padded insertion sequence for the current
 * location in 'p'.  If this is not an insertion site, the string is blank.
 *
 * This variant handles base modifications, but only when "m" is non-NULL.
 *
 * Returns the number of inserted base on success, with string length being
 *        accessable via ins->l;
 *        -1 on failure.
 */
int bam_plp_insertion_mod(const bam_pileup1_t *p,
                          hts_base_mod_state *m,
                          kstring_t *ins, int *del_len) {
    int j, k, indel, nb = 0;
    uint32_t *cigar;

    if (p->indel <= 0) {
        if (ks_resize(ins, 1) < 0)
            return -1;
        ins->l = 0;
        ins->s[0] = '\0';
        return 0;
    }

    if (del_len)
        *del_len = 0;

    // Measure indel length including pads
    indel = 0;
    k = p->cigar_ind+1;
    cigar = bam_get_cigar(p->b);
    while (k < p->b->core.n_cigar) {
        switch (cigar[k] & BAM_CIGAR_MASK) {
        case BAM_CPAD:
        case BAM_CINS:
            indel += (cigar[k] >> BAM_CIGAR_SHIFT);
            break;
        default:
            k = p->b->core.n_cigar;
            break;
        }
        k++;
    }
    nb = ins->l = indel;

    // Produce sequence
    if (ks_resize(ins, indel+1) < 0)
        return -1;
    indel = 0;
    k = p->cigar_ind+1;
    j = 1;
    while (k < p->b->core.n_cigar) {
        int l, c;
        switch (cigar[k] & BAM_CIGAR_MASK) {
        case BAM_CPAD:
            for (l = 0; l < (cigar[k]>>BAM_CIGAR_SHIFT); l++)
                ins->s[indel++] = '*';
            break;
        case BAM_CINS:
            for (l = 0; l < (cigar[k]>>BAM_CIGAR_SHIFT); l++, j++) {
                c = p->qpos + j - p->is_del < p->b->core.l_qseq
                    ? seq_nt16_str[bam_seqi(bam_get_seq(p->b),
                                            p->qpos + j - p->is_del)]
                    : 'N';
                ins->s[indel++] = c;
                int nm;
                hts_base_mod mod[256];
                if (m && (nm = bam_mods_at_qpos(p->b, p->qpos + j - p->is_del,
                                                m, mod, 256)) > 0) {
                    int o_indel = indel;
                    if (ks_resize(ins, ins->l + nm*16+3) < 0)
                        return -1;
                    ins->s[indel++] = '[';
                    int j;
                    for (j = 0; j < nm; j++) {
                        char qual[20];
                        if (mod[j].qual >= 0)
                            snprintf(qual, sizeof(qual), "%d", mod[j].qual);
                        else
                            *qual=0;
                        if (mod[j].modified_base < 0)
                            // ChEBI
                            indel += snprintf(&ins->s[indel], ins->m - indel,
                                              "%c(%d)%s",
                                              "+-"[mod[j].strand],
                                              -mod[j].modified_base,
                                              qual);
                        else
                            indel += snprintf(&ins->s[indel], ins->m - indel,
                                              "%c%c%s",
                                              "+-"[mod[j].strand],
                                              mod[j].modified_base,
                                              qual);
                    }
                    ins->s[indel++] = ']';
                    ins->l += indel - o_indel; // grow by amount we used
                }
            }
            break;
        case BAM_CDEL:
            // eg cigar 1M2I1D gives mpileup output in T+2AA-1C style
            if (del_len)
                *del_len = cigar[k]>>BAM_CIGAR_SHIFT;
            // fall through
        default:
            k = p->b->core.n_cigar;
            break;
        }
        k++;
    }
    ins->s[indel] = '\0';
    ins->l = indel; // string length

    return nb;      // base length
}

/*
 * Fills out the kstring with the padded insertion sequence for the current
 * location in 'p'.  If this is not an insertion site, the string is blank.
 *
 * This is the original interface with no capability for reporting base
 * modifications.
 *
 * Returns the length of insertion string on success;
 *        -1 on failure.
 */
int bam_plp_insertion(const bam_pileup1_t *p, kstring_t *ins, int *del_len) {
    return bam_plp_insertion_mod(p, NULL, ins, del_len);
}

/***********************
 *** Pileup iterator ***
 ***********************/

// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct bam_plp_s {
    mempool_t *mp;
    lbnode_t *head, *tail;
    int32_t tid, max_tid;
    hts_pos_t pos, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;

    // For notification of creation and destruction events
    // and associated client-owned pointer.
    int (*plp_construct)(void *data, const bam1_t *b, bam_pileup_cd *cd);
    int (*plp_destruct )(void *data, const bam1_t *b, bam_pileup_cd *cd);
};

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
{
    bam_plp_t iter;
    iter = (bam_plp_t)calloc(1, sizeof(struct bam_plp_s));
    iter->mp = mp_init();
    iter->head = iter->tail = mp_alloc(iter->mp);
    iter->max_tid = iter->max_pos = -1;
    iter->maxcnt = 8000;
    if (func) {
        iter->func = func;
        iter->data = data;
        iter->b = bam_init1();
    }
    return iter;
}

int bam_plp_init_overlaps(bam_plp_t iter)
{
    iter->overlaps = kh_init(olap_hash);  // hash for tweaking quality of bases in overlapping reads
    return iter->overlaps ? 0 : -1;
}

void bam_plp_destroy(bam_plp_t iter)
{
    lbnode_t *p, *pnext;
    if ( iter->overlaps ) kh_destroy(olap_hash, iter->overlaps);
    for (p = iter->head; p != NULL; p = pnext) {
        if (iter->plp_destruct && p != iter->tail)
            iter->plp_destruct(iter->data, &p->b, &p->cd);
        pnext = p->next;
        mp_free(iter->mp, p);
    }
    mp_destroy(iter->mp);
    if (iter->b) bam_destroy1(iter->b);
    free(iter->plp);
    free(iter);
}

void bam_plp_constructor(bam_plp_t plp,
                         int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd)) {
    plp->plp_construct = func;
}

void bam_plp_destructor(bam_plp_t plp,
                        int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd)) {
    plp->plp_destruct = func;
}

//---------------------------------
//---  Tweak overlapping reads
//---------------------------------

/**
 *  cigar_iref2iseq_set()  - find the first CMATCH setting the ref and the read index
 *  cigar_iref2iseq_next() - get the next CMATCH base
 *  @cigar:       pointer to current cigar block (rw)
 *  @cigar_max:   pointer just beyond the last cigar block
 *  @icig:        position within the current cigar block (rw)
 *  @iseq:        position in the sequence (rw)
 *  @iref:        position with respect to the beginning of the read (iref_pos - b->core.pos) (rw)
 *
 *  Returns BAM_CMATCH, -1 when there is no more cigar to process or the requested position is not covered,
 *  or -2 on error.
 */
static inline int cigar_iref2iseq_set(const uint32_t **cigar,
                                      const uint32_t *cigar_max,
                                      hts_pos_t *icig,
                                      hts_pos_t *iseq,
                                      hts_pos_t *iref)
{
    hts_pos_t pos = *iref;
    if ( pos < 0 ) return -1;
    *icig = 0;
    *iseq = 0;
    *iref = 0;
    while ( *cigar<cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            pos -= ncig;
            if ( pos < 0 ) { *icig = ncig + pos; *iseq += *icig; *iref += *icig; return BAM_CMATCH; }
            (*cigar)++; *iseq += ncig; *icig = 0; *iref += ncig;
            continue;
        }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP )
        {
            pos -= ncig;
            if ( pos<0 ) pos = 0;
            (*cigar)++; *icig = 0; *iref += ncig;
            continue;
        }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    return -1;
}
static inline int cigar_iref2iseq_next(const uint32_t **cigar,
                                       const uint32_t *cigar_max,
                                       hts_pos_t *icig,
                                       hts_pos_t *iseq,
                                       hts_pos_t *iref)
{
    while ( *cigar < cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            if ( *icig >= ncig - 1 ) { *icig = -1;  (*cigar)++; continue; }
            (*iseq)++; (*icig)++; (*iref)++;
            return BAM_CMATCH;
        }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP ) { (*cigar)++; (*iref) += ncig; *icig = -1; continue; }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = -1; continue; }
        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = -1; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = -1; continue; }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    *iref = -1;
    return -1;
}

// Given overlapping read 'a' (left) and 'b' (right) on the same
// template, adjust quality values to zero for either a or b.
// Note versions 1.12 and earlier always removed quality from 'b' for
// matching bases.  Now we select a or b semi-randomly based on name hash.
// Returns 0 on success,
//        -1 on failure
static int tweak_overlap_quality(bam1_t *a, bam1_t *b)
{
    const uint32_t *a_cigar = bam_get_cigar(a),
        *a_cigar_max = a_cigar + a->core.n_cigar;
    const uint32_t *b_cigar = bam_get_cigar(b),
        *b_cigar_max = b_cigar + b->core.n_cigar;
    hts_pos_t a_icig = 0, a_iseq = 0;
    hts_pos_t b_icig = 0, b_iseq = 0;
    uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
    uint8_t *a_seq  = bam_get_seq(a), *b_seq = bam_get_seq(b);

    hts_pos_t iref   = b->core.pos;
    hts_pos_t a_iref = iref - a->core.pos;
    hts_pos_t b_iref = iref - b->core.pos;

    int a_ret = cigar_iref2iseq_set(&a_cigar, a_cigar_max,
                                    &a_icig, &a_iseq, &a_iref);
    if ( a_ret<0 )
        // no overlap or error
        return a_ret<-1 ? -1:0;

    int b_ret = cigar_iref2iseq_set(&b_cigar, b_cigar_max,
                                    &b_icig, &b_iseq, &b_iref);
    if ( b_ret<0 )
        // no overlap or error
        return b_ret<-1 ? -1:0;

    // Determine which seq is the one getting modified qualities.
    uint8_t amul, bmul;
    if (__ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(a))) & 1) {
        amul = 1;
        bmul = 0;
    } else {
        amul = 0;
        bmul = 1;
    }

    // Loop over the overlapping region nulling qualities in either
    // seq a or b.
    int err = 0;
    while ( 1 ) {
        // Step to next matching reference position in a and b
        while ( a_ret >= 0 && a_iref>=0 && a_iref < iref - a->core.pos )
            a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max,
                                         &a_icig, &a_iseq, &a_iref);
        if ( a_ret<0 ) { // done
            err = a_ret<-1?-1:0;
            break;
        }

        while ( b_ret >= 0 && b_iref>=0 && b_iref < iref - b->core.pos )
            b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig,
                                         &b_iseq, &b_iref);
        if ( b_ret<0 ) { // done
            err = b_ret<-1?-1:0;
            break;
        }

        if ( iref < a_iref + a->core.pos )
            iref = a_iref + a->core.pos;

        if ( iref < b_iref + b->core.pos )
            iref = b_iref + b->core.pos;

        iref++;

        // If A or B has a deletion then we catch up the other to this point.
        // We also amend quality values using the same rules for mismatch.
        if (a_iref+a->core.pos != b_iref+b->core.pos) {
            if (a_iref+a->core.pos < b_iref+b->core.pos
                && b_cigar > bam_get_cigar(b)
                && bam_cigar_op(b_cigar[-1]) == BAM_CDEL) {
                // Del in B means it's moved on further than A
                do {
                    a_qual[a_iseq] = amul
                        ? a_qual[a_iseq]*0.8
                        : 0;
                    a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max,
                                                 &a_icig, &a_iseq, &a_iref);
                    if (a_ret < 0)
                        return -(a_ret<-1); // 0 or -1
                } while (a_iref + a->core.pos < b_iref+b->core.pos);
            } else if (a_cigar > bam_get_cigar(a)
                       && bam_cigar_op(a_cigar[-1]) == BAM_CDEL) {
                // Del in A means it's moved on further than B
                do {
                    b_qual[b_iseq] = bmul
                        ? b_qual[b_iseq]*0.8
                        : 0;
                    b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max,
                                                 &b_icig, &b_iseq, &b_iref);
                    if (b_ret < 0)
                        return -(b_ret<-1); // 0 or -1
                } while (b_iref + b->core.pos < a_iref+a->core.pos);
            } else {
                // Anything else, eg ref-skip, we don't support here
                continue;
            }
        }

        // fprintf(stderr, "a_cig=%ld,%ld b_cig=%ld,%ld iref=%ld "
        //         "a_iref=%ld b_iref=%ld a_iseq=%ld b_iseq=%ld\n",
        //         a_cigar-bam_get_cigar(a), a_icig,
        //         b_cigar-bam_get_cigar(b), b_icig,
        //         iref, a_iref+a->core.pos+1, b_iref+b->core.pos+1,
        //         a_iseq, b_iseq);

        if (a_iseq > a->core.l_qseq || b_iseq > b->core.l_qseq)
            // Fell off end of sequence, bad CIGAR?
            return -1;

        // We're finally at the same ref base in both a and b.
        // Check if the bases match (confident) or mismatch
        // (not so confident).
        if ( bam_seqi(a_seq,a_iseq) == bam_seqi(b_seq,b_iseq) ) {
            // We are very confident about this base.  Use sum of quals
            int qual = a_qual[a_iseq] + b_qual[b_iseq];
            a_qual[a_iseq] = amul * (qual>200 ? 200 : qual);
            b_qual[b_iseq] = bmul * (qual>200 ? 200 : qual);;
        } else {
            // Not so confident about anymore given the mismatch.
            // Reduce qual for lowest quality base.
            if ( a_qual[a_iseq] > b_qual[b_iseq] ) {
                // A highest qual base; keep
                a_qual[a_iseq] = 0.8 * a_qual[a_iseq];
                b_qual[b_iseq] = 0;
            } else if (a_qual[a_iseq] < b_qual[b_iseq] ) {
                // B highest qual base; keep
                b_qual[b_iseq] = 0.8 * b_qual[b_iseq];
                a_qual[a_iseq] = 0;
            } else {
                // Both equal, so pick randomly
                a_qual[a_iseq] = amul * 0.8 * a_qual[a_iseq];
                b_qual[b_iseq] = bmul * 0.8 * b_qual[b_iseq];
            }
        }
    }

    return err;
}

// Fix overlapping reads. Simple soft-clipping did not give good results.
// Lowering qualities of unwanted bases is more selective and works better.
//
// Returns 0 on success, -1 on failure
static int overlap_push(bam_plp_t iter, lbnode_t *node)
{
    if ( !iter->overlaps ) return 0;

    // mapped mates and paired reads only
    if ( node->b.core.flag&BAM_FMUNMAP || !(node->b.core.flag&BAM_FPROPER_PAIR) ) return 0;

    // no overlap possible, unless some wild cigar
    if ( (node->b.core.mtid >= 0 && node->b.core.tid != node->b.core.mtid)
         || (llabs(node->b.core.isize) >= 2*node->b.core.l_qseq
         && node->b.core.mpos >= node->end) // for those wild cigars
       ) return 0;

    khiter_t kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(&node->b));
    if ( kitr==kh_end(iter->overlaps) )
    {
        // Only add reads where the mate is still to arrive
        if (node->b.core.mpos >= node->b.core.pos ||
            ((node->b.core.flag & BAM_FPAIRED) && node->b.core.mpos == -1)) {
            int ret;
            kitr = kh_put(olap_hash, iter->overlaps, bam_get_qname(&node->b), &ret);
            if (ret < 0) return -1;
            kh_value(iter->overlaps, kitr) = node;
        }
    }
    else
    {
        lbnode_t *a = kh_value(iter->overlaps, kitr);
        int err = tweak_overlap_quality(&a->b, &node->b);
        kh_del(olap_hash, iter->overlaps, kitr);
        assert(a->end-1 == a->s.end);
        return err;
    }
    return 0;
}

static void overlap_remove(bam_plp_t iter, const bam1_t *b)
{
    if ( !iter->overlaps ) return;

    khiter_t kitr;
    if ( b )
    {
        if ( b->core.flag&BAM_FUNMAP || !(b->core.flag&BAM_FPROPER_PAIR) ) //no need
            return;

        kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(b));
        if ( kitr!=kh_end(iter->overlaps) )
            kh_del(olap_hash, iter->overlaps, kitr);
    }
    else
    {
        // remove all
        for (kitr = kh_begin(iter->overlaps); kitr<kh_end(iter->overlaps); kitr++)
            if ( kh_exist(iter->overlaps, kitr) ) kh_del(olap_hash, iter->overlaps, kitr);
    }
}



// Prepares next pileup position in bam records collected by bam_plp_auto -> user func -> bam_plp_push. Returns
// pointer to the piled records if next position is ready or NULL if there is not enough records in the
// buffer yet (the current position is still the maximum position across all buffered reads).
const bam_pileup1_t *bam_plp64_next(bam_plp_t iter, int *_tid, hts_pos_t *_pos, int *_n_plp)
{
    if (iter->error) { *_n_plp = -1; return NULL; }
    *_n_plp = 0;
    if (iter->is_eof && iter->head == iter->tail) return NULL;
    while (iter->is_eof || iter->max_tid > iter->tid || (iter->max_tid == iter->tid && iter->max_pos > iter->pos)) {
        int n_plp = 0;
        // write iter->plp at iter->pos
        lbnode_t **pptr = &iter->head;
        while (*pptr != iter->tail) {
            if ((*pptr)->next)
                hts_prefetch((*pptr)->next);
            lbnode_t *p = *pptr;
            if (p->b.core.tid < iter->tid || (p->b.core.tid == iter->tid && p->end <= iter->pos)) { // then remove
                overlap_remove(iter, &p->b);
                if (iter->plp_destruct)
                    iter->plp_destruct(iter->data, &p->b, &p->cd);
                *pptr = p->next; mp_free(iter->mp, p);
            }
            else {
                if (p->b.core.tid == iter->tid && p->beg <= iter->pos) { // here: p->end > pos; then add to pileup
                    if (n_plp == iter->max_plp) { // then double the capacity
                        iter->max_plp = iter->max_plp? iter->max_plp<<1 : 256;
                        iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
                    }
                    iter->plp[n_plp].b = &p->b;
                    iter->plp[n_plp].cd = p->cd;
                    if (resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s)) ++n_plp; // actually always true...
                }
                pptr = &(*pptr)->next;
            }
        }
        *_n_plp = n_plp; *_tid = iter->tid; *_pos = iter->pos;
        // update iter->tid and iter->pos
        if (iter->head != iter->tail) {
            if (iter->tid > iter->head->b.core.tid) {
                hts_log_error("Unsorted input. Pileup aborts");
                iter->error = 1;
                *_n_plp = -1;
                return NULL;
            }
        }
        if (iter->tid < iter->head->b.core.tid) { // come to a new reference sequence
            iter->tid = iter->head->b.core.tid; iter->pos = iter->head->beg; // jump to the next reference
        } else if (iter->pos < iter->head->beg) { // here: tid == head->b.core.tid
            iter->pos = iter->head->beg; // jump to the next position
        } else ++iter->pos; // scan contiguously
        // return
        if (n_plp) return iter->plp;
        if (iter->is_eof && iter->head == iter->tail) break;
    }
    return NULL;
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
    hts_pos_t pos64 = 0;
    const bam_pileup1_t *p = bam_plp64_next(iter, _tid, &pos64, _n_plp);
    if (pos64 < INT_MAX) {
        *_pos = pos64;
    } else {
        hts_log_error("Position %"PRId64" too large", pos64);
        *_pos = INT_MAX;
        iter->error = 1;
        *_n_plp = -1;
        return NULL;
    }
    return p;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b)
{
    if (iter->error) return -1;
    if (b) {
        if (b->core.tid < 0) { overlap_remove(iter, b); return 0; }
        // Skip only unmapped reads here, any additional filtering must be done in iter->func
        if (b->core.flag & BAM_FUNMAP) { overlap_remove(iter, b); return 0; }
        if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt)
        {
            overlap_remove(iter, b);
            return 0;
        }
        if (bam_copy1(&iter->tail->b, b) == NULL)
            return -1;
        iter->tail->b.id = iter->id++;
        iter->tail->beg = b->core.pos;
        // Use raw rlen rather than bam_endpos() which adjusts rlen=0 to rlen=1
        iter->tail->end = b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        iter->tail->s = g_cstate_null; iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
        if (b->core.tid < iter->max_tid) {
            hts_log_error("The input is not sorted (chromosomes out of order)");
            iter->error = 1;
            return -1;
        }
        if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
            hts_log_error("The input is not sorted (reads out of order)");
            iter->error = 1;
            return -1;
        }
        iter->max_tid = b->core.tid; iter->max_pos = iter->tail->beg;
        if (iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
            lbnode_t *next = mp_alloc(iter->mp);
            if (!next) {
                iter->error = 1;
                return -1;
            }
            if (iter->plp_construct) {
                if (iter->plp_construct(iter->data, &iter->tail->b,
                                        &iter->tail->cd) < 0) {
                    mp_free(iter->mp, next);
                    iter->error = 1;
                    return -1;
                }
            }
            if (overlap_push(iter, iter->tail) < 0) {
                mp_free(iter->mp, next);
                iter->error = 1;
                return -1;
            }
            iter->tail->next = next;
            iter->tail = iter->tail->next;
        }
    } else iter->is_eof = 1;
    return 0;
}

const bam_pileup1_t *bam_plp64_auto(bam_plp_t iter, int *_tid, hts_pos_t *_pos, int *_n_plp)
{
    const bam_pileup1_t *plp;
    if (iter->func == 0 || iter->error) { *_n_plp = -1; return 0; }
    if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
    else { // no pileup line can be obtained; read alignments
        *_n_plp = 0;
        if (iter->is_eof) return 0;
        int ret;
        while ( (ret=iter->func(iter->data, iter->b)) >= 0) {
            if (bam_plp_push(iter, iter->b) < 0) {
                *_n_plp = -1;
                return 0;
            }
            if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
            // otherwise no pileup line can be returned; read the next alignment.
        }
        if ( ret < -1 ) { iter->error = ret; *_n_plp = -1; return 0; }
        if (bam_plp_push(iter, 0) < 0) {
            *_n_plp = -1;
            return 0;
        }
        if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
        return 0;
    }
}

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
    hts_pos_t pos64 = 0;
    const bam_pileup1_t *p = bam_plp64_auto(iter, _tid, &pos64, _n_plp);
    if (pos64 < INT_MAX) {
        *_pos = pos64;
    } else {
        hts_log_error("Position %"PRId64" too large", pos64);
        *_pos = INT_MAX;
        iter->error = 1;
        *_n_plp = -1;
        return NULL;
    }
    return p;
}

void bam_plp_reset(bam_plp_t iter)
{
    overlap_remove(iter, NULL);
    iter->max_tid = iter->max_pos = -1;
    iter->tid = iter->pos = 0;
    iter->is_eof = 0;
    while (iter->head != iter->tail) {
        lbnode_t *p = iter->head;
        iter->head = p->next;
        mp_free(iter->mp, p);
    }
}

void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
{
    iter->maxcnt = maxcnt;
}

/************************
 *** Mpileup iterator ***
 ************************/

struct bam_mplp_s {
    int n;
    int32_t min_tid, *tid;
    hts_pos_t min_pos, *pos;
    bam_plp_t *iter;
    int *n_plp;
    const bam_pileup1_t **plp;
};

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
{
    int i;
    bam_mplp_t iter;
    iter = (bam_mplp_t)calloc(1, sizeof(struct bam_mplp_s));
    iter->pos = (hts_pos_t*)calloc(n, sizeof(hts_pos_t));
    iter->tid = (int32_t*)calloc(n, sizeof(int32_t));
    iter->n_plp = (int*)calloc(n, sizeof(int));
    iter->plp = (const bam_pileup1_t**)calloc(n, sizeof(bam_pileup1_t*));
    iter->iter = (bam_plp_t*)calloc(n, sizeof(bam_plp_t));
    iter->n = n;
    iter->min_pos = HTS_POS_MAX;
    iter->min_tid = (uint32_t)-1;
    for (i = 0; i < n; ++i) {
        iter->iter[i] = bam_plp_init(func, data[i]);
        iter->pos[i] = iter->min_pos;
        iter->tid[i] = iter->min_tid;
    }
    return iter;
}

int bam_mplp_init_overlaps(bam_mplp_t iter)
{
    int i, r = 0;
    for (i = 0; i < iter->n; ++i)
        r |= bam_plp_init_overlaps(iter->iter[i]);
    return r == 0 ? 0 : -1;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
{
    int i;
    for (i = 0; i < iter->n; ++i)
        iter->iter[i]->maxcnt = maxcnt;
}

void bam_mplp_destroy(bam_mplp_t iter)
{
    int i;
    for (i = 0; i < iter->n; ++i) bam_plp_destroy(iter->iter[i]);
    free(iter->iter); free(iter->pos); free(iter->tid);
    free(iter->n_plp); free(iter->plp);
    free(iter);
}

int bam_mplp64_auto(bam_mplp_t iter, int *_tid, hts_pos_t *_pos, int *n_plp, const bam_pileup1_t **plp)
{
    int i, ret = 0;
    hts_pos_t new_min_pos = HTS_POS_MAX;
    uint32_t new_min_tid = (uint32_t)-1;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
            int tid;
            hts_pos_t pos;
            iter->plp[i] = bam_plp64_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
            if ( iter->iter[i]->error ) return -1;
            if (iter->plp[i]) {
                iter->tid[i] = tid;
                iter->pos[i] = pos;
            } else {
                iter->tid[i] = 0;
                iter->pos[i] = 0;
            }
        }
        if (iter->plp[i]) {
            if (iter->tid[i] < new_min_tid) {
                new_min_tid = iter->tid[i];
                new_min_pos = iter->pos[i];
            } else if (iter->tid[i] == new_min_tid && iter->pos[i] < new_min_pos) {
                new_min_pos = iter->pos[i];
            }
        }
    }
    iter->min_pos = new_min_pos;
    iter->min_tid = new_min_tid;
    if (new_min_pos == HTS_POS_MAX) return 0;
    *_tid = new_min_tid; *_pos = new_min_pos;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
            n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
            ++ret;
        } else n_plp[i] = 0, plp[i] = 0;
    }
    return ret;
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
    hts_pos_t pos64 = 0;
    int ret = bam_mplp64_auto(iter, _tid, &pos64, n_plp, plp);
    if (ret >= 0) {
        if (pos64 < INT_MAX) {
            *_pos = pos64;
        } else {
            hts_log_error("Position %"PRId64" too large", pos64);
            *_pos = INT_MAX;
            return -1;
        }
    }
    return ret;
}

void bam_mplp_reset(bam_mplp_t iter)
{
    int i;
    iter->min_pos = HTS_POS_MAX;
    iter->min_tid = (uint32_t)-1;
    for (i = 0; i < iter->n; ++i) {
        bam_plp_reset(iter->iter[i]);
        iter->pos[i] = HTS_POS_MAX;
        iter->tid[i] = (uint32_t)-1;
        iter->n_plp[i] = 0;
        iter->plp[i] = NULL;
    }
}

void bam_mplp_constructor(bam_mplp_t iter,
                          int (*func)(void *arg, const bam1_t *b, bam_pileup_cd *cd)) {
    int i;
    for (i = 0; i < iter->n; ++i)
        bam_plp_constructor(iter->iter[i], func);
}

void bam_mplp_destructor(bam_mplp_t iter,
                         int (*func)(void *arg, const bam1_t *b, bam_pileup_cd *cd)) {
    int i;
    for (i = 0; i < iter->n; ++i)
        bam_plp_destructor(iter->iter[i], func);
}

#endif // ~!defined(BAM_NO_PILEUP)
