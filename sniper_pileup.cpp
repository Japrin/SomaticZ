/*
 * sniper_pileup.cpp
 *
 * Code from SomaticSniper
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sam.h"
#include <assert.h>
#include "sniper_pileup.h"

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg, end;
	struct __linkbuf_t *next;
} lbnode_t;

/* --- BEGIN: Memory pool */

typedef struct {
	int cnt, n, max;
	lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init()
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

/* --- END: Memory pool */

/* --- BEGIN: Auxiliary functions */

static inline int resolve_cigar(bam_pileup1_t *p, uint32_t pos)
{
	unsigned k;
	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t x = c->pos, y = 0;
	int ret = 1, is_restart = 1;

	if (c->flag&BAM_FUNMAP) return 0; // unmapped read
	assert(x <= pos);
	p->qpos = -1; p->indel = 0; p->is_del = p->is_head = p->is_tail = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK; // operation
		int l = bam1_cigar(b)[k] >> BAM_CIGAR_SHIFT; // length
		if (op == BAM_CMATCH) { // NOTE: this assumes the first and the last operation MUST BE a match or a clip
			if (x + l > pos) { // overlap with pos
				p->indel = p->is_del = 0;
				p->qpos = y + (pos - x);
				if (x == pos && is_restart) p->is_head = 1;
				if (x + l - 1 == pos) { // come to the end of a match
					if (k < c->n_cigar - 1) { // there are additional operation(s)
						uint32_t cigar = bam1_cigar(b)[k+1]; // next CIGAR
						int op_next = cigar&BAM_CIGAR_MASK; // next CIGAR operation
						if (op_next == BAM_CDEL) p->indel = -(int32_t)(cigar>>BAM_CIGAR_SHIFT); // del
						else if (op_next == BAM_CINS) p->indel = cigar>>BAM_CIGAR_SHIFT; // ins
						if (op_next == BAM_CSOFT_CLIP || op_next == BAM_CREF_SKIP || op_next == BAM_CHARD_CLIP)
							p->is_tail = 1; // tail
					} else p->is_tail = 1; // this is the last operation; set tail
				}
			}
			x += l; y += l;
		} else if (op == BAM_CDEL) { // then set ->is_del
			if (x + l > pos) {
				p->indel = 0; p->is_del = 1;
				p->qpos = y + (pos - x);
			}
			x += l;
		} else if (op == BAM_CREF_SKIP) x += l;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		is_restart = (op == BAM_CREF_SKIP || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP);
		if (x > pos) {
			if (op == BAM_CREF_SKIP) ret = 0; // then do not put it into pileup at all
			break;
		}
	}
	assert(x > pos);
	return ret;
}

/* --- END: Auxiliary functions */

struct __bam_plbuf_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	bam_pileup_f func;
	void *func_data;
	int32_t tid, pos, max_tid, max_pos;
	int max_pu, is_eof;
	bam_pileup1_t *pu;
	int flag_mask;
    int mapq_thresh;

};

void bam_plbuf_reset(bam_plbuf_t *buf)
{
	lbnode_t *p, *q;
	buf->max_tid = buf->max_pos = -1;
	buf->tid = buf->pos = 0;
	buf->is_eof = 0;
	for (p = buf->head; p->next;) {
		q = p->next;
		mp_free(buf->mp, p);
		p = q;
	}
	buf->head = buf->tail;
}

void bam_plbuf_set_mapq_thresh(bam_plbuf_t *buf, int thresh)
{
        if (thresh < 0) buf->mapq_thresh = 0;
            else buf->mapq_thresh = thresh;
}


void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask)
{
	if (mask < 0) buf->flag_mask = BAM_DEF_MASK;
	else buf->flag_mask = BAM_FUNMAP | mask;
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
	bam_plbuf_t *buf;
	buf = (bam_plbuf_t*)calloc(1, sizeof(bam_plbuf_t));
	buf->func = func; buf->func_data = data;
	buf->mp = mp_init();
	buf->head = buf->tail = mp_alloc(buf->mp);
	buf->dummy = mp_alloc(buf->mp);
	buf->max_tid = buf->max_pos = -1;
	buf->flag_mask = BAM_DEF_MASK;
	return buf;
}

void bam_plbuf_destroy(bam_plbuf_t *buf)
{
	mp_free(buf->mp, buf->dummy);
	mp_free(buf->mp, buf->head);
	if (buf->mp->cnt != 0)
		fprintf(stderr, "[bam_plbuf_destroy] memory leak: %d. Continue anyway.\n", buf->mp->cnt);
	mp_destroy(buf->mp);
	free(buf->pu);
	free(buf);
}

int get_next_pos(bam_plbuf_t *buf,bamFile fp) {

    if(buf->max_pos != -1){ //this is only true on the first go throught of the loop
        if (buf->head->next) assert(buf->tid <= buf->head->b.core.tid); // otherwise, not sorted
        if (buf->tid < buf->head->b.core.tid) { // come to a new reference sequence
            buf->tid = buf->head->b.core.tid; buf->pos = 0; // jump to the next reference
        }
        else ++buf->pos; // scan contiguously
    }
    if (buf->is_eof && buf->head->next == 0) {
        return -1;
    }

    while(1) {
        while (buf->is_eof || buf->max_tid > buf->tid || (buf->max_tid == buf->tid && buf->max_pos > buf->pos)) {
            int n_pu = 0;
            lbnode_t *p, *q;
            buf->dummy->next = buf->head;
            for (p = buf->head, q = buf->dummy; p->next; q = p, p = p->next) {   //iterate through buf from head to tail in variable p.
                if (p->b.core.tid < buf->tid || (p->b.core.tid == buf->tid && p->end <= buf->pos)) { // then remove from the list
                    q->next = p->next; mp_free(buf->mp, p); p = q;  //copy p's pointer into q, delete p,
                } else if (p->b.core.tid == buf->tid && p->beg <= buf->pos) { // here: p->end > pos; then add to pileup
                    if (n_pu == buf->max_pu) { // then double the capacity
                        buf->max_pu = buf->max_pu? buf->max_pu<<1 : 256;
                        buf->pu = (bam_pileup1_t*)realloc(buf->pu, sizeof(bam_pileup1_t) * buf->max_pu);
                    }
                    buf->pu[n_pu].b = &p->b;
                    if (resolve_cigar(buf->pu + n_pu, buf->pos)) ++n_pu; // skip the read if we are looking at BAM_CREF_SKIP
                }
            }
            buf->head = buf->dummy->next; // dummy->next may be changed
            return n_pu;
            // update tid and pos
        }
        bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
        if (bam_read1(fp, b) >= 0) {
            //if (!(b->core.flag & buf->flag_mask) && !(b->core.qual < buf->mapq_thresh)) { //skip these reads
            if (!(b->core.flag & buf->flag_mask)) { //skip these reads
                bam_copy1(&buf->tail->b, b);
                buf->tail->beg = b->core.pos; buf->tail->end = bam_calend(&b->core, bam1_cigar(b));
                if (!(b->core.tid >= buf->max_tid || (b->core.tid == buf->max_tid && buf->tail->beg >= buf->max_pos))) {
                    fprintf(stderr, "[bam_pileup_core] the input is not sorted. Abort!\n");
                    abort();
                }
                buf->max_tid = b->core.tid; buf->max_pos = buf->tail->beg;
                if (buf->tail->end > buf->pos) { //if the newly added read ends after the current position, make space for another read
                    buf->tail->next = mp_alloc(buf->mp);
                    buf->tail = buf->tail->next;
                }
            }
        } else buf->is_eof = 1;
        free(b->data); free(b);
    }
}

int bam_sspileup_file(bamFile fp1, bamFile fp2, int mask, int thresh, bam_sspileup_f func, void *func_data, FILE *snp_fh)
{
    bam_plbuf_t *buf1, *buf2;
    int ret1 = -1, ret2 = -1;
    buf1 = bam_plbuf_init(NULL, func_data);
    buf2 = bam_plbuf_init(NULL, func_data);

    bam_plbuf_set_mask(buf1, mask);
    bam_plbuf_set_mask(buf2, mask);
    bam_plbuf_set_mapq_thresh(buf1, thresh);
    bam_plbuf_set_mapq_thresh(buf2, thresh);

    while ((ret1 = get_next_pos(buf1, fp1)) >= 0 && (ret2 = get_next_pos(buf2, fp2)) >= 0) {
        do {
            if(buf2->tid > buf1->tid) {
                //no more normal data. Need to catch up
                while(buf2->tid > buf1->tid && (ret1 = get_next_pos(buf1, fp1)) >= 0);
            }
            if(buf1->tid > buf2->tid) {
                //no more tumor data
                while(buf2->tid < buf1->tid && (ret2 = get_next_pos(buf2, fp2)) >= 0);
            }
        } while(buf2->tid != buf1->tid && buf2->pos != buf1->pos);  //removed assert and added this to hopefully deal with tids with no data
        if(ret1 && ret2) {
            func(buf1->tid, buf1->pos, ret1, ret2, buf1->pu, buf2->pu, buf1->func_data, snp_fh);
        }
    }
    bam_plbuf_reset(buf1);  //clear out any remaining data (I hope)
    bam_plbuf_reset(buf2);
    bam_plbuf_destroy(buf1);
    bam_plbuf_destroy(buf2);
    return 0;
}




