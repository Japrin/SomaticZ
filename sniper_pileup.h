/*
 * sniper_pileup.h
 *
 * Code from SomaticSniper
 *
 */

#ifndef SNIPER_PILEUP_H_
#define SNIPER_PILEUP_H_

typedef int (*bam_sspileup_f)(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE* snp_fh);

int get_next_pos(bam_plbuf_t *buf,bamFile fp);
int bam_sspileup_file(bamFile fp1, bamFile fp2, int mask, int thresh, bam_sspileup_f func, void *func_data, FILE *snp_fh);


#endif /* SNIPER_PILEUP_H_ */
