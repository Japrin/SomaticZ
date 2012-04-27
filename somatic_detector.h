/*
 * somatic_detector.h
 *
 *  Created on: Apr 18, 2012
 *      Author: japrin (tao2013@gmail.com)
 */

#ifndef SOMATIC_DETECTOR_H_
#define SOMATIC_DETECTOR_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
using namespace std;

#include "sam.h"
#include "bam.h"
#include "kstring.h"
#include "faidx.h"

#define PHRED_CONST (4.343)
#define expPhred(x) (double)exp((double)(-(x))/PHRED_CONST)
#define logPhred(x) (int)((x) < 1 ? (0.5-PHRED_CONST*log(x)) : (-0.5-PHRED_CONST*log(x)))

#include "ksort.h"
//KSORT_INIT_GENERIC(uint32_t)

//static int aa_bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
//static int aa_bam_nt4_nt16_table[] = { 1, 2, 4, 8 };

#define DEFAULT_PLOIDY 2
#define DEFAULT_DOWNSAMPLE 1000

class CNV_aux
{
public:
	uint32_t curPos;
	string chr_nextRegion;
	uint32_t beg_nextRegion;
	uint32_t end_nextRegion;
	int CN_nextRegion;
	string CNV_filename;
	ifstream CNV_in;
	bool isEnd;

public:
	CNV_aux(string CNV_file);
	int getCN(string chr, uint32_t pos);
	void readRegion();
};

class pileup_aux
{
public:
    bam_header_t *h1;
    bam_header_t *h2;
    faidx_t *fai;
    uint32_t format;
    int tid, len, last_pos;
    int mask;
    int mapQ;
    int min_somatic_qual;
    char *ref;
    //copy number file: cnv_fh1 for tumor, cnv_fh2 for normal
	CNV_aux *cnv_h1;
	CNV_aux *cnv_h2;
	int min_depth;
	int total_len;
	int covered_len;
public:
	pileup_aux()
	{
		h1 = NULL;
		h2 = NULL;
		fai = NULL;
		cnv_h1 = NULL;
		cnv_h2 = NULL;
		ref = NULL;
	}
	~pileup_aux()
	{
		if(cnv_h1 != NULL) delete cnv_h1;
	    if(cnv_h2 != NULL) delete cnv_h2;
	    if(h1 != NULL) bam_header_destroy(h1);
		if(h1 != NULL) bam_header_destroy(h2);
		if(fai != NULL) fai_destroy(fai);
		if(ref !=NULL) free(ref);
	}
};

class likelihood_aux
{
public:
	double * llk;	//log-likelihood of every genotype
	int ploidy;
	int g1;	//most likely genotype;
	int g2; //second most likely genotype;
	double g1_llk;
	double g2_llk;
	int c[4];	//counts of 4 bases
	int c_sum;
	int _ref_base;
	int _non_ref_base;
	double rms_mapQ;
	char _non_ref_base_char;
	int MQ0;			//number of reads with MQ 0
	int c_indel;		//number of reads with indel
	int c_gap;			//number of gap covered this site
	int strandRef[2];	//number of Ref allele supporting reads in both direction
	int strandNonRef[2];//number of NonRef allele supporting reads in both direction
	int c_total;		//number of all covered reads, include low MQ, with gap cover this site etc
	uint32_t *info;

public:
	likelihood_aux(int m);
	~likelihood_aux();
	void setNonRef(int nr);
	void collect_baseInfo(int n, const bam_pileup1_t *pl, uint8_t ref_base, pileup_aux *d);
	void calculate_likelihood();
	double getLogLikelihood(double af);
	void maximize();
};
double calculate_score(likelihood_aux *lk1, likelihood_aux *lk2);
int somatic_core(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh);

#endif /* SOMATIC_DETECTOR_H_ */
