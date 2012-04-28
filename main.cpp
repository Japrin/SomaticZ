/*The MIT License

   Copyright (c) 2012-2012 Zhejiang University Cancer Institute (ZUCI).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
/* Contact: Zheng Liangtao <tao2013@gmail.com> */


#include "faidx.h"
#include "khash.h"
#include "kstring.h"
#include "sam.h"
#include "bam.h"

#include "sniper_pileup.h"
#include "somatic_detector.h"

#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
using namespace std;

//static const char *_default_normal_sample_id = "NORMAL";
//static const char *_default_tumor_sample_id = "TUMOR";
//static const char *_default_output_format = "bed";

void usage(const char* progname) {
    /* we dont like basename(3) */
    const char* pn = strrchr(progname, '/');

    if (pn == NULL)
        pn = progname;
    else
        ++pn;

    fprintf(stderr, "\n");
    fprintf(stderr, "%s [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snp_output_file>\n\n", pn);
    fprintf(stderr, "Required Option: \n");
    fprintf(stderr, "        -f FILE   REQUIRED reference sequence in the FASTA format\n\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "        -a INT    cnv file 1 (for tumor)\n");
    fprintf(stderr, "        -b INT    cnv file 2 (for normal)\n");
    fprintf(stderr, "        -c INT    minimum depth required (for both tumor and normal) [10]\n");
    fprintf(stderr, "        -p INT    only bases with quality not less than INT were used for calling [15]\n");
    fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [1]\n");
    fprintf(stderr, "        -Q INT    filtering somatic snv output with somatic quality less than  INT [1]\n");
    fprintf(stderr, "        -u INT    tumor purity [1.00]\n");
    fprintf(stderr, "        -v INT    normal purity [1.00]\n");
    fprintf(stderr, "\n");
}

int main(int argc, char * argv[]) {
	const char *fn_fa = 0;
	pileup_aux *d = new pileup_aux();
	d->tid = -1;
	d->mask = BAM_DEF_MASK;
	d->baseQ = 15;
	d->mapQ = 1;
	d->min_depth = 10;
	d->min_somatic_qual = 1;
	d->total_len = 0;
	d->covered_len = 0;

	int c;
	while ((c = getopt(argc, argv, "f:p:q:Q:a:b:c:u:v:")) >= 0) {
		switch (c) {
			case 'f': fn_fa = optarg; break;
			case 'p': d->baseQ = atoi(optarg);break;
			case 'q': d->mapQ = atoi(optarg);break;
			case 'Q': d->min_somatic_qual = atoi(optarg); break;
			case 'a': d->cnv_h1 = new CNV_aux(optarg); break;
			case 'b': d->cnv_h2 = new CNV_aux(optarg); break;
			case 'c': d->min_depth = atoi(optarg); break;
			case 'u': break;
			case 'v': break;
			default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
		}
	}
	if (optind == argc) {
	        usage(argv[0]);
	        return 1;
	    }
	if (fn_fa) {
		d->fai = fai_load(fn_fa);
	}else {
        fprintf(stderr, "You MUST specify a reference sequence. It isn't optional.\n");
        delete d;
        exit(1);
    }
	bamFile fp1, fp2;
	fp1 = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
	fprintf(stderr, "Normal bam is %s\n", argv[optind+1]);
	fprintf(stderr, "Tumor bam is %s\n", argv[optind]);
	fp2 = bam_open(argv[optind+1], "r");

	d->h1 = bam_header_read(fp1);
	sam_header_parse_rg(d->h1);
	d->h2 = bam_header_read(fp2);
	sam_header_parse_rg(d->h2);
	FILE* snp_fh = fopen(argv[optind+2], "w");

	if(snp_fh)
	{
		bam_sspileup_file(fp1, fp2, d->mask, d->mapQ, somatic_core, d, snp_fh);
	}
	else {
		fprintf(stderr, "Unable to open snp file!!!!!!!!!\n");
		exit(1);
	}
	//output summary
	cerr<<"Total len:\t"<<d->total_len<<endl;
	cerr<<"Covered len (at least "<<d->min_depth<<"):\t"<<d->covered_len<<endl;

    bam_close(fp1);
    bam_close(fp2);
    fclose(snp_fh);
    delete d;
	return 0;
}
