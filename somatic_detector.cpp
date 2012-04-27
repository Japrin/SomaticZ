/*
 * somatic_detector.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: japrin (tao2013@gmail.com)
 */

#include "somatic_detector.h"

static int aa_bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
static int aa_bam_nt4_nt16_table[] = { 1, 2, 4, 8 };

CNV_aux::CNV_aux(string CNV_file)
{
	chr_nextRegion = "";
	beg_nextRegion = 0;
	end_nextRegion = 0;
	CN_nextRegion = 0;
	curPos = 0;
	CNV_filename=CNV_file;
	CNV_in.open(CNV_file.c_str());
	isEnd = false;
}
int CNV_aux::getCN(string chr, uint32_t pos)
{
	curPos = pos;
	while(!isEnd && (chr != chr_nextRegion || curPos > end_nextRegion ))
	{
		readRegion();
	}
	if(curPos < beg_nextRegion || curPos > end_nextRegion) return DEFAULT_PLOIDY;
	else return CN_nextRegion;
}
void CNV_aux::readRegion()
{
	string buff;
	string _chr;
	uint32_t _beg,_end,_CN;
	if(isEnd) return;
	if(getline(CNV_in,buff))
	{
		istringstream s(buff);
		if(s>>_chr>>_beg>>_end>>_CN)
		{
			_beg++;
			chr_nextRegion = _chr;
			beg_nextRegion = _beg;
			end_nextRegion = _end;
			CN_nextRegion = _CN;
			//cerr<<"readRegion, curPos("<<curPos<<")\t"<<chr_nextRegion<<"\t"<<beg_nextRegion<<"\t"<<end_nextRegion<<"\t"<<CN_nextRegion<<endl;
		}else
		{
			cerr<<"wrong format in CNV file("<<CNV_filename<<")";
		}
	}else
	{
		isEnd = true;
		cerr<<"end of file("<<CNV_filename<<"), curPos:"<<curPos<<"stored region is "<<chr_nextRegion<<":"<<beg_nextRegion<<"-"<<end_nextRegion<<endl;
	}
}


likelihood_aux::likelihood_aux(int m)
{
	ploidy = m;
	llk = new double[m+1];
	g1 = -1;
	g2 = -1;
	g1_llk = - numeric_limits<double>::infinity();
	g2_llk = - numeric_limits<double>::infinity();
	for(int i=0; i<4; i++) c[i]=0;
	c_sum = 0;
	_ref_base = -1;
	_non_ref_base = -1;
	rms_mapQ = 0;
	info = NULL;
	MQ0 = 0;
	c_indel = 0;
	c_gap = 0;
	c_total = 0;
	for(int i=0; i<2; i++) strandRef[i] = 0;
	for(int i=0; i<2; i++) strandNonRef[i] = 0;
}
likelihood_aux::~likelihood_aux()
{
	delete [] llk;
	if(info != NULL) delete [] info;
}
void likelihood_aux::setNonRef(int nr)
{
	_non_ref_base = nr;
}
void likelihood_aux::collect_baseInfo(int n, const bam_pileup1_t *pl, uint8_t ref_base, pileup_aux *d)
{
	info = new uint32_t[n];
	int most_base = -1;
	int second_base = -1;
	_ref_base = aa_bam_nt16_nt4_table[ref_base];

	double rms = 0;
	for (int i = 0; i < n; i++)
	{
		const bam_pileup1_t *p = pl + i;
		uint32_t q, x = 0, qq;
		//if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue;
		bool isSkip = false;
		c_total++;
		if(p->is_del) { c_gap++; isSkip = true; }
		if(p->b->core.qual == 0) { MQ0++; }
		if(p->b->core.qual < d->mapQ) { isSkip = true; }
		for(int k=0; k < p->b->core.n_cigar; k++)
		{
			uint32_t cigar = bam1_cigar(p->b)[k];
			int op = cigar&BAM_CIGAR_MASK;
			if(op == BAM_CDEL || op == BAM_CINS) { c_indel++; break; }
		}
		if(isSkip) continue;

		q = (uint32_t)bam1_qual(p->b)[p->qpos];
		x |= (uint32_t)bam1_strand(p->b) << 18 | q << 8 | p->b->core.qual;
		if (p->b->core.qual < q) q = p->b->core.qual;
		x |= q << 24;
		qq = bam1_seqi(bam1_seq(p->b), p->qpos);
		q = aa_bam_nt16_nt4_table[qq? qq : ref_base];
		if (!p->is_del && q < 4) x |= 1 << 21 | q << 16;
		//if (x>>24 < 4 && (x>>8&0x3f) != 0) x = 4<<24 | (x&0xffffff);

		info[c_sum] = x;
		c_sum++;

		//base on read
		int _b = (x >> 16) & 3;
		//base count
		c[_b]++;
		//find out two most frequent bases
		if(most_base == -1)	most_base = _b;
		if(second_base == -1 && _b != most_base) second_base = _b;
		if(c[_b] > c[most_base])
		{
			second_base = most_base;
			most_base = _b;
		}else if(_b != most_base && c[_b] > c[second_base])
		{
			second_base = _b;
		}
		//mapping quality
		rms += p->b->core.qual * p->b->core.qual;
	}
	rms_mapQ = (sqrt((double)rms / c_sum) + .499);
	//ks_introsort(uint32_t, n, info);
	std::random_shuffle(info,info + c_sum);

	//ref and non-ref
	if(_non_ref_base == -1)
	{
		if(_ref_base == most_base) _non_ref_base = second_base;
		else if(_ref_base == second_base) _non_ref_base = most_base;
		else
		{
			_non_ref_base = most_base;
			//cerr<<"ref_base not present in the two most frequent bases"<<endl;
		}
	}
	_non_ref_base_char = _non_ref_base != -1 ? bam_nt16_rev_table[aa_bam_nt4_nt16_table[_non_ref_base]] : 'N';
	for(int i=0; i<c_sum; i++)
	{
		uint32_t x = info[i];
		int _b = (x >> 16) & 3;
		int _s = (x >> 18) & 1;
		if(_b == _ref_base) strandRef[_s]++;
		else if(_b == _non_ref_base) strandNonRef[_s]++;
	}
}
void likelihood_aux::calculate_likelihood()
{
	//calculate likelihood
	for(int g = 0; g <= ploidy; g++)
	{
		double af = double(g) / ploidy;
		llk[g] = getLogLikelihood(af);
	}
	//get most and second most likely genotype g1 and g2
	maximize();
}
double likelihood_aux::getLogLikelihood(double af)
{
	//af, allele frequency of NonRef-allele
	double ret = 0;
	int l = c_sum > DEFAULT_DOWNSAMPLE ? DEFAULT_DOWNSAMPLE : c_sum;
	for (int i=0; i < l; i++)
	{
		int _b = (info[i] >> 16) & 3;
		int _q = info[i] >> 24;
		if(_b == _ref_base)
		{
			double e = pow(10.0, -_q/10.0);
			ret += log((1 - af) * (1- e) + af * e);
		}else if(_b == _non_ref_base)
		{
			double e = pow(10.0, -_q/10.0);
			ret += log((1 - af) * e + af * (1 - e));
		}
	}
	return ret;
}
void likelihood_aux::maximize()
{
	for(int g = 0; g <= ploidy; g++)
	{
		if(llk[g] > g1_llk)
		{
			g2_llk = g1_llk;
			g2 = g1;
			g1_llk = llk[g];
			g1 = g;
		}else if(llk[g] > g2_llk)
		{
			g2_llk = llk[g];
			g2 = g;
		}
	}
}

double calculate_score(likelihood_aux *lk1, likelihood_aux *lk2 )
{
	//maximize L[1](g^)L[2](g^)
	int CN1 = lk1->ploidy;
	int CN2 = lk2->ploidy;
	double max_llk_sum = - numeric_limits<double>::infinity();
	int max_i=-1;
	int max_j=-1;
	for(int i=0;i<=CN1;i++)
	{
		for(int j=0;j<=CN2;j++)
		{
			if(fabs(double(i)/double(CN1) - double(j)/double(CN2)) > 1e-6) continue;
			double tmp = lk1->llk[i] + lk2->llk[j];
			if(tmp > max_llk_sum)
			{
				max_llk_sum = tmp;
				max_i = i;
				max_j = j;
			}
		}
	}
	//calculate score
	if(max_i == -1 || max_j == -1)
	{
		cerr<<"calculate_score error:"<<endl;
		return 0;
	}
	double D = -2 * (max_llk_sum - lk1->g1_llk - lk2->g1_llk);
	return D;
}

int somatic_core(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh)
{
		pileup_aux *d = (pileup_aux*)data;
		if (d->fai && (int)tid != d->tid) {
			free(d->ref);
			d->ref = fai_fetch(d->fai, d->h1->target_name[tid], &d->len);
			d->tid = tid;
		}
		int rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
		string chr(d->h1->target_name[tid]);
		int CN1 = DEFAULT_PLOIDY;
		if(d->cnv_h1 != NULL)
			CN1 = d->cnv_h1->getCN(chr,pos+1);
		int CN2 = DEFAULT_PLOIDY;
		if(d->cnv_h2 != NULL)
			CN2 = d->cnv_h2->getCN(chr,pos+1);
		d->total_len++;
		if(CN1 == 0 || CN2 == 0) return 0;
		if(n1 < d->min_depth || n2 < d->min_depth || rb - 'N' == 0) { return 0; }
		d->covered_len++;

		//calculate likelihood ratio score
		likelihood_aux * lk1 = new likelihood_aux(CN1);
		lk1->collect_baseInfo(n1,pl1,bam_nt16_table[rb],d);
		likelihood_aux * lk2 = new likelihood_aux(CN2);
		lk2->collect_baseInfo(n2,pl2,bam_nt16_table[rb],d);
		if(lk1->_non_ref_base != -1)
			lk2->setNonRef(lk1->_non_ref_base);
		else
			lk1->setNonRef(lk2->_non_ref_base);
		lk1->calculate_likelihood();
		lk2->calculate_likelihood();

		likelihood_aux * lk3 = new likelihood_aux(DEFAULT_PLOIDY);
		lk3->collect_baseInfo(n1,pl1,bam_nt16_table[rb],d);
		likelihood_aux * lk4 = new likelihood_aux(DEFAULT_PLOIDY);
		lk4->collect_baseInfo(n2,pl2,bam_nt16_table[rb],d);
		if(lk3->_non_ref_base != -1)
			lk4->setNonRef(lk3->_non_ref_base);
		else
			lk3->setNonRef(lk4->_non_ref_base);
		lk3->calculate_likelihood();
		lk4->calculate_likelihood();

		if(lk1->g1 == -1 || lk2->g1 == -1)
			cerr<<"likelihood exceed -infinity!\t"
				<<chr<<"\t"<<pos+1<<"\t"<<char(rb)<<"\t"<<lk1->_non_ref_base_char<<"\t"<<lk2->_non_ref_base_char<<"\t"
				<<lk1->c[lk1->_ref_base]<<"\t"<<(lk1->_non_ref_base == -1 ?0:lk1->c[lk1->_non_ref_base])<<"\t"
				<<lk2->c[lk2->_ref_base]<<"\t"<<(lk2->_non_ref_base == -1 ?0:lk2->c[lk2->_non_ref_base])<<"\t"
				<<CN1<<"\t"<<CN2<<endl;
		else
		{
			double Dp1 = calculate_score(lk1, lk2);
			double Dp2 = calculate_score(lk3, lk4);

			char NonRefChar = lk1->_non_ref_base_char != 'N' ? lk1->_non_ref_base_char : lk2->_non_ref_base_char;
			if((Dp1 >= d->min_somatic_qual || Dp2 >= d->min_somatic_qual) && NonRefChar != 'N')
			{
				if(lk1->_non_ref_base != -1 && lk2->_non_ref_base != -1 && lk1->_non_ref_base != lk2->_non_ref_base)
					fprintf(stderr,"Different non-ref base in tumor and normal at (%s:%d) !\t%d\t%d\n",chr.c_str(),pos+1,lk1->_non_ref_base,lk2->_non_ref_base);
				fprintf(snp_fh,"%s\t%d\t%d\t%c\t%c"		//chr, beg, end, ref allele, non-ref allele
								"\t%4.4f\t%d/%d\t%d/%d"	//score, genotype by aneuploidy hypothesis
								"\t%4.4f\t%d/%d\t%d/%d"	//score, genotype by diploidy hypothesis
								"\t%d\t%d\t%d\t%d"		//ref, non-ref, ref, non-ref
								"\t%d\t%d"				//Copy Number
								"\t%4.4f\t%4.4f"		//rms MQ
								"\t%d\t%d\t%d\t%d\t[%d,%d,%d,%d]\t%d\t%d\t%d\t%d\t[%d,%d,%d,%d]"	//MQ0, indel, gap, n, strand
								"\n",
								chr.c_str(),pos,pos+1,rb,NonRefChar,
								Dp1,lk1->g1,lk1->ploidy,lk2->g1,lk2->ploidy,
								Dp2,lk3->g1,lk3->ploidy,lk4->g1,lk4->ploidy,
								lk1->c[lk1->_ref_base],lk1->_non_ref_base == -1 ?0:lk1->c[lk1->_non_ref_base],lk2->c[lk2->_ref_base],lk2->_non_ref_base == -1 ?0:lk2->c[lk2->_non_ref_base],
								//lk1->c_sum,lk2->c_sum,
								CN1,CN2,
								lk1->rms_mapQ,lk2->rms_mapQ,
								lk1->MQ0,lk1->c_indel,lk1->c_gap,n1,lk1->strandRef[0],lk1->strandRef[1],lk1->strandNonRef[0],lk1->strandNonRef[1],
								lk2->MQ0,lk2->c_indel,lk2->c_gap,n2,lk2->strandRef[0],lk2->strandRef[1],lk2->strandNonRef[0],lk2->strandNonRef[1]
				);
			}
		}

		delete lk1;
		delete lk2;
		delete lk3;
		delete lk4;

		return 0;
}
