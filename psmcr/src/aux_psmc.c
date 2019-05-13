/* aux_psmc.c    2019-05-11
   This file is part of the R-package `psmcr'.
   See the file ../DESCRIPTION for licensing issues. */

//#include <stdio.h>
//#include <assert.h>
//#include <math.h>
//#include <time.h>
//#include <string.h>

#include <R.h>
#include "psmc.h"

void psmc_resamp(psmc_par_t *pp)
{
	int i, n_seqs = 0;
	psmc_seq_t *seqs = 0;
	int64_t L = 0, L_ori = 0;
	for (i = 0; i != pp->n_seqs; ++i) L_ori += pp->seqs[i].L;
	while (1) {
		GetRNGstate();
		psmc_seq_t *ns, *s = pp->seqs + (int)(pp->n_seqs * unif_rand());
		PutRNGstate();
		int tmp1 = L_ori - L;
		int tmp2 = L + s->L - L_ori;
		if (tmp2 <= 0 || (tmp2 > 0 && tmp1 > 0 && tmp2 < tmp1)) { // add seq
			if ((n_seqs&0xff) == 0)
				seqs = (psmc_seq_t*)realloc(seqs, sizeof(psmc_seq_t) * (n_seqs + 0x100));
			ns = seqs + n_seqs;
			ns->name = strdup(s->name);
			ns->seq = (char*)malloc(s->L);
			memcpy(ns->seq, s->seq, s->L);
			ns->L = s->L;
			ns->L_e = s->L_e;
			ns->n_e = s->n_e;
			L += ns->L;
			++n_seqs;
		}
		if (tmp1 >= 0 && tmp2 >= 0) break;
	}
	// delete old information
	/* The block below is silenced because it makes R crash
	   (probably, the memory should not be freed before
	   returning the C call to R --- EP 2019-04-19). */
	/*
	for (i = 0; i != pp->n_seqs; ++i) {
		free(pp->seqs[i].name);
		free(pp->seqs[i].seq);
	}
	free(pp->seqs);
	*/
	// fill up new information
	pp->n_seqs = n_seqs;
	pp->seqs = seqs;
	pp->sum_n = pp->sum_L = 0;
	for (i = 0; i != pp->n_seqs; ++i) {
		pp->sum_n += pp->seqs[i].n_e;
		pp->sum_L += pp->seqs[i].L_e;
	}
}

void psmc_print_data(const psmc_par_t *pp, const psmc_data_t *pd)
{
	int k;
	FLOAT n_recomb = pp->sum_L / pd->C_sigma;
	FLOAT /* theta0, rho0, */ *lambda, sum;
	lambda = (FLOAT*)malloc(sizeof(FLOAT) * (pp->n + 1));
	// apparently unused (EP 2019-04-19):
	/* theta0 = pd->params[0]; rho0 = pd->params[1]; */
	for (k = 0; k <= pp->n; ++k)
		lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
	fprintf(pp->fpout, "LK\t%lf\n", pd->lk);
	fprintf(pp->fpout, "QD\t%lf -> %lf\n", pd->Q0, pd->Q1);
	// calculate Relative Information (KL distnace)
	for (k = 0, sum = 0.0; k <= pp->n; ++k)
		sum += pd->sigma[k] * log(pd->sigma[k] / pd->post_sigma[k]);
	fprintf(pp->fpout, "RI\t%.10lf\n", sum);
	// print other parameters
	fprintf(pp->fpout, "TR\t%lf\t%lf\n", pd->params[0], pd->params[1]);
	fprintf(pp->fpout, "MT\t%lf\n", pd->params[2]);
	if (pp->flag & PSMC_F_DIVERG)
		fprintf(pp->fpout, "DT\t%lf\n", pd->params[pd->n_params - 1]);
	//
	fprintf(pp->fpout, "MM\tC_pi: %lf, n_recomb: %lf\n", pd->C_pi, n_recomb);
	for (k = 0; k <= pp->n; ++k)
		fprintf(pp->fpout, "RS\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", k, pd->t[k], lambda[k], n_recomb * pd->hp->a0[k],
				pd->sigma[k], pd->post_sigma[k]);
	fprintf(pp->fpout, "PA\t%s", pp->pattern);
	for (k = 0; k != pd->n_params; ++k)
		fprintf(pp->fpout, " %.9lf", pd->params[k]);
	if (pp->inp_ti)
		for (k = 0; k <= pp->n; ++k)
			fprintf(pp->fpout, " %.9lf", pp->inp_ti[k]);
	fprintf(pp->fpout, "\n//\n");
	fflush(pp->fpout);
	free(lambda);
}

void psmc_decode(const psmc_par_t *pp, const psmc_data_t *pd)
{
	hmm_par_t *hp = pd->hp;
	int i, k, prev, start;
	FLOAT p, q, *t;
	double *cnt = 0, theta = pd->params[0];
	int32_t n_cnt;
	// compute the time intervals and the coalescent average
	t = (FLOAT*)alloca(sizeof(FLOAT) * (pp->n + 1));
	psmc_avg_t(pp, pd, t);
	for (k = 0; k <= pp->n; ++k) {
		if (t[k] < pd->t[k] || t[k] > pd->t[k+1])
			error("ERROR: (%f <= %f <= %f) does not stand. Contact me if you see this.\n", pd->t[k], t[k], pd->t[k+1]);
		//fprintf(stderr, "ERROR: (%f <= %f <= %f) does not stand. Contact me if you see this.\n", pd->t[k], t[k], pd->t[k+1]);
		fprintf(pp->fpout, "TC\t%d\t%lf\t%lf\t%lf\n", k, pd->t[k] * theta, t[k] * theta, pd->t[k+1] * theta);
	}
	if (pp->fpcnt) {
		fread(&n_cnt, 4, 1, pp->fpcnt); // read the number of counts per base
		cnt = (double*)alloca((pp->n + 1) * n_cnt * sizeof(double));
		memset(cnt, 0, (pp->n + 1) * n_cnt * sizeof(double));
	}
	// the core part
	hmm_pre_backward(hp);
	for (i = 0; i != pp->n_seqs; ++i) {
		hmm_data_t *hd;
		psmc_seq_t *s = pp->seqs + i;
		char *seq = (char*)calloc(s->L+1, 1);
		memcpy(seq, s->seq, s->L);
		hd = hmm_new_data(s->L, seq, hp);
		hmm_forward(hp, hd);
		hmm_backward(hp, hd);
		if (!(pp->flag & PSMC_F_FULLDEC) && (pp->flag & PSMC_F_DECODE)) { // posterior decoding
		    int *x;//, kl;
			/* NOTE: kl is declared here and assigned a few lines below,
			   but this is not apparently not needed. (EP 2019-04-10) */
			hmm_post_decode(hp, hd);
			/* show path */
			x = hd->p;
			start = 1; prev = x[1];
			p = hd->f[1][prev] * hd->b[1][prev] * hd->s[1];
			for (k = 2; k <= s->L; ++k) {
				if (prev != x[k]) {
				    //kl = pp->par_map[prev]; see NOTE above
					fprintf(pp->fpout, "DC\t%s\t%d\t%d\t%d\t%lf\t%.3lf\n", s->name, start, k-1, prev, t[prev] * theta, p);
					prev = x[k]; start = k; p = 0.0;
				}
				q = hd->f[k][x[k]] * hd->b[k][x[k]] * hd->s[k];
				if (p < q) p = q;
			}
			fprintf(pp->fpout, "DC\t%s\t%d\t%d\t%d\t%.3lf\t%.2lf\n", s->name, start, k-1, prev, t[prev] * theta, p);
			fflush(pp->fpout);
		} else if (pp->flag & PSMC_F_DECODE) { // full decoding
			FLOAT *prob = (FLOAT*)malloc(sizeof(FLOAT) * hp->n);
			for (k = 1; k <= s->L; ++k) {
				int l;
				FLOAT p, *fu, *bu1, *eu1; // p is the recombination probability?
				if (k < s->L) {
					p = 0.0; fu = hd->f[k]; bu1 = hd->b[k+1]; eu1 = hp->e[(int)hd->seq[k+1]];
					for (l = 0; l < hp->n; ++l)
						p += fu[l] * hp->a[l][l] * bu1[l] * eu1[l];
					p = 1.0 - p;
				} else p = 0.0;
				hmm_post_state(hp, hd, k, prob);
				fprintf(pp->fpout, "DF\t%d\t%lf", k, p);
				for (l = 0; l < hp->n; ++l)
					fprintf(pp->fpout, "\t%.4f", prob[l]);
				fprintf(pp->fpout, "\n");
			}
			free(prob);
		}
		if (pp->fpcnt) { // very similar to full decoding above
			int32_t *cnt1, l, min_l;
			FLOAT *prob = (FLOAT*)malloc(sizeof(FLOAT) * hp->n);
			fread(&l, 4, 1, pp->fpcnt);
			/* NOTE: the test below does not seem necessary since the file with CpG counts does
			   include the chromosome lengths (see code in psmc/utils/cntcpg.c). (EP 2019-04-10) */
			/* if (l != s->L) */
			/*     Rprintf("WARNING: chromosome length difference %d != %d\n", s->L, l); */
			/*     //fprintf(stderr, "WARNING: chromosome length difference %d != %d\n", s->L, l); */
			cnt1 = calloc(l * n_cnt, 4);
			fread(cnt1, n_cnt * l, 4, pp->fpcnt);
			min_l = s->L < l? s->L : l;
			for (k = 1; k <= min_l; ++k) {
				int j, l;
				hmm_post_state(hp, hd, k, prob);
				for (l = 0; l < hp->n; ++l)
					for (j = 0; j < n_cnt; ++j)
						cnt[l*n_cnt + j] += prob[l] * cnt1[(k-1)*n_cnt + j];
			}
			free(cnt1); free(prob);
		}
		/* free */
		hmm_delete_data(hd);
		free(seq);
	}
	if (pp->fpcnt) {
		for (i = 0; i < hp->n; ++i) {
			fprintf(pp->fpout, "CT\t%d", i);
			for (k = 0; k < n_cnt; ++k)
				fprintf(pp->fpout, "\t%f", cnt[i*n_cnt + k]);
			fprintf(pp->fpout, "\n");
		}
	}
}

/*
void psmc_simulate(const psmc_par_t *pp, const psmc_data_t *pd)
{
	const char conv[3] = { 'T', 'K', 'N' };
	int i, k;
	//srand48(time(0));
	for (i = 0; i != pp->n_seqs; ++i) {
		psmc_seq_t *s = pp->seqs + i;
		char *seq;
		seq = hmm_simulate(pd->hp, s->L);
		for (k = 0; k != s->L; ++k)
			seq[k] = (s->seq[k] == 2)? 'N' : conv[(int)seq[k]];
		// print
		fprintf(pp->fpout, "FA\t>%s", s->name);
		for (k = 0; k != s->L; ++k) {
			if (k%60 == 0) fprintf(pp->fpout, "\nFA\t");
			fputc(seq[k], pp->fpout);
		}
		fputc('\n', pp->fpout);
		free(seq);
	}
}
*/
