/* Rpsmc.c    2023-03-16
   This file is part of the R-package `psmcr'.
   See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "psmc.h"

#define PSMC_VERSION "0.6.5-r67" // from cli.c

static char conv_table_DNAbin[256] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    0, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2,
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2,
    2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

/* the original conversion table working with upper-
   and lowercase letters
static char conv_table[256] = {
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    0, 1, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 0, 2, 0,  2, 2, 2, 0,  2, 2, 2, 1,  2, 1, 2, 2,
    2, 2, 1, 1,  0, 2, 2, 1,  2, 1, 2, 2,  2, 2, 2, 2,
    2, 0, 2, 0,  2, 2, 2, 0,  2, 2, 2, 1,  2, 1, 2, 2,
    2, 2, 1, 1,  0, 2, 2, 1,  2, 1, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2
}; */

// parse a pattern like "4+5*3+4"
// the returned array holds which parameters are linked together
// number of parameters and number of free parameters will be also returned

int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars)
{
	char *q, *p, *tmp;
	int top = 0, *stack = (int*)malloc(sizeof(int) * 0x100);
	int *pars_map, k, l, i;
	p = q = tmp = strdup(pattern);
	k = 1;
	while (1) {
	    //assert(isdigit(*p) || *p == '*' || *p == '+' || *p == '\0'); // allowed characters
		if (*p == '+' || *p == '\0') {
			int is_end = (*p == 0)? 1 : 0;
			*p++ = '\0';
			l = atoi(q); q = p;
			for (i = 0; i < k; ++i) {
				stack[top++] = l;
				//assert(top <= 0xff);
			}
			k = 1;
			if (is_end) break;
		} else if (*p == '*') {
			*p = '\0';
			k = atoi(q); // number of repeats
			*p++ = '*'; q = p;
		} else ++p;
	}
	for (k = l = 0; k != top; ++k) l += stack[k];
	*n_pars = l - 1; *n_free = top;
	pars_map = (int*)malloc(sizeof(int) * (*n_pars + 1));
	for (k = i = 0; k != top; ++k)
		for (l = 0; l < stack[k]; ++l)
			pars_map[i++] = k;
	free(tmp); free(stack);
	return pars_map;
}

/*
void psmc_delete_par(psmc_par_t *par)
{
	int i;
	if (par == 0) return;
	Rprintf("...1");
	free(par->par_map); free(par->inp_pa);
	Rprintf("...2");
	free(par->pre_fn); free(par->pattern);
	Rprintf("...3");
	fclose(par->fpout);
	Rprintf("...4");
	for (i = 0; i != par->n_seqs; ++i) {
		free(par->seqs[i].name);
		free(par->seqs[i].seq);
	}
	Rprintf("...5");
	if (par->fpcnt) fclose(par->fpcnt);
	Rprintf("...6");
	free(par->seqs);
	free(par);
}
*/

void count_CpG(SEXP SEQ, int nint, char *cpgfilename)
{
    PROTECT(SEQ = coerceVector(SEQ, VECSXP));
    int NSEQ = LENGTH(SEQ); // nb of chromosomes

    int step = 100, k = 0, i, s, n;
    int z[] = {0, 0, 0, 0, 0};
    unsigned char a, b, *x;

    FILE *f;
    f = fopen(cpgfilename, "w");

    i = 5; // use the address of 'i' to write the number of values
    fwrite(&i, 4, 1, f);

    for (s = 0; s < NSEQ; s++) {
	x = RAW(VECTOR_ELT(SEQ, s));
	n = LENGTH(VECTOR_ELT(SEQ, s));

	step = n/nint;
	i = (n + step - 1)/step; // number of windows
	fwrite(&i, 4, 1, f);

	i = 0;
	while (i < n - 1) {
	    a = x[i];
	    b = x[i + 1];
	    if ((a == 40 || x[i] == 48) && (a == 72 || b == 192)) {
		++z[0];
		if (a == 48) ++z[1]; // TS in CpG
		if (b == 192) ++z[1]; // id.
		// because of the definition of CpG here, there cannot be a TV
		i += 2;
		k += 2;
	    } else {
		if (a == 48 || a == 192) ++z[3]; // TS in non-CpG
		if (a == 80 || a == 96 || a == 144 || a == 160) ++z[4]; // TV
		i++;
		k++;
	    }
	    if (k >= step) {
		z[2] = n - z[0];
		fwrite(z, 4, 5, f);
		memset(z, 0, 5 * sizeof(int));
		k = 0;
	    }
	}
    }
    fclose(f);
    UNPROTECT(1);
}

SEXP Rpsmc_C(SEXP SEQ, SEXP para)
{
    unsigned char *x;
    int i, j, n, decoding, quiet, resample;
    psmc_par_t *pp;
    psmc_data_t *pd;

    PROTECT(SEQ = coerceVector(SEQ, VECSXP));
    PROTECT(para = coerceVector(para, VECSXP));
    int NSEQ = length(SEQ); // nb of chromosomes
    SEXP seqnames = getAttrib(SEQ, R_NamesSymbol);

    pp = (psmc_par_t*)calloc(1, sizeof(psmc_par_t));

    /* set the parameters */
    pp->n_iters = INTEGER(VECTOR_ELT(para, 2))[0];
    pp->max_t = REAL(VECTOR_ELT(para, 1))[0];
    pp->tr_ratio = REAL(VECTOR_ELT(para, 3))[0];
    pp->pattern = (char *)CHAR(STRING_ELT(VECTOR_ELT(para, 0), 0));
    pp->n_seqs = 0; // initialize
    pp->fpout = fopen((char *)CHAR(STRING_ELT(VECTOR_ELT(para, 7), 0)), "w"); // the output file
    pp->alpha = 0.1;
    pp->dt0 = -1.0;
    pp->flag = 0;
    /* pp->inp_pa = 0;
       pp->inp_ti = 0; */

    decoding = INTEGER(VECTOR_ELT(para, 4))[0];
    quiet = INTEGER(VECTOR_ELT(para, 5))[0];
    resample = INTEGER(VECTOR_ELT(para, 6))[0];

    pp->par_map = psmc_parse_pattern(pp->pattern, &pp->n_free, &pp->n);

    if (decoding) {
	char *cpgfilename;
	cpgfilename = (char *)CHAR(STRING_ELT(VECTOR_ELT(para, 8), 0));
	pp->flag |= PSMC_F_DECODE;
	count_CpG(SEQ, pp->n, cpgfilename);
	pp->fpcnt = fopen(cpgfilename, "r");
    }
    /* if (*simulate) pp->flag |= PSMC_F_SIMU; */

    /* input the sequence data */
    char *s;
    for (i = 0; i < NSEQ; i++) {
	x = RAW(VECTOR_ELT(SEQ, i));
	n = LENGTH(VECTOR_ELT(SEQ, i));

	int L_e = 0, n_e = 0, L = 0;

    //if ((pp->n_seqs & 0xff) == 0) {
    /* always need to reallocate memory (EP 2019-04-01) */
	pp->seqs = (psmc_seq_t*)realloc(pp->seqs, sizeof(psmc_seq_t) * (pp->n_seqs + 0x100));
	s = pp->seqs[pp->n_seqs].seq = (char*)calloc(n, 1);

	for (j = 0; j < n; ++j) {
	    char c = conv_table_DNAbin[x[j]];
	    s[L++] = c;
	    if (c < 2) {
		++L_e;
		if (c == 1) ++n_e;
	    }
	}

	pp->seqs[pp->n_seqs].name = (char *)CHAR(STRING_ELT(seqnames, i));

	pp->seqs[pp->n_seqs].L_e = L_e;
	pp->seqs[pp->n_seqs].n_e = n_e;
	pp->seqs[pp->n_seqs++].L = L;
    }

    pp->sum_n = pp->sum_L = 0;
    for (int i = 0; i < pp->n_seqs; i++) {
	pp->sum_n += pp->seqs[i].n_e;
	pp->sum_L += pp->seqs[i].L_e;
    }

    /* set the data */
    pd = psmc_new_data(pp);

    fprintf(pp->fpout, "CC\n");
    fprintf(pp->fpout, "CC\tBrief Description of the file format:\n");
    fprintf(pp->fpout, "CC\t  CC  comments\n");
    fprintf(pp->fpout, "CC\t  MM  useful-messages\n");
    fprintf(pp->fpout, "CC\t  RD  round-of-iterations\n");
    fprintf(pp->fpout, "CC\t  LL  \\log[P(sequence)]\n");
    fprintf(pp->fpout, "CC\t  QD  Q-before-opt Q-after-opt\n");
    fprintf(pp->fpout, "CC\t  TR  \\theta_0 \\rho_0\n");
    fprintf(pp->fpout, "CC\t  RS  k t_k \\lambda_k \\pi_k \\sum_{l\\not=k}A_{kl} A_{kk}\n");
    fprintf(pp->fpout, "CC\t  DC  begin end best-k t_k+\\Delta_k max-prob\n");
    fprintf(pp->fpout, "CC\n");
    fprintf(pp->fpout, "MM\tVersion: %s\n", PSMC_VERSION);
    fprintf(pp->fpout, "MM\tpattern:%s, n:%d, n_free_lambdas:%d\n", pp->pattern, pp->n, pp->n_free);
    fprintf(pp->fpout, "MM\tn_iterations:%d, skip:1, max_t:%g, theta/rho:%g\n", pp->n_iters, pp->max_t, pp->tr_ratio);
    fprintf(pp->fpout, "MM\tis_decoding:%d\n", (pp->flag&PSMC_F_DECODE)? 1 : 0);
    fprintf(pp->fpout, "MM\tn_seqs:%d, sum_L:%lld, sum_n:%d\n", pp->n_seqs, (long long)pp->sum_L, pp->sum_n);

    if (resample) psmc_resamp(pp);

    /* start the iterations */
    fprintf(pp->fpout, "RD\t0\n");
    psmc_print_data(pp, pd);
    for (i = 0; i < pp->n_iters; i++) {
	if (!quiet) Rprintf("\rIteration %d/%d...", i + 1, pp->n_iters);
	psmc_em(pp, pd);
	fprintf(pp->fpout, "RD\t%d\n", i+1);
	psmc_print_data(pp, pd);
	R_CheckUserInterrupt();
    }
    if (!quiet) Rprintf(" Done.\n");
    if ((pp->flag & PSMC_F_DECODE) || pp->fpcnt) {
	if (!quiet) Rprintf("Starting decoding...");
	psmc_decode(pp, pd);
	if (!quiet) Rprintf(" Done.\n");
    }
    /* if (pp->flag & PSMC_F_SIMU) psmc_simulate(pp, pd); */
    free(pd); //psmc_delete_data(pd);
    free(pp); //psmc_delete_par(pp);
    UNPROTECT(2);
    return R_NilValue;
}

SEXP vcf2DNAbinFromRaw(SEXP x, SEXP refgenome, SEXP sizeofheader);
SEXP seqBinning_DNAbin(SEXP x, SEXP binsize);

static R_CallMethodDef Call_entries[] = {
    {"Rpsmc_C", (DL_FUNC) &Rpsmc_C, 2},
    {"vcf2DNAbinFromRaw", (DL_FUNC) &vcf2DNAbinFromRaw, 4},
    {"seqBinning_DNAbin", (DL_FUNC) &seqBinning_DNAbin, 2},
    {NULL, NULL, 0}
};

void R_init_psmcr(DllInfo *info)
{
    R_registerRoutines(info, NULL, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
