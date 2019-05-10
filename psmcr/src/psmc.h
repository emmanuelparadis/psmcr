#ifndef PSMC_H_LH3
#define PSMC_H_LH3

//#include <sys/types.h>
//#include <stdio.h>
//#include <assert.h>

#include <R.h>
#include "khmm.h"

#ifndef FLOAT
#define FLOAT double
#endif

#define PSMC_N_PARAMS 3

#define PSMC_T_INF 1000.0

#define PSMC_F_DECODE   0x1
#define PSMC_F_FULLDEC  0x2
#define PSMC_F_SIMU     0x4
#define PSMC_F_DIVERG   0x8

typedef struct
{
	int L, L_e, n_e;
	char *seq, *name;
} psmc_seq_t;

typedef struct
{
	int n; // $n$ in psmc.tex. number of intervals equals to $n+1$
	int n_free; // number of free lambdas
	int *par_map; // parameter groups

	int n_seqs; // number of sequences
	psmc_seq_t *seqs; // sequences of segregating sites
	int64_t sum_L;
	int sum_n;

	FILE *fpout, *fpcnt;

	FLOAT *inp_pa; // parameters from the input
	FLOAT *inp_ti;

	int flag;
	int n_iters; // number of iterations
	char *pre_fn; // previous results
	char *pattern; // pattern
	FLOAT max_t; // maximum time
	FLOAT tr_ratio; // theta/rho ratio
	FLOAT alpha;
	FLOAT ran_init;
	FLOAT dt0;
} psmc_par_t;

typedef struct
{
	hmm_par_t *hp;
	FLOAT *t; // time boundaries
	int n_params; // equals to (psmc_par_t::n_free + 2)
	FLOAT *params; // all free parameters: \theta_0, \rho_0 and free \lambda_k
	FLOAT *sigma, *post_sigma;

	FLOAT C_pi, C_sigma; // normalization factor
	FLOAT lk; // log-likelihood
	FLOAT lk_lo; // log-likelihood of the left-out sequence
	FLOAT Q0, Q1; // EM-Q before and after maximization
} psmc_data_t;

#ifdef __cplusplus
extern "C" {
#endif

	// functions for parsing command-line options and inputs
	int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars);
    //void psmc_delete_par(psmc_par_t *pp);
    //void psmc_read_segsites(const char *fn, psmc_par_t *pp);
    void psmc_resamp(psmc_par_t *pp);
	// new/delete psmc_data_t
	psmc_data_t *psmc_new_data(psmc_par_t *pp);
    //void psmc_delete_data(psmc_data_t *pd);
	// calculate HMM parameters given population parameters
	void psmc_update_hmm(const psmc_par_t *pp, psmc_data_t *pd);
	void psmc_avg_t(const psmc_par_t *pp, const psmc_data_t *pd, double *avg_t);
	// Baum-Welch algorithm
	FLOAT psmc_em(psmc_par_t *pp, psmc_data_t *data);
	// calculate and print decoding
	void psmc_decode(const psmc_par_t *pp, const psmc_data_t *pd);
	// print parameters
	void psmc_print_data(const psmc_par_t *pp, const psmc_data_t *pd);
	// simulate sequences
    //void psmc_simulate(const psmc_par_t *pp, const psmc_data_t *pd);
	void psmc_update_intv(int n, FLOAT t[], FLOAT max_t, FLOAT alpha, FLOAT *inp_ti);

#ifdef __cplusplus
}
#endif

#endif
