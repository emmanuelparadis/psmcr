/* cor.c    2019-05-11
   This file is part of the R-package `psmcr'.
   See the file ../DESCRIPTION for licensing issues. */

//#include <stdlib.h>
//#include <math.h>
//#include <string.h>

#include <R.h>
#include "psmc.h"

void psmc_update_intv(int n, FLOAT t[], FLOAT max_t, FLOAT alpha, FLOAT *inp_ti)
{
	int k;
	if (inp_ti == 0) {
		FLOAT beta;
		beta = log(1.0 + max_t / alpha) / n; // beta controls the sizes of intervals
		for (k = 0; k < n; ++k)
			t[k] = alpha * (exp(beta * k) - 1);
		t[n] = max_t; t[n+1] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
	} else {
		memcpy(t, inp_ti, (n+1) * sizeof(FLOAT));
		t[n+1] = PSMC_T_INF;
	}
}

psmc_data_t *psmc_new_data(psmc_par_t *pp)
{
	psmc_data_t *pd;
	int k, n = pp->n;
	pd = (psmc_data_t*)calloc(1, sizeof(psmc_data_t));
	pd->n_params = pp->n_free + PSMC_N_PARAMS + ((pp->flag & PSMC_F_DIVERG)? 1 : 0); // one addition parameter for the divergence model
	pd->hp = hmm_new_par(2, n + 1);
	// initialize
	pd->sigma = (FLOAT*)calloc(n+1, sizeof(FLOAT));
	pd->post_sigma = (FLOAT*)calloc(n+1, sizeof(FLOAT));
	pd->t = (FLOAT*)malloc(sizeof(FLOAT) * (n + 2)); // $t_0,\ldots,t_{n+1}$
	pd->params = (FLOAT*)calloc(pd->n_params, sizeof(FLOAT)); // free lambdas + theta + rho
	// initialize t[] and params[]
	if (pp->inp_pa) { // pameters are loaded from a file
	    memcpy(pd->params, pp->inp_pa, sizeof(FLOAT) * pd->n_params); // FIXME: not working for the divergence model
	} else {
		FLOAT theta;
		// initialize psmc_data_t::params[]
		theta = -log(1.0 - (FLOAT)pp->sum_n / pp->sum_L);
		pd->params[0] = theta; // \theta_0
		pd->params[1] = theta / pp->tr_ratio; // \rho_0
		pd->params[2] = pp->max_t;
		for (k = PSMC_N_PARAMS; k != pp->n_free + PSMC_N_PARAMS; ++k) {
			// \lambda_k
			pd->params[k] = 1.0 + (drand48() * 2.0 - 1.0) * pp->ran_init;
			if (pd->params[k] < 0.1) pd->params[k] = 0.1;
		}
		if (pp->flag & PSMC_F_DIVERG) pd->params[pd->n_params - 1] = pp->dt0;
	}
	psmc_update_hmm(pp, pd);
	return pd;
}

/*
void psmc_delete_data(psmc_data_t *pd)
{
	if (pd == 0) return;
	free(pd->sigma); free(pd->post_sigma);
	free(pd->t); free(pd->params);
	hmm_delete_par(pd->hp);
	free(pd);
}
*/

void psmc_update_hmm(const psmc_par_t *pp, psmc_data_t *pd) // calculate the a_{kl} and e_k(b)
{
	FLOAT *q, tmp, sum_t, max_t, *alpha, *beta, *q_aux, *lambda, theta, rho, *t, *tau, dt = 0;
	hmm_par_t *hp = pd->hp;
	int k, l, n = pp->n;
	t = pd->t;
	lambda = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \lambda_k
	alpha = (FLOAT*)malloc(sizeof(FLOAT) * (n + 2)); // \alpha_k
	beta = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \beta_k
	q_aux = (FLOAT*)malloc(sizeof(FLOAT) * n); // for acceleration
	q = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // q_{kl}
	tau = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \tau_k
	// calculate population parameters: \theta_0, \rho_0, \lambda_k and max_t
	theta = pd->params[0]; rho = pd->params[1]; max_t = pd->params[2];

	for (k = 0; k <= n; ++k) {
	    lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
	}
	psmc_update_intv(pp->n, pd->t, max_t, pp->alpha, pp->inp_ti);
	// set the divergence time parameter if necessary
	if (pp->flag & PSMC_F_DIVERG) {
		dt = pd->params[pd->n_params - 1];
		if (dt < 0) dt = 0;
	}
	// calculate \tau_k
	for (k = 0; k <= n; ++k) tau[k] = t[k+1] - t[k];
	// calculate \alpha
	for (k = 1, alpha[0] = 1.0; k <= n; ++k)
		alpha[k] = alpha[k-1] * exp(-tau[k-1] / lambda[k-1]);
	alpha[k] = 0.0;
	// calculate \beta
	for (k = 1, beta[0] = 0.0; k <= n; ++k)
		beta[k] = beta[k-1] + lambda[k-1] * (1.0 / alpha[k] - 1.0 / alpha[k-1]);
	// calculate q_aux
	for (l = 0; l < n; ++l)
		q_aux[l] = (alpha[l] - alpha[l+1]) * (beta[l] - lambda[l] / alpha[l]) + tau[l];
	// calculate C_pi and C_sigma
	for (l = 0, pd->C_pi = 0.0; l <= n; ++l)
		pd->C_pi += lambda[l] * (alpha[l] - alpha[l+1]);
	pd->C_sigma = 1.0 / (pd->C_pi * rho) + 0.5;
	// calculate all the rest
	for (k = 0, sum_t = 0.0; k <= n; ++k) {
		FLOAT *aa, avg_t, ak1, lak, pik, cpik;
		ak1 = alpha[k] - alpha[k+1]; lak = lambda[k]; // just for convenient
		// calculate $\pi_k$, $\sigma_k$ and Lak
		cpik = ak1 * (sum_t + lak) - alpha[k+1] * tau[k];
		pik = cpik / pd->C_pi;
		pd->sigma[k] = (ak1 / (pd->C_pi * rho) + pik / 2.0) / pd->C_sigma;
		// calculate avg_t, the average time point where mutation happens
		avg_t = - log(1.0 - pik / (pd->C_sigma*pd->sigma[k])) / rho;
		if (isnan(avg_t) || avg_t < sum_t || avg_t > sum_t + tau[k]) // in case something bad happens
			avg_t = sum_t + (lak - tau[k] * alpha[k+1] / (alpha[k] - alpha[k+1]));
		// calculate q_{kl}
		tmp = ak1 / cpik;
		for (l = 0; l < k; ++l) q[l] = tmp * q_aux[l]; // q_{kl}, l<k
		q[l++] = (ak1 * ak1 * (beta[k] - lak/alpha[k]) + 2*lak*ak1 - 2*alpha[k+1]*tau[k]) / cpik; // q_{kk}
		if (k < n) {
			tmp = q_aux[k] / cpik;
			for (; l <= n; ++l) q[l] = (alpha[l] - alpha[l+1]) * tmp; // q_{kl}, l>k
		}
		// calculate p_{kl} and e_k(b)
		tmp = pik / (pd->C_sigma * pd->sigma[k]);
		for (aa = hp->a[k], l = 0; l <= n; ++l) aa[l] = tmp * q[l];
		aa[k] = tmp * q[k] + (1.0 - tmp);
		hp->a0[k] = pd->sigma[k];
		hp->e[0][k] = exp(-theta * (avg_t + dt));
		hp->e[1][k] = 1.0 - hp->e[0][k];
		// update sum_lt
		sum_t += tau[k];
		// for (l = 0, tmp = 0.0; l <= n; ++l) tmp += q[l]; fprintf(stderr, "%d\t%lf\n", k, tmp); // for testing only
	}
	// for (l = 0, tmp = 0.0; l <= n; ++l) tmp += hp->a0[l]; fprintf(stderr, "%lf\n", tmp); // for testing only
	// free
	free(q); free(alpha); free(beta); free(q_aux); free(lambda); free(tau);
}
// compute the average time for each interval; a strip-down version of psmc_update_hmm()
void psmc_avg_t(const psmc_par_t *pp, const psmc_data_t *pd, double *avg_t)
{
	FLOAT sum_t, *alpha, *lambda, rho = pd->params[1], *tau, dt = 0;
	int k, n = pp->n;
	lambda = (FLOAT*)alloca(sizeof(FLOAT) * (n + 1)); // \lambda_k
	alpha = (FLOAT*)alloca(sizeof(FLOAT) * (n + 2)); // \alpha_k
	tau = (FLOAT*)alloca(sizeof(FLOAT) * (n + 1)); // \tau_k
	for (k = 0; k <= n; ++k)
		lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
	if (pp->flag & PSMC_F_DIVERG) {
		dt = pd->params[pd->n_params - 1];
		if (dt < 0) dt = 0;
	}
	for (k = 0; k <= n; ++k) tau[k] = pd->t[k+1] - pd->t[k];
	for (k = 1, alpha[0] = 1.0; k <= n; ++k)
		alpha[k] = alpha[k-1] * exp(-tau[k-1] / lambda[k-1]);
	alpha[k] = 0.0;
	for (k = 0, sum_t = 0.0; k <= n; ++k) {
		FLOAT ak1, lak, pik;
		ak1 = alpha[k] - alpha[k+1]; lak = lambda[k]; // just for convenient
		pik = (ak1 * (sum_t + lak) - alpha[k+1] * tau[k]) / pd->C_pi;
		avg_t[k] = - log(1.0 - pik / (pd->C_sigma*pd->sigma[k])) / rho;
		if (isnan(avg_t[k]) || avg_t[k] < sum_t || avg_t[k] > sum_t + tau[k]) // in case something bad happens
			avg_t[k] = sum_t + (lak - tau[k] * alpha[k+1] / (alpha[k] - alpha[k+1]));
		avg_t[k] += dt;
		sum_t += tau[k];
	}
}
