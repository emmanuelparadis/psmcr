/* kmin.h    2019-05-11
   This file is part of the R-package `psmcr'.
   See the file ../DESCRIPTION for licensing issues. */

#ifndef KMIN_H
#define KMIN_H

#define KMIN_RADIUS  0.5
#define KMIN_EPS     1e-7
#define KMIN_MAXCALL 50000

typedef double (*kmin_f)(int, double*, void*);

#ifdef __cplusplus
extern "C" {
#endif

	double kmin_hj(kmin_f func, int n, double *x, void *data, double r, double eps, int max_calls);

#ifdef __cplusplus
}
#endif

#endif
