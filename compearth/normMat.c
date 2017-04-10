#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#ifndef PARMT_USE_ISCL
static double array_norm64f(const int n, const double *__restrict__ x,
                            const enum normType_enum norm,
                            const double p, int *ierr);
#else
#include "iscl/array/array.h"
#endif
/*!
 * @brief Computes the matrix (Frobenius) norm for a matrix
 *
 * @param[in] n       number of matrices
 * @param[in] M       3 x 3 input matrix as a length 9 array [9*n]
 * @param[in] Lnorm   matrix norm
 *                      TWO_NORM (2)
 *                      ONE_NORM (1)
 *                      P_NORM (in this case must set p)
 *                      INFINITY_NORM
 *                      NEGATIVE_INFINITY_NORM
 * @param[in] p       if using a p-norm this is the value for p (> 0)
 *
 * @param[out] mnorm  matrix norms for each matrix
 *
 * @result 0 indicates success 
 *
 */
int compearth_normMat(const int n,
                      const double *M,
                      const enum normType_enum Lnorm,
                      const double p,
                      double *mnorm)
{
    const char *fcnm = "compearth_normMat\0";
    int i, ierr;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        mnorm[i] = array_norm64f(9, &M[9*i], Lnorm, p, &ierr);
        if (ierr != 0)
        {
            printf("%s: Error computing matrix norm!\n", fcnm);
            mnorm[i] = 0.0;
            break;
        }
    }
    return ierr;
}

#ifndef PARMT_USE_ISCL
static double array_norm64f(const int n, const double *__restrict__ x,
                            const enum normType_enum norm,
                            const double p, int *ierr)
{
    const char *fcnm = "array_norm64f\0";
    double xnorm;
    int i;
    *ierr = 0;
    xnorm = 0.0;
    if (n < 1 || x == NULL)
    {
        if (n < 1){printf("%s: Error no data\n", fcnm);}
        if (x == NULL){printf("%s: Error x is NULL\n", fcnm);}
        *ierr = 1;
        return xnorm;
    }
    // 2 norm
    if (norm == TWO_NORM)
    {
        xnorm = cblas_dnrm2(n, x, 1);
    // 1 norm
    }
    else if (norm == ONE_NORM)
    {
        xnorm = cblas_dasum(n, x, 1);
    }
    // p norm
    else if (norm == P_NORM)
    {
        if (p <= 0.0)
        {
            printf("%s: Invalid p value %f\n", fcnm, p);
            *ierr = 1;
            return xnorm;
        }
        #pragma omp simd reduction(+:xnorm)
        for (i=0; i<n; i++)
        {
            xnorm = xnorm + pow(fabs(x[i]), p);
        }
        xnorm = pow(xnorm, 1./p);
    }
    // infinity norm
    else if (norm == INFINITY_NORM)
    {
        xnorm = fabs(x[cblas_idamax(n, x, 1)]);
    }
    // negative infinity norm
    else if (norm == NEGATIVE_INFINITY_NORM)
    {
        xnorm = fabs(x[0]);
        #pragma omp simd reduction(min:xnorm)
        for (i=1; i<n; i++)
        {
            xnorm = fmin(fabs(x[i]), xnorm);
        }
    }
    // not sure - default to 2-norm
    else
    {
        printf("%s: Defaulting to 2-norm\n", fcnm);
        xnorm = cblas_dnrm2(n, x, 1);
    }
    return xnorm;
}
#endif
