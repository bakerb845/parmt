#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "iscl/array/array.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/random/random.h"
#include "iscl/signal/convolve.h"

int parmt_utils_getEffectiveRank(const double relError,
                                 const struct parmtNoiseBasis_struct basis)
{
    double frac;
    int i, rank;
    rank = 1;
    for (i=1; i<basis.nvals; i++)
    {
        // Modification of Demmel 3.9.3 to deal with sqrt(eig)
        if (100.0*pow(basis.sqrtEvals[i], 2)/pow(basis.sqrtEvals[0], 2) < relError)
        {
            rank = i;
            break;
        }
    }
    return rank;
}


int parmt_utils_makeKLNoise(const int rank,
                            const int npts,
                            const struct parmtNoiseBasis_struct basis,
                            double *__restrict__ xn)
{
    const char *fcnm = "parmt_utils_makeKLNoise\0";
    double *r, xscal;
    int ir, ierr, ldz, rankUse;
    if (rank < 1 || npts < 1 || basis.npts != npts || xn == NULL)
    {
        if (rank < 1){log_errorF("%s: rank is too small\n", fcnm);}
        if (npts < 1){log_errorF("%s: no points\n", fcnm);}
        if (basis.npts != npts){log_errorF("%s: npts != basis.npts\n", fcnm);}
        if (xn == NULL){log_errorF("%s: error x is NULL\n", fcnm);}
        return -1;
    }
    rankUse = MAX(1, MIN(basis.nvals, rank));
    ldz = basis.lde;
    r = random_rand64f(rankUse, &ierr);
    array_zeros64f_work(npts, xn);
    for (ir=0; ir<rank; ir++)
    {
        xscal = r[ir]*basis.sqrtEvals[ir];
printf("%f %f\n", r[ir], basis.sqrtEvals[ir]);
        cblas_daxpy(npts, xscal, &basis.evects[ir*ldz], 1, xn, 1);
    }
    memory_free64f(&r);
    return 0; 
}
//============================================================================//
/*!
 * @brief Releases memory on the noise basis struture
 *
 * @param[out] basis   on exit all memory has been released and all variables
 *                     set to NULL or 0.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_utils_freeNoiseBasis(struct parmtNoiseBasis_struct *basis)
{
    memory_free64f(&basis->sqrtEvals);
    memory_free64f(&basis->evects);
    memset(basis, 0, sizeof(struct parmtNoiseBasis_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the eigenvectors and eigenvalues of the noise
 *        autocorrelation matrix which will be further used in the
 *        Karhunen-Loeve expansion.
 *
 * @param[in] npts     number of data points
 * @param[in] data     noise signal from which to compute noise basis.
 *
 * @param[out] basis   on successful exit contains the eigendecomposition 
 *                     of the autocorrelation noise matrix for use in the KL
 *                     expansion.
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_utils_getNoiseBasis64f(const int npts,
                                 const double *__restrict__ data,
                                 struct parmtNoiseBasis_struct *basis)
{
    const char *fcnm = "parmt_utils_getNoiseBasis64f\0";
    double *C, *CtC, *s, xnorm;
    int i, ierr, j, ldc, ldctc, m, n, nrows;
    const enum corrMatrixType_enum type = CORRMTX_AUTOCORRELATION;
    //------------------------------------------------------------------------//
    //  
    // error checking
    ierr = 0;
    memset(basis, 0, sizeof(struct parmtNoiseBasis_struct));
    if (npts < 1 || data == NULL)
    {   
        if (npts < 1){log_errorF("%s: No points in noise signal\n", fcnm);}
        if (data == NULL){log_errorF("%s; Noise signal is NULL\n", fcnm);}
        return -1; 
    }   
    C = NULL;
    CtC = NULL;
    // Normalize the noise energy in the noise signal
    xnorm = cblas_dnrm2(npts, data, 1);
    if (fabs(xnorm) < 1.e-15)
    {
        log_errorF("%s: Error division by zero\n", fcnm);
        return -1;
    }
    xnorm = 1.0/xnorm;
    s = array_copy64f(npts, data, &ierr);
    cblas_dscal(npts, xnorm, s, 1);
    // Figure out the size of the noise correlation matrix
    n = npts;
    m = npts - 1;
    ldc = n;
    if (type == CORRMTX_AUTOCORRELATION)
    {
        nrows = n + m;
        C = memory_calloc64f(ldc*nrows);
    }
    else
    {
        log_errorF("%s: Invalid correlation matrix type\n", fcnm);
        return -1;
    }
    ldctc = n;
    CtC = memory_calloc64f(n*n);
    // Compute the correlation matrix and C'C
    ierr = convolve_corrmtx64f_work(npts, s,
                                    m, ldc, C,
                                    type, true, ldctc, CtC);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to compute correlation matrix\n", fcnm);
        goto ERROR;
    }
    ierr = parmt_utils_getNoiseBasisFromCtC64f(npts, ldctc, CtC, basis);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to compute noise basis\n", fcnm);
    }
ERROR:;
    memory_free64f(&C);
    memory_free64f(&CtC);
    return ierr;
}
 
int parmt_utils_getNoiseBasisFromCtC64f(const int npts, const int ldctc,
                                        const double *__restrict__ CtC,
                                        struct parmtNoiseBasis_struct *basis)
{
    const char *fcnm = "parmt_utils_getNoiseBasisFromCtC64f\0";
    double *CtCw, *w;
    int ierr, i, info, j;
    ierr = 0;
    memset(basis, 0, sizeof(struct parmtNoiseBasis_struct));
    CtCw = memory_calloc64f(npts*ldctc);
    for (i=0; i<npts; i++)
    {
        cblas_dcopy(npts, &CtC[i*ldctc], 1, &CtCw[i*ldctc], 1);
    }
    // Compute the eigendecompsition of the symmetric matrix
    w = memory_calloc64f(npts);
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', npts, CtCw, ldctc, w);
    if (info != 0)
    {
        ierr = 1;
        if (info > 0)
        {
            log_errorF("%s: Failure to converge on %d'th element\n",
                       fcnm, info);
        }
        else
        {
            log_errorF("%s: %d'th parameter is invalid\n", fcnm, info);
        }
        goto ERROR;
    }
    // Set the basis
    basis->lde = npts;
    basis->nvals = npts;
    basis->npts = npts;
    basis->sqrtEvals = memory_calloc64f(basis->nvals);
    basis->evects = memory_calloc64f(basis->lde*basis->nvals);
    // Put into descending order
    for (i=0; i<npts; i++)
    {
        basis->sqrtEvals[i] = sqrt(w[npts-1-i]);
    }
    // Copy corresponding eigenvectors in descending order 
    for (j=0; j<npts; j++)
    {
        for (i=0; i<npts; i++)
        {
            basis->evects[j*basis->lde+i] = CtCw[i*ldctc+(npts-1-j)];
        }
    }
ERROR:;
    memory_free64f(&w);
    memory_free64f(&CtCw);
    return ierr;
}
