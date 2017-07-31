#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include "parmt_config.h"
#include "compearth.h"
#ifdef PARMT_USE_INTEL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#include <ipps.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "iscl/linalg/linalg.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/statistics/statistics.h"
#include "iscl/sorting/sorting.h"


static double median_sorted_array(int n, const double *__restrict x);
static int setG64f(const int npts,
                   const double *__restrict__ Gxx,
                   const double *__restrict__ Gyy,
                   const double *__restrict__ Gzz,
                   const double *__restrict__ Gxy,
                   const double *__restrict__ Gxz,
                   const double *__restrict__ Gyz,
                   double *__restrict__ G);

#define LDG 8

int parmt_computeL1Magnitude64f(const int ldm,
                                const int nobs,
                                const int nmtSolve,
                                const int maxit, const double eps,
                                const double tol, 
                                const int *__restrict__ signalPtr,
                                const int *__restrict__ lags,
                                const int *__restrict__ mtPtr,
                                const double *__restrict__ Gxx,
                                const double *__restrict__ Gyy,
                                const double *__restrict__ Gzz,
                                const double *__restrict__ Gxy,
                                const double *__restrict__ Gxz,
                                const double *__restrict__ Gyz,
                                const double *__restrict__ mts,
                                const double *__restrict__ d,
                                double *__restrict__ mags)
{
    const char *fcnm = "parmt_computeL1Magnitude64f\0";
    double *G, *est, *obs, *wts, xmag;
    double m6[8] __attribute__((aligned(64))); 
    int i, i1, i2, ierr, imt, jmt, k, maxlag, nptsAll, nptsPad;
    bool luseLag;
    const double p = 1.0;
    const double one = 1.0;
    const double zero = 0.0;
    nptsAll = signalPtr[nobs];
    if (nobs < 1)
    {
        printf("%s: Error no observations\n", fcnm);
        return -1;
    }
    // Get the max number of lags
    maxlag = 0;
    luseLag = false;
    if (lags != NULL)
    {
        luseLag = true;
        maxlag = lags[array_absArgmax32i(nobs, lags, &ierr)];
    }
    if (maxlag == 0){luseLag = false;}
    nptsPad = nptsAll + nobs*maxlag;
    G = memory_calloc64f(nptsPad*LDG);
    obs = memory_calloc64f(nptsPad);
    est = memory_calloc64f(nptsPad);
    wts = memory_calloc64f(nptsPad);
    // Insert the Green's functions and data into the moment tensor matrix
    if (!luseLag)
    {
        for (k=0; k<nobs; k++)
        {
            i1 = signalPtr[k];
            i2 = signalPtr[k+1];
            for (i=i1; i<i2; i++)
            {
                G[i*LDG+0] = Gxx[i];
                G[i*LDG+1] = Gyy[i];
                G[i*LDG+2] = Gzz[i];
                G[i*LDG+3] = Gxy[i];
                G[i*LDG+4] = Gxz[i];
                G[i*LDG+5] = Gyz[i];
                obs[i] = d[i];
            }
        }
    }
    else
    {

    }
    for (imt=0; imt<nmtSolve; imt++)
    {
        // Extract the moment tensor
        jmt = mtPtr[imt];
        for (i=0; i<6; i++)
        {
            m6[0] = mts[ldm*jmt+0];
        }
        // Remove the scalar moment
        compearth_CMT2m0(1, 1, m6, &xmag);
        cblas_dscal(6, 1.0/xmag, m6, 1);
        // Create the linear model [G*m]*M0 = [est]*{M0} = {obs}
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    nptsAll, 6, one, G, LDG, m6, 1, zero, est, 1);
        // Solve for magnitude in the L1 norm
        ierr = linalg_irls64f_work(nptsAll, 1,
                                   maxit, p,
                                   eps, tol,
                                   est, obs, &xmag, wts);
        mags[imt] = xmag;
 double xmagMw;
compearth_m02mw(1, KANAMORI_1978, &xmag, &xmagMw);
printf("%f\n", xmagMw);
    }
    memory_free64f(&wts);
    memory_free64f(&G);
    memory_free64f(&obs);
    memory_free64f(&est);
    return 0;
}

int parmt_computeL1StationMagnitude64f(const int ldm, const int nmt,
                                       const int npts, const int nblock,
                                       const double *__restrict__ Gxx,
                                       const double *__restrict__ Gyy,
                                       const double *__restrict__ Gzz,
                                       const double *__restrict__ Gxy,
                                       const double *__restrict__ Gxz,
                                       const double *__restrict__ Gyz,
                                       const double *__restrict__ mts,
                                       const double *__restrict__ d,
                                       double *__restrict__ mags)
{
    const char *fcnm = "parmt_computeL1StationMagnitude64f\0";
    double *G, *M, *R, *Dmat, *res, *resSort, xmag, xopt;
    const double one = 1.0;
    const double zero = 0.0;
    const double negOne =-one;
    int *perm;
    bool lsorted;
    int i, ic, idx, ierr, imt, jdx, jmt, mblock, Mrows, Ncols;
    const int Kcols = 6; // Number of columns of G
    Mrows = npts;
    mblock = nblock;// + computePadding6f(npts);
    G = memory_calloc64f(LDG*npts);
    R = memory_calloc64f(mblock*npts);
    M = memory_calloc64f(mblock*8);
    Dmat = memory_calloc64f(mblock*npts);
    res = memory_calloc64f(npts);
    resSort = memory_calloc64f(npts);
    perm = memory_calloc32i(npts);
    for (i=0; i<npts; i++){perm[i] = i;}
        // Set the observations
        for (i=0; i<npts; i++) 
        {
            for (ic=0; ic<nblock; ic++)
            {
                Dmat[mblock*i+ic] = d[i];
            }
        }
#ifdef __INTEL_COMPILER
    __assume_aligned(G, 64);
#endif
    ierr = setG64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {    
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
int nsorted = 0;
    for (jmt=0; jmt<nmt; jmt=jmt+nblock)
    {
        Ncols = MIN(nblock, nmt - jmt); // Number of columns of M
        // Set the observations
        cblas_dcopy(npts*mblock, Dmat, 1, R, 1);
        for (i=0; i<6; i++)
        {
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic;
                idx = ldm*imt + i;
                jdx = mblock*i;
                M[jdx+ic] = mts[idx];///xmag;
            }
        }
        // Compute R = GM - D = GM - R
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    Mrows, Ncols, Kcols, one, G, LDG,
                    M, mblock, negOne, R, mblock);
        //lsorted = false;
        //for (i=0; i<npts; i++){perm[i] = i;}
        // Sort the columns of the residual matrix
        for (ic=0; ic<Ncols; ic++)
        {
            imt = jmt + ic;
            // Compute the absolute value of the residual
            for (i=0; i<npts; i++)
            {
                res[i] = fabs(R[ic+i*mblock]);
            }
            //cblas_dcopy(npts, &R[ic], mblock, res, 1);
            // Apply the old permutation
            sorting_applyPermutation64f_work(npts, perm, res, resSort);
            lsorted = sorting_issorted64f(npts, resSort, SORT_ASCENDING, &ierr);
            if (!lsorted)
            {
                sorting_argsort64f_work(npts, res, SORT_ASCENDING, perm);
                sorting_applyPermutation64f_work(npts, perm, res, resSort);
            }
            else
            {
                nsorted = nsorted + 1;
            }
            // Compute the median
            xopt = median_sorted_array(npts, resSort);
            //compearth_CMT2mw(1, 1, &mts[8*imt], &xmag);
            compearth_CMT2m0(1, 1, &mts[8*imt], &xmag);
            xopt = xopt*xmag;
            compearth_m02mw(1, KANAMORI_1978, &xopt, &mags[imt]);
        }
    }
printf("%d %d\n", nsorted, nmt);
    memory_free32i(&perm);
    memory_free64f(&G);
    memory_free64f(&R);
    memory_free64f(&M);
    memory_free64f(&Dmat);
    memory_free64f(&res);
    memory_free64f(&resSort);
    return 0;
}

static double median_sorted_array(int n, const double *__restrict x)
{
    double xmed;
    int n2;
    const double half = 0.5;
    n2 = n/2;
    // Even -> average middle two elements
    if (fmod(n, 2) == 0)
    {
        xmed = half*(x[n2-1] + x[n2]);
    }
    // Median is middle element
    else
    {
        xmed = x[n2];
    }
    return xmed;
}

static int setG64f(const int npts,
                   const double *__restrict__ Gxx,
                   const double *__restrict__ Gyy,
                   const double *__restrict__ Gzz,
                   const double *__restrict__ Gxy,
                   const double *__restrict__ Gxz,
                   const double *__restrict__ Gyz,
                   double *__restrict__ G)
{
    const char *fcnm = "setG64f\0";
    int i;
    bool lalign;
    if (Gxx == NULL || Gyy == NULL || Gzz == NULL || Gxy == NULL ||
        Gxz == NULL || G == NULL || npts < 1)
    
    {
        if (Gxx == NULL){printf("%s: Error Gxx is NULL\n", fcnm);}
        if (Gyy == NULL){printf("%s: Error Gyy is NULL\n", fcnm);}
        if (Gzz == NULL){printf("%s: Error Gzz is NULL\n", fcnm);}
        if (Gxy == NULL){printf("%s: Error Gxy is NULL\n", fcnm);}
        if (Gxz == NULL){printf("%s: Error Gxz is NULL\n", fcnm);}
        if (Gyz == NULL){printf("%s: ERror Gyz is NULL\n", fcnm);}
        if (npts < 1)
        {
            printf("%s: No points in Green's functions\n", fcnm);
        }
        return -1;
    }
    lalign = true;
    if (memory_isAligned(Gxx, 64) != 1 || memory_isAligned(Gyy, 64) != 1 ||
        memory_isAligned(Gzz, 64) != 1 || memory_isAligned(Gxy, 64) != 1 ||
        memory_isAligned(Gxz, 64) != 1 || memory_isAligned(Gyz, 64) != 1 ||
        memory_isAligned(G, 64) != 1)
    {
        lalign = false;
    }
    if (lalign)
    {
        #pragma omp simd aligned(G, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz: 64)
        for (i=0; i<npts; i++)
        {
            G[LDG*i+0] = Gxx[i];
            G[LDG*i+1] = Gyy[i];
            G[LDG*i+2] = Gzz[i];
            G[LDG*i+3] = Gxy[i];
            G[LDG*i+4] = Gxz[i];
            G[LDG*i+5] = Gyz[i];
        }
    }
    else
    {
        #pragma omp simd
        for (i=0; i<npts; i++)
        {
            G[LDG*i+0] = Gxx[i];
            G[LDG*i+1] = Gyy[i];
            G[LDG*i+2] = Gzz[i];
            G[LDG*i+3] = Gxy[i];
            G[LDG*i+4] = Gxz[i];
            G[LDG*i+5] = Gyz[i];
        }
    }
    return 0;
}
