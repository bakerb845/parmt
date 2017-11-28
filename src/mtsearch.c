#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include "parmt_config.h"
#ifdef PARMT_USE_INTEL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#include <ipps.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "parmt_mtsearch.h"
#include "iscl/array/array.h"
#include "iscl/linalg/linalg.h"
#include "iscl/memory/memory.h"
#include "iscl/signal/convolve.h"

#define LDG 8
#define LDCG 8
#define SIZEOF_DOUBLE sizeof(double)

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

enum covarianceMatrixType_enum
{
    INVERSE_COVARIANCE_MATRIX = 0,
    EIGEN_BASIS = 1
};

// Set the Green's functions matrix
static int computePadding64f(const int n);
static int computePadding32f(const int n);
static int setG32f(const int npts,
                   const float *__restrict__ Gxx,
                   const float *__restrict__ Gyy,
                   const float *__restrict__ Gzz,
                   const float *__restrict__ Gxy,
                   const float *__restrict__ Gxz,
                   const float *__restrict__ Gyz,
                   float *__restrict__ G);
static int setG64f(const int npts,
                   const double *__restrict__ Gxx,
                   const double *__restrict__ Gyy,
                   const double *__restrict__ Gzz,
                   const double *__restrict__ Gxy,
                   const double *__restrict__ Gxz,
                   const double *__restrict__ Gyz,
                   double *__restrict__ G);
static int setGxc64f(const int npts,
                     const double *__restrict__ Gxx,
                     const double *__restrict__ Gyy, 
                     const double *__restrict__ Gzz, 
                     const double *__restrict__ Gxy, 
                     const double *__restrict__ Gxz, 
                     const double *__restrict__ Gyz, 
                     const double *__restrict__ d,
                     double *__restrict__ Gxc);
// Compute C*d
static int computeCd32f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const float *__restrict__ CeInv,
                        const float *__restrict__ d,
                        float *__restrict__ Cd);
static int computeCd64f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const double *__restrict__ CeInv,
                        const double *__restrict__ d,
                        double *__restrict__ Cd);
// Compute C*G
static int computeCG64f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const double *__restrict__ CeInv,
                        const double *__restrict__ G,
                        double *__restrict__ CG);
static int computeCG32f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const float *__restrict__ CeInv,
                        const float *__restrict__ G,
                        float *__restrict__ CG);
// Compute G^T C G
static int computeGtCG32f(const int npts,
                          const float *__restrict__ G,
                          const float *__restrict__ CG,
                          float *__restrict__ GtCG);
static int computeGtCG64f(const int npts,
                          const double *__restrict__ G,
                          const double *__restrict__ CG, 
                          double *__restrict__ GtCG);
// Compute L1 objective function
static int performL1Search64f(const int nmt, const int ldm,
                              const int npts,
                              const int blockSize, const int mblock,
                              const int Mrows, const int Kcols,
                              const int klag,
                              const bool lminLags, const bool lwantLags,
                              const bool lrescale,
                              const double dnorm,
                              const double *__restrict__ Dmat,
                              const double *__restrict__ CG,
                              const double *__restrict__ Cd,
                              const double *__restrict__ mts,
                              int *__restrict__ lags,
                              double *__restrict__ phi,
                              double *__restrict__ var);
// Compute m^T G^T G m
static double computemtGTGmt64f(const double *__restrict__ m,
                                const double *__restrict__ GtG);

#pragma omp declare simd
static float computeLeastSquaresResidual32f(const float *__restrict__ m,
                                            const float *__restrict__ GtCG,
                                            const float *__restrict__ GtCd,
                                            const float dtCd);
#pragma omp declare simd
static double computeLeastSquaresResidual64f(const double *__restrict__ m,
                                             const double *__restrict__ GtCG,
                                             const double *__restrict__ GtCd,
                                             const double dtCd);

/*!
 * @brief Computes the likelihood function
 */
/*
int parmt_mtsearch64f( )
{
    
    return 0;
}
*/


int parmt_mtSearchKL64f(const int ldm, const int nlags, const int nmt, 
                        const int npts, const int ldz, 
                        const int rank,
                        const double sigma,
                        const double *__restrict__ Gxx,
                        const double *__restrict__ Gyy,
                        const double *__restrict__ Gzz,
                        const double *__restrict__ Gxy,
                        const double *__restrict__ Gxz,
                        const double *__restrict__ Gyz,
                        const double *__restrict__ sqrtEig,
                        const double *__restrict__ Z,
                        const double *__restrict__ mts, 
                        const double *__restrict__ d,
                        double *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearchKL64f\0";
    double *G, *V, *phiCerf, *phiWork, sigma2, twopin,
           xdiv, xdiv1, xdiv2, xdiv3, xdiv4, xscal;
    int i, ierr, ir, ldv;
    const double sqrt2 = sqrt(2.0);
    const double sqrt2pi = sqrt(2.0*M_PI);
    //------------------------------------------------------------------------//
    ldv = ldz + (64 - ldz%64)/SIZEOF_DOUBLE;
    V = memory_calloc64f(ldv*npts);
    G = memory_calloc64f(LDG*npts);
    sigma2 = sigma*sigma;
    phiWork = memory_calloc64f(nmt);
    phiCerf = memory_calloc64f(nmt);
    ierr = setG64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
    ierr = array_set64f_work(nmt, 1.0, phi);
    if (nlags == 0)
    {
        // Extract 
        for (ir=0; ir<rank; ir++)
        {
            // Compute outer product
            ierr = linalg_outerProduct64f_work(LAPACK_ROW_MAJOR,
                                               npts, &Z[ldz*ir],
                                               npts, &Z[ldz*ir],
                                               ldv, V);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Failed to compute outerproduct\n",
                        __func__);
            }
            // Compute objective function
            ierr = parmt_mtSearch_mahalanobis64f(ldm, nmt, npts,
                                                 ldv, false,
                                                 G, V,
                                                 mts, d, phiWork);
            xdiv1 = 1.0/(sqrt(2.0*sqrtEig[ir])*sigma);
            xdiv2 =-sqrt2pi/(sqrtEig[ir]*sigma);
            xdiv3 = 1.0/(2.0*sqrtEig[ir]*sigma2);
            xdiv4 = 2.0/sqrt(sqrtEig[ir]);
            for (i=0; i<nmt; i++)
            {
                phiCerf[i] = xdiv1*sqrt(phiWork[i]);
            }
            // 
            vmdErfc(nmt, phiCerf, phiCerf, VML_HA);
            for (i=0; i<nmt; i++)
            {
                phiCerf[i] = (xdiv2*sqrt(phiWork[i]))*phiCerf[i];
            }
            for (i=0; i<nmt; i++)
            {
                phiWork[i] = xdiv4*exp(-xdiv3*phiWork[i]);
            }
            for (i=0; i<nmt; i++)
            {
                phi[i] = phi[i]*(phiCerf[i] + phiWork[i]);
            }
        }
    }
    // Scale the product by 1/(2pi^N) - recall we do not need to account for
    // the (-1)^rank term because of MKL's Ei convention.
    xdiv = 1.0/(sigma*(pow(sqrt2pi, npts)));
    cblas_dscal(nmt, xdiv, phi, 1);
    // release memory
    memory_free64f(&V);
    memory_free64f(&G);
    memory_free64f(&phiCerf);
    memory_free64f(&phiWork);
    return 0;
}

//============================================================================//
int parmt_mtSearchXC64f(const int ldm, const int nmt,
                        const int npts, const int blockSize,
                        const int maxlag, const bool lwantLag,
                        const double *__restrict__ Gxx, 
                        const double *__restrict__ Gyy, 
                        const double *__restrict__ Gzz, 
                        const double *__restrict__ Gxy, 
                        const double *__restrict__ Gxz, 
                        const double *__restrict__ Gyz,
                        const double *__restrict__ mts, 
                        const double *__restrict__ d,
                        double *__restrict__ phi,
                        int *__restrict__ lags)
{
    const char *fcnm = "parmt_mtSearchXC64f\0";
    double GtG[36] __attribute__ ((aligned (64)));
    double Gmag[8] __attribute__ ((aligned (64)));
    double m8[8] __attribute__ ((aligned (64)));
    double *DUmat, *DUmatT, *Gxc, *G, *Umat, *M, *uN, dNorm;
    int i, ic, idx, ierr, imt, jdx, jmt, l1, l2, lag, lcref, lhalf,
        maxxc, mblock, mpts, Mrows, MrowsXC, Ncols, nlag, pad;
    bool ltrunc;
    const double one = 1.0;
    const double zero = 0.0;
    const int Kcols = 6; // number of columns of G
    //------------------------------------------------------------------------//
    lcref = 2*npts - 1; // Cross-correlation size
    pad = computePadding64f(blockSize);
    mblock = blockSize + pad;
    mpts = lcref + computePadding64f(lcref);
    G = memory_calloc64f(LDG*npts);
    Gxc = memory_calloc64f(LDG*lcref);
    M = memory_calloc64f(mblock*8);
    Umat = memory_calloc64f(mblock*npts);
    uN = memory_calloc64f(mblock);
    DUmat = memory_calloc64f(mblock*lcref);
    DUmatT = memory_calloc64f(mpts*blockSize);
#ifdef __INTEL_COMPILER
    __assume_aligned(G, 64);
    __assume_aligned(Gxc, 64);
    __assume_aligned(M, 64);
    __assume_aligned(uN, 64);
    __assume_aligned(Umat, 64);
    __assume_aligned(DUmat, 64);
    __assume_aligned(DUmatT, 64);
#endif
    // Figure out the sizes for the cross-correlation
    MrowsXC = lcref;
    Mrows = npts;
    ltrunc = false;
    lhalf = lcref/2;
    l1 = lhalf - maxlag;
    l2 = lhalf + maxlag;
    nlag = l2 - l1 + 1;
    if (maxlag > 0){ltrunc = true;}
    if (ltrunc)
    {
        if (maxlag >= npts)
        {
            printf("%s: Error maxlag can't exceed npts\n", fcnm);
            return -1;
        }
    }
    // Set the data
    ierr = setG64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
    ierr = setGxc64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, d, Gxc);
    if (ierr != 0)
    {
        printf("%s: Failed to set Gxc\n", fcnm);
        return -1;
    }
    // Compute G^T G
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                6, 6, npts, 1.0, G, LDG,
                G, LDG, 0.0, GtG, 6);
    for (i=0; i<8; i++){m8[i] = zero;}
    memory_free64f(&G);
    // Compute the normalization
    dNorm = cblas_dnrm2(npts, d, 1);
    // Cross-correlate the data
/*
        #pragma omp parallel for \
         firstprivate(csum, M, R) \
         private (i, ic, idx, imt, jdx, jmt, Ncols) \
         shared (d, CG, Dmat, mts, phi, Mrows, mblock) \
         default (none) reduction(+:var)
*/
    for (jmt=0; jmt<nmt; jmt=jmt+blockSize)
    {
        Ncols = MIN(blockSize, nmt - jmt); // Number of columns of M
        // Set the moment tensor matrix
        for (i=0; i<6; i++)
        {
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic;
                idx = ldm*imt + i;
                jdx = mblock*i;
                M[jdx+ic] = mts[idx];
            }
        }
        // Compute the normalization - Eqn A7 Herrmann 2011 
        for (ic=0; ic<Ncols; ic++)
        {
            imt = jmt + ic;
            for (i=0; i<6; i++)
            {
                m8[i] = mts[ldm*imt+i];
            }
            uN[ic] = computemtGTGmt64f(m8, GtG);
            uN[ic] = sqrt(uN[ic]);
        }
        // Compute D*U
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    MrowsXC, Ncols, Kcols, one, Gxc, LDG,
                    M, mblock, zero, DUmat, mblock);
        // Get DUmat into a more convenient column major form
        for (ic=0; ic<Ncols; ic++)
        {
            cblas_dcopy(lcref, &DUmat[ic], mblock, &DUmatT[mpts*ic], 1);
        }
        // Compute the normalized cross-correlations
        if (!ltrunc)
        {
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic;
                // Compute the lag where the times range from [-npts+1:npts-1]
                maxxc = array_argmax64f(lcref, &DUmatT[ic*mpts], &ierr);
                if (lwantLag){lags[imt] =-npts + 1 + maxxc;}
                // Compute the normalized cross-correlation
                phi[imt] = DUmatT[ic*mpts+maxxc]/(dNorm*uN[ic]);
                // Convert from [-1,1] to [0,1]
                phi[imt] = 0.5*(1.0 + phi[imt]);
            }
        }
        // Compute the normalized cross-correlation with max lag
        else
        {
            for (ic=0; ic<Ncols; ic++)
            {    
                imt = jmt + ic;
                // Compute the lag where the times range from 
                // [lhalf-maxlag:half+maxlag]
                maxxc = l1 + array_argmax64f(nlag, &DUmatT[ic*mpts+l1], &ierr);
                if (lwantLag){lags[imt] =-npts + 1 + maxxc;}
                phi[imt] = DUmatT[ic*mpts+maxxc]/(dNorm*uN[ic]);
                // Convert from [-1,1] to [0,1]
                phi[imt] = 0.5*(1.0 + phi[imt]);
            }
        }
    }
    // Free memory
    memory_free64f(&uN);
    memory_free64f(&Umat);
    memory_free64f(&DUmat);
    memory_free64f(&DUmatT);
    memory_free64f(&Gxc);
    memory_free64f(&M);
    return ierr;
}
//============================================================================//
int parmt_mtSearchL164f(const int ldm, const int nmt,
                        const int npts, const int blockSize,
                        const int nlags, const bool lwantLags,
                        const bool lrescale,
                        const double *__restrict__ Gxx, 
                        const double *__restrict__ Gyy,
                        const double *__restrict__ Gzz,
                        const double *__restrict__ Gxy,
                        const double *__restrict__ Gxz,
                        const double *__restrict__ Gyz,
                        const double *__restrict__ CeInv,
                        const double *__restrict__ mts, 
                        const double *__restrict__ d,
                        double *__restrict__ phi,
                        double *__restrict__ var,
                        int *__restrict__ lags)
{
    const char *fcnm = "parmt_mtSearchL164f\0";
    double *CG, *Cd, *G, *Dmat,
           dnorm, objfn, res, viMin;
    int i, ic, idx, jdx, ierr, imt, jmt, k, kmt, mblock, Mrows,
        Ncols, nmtBlocks, nthreads0;
    int pad;
    const double zero = 0.0;
    const double one = 1.0;
    const double negOne =-one;
    const int Kcols = 6; // number of columns of G
    // Set space and make G matrix
    pad = computePadding64f(blockSize);
    mblock = blockSize + pad;
    CG = memory_calloc64f(LDCG*npts);
    G = memory_calloc64f(LDG*npts);
    Cd = memory_calloc64f(npts);
    Dmat = memory_calloc64f(mblock*npts);
    Mrows = npts; // Number of rows of G
/*
#ifdef PARMT_USE_INTEL
    mkl_set_num_threads(1);
#endif
*/
#ifdef __INTEL_COMPILER
    __assume_aligned(G, 64); 
    __assume_aligned(CG, 64);
    __assume_aligned(Dmat, 64);
    __assume_aligned(Cd, 64);
    __assume_aligned(mts, 64);
#endif
    ierr = setG64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
    nmtBlocks = (int) ((double) (nmt)/(double) (blockSize) + 1.0);
    if (nmtBlocks*blockSize < nmt)
    {
        printf("%s: Internal error - all mts wont be visited\n", fcnm);
        return -1;
    }
    // They want lags but there are none to compute
    if (lwantLags && nlags == 0 && lags != NULL)
    {
        array_zeros32i_work(nmt, lags);
    }
    // Introduce the regularizer as to compare ||C_d (Gm - d)||
    computeCG64f(npts, 1, true, CeInv, G, CG);
    computeCd64f(npts, 1, true, CeInv, d, Cd);
    dnorm = array_norm64f(npts, Cd, ONE_NORM, one, &ierr);
    if (!lwantLags || nlags == 0)
    {
        // Set the observations
        for (i=0; i<npts; i++) 
        {
            for (ic=0; ic<blockSize; ic++)
            {
                Dmat[mblock*i+ic] = Cd[i];
            }
        }
        ierr = performL1Search64f(nmt, ldm,
                                  npts,
                                  blockSize, mblock,
                                  Mrows, Kcols,
                                  0,
                                  false, false,
                                  lrescale,
                                  dnorm,
                                  Dmat, CG, Cd, mts, 
                                  lags, phi, var);
    }
    else
    {
/*
        varLoc = array_copy64f(npts, var, &ierr);
        rnorm = memory_calloc64f(mblock);
        unorm = memory_calloc64f(mblock);
        M = memory_calloc64f(mblock*8);
        R = memory_calloc64f(mblock*npts);
        U = memory_calloc64f(mblock*npts);
*/
        array_set64f_work(nmt, zero, phi);
        for (k=-nlags; k<=nlags; k++)
        {
            // Null out the the data matrix
            memset(Dmat, 0, (size_t) (mblock*npts)*sizeof(double));
            // Time advance in Greens function emulated by data delay 
            if (k < 0)
            {
                for (i=-k; i<npts; i++)
                {
                    for (ic=0; ic<blockSize; ic++)
                    {
                        Dmat[mblock*i+ic] = Cd[i+k];
                    }
                }
            }
            // Time delay in Greens function emulated by data advance
            else if (k > 0)
            {
                for (i=k; i<npts; i++)
                {
                    for (ic=0; ic<blockSize; ic++)
                    {
                        Dmat[mblock*(i-k)+ic] = Cd[i];
                    }
                }
            }
            // Straight copy
            else
            {
                for (i=0; i<npts; i++)
                {
                    for (ic=0; ic<blockSize; ic++)
                    {
                        Dmat[mblock*i+ic] = Cd[i];
                    }
                }
            }
            // Optimize at this lag
            ierr = performL1Search64f(nmt, ldm,
                                      npts,
                                      blockSize, mblock,
                                      Mrows, Kcols,
                                      k,
                                      true, lwantLags,
                                      lrescale,
                                      dnorm,
                                      Dmat, CG, Cd, mts,
                                      lags, phi, var);
/*
            // Loop on moment tensor block
            for (kmt=0; kmt<nmtBlocks; kmt++)
            {
                jmt = kmt*blockSize;
                Ncols = MIN(blockSize, nmt - jmt); // Number of columns of M
                // Set the moment tensor matrix
                for (i=0; i<6; i++) 
                {
                    #pragma omp simd
                    for (ic=0; ic<Ncols; ic++)
                    {
                        imt = jmt + ic;
                        idx = ldm*imt + i;
                        jdx = mblock*i;
                        M[jdx+ic] = mts[idx];
                    }
                }
                // Compute U = GM
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            Mrows, Ncols, Kcols, one, CG, LDCG,
                            M, mblock, zero, U, mblock);
                memset(rnorm, 0, (size_t) Ncols*sizeof(double));
                memset(unorm, 0, (size_t) Ncols*sizeof(double));
                #pragma omp simd
                for (i=0; i<npts*Ncols; i++)
                {
                    R[i] = fabs(Dmat[i] - U[i]);
                }
                // Compute |d - u|_1
                for (i=0; i<npts; i++)
                {
                    viMin = varLoc[i];
                    #pragma omp simd reduction(min:viMin)
                    for (ic=0; ic<Ncols; ic++)
                    {
                        rnorm[ic] = rnorm[ic] + R[mblock*i+ic];
                        viMin = fmin(viMin, R[mblock*i+ic]);
                    }
                    varLoc[i] = viMin;
                }
                // Compute |u|_1
                for (i=0; i<npts; i++)
                {
                    #pragma omp simd
                    for (ic=0; ic<Ncols; ic++)
                    {
                        unorm[ic] = unorm[ic] + fabs(U[mblock*i+ic]);
                    }
                }
                // Apply the triangle inequality
                for (ic=0; ic<Ncols; ic++)
                {
                    rnorm[ic] = 1.0 - rnorm[ic]/(dnorm + unorm[ic]);
                } 
                // Compute the L1 norm and take the minimum
                if (!lwantLags)
                {
                    for (ic=0; ic<Ncols; ic++)
                    {    
                        imt = jmt + ic;
                        if (phi[imt] > rnorm[ic])
                        {
                            phi[imt] = rnorm[ic];
                            lags[imt] = k;
                        }
                    }
                }
                else
                {
                    for (ic=0; ic<Ncols; ic++)
                    {
                        imt = jmt + ic;
                        phi[imt] = fmax(phi[imt], rnorm[ic]);
                    }
                }
            } // Loop on moment tensor block
*/
        } // Loop on lags
/*
        memory_free64f(&varLoc);
        memory_free64f(&rnorm);
        memory_free64f(&unorm);
        memory_free64f(&M);
        memory_free64f(&R);
        memory_free64f(&U);
*/
    }
    memory_free64f(&CG);
    memory_free64f(&G);
    memory_free64f(&Dmat);
    return ierr;
}
 
//============================================================================//
int parmt_mtSearch64f(const int ldm, const int nlags, const int nmt,
                      const int npts, const int ldc, 
                      const bool ldiag,
                      const double *__restrict__ Gxx, 
                      const double *__restrict__ Gyy, 
                      const double *__restrict__ Gzz, 
                      const double *__restrict__ Gxy, 
                      const double *__restrict__ Gxz, 
                      const double *__restrict__ Gyz, 
                      const double *__restrict__ CeInv,
                      const double *__restrict__ mts, 
                      const double *__restrict__ d,
                      double *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearch32f\0";
    double *G;
    int ierr, ir;
    if (nmt < 1 || npts < 1 || mts == NULL || Gxx == NULL || Gyy == NULL ||
        Gzz == NULL || Gxy == NULL || Gxz == NULL || Gxy == NULL ||
        CeInv == NULL || d == NULL || phi == NULL)
    {
        if (nmt < 1){printf("%s: No moment tensors\n", fcnm);}
        if (npts < 1){printf("%s: No data points\n", fcnm);}
        if (Gxx == NULL){printf("%s: Gxx is NULL\n", fcnm);}
        if (Gyy == NULL){printf("%s: Gyy is NULL\n", fcnm);}
        if (Gzz == NULL){printf("%s: Gzz is NULL\n", fcnm);}
        if (Gxy == NULL){printf("%s: Gxy is NULL\n", fcnm);}
        if (Gxz == NULL){printf("%s: Gxz is NULL\n", fcnm);}
        if (CeInv == NULL){printf("%s: CeInv is NULL\n", fcnm);}
        if (mts == NULL){printf("%s: mts is NULL\n", fcnm);}
        if (d == NULL){printf("%s: d is NULL\n", fcnm);}
        if (phi == NULL){printf("%s: phi is NULL\n", fcnm);}
        return -1;
    }
    if (!ldiag && ldc < npts)
    {
        printf("%s: Error ldc must be >= npts\n", fcnm);
        return -1;
    }
    if (nlags > 0 && nlags > npts)
    {
        printf("%s: Error nlags can't exceed npts\n", fcnm);
        return -1;
    }
    // Set space and make G matrix
    G = memory_calloc64f(LDG*npts);
    ierr = setG64f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
    // straight application of Ce^{-1}
    if (nlags == 0)
    {
        // Diagonal
        if (ldiag)
        {
            ierr = parmt_mtSearch_mahalanobis64f(ldm, nmt, npts,
                                                 ldc, true,
                                                 G, CeInv,
                                                 mts, d, phi);
        }
        // General matrix
        else
        {
            ierr = parmt_mtSearch_mahalanobis64f(ldm, nmt, npts,
                                                 ldc, false,
                                                 G, CeInv,
                                                 mts, d, phi);
        }
    }
    else
    {
        // Diagonal weighting matrix
        if (ldiag)
        {
            ierr = parmt_mtSearch_shiftedMahalanobis64f(nlags, nmt, 
                                                        npts, ldc, 
                                                        true,
                                                        G, CeInv,
                                                        mts, d, phi);
        }
        else
        {
            ierr = parmt_mtSearch_shiftedMahalanobis64f(nlags, nmt, 
                                                        npts, ldc, 
                                                        false,
                                                        G, CeInv,
                                                        mts, d, phi);
        }
    }
    if (ierr != 0)
    {
        printf("%s: Error computing search\n", fcnm);
    }
    // Free memory
    memory_free64f(&G);
    return ierr;
}
//============================================================================//
int parmt_mtSearch32f(const int ldm, const int nlags, const int nmt,
                      const int npts, const int ldc,
                      const bool ldiag,
                      const float *__restrict__ Gxx,
                      const float *__restrict__ Gyy, 
                      const float *__restrict__ Gzz,
                      const float *__restrict__ Gxy,
                      const float *__restrict__ Gxz,
                      const float *__restrict__ Gyz,
                      const float *__restrict__ CeInv,
                      const float *__restrict__ mts,
                      const float *__restrict__ d,
                      float *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearch32f\0";
    float *G;
    int ierr, ir;
    if (nmt < 1 || npts < 1 || mts == NULL || Gxx == NULL || Gyy == NULL ||
        Gzz == NULL || Gxy == NULL || Gxz == NULL || Gxy == NULL ||
        CeInv == NULL || d == NULL || phi == NULL)
    {
        if (nmt < 1){printf("%s: No moment tensors\n", fcnm);}
        if (npts < 1){printf("%s: No data points\n", fcnm);}
        if (Gxx == NULL){printf("%s: Gxx is NULL\n", fcnm);}
        if (Gyy == NULL){printf("%s: Gyy is NULL\n", fcnm);}
        if (Gzz == NULL){printf("%s: Gzz is NULL\n", fcnm);}
        if (Gxy == NULL){printf("%s: Gxy is NULL\n", fcnm);}
        if (Gxz == NULL){printf("%s: Gxz is NULL\n", fcnm);}
        if (CeInv == NULL){printf("%s: CeInv is NULL\n", fcnm);}
        if (mts == NULL){printf("%s: mts is NULL\n", fcnm);}
        if (d == NULL){printf("%s: d is NULL\n", fcnm);}
        if (phi == NULL){printf("%s: phi is NULL\n", fcnm);}
        return -1;
    }
    if (!ldiag && ldc < npts)
    {
        printf("%s: Error ldc must be >= npts\n", fcnm);
        return -1;
    }
    if (nlags > 0 && nlags > npts)
    {
        printf("%s: Error nlags can't exceed npts\n", fcnm);
        return -1;
    }
    // Set space and make G matrix
    G = memory_calloc32f(LDG*npts);
    ierr = setG32f(npts, Gxx, Gyy, Gzz, Gxy, Gxz, Gyz, G);
    if (ierr != 0)
    {
        printf("%s: Failed to set G\n", fcnm);
        return -1;
    }
    // straight application of Ce^{-1}
    if (nlags == 0)
    {
        if (ldiag)
        {
            ierr = parmt_mtSearch_mahalanobis32f(ldm, nmt, npts,
                                                 ldc, true,
                                                 G, CeInv,
                                                 mts, d, phi);
        }
        else
        {
            ierr = parmt_mtSearch_mahalanobis32f(ldm, nmt, npts,
                                                 ldc, false,
                                                 G, CeInv,
                                                 mts, d, phi);
        }
    }
    // shifting the data to the green's functions
    else
    {
        // Diagonal weighting matrix
        if (ldiag)
        {
            ierr = parmt_mtSearch_shiftedMahalanobis32f(nlags, nmt,
                                                        npts, ldc,
                                                        true,
                                                        G, CeInv,
                                                        mts, d, phi);
        }
        // Full inverse covariance matrix
        else
        {
            ierr = parmt_mtSearch_shiftedMahalanobis32f(nlags, nmt,
                                                        npts, ldc,
                                                        false,
                                                        G, CeInv,
                                                        mts, d, phi);
        }
    }
    if (ierr != 0)
    {
        printf("%s: Error computing search\n", fcnm);
    }
    // Free memory
    memory_free32f(&G);
    return ierr;
}
//============================================================================//
int parmt_mtSearch_mahalanobis64f(const int ldm, const int nmt, const int npts,
                                  const int ldc, const bool ldiag,
                                  const double *__restrict__ G,
                                  const double *__restrict__ CeInv,
                                  const double *__restrict__ mts,
                                  const double *__restrict__ d,
                                  double *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearch_mahalanobis64f\0";
    double *Cd, *CG;
    double GtCG[36], GtCd[6], m[6] __attribute__ ((aligned (64)));
    double objfn __attribute__ ((aligned (64))) = 0.0;
    double dtCd __attribute__ ((aligned (64))) = 0.0;
    int imt, M;
    //------------------------------------------------------------------------// 
    //
    // Set space
    CG = memory_calloc64f(LDCG*npts);
    Cd = memory_calloc64f(npts);
    // Compute C G
    computeCG64f(npts, ldc, ldiag, CeInv, G, CG);
    // Compute G^T C G
    computeGtCG64f(npts, G, CG, GtCG);
    // Compute C d
    computeCd64f(npts, ldc, ldiag, CeInv, d, Cd);
    // Compute d^T C G = (G^T C d)^T
    M = npts;
    cblas_dgemv(CblasRowMajor, CblasTrans, M, 6, 1.0, CG, LDCG,
                d, 1, 0.0, GtCd, 1);
    // Compute d^T C d
    dtCd = cblas_ddot(npts, d, 1, Cd, 1);
#ifdef __INTEL_COMPILER
    __assume_aligned(mts, 64);
    __assume_aligned(phi, 64);
    #pragma omp parallel for simd \
     firstprivate (GtCG, GtCd, dtCd) \
     private(m) shared(mts, phi) \
     default(none)
#else
/*
    #pragma omp parallel for simd \
     firstprivate (GtCG, GtCd, dtCd) \
     private(m) shared(mts, phi) \
     default(none) aligned(phi, mts: 64)
*/
#endif
    for (imt=0; imt<nmt; imt++)
    {
        m[0] = mts[ldm*imt];
        m[1] = mts[ldm*imt+1];
        m[2] = mts[ldm*imt+2];
        m[3] = mts[ldm*imt+3];
        m[4] = mts[ldm*imt+4];
        m[5] = mts[ldm*imt+5];
        phi[imt] = computeLeastSquaresResidual64f(m, GtCG, GtCd, dtCd);
    }
    // free space
    memory_free64f(&CG);
    memory_free64f(&Cd);
    return 0;
}

int parmt_mtSearch_mahalanobis32f(const int ldm, const int nmt, const int npts,
                                  const int ldc, const bool ldiag,
                                  const float *__restrict__ G,
                                  const float *__restrict__ CeInv,
                                  const float *__restrict__ mts,
                                  const float *__restrict__ d,
                                  float *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearch_mahalanobis32f\0";
    float *Cd, *CG;
    float GtCG[36], GtCd[6], m[6] __attribute__ ((aligned (64)));
    float objfn __attribute__ ((aligned (64))) = 0.0f;
    float dtCd __attribute__ ((aligned (64))) = 0.0f;
    int imt, M;
    //------------------------------------------------------------------------// 
    //
    // Set space
    CG = memory_calloc32f(LDCG*npts);
    Cd = memory_calloc32f(npts);
    // Compute C G
    computeCG32f(npts, ldc, ldiag, CeInv, G, CG);
    // Compute G^T C G
    computeGtCG32f(npts, G, CG, GtCG);
    // Compute C d
    computeCd32f(npts, ldc, ldiag, CeInv, d, Cd);
    // Compute d^T C G = (G^T C d)^T
    M = npts;
    cblas_sgemv(CblasRowMajor, CblasTrans, M, 6, 1.0f, CG, LDCG,
                d, 1, 0.0f, GtCd, 1);
    // Compute d^T C d
    dtCd = cblas_sdot(npts, d, 1, Cd, 1);
#ifdef __INTEL_COMPILER
    __assume_aligned(mts, 64);
    __assume_aligned(phi, 64);
    #pragma omp parallel for simd \
     firstprivate (GtCG, GtCd, dtCd) \
     private(m) shared(mts, phi) \
     default(none)
#else
/*
    #pragma omp parallel for simd \
     firstprivate (GtCG, GtCd, dtCd) \
     private(m) shared(mts, phi) \
     default(none) aligned(phi, mts: 64)
*/
#endif
    for (imt=0; imt<nmt; imt++)
    {
        m[0] = mts[ldm*imt];
        m[1] = mts[ldm*imt+1];
        m[2] = mts[ldm*imt+2];
        m[3] = mts[ldm*imt+3];
        m[4] = mts[ldm*imt+4];
        m[5] = mts[ldm*imt+5];
        phi[imt] = computeLeastSquaresResidual32f(m, GtCG, GtCd, dtCd);
    }
    // free space
    memory_free32f(&CG);
    memory_free32f(&Cd);
    return 0;
}
//============================================================================//
int parmt_mtSearch_shiftedMahalanobis64f(const int nlags, const int nmt, 
                                         const int npts, const int ldc, 
                                         const bool ldiag,
                                         const double *__restrict__ G,
                                         const double *__restrict__ CeInv,
                                         const double *__restrict__ mts, 
                                         const double *__restrict__ d,
                                         double *__restrict__ phi) 
{
    const char *fcnm = "parmt_mtSearch_shiftedMahalanobis64f\0";
    int ncopy;
    double *Cd, *CG, *dwork, *phiLoc;
    double GtCG[36], GtCd[6], m[6] __attribute__ ((aligned (64)));
    double objfn __attribute__ ((aligned (64))) = 0.0;
    double dtCd __attribute__ ((aligned (64))) = 0.0;
    int i, imt, k, M;
    bool linit;
    //------------------------------------------------------------------------//
    //
    // Try to avoid some problems
    if (nmt < 1 || npts < 1)
    {
        printf("%s: Invalid number of points %d %d\n", fcnm, nmt, npts);
        return -1;
    }
    if (G == NULL || CeInv == NULL || mts == NULL || d == NULL || phi == NULL)
    {
        printf("%s: Null pointer detected\n", fcnm);
        return -1;
    }
    if (memory_isAligned(CeInv, 64) != 1 || memory_isAligned(mts, 64) != 1)
    {
        printf("%s: CeInv and mts must be 64 bit aligned\n", fcnm);
        return -1;
    }
    // Set space
    CG = memory_calloc64f(LDCG*npts);
    Cd = memory_calloc64f(npts);
    dwork = memory_calloc64f(npts);
    phiLoc = memory_calloc64f(nmt);
    // Compute C G
    computeCG64f(npts, ldc, ldiag, CeInv, G, CG);
    // Compute G^T C G
    computeGtCG64f(npts, G, CG, GtCG);
#ifdef __INTEL_COMPILER
    __assume_aligned(mts, 64);
    __assume_aligned(phiLoc, 64);
#endif
    // Loop on lags
    linit = false;
    #pragma omp parallel \
     firstprivate (dwork, Cd, dtCd, GtCd, GtCG, linit, phiLoc) \
     private (i, imt, k, m, M, ncopy, objfn) \
     shared (CeInv, CG, d, mts, phi) \
     default (none)
    {
    for (k=-nlags; k<=nlags; k++)
    {
        // shift greens fns left by shifting data right 
        if (k < 0)
        {
            ncopy = npts + k;
            for (i=0; i<-k; i++){dwork[i] = 0.0;}
            cblas_dcopy(ncopy, d, 1, &dwork[-k], 1);
        }
        else
        {
            // shift greens fns right by shifting data left
            if (k > 0)
            {
                ncopy = npts - k;
                for (i=ncopy; i<npts; i++){dwork[i] = 0.0;}
                cblas_dcopy(ncopy, &d[k], 1, dwork, 1);
            }
            else
            {
                cblas_dcopy(npts, d, 1, dwork, 1);
            }
        }
        // Compute C d
        computeCd64f(npts, ldc, ldiag, CeInv, dwork, Cd);
        // Compute d^T C G = (G^T C d)^T
        M = npts;
        cblas_dgemv(CblasRowMajor, CblasTrans, M, 6, 1.0, CG, LDCG,
                    dwork, 1, 0.0, GtCd, 1);
        // Compute d^T C d
        dtCd = cblas_ddot(npts, dwork, 1, Cd, 1);
        if (linit)
        {
            for (imt=0; imt<nmt; imt++)
            {
                m[0] = mts[8*imt];
                m[1] = mts[8*imt+1];
                m[2] = mts[8*imt+2];
                m[3] = mts[8*imt+3];
                m[4] = mts[8*imt+4];
                m[5] = mts[8*imt+5];
                objfn = computeLeastSquaresResidual64f(m, GtCG, GtCd, dtCd);
                phiLoc[imt] = fmin(phiLoc[imt], objfn);
            }
        }
        else
        {
            for (imt=0; imt<nmt; imt++)
            {    
                m[0] = mts[8*imt];
                m[1] = mts[8*imt+1];
                m[2] = mts[8*imt+2];
                m[3] = mts[8*imt+3];
                m[4] = mts[8*imt+4];
                m[5] = mts[8*imt+5];
                objfn = computeLeastSquaresResidual64f(m, GtCG, GtCd, dtCd);
                phiLoc[imt] = objfn;
            }
            linit = true;
        }
    }
    #pragma omp barrier
    // Reduction
    for (imt=0; imt<nmt; imt++)
    {
        phi[imt] = phiLoc[imt];
    }
    // End parallel section
    }
    // free space
    memory_free64f(&CG);
    memory_free64f(&Cd);
    memory_free64f(&dwork);
    memory_free64f(&phiLoc);
    return 0;
}


int parmt_mtSearch_shiftedMahalanobis32f(const int nlags, const int nmt,
                                         const int npts, const int ldc,
                                         const bool ldiag,
                                         const float *__restrict__ G,
                                         const float *__restrict__ CeInv,
                                         const float *__restrict__ mts,
                                         const float *__restrict__ d,
                                         float *__restrict__ phi)
{
    const char *fcnm = "parmt_mtSearch_shiftedMahalanobis32f\0";
    int ncopy;
    float *Cd, *CG, *dwork, *phiLoc;
    float GtCG[36], GtCd[6], m[6] __attribute__ ((aligned (64)));
    float objfn __attribute__ ((aligned (64))) = 0.0f;
    float dtCd __attribute__ ((aligned (64))) = 0.0f;
    int i, imt, k, M;
    bool linit;
    //------------------------------------------------------------------------//
    //
    // Try to avoid some problems
    if (nmt < 1 || npts < 1)
    {
        printf("%s: Invalid number of points %d %d\n", fcnm, nmt, npts);
        return -1;
    }
    if (G == NULL || CeInv == NULL || mts == NULL || d == NULL || phi == NULL)
    {
        printf("%s: Null pointer detected\n", fcnm);
        return -1;
    }
    if (memory_isAligned(CeInv, 64) != 1 || memory_isAligned(mts, 64) != 1)
    {
        printf("%s: CeInv and mts must be 64 bit aligned\n", fcnm);
        return -1;
    }
    // Set space
    CG = memory_calloc32f(LDCG*npts);
    Cd = memory_calloc32f(npts);
    dwork = memory_calloc32f(npts);
    phiLoc = memory_calloc32f(nmt);
    // Compute C G
    computeCG32f(npts, ldc, ldiag, CeInv, G, CG);
    // Compute G^T C G
    computeGtCG32f(npts, G, CG, GtCG);
#ifdef __INTEL_COMPILER
    __assume_aligned(mts, 64);
    __assume_aligned(phiLoc, 64);
#endif
    // Loop on lags
    linit = false;
    #pragma omp parallel \
     firstprivate (dwork, Cd, dtCd, linit, GtCd, GtCG, phiLoc) \
     private (i, imt, k, m, M, ncopy, objfn) \
     shared (CeInv, CG, d, mts, phi) \
     default (none)
    {
    for (k=-nlags; k<=nlags; k++)
    {
        // shift greens fns left by shifting data right 
        if (k < 0)
        {
            ncopy = npts + k;
            for (i=0; i<-k; i++){dwork[i] = 0.0f;}
            cblas_scopy(ncopy, d, 1, &dwork[-k], 1);
        }
        else
        {
            // shift greens fns right by shifting data left
            if (k > 0)
            {
                ncopy = npts - k;
                for (i=ncopy; i<npts; i++){dwork[i] = 0.0f;}
                cblas_scopy(ncopy, &d[k], 1, dwork, 1);  
            }
            else
            {
                cblas_scopy(npts, d, 1, dwork, 1);
            }
        }
        // Compute C d
        computeCd32f(npts, ldc, ldiag, CeInv, dwork, Cd);
        // Compute d^T C G = (G^T C d)^T
        M = npts;
        cblas_sgemv(CblasRowMajor, CblasTrans, M, 6, 1.0f, CG, LDCG,
                    dwork, 1, 0.0f, GtCd, 1);
        // Compute d^T C d
        dtCd = cblas_sdot(npts, dwork, 1, Cd, 1);
        if (linit)
        {
            for (imt=0; imt<nmt; imt++)
            {
                m[0] = mts[8*imt];
                m[1] = mts[8*imt+1];
                m[2] = mts[8*imt+2];
                m[3] = mts[8*imt+3];
                m[4] = mts[8*imt+4];
                m[5] = mts[8*imt+5];
                objfn = computeLeastSquaresResidual32f(m, GtCG, GtCd, dtCd);
                phiLoc[imt] = fminf(phiLoc[imt], objfn);
            }
        }
        else
        {
            for (imt=0; imt<nmt; imt++)
            {
                m[0] = mts[8*imt];
                m[1] = mts[8*imt+1];
                m[2] = mts[8*imt+2];
                m[3] = mts[8*imt+3];
                m[4] = mts[8*imt+4];
                m[5] = mts[8*imt+5];
                objfn = computeLeastSquaresResidual32f(m, GtCG, GtCd, dtCd);
                phiLoc[imt] = objfn;
            }
            linit = true;
        }
    }
    #pragma omp barrier
    // Reduction
    for (imt=0; imt<nmt; imt++)
    {
        phi[imt] = phiLoc[imt];
    }
    // End parallel section
    }
    // free space
    memory_free32f(&CG);
    memory_free32f(&Cd);
    memory_free32f(&dwork);
    memory_free32f(&phiLoc);
    return 0;
}
//============================================================================//
//                                PRIVATE FUNCTIONS                           //
//============================================================================//
/*!
 * @brief Sets the single precision Green's functions matrix
 *
 * @param[in] npts      number of points in Green's functions.
 * @param[in] Gxx       64 bit aligned Green's functions which are scaled by mxx
 *                      moment tensor term [npts]
 * @param[in] Gyy       64 bit aligned Green's functions which are scaled by myy
 *                      moment tensor term [npts]
 * @param[in] Gzz       64 bit aligned Green's functions which are scaled by mzz
 *                      moment tensor term [npts]
 * @param[in] Gxy       64 bit aligned Green's functions which are scaled by mxy
 *                      moment tensor term [npts]
 * @param[in] Gxz       64 bit aligned Green's functions which are scaled by mxz
 *                      moment tensor term [npts]
 * @param[in] Gyz       64 bit aligned Green's functions which are scaled by myz
 *                      moment tensor term [npts
 *
 * @param[out] G        64 bit aligned Green's functions matrix in row major
 *                      format that multiplies moment tensor array
 *                      $\{ m_{xx}, m_{yy}, m_{zz}, m_{xy}, m_{xz}, m_{yz} \}^T$
 *                      [npts x LDG]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int setG32f(const int npts,
                   const float *__restrict__ Gxx,
                   const float *__restrict__ Gyy,
                   const float *__restrict__ Gzz, 
                   const float *__restrict__ Gxy, 
                   const float *__restrict__ Gxz, 
                   const float *__restrict__ Gyz, 
                   float *__restrict__ G)
{
    const char *fcnm = "setG32f\0";
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
//============================================================================//
static int setGxc64f(const int npts,
                     const double *__restrict__ Gxx,
                     const double *__restrict__ Gyy,
                     const double *__restrict__ Gzz,
                     const double *__restrict__ Gxy,
                     const double *__restrict__ Gxz,
                     const double *__restrict__ Gyz,
                     const double *__restrict__ d,
                     double *__restrict__ Gxc)
{
    const char *fcnm = "setGxc64f\0";
    double *gxxWork, *gyyWork, *gzzWork, *gxyWork, *gxzWork, *gyzWork;
    int ierr, lcref;
    // Set the data
    gxxWork = signal_convolve_correlate64f(npts, d, npts, Gxx,
                                           CONVCOR_FULL, &lcref, &ierr);
    gyyWork = signal_convolve_correlate64f(npts, d, npts, Gyy,
                                           CONVCOR_FULL, &lcref, &ierr);
    gzzWork = signal_convolve_correlate64f(npts, d, npts, Gzz,
                                           CONVCOR_FULL, &lcref, &ierr);
    gxyWork = signal_convolve_correlate64f(npts, d, npts, Gxy,
                                           CONVCOR_FULL, &lcref, &ierr);
    gxzWork = signal_convolve_correlate64f(npts, d, npts, Gxz,
                                           CONVCOR_FULL, &lcref, &ierr);
    gyzWork = signal_convolve_correlate64f(npts, d, npts, Gyz,
                                           CONVCOR_FULL, &lcref, &ierr);
    if (lcref != 2*npts - 1)
    {
        fprintf(stderr, "%s: Error invalid number of points\n", __func__);
        ierr = 1;
    }
    cblas_dcopy(lcref, gxxWork, 1, &Gxc[0], LDG);
    cblas_dcopy(lcref, gyyWork, 1, &Gxc[1], LDG);
    cblas_dcopy(lcref, gzzWork, 1, &Gxc[2], LDG);
    cblas_dcopy(lcref, gxyWork, 1, &Gxc[3], LDG);
    cblas_dcopy(lcref, gxzWork, 1, &Gxc[4], LDG);
    cblas_dcopy(lcref, gyzWork, 1, &Gxc[5], LDG);
    memory_free64f(&gxxWork);
    memory_free64f(&gyyWork);
    memory_free64f(&gzzWork);
    memory_free64f(&gxyWork);
    memory_free64f(&gxzWork);
    memory_free64f(&gyzWork);
    return ierr;
}
/*!
 * @brief Sets the double precision Green's functions matrix
 *
 * @param[in] npts      number of points in Green's functions.
 * @param[in] Gxx       64 bit aligned Green's functions which are scaled by mxx
 *                      moment tensor term [npts]
 * @param[in] Gyy       64 bit aligned Green's functions which are scaled by myy
 *                      moment tensor term [npts]
 * @param[in] Gzz       64 bit aligned Green's functions which are scaled by mzz
 *                      moment tensor term [npts]
 * @param[in] Gxy       64 bit aligned Green's functions which are scaled by mxy
 *                      moment tensor term [npts]
 * @param[in] Gxz       64 bit aligned Green's functions which are scaled by mxz
 *                      moment tensor term [npts]
 * @param[in] Gyz       64 bit aligned Green's functions which are scaled by myz
 *                      moment tensor term [npts
 *
 * @param[out] G        64 bit aligned Green's functions matrix in row major
 *                      format that multiplies moment tensor array
 *                      $\{ m_{xx}, m_{yy}, m_{zz}, m_{xy}, m_{xz}, m_{yz} \}^T$
 *                      [npts x LDG]
 *
 * @result 0 indicates success
 *
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
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
//============================================================================//
/*!
 * @brief Computes the squared residual error for the given moment tensor
 *        by expanding
 *        \f[
 *          (G \textbf{m} - \textbf{d})^T C_e^{-1} (G \textbf{m} - \textbf{d} )
 *         = G^T C_e^{-1} G^T
 *         - 2 m^T G^T C \textbf{d}
 *         + \textbf{d}^T C_{e}^{-1} \textbf{d}
 *        \f]
 *
 * @param[in] m     moment tensor terms [6]
 * @param[in] GtCG  normal equation \f$ G^T C_{e}^{-1} G \f$ [8 x 8].
 * @param[in] GtCd  the \f$ G^T C_{e}^{-1} \textbf{d} \f$ term [6 x 1]
 * @param[in] dtCd  the \f$ \textbf{d}^T C_{e}^{-1} \texbf{d} \f$ term
 *
 * @result squared residual error
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
#pragma omp declare simd
static float computeLeastSquaresResidual32f(const float *__restrict__ m,
                                            const float *__restrict__ GtCG,
                                            const float *__restrict__ GtCd,
                                            const float dtCd)
{
    const float __attribute__ ((aligned (64))) mtwo =-2.0f;
    float l2 __attribute__ ((aligned (64))) = 0.0f;
    float mtGtd, mtGtG, two_mtGtd_dtd __attribute__ ((aligned (64)));
    // m^T G^T C d
    mtGtd = m[0]*GtCd[0] + m[1]*GtCd[1] + m[2]*GtCd[2]
          + m[3]*GtCd[3] + m[4]*GtCd[4] + m[5]*GtCd[5];
    // -2 m^T G^T C d + d^T C d
    two_mtGtd_dtd = mtwo*mtGtd + dtCd;
    // m^T G^T C G m 
    mtGtG = m[0]*(m[0]*GtCG[0]  + m[1]*GtCG[1]  + m[2]*GtCG[2]
                + m[3]*GtCG[3]  + m[4]*GtCG[4]  + m[5]*GtCG[5])
          + m[1]*(m[0]*GtCG[6]  + m[1]*GtCG[7]  + m[2]*GtCG[8]
                + m[3]*GtCG[9]  + m[4]*GtCG[10] + m[5]*GtCG[11])
          + m[2]*(m[0]*GtCG[12] + m[1]*GtCG[13] + m[2]*GtCG[14]
                + m[3]*GtCG[15] + m[4]*GtCG[16] + m[5]*GtCG[17])
          + m[3]*(m[0]*GtCG[18] + m[1]*GtCG[19] + m[2]*GtCG[20]
                + m[3]*GtCG[21] + m[4]*GtCG[22] + m[5]*GtCG[23])
          + m[4]*(m[0]*GtCG[24] + m[1]*GtCG[25] + m[2]*GtCG[26]
                + m[3]*GtCG[27] + m[4]*GtCG[28] + m[5]*GtCG[29])
          + m[5]*(m[0]*GtCG[30] + m[1]*GtCG[31] + m[2]*GtCG[32]
                + m[3]*GtCG[33] + m[4]*GtCG[34] + m[5]*GtCG[35]);
    l2 = MAX(0.0f, mtGtG + two_mtGtd_dtd);
    return l2;
}
//============================================================================//
/*!
 * @brief Computes the squared residual error for the given moment tensor
 *        by expanding
 *        \f[
 *          (G \textbf{m} - \textbf{d})^T C_e^{-1} (G \textbf{m} - \textbf{d} )
 *         = G^T C_e^{-1} G^T
 *         - 2 m^T G^T C \textbf{d}
 *         + \textbf{d}^T C_{e}^{-1} \textbf{d}
 *        \f]
 *
 * @param[in] m     moment tensor terms [6]
 * @param[in] GtCG  normal equation \f$ G^T C_{e}^{-1} G \f$ [8 x 8].
 * @param[in] GtCd  the \f$ G^T C_{e}^{-1} \textbf{d} \f$ term [6 x 1]
 * @param[in] dtCd  the \f$ \textbf{d}^T C_{e}^{-1} \texbf{d} \f$ term
 *
 * @result squared residual error
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
#pragma omp declare simd
static double computeLeastSquaresResidual64f(const double *__restrict__ m,
                                             const double *__restrict__ GtCG,
                                             const double *__restrict__ GtCd,
                                             const double dtCd)
{
    const double __attribute__ ((aligned (64))) mtwo =-2.0;
    double l2 __attribute__ ((aligned (64))) = 0.0;
    double mtGtd, mtGtG, two_mtGtd_dtd __attribute__ ((aligned (64)));
    // m^T G^T C d
    mtGtd = m[0]*GtCd[0] + m[1]*GtCd[1] + m[2]*GtCd[2]
          + m[3]*GtCd[3] + m[4]*GtCd[4] + m[5]*GtCd[5];
    // -2 m^T G^T C d + d^T C d
    two_mtGtd_dtd = mtwo*mtGtd + dtCd;
    // m^T G^T C G m
    mtGtG = m[0]*(m[0]*GtCG[0]  + m[1]*GtCG[1]  + m[2]*GtCG[2]
                + m[3]*GtCG[3]  + m[4]*GtCG[4]  + m[5]*GtCG[5])
          + m[1]*(m[0]*GtCG[6]  + m[1]*GtCG[7]  + m[2]*GtCG[8]
                + m[3]*GtCG[9]  + m[4]*GtCG[10] + m[5]*GtCG[11])
          + m[2]*(m[0]*GtCG[12] + m[1]*GtCG[13] + m[2]*GtCG[14]
                + m[3]*GtCG[15] + m[4]*GtCG[16] + m[5]*GtCG[17])
          + m[3]*(m[0]*GtCG[18] + m[1]*GtCG[19] + m[2]*GtCG[20]
                + m[3]*GtCG[21] + m[4]*GtCG[22] + m[5]*GtCG[23])
          + m[4]*(m[0]*GtCG[24] + m[1]*GtCG[25] + m[2]*GtCG[26]
                + m[3]*GtCG[27] + m[4]*GtCG[28] + m[5]*GtCG[29])
          + m[5]*(m[0]*GtCG[30] + m[1]*GtCG[31] + m[2]*GtCG[32]
                + m[3]*GtCG[33] + m[4]*GtCG[34] + m[5]*GtCG[35]);
    l2 = MAX(0.0, mtGtG + two_mtGtd_dtd);
    return l2;
}
//============================================================================//
static double computemtGTGmt64f(const double *__restrict__ m,
                                const double *__restrict__ GtG)
{
    double mtGtG;
    // m^T G^T C G m
    mtGtG = m[0]*(m[0]*GtG[0]  + m[1]*GtG[1]  + m[2]*GtG[2]
                + m[3]*GtG[3]  + m[4]*GtG[4]  + m[5]*GtG[5])
          + m[1]*(m[0]*GtG[6]  + m[1]*GtG[7]  + m[2]*GtG[8]
                + m[3]*GtG[9]  + m[4]*GtG[10] + m[5]*GtG[11])
          + m[2]*(m[0]*GtG[12] + m[1]*GtG[13] + m[2]*GtG[14]
                + m[3]*GtG[15] + m[4]*GtG[16] + m[5]*GtG[17])
          + m[3]*(m[0]*GtG[18] + m[1]*GtG[19] + m[2]*GtG[20]
                + m[3]*GtG[21] + m[4]*GtG[22] + m[5]*GtG[23])
          + m[4]*(m[0]*GtG[24] + m[1]*GtG[25] + m[2]*GtG[26]
                + m[3]*GtG[27] + m[4]*GtG[28] + m[5]*GtG[29])
          + m[5]*(m[0]*GtG[30] + m[1]*GtG[31] + m[2]*GtG[32]
                + m[3]*GtG[33] + m[4]*GtG[34] + m[5]*GtG[35]);
    return mtGtG;
} 
//============================================================================//
static int performL1Search64f(const int nmt, const int ldm,
                              const int npts,
                              const int blockSize, const int mblock,
                              const int Mrows, const int Kcols,
                              const int klag,
                              const bool lminLags, const bool lwantLags,
                              const bool lrescale,
                              const double dnorm,
                              const double *__restrict__ Dmat,
                              const double *__restrict__ CG,
                              const double *__restrict__ Cd,
                              const double *__restrict__ mts,
                              int *__restrict__ lags,
                              double *__restrict__ phi,
                              double *__restrict__ var)
{
    const char *fcnm = "performL1Search64f\0";
    double *M, *R, *rnorm, *U, *unorm, *varLoc, *xden, robj, xnum, viMin;
    int i, ic, idx, jdx, ierr, imt, jmt, kmt, Ncols, nmtBlocks;
    const double one = 1.0;
    const double zero = 0.0;
    ierr = 0;
    nmtBlocks = (int) ((double) (nmt)/(double) (blockSize) + 1.0);
    if (nmtBlocks*blockSize < nmt) 
    {    
        printf("%s: Internal error - all mts wont be visited\n", fcnm);
        return -1;
    }
    // Get the size of the numerator - recall Dmat is simply a copy of itself
    xnum = 1.0;
    if (lrescale)
    {
        ic = 0;
        i = 0;
        xnum = fabs(Dmat[mblock*i+ic]);
        for (i=0; i<npts; i++)
        {
            xnum = fmax(xnum, fabs(Dmat[mblock*i+ic]));
        }
    }
    #pragma omp parallel \
     firstprivate(dnorm, xnum) \
     private (i, ic, idx, imt, jdx, jmt, M, Ncols, R, rnorm, robj, viMin, \
              U, unorm, xden, varLoc) \
     shared (Cd, CG, Dmat, lags, mts, nmtBlocks, lrescale, phi, var) \
     default (none) reduction(+:ierr)
    {
    varLoc = array_copy64f(npts, var, &ierr);
    rnorm = memory_calloc64f(mblock);
    unorm = memory_calloc64f(mblock);
    xden = memory_calloc64f(mblock);
    M = memory_calloc64f(mblock*8);
    R = memory_calloc64f(mblock*npts);
    U = memory_calloc64f(mblock*npts);
#ifdef __INTEL_COMPILER
    __assume_aligned(varLoc, 64);
    __assume_aligned(rnorm, 64);
    __assume_aligned(unorm, 64);
    __assume_aligned(xden, 64);
    __assume_aligned(M, 64);
    __assume_aligned(U, 64);
    __assume_aligned(R, 64);
    __assume_aligned(CG, 64);
    __assume_aligned(Cd, 64);
#endif
    #pragma omp for
    for (kmt=0; kmt<nmtBlocks; kmt++)
    {
        jmt = kmt*blockSize;
        Ncols = MIN(blockSize, nmt - jmt); // Number of columns of M
        // Set the moment tensor matrix 
        for (i=0; i<6; i++)
        {
            #pragma omp simd aligned(M: 64)
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic; 
                idx = ldm*imt + i;
                jdx = mblock*i;
                M[jdx+ic] = mts[idx];
            }
        }
        // Compute U = GM
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    Mrows, Ncols, Kcols, one, CG, LDCG,
                    M, mblock, zero, U, mblock);
        // Rescale the synthetic data
        if (lrescale)
        {
            i = 0;
            for (ic=0; ic<Ncols; ic++)
            {
                xden[ic] = fabs(U[mblock*i+ic]);
            }
            i = 1;
            for (i=0; i<npts; i++)
            {
                for (ic=0; ic<Ncols; ic++)
                {
                    xden[ic] = fmax(xden[ic], fabs(U[mblock*i+ic]));
                }
            }
            // Rescale the synthetics to the observed 
            for (ic=0; ic<Ncols; ic++)
            {
                xden[ic] = xnum/xden[ic];
                cblas_dscal(npts, xden[ic], &U[ic], mblock);
                //int imax = cblas_idamax(npts, &U[ic], mblock);
                //int jmax = cblas_idamax(npts, &Dmat[ic], mblock);
                //printf("%e %e\n", U[imax*mblock+ic], Dmat[jmax*mblock+ic]);
            }
            //getchar();
        }
        memset(rnorm, 0, (size_t) Ncols*sizeof(double));
        memset(unorm, 0, (size_t) Ncols*sizeof(double));
        #pragma omp simd aligned(R, Dmat, U: 64)
        for (i=0; i<mblock*npts; i++)
        {
            R[i] = fabs(Dmat[i] - U[i]);
        }
        // Compute |d - u|_1 
        for (i=0; i<npts; i++)
        {
            viMin = varLoc[i];
            #pragma omp simd aligned(rnorm, R: 64) reduction(min:viMin)
            for (ic=0; ic<Ncols; ic++)
            {
                rnorm[ic] = rnorm[ic] + R[mblock*i+ic];
                viMin = fmin(viMin, R[mblock*i+ic]);
            }
            varLoc[i] = viMin;
        }
        // Compute |u|_1
        for (i=0; i<npts; i++)
        {
            #pragma omp simd aligned(unorm, U: 64)
            for (ic=0; ic<Ncols; ic++)
            {
                unorm[ic] = unorm[ic] + fabs(U[mblock*i+ic]);
            }
        }
        // Apply the triangle inequality
        if (!lminLags)
        {
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic;
                //printf("%d %e %e\n", imt, Dmat[npts/2+ic], unorm[ic]);
                phi[imt] = 1.0 - rnorm[ic]/(dnorm + unorm[ic]);
            }
        }
        else
        {
            if (!lwantLags)
            {
                for (ic=0; ic<Ncols; ic++)
                {
                    imt = jmt + ic;
                    robj = 1.0 - rnorm[ic]/(dnorm + unorm[ic]);
                    phi[imt] = fmax(phi[imt], robj);
                }
            }
            else
            {
                for (ic=0; ic<Ncols; ic++)
                {
                    imt = jmt + ic;
                    robj = 1.0 - rnorm[ic]/(dnorm + unorm[ic]);
                    if (robj > phi[imt])
                    {
                        phi[imt] = robj;
                        lags[imt] = klag;
                    }
                }
            }
        }
    } // Loop on moment tensor block
    #pragma omp critical
    {
    for (i=0; i<npts; i++)
    {
        var[i] = fmin(var[i], varLoc[i]);
    }
    } // End critical section
    memory_free64f(&varLoc);
    memory_free64f(&rnorm);
    memory_free64f(&unorm);
    memory_free64f(&xden);
    memory_free64f(&M);
    memory_free64f(&R);
    memory_free64f(&U);
    } // end parallel section
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes Cd vector
 *
 * @param[in] npts   number of points
 * @param[in] ldc    leading dimension of CeInv.  if ldiag is true then
 *                   this should be 1 otherwise, this should be >= npts.
 * @param[in] CeInv  inverse of data error (covariance) matrix.  if
 *                   ldiag is true then this is diagonal of length [nobs].
 *                   if false then it is a symmetric [nobs x ldc] matrix
 *                   in row major format where ldc >= nobs.
 * @param[in] d      data vector [npts]
 *
 * @param[out] Cd    multiplication of C*d [npts]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeCd32f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const float *__restrict__ CeInv,
                        const float *__restrict__ d,
                        float *__restrict__ Cd)
{
    int i, M, N;
    // Compute C d
    if (CeInv != NULL)
    {
        if (ldiag)
        {
            #pragma omp simd aligned(Cd, CeInv, d: 64)
            for (i=0; i<npts; i++)
            {
                Cd[i] = CeInv[i]*d[i];
            }
        }
        else
        {
            M = npts;
            N = npts;
            cblas_sgemv(CblasRowMajor, CblasNoTrans, M, N, 1.0f, CeInv, ldc,
                        d, 1, 0.0f, Cd, 1);
        }
    }    
    else 
    {    
        // Set C d but now C is assume to be the identity matrix 
        cblas_scopy(npts, d, 1, Cd, 1);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes Cd vector
 *
 * @param[in] npts   number of points
 * @param[in] ldc    leading dimension of CeInv.  if ldiag is true then
 *                   this should be 1 otherwise, this should be >= npts.
 * @param[in] CeInv  inverse of data error (covariance) matrix.  if
 *                   ldiag is true then this is diagonal of length [nobs].
 *                   if false then it is a symmetric [nobs x ldc] matrix
 *                   in row major format where ldc >= nobs.
 * @param[in] d      data vector [npts]
 *
 * @param[out] Cd    multiplication of C*d [npts]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeCd64f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const double *__restrict__ CeInv,
                        const double *__restrict__ d,
                        double *__restrict__ Cd)
{
    int i, M, N;
    // Compute C d
    if (CeInv != NULL)
    {
        if (ldiag)
        {
            #pragma omp simd aligned(Cd, CeInv, d: 64)
            for (i=0; i<npts; i++)
            {
                Cd[i] = CeInv[i]*d[i];
            }
        }
        else
        {
            M = npts;
            N = npts;
            cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, 1.0, CeInv, ldc,
                         d, 1, 0.0, Cd, 1);
        }
    }
    else
    {
        // Set C d but now C is assume to be the identity matrix 
        cblas_dcopy(npts, d, 1, Cd, 1);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes C*G matrix
 *
 * @param[in] npts     number of points in Green's functions
 * @param[in] ldc      leading dimension of CeInv.  if ldiag is true then
 *                     this should be 1 otherwise, this should be >= npts.
 * @param[in] ldiag    if true then CeInv is diagonal.  if false then
 *                     CeInv is a general [npts x npts] matrix.
 * @param[in] CeInv    inverse of data error (covariance) matrix.  if
 *                     ldiag is true then this is diagonal of length [nobs].
 *                     if false then it is a symmetric [nobs x ldc] matrix
 *                     in row major format where ldc >= nobs.
 * @param[in] G        Green's functions matrix in row major order [npts x LDG]
 *
 * @param[out] CG      row major C*G matrix [npts x LDCG]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeCG32f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const float *__restrict__ CeInv,
                        const float *__restrict__ G,
                        float *__restrict__ CG)
{
    int i, j, jndx, K, kndx, M;
    const int N = 6; // Number of columns of G 
    if (CeInv != NULL)
    {
        if (ldiag)
        {
            // Set C G
            for (i=0; i<npts; i++)
            {
                for (j=0; j<6; j++)
                {
                    jndx = LDG*i + j;
                    kndx = LDCG*i + j;
                    CG[kndx] = CeInv[i]*G[jndx];
                }
            }
        }
        else
        {
            M = npts; // rows of CeInv
            K = npts; // columsn of CeInv
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        M, N, K, 1.0f, CeInv, ldc,
                        G, LDG, 0.0f, CG, LDCG);
        }
    }
    else 
    {    
        // Set C G but now C is assume to be the identity matrix 
        for (i=0; i<npts; i++) 
        {
            for (j=0; j<6; j++) 
            {
                jndx = LDG*i + j; 
                kndx = LDCG*i + j; 
                CG[kndx] = G[jndx];
            }
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes C*G matrix
 *
 * @param[in] npts     number of points in Green's functions
 * @param[in] ldc      leading dimension of CeInv.  if ldiag is true then
 *                     this should be 1 otherwise, this should be >= npts.
 * @param[in] ldiag    if true then CeInv is diagonal.  if false then
 *                     CeInv is a general [npts x npts] matrix.
 * @param[in] CeInv    inverse of data error (covariance) matrix.  if
 *                     ldiag is true then this is diagonal of length [nobs].
 *                     if false then it is a symmetric [nobs x ldc] matrix
 *                     in row major format where ldc >= nobs.
 * @param[in] G        Green's functions matrix in row major order [npts x LDG]
 *
 * @param[out] CG      row major C*G matrix [npts x LDCG]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeCG64f(const int npts,
                        const int ldc,
                        const bool ldiag,
                        const double *__restrict__ CeInv,
                        const double *__restrict__ G,
                        double *__restrict__ CG)
{
    int i, j, jndx, K, kndx, M;
    const int N = 6; // Number of columns of G
    if (CeInv != NULL)
    {
        if (ldiag)
        {
            // Set C G
            for (i=0; i<npts; i++)
            {
                for (j=0; j<6; j++)
                {
                    jndx = LDG*i + j;
                    kndx = LDCG*i + j;
                    CG[kndx] = CeInv[i]*G[jndx];
                }
            }
        }
        else
        {
            M = npts; // rows of CeInv
            K = npts; // columns of CeInv
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        M, N, K, 1.0, CeInv, ldc,
                        G, LDG, 0.0, CG, LDCG);
        }
    }
    else
    {
        // Set C G but now C is assume to be the identity matrix 
        for (i=0; i<npts; i++) 
        {
            for (j=0; j<6; j++) 
            {    
                jndx = LDG*i + j; 
                kndx = LDCG*i + j; 
                CG[kndx] = G[jndx];
            }    
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes G^T C G matrix
 *
 * @param[in] npts    number of points in Green's functions
 * @param[in] G       row major Green's functions matrix [npts x LDG]
 * @param[in] CG      regularized Green's function matrix [npts x LDCG]
 *
 * @param[out] GtCG   symmetric row major matrix containing G^T C G [6 x 6]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeGtCG32f(const int npts,
                          const float *__restrict__ G,
                          const float *__restrict__ CG,
                          float *__restrict__ GtCG)
{
    const int M = 6; // rows of trans(G)
    const int N = 6; // columns of CG
    int K = npts;    // columns of trans(G)
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                M, N, K, 1.0, G, LDG,
                CG, LDCG, 0.0, GtCG, 6);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes G^T C G matrix
 *
 * @param[in] npts    number of points in Green's functions
 * @param[in] G       row major Green's functions matrix [npts x LDG]
 * @param[in] CG      regularized Green's function matrix [npts x LDCG]
 *
 * @param[out] GtCG   symmetric row major matrix containing G^T C G [6 x 6]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
static int computeGtCG64f(const int npts,
                          const double *__restrict__ G,
                          const double *__restrict__ CG,
                          double *__restrict__ GtCG)
{
    const int M = 6; // rows of trans(G)
    const int N = 6; // columns of CG
    int K = npts;    // columns of trans(G)
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                M, N, K, 1.0, G, LDG,
                CG, LDCG, 0.0, GtCG, 6);
    return 0;
}
//============================================================================//
/*!
 * @brief A utility function for computing: 
 *        \f$
 *          \log \left ( \sum_i^n exp(-x_i) \right )
 *        \f$.
 *
 * For an explanation of why I do this whacky arithmetic see:
 * http://jblevins.org/log/log-sum-exp
 *
 * @param[in] n     number of points in array x
 * @param[in] x     objective function values at n points.  notice that these
 *                  are all non-negative.
 *
 * @result log of the sum of the exponential terms
 * 
 * @author Ben Baker
 *
 */
double mtsearch_logSumExp64f(const int n, const double *__restrict__ x)
{
    double arg, sum, hatx, lse;
    int i, iopt;
    // term dominant player in the integral is the `best' or smallest term 
    iopt = cblas_idamin(n, x, 1); 
    hatx = x[iopt];
    // have fun with exponentials and find the constant 'c' s.t.:
    //   log((exp(-x_1 + \hat{x}) + exp(-x_2 + \hat{x}) + \cdots )exp(-\hat{x}))
    // = log(exp(-x_1 + \hat{x} + exp(exp(-x_2 + \hat{x})) + \cdots) - \hat{x}
    sum = 0.0;
    #pragma omp simd reduction(+:sum)
    for (i=0; i<n; i++)
    {
        arg = hatx - x[i];
        sum = sum + exp(arg);
    }
    lse = log(sum) - hatx;
    return lse;
}

static int computePadding64f(const int n)
{
    size_t mod, pad;
    int ipad;
    // Set space and make G matrix
    pad = 0; 
    mod = ((size_t) n*sizeof(double))%64;
    if (mod != 0)
    {
        pad = (64 - mod)/sizeof(double);
    }
    ipad = (int) pad;
    return ipad;
}

static int computePadding32f(const int n)
{
    size_t mod, pad;
    int ipad;
    // Set space and make G matrix
    pad = 0; 
    mod = ((size_t) n*sizeof(float))%64;
    if (mod != 0)
    {
        pad = (64 - mod)/sizeof(float);
    }
    ipad = (int) pad;
    return ipad;
}
