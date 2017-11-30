#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "parmt_mtsearch.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/signal/convolve.h"

int parmt_computeStackedCrossCorrelation_MPI(
    const MPI_Comm comm,
    const int npgrns,
    const double *__restrict__ Gxx,
    const double *__restrict__ Gyy,
    const double *__restrict__ Gzz,
    const double *__restrict__ Gxy,
    const double *__restrict__ Gxz,
    const double *__restrict__ Gyz,
    const int ldm, const int nmt, const double *__restrict__ mt, 
    const int npts, const double *__restrict__ data,
    const int lxc, double *__restrict__ xc)
{
    double *xcwork;
    int ierr, lx;
    ierr = 0;
    lx = npts + npgrns - 1;
    xcwork = memory_calloc64f(lxc);
    ierr = parmt_computeStackedCrossCorrelation(npgrns,
                                                Gxx, Gyy, Gzz,
                                                Gxy, Gxz, Gyz,
                                                ldm, nmt, mt,
                                                npts, data,
                                                lx, xcwork);
    MPI_Allreduce(xcwork, xc, lxc, MPI_DOUBLE, MPI_SUM, comm); 
    memory_free64f(&xcwork);
    return ierr; 
}
//============================================================================//
/*!
 * @brief This is an alignment tool intended to optimize the lag by computing
 *
 *        \f[
 *           x_c = \frac{1}{n_{mt}} \sum_{i=1}^{n_{mt}} \textbf{u}_i \textbf{d}
 *        \f]
 *
 *        where the \f$ \star \f$ implies a normalized cross-correlation.
 *
 * @param[in] npgrns    number of points in Green's functions
 * @param[in] Gxx       Green's functions scaled by \f$ m_{xx} \f$ [npgrns] 
 * @param[in] Gyy       Green's functions scaled by \f$ m_{yy} \f$ [npgrns]
 * @param[in] Gzz       Green's functions scaled by \f$ m_{zz} \f$ [npgrns]
 * @param[in] Gxy       Green's functions scaled by \f$ m_{xy} \f$ [npgrns]
 * @param[in] Gxz       Green's functions scaled by \f$ m_{xz} \f$ [npgrns]
 * @param[in] Gyz       Green's functions scaled by \f$ m_{yz} \f$ [npgrns]
 * @param[in] ldm       leading dimension of moment tensors (>= 6)
 * @param[in] nmt       number of moment tensors
 * @param[in] mt        moment tensors.  the i'th moment tensor begins at
 *                      i*ldm and the moment tensor terms are packed 
 *                      \f$ \{ m_{xx}, m_{yy}, m_{zz},
 *                             m_{xy}, m_{xz}, m_{yz} \f$.
 * @param[in] npts      number of points in signal
 * @param[in] data      Observed data.  This is an array of dimension [npts].
 * @param[in] lxc       max size of xc (should be >= npts + npgrns - 1)
 *
 * @param[out] xc       stack of the observed and synthetic cross-correlations
 *                      [lxc].  the C indexed lag is computed:
 *                        npts - 1 + argmax(xc)
 *                      where the argmax takes values [0,lx-1] 
 *
 * @result 0 indicates success
 *        
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_computeStackedCrossCorrelation(
    const int npgrns,
    const double *__restrict__ Gxx,
    const double *__restrict__ Gyy,
    const double *__restrict__ Gzz,
    const double *__restrict__ Gxy,
    const double *__restrict__ Gxz,
    const double *__restrict__ Gyz,
    const int ldm, const int nmt, const double *__restrict__ mt,
    const int npts, const double *__restrict__ data,
    const int lxc, double *__restrict__ xc)
{
    const char *fcnm = "parmt_computedStackedCrossCorrelation\0";
    double *Gmat, *Gxc, *xcorr, *xcw;
    double m6[8] __attribute__((aligned(64)));
    double dscal __attribute__((aligned(64))) = 0.0;
    double gscal __attribute__((aligned(64))) = 0.0;
    double xscal __attribute__((aligned(64))) = 0.0;
    double est __attribute__((aligned(64))) = 0.0; 
    int i, ierr, j, lc;
    lc = npts + npgrns - 1;
    if (npts < 1 || npgrns < 1 ||
        Gxx == NULL || Gyy == NULL || Gzz == NULL ||
        Gxy == NULL || Gxz == NULL || Gyz == NULL ||
        data == NULL)
    {
        if (npts < 1){fprintf(stderr, "%s: no data points\n", __func__);}
        if (npgrns < 1){fprintf(stderr, "%s: no grns fns points\n", __func__);}
        if (data == NULL){fprintf(stderr, "%s: data is NULL\n", __func__);}
        if (Gxx == NULL){fprintf(stderr, "%s: Gxx is NULL\n", __func__);}
        if (Gyy == NULL){fprintf(stderr, "%s: Gyy is NULL\n", __func__);}
        if (Gzz == NULL){fprintf(stderr, "%s: Gzz is NULL\n", __func__);}
        if (Gxy == NULL){fprintf(stderr, "%s: Gxy is NULL\n", __func__);}
        if (Gxz == NULL){fprintf(stderr, "%s: Gxz is NULL\n", __func__);}
        if (Gyz == NULL){fprintf(stderr, "%s: Gyz is NULL\n", __func__);}
        return -1;
    } 
    if (lc > lxc)
    {
        fprintf(stderr, "%s: Error xc is too small\n", __func__);
        return -1;
    }
    if (ldm < 6 || nmt < 1 || mt == NULL)
    {
        if (ldm < 6){fprintf(stderr, "%s: ldm must be >= 6\n", __func__);}
        if (nmt < 1){fprintf(stderr, "%s: no moment tensors\n", __func__);}
        if (mt == NULL){fprintf(stderr, "%s: mt is NULL\n", __func__);}
        return -1;
    }
    // Set space and correlate data with columns of G
    xcorr = memory_calloc64f(lc);
    Gmat = memory_calloc64f(8*npgrns);
    Gxc = memory_calloc64f(8*lc);
    dscal = cblas_dnrm2(npts, data, 1);
    dscal = dscal*dscal;
#ifdef __INTEL_COMPILER
__assume_aligned(Gxc, 64);
__assume_aligned(Gmat, 64);
__assume_aligned(xcorr, 64);
__assume_aligned(xc, 64);
#endif
    for (i=0; i<6; i++)
    {
        if (i == 0)
        {
            cblas_dcopy(npgrns, Gxx, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gxx,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        else if (i == 1)
        {
            cblas_dcopy(npgrns, Gyy, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gyy,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        else if (i == 2)
        {
            cblas_dcopy(npgrns, Gzz, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gzz,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        else if (i == 3)
        {
            cblas_dcopy(npgrns, Gxy, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gxy,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        else if (i == 4)
        {
            cblas_dcopy(npgrns, Gxz, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gxz,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        else if (i == 5)
        {
            cblas_dcopy(npgrns, Gyz, 1, &Gmat[i], 8);
            ierr = signal_convolve_correlate64f_work(npts, data,
                                                     npgrns,
                                                     Gyz,
                                                     CONVCOR_FULL,
                                                     lc, xcorr);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error correlating Greens fn: %d\n",
                    __func__, i+1);
            return -1;
        }
        cblas_dcopy(lc, xcorr, 1, &Gxc[i], 8);
    }
    array_zeros64f_work(lxc, xc);
    // Compute \f$ x_c = \sum_{i=1}^{n_{mt}} G \textbf{m} \star \textbf{d} \f$
#pragma omp parallel shared(Gxc, Gmat, mt, lc, xc) \
        firstprivate (dscal, xcorr) \
        private (est, gscal, xscal, i, j, m6, xcw) default(none)
{
    xcw = memory_calloc64f(lxc);
#ifdef __INTEL_COMILER
__assume_aligned(xcw, 64);
#endif
    #pragma omp for
    for (i=0; i<nmt; i++)
    {
        m6[0] = mt[ldm*i  ];
        m6[1] = mt[ldm*i+1];
        m6[2] = mt[ldm*i+2];
        m6[3] = mt[ldm*i+3];
        m6[4] = mt[ldm*i+4];
        m6[5] = mt[ldm*i+5];
        // compute the energy in the synthetic
        gscal = 0.0;
        for (j=0; j<npgrns; j++)
        {
            est = Gmat[8*j+0]*m6[0] + Gmat[8*j+1]*m6[1] + Gmat[8*j+2]*m6[2]
                + Gmat[8*j+3]*m6[3] + Gmat[8*j+4]*m6[4] + Gmat[8*j+5]*m6[5];
            gscal = gscal + est*est;
        }
        xscal = 1.0/sqrt(dscal*gscal);
        // compute the normalized cross-correlation
        for (j=0; j<lc; j++)
        {
            est = Gxc[8*j+0]*m6[0] + Gxc[8*j+1]*m6[1] + Gxc[8*j+2]*m6[2]
                + Gxc[8*j+3]*m6[3] + Gxc[8*j+4]*m6[4] + Gxc[8*j+5]*m6[5];
            xcw[j] = xcw[j] + xscal*est;
        }
        //cblas_dgemv(CblasRowMajor, CblasNoTrans,
        //            lc, 6, 1.0, G, 8, &mt[8*i], 1, 1.0, xcw, 1);
    }
    #pragma omp critical
    #pragma omp simd
    for (i=0; i<lxc; i++)
    {
        xc[i] = xc[i] + xcw[i];
    }
    memory_free64f(&xcw);
}
    // normalize by the number of moment tensors
    cblas_dscal(lc, 1.0/(double) nmt, xc, 1); 
    // free memory
    memory_free64f(&Gxc);
    memory_free64f(&Gmat);
    memory_free64f(&xcorr);
    return 0;
}
