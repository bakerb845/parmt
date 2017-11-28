#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

/*!
 * @brief Computes the [mrows x 6] forward modeling matrix for synthesizing
 *        a waveform.
 *
 * @param[in] sacGxx   data structure describing the Green's functions to be
 *                     scaled by the mxx Green's functions
 * @param[in] sacGyy   data structure describing the Green's functions to be
 *                     scaled by the myy Green's functions
 * @param[in] sacGzz   data structure describing the Green's functions to be
 *                     scaled by the mzz Green's functions
 * @param[in] sacGxy   data structure describing the Green's functions to be
 *                     scaled by the mxy Green's functions
 * @param[in] sacGxz   data structure describing the Green's functions to be
 *                     scaled by the mxz Green's functions
 * @param[in] sacGyz   data structure describing the Green's functions to be
 *                     scaled by the myz Green's functions
 * 
 * @param[out] mrows   number of rows in forward modeling matrix 
 * @param[out] ierr    0 indicates success
 *
 * @result the [npts x 6] column major forward modeling matirx packed 
 *         \f$ \{G_{xx}, G_{yy}, G_{zz}, G_{xy}, G_{xz}, G_{yz} \f$. 
 *         this must be freed with ISCL's memory_free64f.

 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
double *parmt_utils_sacGrns2GrnsMat(const struct sacData_struct sacGxx,
                                    const struct sacData_struct sacGyy,
                                    const struct sacData_struct sacGzz,
                                    const struct sacData_struct sacGxy,
                                    const struct sacData_struct sacGxz,
                                    const struct sacData_struct sacGyz, 
                                    int *mrows, int *ierr)
{
    const char *fcnm = "parmt_utils_sacGrns2GrnsMat\0";
    double *G;
    int npts, nptsGxx, nptsGyy, nptsGzz, nptsGxy, nptsGxz, nptsGyz;
    G = NULL;
    *mrows = 0;
    // Verify the inputs
    *ierr  = sacio_getIntegerHeader(SAC_INT_NPTS, sacGxx.header, &nptsGxx);
    *ierr += sacio_getIntegerHeader(SAC_INT_NPTS, sacGyy.header, &nptsGyy);
    *ierr += sacio_getIntegerHeader(SAC_INT_NPTS, sacGzz.header, &nptsGzz);
    *ierr += sacio_getIntegerHeader(SAC_INT_NPTS, sacGxy.header, &nptsGxy);
    *ierr += sacio_getIntegerHeader(SAC_INT_NPTS, sacGxz.header, &nptsGxz);
    *ierr += sacio_getIntegerHeader(SAC_INT_NPTS, sacGyz.header, &nptsGyz);
    if (*ierr != 0 || nptsGxx != nptsGyy ||
        nptsGxx != nptsGzz || nptsGxx != nptsGxy ||
        nptsGxx != nptsGxz || nptsGxx != nptsGyz)
    {
        if (*ierr != 0){printf("%s: Error reading headers\n", fcnm);}
        printf("%s: Inconsistent Green's functions trace lenths\n", fcnm);
        return G;
    }
    if (nptsGxx < 1)
    {
        printf("%s: No points in Green's functions\n", fcnm);
        return G;
    }
    if (nptsGxx != sacGxx.npts || nptsGxx != sacGyy.npts ||
        nptsGxx != sacGzz.npts || nptsGxx != sacGxy.npts ||
        nptsGxx != sacGxz.npts || nptsGxx != sacGyz.npts)
    {
        printf("%s: header npts differs from data npts\n", fcnm);
        return G;
    }
    if (sacGxx.data == NULL || sacGyy.data == NULL || sacGzz.data == NULL ||
        sacGxy.data == NULL || sacGxz.data == NULL || sacGyz.data == NULL)
    {
        printf("%s: null grns fns detected\n", fcnm);
        return G;
    }
    // Finally set the column major modeling matrix 
    npts = nptsGxx;
    *mrows = npts;
    G = memory_calloc64f(npts*6);
    cblas_dcopy(npts, sacGxx.data, 1, &G[0*npts], 1);
    cblas_dcopy(npts, sacGyy.data, 1, &G[1*npts], 1);
    cblas_dcopy(npts, sacGzz.data, 1, &G[2*npts], 1);
    cblas_dcopy(npts, sacGxy.data, 1, &G[3*npts], 1);
    cblas_dcopy(npts, sacGxz.data, 1, &G[4*npts], 1);
    cblas_dcopy(npts, sacGyz.data, 1, &G[5*npts], 1);
    return G;
}
//============================================================================//
/*!
 * @brief Computes a synthetic seismogram from Green's functions
 *
 */
int parmt_utils_sacGrnsToEst(const struct sacData_struct data,
                             const struct sacData_struct sacGxx,
                             const struct sacData_struct sacGyy,
                             const struct sacData_struct sacGzz,
                             const struct sacData_struct sacGxy,
                             const struct sacData_struct sacGxz,
                             const struct sacData_struct sacGyz,
                             const double *__restrict__ mt,
                             const bool lrescale,
                             struct sacData_struct *est)
{
    const char *fcnm = "parmt_utils_sacGrnsToEst\0";
    double *G;
    int ierr, npts, nptsGxx, nptsGyy, nptsGzz, nptsGxy, nptsGxz, nptsGyz;
    const int ldg = 6;
    double xnum, xden, xscal;
    // Release memory on the target
    sacio_free(est);
    // Build the modeling matrix
    G = parmt_utils_sacGrns2GrnsMat(sacGxx, sacGyy, sacGzz,
                                    sacGxy, sacGxz, sacGyz,
                                    &npts, &ierr);
    if (ierr != 0)
    {
        printf("%s: Failed to make forward modeling matrix\n", fcnm);
        memory_free64f(&G);
    }
    // Copy the sacGxx to the estimate
    sacio_copy(data, est);
    est->header.isynth = 501; // set a number
    // Set the modeling matrix 
    G = memory_calloc64f(npts*ldg);
    cblas_dcopy(npts, sacGxx.data, 1, &G[0*npts], 1);
    cblas_dcopy(npts, sacGyy.data, 1, &G[1*npts], 1);
    cblas_dcopy(npts, sacGzz.data, 1, &G[2*npts], 1);
    cblas_dcopy(npts, sacGxy.data, 1, &G[3*npts], 1);
    cblas_dcopy(npts, sacGxz.data, 1, &G[4*npts], 1);
    cblas_dcopy(npts, sacGyz.data, 1, &G[5*npts], 1);
    // Compute Gm = u to compute the estimate
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                npts, 6, 1.0, G, npts, mt, 1, 0.0, est->data, 1);
    memory_free64f(&G);
    // Potentially rescale the synthetic
    if (lrescale)
    {
        xnum = array_maxAbs64f(npts, data.data, &ierr);
        xden = array_maxAbs64f(npts, est->data, &ierr);
        xscal = xnum/xden;
        cblas_dscal(npts, xscal, est->data, 1); 
    }
    return 0;
}

