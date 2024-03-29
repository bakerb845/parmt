#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_invert.h"
#ifdef PARMT_USE_INTEL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/linalg/linalg.h"
#include "iscl/memory/memory.h"

/*!
 * @brief Inverts for a moment tensor at each location in the grid-search.
 *
 * @param[in] data     parMT data structure with observations and 
 *                     corresponding green's functions at all locations
 *                     for all waveforms.
 * @param[in] ldm      Leading dimension of mts.  This must be at least 6.
 * @param[in] lcDev    If true then apply the deviatoric constraint.
 * @param[in] lcDev    Otherwise invert for all moment tensors terms.
 *
 * @param[out] mts     Moment tensors tabulated at each location.  This has
 *                     dimension [data.nlocs x ldm].
 * @param[out] phi     L2 objective function for each moment tensor.  This has
 *                     dimension [data.nlocs].
 *
 * @result 0 indicates success.
 *
 */
int parmt_invertLocationSearch(struct parmtData_struct data,
                               const int ldm,
                               const bool lcDev,
                               double *__restrict__ mts,
                               double *__restrict__ phi)
{
    double *G, *gxxAll, *gyyAll, *gzzAll, *gxyAll, *gxzAll, *gyzAll,
           *obsAll, *wts;
    //double m6[6], phiLoc;
    int *dataPtr, i1, ierr, ierr1, ierr2, iloc, iobs, k, mrows;
    ierr = 0;
    if (ldm < 6 || mts == NULL || phi == NULL)
    {
        if (ldm < 6){fprintf(stderr, "%s: ldm must be at least 6\n", __func__);}
        if (mts == NULL){fprintf(stderr, "%s: mts is NULL\n", __func__);}
        if (phi == NULL){fprintf(stderr, "%s: phi is NULL\n", __func__);}
        return -1;
    }
    G = NULL;
    wts = NULL;
    gxxAll = NULL; 
    gyyAll = NULL;
    gzzAll = NULL;
    gxyAll = NULL;
    gxzAll = NULL;
    gyzAll = NULL;
    // 6 observations required for full MT, 5 for deviatoric constraint
    if (data.nobs < 6 - (int) lcDev)
    {
        fprintf(stderr, "%s: Insufficient number of observations \n",
                __func__);
        return -1;
    }
    // Figure out the workspace
    dataPtr = memory_calloc32i(data.nobs+1); 
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        dataPtr[iobs+1] = dataPtr[iobs] + data.data[iobs].npts;
    }
    if (dataPtr[data.nobs] < 6 - (int) lcDev)
    {
        fprintf(stderr, "%s: Insufficient number of points to invert with\n",
                __func__);
        goto ERROR;
    }
    obsAll = memory_calloc64f(dataPtr[data.nobs]);
    gxxAll = memory_calloc64f(dataPtr[data.nobs]);
    gyyAll = memory_calloc64f(dataPtr[data.nobs]);
    gzzAll = memory_calloc64f(dataPtr[data.nobs]);
    gxyAll = memory_calloc64f(dataPtr[data.nobs]);
    gxzAll = memory_calloc64f(dataPtr[data.nobs]); 
    gyzAll = memory_calloc64f(dataPtr[data.nobs]);
    // Loop on locations
    for (iloc=0; iloc<data.nlocs; iloc++)
    {
        // Loop on observations
        ierr1 = 0;
        ierr2 = 0;
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            k = iobs*data.nlocs + iloc;;
            G = parmt_utils_sacGrns2GrnsMat(data.sacGxx[k], data.sacGyy[k],
                                            data.sacGzz[k], data.sacGxy[k],
                                            data.sacGxz[k], data.sacGyz[k],
                                            &mrows, &ierr1);
            memory_free64f(&G);
            if (ierr1 != 0)
            {
                fprintf(stderr, "%s: Failed to make greens functions matrix\n",
                        __func__);
                ierr2 = ierr2 + 1;
            }
            if (mrows != dataPtr[iobs+1] - dataPtr[iobs])
            {
                fprintf(stderr, "%s: Internal error on mrows and npts\n",
                        __func__);
                ierr2 = ierr2 + 1;
                goto NEXT_OBS;
            }
            // Copy the observations and Green's functions
            i1 = dataPtr[iobs];
            cblas_dcopy(mrows, &G[0*mrows], 1, &gxxAll[i1], 1);
            cblas_dcopy(mrows, &G[1*mrows], 1, &gyyAll[i1], 1);
            cblas_dcopy(mrows, &G[2*mrows], 1, &gzzAll[i1], 1);
            cblas_dcopy(mrows, &G[3*mrows], 1, &gxyAll[i1], 1);
            cblas_dcopy(mrows, &G[4*mrows], 1, &gxzAll[i1], 1);
            cblas_dcopy(mrows, &G[5*mrows], 1, &gyzAll[i1], 1);
            cblas_dcopy(mrows, data.data[iobs].data, 1, &obsAll[i1], 1); 
NEXT_OBS:; 
            memory_free64f(&G); 
             
        } // loop on observations
        // Invert it
        if (ierr2 != 0)
        {
            ierr1 = parmt_invertMT64f(data.nobs,
                                      lcDev, dataPtr,
                                      wts,
                                      gxxAll, gyyAll, gzzAll,
                                      gxyAll, gxzAll, gyzAll,
                                      obsAll, &phi[iloc], &mts[iloc*ldm]);
            //&phiLoc, m6);
            if (ierr1 != 0)
            {
                fprintf(stderr, "%s: Error inverting location %d\n",
                        __func__, iloc);
                ierr = ierr + 1;
            }
        }
        else
        {
            ierr = ierr + 1;
        }
    }
ERROR:;
    memory_free32i(&dataPtr);
    memory_free64f(&gxxAll);
    memory_free64f(&gyyAll);
    memory_free64f(&gzzAll);
    memory_free64f(&gxyAll);
    memory_free64f(&gxzAll);
    memory_free64f(&gyzAll);
    return 0;
}
//============================================================================//
/*!
 * @brief
 *
 * @param[in] nobs     number of waveforms (observations)
 * @param[in] lcDev    if true then force the solution to be purely deviatoric.
 *                     otherwise all moment tensor terms will be inverted for.
 * @param[in] dataPtr  points from iobs'th observation to start index of 
 *                     observation [nobs+1]
 * @param[in] wts      if not NULL then these are the data weights for each
 *                     observations [nobs].
 * @param[in] gxxAll   all Gxx greens functions [dataPtr[nobs]] 
 * @param[in] gyyAll   all Gyy greens functions [dataPtr[nobs]]
 * @param[in] gzzAll   all Gzz greens functions [dataPtr[nobs]]
 * @param[in] gxyAll   all Gxy greens functions [dataPtr[nobs]]
 * @param[in] gxzAll   all Gxz greens functions [dataPtr[nobs]]
 * @param[in] gyzAll   all Gyz greens functions [dataPtr[nobs]]
 * @param[in] obsAll   observations [dataPtr[nobs]]
 *
 * @param[out] phi     L2 difference between observed and estimate seismograms. 
 * @param[out] m       on successful contains the moment tensor elements packed
 *                     \f$\{ m_{xx}, m_{yy}, m_{zz},
 *                           m_{xy}, m_{xz}, m_{yz} \}\f$. 
 *
 * @result 0 indicates success 
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_invertMT64f(const int nobs,
                      const bool lcDev,
                      const int *dataPtr,
                      const double *__restrict__ wts,
                      const double *__restrict__ gxxAll,
                      const double *__restrict__ gyyAll,
                      const double *__restrict__ gzzAll,
                      const double *__restrict__ gxyAll,
                      const double *__restrict__ gxzAll,
                      const double *__restrict__ gyzAll,
                      const double *__restrict__ obsAll,
                      double *phi, double *__restrict__ m)
{
    double *G, *b, *est, m6[6];
    int i1, i2, ierr, iobs, ldg, mrows, ncols, nploc;
    *phi = 0.0;
    ierr = 0;
    // Create the forward modeling matrix
    if (nobs < 5 || dataPtr == NULL ||
        gxxAll == NULL || gyyAll == NULL || gzzAll == NULL ||
        gxyAll == NULL || gxzAll == NULL || gyzAll == NULL ||
        m == NULL)
    {
        if (nobs < 5)
        {
            fprintf(stderr, "%s: At least 6 observations required\n", __func__);
        }
        if (dataPtr == NULL)
        {
            fprintf(stderr, "%s: Error dataPtr is NULL\n", __func__);
        }
        if (gxxAll == NULL)
        {
            fprintf(stderr, "%s: Error gxxAll is NULL\n", __func__);
        }
        if (gyyAll == NULL)
        {
            fprintf(stderr, "%s: Error gyyAll is NULL\n", __func__);
        }
        if (gzzAll == NULL)
        {
            fprintf(stderr, "%s: Error gzzAll is NULL\n", __func__);
        }
        if (gxyAll == NULL)
        {
            fprintf(stderr, "%s: Error gxyAll is NULL\n", __func__);
        }
        if (gxzAll == NULL)
        {
            fprintf(stderr, "%s: Error gxzAll is NULL\n", __func__);
        }
        if (gyzAll == NULL)
        {
            fprintf(stderr, "%s: Error gyzAll is NULL\n", __func__);
        }
        if (m == NULL){fprintf(stderr, "%s: Error m is NULL\n", __func__);}
        return -1;
    }
    memset(m, 0, 6*sizeof(double));
    ncols = 6;
    if (lcDev){ncols = 5;}
    // Get the data size 
    mrows = 0; 
    for (iobs=0; iobs<nobs; iobs++)
    {
        mrows = mrows + dataPtr[iobs-1] - dataPtr[iobs];
    }
    if (mrows < ncols)
    {
        fprintf(stderr, "%s: Error problem is undetermined\n", __func__);
        return -1;
    }
    ldg = mrows;
    G = memory_calloc64f(ldg*ncols);
    b = memory_calloc64f(mrows);
    est = memory_calloc64f(mrows);
    // Set the forward modeling matrix 
    if (wts == NULL)
    {
        if (lcDev)
        {
            cblas_dcopy(mrows, gxxAll, 1, &G[0*ldg], 1);
            cblas_daxpy(mrows, -1.0, gzzAll, 1, &G[0*ldg], 1); // mxx - mzz
            cblas_dcopy(mrows, gyyAll, 1, &G[1*ldg], 1);
            cblas_daxpy(mrows, -1.0, gzzAll, 1, &G[1*ldg], 1); // myy - mzz
            cblas_dcopy(mrows, gxyAll, 1, &G[2*ldg], 1);
            cblas_dcopy(mrows, gxzAll, 1, &G[3*ldg], 1);
            cblas_dcopy(mrows, gyzAll, 1, &G[4*ldg], 1);
        }
        else
        {
            cblas_dcopy(mrows, gxxAll, 1, &G[0*ldg], 1);
            cblas_dcopy(mrows, gyyAll, 1, &G[1*ldg], 1);
            cblas_dcopy(mrows, gyyAll, 1, &G[2*ldg], 1);
            cblas_dcopy(mrows, gxyAll, 1, &G[3*ldg], 1);
            cblas_dcopy(mrows, gxzAll, 1, &G[4*ldg], 1);
            cblas_dcopy(mrows, gyzAll, 1, &G[5*ldg], 1);
        }
        cblas_dcopy(mrows, obsAll, 1, b, 1);
    }
    else
    {
        if (lcDev)
        {
            for (iobs=0; iobs<nobs; iobs++)
            {
                i1 = dataPtr[iobs];
                i2 = dataPtr[iobs];
                nploc = i2 - i1;
                cblas_daxpy(nploc, wts[iobs], &gxxAll[i1], 1, &G[0*ldg+i1], 1);
                cblas_daxpy(nploc,-wts[iobs], &gzzAll[i1], 1, &G[0*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gyyAll[i1], 1, &G[1*ldg+i1], 1);
                cblas_daxpy(nploc,-wts[iobs], &gzzAll[i1], 1, &G[1*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gxyAll[i1], 1, &G[2*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gxzAll[i1], 1, &G[3*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gyzAll[i1], 1, &G[4*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &obsAll[i1], 1, &b[i1], 1);
            }
        }
        else
        {
            for (iobs=0; iobs<nobs; iobs++)
            {
                i1 = dataPtr[iobs];
                i2 = dataPtr[iobs];
                nploc = i2 - i1;
                cblas_daxpy(nploc, wts[iobs], &gxxAll[i1], 1, &G[0*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gyyAll[i1], 1, &G[1*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gzzAll[i1], 1, &G[2*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gxyAll[i1], 1, &G[3*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gxzAll[i1], 1, &G[4*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &gyzAll[i1], 1, &G[5*ldg+i1], 1);
                cblas_daxpy(nploc, wts[iobs], &obsAll[i1], 1, &b[i1], 1);
            }
        }
    }
    // Solve \tilde{G} m = \tilde{b}
    ierr = linalg_lstsq_qr64f_work(LAPACK_COL_MAJOR,
                                   mrows, ncols, 1,
                                   false, G, b,
                                   m6, NULL);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error solving Gm=u\n", __func__);
        ierr = 1;
        goto END;
    }
    if (lcDev) 
    {
        m[0] = m6[0];
        m[1] = m6[1];
        m[2] =-(m6[0] + m6[1]);
        m[3] = m6[2];
        m[4] = m6[3];
        m[5] = m6[4];
    }
    else
    {
        cblas_dcopy(6, m6, 1, m, 1);
    }
    // Apply u = \tilde{G} m
    cblas_dgemv(CblasColMajor, CblasNoTrans, mrows, ncols,
                1.0, G, ncols, m6, 1, 0.0, est, 1);
    // Compute the least-squares residual
    *phi = array_normDiff64f(mrows, est, b, TWO_NORM, 2.0, &ierr);
END:;
    // Release memory
    memory_free64f(&G);
    memory_free64f(&b);
    memory_free64f(&est);
    return ierr;
}
