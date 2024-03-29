#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_mtsearch.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "parmt_utils.h"
#include "iscl/memory/memory.h"

/*!
 * @brief Computes a synthetic seismogram corresponding for the 
 *        given observations, location (depth), and moment tensor.
 *
 * @param[in] iobs      Observations number.
 * @param[in] ilocOpt   Location index.
 * @param[in] imtOpt    Moment tensor number. 
 * @param[in] ldm       Leading dimension of mts.  This must be at leats 6.
 * @param[in] mts       Matrix of moement tensors with dimension [nmt x ldm].
 *
 * @param[in,out] data  On input contains the data.
 * @param[in,out] data  On exit contains the corresponding synthetic.
 * 
 * @result 0 indicates success.
 *
 */
int parmt_computeSynthetic(const int iobs,
                           const int ilocOpt,
                           const int imtOpt,
                           const int ldm,
                           const double *__restrict__ mts,
                           struct parmtData_struct *data)
{
    double *G;
    int k, npts;
    k = iobs*data->nlocs + ilocOpt;
    fprintf(stdout, "%s: Processing observation=%d, location=%d, mt=%d, k=%d\n",
            __func__, iobs, ilocOpt, imtOpt, k);
    npts = data->sacGxx[k].header.npts;
    G = memory_calloc64f(6*npts);
    sacio_copyHeader(data->sacGxx[k].header, &data->est[iobs].header);
    cblas_dcopy(data->sacGxx[k].header.npts,
                data->sacGxx[k].data, 1, &G[0*npts], 1); 
    cblas_dcopy(data->sacGyy[k].header.npts,
                data->sacGyy[k].data, 1, &G[1*npts], 1); 
    cblas_dcopy(data->sacGzz[k].header.npts,
                data->sacGzz[k].data, 1, &G[2*npts], 1); 
    cblas_dcopy(data->sacGxy[k].header.npts,
                data->sacGxy[k].data, 1, &G[3*npts], 1); 
    cblas_dcopy(data->sacGxz[k].header.npts,
                data->sacGxz[k].data, 1, &G[4*npts], 1); 
    cblas_dcopy(data->sacGyz[k].header.npts,
                data->sacGyz[k].data, 1, &G[5*npts], 1);
    data->est[iobs].data = memory_calloc64f(npts);
    data->est[iobs].npts = npts;
    strcpy(data->est[iobs].header.kcmpnm, data->data[iobs].header.kcmpnm);
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                npts, 6, 1.0, G, npts,
                &mts[ldm*imtOpt], 1, 0.0, data->est[iobs].data, 1);
    memory_free64f(&G);
    return 0;
}
