#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include "parmt_config.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "parmt_mtsearch.h"
#include "parmt_utils.h"
#include "compearth.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif


static int getNpMax(const int iobs, const struct parmtData_struct data,
                    const struct localMT_struct mtloc);
/*
static int setDataOnG(const int iobs, const int iloc, const int npmax,
                       const struct parmtData_struct data,
                       double *__restrict__ G);
*/
int parmt_computeL1StationMagnitude64f(const int ldm, const int nmt,
                                       const int npts, const int blockSize,
                                       const double *__restrict__ Gxx,
                                       const double *__restrict__ Gyy,
                                       const double *__restrict__ Gzz,
                                       const double *__restrict__ Gxy,
                                       const double *__restrict__ Gxz,
                                       const double *__restrict__ Gyz,
                                       const double *__restrict__ mts,
                                       const double *__restrict__ d,
                                       double *__restrict__ mags);

int parmt_locSearchL164f(const MPI_Comm locComm,
                         const int iobs, const int blockSize,
                         const int nlags, const bool lwantLags,
                         const bool lrescale,
                         struct localMT_struct mtloc,
                         struct parmtData_struct *data,
                         const double *__restrict__ CeInv,
                         double *__restrict__ phi,
                         double *__restrict__ var,
                         int *__restrict__ lags)
{
    const char *fcnm = "parmt_locSearchL164f\0";
    double *G, *d, *phiLoc, *phiWork, *varLoc, *varWork, xrescaleNum;
    int *lagLoc, *lagWork, ierr, iloc, jndx, jloc, k,
        myloc, npmax, nlocProcs, npts;
    const int master = 0;
    ierr = 0;
    myloc =-1;
    if (mtloc.myid == master)
    {
        MPI_Comm_rank(locComm, &myloc);
        MPI_Comm_size(locComm, &nlocProcs);
    }
    MPI_Bcast(&myloc,     1, MPI_INT, master, mtloc.comm);
    MPI_Bcast(&nlocProcs, 1, MPI_INT, master, mtloc.comm);
    if (iobs < 0 || iobs >= data->nobs)
    {
        printf("%s: Invalid observation number %d\n", fcnm, iobs);
        return -1;
    }
    if (mtloc.ldm <= 6)
    {
        printf("%s: Error leading dimension must be at least 6\n", fcnm);
        return -1;
    }
    if (mtloc.mts == NULL || phi == NULL)
    {
        if (mtloc.mts == NULL)
        {
            fprintf(stderr, "%s: Error mts is NULL\n", __func__);
        }
        if (phi == NULL){fprintf(stderr, "%s: Error phi is NULL\n", __func__);}
        return -1;
    }
    npmax = getNpMax(iobs, *data, mtloc);
    if (npmax < 1)
    {
        printf("%s: No data points\n", fcnm);
        return -1;
    }
    // Set space and copy data
    lagLoc = NULL;
    lagWork = NULL;
    varWork = NULL;
    phiWork = NULL;
    G = memory_calloc64f(6*npmax);
    d = memory_calloc64f(npmax);
    varLoc = array_set64f(npmax, DBL_MAX, &ierr);
    phiLoc = memory_calloc64f(mtloc.nmt);
    if (nlocProcs > 1)
    {
        if (mtloc.myid == master)
        {
            phiWork = memory_calloc64f(data->nlocs*mtloc.nmtAll);
            varWork = memory_calloc64f(npmax);
            if (lwantLags)
            {
                lagWork = memory_calloc32i(data->nlocs*mtloc.nmtAll);
            }
        }
        else
        {
            phiWork = memory_calloc64f(1);
            varWork = memory_calloc64f(1);
            lagWork = memory_calloc32i(1);
        }
    }
    else
    {
        phiWork = phi;
        varWork = var;
        lagWork = lags;
    }
    if (lwantLags)
    {
        lagLoc = memory_calloc32i(mtloc.nmt);
    }
    else
    {
        lagLoc = memory_calloc32i(1);
    }
    npts = data->data[iobs].npts;
    cblas_dcopy(npts, data->data[iobs].data, 1, d, 1);
    for (jloc=0; jloc<data->nlocs; jloc=jloc+nlocProcs)
    {
        iloc = jloc + myloc;
        if (iloc >= data->nlocs){goto NEXT_LOCATION;}
        k = iobs*data->nlocs + iloc;
        // Get the Green's functions onto the matrix G
        ierr = parmt_utils_setDataOnG(iobs, iloc, npmax, *data, G); 
        if (ierr != 0)
        {
            printf("%s: Error setting Greens functions\n", fcnm);
            return -1; 
        }
        // Tally the objective function
        ierr = parmt_mtSearchL164f(mtloc.ldm, mtloc.nmt,
                                   npts, blockSize,
                                   nlags, lwantLags,
                                   lrescale,
                                   &G[0*npmax], &G[1*npmax], &G[2*npmax],
                                   &G[3*npmax], &G[4*npmax], &G[5*npmax],
                                   CeInv, mtloc.mts, d,
                                   phiLoc, varLoc, lagLoc);
        if (ierr != 0)
        {
            printf("%s: Error calling XC64f search %d %d\n",
                   fcnm, myloc, mtloc.myid);
            return -1; 
        }
        // Get the results on the master
        jndx = 0;
        if (mtloc.myid == master){jndx = iloc*mtloc.nmtAll;}
        if (mtloc.commSize == 1)
        {
            array_copy64f_work(mtloc.nmt, phiLoc, &phiWork[jndx]);
            if (lwantLags)
            {
                array_copy32i_work(mtloc.nmt, lagLoc, &lagWork[jndx]); 
            }
        }
        else
        {
            MPI_Gatherv(phiLoc, mtloc.nmt, MPI_DOUBLE,
                        &phiWork[jndx], mtloc.nmtProc, mtloc.offset,
                        MPI_DOUBLE, master, mtloc.comm);
            if (lwantLags)
            {
                MPI_Gatherv(lagLoc, mtloc.nmt, MPI_INT,
                            &lagWork[jndx], mtloc.nmtProc, mtloc.offset,
                            MPI_INT, master, mtloc.comm);
            }
        }
NEXT_LOCATION:;
    }
    // Reduce the variance onto varWork
    MPI_Reduce(varLoc, varWork, npts, MPI_DOUBLE, MPI_MIN,
               master, mtloc.comm);
    // Have the location masters reduce their result onto the master
    if (mtloc.myid == master)
    {
        if (nlocProcs > 1)
        {
            MPI_Reduce(phiWork, phi, data->nlocs*mtloc.nmtAll, 
                       MPI_DOUBLE, MPI_SUM, master, locComm);
            MPI_Reduce(varWork, var, npts,
                       MPI_DOUBLE, MPI_MIN, master, locComm);
            if (lwantLags)
            {
                MPI_Reduce(lagWork, lags, data->nlocs*mtloc.nmtAll,
                           MPI_INT, MPI_SUM, master, locComm);
            }
        }
    }
    cblas_dscal(npts, 1.0/(double) mtloc.nmtAll, varLoc, 1);
    // Free memory 
    memory_free32i(&lagLoc);
    memory_free64f(&phiLoc);
    memory_free64f(&varLoc);
    memory_free64f(&d);
    memory_free64f(&G);
    if (nlocProcs > 1)
    {   
        memory_free64f(&phiWork);
        memory_free64f(&varWork);
        memory_free32i(&lagWork);
    }
    phiWork = NULL;
    lagWork = NULL;
    // Block until everyone is done
    if (mtloc.myid == master){MPI_Barrier(locComm);}
    MPI_Barrier(mtloc.comm);
    return 0;
}

int parmt_locSearchXC64f(const MPI_Comm locComm,
                         const int iobs, const int blockSize,
                         const int nlags, const bool lwantLags,
                         struct localMT_struct mtloc,
                         struct parmtData_struct *data,
                         double *__restrict__ phi, int *__restrict__ lags)
{
    const char *fcnm = "parmt_locSearchXC64f\0";
    double *G, *d, *phiLoc, *phiWork;
    int *lagLoc, *lagWork, ierr, iloc, jndx, jloc, k,
        myloc, npmax, nprocs, npts;
    const int master = 0;
    ierr = 0;
    myloc =-1;
    if (mtloc.myid == master)
    {
        MPI_Comm_rank(locComm, &myloc);
        MPI_Comm_size(locComm, &nprocs);
    }
    MPI_Bcast(&myloc, 1, MPI_INT, master, mtloc.comm); 
    MPI_Bcast(&nprocs, 1, MPI_INT, master, mtloc.comm);
    if (iobs < 0 || iobs >= data->nobs)
    {
        printf("%s: Invalid observation number %d\n", fcnm, iobs);
        return -1; 
    }
    if (mtloc.ldm <= 6)
    {
        printf("%s: Error leading dimension must be at least 6\n", fcnm);
        return -1;
    }
    if (mtloc.mts == NULL || phi == NULL)
    {
        if (mtloc.mts == NULL)
        {
            fprintf(stderr, "%s: Error mts is NULL\n", __func__);
        }
        if (phi == NULL){fprintf(stderr, "%s: Error phi is NULL\n", __func__);}
        return -1;
    }
    npmax = getNpMax(iobs, *data, mtloc);
    if (npmax < 1)
    {
        printf("%s: No data points\n", fcnm);
        return -1;
    }
    // Set space and copy data
    G = memory_calloc64f(6*npmax);
    d = memory_calloc64f(npmax);
    phiLoc = memory_calloc64f(mtloc.nmt);
    if (nprocs > 1)
    {
        if (mtloc.myid == master)
        {
            phiWork = memory_calloc64f(data->nlocs*mtloc.nmtAll);
            if (lwantLags)
            {
                lagWork = memory_calloc32i(data->nlocs*mtloc.nmtAll);
            }
        }
        else
        {
             phiWork = memory_calloc64f(1);
             lagWork = memory_calloc32i(1);
        }
    }
    else
    {
        phiWork = phi;
        lagWork = lags;
    }
    if (lwantLags)
    {
        lagLoc = memory_calloc32i(mtloc.nmt);
    }
    else 
    {
        lagLoc = memory_calloc32i(1);
    }
    npts = data->data[iobs].npts;
    cblas_dcopy(npts, data->data[iobs].data, 1, d, 1);
    for (jloc=0; jloc<data->nlocs; jloc=jloc+nprocs)
    {
        iloc = jloc + myloc;
        if (iloc >= data->nlocs){goto NEXT_LOCATION;}
        k = iobs*data->nlocs + iloc;
        // Get the Green's functions onto the matrix G
        ierr = parmt_utils_setDataOnG(iobs, iloc, npmax, *data, G);
        if (ierr != 0)
        {
            printf("%s: Error setting Greens functions\n", fcnm);
            return -1;
        }
        // Tally the objective function
        ierr = parmt_mtSearchXC64f(mtloc.ldm, mtloc.nmt,
                                   npts, blockSize,
                                   nlags, lwantLags,
                                   &G[0*npmax], &G[1*npmax], &G[2*npmax],
                                   &G[3*npmax], &G[4*npmax], &G[5*npmax],
                                   mtloc.mts, d, phiLoc, lagLoc);
        if (ierr != 0)
        {
            printf("%s: Error calling XC64f search %d %d\n",
                   fcnm, myloc, mtloc.myid);
            return -1;
        }
        // Get the results on the master
        jndx = 0;
        if (mtloc.myid == master){jndx = iloc*mtloc.nmtAll;}
        MPI_Gatherv(phiLoc, mtloc.nmt, MPI_DOUBLE,
                    &phiWork[jndx], mtloc.nmtProc, mtloc.offset,
                    MPI_DOUBLE, master, mtloc.comm);
        if (lwantLags)
        {
            MPI_Gatherv(lagLoc, mtloc.nmt, MPI_INT,
                        &lagWork[jndx], mtloc.nmtProc, mtloc.offset,
                        MPI_INT, master, mtloc.comm);
        }
NEXT_LOCATION:;
        MPI_Barrier(mtloc.comm);
    }
    // Have the location masters reduce their result onto the master
    if (mtloc.myid == master && nprocs > 1)
    {
        MPI_Reduce(phiWork, phi, data->nlocs*mtloc.nmtAll, 
                   MPI_DOUBLE, MPI_SUM, master, locComm); 
        if (lwantLags)
        {
            MPI_Reduce(lagWork, lags, data->nlocs*mtloc.nmtAll, 
                       MPI_INT, MPI_SUM, master, locComm);
        }
    }
    // Free memory 
    memory_free32i(&lagLoc);
    memory_free64f(&phiLoc);
    memory_free64f(&d);
    memory_free64f(&G);
    if (nprocs > 1)
    {
        memory_free64f(&phiWork);
        memory_free32i(&lagWork);
    }
    phiWork = NULL;
    lagWork = NULL;
    return 0;
}

static int getNpMax(const int iobs, const struct parmtData_struct data,
                    const struct localMT_struct mtloc)
{
    const char *fcnm = "getNpMax\0";
    int iloc, k, npmax, npts;
    const int master = 0;
    // Get the workspace size
    npmax = 0;
    npts = data.data[iobs].header.npts;
    if (mtloc.myid == master)
    {
        for (iloc=0; iloc<data.nlocs; iloc++)
        {
            k = iobs*data.nlocs + iloc;
            npmax = MAX(npmax, data.data[iobs].header.npts);
            if (data.data[iobs].header.npts != data.sacGxx[k].header.npts ||
                data.data[iobs].header.npts != data.sacGyy[k].header.npts ||
                data.data[iobs].header.npts != data.sacGzz[k].header.npts ||
                data.data[iobs].header.npts != data.sacGxy[k].header.npts ||
                data.data[iobs].header.npts != data.sacGxz[k].header.npts ||
                data.data[iobs].header.npts != data.sacGyz[k].header.npts)
            {
                printf("%s: Inconsistent sizes\n", fcnm);
                npmax = 0;
                break;
            }
        }
    }
    MPI_Bcast(&npmax, 1, MPI_INTEGER, master, mtloc.comm);
    if (npmax < 1)
    {
        printf("%s: Error no data points\n", fcnm);
        return -1;
    }
    return npmax;
}

/*
static int setDataOnG(const int iobs, const int iloc, const int npmax,
                       const struct parmtData_struct data,
                       double *__restrict__ G)
{
    const char *fcnm = "setDataOnG\0";
    int k;
    if (iobs < 0 || iobs >= data.nobs)
    {
        printf("%s: Error invalid observation\n", fcnm);
        return -1;
    }
    if (iloc < 0 || iloc >= data.nlocs)
    {
        printf("%s: Error invalid location\n", fcnm);
        return -1;
    }
    k = iobs*data.nlocs + iloc;
    if (npmax > data.sacGxx[k].npts || npmax > data.sacGyy[k].npts ||
        npmax > data.sacGzz[k].npts || npmax > data.sacGxy[k].npts ||
        npmax > data.sacGxz[k].npts || npmax > data.sacGyz[k].npts)
    {
        printf("%s: npmax is too small\n", fcnm);
        return -1;
    }
    if (data.sacGxx[k].data == NULL || data.sacGyy[k].data == NULL ||
        data.sacGzz[k].data == NULL || data.sacGxy[k].data == NULL ||
        data.sacGxz[k].data == NULL || data.sacGyz[k].data == NULL)
    {
        printf("%s: Greens function is null\n", fcnm);
        return -1; 
    }
    // Get the Green's functions onto the matrix G
    array_zeros64f_work(6*npmax, G); 
    cblas_dcopy(data.sacGxx[k].header.npts,
                data.sacGxx[k].data, 1, &G[0*npmax], 1); 
    cblas_dcopy(data.sacGyy[k].header.npts,
                data.sacGyy[k].data, 1, &G[1*npmax], 1); 
    cblas_dcopy(data.sacGzz[k].header.npts,
                data.sacGzz[k].data, 1, &G[2*npmax], 1); 
    cblas_dcopy(data.sacGxy[k].header.npts,
                data.sacGxy[k].data, 1, &G[3*npmax], 1); 
    cblas_dcopy(data.sacGxz[k].header.npts,
                data.sacGxz[k].data, 1, &G[4*npmax], 1); 
    cblas_dcopy(data.sacGyz[k].header.npts,
                data.sacGyz[k].data, 1, &G[5*npmax], 1);
    return 0;
}
*/
