#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#include "parmt_mtsearch.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif 
#include "parmt_utils.h"
#include "cps_mpe.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

int parmt_obsSearch64f(const MPI_Comm globalComm,
                       const MPI_Comm obsComm, const MPI_Comm locComm,
                       const bool linObsComm, const bool linLocComm,
                       const struct parmtGeneralParms_struct parms,
                       struct parmtData_struct data,
                       struct localMT_struct mtloc,
                       double *__restrict__ phi)
{
    MPI_Win counterWin;
    const char *archiveName = NULL;
    double lagTime, *CeInv, *phiWork, *phiAll, *var, *varAll, *varOut, *wts,
           defaultWeight, wtNorm;
    int *lagsWork, *nlags, ierr, iobs, myid, nobsProcs,
        npmax, value, value0;
    bool ldefault, lwantLags;
    const double half = 0.5;
    const int master = 0;
    //------------------------------------------------------------------------//
    //
    // Have the root process determine lag information
    ierr = 0;
    MPI_Comm_rank(globalComm, &myid);
    nlags = memory_calloc32i(data.nobs);
    nobsProcs = 0;
    lagsWork = NULL;
    phiWork = NULL;
    phiAll = NULL;
    CeInv = NULL;
    var = NULL;
    varAll = NULL;
    varOut = NULL;
    wts = NULL;
    npmax = 0;
    lwantLags = false;
    if (myid == master)
    {
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            nlags[iobs] = 0;
            if (parms.lwantLags)
            {
                //lagTime = data.data[iobs].header.user0;
                //if (lagTime < zero){lagTime = parms.defaultMaxLagTime;}
                lagTime = parmt_utils_getLagTime(data.data[iobs], 
                                                 parms.defaultMaxLagTime,
                                                 &ldefault);
                nlags[iobs] = (int)
                              (lagTime/data.data[iobs].header.delta + half);
            }
            lwantLags = false;
            if (nlags[iobs] > 0){lwantLags = true;}
            npmax = MAX(npmax, data.data[iobs].npts);
            if (parms.lwantLags)
            {
                printf("%s: Number of lags for waveform %d is %d\n",
                       __func__, iobs+1, nlags[iobs]);
            }
        }
    }
    MPI_Bcast(nlags, data.nobs, MPI_INT,    master, globalComm);
    MPI_Bcast(&npmax,        1, MPI_INT,    master, globalComm);
    MPI_Bcast(&lwantLags,    1, MPI_C_BOOL, master, globalComm);
    // Check for failures
    if (npmax < 1 || mtloc.nmtAll < 1 || parms.objFnType < 1 ||
        parms.objFnType > 3)
    {
        if (npmax < 1){printf("%s: No data points\n", __func__);}
        if (mtloc.nmtAll < 1){printf("%s: Empty mts\n", __func__);}
        if (parms.objFnType < 1 || parms.objFnType > 3)
        {
            printf("%s: Invalid objfn function\n", __func__);
        }
        return -1;
    }
    if (lwantLags)
    {
        if (mtloc.myid == master)
        {
            lagsWork = memory_calloc32i(data.nlocs*mtloc.nmtAll);
        }
        else
        {
            lagsWork = memory_calloc32i(1);
        }
    }
    // Set space
    var = memory_calloc64f(npmax);
    if (linObsComm){MPI_Comm_size(obsComm, &nobsProcs);}
    if (linObsComm)
    {
        if (nobsProcs > 1)
        {
            phiAll = memory_calloc64f(data.nlocs*mtloc.nmtAll);
        }
        else
        {
            if (phi == NULL)
            {
                printf("%s: phi not set on master - segfault\n", __func__);
                return -1; 
            }
            phiAll = phi;
            array_zeros64f_work(data.nlocs*mtloc.nmtAll, phiAll);
        }
        phiWork = memory_calloc64f(data.nlocs*mtloc.nmtAll);
        varAll = memory_calloc64f(data.nobs*npmax);
        defaultWeight = 1.0; ///(double) data.nobs;
        wts = memory_calloc64f(data.nobs); //array_set64f(data.nobs, defaultWeight, &ierr);
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            wts[iobs] = parmt_utils_getWeight(data.data[iobs],
                                              defaultWeight,
                                              &ldefault);
        }
        // Normalize so weighted sum is unity
        wtNorm = 1.0/array_sum64f(data.nobs, wts, &ierr);
        // And rescale the weights
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            wts[iobs] = wts[iobs]*wtNorm; 
        }
    }
    else
    {
        phiWork = memory_calloc64f(1);
        phiAll = memory_calloc64f(1);
    }
    // Set up the MPI RMA window
    value = 0;
    value0 = 0;
    if (linObsComm)
    {
        MPE_setKeyval(MPI_KEYVAL_INVALID);
        MPE_counterCreate(obsComm, 1, &counterWin);
    }
    MPI_Barrier(globalComm);
    // Work-steal the observations to accomodate different sizes/lags
    while (true)
    {
        // Update the counter
        if (linObsComm){MPE_counterNextval(counterWin, 0, &value);}
        // Send the counter to the locations commuincator
        if (linLocComm){MPI_Bcast(&value, 1, MPI_INT, master, locComm);}
        // Trickle the value down to the moment tensor communicator
        MPI_Bcast(&value, 1, MPI_INT, master, mtloc.comm);
        if (value > data.nobs){break;}  // We are done
        if (value == value0){continue;} // Already did this
        value0 = value;   // Update the old value
        iobs = value - 1; // Map to C based numbering to indicate observation
        if (linObsComm)
        {
            printf("%s: Processing observation %d with %d points\n", __func__,
                   iobs+1, data.data[iobs].npts);
        }
        // Figure out the lags
        lwantLags = false;
        if (nlags[iobs] > 0){lwantLags = true;}
        if (parms.objFnType == 1)
        {
            ierr = parmt_locSearchL164f(locComm,
                                        iobs, parms.blockSize,
                                        nlags[iobs], lwantLags,
                                        parms.lrescale,
                                        mtloc, &data,
                                        CeInv, phiWork, var, lagsWork);
            if (ierr != 0)
            {
                printf("%s: L1 search error for observation %d\n",
                       __func__, iobs + 1);
                return -1;
            }
            // Put the variance into the global buffer
            if (linObsComm)
            {
                array_copy64f_work(npmax, var, &varAll[iobs*npmax]);
            }
        }
        else if (parms.objFnType == 2)
        {
            ierr = parmt_locSearchXC64f(locComm, iobs, parms.blockSize,
                                        nlags[iobs], lwantLags,
                                        mtloc, &data, phiWork, lagsWork);
            if (ierr != 0)
            {
                printf("%s: XC search error for observation %d\n",
                       __func__, iobs + 1);
                return -1;
            }
        }
        else if (parms.objFnType == 3)
        {
            printf("%s: need to write this function\n", __func__);
            return -1;
        }
        // Scale by number of observations and stack
        if (linObsComm)
        {
            //printf("%d %f\n", iobs, wts[iobs]);
            cblas_daxpy(data.nlocs*mtloc.nmtAll, wts[iobs], //1.0/(double) data.nobs,
                        phiWork, 1, phiAll, 1);
        }
    }
    // Reduce the objective function
    if (linObsComm)
    {
        if (nobsProcs > 1)
        {
            MPI_Reduce(phiAll, phi, data.nlocs*mtloc.nmtAll, MPI_DOUBLE,
                       MPI_SUM, master, obsComm);
        }
        MPI_Win_free(&counterWin);
    }
    MPI_Barrier(mtloc.comm);
    // Release resources
    if (linObsComm)
    {
        if (nobsProcs > 1){memory_free64f(&phiAll);}
        // Reduce the variance
        varOut = NULL;
        if (myid == master){varOut = memory_calloc64f(npmax*data.nobs);}
        MPI_Reduce(varAll, varOut, npmax*data.nobs, MPI_DOUBLE,
                   MPI_SUM, master, obsComm);
        // Have the master write it
        if (myid == master)
        {
            if (parms.programID == PARMT_ID)
            {
                archiveName = parms.parmtArchive;
            }
            else if (parms.programID == POLARMT_ID)
            {
                archiveName = parms.polarmtArchive;
            }
            for (iobs=0; iobs<data.nobs; iobs++)
            {
                parmt_io_writeVarianceForWaveform64f(
                   //parms.resultsDir, parms.projnm, parms.resultsFileSuffix,
                   archiveName,
                   iobs,
                   data.data[iobs].npts, &varOut[iobs*npmax]);
            }
            archiveName = NULL;
        }
        memory_free64f(&varAll);
        memory_free64f(&varOut);
    }
    else
    {
        memory_free64f(&phiAll);
    }
    memory_free64f(&phiWork);
    memory_free32i(&nlags);
    memory_free32i(&lagsWork);
    memory_free64f(&CeInv);
    memory_free64f(&wts);
    phiAll = NULL;
    MPI_Barrier(globalComm);
    return 0;
}
