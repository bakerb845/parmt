#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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
    const char *fcnm = "parmt_obsSearch64f\0";
    MPI_Win counterWin;
    double lagTime, *CeInv, *phiWork, *phiAll, *var;
    int *lagsWork, *nlags, ierr, iobs, jobs, myid, nobsProcs,
        npmax, value, value0;
    bool lwantLags;
    const double zero = 0.0;
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
    npmax = 0;
    lwantLags = false;
    if (myid == master)
    {
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            nlags[iobs] = 0;
            if (parms.lwantLags)
            {
                lagTime = data.data[iobs].header.user0;
                if (lagTime < zero){lagTime = parms.defaultMaxLagTime;}
                nlags[iobs] = (int)
                              (lagTime/data.data[iobs].header.delta + half);
            }
            lwantLags = false;
            if (nlags[iobs] > 0){lwantLags = true;}
            npmax = MAX(npmax, data.data[iobs].npts);
        }
    }
    MPI_Bcast(nlags, data.nobs, MPI_INT,    master, globalComm);
    MPI_Bcast(&npmax,        1, MPI_INT,    master, globalComm);
    MPI_Bcast(&lwantLags,    1, MPI_C_BOOL, master, globalComm);
    // Check for failures
    if (npmax < 1 || mtloc.nmtAll < 1 || parms.objFnType < 1 ||
        parms.objFnType > 3)
    {
        if (npmax < 1){printf("%s: No data points\n", fcnm);}
        if (mtloc.nmtAll < 1){printf("%s: Empty mts\n", fcnm);}
        if (parms.objFnType < 1 || parms.objFnType > 3)
        {
            printf("%s: Invalid objfn function\n", fcnm);
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
                printf("%s: phi not set on master - segfault\n", fcnm);
            }
            phiAll = phi;
            array_zeros64f_work(data.nlocs*mtloc.nmtAll, phiAll);
        }
        phiWork = memory_calloc64f(data.nlocs*mtloc.nmtAll);
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
        MPI_Comm_size(obsComm, &nobsProcs);
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
        // Figure out the lags
        lwantLags = false;
        if (nlags[iobs] > 0){lwantLags = true;}
        if (parms.objFnType == 1)
        {
            ierr = parmt_locSearchL164f(locComm,
                                        iobs, parms.blockSize,
                                        nlags[iobs], lwantLags,
                                        mtloc, &data,
                                        CeInv, phiWork, var, lagsWork);
            if (ierr != 0)
            {
                printf("%s: L1 search error for observation %d\n",
                       fcnm, iobs + 1);
                return -1;
            }
            if (myid == master)
            {
                parmt_io_writeVarianceForWaveform64f(
                   parms.resultsDir, parms.projnm, parms.resultsFileSuffix,
                   iobs,
                   data.data[iobs].npts, var);
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
                       fcnm, iobs + 1);
                return -1;
            }
        }
        else if (parms.objFnType == 3)
        {
            printf("%s: need to write this function\n", fcnm);
            return -1;
        }
        // Stack objective function
        if (linObsComm)
        {
            cblas_daxpy(data.nlocs*mtloc.nmtAll, 1.0/(double) data.nobs,
                        phiWork, 1, phiAll, 1);
        }
    }
    // Reduce the objective function
    if (linObsComm)
    {
        if (nobsProcs > 1)
        {
            MPI_Reduce(phiAll, phi, data.nobs*mtloc.nmtAll, MPI_DOUBLE,
                       MPI_SUM, master, locComm);
        }
    }
    // Release resources
    if (linObsComm){MPI_Win_free(&counterWin);}
    if (linObsComm)
    {
        if (nobsProcs > 1){memory_free64f(&phiAll);}
        phiAll = NULL;
    }
    else
    {
        memory_free64f(&phiAll);
    }
    memory_free64f(&phiWork);
    memory_free32i(&nlags);
    memory_free32i(&lagsWork);
    memory_free64f(&CeInv);
    // Block until everyone is done
    MPI_Barrier(globalComm);
    return 0;
}
