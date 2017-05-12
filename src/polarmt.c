#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <limits.h>
#include <mpi.h>
#include <omp.h>
#include "parmt_polarity.h"
#ifdef PARMT_USE_INTEL
#include <mkl.h>
#include <ipp.h>
#endif
#include "compearth.h"
#include "iscl/iscl/iscl.h"
#include "iscl/memory/memory.h"

#define PROGRAM_NAME "polarmt"

static int freePolarityData(struct polarityData_struct *polarityData);
int parmt_freeData(struct parmtData_struct *data);
int parmt_freeLocalMTs(struct localMT_struct *mtloc);
static int parseArguments(int argc, char *argv[],
                          const int nprocs,
                          int *npInLocGroups,
                          int *npInMTGroups,
                          char iniFile[PATH_MAX]);
static void printUsage(void);

int main(int argc, char *argv[])
{
    char iniFile[PATH_MAX], fileSuffix[PATH_MAX];
    MPI_Comm obsComm, locComm, mtComm;
    struct parmtData_struct data;
    struct parmtGeneralParms_struct parms;
    struct parmtPolarityParms_struct polarityParms;
    struct parmtMtSearchParms_struct mtsearch;
    struct polarityData_struct polarityData;
    struct localMT_struct mtloc;
    double *betas, *deps, *gammas, *h, *kappas, *M0s, *phi, *sigmas,
           *thetas, *u, *v, du, dv, dh, dk, ds,
           hLower, hUpper, t0, t1, uLower, uUpper, vLower, vUpper;
    bool linLocComm, linMTComm, linObsComm;
    int i, ierr, myid, nmt, nmtAll, npInLocGroups, npInMTGroups, nprocs, provided;
    const int npInObsGroups = 1;
    const int master = 0;
    const int ldm = 8;
    const double M0unit[1] = {1.0/sqrt(2.0)};
    //------------------------------------------------------------------------//
    //
    // initialize MPI, ISCL, and handle threading issues
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    iscl_init();
#ifdef PARMT_USE_INTEL
    //omp_set_num_threads(1);
    mkl_set_num_threads(1);
    ippSetNumThreads(1);
#endif
    // initialize
    ierr = 0;
    betas = NULL;
    gammas = NULL;
    h = NULL;
    kappas = NULL;
    M0s = NULL;
    phi = NULL;
    sigmas = NULL;
    thetas = NULL;
    u = NULL;
    v = NULL;
    memset(&parms, 0, sizeof(struct parmtGeneralParms_struct));
    memset(&data, 0, sizeof(struct parmtData_struct));
    memset(&mtsearch, 0, sizeof(struct parmtMtSearchParms_struct));
    memset(&polarityParms, 0, sizeof(struct parmtPolarityParms_struct));
    memset(&polarityData, 0, sizeof(struct polarityData_struct));
    if (myid == master)
    {
        ierr = parseArguments(argc, argv, nprocs,
                              &npInLocGroups, &npInMTGroups, iniFile);
        if (ierr != 0){goto INIT_ERROR;}
        ierr = 0;
        ierr += parmt_utils_readGeneralParms(iniFile, &parms);
        ierr += parmt_utils_readMtSearch(iniFile, &mtsearch);
        ierr += parmt_utils_readPolarityParms(iniFile, &polarityParms);
        if (ierr != 0)
        {
            printf("%s: Error reading ini file\n", PROGRAM_NAME);
            goto INIT_ERROR;
        }
        printf("%s: Process management:\n", PROGRAM_NAME);
        printf("         Number of processes in observation groups %d\n",
               npInObsGroups);
        printf("         Number of processes in location groups %d\n",
               npInLocGroups);
        printf("         Number of processes in moment tensor groups %d\n",
               npInMTGroups);
        // Load the observations
        printf("%s: Loading waveforms...\n", PROGRAM_NAME);
        ierr = utils_dataArchive_readAllWaveforms(parms.dataFile, &data);
        if (ierr != 0)
        {
            printf("%s: Error reading waveforms\n", PROGRAM_NAME);
            goto INIT_ERROR;
        } 
    }
INIT_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INTEGER, master, MPI_COMM_WORLD);
    if (ierr != 0){goto FINISH;}
    MPI_Bcast(&npInLocGroups, 1, MPI_INT, master, MPI_COMM_WORLD);
    MPI_Bcast(&npInMTGroups,  1, MPI_INT, master, MPI_COMM_WORLD);
    // split the communicators
    ierr = parmt_splitComm(MPI_COMM_WORLD,
                           npInObsGroups, npInLocGroups, npInMTGroups,
                           &linObsComm, &obsComm,
                           &linLocComm, &locComm,
                           &linMTComm,  &mtComm);
    if (ierr != 0)
    {
        printf("%s: Fatal error splitting communicators - have a nice day\n",
               PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    // Send the inputs to everyone
    if (myid == master)
    {
        printf("%s: Broadcasting parameters...\n", PROGRAM_NAME);
    }
    ierr = parmt_broadcast_mtSearchParms(&mtsearch, master, MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("%s: Fatal error broadcasting parameters\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    ierr = parmt_broadcast_generalParms(&parms, master, MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("%s: Fatal error broadcasting general parameters\n",
               PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    } 
    ierr = parmt_broadcast_parmtPolarityParms(&polarityParms,
                                              master, MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("%s: Fatal error broadcasting polarity parms\n",
               PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    ierr = parmt_broadcast_data(&data, master, MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("%s: Fatal error broadcasting waveforms\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    // Discretize the moment tensor space in (u, v, h) space
    ierr = parmt_discretizeCells64f(
               mtsearch.nb, mtsearch.betaLower,  mtsearch.betaUpper,
               mtsearch.ng, mtsearch.gammaLower, mtsearch.gammaUpper,
               mtsearch.nk, mtsearch.kappaLower, mtsearch.kappaUpper,
               mtsearch.ns, mtsearch.sigmaLower, mtsearch.sigmaUpper,
               mtsearch.nt, mtsearch.thetaLower, mtsearch.thetaUpper,
               mtsearch.nm, mtsearch.m0Lower,    mtsearch.m0Upper,
               &betas,  &gammas, &kappas,
               &sigmas, &thetas, &M0s);
    if (ierr != 0)
    {
        printf("%s: Error discretizing MT space\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
/*
    compearth_beta2u(1, &mtsearch.betaLower, &uLower);
    compearth_beta2u(1, &mtsearch.betaUpper, &uUpper);
    compearth_gamma2v(1, &mtsearch.gammaLower, &vLower);
    compearth_gamma2v(1, &mtsearch.gammaUpper, &vUpper);
    compearth_theta2h(1, &mtsearch.thetaLower, &hLower);
    compearth_theta2h(1, &mtsearch.thetaUpper, &hUpper);
    du = (uUpper - uLower)/(double) mtsearch.nb;
    dv = (vUpper - vLower)/(double) mtsearch.ng;
    dh =-(hUpper - hLower)/(double) mtsearch.nt;
    dk = (mtsearch.kappaUpper - mtsearch.kappaLower)/(double) mtsearch.nk;
    ds = (mtsearch.sigmaUpper - mtsearch.sigmaLower)/(double) mtsearch.ns;
    u = array_linspace64f(uLower+du/2, uUpper-du/2, mtsearch.nb, &ierr);
    v = array_linspace64f(vLower+dv/2, vUpper-dv/2, mtsearch.ng, &ierr);
    h = array_linspace64f(hLower-dh/2, hUpper+dh/2, mtsearch.nt, &ierr);
    betas  = memory_calloc64f(mtsearch.nb);
    gammas = memory_calloc64f(mtsearch.ng);
    thetas = memory_calloc64f(mtsearch.nt);
    compearth_u2beta(mtsearch.nb, 20, 2, u, 1.e-6, betas);
    compearth_v2gamma(mtsearch.ng, v, gammas);
    compearth_h2theta(mtsearch.nt, h, thetas);
    M0s    = array_linspace64f(mtsearch.m0Lower,    mtsearch.m0Upper,
                               mtsearch.nm, &ierr);
    kappas = array_linspace64f(mtsearch.kappaLower+dk/2,
                               mtsearch.kappaUpper-dk/2,
                               mtsearch.nk, &ierr);
    sigmas = array_linspace64f(mtsearch.sigmaLower+ds/2,
                               mtsearch.sigmaUpper-ds/2,
                               mtsearch.ns, &ierr);
*/
    // avoid an annoying warning
    for (i=0; i<mtsearch.ns; i++)
    {
        if (fabs(sigmas[i]) < 1.e-7){sigmas[i] = 1.e-7;}
    }
    // discretize the space - note this is a polarity only estimation so M
    // is meaningless and will be incorporated later
    nmtAll = mtsearch.ng*mtsearch.nb*mtsearch.nk
            *mtsearch.ns*mtsearch.nt*mtsearch.nm;
    nmt = mtsearch.ng*mtsearch.nb*mtsearch.nk
         *mtsearch.ns*mtsearch.nt;
    if (myid == master)
    {
        printf("%s: Discretizing moment tensor...\n", PROGRAM_NAME);
    }
    t0 = MPI_Wtime();
    ierr = parmt_discretizeMT64f_MPI(mtComm,
                                     mtsearch.ng, gammas,
                                     mtsearch.nb, betas,
                                     1, M0unit, //mtsearch.nm, M0s,
                                     mtsearch.nk, kappas,
                                     mtsearch.nt, thetas,
                                     mtsearch.ns, sigmas,
                                     ldm, &mtloc);
    if (ierr != 0)
    {
        printf("%s: Error discretizing moment tensor on process %d\n",
               PROGRAM_NAME, myid);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if (myid == master)
    {
        printf("%s: Discretization time %6.2f (seconds)\n",
               PROGRAM_NAME, t1 - t0);
    }
    // Compute the Green's functions
    if (myid == master)
    {
        printf("%s: Computing Green's functions...\n", PROGRAM_NAME);
    }
    t0 = MPI_Wtime();
    ierr = parmt_polarity_computeTTimesGreens(MPI_COMM_WORLD, polarityParms,
                                              data, &polarityData);
    if (ierr != 0)
    {
        printf("%s: Error computing polarity Green's functions\n",
               PROGRAM_NAME);
        ierr = 1;
        goto FINISH;
    }
    t1 = MPI_Wtime();
    if (myid == master)
    {
        printf("%s: Green's functions computation time %6.2f (seconds)\n",
               PROGRAM_NAME, t1 - t0);
    }
    // Perform the grid search
    if (myid == master)
    {
        printf("%s: Performing grid search...\n", PROGRAM_NAME);
    }
    phi = memory_calloc64f(data.nlocs*mtloc.nmtAll);
    t0 = MPI_Wtime();
    ierr = polarity_performLocationSearch64f(
               locComm,
               parms.blockSize,
               polarityData,
               mtloc, 
               phi);
    if (ierr != 0)
    {
        ierr = 1;
        printf("%s: Failed computing grid search on process %d\n",
               PROGRAM_NAME, myid);
        goto FINISH;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == master)
    {
        t1 = MPI_Wtime();
        printf("%s: Grid search time: %6.2f (seconds)\n", PROGRAM_NAME, t1 - t0);
    }
    t0 = MPI_Wtime();
    if (myid == master)
    {
        printf("%s: Initializing output file...\n", PROGRAM_NAME);
        deps = memory_calloc64f(data.nlocs); 
        for (i=0; i<data.nlocs; i++)
        {
            deps[i] = data.sacGxx[i].header.evdp;
        }
        ierr = parmt_io_createObjfnArchive64f(parms.resultsDir, parms.projnm,
                                              parms.polarityFileSuffix,
                                              data.nobs,
                                              data.nlocs, deps,
                                              1, M0unit, //mtsearch.nm, M0s,
                                              mtsearch.nb, betas,
                                              mtsearch.ng, gammas,
                                              mtsearch.nk, kappas,
                                              mtsearch.ns, sigmas,
                                              mtsearch.nt, thetas);
        if (ierr != 0)
        {
            printf("%s: Error initializing file\n", PROGRAM_NAME);
        }
        ierr = parmt_io_writeObjectiveFunction64f(
                   parms.resultsDir, parms.projnm, parms.polarityFileSuffix,
                   mtloc.nmtAll, phi);
        if (ierr != 0)
        {
            printf("%s: Error writing objective function\n", PROGRAM_NAME);
            goto FINISH;
        }

        memory_free64f(&deps); 
    }
    // Destroy communicators, free memory, and shut down
FINISH:;
    if (ierr != 0)
    {
        printf("%s: An error occurred on process %d\n", PROGRAM_NAME, myid);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    if (linObsComm){MPI_Comm_free(&obsComm);}
    if (linLocComm){MPI_Comm_free(&locComm);}
    if (linMTComm){MPI_Comm_free(&mtComm);}
    parmt_freeData(&data);
    freePolarityData(&polarityData);
    parmt_freeLocalMTs(&mtloc);
    memory_free64f(&u);
    memory_free64f(&v);
    memory_free64f(&h);
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&M0s);
    memory_free64f(&phi);
    memory_free64f(&sigmas);
    memory_free64f(&thetas);
    iscl_finalize();
    MPI_Finalize();
    return 0;
}

static int freePolarityData(struct polarityData_struct *polarityData)
{
    memory_free64f(&polarityData->Gxx);
    memory_free64f(&polarityData->Gyy);
    memory_free64f(&polarityData->Gzz);
    memory_free64f(&polarityData->Gxy);
    memory_free64f(&polarityData->Gxz);
    memory_free64f(&polarityData->Gyz);
    memory_free64f(&polarityData->polarity);
    memory_free64f(&polarityData->wts);
    memory_free32i(&polarityData->obsMap);
    memset(polarityData, 0, sizeof(struct polarityData_struct));
    return 0;
}

int parmt_freeLocalMTs(struct localMT_struct *mtloc)
{
    memory_free64f(&mtloc->mts);
    memory_free32i(&mtloc->offset);
    memset(mtloc, 0, sizeof(struct localMT_struct));
    return 0;
}

int parmt_freeData(struct parmtData_struct *data)
{
    int i;
    if (data->nobs > 0 && data->nlocs > 0 &&
        data->sacGxx != NULL && data->sacGyy != NULL && data->sacGzz != NULL &&
        data->sacGxy != NULL && data->sacGxz != NULL && data->sacGyz != NULL)
    {   
        for (i=0; i<data->nobs*data->nlocs; i++)
        {
            sacio_freeData(&data->sacGxx[i]);
            sacio_freeData(&data->sacGyy[i]);
            sacio_freeData(&data->sacGzz[i]); 
            sacio_freeData(&data->sacGxy[i]);
            sacio_freeData(&data->sacGxz[i]);
            sacio_freeData(&data->sacGyz[i]);
        }
        free(data->sacGxx);
        free(data->sacGyy);
        free(data->sacGzz); 
        free(data->sacGxy);
        free(data->sacGxz);
        free(data->sacGyz);
    }   
    if (data->nobs > 0 && data->data != NULL)
    {   
        for (i=0; i<data->nobs; i++)
        {
            sacio_freeData(&data->data[i]);
        }
        free(data->data);
    }   
    if (data->nobs > 0 && data->est != NULL)
    {   
        free(data->est);
    }   
    memset(data, 0, sizeof(struct parmtData_struct));
    return 0;
}

//============================================================================//
static int parseArguments(int argc, char *argv[],
                          const int nprocs,
                          int *npInLocGroups,
                          int *npInMTGroups,
                          char iniFile[PATH_MAX])
{
    bool linFile;
    int prod;
    linFile = false;
    *npInLocGroups = 1;
    *npInMTGroups = 1;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"ini_file", required_argument, 0, 'i'},
            {"mtGroupSize",  required_argument, 0, 'm'},
            {"locGroupSize", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:m:l:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            linFile = true;
        }
        else if (c == 'l' && nprocs > 1)
        {
            *npInLocGroups = atoi(optarg);
        }
        else if (c == 'm' && nprocs > 1)
        {
            *npInMTGroups = atoi(optarg);
        }
        else if (c == 'h' || c == '?')
        {
            printUsage();
            return -2;
        }
        else
        {
            printf("%s: Unknown options: %s\n",
                   PROGRAM_NAME, argv[optionIndex]);
        }
    }
    if (!linFile)
    {
        printf("%s: Error must specify ini file\n\n", PROGRAM_NAME);
        printUsage();
        return -1;
    }
    prod = (*npInLocGroups)*(*npInMTGroups);
    if (prod != nprocs)
    {
        if (prod != 1)
        {
            printf("%s: Invalid process layout - defaulting to (%d,%d)\n",
                   PROGRAM_NAME, 1, nprocs);
        }
        *npInLocGroups = 1;
        *npInMTGroups = nprocs;
    }
    return 0;
}
static void printUsage(void)
{
    printf("Usage:\n   parmt -i input_file\n\n");
    printf("Required arguments:\n");
    printf("   -i input_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("   -m number of processes in moment-tensor groups\n");
    printf("   -l number of processes in location groups\n");
    printf("   -h displays this message\n");
    printf("   Note if m*l must equal the number of processes\n");
    printf("   If they are not set the program will set -m=nprocs\n");
    return;
}

