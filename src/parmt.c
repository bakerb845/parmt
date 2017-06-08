#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <mpi.h>
#include "parmt_config.h"
#ifdef PARMT_USE_INTEL
#include <ippcore.h>
#include <mkl.h>
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "compearth.h"
#include "parmt_mtsearch.h"
#include "parmt_postProcess.h"
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

#define PROGRAM_NAME "parmt"
static int parseArguments(int argc, char *argv[],
                          const int nprocs,
                          int *npInObsGroups,
                          int *npInLocGroups,
                          int *npInMTGroups,
                          char iniFile[PATH_MAX]);
static void printUsage(void);
int parmt_freeData(struct parmtData_struct *data);
int parmt_freeLocalMTs(struct localMT_struct *mts);
/*!
 *
 */
int main(int argc, char *argv[])
{
    char iniFile[PATH_MAX];
    double *betas, *h, *kappas, *gammas, *M0s, *sigmas, *thetas, *u, *v;
    double *deps, *luneMPDF, *phi, t0, t1, du, dv, dh, dk, ds;
    int64_t ngridSearch;
    double hLower, hUpper, uLower, uUpper, vLower, vUpper, xnorm;
    int *lags, i, ierr, iobs, myid, npInLocGroups, nmt, npInMTGroups,
        npInObsGroups, nprocs, provided;
    int ix, iy, k;
    MPI_Comm mtComm, locComm, obsComm;
    bool linMTComm, linLocComm, linObsComm;
    const int master = 0;
    const int ldm = 8;
struct parmtGeneralParms_struct parms;
struct parmtData_struct data;
struct parmtMtSearchParms_struct mtsearch;
    struct localMT_struct mtloc;
    //------------------------------------------------------------------------//
    //
    // initialize MPI, ISCL, and handle threading issues
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    iscl_init();
#ifdef PARMT_USE_INTEL
    //omp_set_num_threads(2);
    mkl_set_num_threads(1);
    ippSetNumThreads(1);
#endif
    // initialize
    phi = NULL;
    M0s = NULL;
    betas = NULL;
    gammas = NULL;
    sigmas = NULL;
    kappas = NULL;
    thetas = NULL;
    h = NULL;
    u = NULL;
    v = NULL;
    lags = NULL;
    memset(&parms, 0, sizeof(struct parmtGeneralParms_struct));
    memset(&data, 0, sizeof(struct parmtData_struct));
    memset(&mtsearch, 0, sizeof(struct parmtMtSearchParms_struct));
    memset(&mtloc, 0, sizeof(struct localMT_struct));
    if (myid == master)
    {
        ierr = parseArguments(argc, argv, nprocs,
                              &npInObsGroups, &npInLocGroups, &npInMTGroups,
                              iniFile);
        if (ierr != 0){goto INIT_ERROR;}
        // read the ini file
        printf("%s: Parsing ini file %s...\n", PROGRAM_NAME, iniFile);
        ierr = 0;
        ierr += parmt_utils_readGeneralParms(iniFile, &parms);
        strcpy(parms.programName, PROGRAM_NAME);
        parms.programID = PARMT_ID;
        ierr += parmt_utils_readMtSearch(iniFile, &mtsearch);
        if (ierr != 0)
        {
            printf("%s: Error reading ini file\n", PROGRAM_NAME);
            goto INIT_ERROR;
        }
        printf("%s: Process management:\n", PROGRAM_NAME);
        printf("        Number of processes in observation groups %d\n",
               npInObsGroups);
        printf("        Number of processes in location groups %d\n",
               npInLocGroups);
        printf("        Number of processes in moment tensor groups %d\n",
               npInMTGroups);
//        npInObsGroups = 1;
//        npInLocGroups = 1;
//        npInMTGroups = nprocs;
printf("verify greens fns are correct\n");
/*
        double pAxis[3], nAxis[3], tAxis[3], beta, gamma, u,v;
        u= 3.0*M_PI/8.0;
        v=-1.0/9.0;
        compearth_u2beta(1, 10, 2, &u, 1.e-10, &beta);
        compearth_v2gamma(1, &v, &gamma);
        postprocess_tt2tnp(beta, gamma, 4.0*M_PI/5.0, -M_PI/2.0, 0.723,
                           pAxis, nAxis, tAxis);
        int nxp = 51;
        double *xw1 = memory_calloc64f(nxp*nxp);
        double *yw1 = memory_calloc64f(nxp*nxp);
        int8_t *pn1 = (int8_t *) calloc(nxp*nxp, sizeof(int8_t));
        const double xc = 1.5; const double yc = 1.5; const double rad = 1.0;
        // Compute the corresponding moment tensor
        double M0loc = 1.0/sqrt(2.0), Muse[6], lamUse[9], Uuse[9];
        double gammaDeg = gamma*180.0/M_PI;
        double deltaDeg = (M_PI_2 - beta)*180.0/M_PI;
        double kappaDeg = (4.0*M_PI/5.0)*180.0/M_PI;
        double sigmaDeg = (-M_PI_2)*180.0/M_PI;
        double thetaDeg = 0.723*180.0/M_PI;
        compearth_tt2cmt(gammaDeg, deltaDeg, M0loc, 
                         kappaDeg, thetaDeg,
                         sigmaDeg,
                         Muse, lamUse, Uuse);
        printf("mt = [%f,%f,%f,%f,%f,%f]\n", Muse[0], Muse[1], Muse[2], Muse[3], Muse[4], Muse[5]);
printf("drawing\n");
for (int i=0; i<10000; i++)
{
        postprocess_tnp2beachballPolarity(nxp, xc, yc, rad,
                                          pAxis, nAxis, tAxis,
                                          xw1, yw1, pn1);
}
        FILE *fwork = fopen("depmag/beachball.txt", "w");
        for (int iy=0; iy<nxp; iy++)
        {
            for (int ix=0; ix<nxp; ix++)
            {
                int k = iy*nxp + ix;
                fprintf(fwork, "%e %e %d\n", xw1[k], yw1[k], pn1[k]);
            }
            fprintf(fwork, "\n");
        }
        fclose(fwork);
        memory_free64f(&xw1);
        memory_free64f(&yw1);
        free(pn1);
printf("back\n");
//getchar();
*/
    }
INIT_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (ierr != 0){goto FINISH;}
    MPI_Bcast(&npInObsGroups, 1, MPI_INT, master, MPI_COMM_WORLD);
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
    // load the observations and greens functions then broadcast
    if (myid == master)
    {   
        printf("%s: Loading data...\n", PROGRAM_NAME);
        memset(&data, 0, sizeof(struct parmtData_struct));
        ierr = utils_dataArchive_readAllWaveforms(parms.dataFile, &data);
        data.est = (struct sacData_struct *)
                   calloc((size_t) data.nobs, sizeof(struct sacData_struct));
        ngridSearch = (int64_t) (mtsearch.nm*mtsearch.nb*mtsearch.ng)
                     *(int64_t) (mtsearch.nk*mtsearch.ns*mtsearch.nt)
                     *(int64_t) data.nlocs;
        if (ngridSearch > INT_MAX)
        {
            printf("%s: Error - grid search too big - integer overflow",
                   PROGRAM_NAME);
            ierr = 1;
        }
    }
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (ierr != 0){goto FINISH;}
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
    // Send the waveforms to everyone
    ierr = parmt_broadcast_data(&data, master, MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("%s: Fatal error broadcasting waveforms\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    // Discretize the moment tensor space in (u, v, h) space
    if (myid == master)
    {
        printf("%s: Discretizing MT space...\n", PROGRAM_NAME);
    }
    // Discretize the moment tensor space in (u, v, h) space
    ierr = parmt_discretizeCells64f(
               mtsearch.nb, mtsearch.betaLower,  mtsearch.betaUpper,
               mtsearch.ng, mtsearch.gammaLower, mtsearch.gammaUpper,
               mtsearch.nk, mtsearch.kappaLower, mtsearch.kappaUpper,
               mtsearch.ns, mtsearch.sigmaLower, mtsearch.sigmaUpper,
               mtsearch.nt, mtsearch.thetaLower, mtsearch.thetaUpper,
               mtsearch.nm, mtsearch.m0Lower,    mtsearch.m0Upper,
               mtsearch.luseLog,
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
//if (myid == master){printf("%f %f %f\n", du, dv, dh);}
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
    kappas = array_linspace64f(mtsearch.kappaLower+dk/2, mtsearch.kappaUpper-dk/2,
                               mtsearch.nk, &ierr);
    sigmas = array_linspace64f(mtsearch.sigmaLower+ds/2, mtsearch.sigmaUpper-ds/2,
                               mtsearch.ns, &ierr);
*/
//if (myid != master){for (int i=0; i<mtsearch.nb; i++){printf("%f %f\n", u[i],betas[i]);}}
    // avoid an annoying warning
    for (i=0; i<mtsearch.ns; i++)
    {
        if (fabs(sigmas[i]) < 1.e-7){sigmas[i] = 1.e-7;}
    }
    nmt = mtsearch.ng*mtsearch.nb*mtsearch.nk
         *mtsearch.ns*mtsearch.nt*mtsearch.nm;
    if (myid == master)
    {
        printf("%s: Discretizing moment tensor...\n", PROGRAM_NAME);
    }
    t0 = MPI_Wtime();
    ierr = parmt_discretizeMT64f_MPI(mtComm,
                                     mtsearch.ng, gammas,
                                     mtsearch.nb, betas,
                                     mtsearch.nm, M0s,
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
    // Initialize the  output file
    t0 = MPI_Wtime();
    if (myid == master)
    {
        printf("%s: Initializing output file...\n", PROGRAM_NAME);
        deps = memory_calloc64f(data.nlocs); 
        for (i=0; i<data.nlocs; i++)
        {
            deps[i] = data.sacGxx[i].header.evdp;
        }
        ierr = parmt_io_createObjfnArchive64f(//parms.resultsDir, parms.projnm,
                                              //parms.resultsFileSuffix,
                                              PROGRAM_NAME, parms.parmtArchive, 
                                              data.nobs,
                                              data.nlocs, deps,
                                              mtsearch.nm, M0s,
                                              mtsearch.nb, betas,
                                              mtsearch.ng, gammas,
                                              mtsearch.nk, kappas,
                                              mtsearch.ns, sigmas,
                                              mtsearch.nt, thetas);
        memory_free64f(&deps);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Compute the objective functions
    t0 = MPI_Wtime();
    phi = NULL;
    lags = NULL;
/*
    nlags = 0;
double lagTime;
bool ldefault, lwantLags;
int iobs = 0;
int nlags = 0;
    phi = NULL;
    lags = NULL;
    nlags = 0;
    if (parms.lwantLags)
    {
        //lagTime = data.data[iobs].header.user0;
        //if (lagTime < 0.0){lagTime = parms.defaultMaxLagTime;}
        lagTime = parmt_utils_getLagTime(data.data[iobs], 
                                         parms.defaultMaxLagTime,
                                         &ldefault);
        nlags = (int) (lagTime/data.data[iobs].header.delta + 0.5);
    }
    lags = NULL;
    lwantLags = false;
    if (nlags > 0){lwantLags = true;}
*/
    if (myid == master)
    {
        phi = memory_calloc64f(data.nlocs*mtloc.nmtAll);
    }
    else
    {
        phi = memory_calloc64f(1);
    }
    // Perform the grid search
    MPI_Wtime();
    ierr = parmt_obsSearch64f(MPI_COMM_WORLD,
                              obsComm, locComm,
                              linObsComm, linLocComm,
                              parms, data, mtloc, phi);
    if (ierr != 0)
    {
        printf("%s: Error calling obsSearch64f\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == master)
    {
        MPI_Wtime();
        printf("%s: Objective function computation time: %f\n",
               PROGRAM_NAME, MPI_Wtime() - t0);
    }
    if (ierr != 0)
    {
        printf("%s: Error calling locSearchXC64f\n", PROGRAM_NAME);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    if (myid == master)
    {
        printf("%s: Writing results...\n", PROGRAM_NAME);
        ierr = parmt_io_writeObjectiveFunction64f(
                   //parms.resultsDir, parms.projnm, parms.resultsFileSuffix,
                   parms.parmtArchive, nmt, phi);
        if (ierr != 0)
        {
            printf("%s: Error writing objective function\n", PROGRAM_NAME);
            goto FINISH;
        }
        int jloc, jm, jb, jg, jk, js, jt; 
        int imtopt = array_argmax64f(data.nlocs*mtloc.nmtAll, phi);
        marginal_getOptimum(data.nlocs, mtsearch.nm, mtsearch.nb,
                            mtsearch.ng, mtsearch.nk, mtsearch.ns,
                            mtsearch.nt, phi,
                            &jloc, &jm, &jb, &jg,
                            &jk, &js, &jt);
        printf("%d %f %f %f %f %f %f\n", jloc,
               gammas[jg]*180.0/M_PI, 90.0-betas[jb]*180.0/M_PI, 1.0,
               kappas[jk]*180.0/M_PI, thetas[jt]*180.0/M_PI,
               sigmas[js]*180.0/M_PI);
        double Muse[6], Mned[6], lam[3], U[9];
        printf("phi opt: %d %e\n", imtopt, phi[imtopt]);
        printf("%d %d %d %d %d %d %d\n", jloc, jg, jb, jm, jk, jt, js);
        compearth_tt2cmt(gammas[jg]*180.0/M_PI,
                         90.0-betas[jb]*180.0/M_PI,
                         M0s[jm],
                         kappas[jk]*180.0/M_PI,
                         thetas[jt]*180.0/M_PI,
                         sigmas[js]*180.0/M_PI,
                         Muse, lam, U);
        parmt_discretizeMT64f(1, &gammas[jg],
                              1, &betas[jb],
                              1, &M0s[jm],
                              1, &kappas[jk],
                              1, &thetas[jt],
                              1, &sigmas[js],
                              6, 1, Mned);
        //printf("%e %e %e\n", phi[0], phi[1000], phi[500000]); 
        printf("mtUSE =[%f,%f,%f,%f,%f,%f]\n",
               Muse[0], Muse[1], Muse[2], Muse[3], Muse[4], Muse[5]);
        printf("mtNED =[%f,%f,%f,%f,%f,%f]\n",
               Mned[0], Mned[1], Mned[2], Mned[3], Mned[4], Mned[5]);
    }
FINISH:;
    if (linObsComm){MPI_Comm_free(&obsComm);}
    if (linLocComm){MPI_Comm_free(&locComm);}
    if (linMTComm){MPI_Comm_free(&mtComm);}
    memory_free64f(&M0s);
    memory_free64f(&h);
    memory_free64f(&u);
    memory_free64f(&v);
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&sigmas);
    memory_free64f(&thetas);
    memory_free64f(&phi);
    memory_free32i(&lags);
    parmt_freeData(&data);
    parmt_freeLocalMTs(&mtloc);
    iscl_finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}

static int parseArguments(int argc, char *argv[],
                          const int nprocs,
                          int *npInObsGroups,
                          int *npInLocGroups,
                          int *npInMTGroups,
                          char iniFile[PATH_MAX])
{
    bool linFile;
    int prod;
    linFile = false;
    *npInObsGroups = 1;
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
            {"obsGroupSize", required_argument, 0, 'w'},
            {"mtGroupSize",  required_argument, 0, 'm'},
            {"locGroupSize", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        }; 
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:w:m:l:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            linFile = true; 
        }
        else if (c == 'w' && nprocs > 1)
        {
            *npInObsGroups = atoi(optarg);
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
    prod = (*npInObsGroups)*(*npInLocGroups)*(*npInMTGroups);
    if (prod != nprocs)
    {
        if (prod != 1)
        {
            printf("%s: Invalid process layout - defaulting to (%d,%d,%d)\n", 
                   PROGRAM_NAME, 1, 1, nprocs); 
        }
        *npInObsGroups = 1;
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
    printf("   -w number of processes in waveform (observation) groups\n");
    printf("   -l number of processes in location groups\n");
    printf("   -h displays this message\n");
    printf("   Note if m*w*l must equal the number of processes\n");
    printf("   If they are not set the program will set -m=nprocs\n");
    return;
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

int parmt_freeLocalMTs(struct localMT_struct *mtloc)
{
    memory_free64f(&mtloc->mts);
    //memory_free32i(&mtloc->l2g); // TODO remove
    memory_free32i(&mtloc->offset);
    memset(mtloc, 0, sizeof(struct localMT_struct));
    return 0;
}
