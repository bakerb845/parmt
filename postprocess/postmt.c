#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <mpi.h>
#include "compearth.h"
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "parmt_polarity.h"
#include "parmt_postProcess.h"
#include "parmt_mtsearch.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

#define PROGRAM_NAME "postmt"
#define OUTDIR "postprocess"
#define WAVOUT_DIR "obsest"

static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX]);
static void printUsage(void);
int parmt_freeData(struct parmtData_struct *data);

int main(int argc, char *argv[])
{
    struct parmtGeneralParms_struct parms;
    struct parmtData_struct data;
    struct parmtPolarityParms_struct polarityParms;
    struct polarityData_struct polarityData;
    FILE *ofl;
    char fname[PATH_MAX];
    char iniFile[PATH_MAX];
    char programNameIn[256];
    bool lpol;
    double U[9], Muse[6], Mned[6], lam[3], *betas, *deps, *depMPDF, *depMagMPDF,
           *G, *gammas, *kappas,
           *sigmas, *thetas, *M0s, *phi, *var, dip, epoch,
           lagTime, Mw, phiLoc, xnorm, xsum;
    int ierr, iobs, imtopt, jb, jg, jk, jloc, jm, joptLoc, 
        js, jt, k, lag, myid, nb, ng, nk, nlags, nlocs,
        nm, nmt, npmax, npts, nprocs, ns, nt,
        provided;
    bool ldefault;
    // Start MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // Initialize
    depMPDF = NULL; 
    depMagMPDF = NULL;
    betas = NULL;
    deps = NULL;
    gammas = NULL;
    kappas = NULL;
    sigmas = NULL;
    thetas = NULL;
    M0s = NULL;
    phi = NULL;
    memset(&parms, 0, sizeof(struct parmtGeneralParms_struct));
    memset(&data, 0, sizeof(struct parmtData_struct));
    memset(&polarityParms, 0, sizeof(struct parmtPolarityParms_struct));
    memset(&polarityData, 0, sizeof(struct polarityData_struct));
    // Parse the input arguments 
    ierr = parseArguments(argc, argv, iniFile);
    if (ierr != 0)
    {
        if (ierr ==-2){return 0;}
        printf("%s: Error parsing arguments\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Load the ini file - from this we should be able to deduce the archive
    ierr = parmt_utils_readGeneralParms(iniFile, &parms);
    if (ierr != 0)
    {
        printf("%s: Error reading general parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    if (!os_path_isfile(parms.postmtFile))
    {
        printf("%s: Archive file %s doesn't exist\n",
               PROGRAM_NAME, parms.postmtFile);
        return EXIT_FAILURE;
    }
    ierr += parmt_utils_readPolarityParms(iniFile, &polarityParms);

    // TODO make this a config
    if (!os_path_isdir(OUTDIR))
    {
        os_makedirs(OUTDIR);
    }
    if (!os_path_isdir(WAVOUT_DIR))
    {
        os_makedirs(WAVOUT_DIR);
    }
    // Load the data
    printf("%s: Reading data...\n", PROGRAM_NAME);
    ierr = utils_dataArchive_readAllWaveforms(parms.dataFile, &data);
    if (ierr != 0)
    {
        printf("%s: Error reading data\n", PROGRAM_NAME);
        goto ERROR;
    }
    data.est = (struct sacData_struct *)
               calloc((size_t) data.nobs, sizeof(struct sacData_struct));
    // Read the archive
    printf("%s: Reading archive %s...\n", PROGRAM_NAME, parms.postmtFile); 
printf("%s %s %s\n", parms.resultsDir, parms.projnm, parms.resultsFileSuffix);
    ierr = parmt_io_readObjfnArchive64f(
                   //parms.resultsDir, parms.projnm, parms.resultsFileSuffix,
                   parms.postmtFile, //parms.parmtArchive,
                   programNameIn,
                   &nlocs, &deps,
                   &nm, &M0s,
                   &nb, &betas,
                   &ng, &gammas,
                   &nk, &kappas,
                   &ns, &sigmas,
                   &nt, &thetas,
                   &nmt, &phi);
    if (ierr != 0)
    {
        printf("%s: Error loading archive\n", PROGRAM_NAME);
        goto ERROR;
    }
/*
double *phi1; 
    ierr = parmt_io_readObjfnArchive64f(
                   "bw", parms.projnm, "bodyWaves",
                   &nlocs, &deps,
                   &nm, &M0s,
                   &nb, &betas,
                   &ng, &gammas,
                   &nk, &kappas,
                   &ns, &sigmas,
                   &nt, &thetas,
                   &nmt, &phi1);
for (int imt=0; imt<nmt; imt++)
{
  phi[imt] = (12.0*phi1[imt] + 11*phi[imt]);
}
        ierr = parmt_io_createObjfnArchive64f("joint", parms.projnm,
                                              "joint",
                                              25, //data.nobs,
                                              nlocs, deps,
                                              nm, M0s,
                                              nb, betas,
                                              ng, gammas,
                                              nk, kappas,
                                              ns, sigmas,
                                              nt, thetas);
        ierr = parmt_io_writeObjectiveFunction64f(
                   "joint", parms.projnm, "joint",
                   nmt, phi);
return 0;
*/
    // Get the optimum moment tensor for the waveforms
    ierr = marginal_getOptimum(nlocs, nm, nb,
                               ng, nk, ns, nt,
                               phi,
                               &jloc, &jm, &jb, &jg,
                               &jk, &js, &jt);
    if (ierr != 0)
    {
        printf("%s: Failed to get optimum\n", PROGRAM_NAME);
        goto ERROR;
    }
    imtopt = jloc*nm*nb*ng*nk*ns*nt
           +      jm*nb*ng*nk*ns*nt
           +         jb*ng*nk*ns*nt
           +            jg*nk*ns*nt
           +               jk*ns*nt
           +                  js*nt
           +                     jt;
//printf("%d %e\n", jm, M0s[jm]);
//getchar();
    compearth_m02mw(1, CE_KANAMORI_1978, &M0s[jm], &Mw);
    joptLoc = jloc;
//printf("%d %d %d %d %d %d %d\n", jloc, jm, jg, jb, jk, js, jt);
//jloc =4; jm = 5; jg =1; jb=0; jk=11; js=2; jt=4;

    printf("%s: Optimum information:\n", PROGRAM_NAME);
    printf("        Value: %f\n", phi[imtopt]);
    printf("        Depth: %f (km)\n", deps[jloc]);
    printf("        Magnitude: %f\n", Mw);
    printf("        Lune longitude: %f (deg)\n", gammas[jg]*180.0/M_PI);
    printf("        Lune latitude: %f (deg)\n", (M_PI_2 - betas[jb])*180.0/M_PI);
    printf("        Strike: %f (deg)\n", kappas[jk]*180.0/M_PI);
    printf("        Slip: %f (deg)\n", sigmas[js]*180.0/M_PI); 
    printf("        Dip: %f (deg)\n", thetas[jt]*180.0/M_PI);
    // Display the moment tensor in USE format which is useful for obspy and gmt
    compearth_tt2cmt(gammas[jg]*180.0/M_PI,
                     (M_PI/2.0 - betas[jb])*180.0/M_PI,
                     M0s[jm],
                     kappas[jk]*180.0/M_PI,
                     thetas[jt]*180.0/M_PI,
                     sigmas[js]*180.0/M_PI,
                     Muse, lam, U);
    ierr = parmt_discretizeMT64f(1, &gammas[jg],
                                 1, &betas[jb],
                                 1, &M0s[jm],
                                 1, &kappas[jk],
                                 1, &thetas[jt],
                                 1, &sigmas[js],
                                 6, 1, Mned);
    printf("mtUSE =[%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n",
           Muse[0], Muse[1], Muse[2], Muse[3], Muse[4], Muse[5]);
    printf("mtNED =[%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n",
           Mned[0], Mned[1], Mned[2], Mned[3], Mned[4], Mned[5]);
//printf("%e\n", M0s[jm]);
//compearth_CMT2m0(1, 1, Mned, &M0s[jm]);
//printf("%e\n", M0s[jm]);
//getchar();
    // Compute the polarities
    lpol = true;
    ierr = parmt_polarity_computeTTimesGreens(MPI_COMM_WORLD, polarityParms,
                                              data, &polarityData);
    if (ierr != 0)
    {   
        printf("%s: Error computing polarity Green's functions\n",
               PROGRAM_NAME);
        ierr = 1;
        lpol = false;
    }
    if (lpol && polarityData.nPolarity > 0)
    {
        double *Gpol = memory_calloc64f(6*polarityData.nPolarity);
        double *est = memory_calloc64f(polarityData.nPolarity);
        double phiPol;
        int ipol, kt;
        kt = jloc*polarityData.nPolarity;
printf("G\n");
        for (ipol=0; ipol<polarityData.nPolarity; ipol++)
        {
            Gpol[6*ipol+0] = polarityData.Gxx[kt+ipol];
            Gpol[6*ipol+1] = polarityData.Gyy[kt+ipol];
            Gpol[6*ipol+2] = polarityData.Gzz[kt+ipol];
            Gpol[6*ipol+3] = polarityData.Gxy[kt+ipol];
            Gpol[6*ipol+4] = polarityData.Gxz[kt+ipol];
            Gpol[6*ipol+5] = polarityData.Gyz[kt+ipol];
printf("%f %f %f %f %f %f\n", Gpol[6*ipol+0], Gpol[6*ipol+1], Gpol[6*ipol+2],
                     Gpol[6*ipol+3], Gpol[6*ipol+4], Gpol[6*ipol+5]);
        }
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    polarityData.nPolarity, 6, 1.0, Gpol, 6,
                    Mned, 1, 0.0, est, 1);
        for (ipol=0; ipol<polarityData.nPolarity; ipol++)
        {
            est[ipol] = copysign(1.0, est[ipol]);
            printf("Polarity: obs %f; est %f\n",
                   polarityData.polarity[ipol], est[ipol]);
        }
        
    }
    // Write the output file

    struct globalMapOpts_struct globalMap;
    memset(&globalMap, 0, sizeof(struct globalMapOpts_struct));
    strcpy(globalMap.outputScript, "gmtScripts/globalMap.sh");
    strcpy(globalMap.psFile, "globalMap.ps");
    array_copy64f_work(6, Muse, globalMap.mts);
    globalMap.basis = CE_USE;
    globalMap.evla = data.data[0].header.evla;
    globalMap.evlo = data.data[0].header.evlo;
    globalMap.evdp = deps[jloc]; 
    globalMap.lwantMT = true;
    globalMap.lwantPolarity = true;
printf("%f %f\n", globalMap.mts[0], globalMap.mts[1]);
    ierr = postmt_gmtHelper_writeGlobalMap(globalMap, data.nobs, data.data);
/*
printf("overriding ned and joptloc\n");
joptLoc=5;
Mned[0] = 8.526325e+15;
Mned[1] =-8.707938e+16;
Mned[2] = 2.226786e+17;
Mned[3] =-5.768950e+17;
Mned[4] =-3.469420e+17;
Mned[5] =-2.756277e+17;
*/
    // Compute the scaling factor
    xnorm = marginal_computeNormalization(nlocs, deps,
                                          nm, M0s,
                                          nb, betas,
                                          ng, gammas,
                                          nk, kappas,
                                          ns, sigmas,
                                          nt, thetas,
                                          phi, &ierr);
    if (xnorm <= 0.0 || ierr != 0)
    {
        printf("%s: Failed to compute normalization factor - setting to 1\n",
               PROGRAM_NAME);
        ierr = 0;
        xnorm = 1.0;
    }
    printf("%s: Scaling factor: %e\n", PROGRAM_NAME, xnorm);
    cblas_dscal(nmt, 1.0/xnorm, phi, 1);
    printf("%s: Computing pure histograms\n", PROGRAM_NAME);

    double *locHist, *magHist, *betaHist, *gammaHist, *kappaHist, *sigmaHist, *thetaHist;
    ierr = postmt_gmtHelper_makeRegularHistograms(nlocs, nm,
                                                  nb, ng, nk, ns, nt, 
                                                  nmt, phi,
                                                  &locHist, &magHist, &betaHist,
                                                  &gammaHist, &kappaHist, &sigmaHist,
                                                  &thetaHist);
    int i;
/*
    printf("Depths\n");
    for (i=0; i<nlocs; i++)
    {
        printf("%f %f\n", deps[i], locHist[i]);
    }
    printf("Magnitudes\n");
    for (i=0; i<nm; i++)
    {
        compearth_m02mw(1, CE_KANAMORI_1978, &M0s[i], &Mw);
        printf("%f %f\n", Mw, magHist[i]);
    }
    printf("Gammas\n");
    for (i=0; i<ng; i++)
    {
        printf("%f %f\n", gammas[i]*180.0/M_PI, gammaHist[i]);
    }
*/
    printf("Betas\n");
    for (i=0; i<nb; i++)
    {
        printf("%f %f\n", 90.0 - betas[i]*180.0/M_PI, betaHist[i]);
    }
/*
    printf("Kappas\n");
    for (i=0; i<nk; i++)
    {
        printf("%f %f\n", kappas[i]*180.0/M_PI, kappaHist[i]);
    }
    printf("Sigmas\n");
    for (i=0; i<ns; i++)
    {
        printf("%f %f\n", sigmas[i]*180.0/M_PI, sigmaHist[i]);
    }
    printf("Thetas\n");
    for (i=0; i<nt; i++)
    {
        printf("%f %f\n", thetas[i]*180.0/M_PI, thetaHist[i]);
    }
*/
const char *scriptFile = "gmtScripts/boxes.sh";
const char *psFile = "boxes.ps";
ierr = postmt_gmtHelper_writeBetaBoxes(false, false, scriptFile, psFile,
                                       nb, jb, betas, betaHist);
ierr = postmt_gmtHelper_writeGammaBoxes(true, false, scriptFile, psFile,
                                        ng, jg, gammas, gammaHist);
ierr = postmt_gmtHelper_writeKappaBoxes(true, false, scriptFile, psFile,
                                        nk, jk, kappas, kappaHist);
ierr = postmt_gmtHelper_writeSigmaBoxes(true, false, scriptFile, psFile,
                                        ns, js, sigmas, sigmaHist);
ierr = postmt_gmtHelper_writeThetaBoxes(true, false, scriptFile, psFile,
                                        nt, jt, thetas, thetaHist);
ierr = postmt_gmtHelper_writeMagnitudeBoxes(true, false, scriptFile, psFile,
                                            nm, jm, M0s, magHist);
ierr = postmt_gmtHelper_writeDepthBoxes(true, true, scriptFile, psFile,
                                        nlocs, joptLoc, deps, locHist);

    // Start the post-processing - compute the depgh magnitude mPDFs
    printf("%s: Computing the depth MPDF's\n", PROGRAM_NAME);
    depMagMPDF = memory_calloc64f(nm*nlocs);
    depMPDF = memory_calloc64f(nlocs);
    ierr = marginal_computeDepthMPDF(nlocs,
                                     nm, M0s,
                                     nb, betas,
                                     ng, gammas,
                                     nk, kappas,
                                     ns, sigmas,
                                     nt, thetas,
                                     phi, depMagMPDF, depMPDF);
    if (ierr != 0)
    {
        printf("%s: Error computing depthMPDF\n", PROGRAM_NAME);
        goto ERROR;
    }
    // Print the station list for my edification
    for (k=0; k<data.nobs; k++)
    {
        printf("%8.3f %8.3f %6s\n", data.data[k].header.stlo,
                                    data.data[k].header.stla,
                                    data.data[k].header.kstnm);
    }
int im, iloc;
    for (im=0; im<nm; im++)
    {
        compearth_m02mw(1, CE_KANAMORI_1978, &M0s[im], &Mw);
        printf("Writing magnitude: %f\n", Mw);
        memset(fname, 0, PATH_MAX*sizeof(char));
        sprintf(fname, "%s/%s_dep.%d.txt", OUTDIR, parms.projnm, im+1);
        ofl = fopen(fname, "w");
        for (iloc=0; iloc<nlocs; iloc++)
        {
            jloc = im*nlocs + iloc;
            fprintf(ofl, "%e %e\n", deps[iloc], depMagMPDF[jloc]);
        }
        fclose(ofl);
    }
    memset(fname, 0, PATH_MAX*sizeof(char));
    sprintf(fname, "%s/%s_dep.txt", OUTDIR, parms.projnm);
    ofl = fopen(fname, "w");
    for (iloc=0; iloc<nlocs; iloc++)
    {
        fprintf(ofl, "%e %e\n", deps[iloc], depMPDF[iloc]);
    }
    fclose(ofl);
    // Compute the optimal synthetics
    xsum = 0.0;
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        // Extract the depth waveform
        npts = data.data[iobs].npts;
        npmax = npts;
        G = memory_calloc64f(6*npmax);
        var = memory_calloc64f(npmax);
        ierr = parmt_utils_setDataOnG(iobs, joptLoc, npmax, data, G);
        if (ierr != 0)
        {
            printf("%s: Failed to set data on G\n", PROGRAM_NAME);
            goto ERROR;
        }
        // Compute the lag time
        nlags = 0;
        lag = 0;
        if (parms.lwantLags)
        {
            //lagTime = data.data[iobs].header.user0;
            //if (lagTime < 0.0){lagTime = parms.defaultMaxLagTime;}
            lagTime = parmt_utils_getLagTime(data.data[iobs], 
                                             parms.defaultMaxLagTime,
                                             &ldefault);
            nlags = (int) (lagTime/data.data[iobs].header.delta + 0.5);
        }
        ierr = parmt_mtSearchL164f(6, 1,
                                   npts, 1,
                                   nlags, parms.lwantLags,
                                   &G[0*npmax], &G[1*npmax], &G[2*npmax],
                                   &G[3*npmax], &G[4*npmax], &G[5*npmax],
                                   NULL, Mned, data.data[iobs].data,
                                   &phiLoc, var, &lag);
        if (ierr != 0)
        {
            printf("%s: Error computing lags for waveform %d\n",
                   PROGRAM_NAME, iobs);
            goto ERROR;
        }
        // Compute the synthetic
        k = iobs*data.nlocs + joptLoc;
        parmt_utils_sacGrnsToEst(data.data[iobs],
                                     data.sacGxx[k], data.sacGyy[k],
                                     data.sacGzz[k], data.sacGxy[k],
                                     data.sacGxz[k], data.sacGyz[k],
                                     Mned,
                                     &data.est[iobs]);
        // Fix the timing
        ierr = sacio_getEpochalStartTime(data.est[iobs].header, &epoch);
        epoch = epoch + (double) lag*data.data[iobs].header.delta;
        sacio_setEpochalStartTime(epoch, &data.est[iobs].header);
        // Write the file
        memset(fname, 0, PATH_MAX*sizeof(char));
        sprintf(fname, "%s/%s.%s.%s.%s.SAC",
                WAVOUT_DIR,
                data.data[iobs].header.knetwk,
                data.data[iobs].header.kstnm,
                data.data[iobs].header.kcmpnm,
                data.data[iobs].header.khole);
        sacio_writeTimeSeriesFile(fname, data.data[iobs]);
        memset(fname, 0, PATH_MAX*sizeof(char));
        sprintf(fname, "%s/%s.%s.%s.%s.EST.SAC",
                WAVOUT_DIR,
                data.data[iobs].header.knetwk,
                data.data[iobs].header.kstnm,
                data.data[iobs].header.kcmpnm,
                data.data[iobs].header.khole);
        sacio_writeTimeSeriesFile(fname, data.est[iobs]);
        printf("Waveform %6s has lag %+4d=%+8.3e (s) and a fit of %f\n",
               data.data[iobs].header.kstnm,
               lag, lag*data.data[iobs].header.delta, phiLoc);
        memory_free64f(&G); 
        memory_free64f(&var);
        xsum = xsum + phiLoc;
    } 
    printf("Average %f\n", xsum/(double) data.nobs );
/*
double *db = memory_calloc64f(nb);
double *dg = memory_calloc64f(ng);
double *dt = memory_calloc64f(nt);
postprocess_computeBetaCellSpacing(nb, betas, db); 
 printf("\n");
postprocess_computeGammaCellSpacing(ng, gammas, dg);
printf("\n");
postprocess_computeThetaCellSpacing(nt, thetas, dt);
*/
    // Free memory
ERROR:;
    parmt_freeData(&data);
    memory_free64f(&depMagMPDF);
    memory_free64f(&depMPDF);
    memory_free64f(&deps);
    memory_free64f(&M0s); 
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&sigmas);
    memory_free64f(&thetas);
    memory_free64f(&phi);
    iscl_finalize();
    MPI_Finalize();
    return ierr;
}
//============================================================================//
/*!
 * @brief Parses the input arguments to get the ini file anme
 */
static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX])
{
    bool linFile;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"ini_file", required_argument, 0, 'i'},
            {0, 0, 0, 0}
        };  
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {   
            strcpy(iniFile, (const char *) optarg);
            linFile = true; 
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
    return 0;
}

static void printUsage(void)
{
    printf("Usage:\n   postmt -i input_file\n\n");
    printf("Required arguments:\n");
    printf("   -i input_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("   -h displays this message\n");
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
