#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include "compearth.h"
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "parmt_postProcess.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

#define PROGRAM_NAME "postmt"
#define OUTDIR "postprocess"

static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX]);
static void printUsage(void);


int main(int argc, char *argv[])
{
    struct parmtGeneralParms_struct parms;
    FILE *ofl;
    char fname[PATH_MAX];
    char iniFile[PATH_MAX];
    double U[9], Muse[6], lam[3], *betas, *deps, *depMPDF, *depMagMPDF,
           *gammas, *kappas,
           *sigmas, *thetas, *M0s, *phi, dip, Mw, xnorm;
    int ierr, imtopt, jb, jg, jk, jloc, jm, js, jt, nb, ng, nk, nlocs, nm, nmt, ns, nt;
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
    // TODO make this a config
    if (!os_path_isdir(OUTDIR))
    {
        os_makedirs(OUTDIR);
    }
    // Read the archive
    printf("%s: Reading archive...\n", PROGRAM_NAME); 
    ierr = parmt_io_readObjfnArchive64f(
                   parms.resultsDir, parms.projnm, parms.resultsFileSuffix,
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
    compearth_m02mw(1, 1, &M0s[jm], &Mw);
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
                     90.0 - betas[jb]*180.0/M_PI,
                     M0s[jm],
                     kappas[jk]*180.0/M_PI,
                     thetas[jt]*180.0/M_PI,
                     sigmas[js]*180.0/M_PI,
                     Muse, lam, U);
    printf("mtUSE =[%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n",
           Muse[0], Muse[1], Muse[2], Muse[3], Muse[4], Muse[5]);
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
    for (int im=0; im<nm; im++)
    {
        compearth_m02mw(1, 1, &M0s[im], &Mw);
        printf("Writing magnitude: %f\n", Mw);
        memset(fname, 0, PATH_MAX*sizeof(char));
        sprintf(fname, "%s/%s_dep.%d.txt", OUTDIR, parms.projnm, im+1);
        ofl = fopen(fname, "w");
        for (int iloc=0; iloc<nlocs; iloc++)
        {
            int jloc = im*nlocs + iloc;
            fprintf(ofl, "%e %e\n", deps[iloc], depMagMPDF[jloc]);
        }
        fclose(ofl);
    }
    memset(fname, 0, PATH_MAX*sizeof(char));
    sprintf(fname, "%s/%s_dep.txt", OUTDIR, parms.projnm);
    ofl = fopen(fname, "w");
    for (int iloc=0; iloc<nlocs; iloc++)
    {
        fprintf(ofl, "%e %e\n", deps[iloc], depMPDF[iloc]);
    }
    fclose(ofl);
    
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

