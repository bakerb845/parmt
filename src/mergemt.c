#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <iniparser.h>
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#define PROGRAM_NAME "mergemt"

static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX]);

struct mergemtParms_struct
{
    char **files;
    double *wts;
    char archiveFile[PATH_MAX];
    int nfiles;
};

int mergemt_readIni(const char *iniFile, struct mergemtParms_struct *parms);


int main(int argc, char *argv[])
{
    struct mergemtParms_struct parms;
    char iniFile[PATH_MAX], programNameIn[256]; 
    double *beta0, *beta1, *betaLoc, *dep0, *dep1, *depLoc,
           *gamma0, *gamma1, *gammaLoc, *kappa0, *kappa1, *kappaLoc,
           *M0, *M1, *M0loc, *phi0, *phi1,
           *phiLoc, *sigma0, *sigma1, *sigmaLoc, *theta0, *theta1,
           *thetaLoc, xsum;
    int i, ierr, il, im, indx, jndx, nb, nb0, nb1, ncopy, ng, ng0, ng1,
        nk, nk0, nk1, nlocs, nlocs0, nlocs1, nm, nm0, nm1, nmt, nmt0, nmt1,
        ns, ns0, ns1, nt, nt0, nt1;
    const int jm = 0;
    ierr = parseArguments(argc, argv, iniFile);
    if (ierr != 0)
    {
        if (ierr ==-2){return EXIT_SUCCESS;}
        printf("%s: Error parsing arguments\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    } 
    if (!os_path_isfile(iniFile))
    {
        printf("%s: Error ini file doesn't exist\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = mergemt_readIni(iniFile, &parms);
    if (ierr != 0)
    {
        printf("%s: Error reading ini file\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    xsum = array_sum64f(parms.nfiles, parms.wts, &ierr);
    if (fabs(xsum) < 1.e-14)
    {
        printf("%s: Division by zero - setting to 1\n", PROGRAM_NAME);
        xsum = 1.0;
    }
    cblas_dscal(parms.nfiles, 1.0/xsum, parms.wts, 1);
    // Load the files
    nlocs0 =-1;
    nlocs1 =-1;
    nm0 =-1;
    nm1 =-1;
    nb0 =-1;
    nb1 =-1;
    ng0 =-1;
    ng1 =-1;
    nk0 =-1;
    nk1 =-1;
    ns0 =-1;
    ns1 =-1;
    nt0 =-1;
    nt1 =-1;
    nmt0 =-1;
    nmt1 =-1;
    dep0 = NULL;
    dep1 = NULL;
    M0 = NULL;
    M1 = NULL;
    beta0 = NULL;
    beta1 = NULL;
    gamma0 = NULL;
    gamma1 = NULL;
    kappa0 = NULL;
    kappa1 = NULL;
    sigma0 = NULL;
    sigma1 = NULL;
    theta0 = NULL;
    theta1 = NULL;
    phi0 = NULL;
    phi1 = NULL;
    for (i=0; i<parms.nfiles; i++)
    {
        if (!os_path_isfile(parms.files[i])){continue;}
        printf("%s: Loading archive: %s\n", PROGRAM_NAME, parms.files[i]); 
        ierr = parmt_io_readObjfnArchive64f(parms.files[i],
                                            programNameIn,
                                            &nlocs, &depLoc,
                                            &nm, &M0loc,
                                            &nb, &betaLoc,
                                            &ng, &gammaLoc,
                                            &nk, &kappaLoc,
                                            &ns, &sigmaLoc,
                                            &nt, &thetaLoc,
                                            &nmt, &phiLoc);
        if (ierr != 0)
        {
            printf("%s: Error loading %s; skipping...\n",
                   PROGRAM_NAME, parms.files[i]);
            continue;
        }
        if (strcasecmp("parmt", programNameIn) == 0)
        {
            if (nm0 ==-1)
            {   
                nlocs0 = nlocs;
                nm0 = nm; 
                nb0 = nb; 
                ng0 = ng; 
                nk0 = nk; 
                ns0 = ns; 
                nt0 = nt; 
                nmt0 = nmt;
                phi0 = array_zeros64f(nmt, &ierr);
                dep0 = array_copy64f(nlocs, depLoc, &ierr);
                M0 = array_copy64f(nm, M0loc, &ierr);
                beta0 = array_copy64f(nb, betaLoc, &ierr);
                gamma0 = array_copy64f(ng, gammaLoc, &ierr);
                kappa0 = array_copy64f(nk, kappaLoc, &ierr);
                sigma0 = array_copy64f(ns, sigmaLoc, &ierr);
                theta0 = array_copy64f(nt, thetaLoc, &ierr);
            }
            if (nm0 != nm || nb0 != nb || ng0 != ng || nk0 != nk ||
                ns0 != ns || nt0 != nt || nmt0 != nmt || nlocs0 != nlocs)
            {
                printf("%s: Size mismatch in parmt\n", PROGRAM_NAME);
                continue;
            }
            cblas_daxpy(nmt, parms.wts[i], phiLoc, 1, phi0, 1); 
        }
        else if (strcasecmp("polarmt", programNameIn) == 0)
        {
            if (nm1 ==-1)
            {
                nlocs1 = nlocs;
                nm1 = nm;
                nb1 = nb;
                ng1 = ng;
                nk1 = nk;
                ns1 = ns;
                nt1 = nt;
                nmt1 = nmt;
                phi1 = array_zeros64f(nmt, &ierr);
                dep1 = array_copy64f(nlocs, depLoc, &ierr);
                M1 = array_copy64f(nm, M0loc, &ierr);
                beta1 = array_copy64f(nb, betaLoc, &ierr);
                gamma1 = array_copy64f(ng, gammaLoc, &ierr);
                kappa1 = array_copy64f(nk, kappaLoc, &ierr);
                sigma1 = array_copy64f(ns, sigmaLoc, &ierr);
                theta1 = array_copy64f(nt, thetaLoc, &ierr);
            }
printf("%f\n",parms.wts[i]);
            cblas_daxpy(nmt, parms.wts[i], phiLoc, 1, phi1, 1);
        }
        else
        {
            printf("%s: Unkown origin: %s\n", PROGRAM_NAME, programNameIn); 
            goto NEXT; 
        }

NEXT:;
        memory_free64f(&M0loc);
        memory_free64f(&betaLoc);
        memory_free64f(&gammaLoc);
        memory_free64f(&kappaLoc);
        memory_free64f(&sigmaLoc);
        memory_free64f(&thetaLoc);
        memory_free64f(&phiLoc);
    }
    if (phi0 != NULL && phi1 != NULL)
    {
        if (nb0 != nb1 || ng0 != ng1 || nk0 != nk1 ||
            ns0 != ns1 || nt0 != nt1 || nlocs0 != nlocs)
        {
            printf("%s: Inconsistent grid search sizes\n", PROGRAM_NAME);
            return -1;
        }
        printf("%f\n", array_max64f(nmt0, phi0));
        printf("%f\n", array_max64f(nmt1, phi1)); 
        printf("%s: Stacking polarities into waveforms...\n", PROGRAM_NAME);
/*
        for (i=0; i<nm0; i++)
        {
            //printf("%d %d\n", nmt, nmt0);
            cblas_daxpy(nmt1, 1.0, phi1, 1, &phi0[nmt1*i], 1);
        }
*/
        for (il=0; il<nlocs0; il++)
        {
            for (im=0; im<nm0; im++)
            {
                indx = il*nm0*nb0*ng0*nk0*ns0*nt0
                     + im*nb0*ng0*nk0*ns0*nt0;
                jndx = il*nm1*nb1*ng1*nk1*ns1*nt1
                     + jm*nb1*ng1*nk1*ns1*nt1;
                ncopy = nb0*ng0*nk0*ns0*nt0;
                cblas_daxpy(ncopy,  1.0, &phi1[jndx], 1, &phi0[indx], 1); 
            }
        }
printf("%f\n", array_max64f(nmt0, phi0));
        printf("%s: Writing joint archive\n", PROGRAM_NAME);
        ierr = parmt_io_createObjfnArchive64f(PROGRAM_NAME, parms.archiveFile,
                                              1, //data.nobs, TODO - read and fix
                                              nlocs0, dep0,
                                              nm0, M0,
                                              nb0, beta0,
                                              ng0, gamma0,
                                              nk0, kappa0,
                                              ns0, sigma0,
                                              nt0, theta0);
        ierr = parmt_io_writeObjectiveFunction64f(parms.archiveFile,
                                                  nmt0, phi0);
    }
    else
    {
        if (phi0 != NULL)
        {
            printf("%s: Writing joint waveform archive\n", PROGRAM_NAME);
            ierr = parmt_io_createObjfnArchive64f(PROGRAM_NAME, parms.archiveFile,
                                              1, //data.nobs, TODO - read and fix
                                              nlocs0, dep0,
                                              nm0, M0, 
                                              nb0, beta0,
                                              ng0, gamma0,
                                              nk0, kappa0,
                                              ns0, sigma0,
                                              nt0, theta0);
            ierr = parmt_io_writeObjectiveFunction64f(parms.archiveFile,
                                                      nmt0, phi0);
        }
        else if (phi1 != NULL)
        {
            printf("%s: Writing joint waveform archive\n", PROGRAM_NAME);
            ierr = parmt_io_createObjfnArchive64f(PROGRAM_NAME, parms.archiveFile,
                                              1, //data.nobs, TODO - read and fix
                                              nlocs1, dep1,
                                              nm1, M1, 
                                              nb1, beta1,
                                              ng1, gamma1,
                                              nk1, kappa1,
                                              ns1, sigma1,
                                              nt1, theta1);
        }
    }
    memory_free64f(&dep0);
    memory_free64f(&dep1);
    memory_free64f(&M0);
    memory_free64f(&M1);
    memory_free64f(&beta0);
    memory_free64f(&beta1);
    memory_free64f(&gamma0);
    memory_free64f(&gamma1);
    memory_free64f(&kappa0);
    memory_free64f(&kappa1);
    memory_free64f(&sigma0);
    memory_free64f(&sigma1);
    memory_free64f(&theta0);
    memory_free64f(&theta1);
    memory_free64f(&phi0);
    memory_free64f(&phi1);
    return EXIT_SUCCESS;
}
//============================================================================//
int mergemt_readIni(const char *iniFile,
                    struct mergemtParms_struct *parms) 
{
    const char *fcnm = "mergemt_readIni\0";
    const char *s;
    char vname[256];
    dictionary *ini;
    int i, ierr;
    memset(parms, 0, sizeof(struct mergemtParms_struct));
    if (!os_path_isfile(iniFile))
    {
        printf("%s: Error ini file doesn't exist\n", fcnm); 
        return -1;
    }
    ini = iniparser_load(iniFile);

    parms->nfiles = iniparser_getint(ini, "general:nmerge\0", 0);
    if (parms->nfiles < 1)
    {
        printf("%s: No files to merge\n", fcnm);
        return -1;
    }

    parms->files = (char **) calloc((size_t) parms->nfiles, sizeof(char *));
    parms->wts = array_set64f(parms->nfiles, 1.0, &ierr);
    for (i=0; i<parms->nfiles; i++)
    {
        parms->files[i] = (char *) calloc(PATH_MAX, sizeof(char));
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "general:mergeFile_%d", i+1);
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            printf("%s: File %s doesn't exist\n", fcnm, s);
            continue;
        }
        strcpy(parms->files[i], s); 

        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "general:wtFile_%d", i+1);
        parms->wts[i] = iniparser_getdouble(ini, vname, 1.0);
        if (parms->wts[i] < 0.0)
        {
            printf("%s: Weight can't be negative: %f\n", fcnm, parms->wts[i]);
            return -1;
        }
    } 

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "general:mergeArchive");
    s = iniparser_getstring(ini, vname, "mergemt.h5");
    strcpy(parms->archiveFile, s);

    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX])
{
    bool linFile;
    int prod;
    linFile = false;
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
        c = getopt_long(argc, argv, "?hi:m:l:",
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
//============================================================================//
static void printUsage(void)
{
    printf("Usage:\n   mergemt -i input_file\n\n");
    printf("Required arguments:\n");
    printf("   -i input_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("   -h displays this message\n");
    return;
}
