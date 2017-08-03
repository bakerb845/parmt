#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <iniparser.h>
#include "prepmt/prepmt.h"
#include "iscl/iscl/iscl.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

#define PROGRAM_NAME "xgrnsRegSW" 

struct ffGreensTable_struct
{
    double *depths; /*!< Event depths (km) [ndepths] */
    double *dists;  /*!< Great circle distances (degrees) [ndists] */
    int ndepths;    /*!< Number of depths. */
    int ndists;     /*!< Number of distances. */
};

int prepmt_grnsRegSW_readParameters(const char *iniFile, const char *section,
                                    char archiveFile[PATH_MAX],
                                    char hspecFile[PATH_MAX],
                                    char parmtDataFile[PATH_MAX],
                                    char modelName[64],
                                    double *surfaceWaveVel,
                                    bool *lalignXC,
                                    bool *luseEnvelope,
                                    bool *lnormalizeXC,
                                    double *maxXCtimeLag,
                                    int *ndepth, double **depths);

static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256]);
static void printUsage(void);
int prepmt_grnsRegSW_freeFFGreensTable(
    struct ffGreensTable_struct *greensTable);
int prepmt_grnsRegSW_loadFFGreensTable(const char *grnsfl, const char *model,
                                       struct ffGreensTable_struct *greensTable);
int prepmt_grnsRegSW_getDistanceIndex(
    const double distance, const struct ffGreensTable_struct greensTable);
int prepmt_grnsRegSW_loadFFGreensForModelDepthDistance(
    const char *grnsfl, const char *modelName,
    const int idep, const int idist,
    struct sacData_struct *sacGrns);
int prepmt_grnsRegSW_getDepthIndex(
    const double depth, const struct ffGreensTable_struct greensTable);

int main(int argc, char **argv)
{
    char iniFile[PATH_MAX], archiveFile[PATH_MAX],
         hspecFile[PATH_MAX], parmtDataFile[PATH_MAX],
         section[256], modelName[64];
    struct sacData_struct *ffGrns, *grns, *sacData;
    struct prepmtEventParms_struct event;
    struct prepmtCommands_struct cmds;
    struct prepmtModifyCommands_struct options;
    struct ffGreensTable_struct ffGrnsTable;
    double *depths, gcarc, maxXCtimeLag, surfaceWaveVel, ptime, o;
    int ierr, id, indx, k, kndx, iobs, jdep, jdist, ndepth, nobs;
    bool lalignXC, luseEnvelope, lnormalizeXC;
    struct hpulse96_parms_struct hpulse96Parms;
    const char *hpulseSection = "hpulse96\0";
    //------------------------------------------------------------------------//
    iscl_init();
    depths = NULL;
    memset(&options, 0, sizeof(struct prepmtModifyCommands_struct));
    ierr = parseArguments(argc, argv, iniFile, section);
    if (ierr != 0)
    {
        if (ierr ==-2){return EXIT_FAILURE;}
        return EXIT_SUCCESS;
    }
    // Load the event
    ierr = prepmt_event_initializeFromIniFile(iniFile, &event);
    if (ierr != 0)
    {
        printf("%s: Error reading event info\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Read the program basics
    ierr = prepmt_grnsRegSW_readParameters(iniFile, section,
                                           archiveFile, hspecFile,
                                           parmtDataFile, modelName,
                                           &surfaceWaveVel,
                                           &lalignXC, &luseEnvelope,
                                           &lnormalizeXC, &maxXCtimeLag,
                                           &ndepth, &depths);
    if (ierr != 0)
    {
        printf("%s: Failed to read ini file\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Need inversion units
    ierr = prepmt_hpulse96_readHpulse96Parameters(iniFile, hpulseSection,
                                                  &hpulse96Parms);
    if (ierr != 0)
    {
        printf("%s: Failed to read hpulse parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Read how the default commands will be modified
    options.iodva = hpulse96Parms.iodva;
    options.ldeconvolution = false; // We are working with the greens fns 
    ierr = prepmt_prepData_getDefaultDTAndWindowFromIniFile(iniFile, section,
                                                            &options.targetDt,
                                                            &options.cut0,
                                                            &options.cut1);
    if (ierr != 0)
    {
        printf("%s: Failed to read default command info\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Load the data
    printf("%s: Loading data...\n", PROGRAM_NAME);
    sacData = prepmt_prepData_readArchivedWaveforms(archiveFile, &nobs, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error loading data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Read the processing list
    printf("%s: Reading pre-processing commands...\n", PROGRAM_NAME);
    cmds = prepmt_commands_readFromIniFile(iniFile, section, nobs,
                                           sacData, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error loading grns prep commands\n", PROGRAM_NAME);
        return -1;
    }
    // Load the fundamental faults table - then do a lookup
    printf("%s: Loading fundamental faults from model %s in %s...\n",
           PROGRAM_NAME, modelName, hspecFile);
    ierr = prepmt_grnsRegSW_loadFFGreensTable(hspecFile, modelName,
                                              &ffGrnsTable);
    if (ierr != 0)
    {
        printf("%s: Error loading ff table\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ffGrns = (struct sacData_struct *)
             calloc( (size_t) (nobs*ndepth*10), sizeof(struct sacData_struct));
    for (iobs=0; iobs<nobs; iobs++)
    {
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                    sacData[iobs].header, &gcarc); 
        if (ierr != 0)
        {
            printf("%s: Error gcarc not defined on SAC header\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
        jdist = prepmt_grnsRegSW_getDistanceIndex(gcarc, ffGrnsTable);
        for (id=0; id<ndepth; id++)
        {
            indx = prepmt_hudson96_observationDepthTstarToIndex(
                       ndepth, 1, iobs, id, 0);
            kndx = indx*10;                                              
            jdep = prepmt_grnsRegSW_getDepthIndex(depths[id], ffGrnsTable);
            //printf("%d %d %d %d\n", jdist, jdep, indx, kndx);
            ierr =  prepmt_grnsRegSW_loadFFGreensForModelDepthDistance(
                                 hspecFile, modelName, jdep, jdist,
                                 &ffGrns[kndx]);
            if (ierr != 0)
            {
                printf("%s: Error loading ffGreen's functions\n", PROGRAM_NAME);
                return EXIT_FAILURE;
            }
            // Need to scrub the arrivals because this is for surface waves
            for (k=0; k<10; k++) 
            {
                sacio_setFloatHeader(SAC_FLOAT_A,  -12345.0,
                                     &ffGrns[kndx+k].header);
                sacio_setFloatHeader(SAC_FLOAT_T0, -12345.0,
                                     &ffGrns[kndx+k].header);
                sacio_setFloatHeader(SAC_FLOAT_T1, -12345.0,
                                     &ffGrns[kndx+k].header);
                sacio_setCharacterHeader(SAC_CHAR_KA, "-12345",
                                         &ffGrns[kndx+k].header);
                sacio_setCharacterHeader(SAC_CHAR_KT0, "-12345",
                                         &ffGrns[kndx+k].header);
                sacio_setCharacterHeader(SAC_CHAR_KT1, "-12345",
                                         &ffGrns[kndx+k].header);
                ierr = sacio_getFloatHeader(SAC_FLOAT_O,
                                            ffGrns[kndx+k].header, &o);
                if (ierr != 0)
                {
                    printf("%s: o not set\n", PROGRAM_NAME);
                    o = 0.0;
                }
                ptime = o + gcarc*111.195/surfaceWaveVel;
                //printf("%f %f\n", o, ptime);
                sacio_setFloatHeader(SAC_FLOAT_A, ptime,
                                     &ffGrns[kndx+k].header);
                sacio_setCharacterHeader(SAC_CHAR_KA, "R",
                                         &ffGrns[kndx+k].header);
            }
        }
    }
    printf("%s: Contextualizing Green's functions...\n", PROGRAM_NAME);
    grns = (struct sacData_struct *)
           calloc((size_t) (6*ndepth*nobs),
                  sizeof(struct sacData_struct));
    ierr = prepmt_greens_ffGreensToGreens(nobs, sacData,
                                          ndepth, 1, ffGrns, grns);
    if (ierr != 0)
    {
        printf("%s: Error contextualizing Greens functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    for (k=0; k<10*nobs*ndepth; k++)
    {
        sacio_free(&ffGrns[k]);
    }
    free(ffGrns);
    // Process the Green's functions
    printf("%s: Modifying Green's functions processing commands...\n",
           PROGRAM_NAME);
    for (k=0; k<cmds.nobs; k++)
    {
        kndx = k*ndepth*6;
        options.targetDt = sacData[k].header.delta;
        if (fabs(sacData[k].header.delta - grns[kndx].header.delta) > 1.e-5)
        {
            printf("%s: Resample not yet done\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
        ierr = prepmt_commands_modifyCommandsCharsStruct(options,
                                                         sacData[k],
                                                         &cmds.cmds[k]);
        if (ierr != 0)
        {
            printf("%s: Failed to modify processing commands\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
    }
    printf("%s: Processing Green's functions...\n", PROGRAM_NAME);
    ierr = prepmt_greens_processHudson96Greens(nobs, 1, ndepth,
                                               cmds, grns);
    if (ierr != 0)
    {
        printf("%s: Error processing Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Align with cross-correlation
    if (lalignXC)
    {
        printf("%s: Refining waveform alignment with cross-correlation\n",
               PROGRAM_NAME);
        for (iobs=0; iobs<nobs; iobs++)
        {
            for (k=0; k<ndepth; k++)
            {
                kndx = iobs*ndepth*6 + k*6;
                ierr = prepmt_greens_xcAlignGreensToData(sacData[iobs],
                                                         luseEnvelope,
                                                         lnormalizeXC,
                                                         maxXCtimeLag,
                                                         &grns[kndx]);
                if (ierr != 0)
                {
                    printf("%s: Error aligning with XC\n", PROGRAM_NAME);
                    return EXIT_FAILURE;
                }
            }
        }
    }
    // Trim
    printf("%s: Windowing Green's functions to data\n", PROGRAM_NAME);
    ierr = prepmt_greens_cutHudson96FromData(nobs, sacData,
                                             ndepth, 1, grns);
    if (ierr != 0)
    {
        printf("%s: Error trimming Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
/*
sacio_writeTimeSeriesFile("test_gxx.sac", grns[0*ndepth*6+6]); 
sacio_writeTimeSeriesFile("test_gyy.sac", grns[0*ndepth*6+7]);
sacio_writeTimeSeriesFile("test_gzz.sac", grns[0*ndepth*6+8]);
sacio_writeTimeSeriesFile("test_gxy.sac", grns[0*ndepth*6+9]);
sacio_writeTimeSeriesFile("test_gxz.sac", grns[0*ndepth*6+10]);
sacio_writeTimeSeriesFile("test_gyz.sac", grns[0*ndepth*6+11]);
*/
    printf("%s: Writing archive: %s\n", PROGRAM_NAME, parmtDataFile);
    ierr = prepmt_greens_writeArchive(parmtDataFile,
                                      nobs, ndepth,
                                      event.latitude, event.longitude,
                                      depths, sacData, grns);
    if (ierr != 0)
    {
        printf("%s: Failed to write archive\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Clean up
    for (k=0; k<6*ndepth*nobs; k++){sacio_freeData(&grns[k]);}
    free(grns);
    for (k=0; k<nobs; k++){sacio_freeData(&sacData[k]);}
    free(sacData); 
    prepmt_grnsRegSW_freeFFGreensTable(&ffGrnsTable);
    memory_free64f(&depths);
    iscl_finalize();
    return 0;
}
//============================================================================//
int prepmt_grnsRegSW_readParameters(const char *iniFile, const char *section,
                                    char archiveFile[PATH_MAX],
                                    char hspecFile[PATH_MAX],
                                    char parmtDataFile[PATH_MAX],
                                    char modelName[64],
                                    double *surfaceWaveVel,
                                    bool *lalignXC,
                                    bool *luseEnvelope,
                                    bool *lnormalizeXC,
                                    double *maxXCtimeLag,
                                    int *ndepth, double **depths)
{
    const char *fcnm = "prepmt_grnsRegSW_readParameters\0";
    char vname[256];
    const char *s;
    char *dirName;
    double *deps, dmax, dmin;
    int ierr;

    deps = NULL; 
    *luseEnvelope = false;
    *lalignXC = false;
    *lnormalizeXC = false;
    *maxXCtimeLag =-1.0;
    *ndepth = 0;
    memset(archiveFile, 0, PATH_MAX*sizeof(char));
    memset(hspecFile, 0, PATH_MAX*sizeof(char));
    memset(modelName, 0, 64*sizeof(char));
    memset(parmtDataFile, 0, PATH_MAX*sizeof(char));

    dictionary *ini;
    if (!os_path_isfile(iniFile))
    {
        printf("%s: Error ini file %s doesn't exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dataH5File", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (!os_path_isfile(s))
    {
        printf("%s: Data file %s does not exist\n", fcnm, s);
        ierr = 1;
        goto ERROR;
    }
    strcpy(archiveFile, s);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:parmtDataFile", section);
    s = iniparser_getstring(ini, vname, "bodyWave.h5"); 
    strcpy(parmtDataFile, s);

    dirName = os_dirname(parmtDataFile);
    if (!os_path_isdir(dirName))
    {    
        ierr = os_makedirs(dirName);
        if (ierr != 0)
        {
            printf("%s: Failed to make directory: %s\n", fcnm, dirName);
            goto ERROR;
        }
    }    
    memory_free8c(&dirName);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:hspecArchiveFile", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (!os_path_isfile(s))
    {
        printf("%s: hspec archive file %s does not exist\n", fcnm, s);
        ierr = 1;
        goto ERROR;
    }
    strcpy(hspecFile, s);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:modelName", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (s == NULL)
    {
        printf("%s: model name is NULL - weirdness may happen\n", fcnm);
    }
    else
    {
        strcpy(modelName, s);
    }

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:surfaceWaveVel", section);
    *surfaceWaveVel = iniparser_getdouble(ini, vname, 3.2);
    if (*surfaceWaveVel <= 0.0)
    {
        printf("%s: Invalid surface wave velocity: %f\n",
               fcnm, *surfaceWaveVel);
        ierr = 1;
        goto ERROR;
    }

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:ndepths\0");
    *ndepth = iniparser_getint(ini, vname, 0);
    if (*ndepth < 1) 
    {    
        printf("%s: Inadequate number of depths %d\n", fcnm, *ndepth);
        ierr = 1; 
        goto ERROR;
    }    

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:depthMin\0");
    dmin = iniparser_getdouble(ini, vname, -1.0);

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:depthMax\0");
    dmax = iniparser_getdouble(ini, vname, -1.0);
    if (dmin < 0.0 || dmin > dmax)
    {    
        printf("%s: Invalid dmin/dmax %f %f\n", fcnm, dmin, dmax);
        ierr = 1; 
        goto ERROR;
    }    
    deps = array_linspace64f(dmin, dmax, *ndepth, &ierr);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lalignXC", section);
    *lalignXC = iniparser_getboolean(ini, vname, false);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:luseEnvelope", section);
    *luseEnvelope = iniparser_getboolean(ini, vname, false);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lnormalizeXC", section);
    *lnormalizeXC = iniparser_getboolean(ini, vname, false);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:maxXCtimeLag", section);
    *maxXCtimeLag = iniparser_getdouble(ini, vname, -1.0);

ERROR:;
    *depths = deps; 
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Parses the command line arguments for the ini file.
 *
 * @param[in] argc      Number of arguments input to the program.
 * @param[in] argv      Command line arguments.
 *
 * @param[out] iniFile  If the result is 0 then this is the ini file.
 * @param[out] section  Section of ini file to read.  The default is
 *                      prepmt:telePGrns. 
 *
 * @result 0 indicates success. \n
 *        -1 indicates the user inquired about the usage. \n
 *        -2 indicates the command line argument is invalid.
 *
 * @author Ben Baker, ISTI
 *
 */
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256])
{
    bool linFile, lsection;
    linFile = false;
    lsection = false;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    memset(section, 0, 256*sizeof(char));
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"ini_file", required_argument, 0, 'i'},
            {"section", required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:s:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            linFile = true;
        }
        else if (c == 's')
        {
            strcpy(section, (const char *) optarg);
            lsection = true;
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
    else
    {
        if (!os_path_isfile(iniFile))
        {
            printf("%s: Error - ini file: %s does not exist\n",
                   PROGRAM_NAME, iniFile);
            return EXIT_FAILURE;
        }
    }
    if (!lsection){strcpy(section, "prepmt:grnsRegSW\0");}
    return 0;
}
//============================================================================//
static void printUsage(void)
{
    printf("Usage:\n   %s -i ini_file\n\n", PROGRAM_NAME);
    printf("Required arguments:\n");
    printf("    -i ini_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("    -h displays this message\n");
    printf("    -s section of ini file to read\n");
    return;
}
//============================================================================//
/*!
 * @brief Loads the fundamental fault Green's functions for the given
 *        model, depth, and distance.
 *
 * @param[in] grnsfl       Name of HDF5 finite fault Green's functions archive.
 * @param[in] modelName    Name of earth model.
 * @param[in] idep         Depth index.
 * @param[in] idist        Distance index.
 *
 * @param[in,out] sacGrns  On input this is an array of length 10.
 *                         On output the 10 elements of sacGrns have been
 *                         filled with the ZDS, ZSS, ZDD, ZEX, RDS, RSS,
 *                         RDD, REX, TDS, and TSS fundamental fault Green's
 *                         functions.
 *
 * @author Ben Baker
 * 
 */
int prepmt_grnsRegSW_loadFFGreensForModelDepthDistance(
    const char *grnsfl, const char *modelName,
    const int idep, const int idist,
    struct sacData_struct *sacGrns)
{
    const char *fcnm = "prepmt_grnsRegSW_loadFFGreensForModelDepthDistance\0";
    char greensName[512], greensGroup[512];
    hid_t fid, groupID;
    int ierr, k;
    const int ngrns = 10;
    const char *cfaults[10] = {"ZDS\0", "ZSS\0", "ZDD", "ZEX\0",
                               "RDS\0", "RSS\0", "RDD", "REX\0",
                               "TDS\0", "TSS\0"};
    fid = H5Fopen(grnsfl, H5F_ACC_RDONLY, H5P_DEFAULT);
    for (k=0; k<ngrns; k++)
    {
        memset(greensGroup, 0, 512*sizeof(char));
        memset(greensName, 0, 512*sizeof(char));
        sprintf(greensGroup, "/%s/Depth_%d", modelName, idep);
        sprintf(greensName, "Greens_%s_%d", cfaults[k], idist);
        if (!H5Lexists(fid, greensGroup, H5P_DEFAULT))
        {
            printf("%s: Green's functions group %s does not exist\n",
                   fcnm, greensGroup);
            return -1;
        }
        groupID = H5Gopen(fid, greensGroup, H5P_DEFAULT);
        if (!H5Lexists(groupID, greensName, H5P_DEFAULT))
        {
            printf("%s: Green's function %s does not exist\n",
                   fcnm, greensName);
            return -1;
        }
        ierr = sacioh5_readTimeSeries2(greensName, groupID, &sacGrns[k]);
        if (ierr != 0)
        {
            printf("%s: Failed to load greens function: %s\n",
                   fcnm, greensName);
            return -1;
        }
        H5Gclose(groupID);
    }
    H5Fclose(fid);
    return 0;
}
//============================================================================//
/*!
 * @brief Utility function which returns the appropriate distance index in
 *        H5 archive.
 *
 * @param[in] distance     Great circle distance (degrees).
 * @param[in] greensTable  Table containing a descriptor of the H5 archive's
 *                         depth and distance pre-tabulation distribution.
 *
 * @result On successful exit this is the distance index in the Green's 
 *         function archive for loading Greens_???_distanceIndex.
 *
 * @author Ben Baker
 *
 */
int prepmt_grnsRegSW_getDistanceIndex(
    const double distance, const struct ffGreensTable_struct greensTable)
{
    const char *fcnm = "greensTable_getDistanceIndex\0";
    double diff, dmin;
    int i, imin; 
    dmin = fabs(greensTable.dists[0] - distance);
    imin = 0;
    if (greensTable.dists[0] > distance)
    {
        printf("%s: Too close - collocating to nearest greens distance\n",
               fcnm);
        imin = 0;
        return imin;
    }
    if (greensTable.dists[greensTable.ndists-1] < distance)
    {
        // less than a kilometer - i don't really care
        if (fabs(greensTable.dists[greensTable.ndists-1] - distance) > 0.01)
        {
            printf("%s: Too far - collocating to furthest distance\n", fcnm);
        }
        imin = greensTable.ndists - 1;
        return imin;
    }
    // Find the smallest difference
    for (i=1; i<greensTable.ndists; i++)
    {
        diff = fabs(greensTable.dists[i] - distance);
        if (diff < dmin)
        {
            dmin = diff;
            imin = i;
        }
    }
    return imin;
}
//============================================================================//
/*!
 * @brief Utility function which returns the appropriate depth index in H5
 *        archive.
 *
 * @param[in] depth        Source depth (km).
 * @param[in] greensTable  Table containing a descriptor of the H5 archive's
 *                         depth and distance pre-tabulation distribution.
 *
 * @result On successful exit this is the depth index in the Green's 
 *         function archive for loading directory Depth_depthIndex.
 *
 * @author Ben Baker
 *
 */
int prepmt_grnsRegSW_getDepthIndex(
    const double depth, const struct ffGreensTable_struct greensTable)
{
    const char *fcnm = "prepmt_grnsRegSW_getDepthIndex\0";
    double diff, dmin;
    int i, imin;
    dmin = fabs(greensTable.depths[0] - depth);
    imin = 0;  
    if (greensTable.depths[0] > depth)
    {
        printf("%s: Too shallow - collocating to nearest greens depth\n",
               fcnm);
        imin = 0;
        return imin;
    }
    if (greensTable.depths[greensTable.ndepths-1] < depth)
    {   
        printf("%s: Too deep - collocating to deepest depth\n", fcnm);
        imin = greensTable.ndepths - 1;
        return imin;
    }
    // Find the smallest difference
    for (i=1; i<greensTable.ndepths; i++)
    {
        diff = fabs(greensTable.depths[i] - depth);
        if (diff < dmin)
        {
            dmin = diff;
            imin = i;
        }
    }
    return imin;
}
//============================================================================//
/*!
 * @brief Loads the 1D fundamental fault Green's function distance/depth table
 *        for the given earth model.
 *
 * @param[in] grnsfl        hspec96 archive file containing fundamental fault
 *                          Green's functions.
 * @param[in] model         Name of earth model.
 *
 * @param[out] greensTable  On successful output contains the Green's function
 *                          distance/depth table. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */
int prepmt_grnsRegSW_loadFFGreensTable(const char *grnsfl, const char *model,
                                       struct ffGreensTable_struct *greensTable)
{
    const char *fcnm = "prepmt_grnsRegSW_loadFFGreensTable\0";
    char cdepths[512], cdists[512];
    int rank;
    hsize_t dims[1];
    hid_t dataSet, dataSpace, fid, memSpace;
    //------------------------------------------------------------------------//
    memset(greensTable, 0, sizeof(struct ffGreensTable_struct));
    if (!os_path_isfile(grnsfl))
    {
        printf("%s: Greens function file %s doesnt exist\n", fcnm, grnsfl);
        return -1;
    }
    fid = H5Fopen(grnsfl, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (!H5Lexists(fid, model, H5P_DEFAULT))
    {
        printf("%s: Model %s does not exist\n", fcnm, model);
        return -1;
    }
    memset(cdepths, 0, 512*sizeof(char));
    memset(cdists, 0, 512*sizeof(char));
    sprintf(cdepths, "/%s/Depths", model);
    sprintf(cdists, "/%s/Distances", model);
    if (!H5Lexists(fid, cdepths, H5P_DEFAULT))
    {
        printf("%s: Depths does not exist\n", fcnm);
        return -1;
    }
    // Read the distances
    dataSet = H5Dopen(fid, cdepths, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank != 1){printf("this will be bad\n");}
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    greensTable->ndepths = (int) dims[0];
    greensTable->depths = memory_calloc64f((int) dims[0]);
    memSpace = H5Screate_simple(rank, dims, NULL);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
            H5P_DEFAULT, greensTable->depths);
    H5Sclose(dataSpace);
    H5Sclose(memSpace);
    H5Dclose(dataSet);
    // Read the depths
    dataSet = H5Dopen(fid, cdists, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank != 1){printf("this will be bad\n");}
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    greensTable->ndists = (int) dims[0];
    greensTable->dists = memory_calloc64f((int) dims[0]);
    memSpace = H5Screate_simple(rank, dims, NULL);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
            H5P_DEFAULT, greensTable->dists);
    H5Sclose(dataSpace);
    H5Sclose(memSpace);
    H5Dclose(dataSet);
    H5Fclose(fid);
    return 0;
}
//============================================================================//
/*!
 * @brief Releases memory on the fundamental fault Green's function table
 *        structure.
 *
 * @param[out] greensTable   On exit all memory has been freed and variables
 *                           set to 0 or NULL.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */
int prepmt_grnsRegSW_freeFFGreensTable(struct ffGreensTable_struct *greensTable)
{
    memory_free64f(&greensTable->dists);
    memory_free64f(&greensTable->depths);
    memset(greensTable, 0, sizeof(struct ffGreensTable_struct));
    return 0;
}

