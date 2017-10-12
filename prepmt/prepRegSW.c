#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <limits.h>
#include <iniparser.h>
#include "prepmt/prepmt.h"
#include "iscl/iscl/iscl.h"
#include "iscl/os/os.h"
#include "iscl/string/string.h"
#include "ttimes.h"
 
#define PROGRAM_NAME "xprepRegSW"
#define DO_RAYLEIGH_PICKS true
static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256]);
int prepmt_telep_intermediateFileOptions(const char *iniFile,
                                        const char *section,
                                        bool *lwriteIntermediateFiles,
                                        char wfDir[PATH_MAX],
                                        char wfSuffix[128],
                                        char archiveFile[PATH_MAX]);
int prepmt_telep_archiveWaveforms(const char *archiveFile,
                                  const int nobs,
                                  const struct sacData_struct *data);
int prepmt_regsw_readPickModel(const char *iniFile,
                               const char *section,
                               bool *lsetNewPicks,
                               bool *lusePickFile,
                               char pickFile[PATH_MAX],
                               double *surfaceWaveVelocity);

int main(int argc, char **argv)
{
    char **sacFiles, **sacpzFiles;
    char iniFile[PATH_MAX], pickFile[PATH_MAX], wfDir[PATH_MAX],
         fname[PATH_MAX], archiveFile[PATH_MAX], section[256], wfSuffix[128];
    struct prepmtEventParms_struct event;
    struct sacData_struct *sacData;
    struct prepmtCommands_struct cmds;
    struct prepmtModifyCommands_struct options;
    struct hpulse96_parms_struct hpulse96Parms;
    double surfaceWaveVel;
    int ierr, k, nfiles;
    bool lsetNewPicks, lusePickFile, lwriteIntermediateFiles;
    char khole[8];
    const double dmin = PREPMT_MIN_REGIONAL_DIST;
    const double dmax = PREPMT_MAX_REGIONAL_DIST;
    const char *hpulseSection = "hpulse96\0";
    //------------------------------------------------------------------------// 
    // Initialize
    memset(&event, 0, sizeof(struct prepmtEventParms_struct));
    iscl_init();
    ierr = parseArguments(argc, argv, iniFile, section);
    if (ierr != 0){return EXIT_FAILURE;}
    // Load the event
    ierr = prepmt_event_initializeFromIniFile(iniFile, &event);
    if (ierr != 0){return EXIT_FAILURE;}
    // Load the pick model
    ierr = prepmt_regsw_readPickModel(iniFile,
                                      section,
                                      &lsetNewPicks, &lusePickFile,
                                      pickFile, &surfaceWaveVel);
    if (ierr != 0){return EXIT_FAILURE;}
    // Figure out some sampling period and windowing info
    memset(&options, 0, sizeof(struct prepmtModifyCommands_struct));
    prepmt_hpulse96_readHpulse96Parameters(iniFile, hpulseSection,
                                           &hpulse96Parms);
    options.iodva = hpulse96Parms.iodva;
    options.ldeconvolution = true; // We are working the data
    ierr = prepmt_prepData_getDefaultDTAndWindowFromIniFile(iniFile, section,
                                                            &options.targetDt,
                                                            &options.cut0,
                                                            &options.cut1);
    if (ierr != 0){return EXIT_FAILURE;}
    // Load the scratch output data info
    ierr = prepmt_telep_intermediateFileOptions(iniFile, section,
                                                &lwriteIntermediateFiles,
                                                wfDir, wfSuffix, archiveFile);
    // Read the data
    printf("%s: Loading data...\n", PROGRAM_NAME);
    ierr = prepmt_prepData_readDataListFromIniFile(iniFile,
                                               section,
                                               &nfiles, &sacFiles, &sacpzFiles);
    if (ierr != 0 || nfiles < 1)
    {
        printf("%s: Failed to find data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    sacData = (struct sacData_struct *)
              calloc((size_t) nfiles, sizeof(struct sacData_struct));
    for (k=0; k<nfiles; k++)
    {
        // Read the data
        ierr = sacio_readTimeSeriesFile(sacFiles[k], &sacData[k]); 
        if (ierr != 0)
        {
            printf("%s: Failed to read sacFile: %s\n",
                   PROGRAM_NAME, sacFiles[k]);
            return EXIT_FAILURE;
        }
        // Fix the location code
        ierr = sacio_getCharacterHeader(SAC_CHAR_KHOLE, sacData[k].header,
                                        khole);
        if (ierr != 0 || strcasecmp(khole, "-12345\0") == 0)
        {
            sacio_setCharacterHeader(SAC_CHAR_KHOLE, "--\0",
                                     &sacData[k].header);
        }
        // Get the metadata
        if (os_path_isfile(sacpzFiles[k]))
        {
            ierr = sacio_readPoleZeroFile(sacpzFiles[k], &sacData[k].pz);
            if (ierr != 0)
            {
                printf("%s: Failed to read pole-zero information\n",
                       PROGRAM_NAME);
            }
        }
        free(sacpzFiles[k]);
        free(sacFiles[k]);
    }
    free(sacpzFiles);
    free(sacFiles);
    // Read the processing commands
    cmds = prepmt_commands_readFromIniFile(iniFile, section,
                                           nfiles, sacData, &ierr);
    if (ierr != 0)
    {
        printf("%s: Failed to get processing commands\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Attach the event to the SAC data
    ierr = prepmt_prepData_setEventInformation(event.latitude,
                                               event.longitude,
                                               event.depth,
                                               event.time,
                                               nfiles, sacData);
    if (ierr != 0)
    {
        printf("%s: Failed to set event information\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Verify the data ranges
    printf("%s: Verifying data ranges...\n", PROGRAM_NAME);
    for (k=0; k<nfiles; k++)
    {
        ierr = prepmt_prepData_verifyRegionalDistance(&dmin, &dmax,
                                                     sacData[k]);
        if (ierr != 0)
        {
            printf("%s: Cannot verify distance - skipping\n", PROGRAM_NAME);
            continue;
        } 
    }
    // Attach picks? 
    if (lsetNewPicks)
    {
        if (!lusePickFile)
        {
            printf("%s: Setting theoretical surface wave picks\n",
                   PROGRAM_NAME);
            ierr = prepmt_prepData_setTheoreticalSurfaceWaveArrivalTime(
                         surfaceWaveVel, DO_RAYLEIGH_PICKS,
                         SAC_FLOAT_A, SAC_CHAR_KA, nfiles, sacData);
            if (ierr != 0)
            {
                printf("%s: Failed to set theoretical surface wave picks\n",
                       PROGRAM_NAME);
                return EXIT_FAILURE;
            }
        }
        else
        {
            printf("%s: Attaching picks from NLL file: %s\n",
                   PROGRAM_NAME, pickFile);
            ierr = prepmt_pickFile_nonLinLoc2sac(pickFile,
                                                 SAC_FLOAT_A, SAC_CHAR_KA,
                                                 nfiles, sacData);
            if (ierr != 0)
            {
                printf("%s: Failed to attached picks to data\n",
                       PROGRAM_NAME);
                return EXIT_FAILURE;
            }
        }
    }
    // Modify the processing commands
    for (k=0; k<cmds.nobs; k++)
    {
        ierr = prepmt_commands_modifyCommandsCharsStruct(options,
                                                         sacData[k],
                                                         &cmds.cmds[k]);
        if (ierr != 0)
        {
            printf("%s: Failed to modify processing commands\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
    }
    // Process the data
    printf("%s: Processing data...\n", PROGRAM_NAME);
    ierr = prepmt_prepData_process(cmds, nfiles, sacData);

    if (lwriteIntermediateFiles)
    {
        printf("%s: Writing intermediate files to directory: %s\n",
               PROGRAM_NAME, wfDir);
        for (k=0; k<nfiles; k++)
        {
            memset(fname, 0, PATH_MAX*sizeof(char));
            if (strlen(wfSuffix) > 0)
            {
                sprintf(fname, "%s/%s.%s.%s.%s.%s.SAC",
                        wfDir,
                        sacData[k].header.knetwk,
                        sacData[k].header.kstnm,
                        sacData[k].header.kcmpnm,
                        sacData[k].header.khole,
                        wfSuffix);
            }
            else
            {
                sprintf(fname, "%s/%s.%s.%s.%s.SAC",
                        wfDir,
                        sacData[k].header.knetwk,
                        sacData[k].header.kstnm,
                        sacData[k].header.kcmpnm,
                        sacData[k].header.khole);
            }
            sacio_writeTimeSeriesFile(fname, sacData[k]);
        }
    }
    // Write the results to an h5 archive
    ierr = prepmt_prepData_archiveWaveforms(archiveFile, nfiles, sacData);
    if (ierr != 0)
    {
        printf("%s: Error writing archive!\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Free data
    prepmt_commands_freePrepmtCommands(&cmds);
    for (k=0; k<nfiles; k++) 
    {    
        sacio_free(&sacData[k]);
    }    
    free(sacData);
    iscl_finalize();
    return 0;
}
//============================================================================//
int prepmt_telep_intermediateFileOptions(const char *iniFile,
                                        const char *section,
                                        bool *lwriteIntermediateFiles,
                                        char wfDir[PATH_MAX],
                                        char wfSuffix[128],
                                        char archiveFile[PATH_MAX])
{
    const char *fcnm = "prepmt_telep_intermediateFileOptions\0";
    const char *s;
    char defaultName[PATH_MAX], vname[256];
    dictionary *ini;
    int ierr;
    size_t lenos;
    //------------------------------------------------------------------------//
    ierr = 0;
    *lwriteIntermediateFiles = false;
    memset(wfDir, 0, PATH_MAX*sizeof(char));
    memset(wfSuffix, 0, 128*sizeof(char));
    memset(archiveFile, 0, PATH_MAX*sizeof(char));
    if (!os_path_isfile(iniFile))
    {
        printf("%s: Error - ini file %s does not exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lwriteIntermediateFiles", section);
    *lwriteIntermediateFiles = iniparser_getboolean(ini, vname, true);
    if (*lwriteIntermediateFiles)
    {
        memset(vname, 0, 256*sizeof(char)); 
        sprintf(vname, "%s:intermediateWaveformDir", section);
        s = iniparser_getstring(ini, vname, "./\0");
        if (!os_path_isdir(s))
        {
            ierr = os_makedirs(s);
            if (ierr != 0)
            {
                printf("%s: Couldn't make output directory %s\n", fcnm, s);
                ierr = 1;
                goto END;
            }
        }
        strcpy(wfDir, s);
        lenos = strlen(wfDir);
        if (lenos > 0)
        {
            if (wfDir[lenos-1] != '/'){wfDir[lenos] = '/';}
        }

        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:outputSuffix", section);
        s = iniparser_getstring(ini, vname, NULL);
        if (s != NULL){strcpy(wfSuffix, s);}
    }
    memset(vname, 0, 256*sizeof(char));
    memset(defaultName, 0, PATH_MAX*sizeof(char));
    sprintf(vname, "%s:archivedWaveforms", section);
    sprintf(defaultName, "%sobservedWaveforms.h5", wfDir);
    s = iniparser_getstring(ini, vname, defaultName);
    strcpy(archiveFile, s);
END:;
    iniparser_freedict(ini); 
    return ierr;
}
//============================================================================//
/*!
 * @brief Defines the way picks will be added to the data.
 *
 * @param[out] lsetNewPicks      If true then the program will update the picks.
 * @param[out] lusePickFile      If true then picks are derived from a pick
 *                               file.  This is false if lsetNewPicks is false.
 * @param[out] pickFile          If lusePickFile is true then this is the pick
 * @param[out] surfaceWaveVel    Default surface-wave velocity for making
 *                               theoretical arrival time estimates.
 *
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_regsw_readPickModel(const char *iniFile,
                               const char *prefix,
                               bool *lsetNewPicks,
                               bool *lusePickFile,
                               char pickFile[PATH_MAX],
                               double *surfaceWaveVel)
{
    const char *fcnm = "prepmt_telep_readPickModel\0";
    const char *s;
    char vname[256];
    int ierr;
    dictionary *ini;
    //------------------------------------------------------------------------//
    ierr = 0;
    *lsetNewPicks = false;
    *lusePickFile = false;
    *surfaceWaveVel = 0.9*3.2;
    memset(pickFile, 0, PATH_MAX*sizeof(char));

    if (!os_path_isfile(iniFile))
    {   
        printf("%s: Error - ini file %s does not exist\n", fcnm, iniFile); 
        return -1;
    }
    ini = iniparser_load(iniFile);
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lsetNewPicks", prefix);
    *lsetNewPicks = iniparser_getboolean(ini, vname, false);
    if (!*lsetNewPicks){goto END;} 

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lusePickFile", prefix);
    *lusePickFile = iniparser_getboolean(ini, vname, false);

    if (*lusePickFile)
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:pickFile", prefix);
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            printf("%s: Pick file %s does not exist\n", fcnm, s);
            ierr = 1;
            goto END;
        }
        strcpy(pickFile, s); 
    }
    else
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:surfaceWaveVel", prefix);
        *surfaceWaveVel = iniparser_getdouble(ini, vname, 3.2*0.9);
        if (*surfaceWaveVel <= 0.0)
        {
            printf("%s: surfaceWaveVel must be positive\n", fcnm);
            ierr = 1;
            goto END;
        }
    }
END:;
    iniparser_freedict(ini);
    return ierr;
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
 *                      prepmt:prepSurfaceWaveData
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
    if (!lsection){strcpy(section, "prepmt:prepSurfaceWaveData\0");} 
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
