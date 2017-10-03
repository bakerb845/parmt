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
 
#define PROGRAM_NAME "xprepTeleP"
static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256]);

int main(int argc, char **argv)
{
    char iniFile[PATH_MAX], section[256];
    const char *hpulseSection = "hpulse96\0";
    const bool lisP = true;
    int ierr;
    const double dmin = PREPMT_P_MIN_TELESEISMIC_DIST;
    const double dmax = PREPMT_P_MAX_TELESEISMIC_DIST;
    iscl_init();
    ierr = parseArguments(argc, argv, iniFile, section);
    if (ierr != 0){return EXIT_FAILURE;}
    ierr = prepmt_prepData_prepTeleseismicBodyWaves(iniFile, section,
                                                    hpulseSection,lisP,
                                                    dmin, dmax);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Error preparing data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    iscl_finalize();
    return EXIT_SUCCESS;
}
/*
int main_original(int argc, char **argv)
{
    char **sacFiles, **sacpzFiles;
    char iniFile[PATH_MAX], pickFile[PATH_MAX], wfDir[PATH_MAX],
         ttimesTableDir[PATH_MAX], fname[PATH_MAX], archiveFile[PATH_MAX],
         section[256], ttimesModel[128], wfSuffix[128];
    struct prepmtEventParms_struct event;
    struct sacData_struct *sacData;
    struct prepmtCommands_struct cmds;
    struct prepmtModifyCommands_struct options;
    struct hpulse96_parms_struct hpulse96Parms;
    int ierr, k, nfiles;
    bool lsetNewPicks, lusePickFile, lwrtIntFiles;
    char khole[8];
    const double dmin = PREPMT_P_MIN_TELESEISMIC_DIST;
    const double dmax = PREPMT_P_MAX_TELESEISMIC_DIST;
    const char *hpulseSection = "hpulse96\0";
    const bool lisP = true;
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
    ierr = prepmt_prepData_readPickModel(iniFile,
                                         section,
                                         &lsetNewPicks, &lusePickFile,
                                         pickFile, ttimesTableDir, ttimesModel);
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
    ierr = prepmt_prepData_intermediateFileOptions(iniFile, section,
                                                   &lwrtIntFiles,
                                                   wfDir, wfSuffix,
                                                   archiveFile);
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
        // Get the pole-zero file
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
        ierr = prepmt_prepData_verifyTeleseismicDistance(lisP, &dmin, &dmax,
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
            printf("%s: Setting theoretical P picks\n", PROGRAM_NAME);
            ierr = prepmt_prepData_setPrimaryPorSPickFromTheoreticalTime(
                                                      NULL, "ak135", true,
                                                      SAC_FLOAT_A, SAC_CHAR_KA,
                                                      nfiles, sacData);
            if (ierr != 0)
            {
                printf("%s: Failed to set theoretical primary picks\n",
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

    if (lwrtIntFiles)
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
    free(sacData);
    iscl_finalize();
    return 0;
}
*/
//============================================================================//
/*!
 * @brief Parses the command line arguments for the ini file.
 *
 * @param[in] argc      Number of arguments input to the program.
 * @param[in] argv      Command line arguments.
 *
 * @param[out] iniFile  If the result is 0 then this is the ini file.
 * @param[out] section  Section of ini file to read.  The default is
 *                      prepmt:prepTeleseismicData.
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
    if (!lsection){strcpy(section, "prepmt:prepTeleseismicData\0");} 
    return 0;
}
//============================================================================//
/*!
 * @brief Utility routine to demonstrate how to call the program.
 */
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
