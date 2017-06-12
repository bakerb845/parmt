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
#define DO_P_PICKS true
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
int prepmt_telep_readPickModel(const char *iniFile,
                               const char *section,
                               bool *lsetNewPicks,
                               bool *lusePickFile,
                               char pickFile[PATH_MAX],
                               char ttimesTableDir[PATH_MAX],
                               char ttimesModel[128]);

int main(int argc, char **argv)
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
    bool lsetNewPicks, lusePickFile, lwriteIntermediateFiles;
    const double dmin = PREPMT_MIN_TELESEISMIC_DIST;
    const double dmax = PREPMT_MAX_TELESEISMIC_DIST;
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
    ierr = prepmt_telep_readPickModel(iniFile,
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
        ierr = sacio_readTimeSeriesFile(sacFiles[k], &sacData[k]); 
        if (ierr != 0)
        {
            printf("%s: Failed to read sacFile: %s\n",
                   PROGRAM_NAME, sacFiles[k]);
            return EXIT_FAILURE;
        }
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
        ierr = prepmt_prepData_verifyTeleseismicDistance(&dmin, &dmax,
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
                                                      NULL, "ak135", DO_P_PICKS,
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
    ierr = prepmt_telep_archiveWaveforms(archiveFile, nfiles, sacData);
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
//============================================================================//
int prepmt_telep_archiveWaveforms(const char *archiveFile,
                                  const int nobs,
                                  const struct sacData_struct *data)
{
    const char *fcnm = "prepmt_telep_archiveWaveforms\0";
    char dir[PATH_MAX];
    char **csplit;
    char objName[256];
    int ierr, k, nsplit;
    hid_t attr, dataSpace, fileID, groupID;
    size_t lc;
    ierr = 0;
    // Check inputs
    if (archiveFile == NULL)
    {
        printf("%s: Error archive file is NULL\n", fcnm);
        return -1;
    }
    if (strlen(archiveFile) == 0)
    {
        printf("%s: Error archive file is blank\n", fcnm);
        return -1;
    }
    // Ensure the archive directory exists
    memset(dir, 0, PATH_MAX*sizeof(char));
    csplit = string_rsplit("/", archiveFile, &nsplit);
    if (nsplit > 1 && csplit != NULL)
    {
        // Get the directory
        lc = strlen(archiveFile) - strlen(csplit[nsplit-1]) + 1;
        strncpy(dir, archiveFile, lc);
        if (!os_path_isdir(dir))
        {
            ierr = os_makedirs(dir);
            if (ierr != 0)
            {
                printf("%s: Failed to make output dircectory: %s\n", fcnm, dir);
                return -1;
            }
        }
        for (k=0; k<nsplit; k++)
        {
            if (csplit[k] != NULL){free(csplit[k]);}
        }
        free(csplit);
    }
    if (os_path_isfile(archiveFile))
    {
        printf("%s: Clobbering file %s\n", fcnm, archiveFile);
    } 
    // Create a brand new file to hold the data
    fileID = H5Fcreate(archiveFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Write to /ObservedWaveforms directory
    groupID = H5Gcreate2(fileID, "/ObservedWaveforms", 
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataSpace = H5Screate(H5S_SCALAR);
    attr = H5Acreate(groupID, "NumberOfWaveforms", H5T_NATIVE_INT,
                     dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &nobs);
    H5Aclose(attr);
    H5Sclose(dataSpace);
    // Write each waveform
    for (k=0; k<nobs; k++)
    {
        memset(objName, 0, 256*sizeof(char));
        sprintf(objName, "%s.%s.%s.%s.SAC",
                data[k].header.knetwk,
                data[k].header.kstnm,
                data[k].header.kcmpnm,
                data[k].header.khole);
        ierr = sacioh5_writeTimeSeries2(objName, groupID, data[k]);
        if (ierr != 0)
        {
            printf("%s: Error writing waveform: %d\n", fcnm, k + 1);
            goto ERROR;
        }
    }
/*
 int nfilesRef;
 char **files;
files = sacioh5_getFilesInGroup(groupID, &nfilesRef, &ierr);
for (k=0; k<nfilesRef; k++)
{
 printf("%s\n", files[k]);
}
*/
ERROR:;
    H5Gclose(groupID);
    H5Fclose(fileID);
    return ierr;
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
 *                               file from which to set the picks.
 * @param[out] ttimesTableDir    Directory where the ttimes precomputed
 *                               models reside.  This is not defined if
 *                               lusePickFile is true.
 * @param[out] ttimesModel       Name of ttimes model (e.g., ak135 or iasp91).
 *                               This is not defined if lusePickFile is true.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_telep_readPickModel(const char *iniFile,
                               const char *prefix,
                               bool *lsetNewPicks,
                               bool *lusePickFile,
                               char pickFile[PATH_MAX],
                               char ttimesTableDir[PATH_MAX],
                               char ttimesModel[128])
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
    memset(pickFile, 0, PATH_MAX*sizeof(char));
    memset(ttimesModel, 0, 128*sizeof(char)); 
    memset(ttimesTableDir, 0, PATH_MAX*sizeof(char));

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
        sprintf(vname, "%s:ttimesModel", prefix);
        s = iniparser_getstring(ini, vname, TTIMES_DEFAULT_MODEL);
        strcpy(ttimesModel, s);

        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:ttimesTableDir", prefix);
        s = iniparser_getstring(ini, vname, TTIMES_DEFAULT_TABLE_DIRECTORY);
        if (!os_path_isdir(s))
        {
            printf("%s: ttimes table directory %s doesn't exist\n", fcnm, s);
            ierr = 1;
            goto END;
        }
        strcpy(ttimesTableDir, s);
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
