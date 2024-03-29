#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iniparser.h>
#include "sacio.h"
#include "ttimes.h"
#include "prepmt/prepmt.h"
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/geodetic/geodetic.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/string/string.h"
#include "iscl/time/time.h"

#define DO_P_PICKS true
#define DO_S_PICKS false

/*!
 * @brief Preprocesses the teleseismic P or S body waves.
 *
 * @param[in] iniFile        Name of initialization file.
 * @param[in] section        Section of initialization file to read containing
 *                           the data pre-processing parameters.
 * @param[in] hpulseSection  Section of initialization file with hpulse
 *                           parameters. 
 * @param[in] lisP           If true then this is for P waves. \n
 *                           Otherwise, this is for S waves.
 * @param[in] dmin           Minimum teleseismic distance (degrees).
 * @param[in] dmax           Maximum teleseismic distance (degrees).
 *
 * @result EXIT_SUCCESS indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_prepTeleseismicBodyWaves(const char *iniFile,
                                             const char *section,
                                             const char *hpulseSection,
                                             const bool lisP,
                                             const double dmin,
                                             const double dmax)
{
    char **sacFiles, **sacpzFiles;
    char pickFile[PATH_MAX], wfDir[PATH_MAX],
         ttimesTableDir[PATH_MAX], fname[PATH_MAX], archiveFile[PATH_MAX],
         ttimesModel[128], wfSuffix[128];
    struct prepmtEventParms_struct event;
    struct sacData_struct *sacData;
    struct prepmtCommands_struct cmds;
    struct prepmtModifyCommands_struct options;
    struct hpulse96_parms_struct hpulse96Parms;
    int ierr, k, nfiles;
    bool lsetNewPicks, lusePickFile, lwrtIntFiles;
    char khole[8];
    //------------------------------------------------------------------------// 
    // Initialize
    memset(&event, 0, sizeof(struct prepmtEventParms_struct));
    // Load the event
    ierr = prepmt_event_initializeFromIniFile(iniFile, &event);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error loading event\n", __func__);
        return EXIT_FAILURE;
    }
    // Load the pick model
    ierr = prepmt_prepData_readPickModel(iniFile,
                                         section,
                                         &lsetNewPicks, &lusePickFile,
                                         pickFile, ttimesTableDir, ttimesModel);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error reading pick model\n", __func__);
        return EXIT_FAILURE;
    }
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
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error getting default dt\n", __func__);
        return EXIT_FAILURE;
    }
    // Load the scratch output data info
    ierr = prepmt_prepData_intermediateFileOptions(iniFile, section,
                                                   &lwrtIntFiles,
                                                   wfDir, wfSuffix,
                                                   archiveFile);
    // Read the data
    fprintf(stdout, "%s: Loading data...\n", __func__);
    ierr = prepmt_prepData_readDataListFromIniFile(iniFile,
                                               section,
                                               &nfiles, &sacFiles, &sacpzFiles);
    if (ierr != 0 || nfiles < 1)
    {
        fprintf(stderr, "%s: Failed to find data\n", __func__);
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
            fprintf(stderr, "%s: Failed to read sacFile: %s\n",
                    __func__, sacFiles[k]);
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
                fprintf(stderr, "%s: Failed to read pole-zero information\n",
                        __func__);
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
        fprintf(stderr, "%s: Failed to get processing commands\n", __func__);
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
        fprintf(stderr, "%s: Failed to set event information\n", __func__);
        return EXIT_FAILURE;
    }
    // Verify the data ranges
    fprintf(stdout, "%s: Verifying data ranges...\n", __func__);
    for (k=0; k<nfiles; k++)
    {
        ierr = prepmt_prepData_verifyTeleseismicDistance(lisP, &dmin, &dmax,
                                                         sacData[k]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Cannot verify distance - skipping\n", __func__);
            continue;
        } 
    }
    // Attach picks? 
    if (lsetNewPicks)
    {
        if (!lusePickFile)
        {
            if (lisP)
            {
                fprintf(stderr, "%s: Setting theoretical P picks\n", __func__);
                ierr = prepmt_prepData_setPrimaryPorSPickFromTheoreticalTime(
                                                      NULL, "ak135", DO_P_PICKS,
                                                      SAC_FLOAT_A, SAC_CHAR_KA,
                                                      nfiles, sacData);
            }
            else
            {
                fprintf(stderr, "%s: Setting theoretical S picks\n", __func__);
                ierr = prepmt_prepData_setPrimaryPorSPickFromTheoreticalTime(
                                                      NULL, "ak135", DO_S_PICKS,
                                                      SAC_FLOAT_A, SAC_CHAR_KA,
                                                      nfiles, sacData);
            }
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Failed to set theoretical primary picks\n",
                        __func__);
                return EXIT_FAILURE;
            }
        }
        else
        {
            fprintf(stderr, "%s: Attaching picks from NLL file: %s\n",
                    __func__, pickFile);
            ierr = prepmt_pickFile_nonLinLoc2sac(pickFile, 
                                                 SAC_FLOAT_A, SAC_CHAR_KA,
                                                 nfiles, sacData);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Failed to attached picks to data\n",
                        __func__);
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
            fprintf(stderr, "%s: Failed to modify processing commands\n",
                    __func__);
            return EXIT_FAILURE;
        }
    }
    // Process the data
    fprintf(stdout, "%s: Processing data...\n", __func__);
    time_tic();
    ierr = prepmt_prepData_process(cmds, nfiles, sacData);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error processing data\n", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Processing time: %lf (s)\n", __func__, time_toc());
    if (lwrtIntFiles)
    {
        fprintf(stdout, "%s: Writing intermediate files to directory: %s\n",
                __func__, wfDir);
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
    fprintf(stdout, "%s: Writing H5 archive...\n", __func__);
    ierr = prepmt_prepData_archiveWaveforms(archiveFile, nfiles, sacData);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error writing archive!\n", __func__);
        return EXIT_FAILURE;
    }
    // Free data
    prepmt_commands_freePrepmtCommands(&cmds);
    for (k=0; k<nfiles; k++)
    {
        sacio_free(&sacData[k]);
    } 
    free(sacData);
    return EXIT_SUCCESS;
}
/*!
 * @brief Convenience function to check if SAC data is at a valid
 *        teleseismic distance.
 *
 * @param[in] lisP     If true then this is for P waves. \n
 *                     Otherwise, it is for S waves.
 * @param[in] dminIn   Minimum override distance (degrees).  If NULL then
 *                     this will be set to 30.0.  Otherwise, MAX(30, dminIn)
 *                     will be used.
 * @param[in] dmaxIn   Maximum override distance (degrees).  If NULL then
 *                     this will be set to 95.0.  Otherwise, MIN(95, dmaxIn)
 *                     will be used.
 * @param[in] data     Data with SAC gcarc header variable defined.
 *
 * @result -2 An internal error occurred and the distance couldn't be
 *            verified. \n
 *         -1 The data is out of range. \n
 *          0 The data is in range. 
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_verifyTeleseismicDistance(const bool lisP,
                                              const double *dminIn,
                                              const double *dmaxIn,
                                              const struct sacData_struct data)
{
    double az, baz, dmax, dmin, dist, evla, evlo, gcarc, stla, stlo;
    char stat[8], netw[8], chan[8], loc[8];
    int ierr;
    if (lisP)
    {
        dmin = PREPMT_P_MIN_TELESEISMIC_DIST;
        dmax = PREPMT_P_MAX_TELESEISMIC_DIST;
    }
    else
    {
        dmin = PREPMT_S_MIN_TELESEISMIC_DIST;
        dmax = PREPMT_S_MAX_TELESEISMIC_DIST;
    }
    if (dminIn != NULL){dmin = fmax(*dminIn, dmin);}
    if (dmaxIn != NULL){dmax = fmin(*dmaxIn, dmax);}
    if (dmin > dmax)
    {
        fprintf(stderr, "%s: Error - dmin %f > dmax %f\n",
                 __func__, dmin, dmax);
        return -2;
    }
    ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC, data.header, &gcarc);
    if (ierr != 0)
    {
        ierr = 0;
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA, data.header, &evla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO, data.header, &evlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLA, data.header, &stla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLO, data.header, &stlo);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Unable to compute gcarc and not in header\n",
                    __func__);
            return -2;
        }
        // Compute the azimuth, back-azimuth, distances 
        geodetic_gps2distanceAzimuth(evla, evlo, stla, stlo,
                                     &dist, &gcarc, &az, &baz);
    }
    if (gcarc >= dmin && gcarc <= dmax){return 0;}
    memset(netw, 0, 8*sizeof(char));
    memset(stat, 0, 8*sizeof(char));
    memset(chan, 0, 8*sizeof(char));
    memset(loc,  0, 8*sizeof(char));
    sacio_getCharacterHeader(SAC_CHAR_KNETWK, data.header, netw);
    sacio_getCharacterHeader(SAC_CHAR_KSTNM,  data.header, stat);
    sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data.header, chan);
    sacio_getCharacterHeader(SAC_CHAR_KHOLE,  data.header, loc);
    fprintf(stderr, "%s: %s.%s.%s.%s with distance %f out of bounds [%f,%f]\n",
            __func__, netw, stat, chan, loc, gcarc, dmin, dmax);
    return -1;
}
//============================================================================//
/*!
 * @brief Convenience function to check if SAC data is at a valid
 *        regional distance.
 *
 * @param[in] dminIn   Minimum override distance (degrees).  If NULL then
 *                     this will be set to 1.0.  Otherwise, MAX(1, dminIn)
 *                     will be used.
 * @param[in] dmaxIn   Maximum override distance (degrees).  If NULL then
 *                     this will be set to 15.0.  Otherwise, MIN(15, dmaxIn)
 *                     will be used.
 * @param[in] data     Data with SAC gcarc header variable defined.
 *
 * @result -2 An internal error occurred and the distance couldn't be
 *            verified. \n
 *         -1 The data is out of range. \n
 *          0 The data is in range. 
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_verifyRegionalDistance(const double *dminIn,
                                           const double *dmaxIn,
                                           const struct sacData_struct data)
{
    double az, baz, dmax, dmin, dist, evla, evlo, gcarc, stla, stlo;
    char stat[8], netw[8], chan[8], loc[8];
    int ierr;
    dmin = PREPMT_MIN_REGIONAL_DIST;
    dmax = PREPMT_MAX_REGIONAL_DIST;
    if (dminIn != NULL){dmin = fmax(*dminIn, dmin);}
    if (dmaxIn != NULL){dmax = fmin(*dmaxIn, dmax);}
    if (dmin > dmax)
    {
        fprintf(stderr, "%s: Error - dmin %f > dmax %f\n",
                __func__, dmin, dmax);
        return -2;
    }
    ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC, data.header, &gcarc);
    if (ierr != 0)
    {
        ierr = 0;
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA, data.header, &evla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO, data.header, &evlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLA, data.header, &stla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLO, data.header, &stlo);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Unable to compute gcarc not in header\n",
                    __func__);
            return -2;
        }
        // Compute the azimuth, back-azimuth, distances 
        geodetic_gps2distanceAzimuth(evla, evlo, stla, stlo,
                                     &dist, &gcarc, &az, &baz);
    }
    if (gcarc >= dmin && gcarc <= dmax){return 0;}
    memset(netw, 0, 8*sizeof(char));
    memset(stat, 0, 8*sizeof(char));
    memset(chan, 0, 8*sizeof(char));
    memset(loc,  0, 8*sizeof(char));
    sacio_getCharacterHeader(SAC_CHAR_KNETWK, data.header, netw);
    sacio_getCharacterHeader(SAC_CHAR_KSTNM,  data.header, stat);
    sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data.header, chan);
    sacio_getCharacterHeader(SAC_CHAR_KHOLE,  data.header, loc); 
    fprintf(stderr, "%s: %s.%s.%s.%s with distance %f out of bounds [%f,%f]\n",
            __func__, netw, stat, chan, loc, gcarc, dmin, dmax);
    return -1;
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
int prepmt_prepData_readPickModel(const char *iniFile,
                                  const char *prefix,
                                  bool *lsetNewPicks,
                                  bool *lusePickFile,
                                  char pickFile[PATH_MAX],
                                  char ttimesTableDir[PATH_MAX],
                                  char ttimesModel[128])
{
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
        fprintf(stderr, "%s: Error - ini file %s does not exist\n",
                __func__, iniFile);
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
            fprintf(stderr, "%s: Pick file %s does not exist\n", __func__, s);
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
            fprintf(stderr, "%s: ttimes table directory %s doesn't exist\n",
                    __func__, s);
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
 * @brief Extracts the file/directory information for file writing 
 *        and whether or not to  write individual SAC files.
 *
 * @param[in] iniFile        Name of initialization file.
 * @param[in] section        Section of ini file to read.
 *
 * @param[out] lwrtIntFiles  If true then the intermediate SAC files will be
 *                           written.
 * @param[out] wfDir         Directory to which the SAC files will be written.
 * @param[out] wfSuffix      An identifier to append to the waveform names.
 * @param[out] archiveFile   Name of HDF5 archive file.
 *
 * @result 0 indicates success.
 * 
 */
int prepmt_prepData_intermediateFileOptions(const char *iniFile,
                                            const char *section,
                                            bool *lwrtIntFiles,
                                            char wfDir[PATH_MAX],
                                            char wfSuffix[128],
                                            char archiveFile[PATH_MAX])
{
    const char *s;
    char defaultName[PATH_MAX], vname[256];
    dictionary *ini;
    int ierr;
    size_t lenos;
    //------------------------------------------------------------------------//
    ierr = 0;
    *lwrtIntFiles = false;
    memset(wfDir, 0, PATH_MAX*sizeof(char));
    memset(wfSuffix, 0, 128*sizeof(char));
    memset(archiveFile, 0, PATH_MAX*sizeof(char));
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error - ini file %s does not exist\n",
                __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lwrtIntFiles", section);
    *lwrtIntFiles = iniparser_getboolean(ini, vname, true);
    if (*lwrtIntFiles)
    {
        memset(vname, 0, 256*sizeof(char)); 
        sprintf(vname, "%s:intermediateWaveformDir", section);
        s = iniparser_getstring(ini, vname, "./\0");
        if (!os_path_isdir(s))
        {
            ierr = os_makedirs(s);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Couldn't make output directory %s\n",
                        __func__, s);
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
 * @brief Reads the SAC data files and SAC pole-zero files from the
 *        the iniFile.
 *
 * @param[in] iniFile     Name of ini file.
 * @param[in] prefix      This is the `directory' in the ini file where
 *                        the data list exists.  For example:
 *                        parmt:prepTeleseismicData.
 *
 * @param[out] nfiles     Number of data files read.
 * @param[out] sacFiles   The nfiles SAC files.
 * @param[out] sacpzFiles The nfiles SAC pole-zero files.  If the ifile'th
 *                        file is empty then no pole-zero file was specified
 *                        for this SAC file.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_readDataListFromIniFile(const char *iniFile,
                                            const char *prefix,
                                            int *nfiles,
                                            char ***sacFiles,
                                            char ***sacpzFiles)
{
    FILE *flist;
    const char *s;
    char **csplit, **sFiles, **spzFiles, vname[256], cline[2*PATH_MAX];
    dictionary *ini;
    int ierr, item, k, nitems, nread;
    bool luseDataList;
    //------------------------------------------------------------------------//
    ierr = 0;
    *nfiles = 0;
    sFiles = NULL;
    spzFiles = NULL;
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file %s doesn't exist\n",
                __func__, iniFile);
        return -1;
    }    
    // Load the ini file
    ini = iniparser_load(iniFile);
    // Determine if i'm using a data list or not
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:luseDataList", prefix);
    luseDataList = iniparser_getboolean(ini, vname, false);
    // Read from external file list
    if (luseDataList)
    {
        // Get the file list name and verify it exists
        memset(vname, 0, 256*sizeof(char)); 
        sprintf(vname, "%s:dataList", prefix);
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            fprintf(stderr, "%s: Error file list %s does not exist\n",
                    __func__, s);
            return -1;
        }
        // Count the number of files
        flist = fopen(s, "r");
        while(fgets(cline, 128, flist) != NULL){*nfiles = *nfiles + 1;}
        if (*nfiles < 1)
        {
            fprintf(stderr, "%s: Error no files in file list\n", __func__);
            fclose(flist);
        }
        // Set space
        nread = 0;
        sFiles = (char **) calloc((size_t) *nfiles, sizeof(char *));
        spzFiles = (char **) calloc((size_t) *nfiles, sizeof(char *));
        rewind(flist);
        // Load the number of files
        for (k=0; k<*nfiles; k++)
        {
            memset(cline, 0, 2*PATH_MAX*sizeof(char));
            fgets(cline, 2*PATH_MAX, flist);
            csplit = string_rsplit(NULL, cline, &nitems);
            if (csplit != NULL)
            {
                if (nitems > 0)
                {
                    if (!os_path_isfile(csplit[0]))
                    {   
                        fprintf(stderr, "%s: Error data file %s doesnt exist\n",
                                __func__, csplit[0]);
                        goto NEXT_LINE;
                    }
                    else
                    {   
                        sFiles[nread] = (char *)
                                        calloc(PATH_MAX, sizeof(char));
                        strcpy(sFiles[nread], csplit[0]);
                        nread = nread + 1;
                    }
                }
                if (nitems == 2)
                {
                    if (!os_path_isfile(csplit[1]))
                    {
                        fprintf(stdout, "%s: Warning pz file %s doesnt exist\n",
                                __func__, csplit[1]);
                    }
                    else
                    {   
                        spzFiles[nread-1] = (char *)
                                            calloc(PATH_MAX, sizeof(char));
                        strcpy(spzFiles[nread-1], csplit[1]);
                    }
                }
                NEXT_LINE:;
                for (item=0; item<nitems; item++){free(csplit[item]);}
                free(csplit);
            }
        } // Loop on files in file list
        fclose(flist);
        if (nread != *nfiles)
        {
            fprintf(stdout, "%s: Only %d files from %s read\n",
                    __func__, nread, s);
        }
        *nfiles = nread;
    }
     // Read from internal file list
    else
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:nobs", prefix);
        *nfiles = iniparser_getint(ini, vname, 0);
        if (*nfiles < 1)
        {
            fprintf(stderr, "%s: No files to read!\n", __func__);
            ierr = 1;
            goto ERROR;
        }
        nread = 0;
        sFiles = (char **) calloc((size_t) *nfiles, sizeof(char *));
        spzFiles = (char **) calloc((size_t) *nfiles, sizeof(char *));
        for (k=0; k<*nfiles; k++)
        {
            memset(vname, 0, 256*sizeof(char));
            sprintf(vname, "%s:sacFile_%d", prefix, k+1);
            sFiles[k] = (char *) calloc(PATH_MAX, sizeof(char));
            spzFiles[k] = (char *) calloc(PATH_MAX, sizeof(char));
            s = iniparser_getstring(ini, vname, NULL);
            if (!os_path_isfile(s))
            {
                fprintf(stderr, "%s: Error sac file %s doesn't exist\n", 
                        __func__, s);
                continue;
            }
            strcpy(sFiles[nread], s);
            memset(vname, 0, 256*sizeof(char));
            sprintf(vname, "%s:sacpzFile_%d", prefix, k+1);
            s = iniparser_getstring(ini, vname, NULL);
            if (os_path_isfile(s))
            {
                strcpy(spzFiles[nread], s);
            }
            nread = nread + 1;
        }
    }
    if (*nfiles < 1)
    {
        fprintf(stderr, "%s: Error - no data files read\n", __func__);
        ierr = 1;
    }
ERROR:;
    *sacFiles = sFiles; 
    *sacpzFiles = spzFiles; 
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Reads the default sampling period and window information for the
 *        estimation.
 *
 * @param[in] iniFile     Name of ini file.
 * @param[in] section     Section of ini file to read.
 *
 * @param[out] targetDt   Target sampling period (s) for all data in inversion.
 * @param[out] cutStart   Default start window time (s) relative to first
 *                        arrival.
 * @param[out] cutEnd     Default end window time (s) relative to first arrival.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_getDefaultDTAndWindowFromIniFile(const char *iniFile,
                                                     const char *section,
                                                     double *targetDt,
                                                     double *cutStart,
                                                     double *cutEnd)
{
    char vname[256];
    dictionary *ini;
    *cutStart =-2.0;
    *cutEnd = 3.0;  
    *targetDt = 0.1;
    if (!os_path_isfile(iniFile))
    {    
        fprintf(stderr, "%s: Error ini file %s doesn't exist\n",
                __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:cutStart", section);
    *cutStart = iniparser_getdouble(ini, vname, -2.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:cutEnd", section);
    *cutEnd   = iniparser_getdouble(ini, vname, 3.0);
    if (*cutEnd <= *cutStart)
    {    
        fprintf(stderr, "%s: Error cutStart: %f must be less than cutEnd: %f\n",
                __func__, *cutStart, *cutEnd);
        return -1; 
    }

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:targetDt", section);
    *targetDt = iniparser_getdouble(ini, vname, 0.1);
    if (*targetDt <= 0.0)
    {    
        fprintf(stderr, "%s: Error targetDt %f must be positive\n",
                __func__, *targetDt);
        return -1; 
    }
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the event latitude, longitude, depth, and event time on the
 *        SAC headers.  Given this
 *
 * @param[in] evla         Event latitude (degrees).
 * @param[in] evlo         Event longitude (degrees).
 * @param[in] evdp         Event depth (km).
 * @param[in] evtime       Event time (UTC-seconds).
 * @param[in] nobs         Number of observations.
 *
 * @param[in,out] data     On input contains the SAC data and the station
 *                         latitude and longitude.
 *                         On output contains the event latitude, longitude,
 *                         depth, and derivative information such as the
 *                         great circle distance, distance, azimuth,
 *                         and back-azimuth.
 *                         This is an array of dimension [nobs].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_setEventInformation(const double evla,
                                        const double evlo,
                                        const double evdp,
                                        const double evtime,
                                        const int nobs,
                                        struct sacData_struct *data)
{
    double az, baz, dist, epoch, gcarc, o, stla, stlo, var;
    int ierr, ierr1, k;
    ierr = 0;
    if (data == NULL || nobs < 1)
    {
        fprintf(stderr, "%s: Error - no data\n", __func__);
        return -1; 
    }
    for (k=0; k<nobs; k++)
    {
        // Figure out the origin time
        ierr = sacio_getEpochalStartTime(data[k].header, &epoch);
        o = evtime - epoch;
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get start time\n", __func__);
            break;
        }
        // Get the station location
        ierr = sacio_getFloatHeader(SAC_FLOAT_STLA, data[k].header, &stla);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get station latitude\n", __func__);
            break;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_STLO, data[k].header, &stlo);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get station longitude\n", __func__);
            break;
        }
        // Compute the azimuth, back-azimuth, distances 
        geodetic_gps2distanceAzimuth(evla, evlo, stla, stlo,
                                     &dist, &gcarc, &az, &baz);
        // Let user know what's going to be overwritten
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_EVLA, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - evla) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event latitude: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_EVLO, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - evlo) > 1.e-3)
        {   
            fprintf(stdout, "%s: Overwriting event longitude: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_EVDP, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - evdp) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event depth: %f\n",
                     __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_GCARC, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - gcarc) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event great-circle dist: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_DIST, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - dist) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event distance: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_AZ, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - az) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event azimuhh: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_BAZ, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - baz) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event back-azimuth: %f\n",
                    __func__, var);
        }
        ierr1 = sacio_getFloatHeader(SAC_FLOAT_O, data[k].header, &var);
        if (ierr1 == 0 && fabs(var - o) > 1.e-3)
        {
            fprintf(stdout, "%s: Overwriting event origin time: %f\n",
                    __func__, var);
        }
        // Set the variables
        sacio_setFloatHeader(SAC_FLOAT_EVLA,  evla,  &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_EVLO,  evlo,  &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_EVDP,  evdp,  &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_GCARC, gcarc, &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_DIST,  dist,  &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_AZ,    az,    &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_BAZ,   baz,   &data[k].header);
        sacio_setFloatHeader(SAC_FLOAT_O,     o,     &data[k].header);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Loads the H5 archived waveforms.
 *
 * @param[in] archiveFile    Name of archive file.
 * @param[out] nobs          Number of observations read.
 * @param[out] ierr          0 indicates success.
 *
 * @result An array of SAC structures with the input waveforms from the
 *         archive file's /ObservedWaveforms directory.  The array is of
 *         dimension [nobs] and can be freed with free.
 *
 * @author Ben Baker
 *
 */
struct sacData_struct *prepmt_prepData_readArchivedWaveforms(
    const char *archiveFile, int *nobs, int *ierr)
{
    char **sacFiles;
    hid_t groupID, fileID;
    int i, nfiles;
    struct sacData_struct *sacData;
    *nobs = 0; 
    sacData = NULL;
    if (!os_path_isfile(archiveFile))
    {
        fprintf(stderr, "%s: Error archive file %s doesn't exist\n",
                __func__, archiveFile);
        *ierr = 1; 
        return sacData;
    }    
    fileID = H5Fopen(archiveFile, H5F_ACC_RDONLY, H5P_DEFAULT); 
    groupID = H5Gopen2(fileID, "/ObservedWaveforms", H5P_DEFAULT);
    sacFiles = sacioh5_getFilesInGroup(groupID, &nfiles, ierr);
    if (*ierr != 0 || sacFiles == NULL)
    {
        fprintf(stderr, "%s: Error getting names of SAC flies\n", __func__);
        *ierr = 1; 
        return sacData;
    }    
    sacData = sacioh5_readTimeSeriesList(nfiles, (const char **) sacFiles,
                                         groupID, nobs, ierr);
    if (*ierr != 0)
    {    
        fprintf(stderr, "%s: Errors while reading SAC files\n", __func__);
    }    
    if (*nobs != nfiles)
    {
        fprintf(stdout, "%s: Warning - subset of data was read\n", __func__);
    }
    // Clean up and close archive file
    for (i=0; i<nfiles; i++)
    {
        if (sacFiles[i] != NULL){free(sacFiles[i]);}
    }
    free(sacFiles);
    H5Gclose(groupID);
    H5Fclose(fileID);
    return sacData;
}
//============================================================================//
/*!
 * @brief Writes the metadata and processed data to an H5 archive.
 *
 * @param[in] archiveFile   Name of HDF5 archive file.
 * @param[in] nobs          Number of observation.
 * @param[in] data          Observed SAC waveforms.  This is an array of
 *                          dimension [nobs].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */
int prepmt_prepData_archiveWaveforms(const char *archiveFile,
                                     const int nobs,
                                     const struct sacData_struct *data)
{
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
        fprintf(stderr, "%s: Error archive file is NULL\n", __func__);
        return -1; 
    }   
    if (strlen(archiveFile) == 0)
    {
        fprintf(stderr, "%s: Error archive file is blank\n", __func__);
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
                fprintf(stderr, "%s: Failed to make output dircectory: %s\n",
                        __func__, dir);
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
        fprintf(stderr, "%s: Clobbering file %s\n", __func__, archiveFile);
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
            fprintf(stderr, "%s: Error writing waveform %d, object %sn",
                    __func__, k + 1, objName);
            goto ERROR;
        }
    }
ERROR:;
    H5Gclose(groupID);
    H5Fclose(fileID);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the theoretical first arriving P pick time or S pick time.
 *
 * @param[in] dirnm     Directory containing the ttimes models.  If NULL this
 *                      will use the default as specified in the ttimes config
 *                      header file.
 * @param[in] model     ttimes model.  Likely "ak135" or "iasp91".
 * @param[in] ldop      If true then compute the primary P pick.
 *                      Otherwise compute the primary S pick.
 * @param[in] data      Data with event depth, event origin time, and the
 *                      source-to-receiver great circle distance for each 
 *                      observation.
 *
 * @param[out] ptimes   Epochal theoretical primary P or S arrival times
 *                      (UTC-seconds) for each observation.  This has dimension
 *                      [nobs].
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_computeTheoreticalPorSPickTimes(
    const char *dirnm, const char *model,
    const bool ldop,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes)
{
    struct ttimesTravelTime_struct ppick;
    double delta, depth, epoch, o;
    int ierr, k;
    if (nobs < 1 || ptimes == NULL)
    {
        fprintf(stderr, "%s: Insufficient space for output ptimes\n", __func__);
        return -1;
    }
    memset(&ppick, 0, sizeof(struct ttimesTravelTime_struct));
    for (k=0; k<nobs; k++)
    {
        ierr = sacio_getEpochalStartTime(data[k].header, &epoch);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get start time\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_O, data[k].header, &o);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Origin time not set\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                    data[k].header, &delta);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error event distance not set\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVDP,
                                    data[k].header, &depth);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error event depth not set\n", __func__);
            return -1;
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting origin time\n", __func__);
            return -1;
        }
        if (ldop)
        {
            ierr = ttimes_getFirstPPhase(delta, depth, dirnm, model, &ppick);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error computing P pick time\n", __func__);
                return -1;
            } 
        }
        else
        {
            ierr = ttimes_getFirstSPhase(delta, depth, dirnm, model, &ppick);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error computing S pick time\n", __func__);
                return -1;
            }
        }
        ptimes[k] = epoch + o + ppick.tt;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the theoretical first arriving P pick times.
 *
 * @param[in] dirnm     Directory containing the ttimes models.  If NULL this
 *                      will use the default as specified in the ttimes config
 *                      header file.
 * @param[in] model     ttimes model.  Likely "ak135" or "iasp91".
 * @parma[in] nobs      Number of observations.
 * @param[in] data      Data with event depth, event origin time, and the
 *                      source-to-receiver great circle distance for each 
 *                      observation.  This is an array with dimension [nobs].
 *
 * @param[out] ptimes   Epochal theoretical primary P arrival times
 *                      (UTC-seconds) for each observation.  This has dimension
 *                      [nobs].
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_computeTheoreticalPPickTimes(
    const char *dirnm, const char *model,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes)
{
    int ierr;
    ierr = prepmt_prepData_computeTheoreticalPorSPickTimes(dirnm, model, true,
                                                           nobs, data, ptimes);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing P pick times\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the theoretical first arriving S pick times.
 *
 * @param[in] dirnm     Directory containing the ttimes models.  If NULL this
 *                      will use the default as specified in the ttimes config
 *                      header file.
 * @param[in] model     ttimes model.  Likely "ak135" or "iasp91".
 * @param[in] data      Data with event depth, event origin time, and the
 *                      source-to-receiver great circle distance for each 
 *                      observation.
 *
 * @param[out] ptimes   Epochal theoretical primary S arrival times
 *                      (UTC-seconds) for each observation.  This has dimension
 *                      [nobs]
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_computeTheoreticalSPickTimes(
    const char *dirnm, const char *model,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes)
{
    int ierr;
    ierr = prepmt_prepData_computeTheoreticalPorSPickTimes(dirnm, model, false,
                                                          nobs, data, ptimes);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error computing S pick times\n", __func__);
        return -1; 
    }   
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the primary pick time on the header from the theoretical
 *        first arriving P or S travel time.
 *
 * @param[in] dirnm           Directory name where ttimes models reside.
 *                            If NULL then the default directory will be used.
 * @param[in] model           Earth model to be used in travel time computation.
 *                            If NULL then ak135 will be used.
 * @param[in] ldop            If true then set the first arriving P travel time.
 *                            Otherwise set the first arriving S travel time. 
 * @param[in] pickHeaderTime  Header variable identifier to hold the pick
 *                            time.  Likely SAC_FLOAT_A for a primary arrival,
 *                            but also may be SAC_FLOAT_T0-SAC_FLOAT_T9
 *                            for other phases.
 * @param[in] pickHeaderName  Header variable identifer to hold the phase name.
 *                            Likely SAC_CHAR_KA for a primary arrival,
 *                            but also may be SAC_CHAR_KT0-SAC_CHAR_KT9
 *                            for other phases.
 * @param[in] nobs            Number of observations.
 *
 * @param[in,out] data        On input holds the data with the event depth,
 *                            time, and great circle distance for all 
 *                            observations.
 *                            On output holds the theoretical first arriving
 *                            P or S travel time.
 *                            This is an array of dimension [nobs].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_setPrimaryPorSPickFromTheoreticalTime(
    const char *dirnm, const char *model,
    const bool ldop,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    const int nobs, struct sacData_struct *data)
{
    double *ptimes;
    char phaseName[8];
    int ierr, k;
    //------------------------------------------------------------------------//
    //   
    // Initialize and get number of observations
    ierr = 0; 
    ptimes = NULL;
    memset(phaseName, 0, 8*sizeof(char));
    phaseName[0] = 'P';
    if (nobs < 1 || data == NULL) 
    {    
        fprintf(stderr, "%s: Error - no data\n", __func__);
        return -1;
    }
    if (!ldop){phaseName[0] = 'S';}
    // Compute the theoretical primary arrival times
    ptimes = memory_calloc64f(nobs);
    ierr = prepmt_prepData_computeTheoreticalPorSPickTimes(dirnm, model, ldop,
                                                          nobs, data,
                                                          ptimes);
    if (ierr != 0)
    {
        fprintf(stderr,
                "%s: Error computing theorietical primary arrival times\n",
                __func__);
        goto ERROR;
    }
    // Set the theoretical arrival times on the header
    for (k=0; k<nobs; k++)
    {
        ierr = sacio_setEpochalPickOnHeader(ptimes[k], phaseName,
                                            pickHeaderTime,
                                            pickHeaderName,
                                            &data[k].header);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set pick time on header\n",
                    __func__);
            goto ERROR;
        }
    }
ERROR:;
    memory_free64f(&ptimes);
    return ierr;
}
//============================================================================//
int prepmt_prepData_setTheoreticalSurfaceWaveArrivalTime(
    const double vel, const bool lr,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    const int nobs, struct sacData_struct *data)
{
    char phaseName[8];
    double dist, epoch, o, ptime;
    int ierr, k;
    // Initialize and get number of observations
    ierr = 0;
    memset(phaseName, 0, 8*sizeof(char));
    phaseName[0] = 'R';
    if (nobs < 1 || data == NULL)
    {
        fprintf(stderr, "%s: Error - no data\n", __func__);
        return -1;
    }
    if (vel <= 0.0)
    {
        fprintf(stderr, "%s: Error - velocity %f must be positive\n",
                __func__, vel);
        return -1;
    }
    if (!lr){phaseName[0] = 'L';}
    for (k=0; k<nobs; k++)
    {
        ierr = sacio_getEpochalStartTime(data[k].header, &epoch);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get start time\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_O, data[k].header, &o);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Origin time not set\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_DIST,
                                    data[k].header, &dist);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error event distance not set\n", __func__);
            return -1;
        }
        ptime = epoch + o + dist/vel;
        ierr = sacio_setEpochalPickOnHeader(ptime, phaseName,
                                            pickHeaderTime,
                                            pickHeaderName,
                                            &data[k].header);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set pick time on header\n",
                     __func__);
            goto ERROR;
        }
    }
ERROR:;
    return ierr;
}
//============================================================================//
/*!
 * @brief Processes the data.
 *
 * @param[in] cmds      String based commands for data processing corresponding
 *                      to each observation.
 * @param[in] nobs      Number of observations.
 *
 * @param[in,out] data  On input contains the data to process.
 *                      On output contains the
 *                      This is an array of dimension [nobs].
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_prepData_process(const struct prepmtCommands_struct cmds,
                            const int nobs, struct sacData_struct *data)
{
    struct serialCommands_struct *commands;
    double dt, dt0, epoch, epoch0, time, *ycopy;
    int *nyAll, i, i0, ierr, k, npts0, nq, nwork;
    bool lnewDt, lnewStartTime;
    const int nTimeVars = 12;
    const enum sacHeader_enum timeVars[12]
       = {SAC_FLOAT_A, SAC_FLOAT_O, 
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    // No data
    ycopy = NULL;
    nyAll = NULL;
    if (nobs < 1){return 0;}
    if (cmds.nobs != nobs)
    {
        fprintf(stderr, "%s: Inconsistent sizes between cmds and nobs\n",
                __func__);
        return -1;
    }
    if (data == NULL)
    {
        fprintf(stderr, "%s: Error - data is NULL\n", __func__);
        return -1;
    }
    commands = (struct serialCommands_struct *)
               calloc((size_t) nobs, sizeof(struct serialCommands_struct));
    // Set the processing commands and the data
    for (k=0; k<nobs; k++) 
    {
        ierr = process_stringsToSerialCommandsOptions(cmds.cmds[k].ncmds,
                                      (const char **) cmds.cmds[k].cmds,
                                      &commands[k]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting serial command string\n",
                    __func__);
            goto ERROR;
        }
        process_setSerialCommandsData64f(data[k].npts,
                                         data[k].data, 
                                         &commands[k]);
    }
    // Process the data
    ierr = process_applyMultipleSerialCommands(nobs, commands);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error applying serial commands chains to data\n",
                __func__);
        goto ERROR;
    }    
    // Get the workspace
    nwork =-1; 
    nyAll = memory_calloc32i(nobs);
    for (k=0; k<nobs; k++)
    {
        process_getSerialCommandsData64f(commands[k], -1, &nyAll[k], NULL);
        nwork = MAX(nwork, nyAll[k]);
    }
    ycopy = memory_calloc64f(nwork);
    // Move the results back onto data
    for (k=0; k<nobs; k++)
    {
        lnewDt = false;
        lnewStartTime = false;
        // Loop through commands and check if dt or npts changed
        ierr = sacio_getEpochalStartTime(data[k].header, &epoch0);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get start time of trace\n",
                    __func__);
            goto ERROR;
        }
        epoch = epoch0;
        sacio_getFloatHeader(SAC_FLOAT_DELTA, data[k].header, &dt0);
        sacio_getIntegerHeader(SAC_INT_NPTS, data[k].header, &npts0);
        dt = dt0;
        for (i=0; i<commands[k].ncmds; i++)
        {
            if (commands[k].commands[i].type == CUT_COMMAND)
            {
                i0 = commands[k].commands[i].cut.i0;
                epoch = epoch + (double) i0*dt;
                lnewStartTime = true;
            }
            if (commands[k].commands[i].type == DOWNSAMPLE_COMMAND)
            {
                nq = commands[k].commands[i].downsample.nq;
                dt = dt*(double) nq;
                lnewDt = true;
            }
            if (commands[k].commands[i].type == DECIMATE_COMMAND)
            {
                nq = commands[k].commands[i].decimate.nqAll;
                dt = dt*(double) nq;
                lnewDt = true;
            }
        }
        // Extract the data - here a resize is required
        if (nyAll[k] != npts0)
        {
            ierr = process_getSerialCommandsData64f(commands[k], nwork,
                                                    &nyAll[k], ycopy);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error extracting data onto ycopy\n",
                        __func__);
                goto ERROR;
            }
            sacio_freeData(&data[k]);
            if (nyAll[k] > 0)
            {
                sacio_freeData(&data[k]);
                data[k].data = sacio_malloc64f(nyAll[k]);
                data[k].npts = nyAll[k];
                sacio_setIntegerHeader(SAC_INT_NPTS, nyAll[k],
                                       &data[k].header);
                ierr = array_copy64f_work(nyAll[k], ycopy, data[k].data);
            }
        }
        // Otherwise just copy back onto the data structure
        else
        {
            ierr = process_getSerialCommandsData64f(commands[k], nwork,
                                                    &nyAll[k],
                                                    data[k].data);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error extracting data\n", __func__);
                goto ERROR;
            }
        }
        // Update the times
        if (lnewStartTime)
        {
            // Update the picks
            for (i=0; i<nTimeVars; i++)
            {
                // Change all times to be relative to epoch0
                ierr = sacio_getFloatHeader(timeVars[i], data[k].header,
                                            &time);
                if (ierr == 0)
                {
                    time = time + epoch0; // Turn to real time
                    time = time - epoch;  // Make relative to new time 
                    //printf("new d %d %f %e %e\n", i, time, epoch, epoch0);
                    //if (epoch < epoch0){printf("this is broken\n");}
                    sacio_setFloatHeader(timeVars[i], time,
                                         &data[k].header);
                }
                ierr = 0;
            }
            sacio_setEpochalStartTime(epoch, &data[k].header);
        }
        // Update the sampling period
        if (lnewDt)
        {
            sacio_setFloatHeader(SAC_FLOAT_DELTA, dt, &data[k].header);
        }
     }
ERROR:;
    for (k=0; k<nobs; k++)
    {
        process_freeSerialCommands(&commands[k]);
    }
    memory_free64f(&ycopy);
    memory_free32i(&nyAll);
    free(commands);
    return ierr;
}
