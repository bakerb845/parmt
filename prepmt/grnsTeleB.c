#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <iniparser.h>
#include <limits.h>
#include "prepmt/prepmt.h"
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#define PROGRAM_NAME "xgrnsTeleB"

static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256]);
struct sacData_struct *
    prepmt_grnsTeleB_readData(const char *archiveFile, int *nobs, int *ierr);

int prepmt_grnsTeleB_readParameters(const char *iniFile,
                                    const char *section,
                                    char archiveFile[PATH_MAX],
                                    char parmtDataFile[PATH_MAX],
                                    bool *luseCrust1,
                                    char crustDir[PATH_MAX],
                                    bool *luseSourceModel,
                                    char sourceModel[PATH_MAX],
                                    bool *luseTstarTable,
                                    double *defaultTstar,
                                    char tstarTable[PATH_MAX],
                                    bool *lrepickGrns,
                                    double *staWin, double *ltaWin,
                                    double *staltaThreshPct,
                                    bool *lalignXC,
                                    bool *luseEnvelope,
                                    bool *lnormalizeXC,
                                    double *maxXCtimeLag,
                                    int *ndepth, double **depths);
int prepmt_grnsTeleB_loadTstarTable(const char *tstarTable,
                                    const double defaultTstar,
                                    const int nobs,
                                    struct sacData_struct *data,
                                    double **tstars);

int main(int argc, char **argv)
{
    char iniFile[PATH_MAX], archiveFile[PATH_MAX], tstarTable[PATH_MAX],
         parmtDataFile[PATH_MAX], crustDir[PATH_MAX],
         sourceModel[PATH_MAX], section[256];
    struct prepmtCommands_struct cmds;
    struct vmodel_struct *recmod, telmod, srcmod;
    struct prepmtEventParms_struct event;
    struct hudson96_parms_struct hudson96Parms;
    struct hpulse96_parms_struct hpulse96Parms;
    struct prepmtModifyCommands_struct options;
    struct sacData_struct *sacData, *grns, *ffGrns, *locFF;
    const char *hudsonSection = "hudson96\0";
    const char *hpulseSection = "hpulse96\0";
    const int ntstar = 1;
    double *depths, *tstars, defaultTstar, ltaWin, maxXCtimeLag,
           staltaThreshPct, staWin;
    int i, ierr, iobs, k, kndx, ndepth, nobs;
    bool lalignXC, lnormalizeXC, lrepickGrns, luseCrust1,
         luseSourceModel, luseEnvelope, luseTstarTable;

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
    // Read the forward modeling parameters
    ierr = prepmt_hudson96_readHudson96Parameters(iniFile, hudsonSection,
                                                  &hudson96Parms);
    if (ierr != 0)
    {
        printf("%s: Failed to read hudson parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = prepmt_hpulse96_readHpulse96Parameters(iniFile, hpulseSection,
                                                  &hpulse96Parms);
    if (ierr != 0)
    {
        printf("%s: Failed to read hpulse parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = prepmt_grnsTeleB_readParameters(iniFile, section,
                                           archiveFile,
                                           parmtDataFile,
                                           &luseCrust1, crustDir,
                                           &luseSourceModel, sourceModel,
                                           &luseTstarTable,
                                           &defaultTstar,
                                           tstarTable,
                                           &lrepickGrns,
                                           &staWin, &ltaWin, &staltaThreshPct,
                                           &lalignXC, &luseEnvelope,
                                           &lnormalizeXC, &maxXCtimeLag,
                                           &ndepth, &depths);
    if (ierr != 0)
    {
        printf("%s: Failed to read modeling parameters\n", PROGRAM_NAME);
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
//memset(archiveFile, 0, PATH_MAX*sizeof(char));
//strcpy(archiveFile, "windowedPData/observedWaveforms.h5\0");
    printf("%s: Loading data...\n", PROGRAM_NAME);
    sacData = prepmt_prepData_readArchivedWaveforms(archiveFile, &nobs, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error loading data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    if (luseTstarTable)
    {
        printf("%s: Loading tstar table...\n", PROGRAM_NAME);
        ierr = prepmt_grnsTeleB_loadTstarTable(tstarTable, defaultTstar,
                                               nobs, sacData, &tstars);
        if (ierr != 0)
        {
            printf("%s: Failed to load t* table\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
    }
    else
    {
        printf("%s: Using default tstar: %f\n", PROGRAM_NAME, defaultTstar); 
        tstars = array_set64f(nobs, defaultTstar, &ierr);
    }
    printf("%s: Reading pre-processing commands...\n", PROGRAM_NAME);
    cmds = prepmt_commands_readFromIniFile(iniFile, section, nobs,
                                           sacData, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error loading grns prep commands\n", PROGRAM_NAME);
        return -1;
    }
    printf("%s: Loading velocity models...\n", PROGRAM_NAME);
    recmod = (struct vmodel_struct *)
             calloc((size_t) nobs, sizeof(struct vmodel_struct));
    ierr = hudson96_getModels(nobs, sacData,
                              luseCrust1, crustDir,
                              luseSourceModel, sourceModel,
                              &telmod, &srcmod, recmod);
    if (ierr != 0)
    {
        printf("%s: Failed to load models\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    printf("%s: Computing fundamental fault solutions...\n", PROGRAM_NAME);
    ffGrns = (struct sacData_struct *)
             calloc((size_t) (ndepth*nobs*ntstar*10),
                    sizeof(struct sacData_struct));
    for (k=0; k<nobs; k++)
    {
        locFF = prepmt_hudson96_computeGreensFF(hudson96Parms,
                                                hpulse96Parms,
                                                telmod, srcmod, recmod,
                                                ntstar, &tstars[1],
                                                ndepth, depths,
                                                1, &sacData[k], &ierr);
        for (i=0; i<10*ndepth; i++)
        {
            kndx = k*10*ndepth + i;
            sacio_copy(locFF[i], &ffGrns[kndx]);
            sacio_freeData(&locFF[i]);
        }
        free(locFF);
    }
    printf("%s: Contextualing Green's functions...\n", PROGRAM_NAME);
    grns = (struct sacData_struct *)
           calloc((size_t) (6*ndepth*nobs*ntstar),
                  sizeof(struct sacData_struct));
    ierr = prepmt_greens_ffGreensToGreens(nobs, sacData,
                                          ndepth, ntstar, ffGrns, grns);
    if (ierr != 0)
    {
        printf("%s: Error contextualizing Greens functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    for (k=0; k<10*nobs*ndepth*ntstar; k++){sacio_freeData(&ffGrns[k]);}
    free(ffGrns);
    // Release memory
    for (k=0; k<nobs; k++)
    {
        cps_utils_freeVmodelStruct(&recmod[k]);
    }
    free(recmod);
    cps_utils_freeVmodelStruct(&srcmod);
    cps_utils_freeVmodelStruct(&telmod);
    // Process the Green's functions
    printf("%s: Modifying Green's functions processing commands...\n",
           PROGRAM_NAME);
    for (k=0; k<cmds.nobs; k++)
    {
        ierr = prepmt_commands_modifyCommandsCharsStruct(options,
                                                         sacData[k],
                                                         &cmds.cmds[k]);
        /*
        for (int l=0; l<cmds.cmds[k].ncmds; l++)
        {
             printf("%s\n", cmds.cmds[k].cmds[l]);
        }
        printf("\n");
        */
        if (ierr != 0)
        {
            printf("%s: Failed to modify processing commands\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
    }
    printf("%s: Processing Green's functions...\n", PROGRAM_NAME);
    ierr = prepmt_greens_processHudson96Greens(nobs, ntstar, ndepth,
                                               cmds, grns);
    if (ierr != 0)
    {
        printf("%s: Error processing Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Repick the Green's functions with an STA/LTA
    if (lrepickGrns)
    {
        printf("%s: Repicking Green's functions...\n", PROGRAM_NAME);
        for (k=0; k<ndepth*ntstar*nobs*6; k=k+6)
        {
            ierr = prepmt_greens_repickGreensWithSTALTA(staWin, ltaWin,
                                                        staltaThreshPct,
                                                        &grns[k]);
            if (ierr != 0)
            {
                printf("%s: Error repicking Green's functions\n", PROGRAM_NAME);
                return EXIT_FAILURE;
            }
        }
    }
    // Shift 
    if (lalignXC)
    {
        printf("%s: Refining waveform alignment with cross-correlation\n",
               PROGRAM_NAME);
        for (iobs=0; iobs<nobs; iobs++)
        {
            for (k=0; k<ndepth*ntstar; k++)
            {
                kndx = iobs*ndepth*ntstar*6 + k*6;
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
                                             ndepth, ntstar, grns);
    if (ierr != 0)
    {
        printf("%s: Error trimming Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Dump archive
/*
sacio_writeTimeSeriesFile("test_gxx.sac", grns[11*ndepth*6+6]);
sacio_writeTimeSeriesFile("test_gyy.sac", grns[11*ndepth*6+7]);
sacio_writeTimeSeriesFile("test_gzz.sac", grns[11*ndepth*6+8]);
sacio_writeTimeSeriesFile("test_gxy.sac", grns[11*ndepth*6+9]);
sacio_writeTimeSeriesFile("test_gxz.sac", grns[11*ndepth*6+10]);
sacio_writeTimeSeriesFile("test_gyz.sac", grns[11*ndepth*6+11]);
*/
    printf("%s: Writing archive: %s\n", PROGRAM_NAME, parmtDataFile);
    ierr = prepmt_greens_writeArchive(parmtDataFile, //"./", "pwave",
                                      nobs, ndepth,
                                      event.latitude, event.longitude,
                                      depths, sacData, grns);
    if (ierr != 0)
    {
        printf("%s: Failed to write archive file\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // release rest of memory
    for (k=0; k<nobs; k++){sacio_freeData(&sacData[k]);}
    for (k=0; k<ndepth*ntstar*nobs*6; k++){sacio_freeData(&grns[k]);}
    memory_free64f(&tstars);
    free(grns);
    free(sacData);
    iscl_finalize();

    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Obtains the forward modeling t* values for each teleseismic
 *        waveform observation from a text file.
 *
 * @param[in] tstarTable    Name of tstar table.
 * @param[in] defaultTstar  Default t* if a piece of data cannot be found in
 *                          the table.
 * @param[in] nobs          Number of observations.
 * @param[in] data          List of teleseismic waveforms to map t* values to.
 *                          This is an array of dimension [nobs].
 *
 * @param[out] tstars       tstar values corresponding to each station.
 *                          If a SNCL in the data list cannot be found in the
 *                          tstarTable file then this tstar will be set to the
 *                          default value.
 *                          This is an array of length [nobs].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */
int prepmt_grnsTeleB_loadTstarTable(const char *tstarTable,
                                    const double defaultTstar,
                                    const int nobs,
                                    struct sacData_struct *data,
                                    double **tstars)
{
    FILE *tsf;
    char cline[64], search[64], item[64];
    double *tstar, tin;
    int i, ierr, k, nlines;
    bool lfound;
    *tstars = NULL;
    tstar = NULL;
    if (!os_path_isfile(tstarTable))
    {
        fprintf(stderr, "%s: Error tstar table %s does not exist\n",
                 __func__, tstarTable);
        return -1;
    }
    if (nobs < 1 || data == NULL)
    {
        fprintf(stderr, "%s: No data\n", __func__);
        return -1;
    }
    // Initialize to default
    tstar = array_set64f(nobs, defaultTstar, &ierr);
    // Get number of lines in text file
    tsf = fopen(tstarTable, "r");
    nlines = 0;
    while (fgets(cline, 64, tsf) != NULL){nlines = nlines + 1;}  
    rewind(tsf);
    // Loop through and match observations to tstar's in table
    for (k=0; k<nobs; k++)
    {
        memset(search, 0, 64*sizeof(char));
        sprintf(search, "%s.%s.%s.%s",
                data[k].header.knetwk, data[k].header.kstnm,
                data[k].header.kcmpnm, data[k].header.khole);
        lfound = false;
        for (i=0; i<nlines; i++)
        {
            memset(cline, 0, 64*sizeof(char));
            memset(item,  0, 64*sizeof(char));
            fgets(cline, 64, tsf);
            sscanf(cline, "%s %lf\n", item, &tin);
            if (strcasecmp(search, item) == 0)
            {
                tstar[k] = defaultTstar; 
                lfound = true;
                break;
            }
        }
        if (!lfound)
        {
            fprintf(stderr, "%s: Setting %s tstar to: %f\n",
                    __func__, search, defaultTstar);
        }
        rewind(tsf);
    }
    fclose(tsf);
    *tstars = tstar;
    return 0;
}
//============================================================================//
/*!
 * @brief TODO: delete this function. 
 */
int prepmt_grnsTeleB_windowHudson96(
    const int nobs, const int ntstar, const int ndepth,
    const struct sacData_struct *data,
    struct sacData_struct *grns)
{
    double epoch, epochGrns;
    int indices[6], idep, ierr, iobs, it;
    for (iobs=0; iobs<nobs; iobs++)
    {
        // Figure out the origin time
        ierr = sacio_getEpochalStartTime(data[iobs].header, &epoch);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to get start time\n", __func__);
            break;
        }
        for (idep=0; idep<ndepth; idep++)
        {
            for (it=0; it<ntstar; it++)
            {
                ierr = prepmt_greens_getHudson96GreensFunctionsIndices(
                                        nobs, ntstar, ndepth,
                                        iobs, it, idep, indices);
                // Figure out the origin time
                ierr = sacio_getEpochalStartTime(grns[indices[0]].header,
                                                 &epochGrns);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Failed to get start time\n", __func__);
                    break;
                }
                // Fig
            }
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Processes the Green's functions.
 *
 * @param[in,out] grns    On input contains the Green's functions for all
 *                        observations, depths, and t*'s as well as the
 *                        processing chains. \n
 *                        On output contains the filtered Green's functions
 *                        for all observations, depths, and t*'s.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_grnsTeleB_processHudson96Greens(
    const int nobs, const int ntstar, const int ndepth,
    const struct prepmtCommands_struct cmds,
    struct sacData_struct *grns)
{
    struct serialCommands_struct commands;
    struct parallelCommands_struct parallelCommands;
    double *G, dt, dt0, epoch, epoch0, time;
    int *dataPtr, indices[6], 
        i, i0, i1, i2, ierr, iobs, idep, it, kndx, l, npts, npts0, nq,
        nwork, ny, ns, nsuse;
    bool lnewDt, lnewStartTime;
    const int nTimeVars = 12;
    const enum sacHeader_enum timeVars[12]
       = {SAC_FLOAT_A, SAC_FLOAT_O, 
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    const int spaceInquiry =-1;
    // Loop on the observations
    ns = ndepth*ntstar*6;
    dataPtr = memory_calloc32i(ns);
    for (iobs=0; iobs<nobs; iobs++)
    {
        // Set the processing structure
        memset(&parallelCommands, 0,
               sizeof(struct parallelCommands_struct)); 
               memset(&commands, 0, sizeof(struct serialCommands_struct));
        kndx = prepmt_greens_getHudson96GreensFunctionIndex(G11_GRNS,
                                                      nobs, ntstar, ndepth,
                                                      iobs, 0, 0);
        // Parse the commands
        ierr = process_stringsToSerialCommandsOptions(
                                      cmds.cmds[iobs].ncmds,
                                      (const char **) cmds.cmds[iobs].cmds,
                                      &commands);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting serial command string\n",
                    __func__);
            goto ERROR;
        }
        // Determine some characteristics of the processing
        sacio_getEpochalStartTime(grns[kndx].header, &epoch0);
        sacio_getFloatHeader(SAC_FLOAT_DELTA, grns[kndx].header, &dt0);
        lnewDt = false;
        lnewStartTime = false;
        epoch = epoch0;
        dt = dt0;
        for (i=0; i<commands.ncmds; i++)
        {
            if (commands.commands[i].type == CUT_COMMAND)
            {
                i0 = commands.commands[i].cut.i0;
                epoch = epoch + (double) i0*dt;
                lnewStartTime = true;
            }
            if (commands.commands[i].type == DOWNSAMPLE_COMMAND)
            {
                nq = commands.commands[i].downsample.nq;
                dt = dt*(double) nq;
                lnewDt = true;
            }
            if (commands.commands[i].type == DECIMATE_COMMAND)
            {
                nq = commands.commands[i].decimate.nqAll;
                dt = dt*(double) nq;
                lnewDt = true;
            }
        }
        // Set the commands on the parallel processing structure
        ierr = process_setCommandOnAllParallelCommands(ns, commands,
                                                       &parallelCommands);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting the parallel commands\n",
                    __func__);
            goto ERROR;
        }
        // Get the data
        dataPtr[0] = 0;
        for (idep=0; idep<ndepth; idep++)
        {
            for (it=0; it<ntstar; it++)
            {
                kndx = prepmt_greens_getHudson96GreensFunctionIndex(G11_GRNS,
                                                          nobs, ntstar, ndepth,
                                                          iobs, 0, 0);
                ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                              grns[kndx].header, &npts);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error getting npts\n", __func__);
                    goto ERROR;
                }
                i1 = idep*6*ntstar + it*6 + 0;
                i2 = i1 + 6;
                for (i=i1; i<i2; i++)
                {
                    dataPtr[i+1] = dataPtr[i] + npts;
                }
            }
        }
        nwork = dataPtr[6*ndepth*ntstar];
        if (nwork < 1)
        {
            fprintf(stderr, "%s: Invalid workspace size: %d\n",
                    __func__, nwork);
            ierr = 1;
            goto ERROR;
        } 
        G = memory_calloc64f(nwork);
        for (idep=0; idep<ndepth; idep++)
        {
            for (it=0; it<ntstar; it++)
            {
                ierr = prepmt_greens_getHudson96GreensFunctionsIndices(
                                        nobs, ntstar, ndepth,
                                        iobs, it, idep, indices);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error getting index\n", __func__);
                    goto ERROR;
                }
                for (i=0; i<6; i++)
                {
                    i1 = idep*6*ntstar + it*6 + i;
                    cblas_dcopy(grns[indices[i]].npts, 
                                grns[indices[i]].data, 1, &G[dataPtr[i1]], 1);
                }
            }
        }
        // Set the data
        ierr =  process_setParallelCommandsData64f(ns, dataPtr,
                                                   G, &parallelCommands);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting data\n", __func__);
            goto ERROR;
        }
        // Apply the commands
        ierr = process_applyParallelCommands(&parallelCommands);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error processing data\n", __func__);
            goto ERROR;
        }
        // Get the data
        ierr = process_getParallelCommandsData64f(parallelCommands,
                                                  spaceInquiry, spaceInquiry,
                                                  &ny, &nsuse,
                                                  dataPtr, G);
        if (ny < nwork)
        {
            memory_free64f(&G);
            G = memory_calloc64f(ny);
        }
        ny = nwork;
        ierr = process_getParallelCommandsData64f(parallelCommands,
                                                  nwork, ns,
                                                  &ny, &nsuse, dataPtr, G); 
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting data\n", __func__);
            goto ERROR;
        }
        // Unpack the data
        for (idep=0; idep<ndepth; idep++)
        {
            for (it=0; it<ntstar; it++)
            {
                ierr = prepmt_greens_getHudson96GreensFunctionsIndices(
                                        nobs, ntstar, ndepth,
                                        iobs, it, idep, indices);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error getting index\n", __func__);
                    goto ERROR;
                }
                for (i=0; i<6; i++)
                {
                    i1 = idep*6*ntstar + it*6 + i;
                    sacio_getIntegerHeader(SAC_INT_NPTS,  
                                           grns[indices[i]].header, &npts0);
                    npts = dataPtr[i1+1] - dataPtr[i1];
                    // Resize event
                    if (npts != npts0)
                    {
                        sacio_freeData(&grns[indices[i]]);
                        grns[indices[i]].data = sacio_malloc64f(npts);
                        grns[indices[i]].npts = npts;
                        sacio_setIntegerHeader(SAC_INT_NPTS, npts,
                                               &grns[indices[i]].header);
                        ierr = array_copy64f_work(npts,
                                                  &G[dataPtr[i]],
                                                  grns[indices[i]].data);
                    }
                    else
                    {
                        ierr = array_copy64f_work(npts,
                                                  &G[dataPtr[i]],
                                                  grns[indices[i]].data);
                    }
                    // Update the times
                    if (lnewStartTime)
                    {
                        for (i=0; i<6; i++)
                        {
                            // Update the picks
                            for (l=0; l<nTimeVars; l++)
                            {
                                ierr = sacio_getFloatHeader(timeVars[l],
                                              grns[indices[i]].header, &time);
                                if (ierr == 0)
                                {
                                    time = time + epoch0; // Turn to real time
                                    time = time - epoch;  // Relative to new time
                                    sacio_setFloatHeader(timeVars[l], time,
                                                    &grns[indices[i]].header);
                                }
                            } // Loop on picks
                            sacio_setEpochalStartTime(epoch,
                                                      &grns[indices[i]].header);
                        } // Loop on signals
                    }
                    // Update the sampling period
                    if (lnewDt)
                    {
                        for (i=0; i<6; i++)
                        {
                            sacio_setFloatHeader(SAC_FLOAT_DELTA, dt,
                                                 &grns[indices[i]].header);
                        }
                    }
                } // Loop on Green's functions gxx, gyy, ...
            } // Loop on t*
        } // Loop on depths
        process_freeSerialCommands(&commands);
        process_freeParallelCommands(&parallelCommands);
        memory_free64f(&G);
    }
ERROR:;
    memory_free32i(&dataPtr);
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the teleseismic body wave Green's functions modeling parameters.
 *
 * @param[in] iniFile           Name of initialization file.
 * @param[in] section           Section of initalization file to read.  For
 *                              example, this may be prepmt:telePGrns.
 * @param[out] archiveFile      This is the HDF5 file with the observed waveforms.
 * @param[out] parmtDataFile    This is the name of the data archive file to
 *                              be used by parmt.
 * @param[out] luseCrust1       If true then the modeling may use crust1 
 *                              at the source and receiver (unless superseded
 *                              by a specified local model). \n
 *                              Otherwise, use a local model or the global 
 *                              ak135 earth model.
 * @param[out] crustDir         The crust1.0 directory.  This is only relevant
 *                              if luseCrust1 is true.
 * @param[out] luseSourceModel  If true then use a local source model.
 * @param[out] sourceModel      If luseSourceModel is true then this is the
 *                              filename with the local source model.
 * @param[out] luseTstarTable   If true then read the t* attenuation factors
 *                              from a table.
 * @param[out] defaultTstar     This is the default t* (e.g., 0.4).
 * @param[out] tstarTable       If luseTstarTable is true then this is the file
 *                              with the station/t* table pairings.
 * @param[out] lreprickGrns     If true then repick the Green's functions with
 *                              an STA/LTA function.  This is in the context of
 *                              severe attenutation which noticeably effect
 *                              theoretical arrival times.
 * @param[out] staWin           If repicking Green's functions then this
 *                              is the short term window (seconds).
 * @param[out] ltaWin           If repicking Green's functions then this
 *                              is the long term window (seconds).  
 * @param[out] staltaThreshPct  This is the fraction (0,1) after which an 
 *                              arrival is declared in the STA/LTA repicking.
 * @param[out] lalignXC         If true the the Green's functions will be 
 *                              re-aligned with cross-correlation.
 * @param[out] luseEnvelope     If true and lalignXC is true then the 
 *                              Green's functions will be aligned based on 
 *                              the cross-correlation of their envelopes.
 * @param[out] lnormalizeXC     If true then each Green's function is normalized
 *                              in the Green's function cross-correlation thus
 *                              Green's functions with larger amplitudes may 
 *                              bias the result.
 * @param[out] maxXCtimeLag     This is the max time lag (seconds) available to
 *                              the waveform cross-correlation alignment.
 * @param[out] ndepth           Number of depths in modeling.
 * @param[in,out] depths        On input this is a NULL pointer.
 *                              On output this is a pointer to an array 
 *                              of dimension [ndepth] with the modeling depths
 *                              (km).
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_grnsTeleB_readParameters(const char *iniFile,
                                    const char *section,
                                    char archiveFile[PATH_MAX],
                                    char parmtDataFile[PATH_MAX],
                                    bool *luseCrust1,
                                    char crustDir[PATH_MAX],
                                    bool *luseSourceModel,
                                    char sourceModel[PATH_MAX],
                                    bool *luseTstarTable,
                                    double *defaultTstar,
                                    char tstarTable[PATH_MAX],
                                    bool *lrepickGrns,
                                    double *staWin, double *ltaWin,
                                    double *staltaThreshPct,
                                    bool *lalignXC,
                                    bool *luseEnvelope,
                                    bool *lnormalizeXC,
                                    double *maxXCtimeLag,
                                    int *ndepth, double **depths)
{
    const char *s;
    char vname[256];
    dictionary *ini;
    char *dirName;
    double *deps, dmin, dmax;
    int ierr;

    ierr = 0;
    deps = NULL;
    memset(archiveFile, 0, PATH_MAX*sizeof(char));
    memset(crustDir, 0, PATH_MAX*sizeof(char));
    memset(sourceModel, 0, PATH_MAX*sizeof(char));
    memset(tstarTable, 0, PATH_MAX*sizeof(char));
    memset(parmtDataFile, 0, PATH_MAX*sizeof(char));
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file does not exist\n", __func__);
        return -1;
    }
    ini = iniparser_load(iniFile);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dataH5File", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (!os_path_isfile(s))
    {
        fprintf(stderr, "%s: Data file %s does not exist\n", __func__, s);
        ierr = 1;
        goto ERROR;
    }
    strcpy(archiveFile, s);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:parmtDataFile", section);
    s = iniparser_getstring(ini, vname, "bodyWave.h5"); 
    strcpy(parmtDataFile, s);

    dirName = os_dirname(parmtDataFile, &ierr);
    if (!os_path_isdir(dirName))
    {
        ierr = os_makedirs(dirName);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to make directory: %s\n",
                    __func__, dirName);
            goto ERROR;
        }
    }
    memory_free8c(&dirName);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:lrepickGrns", section);
    *lrepickGrns = iniparser_getboolean(ini, vname, false);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:staWin", section);
    *staWin = iniparser_getdouble(ini, vname, 0.2);
    if (*staWin <= 0.0)
    {
        fprintf(stderr, "%s: Invalid STA length %f\n", __func__, *staWin);
        ierr = 1;
        goto ERROR;
    }

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:ltaWin", section);
    *ltaWin = iniparser_getdouble(ini, vname, 1.2);
    if (*ltaWin <= *staWin)
    {
        fprintf(stderr, "%s: Invalid LTA/STA lengths %f %f\n",
                __func__, *staWin, *ltaWin);
        ierr = 1;
        goto ERROR;
    }

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:staltaThreshPct", section);
    *staltaThreshPct = iniparser_getdouble(ini, vname, 0.8);
    if (*staltaThreshPct <= 0.0 || *staltaThreshPct > 1.0)
    {
        fprintf(stderr, "%s: Invalid STA/LTA thresh pct %f\n",
                 __func__, *staltaThreshPct);
        ierr = 1;
        goto ERROR;
    }

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:luseTstarTable", section);
    *luseTstarTable = iniparser_getboolean(ini, vname, false);
    if (*luseTstarTable)
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:tstarTable", section);
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            fprintf(stderr, "%s: tstar table %s doesn't exist\n",
                    __func__, s);
            ierr = 1;
            goto ERROR;
        }
        strcpy(tstarTable, s);
    }

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:ndepths\0");
    *ndepth = iniparser_getint(ini, vname, 0);
    if (*ndepth < 1)
    {
        fprintf(stderr, "%s: Inadequate number of depths %d\n",
                __func__, *ndepth);
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
        fprintf(stderr, "%s: Invalid dmin/dmax %f %f\n", __func__, dmin, dmax);
        ierr = 1;
        goto ERROR;
    }
    deps = array_linspace64f(dmin, dmax, *ndepth, &ierr);

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:luseCrust1\0");
    *luseCrust1 = iniparser_getboolean(ini, vname, true);
    if (*luseCrust1)
    {
        memset(vname, 0, 256*sizeof(char));
        strcpy(vname, "precompute:crustDir\0");
        s = iniparser_getstring(ini, vname, CPS_DEFAULT_CRUST1_DIRECTORY);
        if (!os_path_isdir(s))
        {
            fprintf(stderr, "%s: crust1.0 directory %s doesn't exist\n",
                    __func__, s);
            ierr = 1;
            goto ERROR;
        }
        strcpy(crustDir, s);
    }

    memset(vname, 0, 256*sizeof(char));
    strcpy(vname, "precompute:luseSourceModel\0");
    *luseSourceModel = iniparser_getboolean(ini, vname, false);
    if (*luseSourceModel)
    {
        memset(vname, 0, 256*sizeof(char));
        strcpy(vname, "precompute:sourceModel\0");
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            fprintf(stderr, "%s: Source model %s does not exist\n",
                    __func__, s);
            ierr = 1;
            goto ERROR;
        }
        strcpy(sourceModel, s);
    }
       
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:defaultTstar", section); 
    *defaultTstar = iniparser_getdouble(ini, vname, 0.4);

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
    return ierr;
}
//============================================================================//
/*!
 * @brief Utility function for reading the data for which the Green's functions
 *        will be generated.
 *
 * @param[in] archiveFile  Name of HDF5 archive file containing the
 *                         observations.
 *
 * @param[out] nobs        Number of observations.
 * @param[out] ierr        0 indicates success.
 *
 * @result The SAC data from the H5 archive for which the Green's functions
 *         will be created.  This is an array of dimension [nobs].
 *
 * @author Ben Baker, ISTI
 *
 */
struct sacData_struct *
    prepmt_grnsTeleB_readData(const char *archiveFile, int *nobs, int *ierr)
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
        fprintf(stderr, "%s: Warning - subset of data was read\n", __func__);
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
    if (!lsection){strcpy(section, "prepmt:telePGrns\0");}
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
