#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <iniparser.h>
#include "prepmt/prepmt_commands.h"
#include "sacio.h"
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/string/string.h"

int prepmt_commands_freePrepmtCommands(struct prepmtCommands_struct *cmds)
{
    int i, k;
    for (k=0; k<cmds->nobs; k++)
    {
        if (cmds->cmds[k].cmds != NULL)
        {
            for (i=0; i<cmds->cmds[k].ncmds; i++)
            {
                if (cmds->cmds[k].cmds[i] != NULL)
                {
                    free(cmds->cmds[k].cmds[i]);
                }
            }
            free(cmds->cmds[k].cmds);
        }
        memset(&cmds->cmds[k], 0, sizeof(struct prepmtCommandsChars_struct));
    }
    memset(cmds, 0, sizeof(struct prepmtCommands_struct));
    return 0;
}

struct prepmtCommands_struct
    prepmt_commands_readFromIniFile(const char *iniFile,
                                     const char *section,
                                     const int nobs,
                                     const struct sacData_struct *data,
                                     int *ierr) 
{
    const char *fcnm = "prepmt_commands_readFromIniFile\0";
    FILE *pfl;
    const char *s;
    char **cwork, vname[256], cline[MAX_CMD_LEN];
    char cdum[64], chan[64], loc[64], netw[64], stat[64],
         kcmpnm[64], khole[64], knetwk[64], kstnm[64];
    struct prepmtCommands_struct cmds;
    int *indx, i, icmd, j, k, ncmds, ncmdsWork, nindex, nlines;
    dictionary *ini;
    size_t lenos;
    bool lfound, luseProcessingList;
    //------------------------------------------------------------------------//
    *ierr = 0;
    memset(&cmds, 0, sizeof(struct prepmtCommands_struct));
    if (!os_path_isfile(iniFile))
    {
        printf("%s: Error - ini file %s does not exist\n", fcnm, iniFile);
        *ierr = 1;
        return cmds;
    }
    ini = iniparser_load(iniFile);
 
    cmds.nobs = MAX(1, nobs);
    cmds.cmds = (struct prepmtCommandsChars_struct *)
                calloc((size_t) cmds.nobs,
                       sizeof(struct prepmtCommandsChars_struct));

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:useProcessingList", section);
    luseProcessingList = iniparser_getboolean(ini, vname, false); 
    if (!luseProcessingList)
    {
        cmds.cmds[0].lgeneric = true;
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:ncmds", section);
        ncmds = iniparser_getint(ini, vname, 0);
        if (ncmds > 0)
        {
            ncmdsWork = ncmds;
            ncmds = 0;
            cwork = (char **) calloc((size_t) ncmdsWork, sizeof(char *));
            for (k=0; k<ncmdsWork; k++)
            {
                memset(vname, 0, 256*sizeof(char));
                sprintf(vname, "%s:command_%d", section, k+1);
                s = iniparser_getstring(ini, vname, NULL);
                if (s == NULL){continue;}
                lenos = strlen(s);
                cwork[ncmds] = (char *) calloc(lenos+1, sizeof(char));
                strcpy(cwork[ncmds], s);
                //cwork[ncmds] = strdup(s); 
                ncmds = ncmds + 1;
            }
            cmds.cmds[0].lgeneric = true;
            cmds.cmds[0].ncmds = ncmdsWork;
            cmds.cmds[0].cmds = cwork;
        }
        // Copy the results to all other observations 
        for (k=1; k<nobs; k++)
        {
            cmds.cmds[k].lgeneric = cmds.cmds[0].lgeneric;
            cmds.cmds[k].ncmds = cmds.cmds[0].ncmds;
            if (cmds.cmds[k].ncmds == 0){continue;}
            cmds.cmds[k].cmds =
                  (char **) calloc((size_t) cmds.cmds[k].ncmds, sizeof(char *));
            for (i=0; i<cmds.cmds[0].ncmds; i++)
            {
                lenos = strlen(cmds.cmds[0].cmds[i]);
                cmds.cmds[k].cmds[i] = (char *) calloc(lenos, sizeof(char));
                strcpy(cmds.cmds[k].cmds[i], cmds.cmds[0].cmds[i]);
                //cmds.cmds[k].cmds[i] = strdup(cmds.cmds[0].cmds[i]);
            }
        }
    }
    // Read a processing list
    else
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:processingList", section);
        if (!os_path_isfile(s))
        {
            printf("%s: Processing list file %s does not exist\n", fcnm, s);
            *ierr = 1;
            return cmds;
        }
        pfl = fopen(s, "r");
        // Count the number of lines
        nlines = 0;
        while (fgets(cline, MAX_CMD_LEN, pfl) != NULL){nlines = nlines + 1;}
        // Index the file
        rewind(pfl);
        nindex = 0;
        indx = (int *) calloc((size_t) nlines+1, sizeof(int));
        for (i=0; i<nlines; i++)
        {
            memset(cline, 0, MAX_CMD_LEN*sizeof(char));
            fgets(cline, MAX_CMD_LEN, pfl);
            if (strcasecmp("SNCL:", cline) >= 0)
            {
                nindex = nindex + 1;
                indx[nindex] = i;
            }
        }
        // Hunt through file
        rewind(pfl);
        for (k=0; k<nobs; k++)
        {
            lfound = false; 
            sacio_getCharacterHeader(SAC_CHAR_KNETWK, data[k].header, knetwk);
            sacio_getCharacterHeader(SAC_CHAR_KSTNM,  data[k].header, kstnm);
            sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data[k].header, kcmpnm);
            sacio_getCharacterHeader(SAC_CHAR_KHOLE,  data[k].header, khole); 
            for (i=0; i<nindex; i++)
            {
                for (j=indx[i]; j<indx[i+1]; j++)
                {
                    memset(cline, 0, MAX_CMD_LEN*sizeof(char));
                    fgets(cline, MAX_CMD_LEN, pfl);
                    if (j == indx[i])
                    {
                        memset(cdum, 0, 64*sizeof(char));
                        memset(netw, 0, 64*sizeof(char));
                        memset(stat, 0, 64*sizeof(char));
                        memset(chan, 0, 64*sizeof(char));
                        memset(loc,  0, 64*sizeof(char)); 
                        sscanf(cline, "%s:%s.%s.%s.%s\n",
                                      cdum, netw, stat, chan, loc); 
                        if (strcasecmp(knetwk, netw) == 0 &&
                            strcasecmp(kstnm,  stat) == 0 &&
                            strcasecmp(kcmpnm, chan) == 0)
                        {
                            if (strcasecmp(khole, loc) != 0)
                            {
                                printf("%s: Location code mismatch\n", fcnm);
                            }
                            ncmds = indx[i+1] - indx[i];
                            cmds.cmds[k].lgeneric = false;
                            cmds.cmds[k].cmds
                             = (char **) calloc((size_t) ncmds, sizeof(char *));
                            lfound = true;
                        }
                    }
                    else
                    {
                        if (lfound)
                        {
                            icmd = j - (indx[i] + 1); // + 1 removes header line
                            lenos = strlen(cline);
                            cmds.cmds[k].cmds[icmd]
                               = (char *) calloc(lenos+1, sizeof(char));
                            strcpy(cmds.cmds[k].cmds[icmd], cline);
                        }
                    }
                }
                if (lfound){break;}
            }
            if (!lfound)
            {
                printf("%s: Warning - no commands for: %s.%s.%s.%s\n",
                       fcnm, knetwk, kstnm, kcmpnm, khole);
            }
            rewind(pfl);
        } // Loop on stations
        fclose(pfl);
        free(indx);
    }
    iniparser_freedict(ini);
    return cmds;
}

int prepmt_commands_modifyCommandsCharsStruct(
    const struct prepmtModifyCommands_struct options,
    const struct sacData_struct data,
    struct prepmtCommandsChars_struct *cmds)
{
    const char *fcnm = "prepmt_commands_modifyCommandsCharsStruct\0";
    char **cwork;
    int i, ierr;
    size_t lenos;
    if (cmds->ncmds < 1){return 0;}
    if (cmds->cmds == NULL)
    {
        printf("%s: Commands are NULL\n", fcnm);
        return -1;
    }
    if (!cmds->lgeneric){return 0;}
    cwork = prepmt_commands_modifyCommands(cmds->ncmds,
                                           (const char **) cmds->cmds,
                                           options, data, &ierr); 
    if (ierr != 0)
    {
        printf("%s: Failed to modify commands\n", fcnm);
        return -1;
    }
    for (i=0; i<cmds->ncmds; i++)
    {
        free(cmds->cmds[i]);
        if (cwork[i] != NULL)
        {
            lenos = strlen(cwork[i]);
            cmds->cmds[i] = (char *) calloc(lenos+1, sizeof(char));
            strcpy(cmds->cmds[i], cwork[i]);
            free(cwork[i]);
        }
        else
        {
            cmds->cmds[i] = (char *) calloc(1, sizeof(char));
        }
    }
    free(cwork);
    return 0;
}
                                    
char **prepmt_commands_modifyCommands(
    const int ncmds, const char **cmds,
    const struct prepmtModifyCommands_struct options,
    const struct sacData_struct data,
    int *ierr)
{
    const char *fcnm = "prepmt_commands_modifyCommands\0";
    char **newCmds, **cmdSplit, cwork[MAX_CMD_LEN], c64[64], cmd1[64], cmd2[64];
    double *freqs, cut0, cut1, dt0, epoch, gainFix, ptime, t0, t1, targetDt;
    struct signalZPK_struct zpkFrom, zpkTo;
    int i, iodva, l, npolesAdd, nsplit, nzerosAdd;
    size_t lenos;
    bool ldeconvolution, oneCorner;
    const bool lisDigital = false; // SAC-PZ's aren't digital
    const bool laaFilter = true; // anti-alias filter in decimation
    const bool lfixPhase = true; // don't let anti-alias filter mess up picks
    *ierr = 0;
    freqs = NULL;
    cut0 = options.cut0;
    cut1 = options.cut1;
    targetDt = options.targetDt;
    ldeconvolution = options.ldeconvolution;
    iodva = options.iodva;
    if (iodva < 0 || iodva > 3)
    {
        printf("%s: No idea what iodva=%d means\n", fcnm, iodva);
    }
    newCmds = NULL;
    if (ncmds < 1){return 0;} // Nothing to do
    // Modify the problematic commands
    *ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA, data.header, &dt0);
    newCmds = (char **) calloc((size_t) ncmds, sizeof(char *));
    for (i=0; i<ncmds; i++)
    {
        newCmds[i] = (char *) calloc(MAX_CMD_LEN, sizeof(char));
        lenos = strlen(cmds[i]);
        memset(cwork, 0, MAX_CMD_LEN*sizeof(char));
        memset(cmd1, 0, 64*sizeof(char));
        memset(cmd2, 0, 64*sizeof(char));
        strncpy(cmd1, cmds[i], MIN(3, lenos));
        strncpy(cmd2, cmds[i], MIN(3, lenos));
        npolesAdd = 0;
        nzerosAdd = 0;
        if (strcasecmp(cmds[i], "transfer\0") == 0)
        {
            memset(&zpkTo, 0, sizeof(struct signalZPK_struct));
            memset(&zpkFrom, 0, sizeof(struct signalZPK_struct));
            gainFix = 1.0;
            //printf("%d\n", (int) data.pz.lhavePZ);
            if (data.pz.lhavePZ)
            {
                // 1/From(w) -> need to eliminate a zero from To(w)
                // by adding a zero
                if (data.pz.inputUnits == SAC_METERS ||
                    data.pz.inputUnits == SAC_NANOMETERS)
                {
                    // Displacement to displacement
                    if (iodva == 0)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                    // Velocity to displacement
                    else if (iodva == 1)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 1;
                    }
                    // Acceleration to displacement
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 2;
                    }
                    // Make input unit meters
                    if (data.pz.inputUnits == SAC_NANOMETERS)
                    {
                        gainFix = 1.e-9;
                    }
                }
                // 1/From(w) = 1 -> do nothing
                else if (data.pz.inputUnits == SAC_METERS_SECOND ||
                         data.pz.inputUnits == SAC_NANOMETERS_SECOND)
                {
                    // Displacement to velocity
                    if (iodva == 0)
                    {
                        npolesAdd = 1;
                        nzerosAdd = 0;
                    }
                    // Velocity to velocity
                    else if (iodva == 1)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                    // Acceleration to velocity
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 1;
                    }
                    // Make input unit meters/second
                    if (data.pz.inputUnits == SAC_NANOMETERS_SECOND)
                    {   
                        gainFix = 1.e-9;
                    }
                }
                // 1/From(w) -> add a zero to To(w) by adding a pole
                else if (data.pz.inputUnits == SAC_METERS_SECOND_SECOND ||
                         data.pz.inputUnits == SAC_NANOMETERS_SECOND_SECOND)
                {
                    // Displacement to acceleration
                    if (iodva == 0)
                    {
                        npolesAdd = 2;
                        nzerosAdd = 0;
                    }
                    // Velocity to acceleration
                    else if (iodva == 1)
                    {
                        npolesAdd = 1;
                        nzerosAdd = 0;
                    }
                    // Acceleration to acceleration
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                    // make input unit meters/second**2
                    if (data.pz.inputUnits == SAC_NANOMETERS_SECOND_SECOND)
                    {
                        gainFix = 1.e-9;
                    }
                }
                else if (data.pz.inputUnits == SAC_UNKNOWN_UNITS)
                {
                    printf("%s: Unknown input units\n", fcnm);
                    *ierr = 1;
                    break;
                }
            }
            // Set the poles and zeros
            zpkFrom.npoles = npolesAdd;
            zpkFrom.nzeros = nzerosAdd; 
            zpkFrom.k = 1.0;
            if (zpkFrom.npoles > 0)
            {
                zpkFrom.p = memory_calloc64z(zpkFrom.npoles);
            }
            if (zpkFrom.nzeros > 0)
            {
                zpkFrom.z = memory_calloc64z(zpkFrom.nzeros);
            }
            zpkTo.npoles = data.pz.npoles;
            zpkTo.nzeros = data.pz.nzeros;
            zpkTo.k = data.pz.constant*gainFix;
            if (zpkTo.npoles > 0)
            {
                zpkTo.p = array_copy64z(zpkTo.npoles, data.pz.poles, ierr);
            }
            if (zpkTo.nzeros > 0)
            {
                zpkTo.z = array_copy64z(zpkTo.nzeros, data.pz.zeros, ierr);
            }
            if (!ldeconvolution)
            {
                *ierr = transfer_poleZeroInstrumentResponsesToString(dt0,
                                                             freqs,
                                                             lisDigital,
                                                             true, zpkFrom,
                                                             true, zpkTo,
                                                             cwork);
            }
            else
            {
                *ierr = transfer_poleZeroInstrumentResponsesToString(dt0,
                                                             freqs,
                                                             lisDigital,
                                                             true, zpkTo,
                                                             true, zpkFrom,
                                                             cwork);
            }
            memory_free64z(&zpkTo.p);
            memory_free64z(&zpkTo.z);
            memory_free64z(&zpkFrom.p);
            memory_free64z(&zpkFrom.z);
            if (*ierr != 0)
            {
                log_errorF("%s: Error setting transfer command\n", fcnm);
                goto ERROR;
            }
/*
transfer_poleZeroInstrumentResponsesToString(
    dt0,
    const double *freqs,
    const bool lisDigital,
    const bool lhaveFrom, const struct signalZPK_struct zpkFrom,
    const bool lhaveTo, const  struct signalZPK_struct zpkTo,
    cwork);
*/
//printf("%f\n", gainFix);
//printf("%s\n", cwork);
//getchar();
        }
        else if (strcasecmp(cmd1, "lowpass\0") == 0 ||
                 strcasecmp(cmd1, "highpas\0") == 0)
        {
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " "); 
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corner\0") == 0)
                {
                    memset(c64, 0, 64*sizeof(char));
                    sprintf(c64, "dt %f corner %s", dt0, cmdSplit[l+1]);
                    strcat(cwork, c64);
                    l = l + 1;
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmd1, "bandpas\0") == 0 ||
                 strcasecmp(cmd1, "bandrej\0") == 0)
        {
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " ");
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corners\0") == 0)
                {
                    memset(c64, 0, 64*sizeof(char));
                    sprintf(c64, "dt %f corners %s %s",
                            dt0, cmdSplit[l+1], cmdSplit[l+2]);
                    strcat(cwork, c64);
                    l = l + 2;
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmd2, "sos\0") == 0)
        {
            oneCorner = false;
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " ");
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "lowpass\0") == 0 ||
                    strcasecmp(cmdSplit[l], "highpas\0") == 0)
                {
                    oneCorner = true;
                }
                if (strcasecmp(cmdSplit[l], "bandpas\0") == 0 ||
                    strcasecmp(cmdSplit[l], "bandrej\0") == 0)
                {
                    oneCorner = false;
                }
            }
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corner\0") == 0 ||
                    strcasecmp(cmdSplit[l], "corners\0") == 0)
                {
                    if (oneCorner)
                    {
                        memset(c64, 0, 64*sizeof(char));
                        sprintf(c64, "dt %f corner %s", dt0, cmdSplit[l+1]);
                        strcat(cwork, c64);
                        l = l + 1;
                    }
                    else
                    {
                        memset(c64, 0, 64*sizeof(char));
                        sprintf(c64, "dt %f corners %s %s",
                                dt0, cmdSplit[l+1], cmdSplit[l+2]);
                        strcat(cwork, c64);
                        l = l + 2;
                    }
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmds[i], "cut\0") == 0)
        {
            *ierr = sacio_getEpochalStartTime(data.header, &epoch);
            if (*ierr != 0)
            {
                log_errorF("%s: Failed to get start time\n", fcnm);
                goto ERROR;
            }
            *ierr = sacio_getFloatHeader(SAC_FLOAT_A, data.header,
                                         &ptime);
            if (*ierr != 0)
            {
                log_errorF("%s: Failed to get pick time\n", fcnm);
                goto ERROR;
            }
            t0 = epoch + ptime + cut0; // superfluous; epoch will be removed
            t1 = epoch + ptime + cut1; // superfluous; epoch will be removed
            *ierr = cut_cutEpochalTimesToString(dt0, epoch, t0, t1, cwork);
            if (*ierr != 0)
            {
                log_errorF("%s: Failed to modify cut command\n", fcnm);
                goto ERROR;
            }
        }
        else if (strcasecmp(cmds[i], "decimate\0") == 0)
        {
            *ierr = decimate_createDesignCommandsWithDT(dt0, targetDt,
                                                        laaFilter, lfixPhase,
                                                        cwork);
            if (*ierr != 0)
            {
                log_errorF("%s: Couldn't modify the decimate command\n",
                           fcnm);
                goto ERROR;
            }
            dt0 = targetDt;
        }
        else if (strcasecmp(cmds[i], "downsample\0") == 0)
        {
            *ierr = downsample_downsampleTargetSamplingPeriodToString(
                        dt0, targetDt, cwork);
            if (*ierr != 0)
            {
                log_errorF("%s: Couldn't modify downsample command\n",
                           fcnm);
                goto ERROR;
            }
            dt0 = targetDt;
        }
        else
        {
            strcpy(cwork, cmds[i]);
        }
        // Update the command
        strcpy(newCmds[i], cwork);
    } // Loop on commands
ERROR:;
    return newCmds;
}
