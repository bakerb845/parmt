#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "prepmt/prepmt_pickFile.h"
#include "sacio.h"
#include "iscl/os/os.h"
#include "iscl/time/time.h"

/*!
 * @brief Reads a NonLinLoc hyp pick file and populates the SAC headers.
 *
 * @param[in] pickFile     Name of file with picks.
 * @param[in] pickVar      SAC float header variable to which the pick will be
 *                         written.  For example, this could be SAC_FLOAT_A.
 * @param[in] pickNameVar  SAC character header variable to which the phase
 *                         name, and possibly, polarity will be written. 
 *                         For example, this could be SAC_CHAR_KA.
 * @param[in] nobs         Number of waveforms.
 * @param[out] data        Contains the first arrival picks corresponding to
 *                         the picks in the pick file.  The waveforms are
 *                         matched by SNCL.  This is an array of dimension
 *                         [nobs].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_pickFile_nonLinLoc2sac(const char *pickFile,
                                  const enum sacHeader_enum pickVar, 
                                  const enum sacHeader_enum pickNameVar,
                                  const int nobs, struct sacData_struct *data)
{
    const char *fcnm = "prepmt_pickFile_nonLinLoc2sac\0";
    double *epochs, epoch, errMag, codaDuration, amplitude, period,
           pick, weight, second;
    char **netStat, **chans, **locs, **phases, cline[256];
    char knetst[64], knetw[64], kstnm[64], kcmpnm[64], khole[64];
    char phase[8], fm[2], pid[2], cdum[4];
    FILE *hfl;
    long ymd, hm, work;
    int dom, i, ierr, k, month, nlines, nzhour, nzmin, nzmusec, nzsec, nzyear;
    bool lfound;
    ierr = 0;
    if (!os_path_isfile(pickFile))
    {
        printf("%s: Error pick file %s does not exist\n", fcnm, pickFile);
        return -1;
    }
    // Read the pick file
    hfl = fopen(pickFile, "r");
    nlines = 0;
    while (fgets(cline, 256, hfl) != NULL){nlines = nlines + 1;} 
    if (nlines < 1)
    {
        printf("%s: No picks in %s\n", fcnm, pickFile); 
        fclose(hfl);
        return 0;
    }
    rewind(hfl);
    netStat = (char **) calloc((size_t) nlines, sizeof(char *));
    chans  = (char **) calloc((size_t) nlines, sizeof(char *));
    locs   = (char **) calloc((size_t) nlines, sizeof(char *));
    phases = (char **) calloc((size_t) nlines, sizeof(char *));
    epochs = (double *) calloc((size_t) nlines, sizeof(double));
    for (i=0; i<nlines; i++)
    {
        memset(cline, 0, 256*sizeof(char));
        fgets(cline, 256, hfl);
        netStat[i] = (char *) calloc(16, sizeof(char));
        chans[i] = (char *) calloc(8, sizeof(char));
        locs[i]  = (char *) calloc(8, sizeof(char));
        memset(pid, 0, 2*sizeof(char));
        memset(phase, 0, 8*sizeof(char));
        memset(fm, 0, 2*sizeof(char));
        // Here i'm assuming the network code is 2 characters but allowing
        // the station lenght to be somewhat longer
        sscanf(cline,
              "%s %4s %4s %1s %6s %1s %ld %ld %lf %s %lf %lf %lf %lf %lf\n",
               netStat[i], locs[i], chans[i],
               pid, phase, fm, 
               &ymd, &hm, &second,
               cdum,
               &errMag, &codaDuration, &amplitude, &period, &weight);
        // Convert YYYYMMDD to year, day of month, hour
        nzyear = (int) ((double) ymd*1.e-4);
        work   = ymd - (long) nzyear*10000;
        month  = (int) ((double) work*1.e-2);
        work   = ymd - (long) nzyear*10000 - (long) month*100;
        dom    = (int) ((double) work);

        nzhour = (int) ((double) hm*1.e-2);
        work   = hm - ((long) nzhour*100); 
        nzmin  = (int) ((double) work);

        nzsec = (int) second;
        nzmusec = (int) ((second - (double) nzsec)*1.e6);
        //printf("%d %d %d %d %d %d %d\n", nzyear, month,
        //       dom, nzhour, nzhour, nzmin, nzmusec);
        phases[i] = (char *) calloc(8, sizeof(char));
        strcpy(phases[i], phase);
        if (fm[0] != '?')
        {
            if (fm[0] == '+'){strcat(phases[i], "+");}
            if (fm[0] == '-'){strcat(phases[i], "-");}
        }
        epochs[i] = time_calendar2epoch2(nzyear, month, dom,
                                         nzhour, nzmin, nzsec,
                                         nzmusec);
    }
    fclose(hfl);
    // Match the picks
    for (k=0; k<nobs; k++)
    {
        memset(knetst, 0, 64*sizeof(char));
        memset(knetw,  0, 64*sizeof(char));
        memset(kstnm,  0, 64*sizeof(char));
        memset(kcmpnm, 0, 64*sizeof(char));
        memset(khole,  0, 64*sizeof(char));
        ierr  = sacio_getCharacterHeader(SAC_CHAR_KNETWK,
                                         data[k].header, knetw);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KSTNM,
                                         data[k].header, kstnm);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM,
                                         data[k].header, kcmpnm);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KHOLE,
                                         data[k].header, khole);
        if (ierr != 0)
        {
            printf("%s: Couln't get SNCL\n", fcnm);
            goto ERROR;
        }
        sprintf(knetst, "%-s%-s", knetw, kstnm);
        // Figure out pick epochal time
        ierr = sacio_getEpochalStartTime(data[k].header, &epoch);
        if (ierr != 0)
        {
            printf("%s: Couldn't get start time of trace\n", fcnm);
            goto ERROR; 
        }
        lfound = false;
        for (i=0; i<nlines; i++)
        {
              if (strcasecmp(netStat[i], knetst) == 0 &&
                  strcasecmp(chans[i], kcmpnm) == 0 &&
                  strcasecmp(locs[i],  khole) == 0)
            {

                pick = epochs[i] - epoch;
                //sacio_setFloatHeader(SAC_FLOAT_A, pick, &data[k].header);
                //sacio_setCharacterHeader(SAC_CHAR_KA, phases[i],
                //                         &data[k].header);
                sacio_setFloatHeader(pickVar, pick, &data[k].header);
                sacio_setCharacterHeader(pickNameVar, phases[i], 
                                         &data[k].header);
                lfound = true;
                break;
            }
        }
        if (!lfound)
        {
            printf("%s: Couldn't find pick for %s.%s.%s\n",
                   fcnm, netStat[i], chans[i], locs[i]);
        }
    }

ERROR:;
    // Release memory
    for (i=0; i<nlines; i++)
    {
        free(netStat[i]);
        free(chans[i]);
        free(locs[i]);
        free(phases[i]);
    }
    free(netStat);
    free(chans);
    free(locs);
    free(phases);
    free(epochs);
    return ierr;
}
//============================================================================//
/*!
 * @brief Writes a line which can be written to a NonLinLoc *hyp pick file.
 *
 * @param[in] data         SAC data from which to extract pick information.
 * @param[in] pickVar      SAC float header variable to extract the pick time.
 *                         For example this could be SAC_FLOAT_A.  This
 *                         must be defined.
 * @param[in] pickNameVar  SAC float header variables to extract the pick
 *                         phase.  For example this could be SAC_CHAR_KA.
 *                         This must be specified.  If a polarity is associated
 *                         with the phase then this should end in + or - for
 *                         positive or negative polarity, respectively.
 * @param[in] errMagVar    SAC float header variable containing the 
 *                         magnitude error in seconds.  For example
 *                         this could be SAC_FLOAT_USER0.  This variable
 *                         does not need to be defined in the SAC header.
 * @param[in] codaDurVar   SAC float header variable containing the 
 *                         coda duration.  For example, this could be
 *                         SAC_FLOAT_USER1.  This variable does not need
 *                         to be defined in the SAC header.
 * @param[in] ampVar       SAC float header variable with the maximum 
 *                         peak-to-peak amplitude. For example, this could be
 *                         SAC_FLOAT_USER2.  This variable does not need
 *                         to be defined in the SAC header.
 * @param[in] periodVar    SAC float header variable with the period of the
 *                         amplitude reading.   For example, this could be
 *                         SAC_FLOAT_USER3.  This variable does not need to
 *                         be defined.
 * @param[in] weightVar    SAC float header with the prior weight for this
 *                         observation.  This variable does not need to be
 *                         defined.
 *
 * @param[out] line        On successful exit this line can be safely written
 *                         to a hyp file.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_pickFile_sacToNonLinLocLine(const struct sacData_struct data,
                                       const enum sacHeader_enum pickVar,
                                       const enum sacHeader_enum pickNameVar,
                                       const enum sacHeader_enum errMagVar,
                                       const enum sacHeader_enum codaDurVar,
                                       const enum sacHeader_enum ampVar,
                                       const enum sacHeader_enum periodVar,
                                       const enum sacHeader_enum weightVar,
                                       char line[256])
{
    const char *fcnm = "prepmt_pickFile_sacToNonLinLocLine\0";
    char phase[8], firstMotion[1], pid[1];
    int ierr, nzyear, nzjday, month, mday, nzhour, nzmin, nzsec, nzmusec;
    double epoch, pick, second, temp;
    size_t lenos;
    double errMag = 0.0;
    double codaDuration = 0.0;
    double amplitude = 0.0;
    double period = 0.0;
    double weight = 1.0;
    // Null out result
    memset(line, 0, 256*sizeof(char));
    memset(phase, 0, 8*sizeof(char));
    // Figure out pick epochal time
    ierr = sacio_getEpochalStartTime(data.header, &epoch);
    if (ierr != 0)
    {
        printf("%s: Couldn't get start time of trace\n", fcnm);
        return -1;
    }
    ierr = sacio_getFloatHeader(pickVar, data.header, &pick);
    if (ierr != 0)
    {
        printf("%s: Couldn't get pick time\n", fcnm);
        return -1;
    }
    pick = pick + epoch;
    lenos = strlen(phase); 
    // Convert the pick to calendrical time
    ierr = time_epoch2calendar(epoch,
                               &nzyear, &nzjday, &month, &mday,
                               &nzhour, &nzmin, &nzsec, &nzmusec);
    second = (double) nzsec + (double) nzmusec*1.e-6;
    // Figure out phase
    ierr = sacio_getCharacterHeader(pickNameVar, data.header, phase);
    if (ierr != 0)
    {
        printf("%s: Error - could not get phase\n", fcnm);
        return -1;
    }
    if (phase[0] != 'P' || phase[0] != 'p' ||
        phase[0] != 'S' || phase[0] != 's')
    {
        printf("%s: Warning - I don't know what type of phase this is\n", fcnm);
    }
    pid[0] = '?'; //phase[0];
    // Figure out first motion and remove it from string
    firstMotion[0] = '?';
    if (lenos > 0)
    {
        if (phase[lenos-1] == '+')
        {
            firstMotion[0] = '+';
            phase[lenos-1] = '\0';
        }
        else if (phase[lenos-1] == '-')
        {
            firstMotion[0] = '-';
            phase[lenos-1] = '\0';
        }
    }
    // Extract some extra info from header if it is defined
    errMag = 0.0;
    if (errMagVar != SAC_UNKNOWN_HDRVAR)
    {
        ierr = sacio_getFloatHeader(errMagVar, data.header, &temp);
        if (ierr != 0 && temp > 0.0){errMag = temp;}
    }

    codaDuration = 0.0;
    if (codaDurVar != SAC_UNKNOWN_HDRVAR)
    {
        ierr = sacio_getFloatHeader(codaDurVar, data.header, &temp);
        if (ierr != 0 && temp > 0.0){codaDuration = temp;}
    }

    amplitude = 0.0;
    if (ampVar != SAC_UNKNOWN_HDRVAR)
    {
        ierr = sacio_getFloatHeader(ampVar, data.header, &temp);
        if (ierr != 0 && amplitude > 0.0){amplitude = temp;}
    }

    period = 0.0;
    if (periodVar != SAC_UNKNOWN_HDRVAR)
    {
        ierr = sacio_getFloatHeader(periodVar, data.header, &temp);
        if (ierr != 0 && temp > 0.0){period = temp;}
    }

    weight = 1.0;
    if (weightVar != SAC_UNKNOWN_HDRVAR)
    {
        ierr = sacio_getFloatHeader(weightVar, data.header, &temp);
        if (ierr != 0 && temp >= 0.0){weight = temp;}
    }

    // Get the epochal time.  Note in GridLib.c Lomax reads with %s 
    sprintf(line,
            "%-4s%-6s %-4s %-4s %1s %-6s %1s %4d%02d%02d %02d%02d %7.4f GAU %9.2e %9.2e %9.2e %9.2e %9.2e",
            data.header.knetwk, data.header.kstnm,
            data.header.khole,  data.header.kcmpnm,
            pid, phase, firstMotion,
            nzyear, month, mday, nzhour, nzmin, second,
            errMag, codaDuration, amplitude, period, weight);
    return 0;
}
