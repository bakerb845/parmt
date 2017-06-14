#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "parmt_utils.h"
#include "prepmt/prepmt_dataArchive.h"
#include "prepmt/prepmt_greens.h"
#include "prepmt/prepmt_hudson96.h"
#include "prepmt/prepmt_commands.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/log/log.h"
#include "iscl/fft/fft.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/signal/convolve.h"
#include "iscl/signal/filter.h"

static int getPrimaryArrival(const struct sacHeader_struct hdr,
                             double *time, char phaseName[8]);
static int getPrimaryArrivalNoPolarity(const struct sacHeader_struct hdr,
                                       double *time, char phaseName[8]);
static void shiftGreens(const int npgrns, const int lag,
                        double *__restrict__ G, double *__restrict__ work);


//============================================================================//
/*!
 * @brief Converts the fundamental faults Green's functions to Green's
 *        functions that can be used by parmt.
 *
 * @param[in] nobs      Number of observations.
 * @param[in] obs       Observed waveforms to contextualize.  This is an 
 *                      array of length [nobs].
 * @param[in] ndepth    Number of depths.
 * @param[in] ntstar    Number of t*s'.
 * @param[in] ffGrns    Fundamental fault Green's functions for every t* and
 *                      depth in the grid search for each observation.  This
 *                      is an array of dimension [nobs x ndepth x ntstar x 10]. 
 *
 * @param[in,out] grns  Holds space for Green's functions.
 *                      Contains the Green's functions that can be applied to
 *                      a moment tensor to produce a synthetic for every t*
 *                      and depth in the grid search for each observation.
 *                      This is an array of length [nobs x ndepth x ntstar x 6].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_greens_ffGreensToGreens(const int nobs,
                                   const struct sacData_struct *obs,
                                   const int ndepth,
                                   const int ntstar,
                                   const struct sacData_struct *ffGrns,
                                   struct sacData_struct *grns)
{
    const char *fcnm = "prepmt_greens_ffGreenToGreens\0";
    char knetwk[8], kstnm[8], kcmpnm[8], khole[8], phaseName[8],
         phaseNameGrns[8];
    double az, baz, cmpaz, cmpinc, cmpincSEED, epoch, epochNew,
           evla, evlo, o, pick, pickTime, pickTimeGrns, stel, stla, stlo;
    int indices[6], i, icomp, id, ierr, idx, iobs, it, kndx, l, npts;
    const char *kcmpnms[6] = {"GXX\0", "GYY\0", "GZZ\0",
                              "GXY\0", "GXZ\0", "GYZ\0"};
/*
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
*/
    const int nTimeVars = 11; 
    const enum sacHeader_enum pickVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    //memset(grns, 0, sizeof(struct tdSearchGreens_struct));
    for (iobs=0; iobs<nobs; iobs++)
    {
        ierr = 0;
        ierr += sacio_getFloatHeader(SAC_FLOAT_AZ,
                                     obs[iobs].header, &az);
        ierr += sacio_getFloatHeader(SAC_FLOAT_BAZ,
                                     obs[iobs].header, &baz);
        ierr += sacio_getFloatHeader(SAC_FLOAT_CMPINC,
                                     obs[iobs].header, &cmpinc);
        ierr += sacio_getFloatHeader(SAC_FLOAT_CMPAZ,
                                     obs[iobs].header, &cmpaz);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA,
                                     obs[iobs].header, &evla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO,
                                     obs[iobs].header, &evlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLA,
                                     obs[iobs].header, &stla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLO,
                                     obs[iobs].header, &stlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STEL,
                                     obs[iobs].header, &stel);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KNETWK,
                                         obs[iobs].header, knetwk);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KSTNM,
                                         obs[iobs].header, kstnm);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM,
                                         obs[iobs].header, kcmpnm);
        if (ierr != 0)
        {
            log_errorF("%s: Error reading header variables\n", fcnm);
            break;
        }
        // Station location code is not terribly important
        sacio_getCharacterHeader(SAC_CHAR_KHOLE, obs[iobs].header, khole);
        cmpincSEED = cmpinc - 90.0; // SAC to SEED convention
        // Get the primary arrival
        ierr = sacio_getEpochalStartTime(obs[iobs].header, &epoch);
        if (ierr != 0)
        {
            log_errorF("%s: Error getting start time\n", fcnm);
            break;
        }
        ierr += getPrimaryArrivalNoPolarity(obs[iobs].header, &pickTime,
                                            phaseName);
        //ierr += getPrimaryArrival(obs[iobs].header, &pickTime, phaseName);
        if (ierr != 0)
        {
            log_errorF("%s: Error getting primary pick\n", fcnm);
            break;
        }
        // Need to figure out the component
        icomp = 1;
        if (kcmpnm[2] == 'Z' || kcmpnm[2] == 'z' || kcmpnm[2] == '1')
        {
            icomp = 1;
        }
        else if (kcmpnm[2] == 'N' || kcmpnm[2] == 'n' || kcmpnm[2] == '2')
        {
            icomp = 2;
        }
        else if (kcmpnm[2] == 'E' || kcmpnm[2] == 'e' || kcmpnm[2] == '3')
        {
            icomp = 3;
        }
        else
        {
            log_errorF("%s: Can't classify component: %s\n", fcnm, kcmpnm);
        }
        // Process all Green's functions in this block
        for (id=0; id<ndepth; id++)
        {
            for (it=0; it<ntstar; it++)
            {
                idx = prepmt_hudson96_observationDepthTstarToIndex(
                             ndepth, ntstar, iobs, id, it);
                kndx = 10*idx;
                sacio_getFloatHeader(SAC_FLOAT_O, ffGrns[kndx].header, &o);
                getPrimaryArrival(ffGrns[kndx].header,
                                  &pickTimeGrns, phaseNameGrns);
                if (strcasecmp(phaseNameGrns, phaseName) != 0)
                {
                    log_warnF("%s: Phase name mismatch %s %s\n",
                              fcnm, phaseName, phaseNameGrns);
                }
                npts = ffGrns[kndx].npts;
                ierr = prepmt_greens_getHudson96GreensFunctionsIndices(
                                 nobs, ntstar, ndepth,
                                 iobs, it, id,
                                 indices);
                /*
                indx = prepmt_greens_getHudson96GreensFunctionIndex(G11_GRNS,
                                                            nobs, ntstar, ndepth,
                                                            iobs, it, id);
                */
                for (i=0; i<6; i++)
                {
                    sacio_copy(ffGrns[kndx], &grns[indices[i]]);
                    if (obs[iobs].pz.lhavePZ)
                    {
                        sacio_copyPolesAndZeros(obs[iobs].pz,
                                                &grns[indices[i]].pz);
                    }
                    sacio_setFloatHeader(SAC_FLOAT_AZ, az,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_BAZ, baz,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPAZ, cmpaz,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPINC, cmpinc,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_EVLA, evla,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_EVLO, evlo,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_STLA, stla,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_STLO, stlo,
                                         &grns[indices[i]].header);
                    sacio_setFloatHeader(SAC_FLOAT_STEL, stel,
                                         &grns[indices[i]].header); 
                    sacio_setCharacterHeader(SAC_CHAR_KNETWK, knetwk,
                                             &grns[indices[i]].header);
                    sacio_setCharacterHeader(SAC_CHAR_KSTNM, kstnm,
                                             &grns[indices[i]].header);
                    sacio_setCharacterHeader(SAC_CHAR_KHOLE, khole,
                                             &grns[indices[i]].header);
                    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, kcmpnms[i],
                                             &grns[indices[i]].header);
                    sacio_setCharacterHeader(SAC_CHAR_KEVNM, "SYNTHETIC\0",
                                             &grns[indices[i]].header);
                    // Set the start time by aligning on the arrival
                    epochNew = epoch + (pickTime - o) - pickTimeGrns;
                    sacio_setEpochalStartTime(epochNew,
                                              &grns[indices[i]].header);
                    // Update the pick times
                    for (l=0; l<nTimeVars; l++)
                    {
                        ierr = sacio_getFloatHeader(pickVars[l],
                                                    grns[indices[i]].header,
                                                    &pick);
                        if (ierr == 0)
                        {
                            pick = pick + o;
                            sacio_setFloatHeader(pickVars[l],
                                                 pick,
                                                 &grns[indices[i]].header);
                        }
                    } 
                }
                //printf("%d %d\n", kndx, indices[0]);
                //printf("%e\n", array_max64f(npts, ffGrns[kndx+0].data));
                ierr = parmt_utils_ff2mtGreens64f(npts, icomp,
                                                  az, baz,
                                                  cmpaz, cmpincSEED,
                                                  ffGrns[kndx+0].data,
                                                  ffGrns[kndx+1].data,
                                                  ffGrns[kndx+2].data,
                                                  ffGrns[kndx+3].data,
                                                  ffGrns[kndx+4].data,
                                                  ffGrns[kndx+5].data,
                                                  ffGrns[kndx+6].data,
                                                  ffGrns[kndx+7].data,
                                                  ffGrns[kndx+8].data,
                                                  ffGrns[kndx+9].data,
                                                  grns[indices[0]].data,
                                                  grns[indices[1]].data,
                                                  grns[indices[2]].data,
                                                  grns[indices[3]].data,
                                                  grns[indices[4]].data,
                                                  grns[indices[5]].data);
                if (ierr != 0)
                {
                    log_errorF("%s: Failed to rotate Greens functions\n", fcnm);
                }
                // Fix the characteristic magnitude scaling in CPS 
                // N.B. this is done early now
                /*
                for (i=0; i<6; i++)
                {
                    cblas_dscal(npts, xscal, grns[indices[i]].data, 1); 
                    //printf("%d %d %e\n", kndx, indices[i],
                    //          array_max64f(npts,grns[indices[i]].data));
                }
                */
            }
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Convenience function which returns the index of the Green's
 *        function on the Green's function structure.
 *
 * @param[in] GMT_TERM    Name of the desired Green's function:
 *                        (G11_TERM, G22_TERM, ..., G23_TERM).
 * @param[in] nobs        Number of observations.
 * @param[in] ntstar      Number of t*'s.
 * @param[in] ndepth      Number of depths.
 * @param[in] iobs        Desired observation number (C numbering).
 * @param[in] itstar      Desired t* (C numbering).
 * @param[in] idepth      Desired depth (C numbering).
 *
 * @result Negative indicates failure.  Otherwise, this is the index in 
 *         grns.grns corresponding to the desired 
 *         (iobs, idepth, itstar, G??_GRNS) coordinate.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_greens_getHudson96GreensFunctionIndex(
    const enum prepmtGreens_enum GMT_TERM,
    const int nobs, const int ntstar, const int ndepth,
    const int iobs, const int itstar, const int idepth)
{
    const char *fcnm = "prepmt_greens_getGreensHudson96FunctionIndex\0";
    int igx, indx, ngrns;
    ngrns = 6*nobs*ntstar*ndepth;
    indx =-1;
    igx = (int) GMT_TERM - 1;
    if (igx < 0 || igx > 5)
    {
        log_errorF("%s: Can't classify Green's functions index\n", fcnm);
        return indx;
    }
    indx = iobs*(6*ntstar*ndepth)
         + idepth*(6*ntstar)
         + itstar*6
         + igx;
    if (indx < 0 || indx >= ngrns)
    {
        log_warnF("%s: indx out of bounds - segfault is coming\n", fcnm);
        return -1;
    }
    return indx;
}
//============================================================================//
/*!
 * @brief Repicks the Green's functions onset time with an STA/LTA picker.
 *
 * @param[in] sta        Short term average window length (seconds).
 * @param[in] lta        Long term average window length (seconds).
 * @param[in] threshPct  Percentage of max STA/LTA after which an arrival
 *                       is declared.
 * @param[in] nobs       Number of observations.
 * @param[in] ntstar     Number of t*'s.
 * @param[in] ndepth     Number of depths.
 * @param[in] iobs       C numbered observation index.
 * @param[in] itstar     C numbered t* index.
 * @param[in] idepth     C numbered depth index.
 *
 * @param[in,out] grns   On input contains the Green's functions.
 *                       On output the first defined arrival time has been
 *                       modified with an STA/LTA picker.  The header variable
 *                       to be modified is most likely SAC_FLOAT_A.
 *                       This is an array of dimension [6]. 
 *
 * @brief 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_greens_repickGreensWithSTALTA(
    const double sta, const double lta, const double threshPct,
    struct sacData_struct *grns)
{
    const char *fcnm = "prepmt_greens_repickGreensWithSTALTA\0";
    struct stalta_struct stalta;
    enum sacHeader_enum pickHeader;
    double *charFn, *g, *Gxx, *Gyy, *Gzz, *Gxy, *Gxz, *Gyz,
           *gxxPad, *gyyPad, *gzzPad, *gxyPad, *gxzPad, *gyzPad,
           apick, charMax, dt, tpick;
    int ierr, k, npts, nlta, npad, nsta, nwork, prePad;
    const int nTimeVars = 11;
    const enum sacHeader_enum pickVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    // Check STA/LTA 
    ierr = 0;
    memset(&stalta, 0, sizeof(struct stalta_struct));
    if (lta < sta || sta < 0.0)
    {
        if (lta < sta){log_errorF("%s: Error lta < sta\n", fcnm);}
        if (sta < 0.0){log_errorF("%s: Error sta is negative\n", fcnm);}
        return -1;
    }
    ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                  grns[0].header, &npts);
    if (ierr != 0 || npts < 1)
    {
        if (ierr != 0)
        {
            log_errorF("%s: Error getting number of points from header\n",
                       fcnm);
        }
        else
        {
            log_errorF("%s: Error no data points\n", fcnm);
        }
        return -1;
    }
    ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                grns[0].header, &dt);
    if (ierr != 0 || dt <= 0.0)
    {
        if (ierr != 0){log_errorF("%s: failed to get dt\n", fcnm);}
        if (dt <= 0.0){log_errorF("%s: invalid sampling period\n", fcnm);}
        return -1;
    }
    // Define the windows
    nsta = (int) (sta/dt + 0.5);
    nlta = (int) (lta/dt + 0.5);
    prePad = MAX(64, fft_nextpow2(nlta));
    npad = prePad + npts;
    // Set space
    gxxPad = memory_calloc64f(npad);
    gyyPad = memory_calloc64f(npad);
    gzzPad = memory_calloc64f(npad);
    gxyPad = memory_calloc64f(npad);
    gxzPad = memory_calloc64f(npad);
    gyzPad = memory_calloc64f(npad);
    charFn = memory_calloc64f(npad);
    // Reference pointers
    Gxx = grns[0].data;
    Gyy = grns[1].data;
    Gzz = grns[2].data;
    Gxy = grns[3].data;
    Gxz = grns[4].data;
    Gyz = grns[5].data;
    // Pre-pad signals
    array_set64f_work(prePad, Gxx[0], gxxPad);
    array_set64f_work(prePad, Gyy[0], gyyPad);
    array_set64f_work(prePad, Gzz[0], gzzPad); 
    array_set64f_work(prePad, Gxy[0], gxyPad);
    array_set64f_work(prePad, Gxz[0], gxzPad);
    array_set64f_work(prePad, Gyz[0], gyzPad);
    // Copy rest of array
    array_copy64f_work(npts, Gxx, &gxxPad[prePad]);
    array_copy64f_work(npts, Gyy, &gyyPad[prePad]);
    array_copy64f_work(npts, Gzz, &gzzPad[prePad]);
    array_copy64f_work(npts, Gxy, &gxyPad[prePad]);
    array_copy64f_work(npts, Gxz, &gxzPad[prePad]);
    array_copy64f_work(npts, Gyz, &gyzPad[prePad]);
    // apply the sta/lta
    for (k=0; k<6; k++)
    {
        g = NULL;
        if (k == 0)
        {
            g = gxxPad;
        }
        else if (k == 1)
        {
            g = gyyPad;
        }
        else if (k == 2)
        {
            g = gzzPad;
        }
        else if (k == 3)
        {
            g = gxyPad;
        }
        else if (k == 4)
        {
            g = gxzPad;
        }
        else if (k == 5)
        {
            g = gyzPad;
        }
        ierr = stalta_setShortAndLongTermAverage(nsta, nlta, &stalta);
        if (ierr != 0)
        {
            printf("%s: Error setting STA/LTA\n", fcnm);
            break;
        }
        ierr = stalta_setData64f(npad, g, &stalta);
        if (ierr != 0)
        {
            printf("%s: Error setting data\n", fcnm);
            break;
        }
        ierr = stalta_applySTALTA(&stalta);
        if (ierr != 0)
        {
            printf("%s: Error applying STA/LTA\n", fcnm);
            break;
        }
        ierr = stalta_getData64f(stalta, npad, &nwork, g);
        if (ierr != 0)
        {
            printf("%s: Error getting result\n", fcnm);
            break;
        }
        cblas_daxpy(npad, 1.0, g, 1, charFn, 1);
        stalta_resetInitialConditions(&stalta);
        stalta_resetFinalConditions(&stalta);
        g = NULL;
    }
    // Compute the pick time
    charMax = array_max64f(npts, &charFn[prePad]);
    tpick =-1.0;
    for (k=prePad; k<npad; k++)
    {
        if (charFn[k] > 0.01*threshPct*charMax)
        {
            tpick = (double) (k - prePad)*dt;
            break;
        }
    }
    if (tpick ==-1.0)
    {
        tpick = (double) (array_argmax64f(npad, charFn) - prePad)*dt;
    }
    // Overwrite the pick
    pickHeader = SAC_UNKNOWN_HDRVAR;
    for (k=0; k<nTimeVars; k++)
    {
        if (sacio_getFloatHeader(pickVars[k], grns[k].header, &apick) == 0)
        {
            pickHeader = pickVars[k];
            break;
         }
    } 
    if (pickHeader == SAC_UNKNOWN_HDRVAR)
    {
        log_warnF("%s: Could not locate primary arrival - assume A\n", fcnm);
        pickHeader = SAC_FLOAT_A;
    } 
//printf("%f %f\n", tpick, apick);
    //double apick;
    //sacio_getFloatHeader(SAC_FLOAT_A, grns[0].header, &apick);
    for (k=0; k<6; k++)
    {
        sacio_setFloatHeader(pickHeader, tpick,
                             &grns[k].header);
    }
    // Dereference pointers and free space
    Gxx = NULL;
    Gyy = NULL;
    Gzz = NULL;
    Gxy = NULL;
    Gxz = NULL;
    Gyz = NULL;
    memory_free64f(&gxxPad);
    memory_free64f(&gyyPad);
    memory_free64f(&gzzPad);
    memory_free64f(&gxyPad);
    memory_free64f(&gxzPad);
    memory_free64f(&gyzPad);
    memory_free64f(&charFn);
    return ierr;
}
//============================================================================//
int prepmt_greens_processHudson96Greens(
    const int nobs, const int ntstar, const int ndepth,
    const struct prepmtCommands_struct cmds,
    struct sacData_struct *grns)
{
    const char *fcnm = "prepmt_greens_processHudson96Greens\0";
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
            log_errorF("%s: Error setting serial command string\n", fcnm);
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
            log_errorF("%s: Error setting the parallel commands\n", fcnm);
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
                                                          iobs, idep, it);
                ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                              grns[kndx].header, &npts);
                if (ierr != 0)
                {
                    log_errorF("%s: Error getting npts\n", fcnm);
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
            log_errorF("%s: Invalid workspace size: %d\n", fcnm, nwork);
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
                    log_errorF("%s: Error getting index\n", fcnm);
                    goto ERROR;
                }
                for (i=0; i<6; i++)
                {
                    i1 = idep*6*ntstar + it*6 + i;
                    //if (i == 0){printf("fill: %d %d %d %e\n", i1, indices[i], dataPtr[i1], array_max64f(grns[indices[i]].npts, grns[indices[i]].data));}

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
            log_errorF("%s: Error setting data\n", fcnm);
            goto ERROR;
        }
        // Apply the commands
        ierr = process_applyParallelCommands(&parallelCommands);
        if (ierr != 0)
        {
            log_errorF("%s: Error processing data\n", fcnm);
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
        //printf("%d\n", nsuse);
        if (ierr != 0)
        {
            log_errorF("%s: Error getting data\n", fcnm);
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
                    log_errorF("%s: Error getting index\n", fcnm);
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
                                                  &G[dataPtr[i1]],
                                                  grns[indices[i]].data);
                    }
                    else
                    {
                        ierr = array_copy64f_work(npts,
                                                  &G[dataPtr[i1]],
                                                  grns[indices[i]].data);
                        //if (i == 0){printf("%d %d %e\n", iobs, indices[i], array_max64f(grns[indices[i]].npts, grns[indices[i]].data));}
                    }
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
int prepmt_greens_writeArchive(const char *archiveName, //const char *archiveDir, const char *projnm,
                               const int nwaves, const int ndepths,
                               const double evla, const double evlo,
                               const double *__restrict__ depths,
                               const struct sacData_struct *sac,
                               const struct sacData_struct *sacGrns)
{
    const char *fcnm = "prepgrns_writeArchive\0";
    int id, ierr, indx, k;
    hid_t h5fl;
    // initialize the archive
    ierr = prepmt_dataArchive_createArchive(archiveName, //archiveDir, projnm,
                                            ndepths, evla, evlo, depths);
    if (ierr != 0)
    {    printf("%s: Error creating archive\n", fcnm);
         return -1;
    }
    // open it for writing
    h5fl = prepmt_dataArchive_openArchive(archiveName, &ierr); //archiveDir, projnm, &ierr);
    // write each observation and greens fns
    for (k=0; k<nwaves; k++)
    {
        ierr = prepmt_dataArchive_addObservation(h5fl, sac[k]);
        for (id=0; id<ndepths; id++)
        {
            indx = k*ndepths*6 + id*6;
            ierr = prepmt_dataArchive_addGreensFunctions(h5fl, sac[k],
                                              sacGrns[indx+0], sacGrns[indx+1],
                                              sacGrns[indx+2], sacGrns[indx+3],
                                              sacGrns[indx+4], sacGrns[indx+5]);
            //printf("%f %e\n", sacGrns[indx+0].header.evdp, sacGrns[indx+0].data[12]);
            if (ierr != 0){return -1;}
        }
    }
    prepmt_dataArchive_closeArchive(h5fl);
    return ierr;
}
//============================================================================//
int prepmt_greens_cutHudson96FromData(const int nobs,
                                      const struct sacData_struct *data,
                                      const int ndepth, const int ntstar,
                                      struct sacData_struct *grns)
{
    const char *fcnm = "prepmt_greens_cutHudson96FromData\0";
    char **newCmds;
    struct prepmtModifyCommands_struct options;
    double cut0, dt, epoch, epochData;
    int i, ierr, idep, indx, iobs, it, npts;
    const int ncmds = 2;
    //const char *cmds[2] = {"cut\0", "demean\0"};
    struct prepmtCommands_struct prepMTcmds;
    newCmds = (char **) calloc(2, sizeof(char *));
    for (i=0; i<2; i++)
    {
        newCmds[i] = (char *) calloc(MAX_CMD_LEN, sizeof(char));
    }
    strcpy(newCmds[1], "demean\0");

    prepMTcmds.cmds = (struct prepmtCommandsChars_struct *)
                      calloc(1, sizeof(struct prepmtCommandsChars_struct)); 
    prepMTcmds.cmds[0].ncmds = ncmds;
    memset(&options, 0, sizeof(struct prepmtModifyCommands_struct));
    for (iobs=0; iobs<nobs; iobs++)
    { 
        // Get the primary pick
        ierr = sacio_getEpochalStartTime(data[iobs].header, &epoch);
        if (ierr != 0)
        {
            log_errorF("%s: Error getting start time\n", fcnm);
            break;
        }
        epochData = epoch;
        sacio_getFloatHeader(SAC_FLOAT_DELTA, data[iobs].header, &dt);
        sacio_getIntegerHeader(SAC_INT_NPTS, data[iobs].header, &npts); 
        cut0 = epoch;
        //cut1 = cut0 + dt*(double) (npts - 1);
        for (idep=0; idep<ndepth; idep++)
        {
            for (it=0; it<ntstar; it++)
            {
                indx = prepmt_greens_getHudson96GreensFunctionIndex(G11_GRNS,
                                                           nobs, ntstar, ndepth,
                                                           iobs, it, idep);
//printf("%d %d %d %d\n", iobs, it, idep, indx);
                sacio_getEpochalStartTime(grns[indx].header, &epoch);
                if (cut0 - epoch < 0.0)
                {
                    log_warnF("%s: Cut may be funky\n", fcnm);
                }
                options.cut0 = fmax(0.0, cut0 - epoch);
                options.cut1 = options.cut0 + dt*(double) (npts - 1);
/*
                newCmds = prepmt_commands_modifyCommands(ncmds, cmds,
                                                         options, grns[indx],
                                                         &ierr);
*/
                ierr = cut_cutEpochalTimesToString(grns[indx].header.delta,
                                                   0.0,
                                                   options.cut0,
                                                   options.cut1,
                                                   newCmds[0]);
                //strcpy(newCmds[1], "demean");
                //printf("%s %f\n", newCmds[0], grns[indx].header.delta);
                prepMTcmds.cmds[0].cmds = newCmds;
                ierr = prepmt_greens_processHudson96Greens(1, 1, 1,
                                                           prepMTcmds,
                                                           &grns[indx]);
                // Enforce the epochal time
                for (i=0; i<6; i++)
                {
                    sacio_setEpochalStartTime(epochData, &grns[indx+i].header);
                }
//printf("%d %e\n", indx, array_max64f(grns[indx].npts, grns[indx].data));
                for (i=0; i<ncmds; i++)
                {
               //     free(newCmds[i]);
                }
                //free(newCmds);
            }
        }
    } 
    for (i=0; i<2; i++)
    {
        free(newCmds[i]);
    }
    free(newCmds[i]);
    free(prepMTcmds.cmds);
    return 0; 
}
//============================================================================//
/*!
 * @brief Convenience function for extracting the: 
 *        \$ \{ G_{xx}, G_{yy}, G_{zz}, G_{xy}, G_{xz}, G_{yz} \} \$
 *        Green's functions indices for the observation, t*, and depth.
 *
 * @param[in] nobs      Total number of observations.
 * @param[in] ntstar    Total number of t*.
 * @param[in] ndepth    Total number of depths.
 * @param[in] iobs      Observation number.  This is C indexed.
 * @param[in] itstar    t* index.  This is C indexed.
 * @param[in] idepth    Depth index.  This is C indexed.
 * @param[in] grns      Contains the Green's functions.
 * @param[out] indices  Contains the Green's functions indices defining
 *                      the indices that return the:
 *                       \$ \{ G_{xx}, G_{yy}, G_{zz}, 
 *                             G_{xy}, G_{xz}, G_{yz} \} \$
 *                      for this observation, t*, and depth. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_greens_getHudson96GreensFunctionsIndices(
    const int nobs, const int ntstar, const int ndepth,
    const int iobs, const int itstar, const int idepth,
    int indices[6])
{
    const char *fcnm = "prepmt_greens_getGreensFunctionsIndices\0";
    int i, ierr;
    const enum prepmtGreens_enum mtTerm[6] = 
       {G11_GRNS, G22_GRNS, G33_GRNS, G12_GRNS, G13_GRNS, G23_GRNS};
    ierr = 0;
    for (i=0; i<6; i++)
    {
        indices[i] = prepmt_greens_getHudson96GreensFunctionIndex(mtTerm[i],
                                                          nobs, ntstar, ndepth,
                                                          iobs, itstar, idepth);
        if (indices[i] < 0){ierr = ierr + 1;}
    }
    if (ierr != 0)
    {
        log_errorF("%s: Error getting indices (obs,t*,depth)=(%d,%d,%d)\n",
                   fcnm, iobs, itstar, idepth);
    }
    return ierr;
}
//============================================================================//
int prepmt_greens_xcAlignGreensToData(const struct sacData_struct data,
                                      const bool luseEnvelope, const bool lnorm,
                                      const double maxTimeLag,
                                      struct sacData_struct *grns)
{
    const char *fcnm = "prepmt_greens_xcAlignGreensToData\0";
    double *dataPad, *Gxx, *Gyy, *Gzz, *Gxy, *Gxz, *Gyz, *xc,
           dt, dtGrns, epochData, epochGrns;
    int ierr, insertData, lxc, maxShift, npadData, npts, npgrns;
    ierr = 0;
    sacio_getIntegerHeader(SAC_INT_NPTS, data.header, &npts);
    sacio_getIntegerHeader(SAC_INT_NPTS, grns[0].header, &npgrns);
    if (npts < 1 || npgrns < 1)
    {
        printf("%s: Error no data\n", fcnm);
        printf("%s: Error no grns functions\n", fcnm);
        return -1;
    }
    sacio_getFloatHeader(SAC_FLOAT_DELTA, data.header, &dt);
    sacio_getFloatHeader(SAC_FLOAT_DELTA, data.header, &dtGrns);
    if (fabs(dt - dtGrns) > 1.e-6 || dt <= 0.0)
    {
        if (dt <= 0.0){printf("%s: dt is invalid\n", fcnm);}
        if (fabs(dt - dtGrns) > 1.e-6)
        {
            printf("%s: dt is inconsistent\n", fcnm);
        }
        return -1;
    }
    // Need a common reference - choose the start time
    ierr = sacio_getEpochalStartTime(data.header,    &epochData);
    if (ierr != 0)
    {
        printf("%s; Error getting data start time\n", fcnm);
        return -1;
    }
    ierr = sacio_getEpochalStartTime(grns[0].header, &epochGrns);
    if (ierr != 0)
    {
        printf("%s: Error getting grns start time\n", fcnm);
        return -1;
    }
    if (epochData < epochGrns)
    {
        printf("%s: epochData < epochGrns not programmed\n", fcnm);
        return -1;
    }
    // Insert the data in a epochal time aligned array
    insertData = (int) ((epochData - epochGrns)/dt + 0.5);
    npadData = insertData + npts;
    //printf("%d %d\n", insertData, npts);
    dataPad = memory_calloc64f(npadData);
    array_copy64f_work(npts, data.data, &dataPad[insertData]);
    // Create pointers to Green's functions 
    Gxx = grns[0].data;
    Gyy = grns[1].data;
    Gzz = grns[2].data;
    Gxy = grns[3].data;
    Gxz = grns[4].data;
    Gyz = grns[5].data;
    // Compute hte max lag
    maxShift =-1;
    if (maxTimeLag >= 0.0){maxShift = (int) (maxTimeLag/dt + 0.5);}
//printf("%d %f %f\n", maxShift, maxTimeLag/dt, maxTimeLag);
    // Align
    lxc = npadData + npgrns - 1; //npts + npgrns - 1;
    xc = memory_calloc64f(lxc);
    ierr = prepmt_greens_xcAlignGreensToData_work(npadData, dataPad, //data.data,
                                                  npgrns, luseEnvelope, lnorm,
                                                  maxShift,
                                                  Gxx, Gyy, Gzz,
                                                  Gxy, Gxz, Gyz,
                                                  lxc, xc);
    memory_free64f(&xc);
    memory_free64f(&dataPad);
    // Dereference pointers
    Gxx = NULL;
    Gyy = NULL;
    Gzz = NULL;
    Gxy = NULL;
    Gxz = NULL;
    Gyz = NULL;
    return ierr;
}
//============================================================================//
int prepmt_greens_xcAlignGreensToData_work(
    const int npts, double *__restrict__ data,
    const int npgrns,
    const bool luseEnvelope, const bool lnorm,
    const int maxShift,
    double *__restrict__ Gxx,
    double *__restrict__ Gyy,
    double *__restrict__ Gzz,
    double *__restrict__ Gxy,
    double *__restrict__ Gxz,
    double *__restrict__ Gyz,
    const int lxc, double *__restrict__ xc)
{
    const char *fcnm = "utils_alignGreensToData\0";
    double *dwork, *Gwork, *xcorrWork, egrns, esig, xdiv;
    int i, ierr, ierr1, kmt, l1, l2, lag, lcref, nlag;
    //------------------------------------------------------------------------//
    ierr = 0;
    lcref = npgrns + npts - 1;
    if (lcref > lxc){return -1;}
    xcorrWork = memory_calloc64f(lcref);
    array_zeros64f_work(lxc, xc);
    if (luseEnvelope)
    {
        Gwork = memory_calloc64f(npgrns);
        dwork = memory_calloc64f(npts);
        ierr = signal_filter_envelope64f_work(npts, data, dwork);
    }
    else
    {
        Gwork = NULL;
        dwork = data;
    }
    esig = cblas_dnrm2(npts, data, 1);
    // Loop on the Green's functions
    for (kmt=0; kmt<6; kmt++)
    {
        if (kmt == 0)
        {
            if (!luseEnvelope)
            {
                Gwork = Gxx;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gxx, Gwork);
            }
        }
        else if (kmt == 1)
        {
            if (!luseEnvelope)
            {
                Gwork = Gyy;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gyy, Gwork);
            }
        }
        else if (kmt == 2)
        {
            if (!luseEnvelope)
            {
                Gwork = Gzz;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gzz, Gwork);
            }
        }
        else if (kmt == 3)
        {
            if (!luseEnvelope)
            {
                Gwork = Gxy;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gxy, Gwork);
            }
        }
        else if (kmt == 4)
        {
            if (!luseEnvelope)
            {
                Gwork = Gxz;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gxz, Gwork);
            }
        }
        else
        {
            if (!luseEnvelope)
            {
                Gwork = Gyz;
            }
            else
            {
                signal_filter_envelope64f_work (npgrns, Gyz, Gwork);
            }
        }
        egrns = cblas_dnrm2(npgrns, Gwork, 1);
        ierr1 = signal_convolve_correlate64f_work(npts, dwork,
                                                  npgrns,
                                                  Gwork,
                                                  CONVCOR_FULL,
                                                  lcref, xcorrWork);
        if (ierr1 != 0)
        {
            if (ierr1 != 0)
            {
                printf("%s: Error correlating kmt %d\n", fcnm, kmt);
            }
            ierr = ierr + 1;
        }
        // sum in the result - absolute value handles polarity
        else
        {
            xdiv = 1.0;
            if (lnorm){xdiv = 1.0/(esig*egrns);}
            if (luseEnvelope)
            {
                #pragma omp simd
                for (i=0; i<lcref; i++)
                {
                    xc[i] = xc[i] + xcorrWork[i]*xdiv;
                }
            }
            else
            {
                #pragma omp simd
                for (i=0; i<lcref; i++)
                {
                    xc[i] = xc[i] + fabs(xcorrWork[i]*xdiv);
//if (kmt == 5){printf("%e\n", xc[i]);}
                }
            }
        }
    } // Loop on green's functions
    // Compute the lag where the times range from [-npgrns:npts].  Hence,
    // lag = npgrns is the unlagged cross-correlation. 
    if (maxShift < 0)
    {
        lag = array_argmax64f(lcref, xc);
    }
    else
    {
        l1 = MAX(0, npgrns - maxShift);
        l2 = MIN(lcref - 1, npgrns + maxShift);
        nlag = l2 - l1 + 1;
        lag = l1 + array_argmax64f(nlag, &xc[l1]); 
    }
    lag =-npgrns + lag;
    //printf("%d %d %d\n", lag, npgrns, npts);
/*
for (i=0; i<npts; i++)
{
printf("%e\n", data[i]);
}
for (i=0; i<npgrns; i++)
{
printf("%e %e %e %e %e %e\n", Gxx[i], Gyy[i], Gzz[i], Gxy[i], Gxz[i], Gyz[i]); 
}
for (i=0; i<lcref; i++)
{
printf("%d %e\n",-npgrns+i, xc[i]); 
}
printf("%d\n", lag);
*/
//getchar();
    // Shift the greens' functions
    array_zeros64f_work(lcref, xcorrWork);
    shiftGreens(npgrns, lag, Gxx, xcorrWork);
    shiftGreens(npgrns, lag, Gyy, xcorrWork);
    shiftGreens(npgrns, lag, Gzz, xcorrWork);
    shiftGreens(npgrns, lag, Gxy, xcorrWork);
    shiftGreens(npgrns, lag, Gxz, xcorrWork);
    shiftGreens(npgrns, lag, Gyz, xcorrWork);
    if (luseEnvelope)
    {
        memory_free64f(&Gwork);
        memory_free64f(&dwork);
    }
    else
    {
        dwork = NULL;
        Gwork = NULL;
    }
    memory_free64f(&xcorrWork);
    return 0;
}
//============================================================================//
static int getPrimaryArrival(const struct sacHeader_struct hdr,
                             double *time, char phaseName[8])
{
    const char *fcnm = "getPrimaryArrival\0";
    const enum sacHeader_enum timeVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    const enum sacHeader_enum timeVarNames[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};
    int i, ifound1, ifound2; 
    memset(phaseName, 0, 8*sizeof(char));
    for (i=0; i<11; i++)
    {
        ifound1 = sacio_getFloatHeader(timeVars[i], hdr, time);
        ifound2 = sacio_getCharacterHeader(timeVarNames[i], hdr, phaseName); 
        if (ifound1 == 0 && ifound2 == 0){return 0;}
    }
    printf("%s: Failed to get primary pick\n", fcnm);
    *time =-12345.0;
    memset(phaseName, 0, 8*sizeof(char));
    strcpy(phaseName, "-12345"); 
    return -1;
}
//============================================================================//
static int getPrimaryArrivalNoPolarity(const struct sacHeader_struct hdr,
                                       double *time, char phaseName[8])
{
    int ierr;
    size_t lenos;
    ierr = getPrimaryArrival(hdr, time, phaseName);
    if (ierr == 0)
    {
        lenos = strlen(phaseName);
        if (lenos > 0)
        {
            if (phaseName[lenos-1] == '+' || phaseName[lenos-1] == '-')
            {
                phaseName[lenos-1] = '\0';
            }
        }
    }
    return ierr;
}
//============================================================================//
static void shiftGreens(const int npgrns, const int lag,
                        double *__restrict__ G, double *__restrict__ work)
{
    int ncopy; 
    if (lag == 0){return;}
    array_copy64f_work(npgrns, G, work);
    array_zeros64f_work(npgrns, G);
    if (lag > 0)
    {
        ncopy = npgrns - lag;
        array_copy64f_work(ncopy, work, &G[lag]);
    }
    else
    {
        ncopy = npgrns + lag;
        array_copy64f_work(ncopy, &work[-lag], G);
    }
    return;
}
