#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "prepmt/prepmt_hudson96.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "cps.h"
#include "iscl/array/array.h"
#include "iscl/fft/fft.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

/*!
 * @brief Reads the ini file for Computer Programs in Seismology hudson96
 *        forward modeling variables.
 *
 * @param[in] iniFile    Name of ini file.
 * @param[in] section    Section of ini file to read.
 *
 * @param[out] parms     The hudson96 parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_hudson96_readHudson96Parameters(const char *iniFile,
                                           const char *section,
                                           struct hudson96_parms_struct *parms)
{
    const char *s;
    char vname[256];
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHudson96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: ini file: %s does not exist\n", __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        fprintf(stderr, "%s: Cannot parse ini file\n", __func__);
        return -1;
    }
    // Teleseismic model
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:modeltel", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (s != NULL)
    {
        strcpy(parms->modeltel, s);
    }
    else
    {
        strcpy(parms->modeltel, "tak135sph.mod\0");
    }
    // Receiver model
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:modelrec", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (s != NULL)
    {
        strcpy(parms->modelrec, s);
    }
    else
    {
        strcpy(parms->modelrec, parms->modelrec);
    }
    // Source model
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:modelsrc", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (s != NULL)
    {
        strcpy(parms->modelsrc, s);
    }
    else
    {
        strcpy(parms->modelsrc, parms->modelsrc);
    }
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:hs", section);
    parms->hs = iniparser_getdouble(ini, vname, 0.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dt", section);
    parms->dt = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:npts", section);
    parms->npts = iniparser_getint(ini, vname, 1024);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:gcarc", section);
    parms->gcarc = iniparser_getdouble(ini, vname, 50.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:offset", section);
    parms->offset = iniparser_getdouble(ini, vname, 10.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dosrc", section);
    parms->dosrc = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dorec", section);
    parms->dorec = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dotel", section);
    parms->dotel = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dop", section);
    parms->dop = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dokjar", section);
    parms->dokjar = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:loffsetdefault", section);
    parms->loffsetdefault = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:verbose", section);
    parms->verbose = iniparser_getboolean(ini, vname, 0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:utstar", section);
    parms->utstar = iniparser_getdouble(ini, vname, -12345.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dottonly", section);
    parms->dottonly = iniparser_getboolean(ini, vname, 0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:zsrc", section);
    parms->zsrc = iniparser_getdouble(ini, vname, 100.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:zrec", section);
    parms->zrec = iniparser_getdouble(ini, vname, 60.0);
    // Free ini dictionary
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the fundamental fault Green's functions for the teleseismic
 *        body waves with hudson96.
 *
 * @param[in] nobs            Number of observations.
 * @param[in] obs             Observations with source location and
 *                            receiver location.  This is an array of length
 *                            nobs.
 * @param[in] luseCrust1      If true then use crust1.0 models at source
 *                            and receiver.
 * @param[in] crust1Dir       If luseCrust1 is true then this is the name
 *                            of the crust1.0 directory.  If NULL and 
 *                            luseCrust1 is used then it will default to
 *                            libcps configuration.
 * @param[in] luseSrcModel    If true then use the local source velocity
 *                            at the source.  This supersedes the crust1.0
 *                            source velocity if it is defined.
 * @param[in] srcModelName    If true then this is the name of CPS style
 *                            1D local source velocity model. 
 * @param[in] ntstar          Number of t*'s in grid.
 *
 * @param[out] telmod         Holds the teleseismic model (ak135).
 * @param[out] srcmod         Holds the local source model.
 * @param[in,out] recmod      On input contains sufficient space and should
 *                            be an array of length nobs.
 *                            On output holds the receiver models.
 *
 * @param[out] ierr           0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 * @bug I need better handling of the offset.
 *
 */
int hudson96_getModels(const int nobs, const struct sacData_struct *obs,
                       const bool luseCrust1, const char *crust1Dir,
                       const bool luseSrcModel, const char *sourceModel,
                       struct vmodel_struct *telmod,
                       struct vmodel_struct *srcmod, 
                       struct vmodel_struct *recmod)
{
    double *lats, *lons, evla, evlo;
    int ierr, iobs;
    ierr = 0;
    // Set the teleseismic model
    memset(srcmod, 0, sizeof(struct vmodel_struct));
    memset(telmod, 0, sizeof(struct vmodel_struct));
    cps_globalModel_ak135f(telmod);
    if (luseCrust1)
    {
        fprintf(stdout, "%s: Reading crust1.0...\n", __func__);
        lats = memory_calloc64f(nobs);
        lons = memory_calloc64f(nobs);
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVLA, obs[0].header, &evla);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error - evla not set\n", __func__);
            return ierr;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVLO, obs[0].header, &evlo);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error - evlo not set\n", __func__);
            return ierr;
        }
        for (iobs=0; iobs<nobs; iobs++)
        {
            sacio_getFloatHeader(SAC_FLOAT_STLA, obs[iobs].header,
                                 &lats[iobs]);
            sacio_getFloatHeader(SAC_FLOAT_STLO, obs[iobs].header,
                                 &lons[iobs]);
        }
        ierr = cps_crust1_getCrust1ForHerrmann(crust1Dir, false,
                                               evla, evlo,
                                               lats, lons,
                                               nobs,
                                               srcmod, recmod);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to load crust1.0 model\n", __func__);
            return ierr;
        }
        memory_free64f(&lats);
        memory_free64f(&lons);
    }
    // Use teleseismic model for receiver model
    else
    {
        fprintf(stdout, "%s: Setting receiver models to teleseismic model\n",
                __func__);
        for (iobs=0; iobs<nobs; iobs++)
        {
            cps_utils_copyVmodelStruct(*telmod, &recmod[iobs]);
        }
    }
    // Use source model
    if (luseSrcModel)
    {
        fprintf(stdout, "%s: Using source model: %s\n", __func__, sourceModel);
        cps_utils_freeVmodelStruct(srcmod);
        ierr = cps_getmod(sourceModel, srcmod);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error loading srcmod\n", __func__);
            return ierr;
        }
    }
    fprintf(stdout, "%s: Computing Green's functions...\n", __func__);
    // Otherwise source to teleseismic model
    if (!luseSrcModel && !luseCrust1)
    {
        fprintf(stdout, "%s: Setting source model to teleseismic model\n",
                __func__);
        cps_utils_copyVmodelStruct(*telmod, srcmod);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the fundamental fault Green's functions for the teleseismic
 *        body waves with hudson96.
 *
 * @param[in] hudson96Parms   hudson96 forward modeling parameters.
 * @param[in] hpulse96Parms   hpulse96 forward modeling parameters.
 * @param[in] telmod          Holds the teleseismic model (ak135).
 * @param[in] srcmod          Holds the local source model.
 * @param[in] recmod          Holds each receiver model.
 *                            This is an array of dimension [nobs].
 * @param[in] ntstar          Number of t*'s in grid.
 * @param[in] tstars          Array of dimension [ntstar] with the attenuation
 *                            factors.
 * @param[in] ndepth          Number of source depths in grid.
 * @param[in] depths          Array of dimension [ndepth] with the source
 *                            depths specified in (km).
 * @param[in] nobs            Number of observations.
 * @param[in] obs             Observed waveforms which contain the receiver
 *                            locations, components of motion, first
 *                            arrival time, sampling period, and number of
 *                            points in waveform.
 *
 * @param[out] ierr           0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 * @bug I need better handling of the offset.
 *
 */
struct sacData_struct *prepmt_hudson96_computeGreensFF(
    const struct hudson96_parms_struct hudson96Parms,
    const struct hpulse96_parms_struct hpulse96Parms,
    const struct vmodel_struct telmod,
    const struct vmodel_struct srcmod,
    const struct vmodel_struct *recmod,
    const int ntstar, const double *__restrict__ tstars,
    const int ndepth, const double *__restrict__ depths,
    const int nobs, const struct sacData_struct *obs, int *ierr)
{
    struct hudson96_parms_struct hudson96ParmsWork;
    struct hpulse96_parms_struct hpulse96ParmsWork;
    struct hpulse96_data_struct zresp;
    struct sacData_struct *sacFFGrns;
    struct hwave_greens_struct *ffGrns;
    char sncl[64], phaseName[8];
    double *dptr, cmpaz, cmpinc, dt, gcarc, offset0,
           pickTime, xnorm;
    int i, idep, ierrAll, ierr1, ierr2, indx, iobs, iobs0, ip, it,
        kndx, maxPts, nloop, npts;
    bool lzero[10], lfound, lsh;
    enum isclError_enum isclError;
    const enum sacHeader_enum pickVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    const enum sacHeader_enum pickTypes[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};
    // Here's the idea: I want to input magnitudes with units of N-m so:
    // N-m -> Dyne-dm (1.e+7)
    // CPS internally scales from Dyne-cm to cm (1.e+20)
    // Finally, I want outputs proprotional to m (1.e-2)
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
    const char *cfaults[10] = {"ZDS\0", "ZSS\0", "ZDD", "ZEX\0",
                               "RDS\0", "RSS\0", "RDD", "REX\0",
                               "TDS\0", "TSS\0"};
    const int npTypes = 11;
    const int nrDist = 1;
    const int nrDepths = 1;
    const int nsDepths = 1;
    //------------------------------------------------------------------------//
    //
    // Quick error checks 
    *ierr = 0;
    sacFFGrns = NULL;
    if (nobs < 1 || obs == NULL ||
        ntstar < 1 || tstars == NULL ||
        ndepth < 1 || depths == NULL)
    {
        *ierr = 1;
        if (nobs < 1){fprintf(stderr, "%s: Error no observations\n", __func__);}
        if (obs == NULL){fprintf(stderr, "%s: Error obs is NULL\n", __func__);}
        if (ntstar < 1){fprintf(stderr, "%s: Error no t*'s\n", __func__);}
        if (tstars == NULL)
        {
            fprintf(stderr, "%s: Error tstars is NULL\n", __func__);
        }
        if (ndepth < 1){fprintf(stderr, "%s: Error no depths\n", __func__);}
        if (depths == NULL)
        {
            fprintf(stderr, "%s: Error depths is NULL\n", __func__);
        }
        return sacFFGrns;
    }
    // Set the modeling structures
    memset(&hudson96ParmsWork, 0, sizeof(struct hudson96_parms_struct));
    memset(&hpulse96ParmsWork, 0, sizeof(struct hpulse96_parms_struct));
    cps_setHudson96Defaults(&hudson96ParmsWork);
    cps_setHpulse96Defaults(&hpulse96ParmsWork);
    cps_utils_copyHudson96ParmsStruct(hudson96Parms, &hudson96ParmsWork);
    cps_utils_copyHpulse96ParmsStruct(hpulse96Parms, &hpulse96ParmsWork);
    offset0 = hudson96Parms.offset;
    // Set space for output
    nloop = nobs*ndepth*ntstar;
    sacFFGrns = (struct sacData_struct *)
                 calloc((size_t) (10*nloop), sizeof(struct sacData_struct));
    // Loop on the distances, depths, t*'s and compute greens functions
    iobs0 =-1;
    ierrAll = 0;
    for (indx=0; indx<nloop; indx++)
    {
        ffGrns = NULL;
        // Convert 3D grid into loop indices (i->n1 is fast; k->n3 is slow)
        *ierr = prepmt_hudson96_grd2ijk(indx,
                                        ntstar, ndepth, nobs,
                                        &it, &idep, &iobs);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Failed to convert to grid\n", __func__);
            break;
        }
        // New observation (station) -> may need to update velocity model
        if (iobs != iobs0)
        {
            iobs0 = iobs;
        }
        *ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                     obs[iobs].header, &dt);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Could not get sampling period\n", __func__);
            break;
        }
        *ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                     obs[iobs].header, &gcarc);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Could not get gcarc\n", __func__);
            break;
        }
        *ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                       obs[iobs].header, &npts);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Could not get npts\n", __func__);
            break;
        }
        // Tie waveform modeling to first pick type 
        lfound = false;
        hudson96ParmsWork.dop = true;
        for (ip=0; ip<npTypes; ip++)
        {
            ierr1 = sacio_getCharacterHeader(pickTypes[ip],
                                             obs[iobs].header, phaseName);
            ierr2 = sacio_getFloatHeader(pickVars[ip], 
                                         obs[iobs].header, &pickTime);
            if (ierr1 == 0 && ierr2 == 0)
            {
                lfound = true;
                if (strncasecmp(phaseName, "P", 1) == 0)
                {
                    hudson96ParmsWork.dop = true;
                }
                else if (strncasecmp(phaseName, "S", 1) == 0)
                {
                    hudson96ParmsWork.dop = false;
                }
                else
                {
                    fprintf(stderr, "%s: Can't classify phase: %s\n",
                            __func__, phaseName);
                    lfound = false;
                }
                break;
            }
        }
        if (!lfound)
        {
            fprintf(stdout, "%s: No pick available - won't be able to align\n",
                    __func__);
            continue;
        }
        memset(&zresp, 0, sizeof(struct hpulse96_data_struct));
        hudson96ParmsWork.hs = depths[idep];
        hudson96ParmsWork.dt = dt;
        hudson96ParmsWork.gcarc = gcarc;
        hudson96ParmsWork.utstar = tstars[it];
        maxPts = 256;
        if (hudson96ParmsWork.dt < 0.1){maxPts = 512;}
        hudson96ParmsWork.npts = MAX(maxPts,
                                     MIN(2048,
                                         2*fft_nextpow2(MAX(1, npts), ierr)));
        hudson96ParmsWork.offset = fmin(pickTime, offset0); //TODO: fmin or fmax?
//printf("%s %f\n", phaseName, pickTime);
        // Try to pad this out a little
        hudson96ParmsWork.offset = fmax(hudson96ParmsWork.offset,
                                  (double) (hudson96ParmsWork.npts - 1)*dt*0.2);
/*
        if (pickTime < hudson96ParmsWork.offset)
        {
            hudson96ParmsWork.offset = (double) (int) (pickTime/dt + 0.5)*dt;
        }
*/
        zresp = hudson96_interface(&hudson96ParmsWork,
                                   telmod, recmod[iobs], srcmod, ierr);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Error calling hudson96\n", __func__);
            ierrAll = ierrAll + 1;
            goto NEXT_OBS;
        }
        ffGrns = hpulse96_interface(nrDist, nrDepths, nsDepths,
                                    &hpulse96ParmsWork,
                                    &zresp, ierr);
        if (*ierr != 0)
        {
            fprintf(stderr, "%s: Error calling hpulse96 interface\n", __func__);
            ierrAll = ierrAll + 1;
            goto NEXT_OBS;
        }
        // Set the fundamental fault Green's functions
        for (i=0; i<10; i++)
        {
            dptr = NULL;
            lsh = false;
            lzero[i] = false;
            if (i == 0)
            {
                dptr = ffGrns->zds;
            }
            else if (i == 1)
            {
                dptr = ffGrns->zss;
            }
            else if (i == 2)
            {
                dptr = ffGrns->zdd;
            }
            else if (i == 3)
            {
                dptr = ffGrns->zex;
            }
            else if (i == 4)
            {
                dptr = ffGrns->rds;
                cmpinc = 90.0;
            }
            else if (i == 5)
            {
                dptr = ffGrns->rss;
                cmpinc = 90.0;
            }
            else if (i == 6)
            {
                dptr = ffGrns->rdd;
                cmpinc = 90.0;
            }
            else if (i == 7)
            {
                dptr = ffGrns->rex;
                cmpinc = 90.0;
            }
            else if (i == 8)
            {
                dptr = ffGrns->tds;
                lsh = true;
                cmpaz = 90.0;
                cmpinc = 90.0;
            }
            else if (i == 9)
            {
                dptr = ffGrns->tss;
                lsh = true;
                cmpaz = 90.0;
                cmpinc = 90.0;
            }
            kndx = indx*10 + i;
            sacio_setDefaultHeader(&sacFFGrns[kndx].header);
            sacFFGrns[kndx].header.lhaveHeader = true;
            sacio_setIntegerHeader(SAC_INT_NPTS, ffGrns->npts,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZYEAR, 1970,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZJDAY, 1,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZHOUR, 0,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZMIN, 0,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZSEC, 0,
                                   &sacFFGrns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZMSEC, 0,
                                   &sacFFGrns[kndx].header);
            sacio_setBooleanHeader(SAC_BOOL_LCALDA, false,
                                   &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_DELTA, ffGrns->dt,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_GCARC, ffGrns->dist/111.195,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_DIST, ffGrns->dist,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_EVDP, depths[idep],
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_STEL, ffGrns->stelel,
                                 &sacFFGrns[kndx].header); 
            sacio_setFloatHeader(SAC_FLOAT_AZ,  0.0,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_BAZ, 180.0,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_O, -ffGrns->t0, 
                                 &sacFFGrns[kndx].header); 
            sacio_setFloatHeader(SAC_FLOAT_B, ffGrns->t0,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_CMPAZ, cmpaz,
                                 &sacFFGrns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_CMPINC, cmpinc,
                                 &sacFFGrns[kndx].header);
            if (hudson96ParmsWork.dop)
            {
                sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timep,
                                     &sacFFGrns[kndx].header);
                sacio_setCharacterHeader(SAC_CHAR_KA, "P",
                                         &sacFFGrns[kndx].header);
            }
            else
            {
                if (!lsh)
                {
                    sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timesv,
                                         &sacFFGrns[kndx].header);
                }
                else
                {
                    sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timesh,
                                         &sacFFGrns[kndx].header);
                }
                sacio_setCharacterHeader(SAC_CHAR_KA, "S",
                                         &sacFFGrns[kndx].header);
            }
            sacio_setCharacterHeader(SAC_CHAR_KO, "O\0",
                                     &sacFFGrns[kndx].header);
            sacio_setCharacterHeader(SAC_CHAR_KEVNM, "SYNTHETIC\0",
                                     &sacFFGrns[kndx].header);
            sacio_setCharacterHeader(SAC_CHAR_KCMPNM, cfaults[i],
                                     &sacFFGrns[kndx].header);
            sacFFGrns[kndx].npts = ffGrns->npts;
            sacFFGrns[kndx].data = sacio_malloc64f(ffGrns->npts);
            if (dptr != NULL)
            {
                array_copy64f_work(ffGrns->npts, dptr, sacFFGrns[kndx].data);
            }
            else
            {
                array_set64f_work(ffGrns->npts, 0.0, sacFFGrns[kndx].data);
            }
            xnorm = cblas_dnrm2(sacFFGrns[kndx].npts, sacFFGrns[kndx].data, 1);
            if (xnorm == 0.0){lzero[i] = true;}
            // Handle scaling
            cblas_dscal(sacFFGrns[kndx].npts, xscal, sacFFGrns[kndx].data, 1);
            dptr = NULL;
        } // Loop on fundamnetal faults
        if (array_sum8l(10, lzero, &isclError) == 10)
        {
            memset(sncl, 0, 64*sizeof(char));
            sprintf(sncl, "%s.%s.%s.%s",
                    obs[iobs].header.knetwk, obs[iobs].header.kstnm, 
                    obs[iobs].header.kcmpnm, obs[iobs].header.khole);
            fprintf(stdout,
                    "%s: Warning Greens function %s at gcarc=%e is zero %d\n",
                    __func__, sncl, gcarc, hudson96ParmsWork.npts);
        }
NEXT_OBS:;
        if (ffGrns != NULL)
        {
            cps_utils_freeHwaveGreensStruct(ffGrns);
            free(ffGrns);
        }
        cps_utils_freeHpulse96DataStruct(&zresp);
    }
    return sacFFGrns;
}
//============================================================================//
/*!
 * @brief Converts a 1D index to 3D grid indices:
 *
 *        Changes:
 *        for (igrd=0; igrd<n1*n2*n3; igrd++)
 *        {
 *
 *        }
 *        To: 
 *        for (k=0; k<n3; k++)
 *        {
 *            for (j=0; j<n2; j++)
 *            {
 *                for (i=0; i<n1; i++)
 *                {
 *                    igrd = k*n2*n1 + j*n1 + i;
 *                }
 *            }
 *        }
 *
 * @param[in] igrd  Flattened C numbered grid index: igrd = k*n2*n1 + j*n1 + i.
 * @param[in] n1    This is the fastest loop iterator (innermost loop).
 * @param[in] n2    This is the middle loop iterator (intermediate loop).
 * @param[in] n3    This is the slowest loop iterator (outermost loop).
 * @param[out] i    C index in inner loop [0,n1).
 * @param[out] j    C index in intermediate loop [0,n2).
 * @param[out] k    C index in outermost loop [0,n3).
 *
 * @result 0 indicates success.
 *
 */
int prepmt_hudson96_grd2ijk(const int igrd,
                            const int n1, const int n2, const int n3,
                            int *i, int *j, int *k) 
{
    int ierr, n12;
    ierr = 0;
    n12 = n1*n2;
    *k = (igrd)/n12;
    *j = (igrd - *k*n12)/n1;
    *i =  igrd - *k*n12 - *j*n1;
    if (*i < 0 || *i > n1 - 1){ierr = ierr + 1;}
    if (*j < 0 || *j > n2 - 1){ierr = ierr + 1;} 
    if (*k < 0 || *k > n3 - 1){ierr = ierr + 1;} 
    return ierr; 
}
//============================================================================//
/*!
 * @brief Maps the observation, depth, and t* to the global index.
 *
 * @param[in] ndepth    Number of depths.
 * @param[in] ntstar    Number of t*'s.
 * @param[in] iobs      C indexed observation number [0,nobs-1]
 * @param[in] idep      C indexed depth number [0,ndepth-1]
 * @param[in] it        C indexed t* number [0,ntstar-1]
 *
 * @result If >= 0 then this (iobs, idep, it) index in the Green's functions
 *         table.
 * 
 * @author Ben Baker, ISTI
 *
 */
int prepmt_hudson96_observationDepthTstarToIndex(
    const int ndepth, const int ntstar,
    const int iobs, const int idep, const int it)
{
    int indx;
    indx = iobs*ntstar*ndepth + idep*ntstar + it;
    return indx;
}
