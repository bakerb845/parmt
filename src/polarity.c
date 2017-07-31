#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "parmt_polarity.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "ttimes.h"
#include "iscl/array/array.h"
#include "iscl/geodetic/geodetic.h"
#include "iscl/memory/memory.h"


#define LDG 8
#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif

static void fillBasis(const double i, const double phi,
                      double *__restrict__ gam,
                      double *__restrict__ phat,
                      double *__restrict__ phihat);
static double computeContraction3x3(const double *__restrict__ a,
                                    const double *__restrict__ M,
                                    const double *__restrict__ b);
static void setM3x3(const int k, double *__restrict__ M);
static int computePadding64f(const int n);
static int performPolaritySearch64f(const int nmt, const int ldm, 
                                    const int nobs,
                                    const int blockSize, const int mblock,
                                    const int Mrows, const int Kcols,
                                    const double *__restrict__ Dmat,
                                    const double *__restrict__ G,
                                    const double *__restrict__ Sigma,
                                    const double *__restrict__ mts, 
                                    double *__restrict__ phi);

int parmt_polarity_computeTTimesGreens(
    const MPI_Comm globalComm,
    const struct parmtPolarityParms_struct parms,
    const struct parmtData_struct data,
    struct polarityData_struct *polarityData)
{
    const char *fcnm = "parmt_polarity_computeTTimesGreens64f\0";
    char kt0[8], kcmpnm[8], stat[8];
    double G6[6], *cmpazs, *cmpincs, *deps, *evlas, *evlos,
           *GxxBuf, *GyyBuf, *GzzBuf, *GxyBuf, *GxzBuf, *GyzBuf,
           *stlas, *stlos, cmpinc, cmpaz, stla, stlo;
    int *icomps, *observation, *polarity, *waveType, icomp, ierr, ierrAll,
        iloc, iobs, ipol, it, iwav, jloc, k, kt, myid, nPolarity, nprocs;
    size_t lenos;
    const int master = 0;
    const int nTimeVars = 11;
    const enum sacHeader_enum timeVarNames[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};
    //------------------------------------------------------------------------//
    //
    // Initialize
    ierr = 0;
    cmpazs = NULL;
    cmpincs = NULL;
    deps = NULL;
    evlas = NULL;
    evlos = NULL;
    GxxBuf = NULL;
    GyyBuf = NULL;
    GzzBuf = NULL;
    GxyBuf = NULL;
    GxzBuf = NULL;
    GyzBuf = NULL;
    icomps = NULL;
    observation = NULL;
    polarity = NULL;
    stlas = NULL;
    stlos = NULL;
    waveType = NULL;
    memset(polarityData, 0, sizeof(struct polarityData_struct));
    // ttimes uses read-only file-io and i don't want to fix it
    MPI_Comm_rank(globalComm, &myid);
    MPI_Comm_size(globalComm, &nprocs);
    if (myid != master){goto WAIT;}
    // Verify there is a chance for something to do
    if (data.nobs < 1 || data.nlocs < 1)
    {
        if (data.nobs < 1){printf("%s: No observations\n", fcnm);}
        if (data.nlocs < 1){printf("%s: No locations\n", fcnm);}
        return 0;
    }
    // Extract the depths in the grid-search from the first waveform 
    iobs = 0;
    deps = memory_calloc64f(data.nlocs);
    evlas = memory_calloc64f(data.nlocs);
    evlos = memory_calloc64f(data.nlocs);
    for (iloc=0; iloc<data.nlocs; iloc++)
    {
        k = iobs*data.nlocs + iloc;
        ierr  = sacio_getFloatHeader(SAC_FLOAT_EVDP, data.sacGxx[k].header,
                                     &deps[iloc]);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA, data.sacGxx[k].header,
                                     &evlas[iloc]);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO, data.sacGxx[k].header,
                                     &evlos[iloc]);
        if (ierr != 0)
        {
            printf("%s: Unable to get event coordinates on %d\n", fcnm, myid);
            goto ERROR;
        }
    }
    // Compute the number of usable polarities (i.e. observations with P picks
    // or S picks with first motions)
    nPolarity = 0;
    if (myid == master)
    {
        polarity = memory_calloc32i(data.nobs);
        waveType = memory_calloc32i(data.nobs);
        icomps = memory_calloc32i(data.nobs);
        stlas = memory_calloc64f(data.nobs);
        stlos = memory_calloc64f(data.nobs);
        cmpazs = memory_calloc64f(data.nobs);
        cmpincs = memory_calloc64f(data.nobs);
        observation = array_set32i(data.nobs, -1, &ierr);
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            // Ensure the essential preliminary information is defined
            ierr  = sacio_getFloatHeader(SAC_FLOAT_CMPINC,
                                         data.data[iobs].header, &cmpinc);
            ierr += sacio_getFloatHeader(SAC_FLOAT_CMPAZ,
                                         data.data[iobs].header, &cmpaz);
            ierr += sacio_getFloatHeader(SAC_FLOAT_STLA,
                                         data.data[iobs].header, &stla);
            ierr += sacio_getFloatHeader(SAC_FLOAT_STLO,
                                         data.data[iobs].header, &stlo);
            ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM,
                                             data.data[iobs].header, kcmpnm);
            if (ierr != 0)
            {
                printf("%s: Failed to get header information\n", fcnm);
                continue;
            }
            // SAC to SEED convention
            cmpinc = cmpinc - 90.0;
            // Figure out the component
            lenos = MAX(1, strlen(kcmpnm));
            icomp = 1;
            if (kcmpnm[lenos-1] == 'Z' || kcmpnm[lenos-1] == '1')
            {
                icomp = 1;
            }
            else if (kcmpnm[lenos-1] == 'N' || kcmpnm[lenos-1] == '2')
            {
                icomp = 2;
            }
            else if (kcmpnm[lenos-1] == 'E' || kcmpnm[lenos-1] == '3')
            {
                icomp = 3;
            }
            else
            {
                printf("%s: Cannot classify component %s\n", fcnm, kcmpnm);
                continue; 
            }
            // Get the primary pick
            for (it=0; it<nTimeVars; it++)
            {
                ierr = sacio_getCharacterHeader(timeVarNames[it],
                                                data.data[iobs].header,
                                                kt0);
                if (ierr == 0){break;}
            }
            lenos = strlen(kt0);
            // Ensure the pick is defined
            if (ierr == 0 && lenos > 1)
            {
                iwav = 1;
                if (kt0[0] == 'P' || kt0[0] == 'p')
                {
                    iwav = 1;
                }
                else if (kt0[0] == 'S' || kt0[0] == 's')
                {
                    iwav = 2;
                }
                else
                {
                    // surface waves will commonly be processed and not
                    // have polarities
                    if (kt0[0] == 'R' || kt0[0] == 'r' ||
                        kt0[0] == 'L' || kt0[0] == 'l')
                    {
                        continue;
                    }
                    printf("%s: t0 phase is not a P or S phase %s\n",
                           fcnm, kt0);
                    continue;
                }
                if (kt0[lenos-1] == '+')
                {
                    ipol = 1;
                }
                else if (kt0[lenos-1] == '-')
                {
                    ipol =-1;
                }
                else
                {
                    printf("%s: could not classify polarity %s\n", fcnm, kt0);
                    continue;
                }
                // let user know something happened 
                sacio_getCharacterHeader(SAC_CHAR_KSTNM,
                                         data.data[iobs].header, stat);
                printf("%s: Polarity for %s is %d\n", fcnm, stat, ipol); 
//printf("%f %f %f %f %d %d %d %d\n", stla, stlo, cmpinc, cmpaz, icomp, iobs, iwav, ipol);
                // save information for modeling
                stlas[nPolarity] = stla;
                stlos[nPolarity] = stlo;
                cmpincs[nPolarity] = cmpinc;
                cmpazs[nPolarity] = cmpaz;
                icomps[nPolarity] = icomp;
                observation[nPolarity] = iobs;
                waveType[nPolarity] = iwav;
                polarity[nPolarity] = ipol;
                nPolarity = nPolarity + 1;
                //printf("%d %d %f %f\n", iwav, ipol, cmpinc, cmpaz);
           } // End basic header info check
        } // Loop on observations
    } // End check on myid == master
    // Tell all other processes about the forward modeling information
    MPI_Bcast(&nPolarity, 1, MPI_INT, master, globalComm);
    if (nPolarity == 0)
    {
        if (myid == master){printf("%s: There are no polarities\n", fcnm);}
        ierr = 0;
        goto ERROR;
    }
    if (myid != master)
    {
        stlas = memory_calloc64f(nPolarity);
        stlos = memory_calloc64f(nPolarity);
        cmpincs = memory_calloc64f(nPolarity);
        cmpazs = memory_calloc64f(nPolarity);
        icomps = memory_calloc32i(nPolarity);
        observation = memory_calloc32i(nPolarity);
        waveType = memory_calloc32i(nPolarity);
        polarity = memory_calloc32i(nPolarity);
    }
    else
    {
        printf("%s: Warning - i'm setting wts to unity for now\n", fcnm);
    }
    MPI_Bcast(stlas,   nPolarity, MPI_DOUBLE, master, globalComm);
    MPI_Bcast(stlos,   nPolarity, MPI_DOUBLE, master, globalComm);
    MPI_Bcast(cmpincs, nPolarity, MPI_DOUBLE, master, globalComm);
    MPI_Bcast(cmpazs,  nPolarity, MPI_DOUBLE, master, globalComm);
    MPI_Bcast(icomps,      nPolarity, MPI_INT, master, globalComm);
    MPI_Bcast(observation, nPolarity, MPI_INT, master, globalComm);
    MPI_Bcast(waveType,    nPolarity, MPI_INT, master, globalComm);
    MPI_Bcast(polarity,    nPolarity, MPI_INT, master, globalComm); 
    // Set space
    ierrAll = 0;
    GxxBuf = memory_calloc64f(data.nlocs*nPolarity);
    GyyBuf = memory_calloc64f(data.nlocs*nPolarity);
    GzzBuf = memory_calloc64f(data.nlocs*nPolarity);
    GxyBuf = memory_calloc64f(data.nlocs*nPolarity);
    GxzBuf = memory_calloc64f(data.nlocs*nPolarity);
    GyzBuf = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gxx = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gyy = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gzz = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gxy = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gxz = memory_calloc64f(data.nlocs*nPolarity);
    polarityData->Gyz = memory_calloc64f(data.nlocs*nPolarity);
    // Set the data
    polarityData->polarity = memory_calloc64f(nPolarity);
    for (ipol=0; ipol<nPolarity; ipol++)
    {
        polarityData->polarity[ipol] = (double) polarity[ipol];
    }
    // Set the weights
    // TODO fix me
    polarityData->wts = array_set64f(nPolarity, 1.0, &ierr); 
    // Compute the forward modeling matrix columns
    for (jloc=0; jloc<data.nlocs; jloc=jloc+nprocs)
    {
        iloc = jloc + myid; 
        if (iloc >= data.nlocs){continue;}
        // Compute the forward modeling matrix for this observation group
        for (ipol=0; ipol<nPolarity; ipol++)
        {
            kt = iloc*nPolarity + ipol;
            ierr = parmt_polarity_computeGreensRowFromTtimes(
                                waveType[ipol], icomps[ipol],
                                evlas[iloc], evlos[iloc], deps[iloc],
                                stlas[ipol], stlos[ipol],
                                cmpincs[ipol], cmpazs[ipol],
                                parms.ttimesTablesDir,
                                parms.ttimesModel,
                                G6);
            if (ierr != 0)
            {
                printf("%s: Error computing polarities %d %d on PID %d\n",
                       fcnm, iloc, ipol, myid);
                ierrAll = ierrAll + 1;
                continue;
            }
            // Save it
            GxxBuf[kt] = G6[0];
            GyyBuf[kt] = G6[1];
            GzzBuf[kt] = G6[2];
            GxyBuf[kt] = G6[3];
            GxzBuf[kt] = G6[4];
            GyzBuf[kt] = G6[5];
            //printf("%f %f %f %e %e %e %e %e %e\n", stlas[ipol], stlos[ipol], deps[iloc], G6[0], G6[1], G6[2], G6[3], G6[4], G6[5]);
        } // Loop on the polarities 
        //printf("\n");
    } // Loop on the locations
    ierr = ierrAll;
    if (ierr != 0)
    {
        printf("%s: Error computing polarities from ttimes\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    polarityData->nPolarity = nPolarity;
    polarityData->nlocs = data.nlocs;
    polarityData->nobs = data.nobs;
    polarityData->obsMap = array_copy32i(nPolarity, observation, &ierr);
    MPI_Allreduce(GxxBuf, polarityData->Gxx, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    MPI_Allreduce(GyyBuf, polarityData->Gyy, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    MPI_Allreduce(GzzBuf, polarityData->Gzz, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    MPI_Allreduce(GxyBuf, polarityData->Gxy, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    MPI_Allreduce(GxzBuf, polarityData->Gxz, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    MPI_Allreduce(GyzBuf, polarityData->Gyz, data.nlocs*nPolarity,
                  MPI_DOUBLE, MPI_SUM, globalComm);
    // Finally assemble Green's functions into modeling matrix
WAIT:;
ERROR:;
    MPI_Barrier(globalComm);
    // free memory
    memory_free64f(&GxxBuf);
    memory_free64f(&GyyBuf);
    memory_free64f(&GzzBuf);
    memory_free64f(&GxyBuf);
    memory_free64f(&GxzBuf);
    memory_free64f(&GyzBuf);
    memory_free64f(&deps);
    memory_free64f(&evlas);
    memory_free64f(&evlos);
    memory_free64f(&stlas);
    memory_free64f(&stlos);
    memory_free64f(&cmpincs);
    memory_free64f(&cmpazs);
    memory_free32i(&observation);
    memory_free32i(&polarity);
    memory_free32i(&icomps);
    memory_free32i(&waveType);
    return ierr;
}
//============================================================================//
int parmt_polarity_computeGreensRowFromData(const struct sacData_struct data,
                                            const int wavetype,
                                            const double evdp,
                                            const char *dirnm,
                                            const char *model,
                                            double *__restrict__ G)
{
    const char *fcnm = "parmt_polarity_computeGreensRowFromData\0";
    char kcmpnm[16];
    double cmpaz, cmpinc, evla, evlo, stla, stlo;
    int icomp, ierr;
    size_t lenos;
    memset(kcmpnm, 16, 16*sizeof(char));
    ierr  = sacio_getFloatHeader(SAC_FLOAT_CMPINC, data.header, &cmpinc); 
    ierr += sacio_getFloatHeader(SAC_FLOAT_CMPAZ,  data.header, &cmpaz);
    ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA,   data.header, &evla);
    ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO,   data.header, &evlo);
    ierr += sacio_getFloatHeader(SAC_FLOAT_STLA,   data.header, &stla);
    ierr += sacio_getFloatHeader(SAC_FLOAT_STLO,   data.header, &stlo);
    ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data.header, kcmpnm);
    if (ierr != 0) 
    {
        printf("%s: Failed to get header information\n", fcnm);
        return -1;
    }
    // SAC to SEED convention
    cmpinc = cmpinc - 90.0;
    // Figure out the component
    lenos = MAX(1, strlen(kcmpnm));
    icomp = 1;
    if (kcmpnm[lenos-1] == 'Z' || kcmpnm[lenos-1] == '1')
    {
        icomp = 1;
    }
    else if (kcmpnm[lenos-1] == 'N' || kcmpnm[lenos-1] == '2')
    {
        icomp = 2;
    }
    else if (kcmpnm[lenos-1] == 'E' || kcmpnm[lenos-1] == '3')
    {
        icomp = 3;
    }
    else
    {
        printf("%s: Can't classify component %s\n", fcnm, kcmpnm);
        return -1;
    }
    ierr = parmt_polarity_computeGreensRowFromTtimes(wavetype, icomp,
                                                     evla, evlo, evdp,
                                                     stla, stlo,
                                                     cmpinc, cmpaz,
                                                     dirnm, model,
                                                     G);
    if (ierr != 0)
    {
        printf("%s: Failed to compute G\n", fcnm);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes a row of the Green's functions matrix from ttimes
 *
 * @param[in] wavetype   =1 then this is a P wave.
 *                       =2 then this is an S wave
 * @param[in] icomp      =1 then this is the vertical channel.
 *                       =2 then this is the north (2) channel.
 *                       =3 then this is the east (3) channel.
 * @param[in] evla       event latitude (degrees)
 * @param[in] evlo       event longitude (degrees)
 * @param[in] evdp       event depth (km)
 * @param[in] stla       station latitude (degrees)
 * @param[in] stlo       station longitude (degrees)
 * @param[in] cmpaz      component azimuth (0 north, +90 east)
 * @param[in] cmpinc     component inclinantion (-90 up, 0 east/north, +90 down)
 * @param[in] dirnm      directory containing the ttimes precomputed binary
 *                       files
 * @param[in] model      model (e.g. ak135 or iasp91)
 *
 * @param[out] G         row of matrix s.t. G*m produces estimates the polarity
 *                       at the station.  Here m is packed 
 *                       \f$ \{m_{xx}, m_{yy}, m_{zz},
 *                             m_{xy}, m_{xz}, m_{yz} \} \f$
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_polarity_computeGreensRowFromTtimes(
    const int wavetype, const int icomp,
    const double evla, const double evlo, const double evdp,
    const double stla, const double stlo,
    const double cmpinc, const double cmpaz,
    const char *dirnm, const char *model,
    double *__restrict__ G)
{
    const char *fcnm = "parmt_polarity_computeGreensRowFromCoordinates\0";
    double aoiRec, az, azSrc, baz, bazRec, dist, delta, toaSrc;
    struct ttimesTravelTime_struct ttime;
    int ierr;
    ierr = 0;
    if (G == NULL)
    {
        printf("%s: Error G is NULL\n", fcnm);
        return -1;
    }
    memset(G, 0, 6*sizeof(double));
    if (evdp < 0.0 || evdp > ttimes_getMaxDepth())
    {
        printf("%s: Error depth must be between [0,%f]\n", fcnm,
                ttimes_getMaxDepth()); 
        return -1;
    }
    memset(&ttime, 0, sizeof(struct ttimesTravelTime_struct));
    geodetic_gps2distanceAzimuth(evla, evlo, stla, stlo,
                                 &dist, &delta, &az, &baz);
    if (wavetype == 1)
    {
        ierr = ttimes_getFirstPPhase(delta, evdp, dirnm, model, &ttime);
    }
    else if (wavetype == 2)
    {
        ierr = ttimes_getFirstSPhase(delta, evdp, dirnm, model, &ttime); 
    }
    else
    {
        printf("%s: Invalid phase type - must be 1 (P) or 2 (S)\n", fcnm);
        return -1;
    }
    if (ierr != 0)
    {
        printf("%s: Error computing theoretical traveltime info\n", fcnm);
        return -1;
    }
    // Compute the column in the Green's function matrix
    toaSrc = ttime.toang;
    azSrc  = az;
    bazRec = baz;
    aoiRec = ttime.aoi;
    //printf("%f %f %f %f %f %f %f\n", delta, stla, stlo, toaSrc, azSrc, aoiRec, bazRec);
    ierr = parmt_polarity_computeGreensMatrixRow(wavetype, icomp,
                                                 azSrc, toaSrc, bazRec, aoiRec,
                                                 cmpinc, cmpaz, G);
   
    if (ierr != 0)
    {
        printf("%s: Failed to compute polarity greens row\n", fcnm);
        memset(G, 0, 6*sizeof(double));
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes a column for the Green's function s.t. G*m produces an 
 *        estimate of polarity on the icomp'th component measured in the 
 *        far-field.  The moment tensor which would be applied to this row
 *        is packed: 
 *          \f$ \{m_{xx}, m_{yy}, m_{zz}, m_{xy}, m_{xz}, m_{yz} \} \f$
 *        with convention North, East, Down (e.g., Jost and Herrmann).
 *        For more see Quantitative Seismology - Aki and Richards 2002,
 *        Eqn 4.96 on pg 111 and Source Mechanisms of Earthquakes: Theory
 *        and Practice - Udias et al. pg 100..
 *
 * @param[in] wavetype   =1 -> P wave
 *                       =2 -> S wave
 * @param[in] icomp      receiver component of motion:
 *                         =1 -> Vertical channel
 *                         =2 -> 1 or North channel
 *                         =3 -> 2 or East channel
 * @param[in] azSrc      source to receiver azimuth is measured positive from
 *                       north (degrees)
 * @param[in] toaSrc     take-off angle (measured positive from x3 where x3 
 *                       points down) (degrees)
 * @param[in] bazRec     receiver to source back azimuth measured positive
 *                       from north (degrees)
 * @param[in] aoiRec     angle of incidence at receiver
 * @param[in] cmpaz      component azimuth (0 north, +90 east)
 * @param[in] cmpinc     component inclinantion (-90 up, 0 east/north, +90 down)
 *
 * @param[out] G         row of matrix s.t. G*m produces estimates the polarity
 *                       at the station.  Here m is packed 
 *                       \f$ \{m_{xx}, m_{yy}, m_{zz},
 *                             m_{xy}, m_{xz}, m_{yz} \} \f$
 *
 * @bugs i've really only looked at the vertical p-teleseismic case - no idea
 *        about other phases, wavetypes
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_polarity_computeGreensMatrixRow(const int wavetype,
                                          const int icomp,
                                          const double azSrc,
                                          const double toaSrc,
                                          const double bazRec,
                                          const double aoiRec,
                                          const double cmpinc,
                                          const double cmpaz,
                                          double *__restrict__ G)
{
    const char *fcnm = "parmt_polarity_computeGreensMatrixRow\0"; 
    double M[9], G63[18], up[6], ush[6], usv[6],
           gam[3], lhat[3], phat[3], phihat[3],
           cosba, cos_cmpaz, cost_rec,
           r11, r12, r13, r21, r22, r23, r31, r32, r33,
           sinba, sin_cmpaz,
           sint_rec,
           t1, t2, t3, theta, xsign, u1, u2, ue, un, uz;
    int k;
    const double pi180 = M_PI/180.0;
    const bool lrot = true;//false;
    const int P_WAVE = 1;
    const int S_WAVE = 2;
    //------------------------------------------------------------------------//
    //
    // Error check
    if (icomp < 1 || icomp > 3)
    {
        printf("%s: Invalid component\n", fcnm);
        return -1;
    } 
    if (wavetype < P_WAVE || wavetype > S_WAVE)
    {
        printf("%s: Invalid wavetype\n", fcnm);
        return -1;
    }
    // Fill the basis at the source (Aki and Richards Eqn 4.88)
    fillBasis(toaSrc, azSrc, gam, phat, phihat);
    // Compute the contractions for u_p, u_sv, and u_sh (A&R Eqn 4.96)
    // for all 6 individual moment tensor terms
    for (k=0; k<6; k++)
    {
        // Set the moment tensor with the k'th mt term 1 others zero
        setM3x3(k, M);
        // Compute contraction
        up[k]  = computeContraction3x3(gam,    M, gam);
        usv[k] = computeContraction3x3(phat,   M, gam);
        ush[k] = computeContraction3x3(phihat, M, gam); 
    }
    // Fill the basis at the receiver - notice the basis uses the forward
    // azimuth from the receiver to the source 
    fillBasis(aoiRec, (bazRec + 180.0), lhat, phat, phihat); 
    // Compute geometric factors at receiver 
//printf("%f %f %f %f %f\n", toaSrc, azSrc, bazRec, aoiRec, cmpaz);
    theta = aoiRec*pi180;
    cosba     = cos(bazRec*pi180);
    cost_rec  = cos(theta);
    sinba     = sin(bazRec*pi180);
    sint_rec  = sin(theta);
    cos_cmpaz = cos(cmpaz*pi180);
    sin_cmpaz = sin(cmpaz*pi180);
    // Set 3D rotation matrix
    r11 = cost_rec;
    r21 = sint_rec;
    r31 = 0.0;
    r12 =-sint_rec*sinba;
    r22 = cost_rec*sinba;
    r32 =         -cosba;
    r13 =-sint_rec*cosba;
    r23 = cost_rec*cosba;
    r33 =          sinba;
    // Compute 6 x 3 subforward modeling matrix with row for up (1 channel)
    if (wavetype == P_WAVE)
    {
        // Flip sign for receivers that acquire positive down 
        xsign = 1.0;
        if (fabs(cmpinc - 90.0) < 1.e-4){xsign =-1.0;}
        // Loop on mts terms
        for (k=0; k<6; k++)
        {
            // Extract the (north, east, down) component
            t3 =-up[k]*lhat[0]; // z-down -> z-up
            //t3 = up[k]*lhat[0]; // z-up
            t2 = up[k]*lhat[1]; // east
            t1 = up[k]*lhat[2]; // north 
            // Not sure if i have to rotate LQT -> ZNE
            if (lrot)
            {
                uz = t1*r11 + t2*r21 + t3*r31; // Z
                ue = t1*r12 + t2*r22 + t3*r32; // E
                un = t1*r13 + t2*r23 + t3*r33; // N
            }
            else
            {
                uz = t3;
                ue = t2;
                un = t1;
            }
            // Rotate into (1, 2)
            u1 = un*cos_cmpaz + ue*sin_cmpaz;
            u2 =-un*sin_cmpaz + ue*cos_cmpaz;
            // Finish
            G63[k*3+0] = xsign*uz;
            //G63[k*3+0] =-xsign*uz;
            G63[k*3+1] = u1;
            G63[k*3+2] = u2;
        }
    }
    // SH wave
    else
    {
        // Loop on mts terms
        for (k=0; k<6; k++)
        {
            // Extract the (north, east, down) component
            t3 =-usv[k]*phat[0] - ush[k]*phihat[0]; // z-down -> z-up
            t2 = usv[k]*phat[1] + ush[k]*phihat[1]; // east
            t1 = usv[k]*phat[2] + ush[k]*phihat[2]; // north
            // Not sure if i have to rotate LQT -> ZNE
            if (lrot)
            {
                uz = t1*r11 + t2*r21 + t3*r31; // Z
                ue = t1*r12 + t2*r22 + t3*r32; // E
                un = t1*r13 + t2*r23 + t3*r33; // N
            }
            else
            {
                uz = t3;
                ue = t2;
                un = t1;
            }
            // Rotate into (1, 2)
            u1 = un*cos_cmpaz + ue*sin_cmpaz;
            u2 =-un*sin_cmpaz + ue*cos_cmpaz;
            // Flip sign for receivers that acquire positive down 
            xsign = 1.0;
            if (fabs(cmpinc - 90.0) < 1.e-4){xsign =-1.0;}
            // Finish
            G63[k*3+0] = xsign*uz;
            G63[k*3+1] = u1;
            G63[k*3+2] = u2;
        }
    }
    // Copy the result - vertical
    if (icomp == 1)
    {
        for (k=0; k<6; k++)
        {
            G[k] = G63[k*3+0];
        }
    }
    // 1 component
    else if (icomp == 2)
    {
        for (k=0; k<6; k++)
        {
            G[k] = G63[k*3+1];
        }
    }
    // 2 component
    else if (icomp == 3)
    {
        for (k=0; k<6; k++)
        {
            G[k] = G63[k*3+2];
        }
    }
    return 0;
}

/*!
 * @brief Fills in the basis vectors - Aki and Richards pg 108
 */
static void fillBasis(const double i, const double phi,
                      double *__restrict__ gam,
                      double *__restrict__ phat,
                      double *__restrict__ phihat)
{
    double cosi, sini, cosp, sinp;
    const double pi180 = M_PI/180.0;
    cosi = cos(i*pi180);
    sini = sin(i*pi180);
    cosp = cos(phi*pi180);
    sinp = sin(phi*pi180);
    gam[0] = sini*cosp;
    gam[1] = sini*sinp;
    gam[2] = cosi;

    phat[0] = cosi*cosp;
    phat[1] = cosi*sinp;
    phat[2] =-sini; 

    phihat[0] =-sinp;
    phihat[1] = cosp;
    phihat[2] = 0.0;
}

/*!
 * @brief Sets the moment tensor matrix where the k'th moment tensor
 *        term is 1 and others are zero.  Here moment tensor terms
 *        are counted {0,1,2,3,4,5} = {xx,yy,zz,xy,xz,yz}
 */
static void setM3x3(const int k, double *__restrict__ M)
{
    memset(M, 0, 9*sizeof(double));
    // mxx (fill mtt)
    if (k == 0)
    {
        M[4] = 1.0; //M[0][0] = 1.0;
    }
    // myy (fill mpp)
    else if (k == 1)
    {
        M[8] = 1.0; //M[1][1] = 1.0;
    }
    // mzz (fill mrr)
    else if (k == 2)
    {
        M[0] = 1.0; //M[2][2] = 1.0;
    }
    // mxy and myz (fill mtp)
    else if (k == 3)
    {
        M[5] =-1.0; //M[0][1] = 1.0;
        M[7] =-1.0; //M[1][0] = 1.0;
    }
    // mxz and mzx (fill mrp)
    else if (k == 4)
    {
        M[1] = 1.0; //M[0][2] = 1.0;
        M[3] = 1.0; //M[2][0] = 1.0;
    }
    // myz and mzy (fill mrp)
    else
    {
       M[2] =-1.0; //M[1][2] = 1.0;
       M[6] =-1.0; //M[2][1] = 1.0;
    }
/*
    // mxx
    if (k == 0)
    {
        M[0] = 1.0; //M[0][0] = 1.0;
    }
    // myy
    else if (k == 1)
    {
        M[4] = 1.0; //M[1][1] = 1.0;
    }
    // mzz
    else if (k == 2)
    {
        M[8] = 1.0; //M[2][2] = 1.0;
    }
    // mxy and myz
    else if (k == 3)
    {
        M[1] = 1.0; //M[0][1] = 1.0;
        M[3] = 1.0; //M[1][0] = 1.0;
    }
    // mxz and mzx
    else if (k == 4)
    {
        M[2] = 1.0; //M[0][2] = 1.0;
        M[6] = 1.0; //M[2][0] = 1.0;
    }
    // myz and mzy
    else
    {
       M[5] = 1.0; //M[1][2] = 1.0;
       M[7] = 1.0; //M[2][1] = 1.0;
    }
*/
    return;
}
/*!
 * @brief Computes the contraction in Aki and Richards Eqn 4.96
 *        for the given basis vectors and moment tensor.
 */
static double computeContraction3x3(const double *__restrict__ a,
                                    const double *__restrict__ M,
                                    const double *__restrict__ b)
{ 
    double res;
    int i, ij, j;
    res = 0.0;
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++) 
        {
            ij = 3*i + j;
            res = res + a[i]*M[ij]*b[j];
        }
    }
    return res;
}

int polarity_performLocationSearch64f(const MPI_Comm locComm,
                                      const int blockSize,
                                      struct polarityData_struct polarityData, 
                                      struct localMT_struct mtloc,
                                      double *__restrict__ phi)
{
    const char *fcnm = "polarity_performLocationSearch64f\0"; 
    double *Dmat, *G, *Sigma, *phiLoc, *phiWork;
    int icol, ierr, ierrAll, iloc, imtLoc, ipol, jloc, jndx,
        kt, mylocID, nlocProcs, mblock, pad;
    const int nPolarity = polarityData.nPolarity;
    const int Mrows = nPolarity;
    const int Kcols = 6;
    const int nlocs = polarityData.nlocs;
    const int ldm = mtloc.ldm;
    const int nmt = mtloc.nmt;
    const int master = 0;
    // Initialize
    ierr = 0;
    ierrAll = 0;
    MPI_Comm_size(locComm, &nlocProcs);
    MPI_Comm_rank(locComm, &mylocID);
    // Compute padding for 64 bit data alignment in observations and mt blocks
    pad = computePadding64f(blockSize);
    mblock = blockSize + pad;
    // Set space
    phiLoc = memory_calloc64f(nmt);
    G = memory_calloc64f(nPolarity*LDG);
    Dmat = memory_calloc64f(nPolarity*mblock);
    Sigma = array_set64f(nPolarity, 1.0, &ierr); // Default to identity
    phiWork = memory_calloc64f(nlocs*mtloc.nmtAll);
    // Set the row major data matrix where each row is an observation
    for (ipol=0; ipol<nPolarity; ipol++) 
    {
        for (icol=0; icol<mblock; icol++)
        {
            Dmat[ipol*mblock+icol] = polarityData.polarity[ipol];
            if (icol == 1){printf("%f\n", polarityData.polarity[ipol]);}
        }
    }
    // Input weights are given - copy them to diagonal weight matrix
    if (polarityData.wts != NULL)
    {
        ierr = array_copy64f_work(nPolarity, polarityData.wts, Sigma);
    }
    // Loop on the source locations in the grid search
    for (jloc=0; jloc<nlocs; jloc=jloc+nlocProcs)
    {
        iloc = jloc + mylocID;
        kt = iloc*nPolarity;
        imtLoc = iloc*nmt;
        // Assemble the row major Green's functions matrix
        array_zeros64f_work(LDG*nPolarity, G);
        #pragma omp simd
        for (ipol=0; ipol<nPolarity; ipol++)
        {
            G[LDG*ipol+0] = polarityData.Gxx[kt+ipol];
            G[LDG*ipol+1] = polarityData.Gyy[kt+ipol];
            G[LDG*ipol+2] = polarityData.Gzz[kt+ipol];
            G[LDG*ipol+3] = polarityData.Gxy[kt+ipol];
            G[LDG*ipol+4] = polarityData.Gxz[kt+ipol];
            G[LDG*ipol+5] = polarityData.Gyz[kt+ipol];
        }
        // Perform the polarity search for all mt's at this location
        ierr = performPolaritySearch64f(nmt, ldm,
                                        nPolarity,
                                        blockSize, mblock,
                                        Mrows, Kcols,
                                        Dmat, G, Sigma, mtloc.mts,
                                        phiLoc);
        if (ierr != 0)
        {
            printf("%s: Error in mt polarity search\n", fcnm);
            ierrAll = ierrAll + 1;
        }
        // Gather the moment tensor search results onto the master 
        jndx = 0;
        if (mtloc.myid == master){jndx = iloc*mtloc.nmtAll;}
        if (mtloc.commSize == 1)
        {
            array_copy64f_work(mtloc.nmt, phiLoc, &phiWork[jndx]);
        }
        else
        {
            MPI_Gatherv(phiLoc, mtloc.nmt, MPI_DOUBLE,
                        &phiWork[jndx], mtloc.nmtProc, mtloc.offset,
                        MPI_DOUBLE, master, mtloc.comm);
        }
printf("%f %f\n", array_min64f(nmt, phiLoc, &ierr), array_max64f(nmt, phiLoc, &ierr));

    }
    // Reduce the search onto the master
    //if (linLocComm && mylocID == master)
    {
         MPI_Reduce(phiWork, phi, nlocs*mtloc.nmtAll, MPI_DOUBLE,
                    MPI_SUM, master, locComm);
    }
    // Free memory
    memory_free64f(&phiWork);
    memory_free64f(&phiLoc);
    memory_free64f(&G);
    memory_free64f(&Dmat);
    memory_free64f(&Sigma);
    return 0;
}

//============================================================================//
static int performPolaritySearch64f(const int nmt, const int ldm, 
                                    const int nPolarity,
                                    const int blockSize, const int mblock,
                                    const int Mrows, const int Kcols,
                                    const double *__restrict__ Dmat,
                                    const double *__restrict__ G,
                                    const double *__restrict__ Sigma,
                                    const double *__restrict__ mts, 
                                    double *__restrict__ phi) 
{
    const char *fcnm = "performPolaritySearch64f\0";
    double *U, *Usign, *res2, *sqrtSigmaWt, res, traceSigma;
    int i, ic, ierr, ipol, jmt, kmt, nmtBlocks, Ncols;
    const double one = 1.0; 
    const double zero = 0.0; 
    ierr = 0; 
    nmtBlocks = (int) ((double) (nmt)/(double) (blockSize) + 1.0);
    if (nmtBlocks*blockSize < nmt) 
    {    
        printf("%s: Internal error - all mts wont be visited\n", fcnm);
        return -1;
    }
#ifdef __INTEL_COMPILER
    __assume_aligned(G, 64);
    __assume_aligned(Dmat, 64);
#endif
    // The VR in Chiang 2016 is:
    //   (1 - \frac{ \sum_i w_i (Pol_{i,obs} - Pol_{i,est})^2 }
    //             { \sum_i w_i (Pol_{i,obs}^2) } )
    // Because Pol_{obs} and Pol_{est} take values of + or -1  we can
    // reduce the demoniator so that the VR is
    //   (1 - \frac{ \sum_i w_i (Pol_{i,obs} - Pol_{i,est})^2 }
    //             { \sum_i w_i })
    // Furthermore, we can incorporate this term into the weighting
    // function s.t.  Sigma/trace(Sigma) to obtain
    //   (1 - \frac{ \sum_i \hat{Sigma} (Pol_{i,obs} - Pol_{i,est})^2)
    // Finally, because we want likelihoods we want to rescale the
    // residual squared from [0,4] to [0,1].  To do this we multiply
    // the weight by 0.5 s.t. the residual [-2,2] goes to [-1,1] which,
    // when squared, fits in the desired bounds of [0,1]; i.e. multiply by
    // (1/2)^2 = 0.25
    traceSigma = array_sum64f(nPolarity, Sigma, &ierr);
    if (ierr != 0){traceSigma = 1.0;}
    sqrtSigmaWt = memory_calloc64f(nPolarity);
#ifdef __INTEL_COMPILER
    __assume_aligned(sqrtSigmaWt, ISCL_MEMORY_ALIGN);
#endif
    for (i=0; i<nPolarity; i++){sqrtSigmaWt[i] = sqrt(0.25*Sigma[i]/traceSigma);}
    #pragma omp parallel \
     private (i, ic, ipol, jmt, kmt, Ncols, res, res2, U, Usign) \
     shared (G, Dmat, mts, nmtBlocks, phi, sqrtSigmaWt) \
     default (none) reduction(+:ierr)
    {
    U = memory_calloc64f(mblock*nPolarity);
    res2 = memory_calloc64f(mblock);
    Usign = memory_calloc64f(mblock*nPolarity);
#ifdef __INTEL_COMPILER
    __assume_aligned(res2, ISCL_MEMORY_ALIGN);
    __assume_aligned(Usign, ISCL_MEMORY_ALIGN);
#endif
    #pragma omp for
    for (kmt=0; kmt<nmtBlocks; kmt++)
    {
        jmt = kmt*blockSize;
        Ncols = MIN(blockSize, nmt - jmt); // Number of columns of M
        // Compute U = GM
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    Mrows, Ncols, Kcols, one, G, LDG,
                    &mts[ldm*jmt], ldm, zero, U, mblock);
        // Make the theoretical polarity +1 or -1 to match the data
        array_copysign64f_work(mblock*nPolarity, U, Usign);
        memset(res2, 0, (size_t) mblock*sizeof(double));
        // Compute the weighted residual
        for (ipol=0; ipol<nPolarity; ipol++)
        {
            #pragma omp simd aligned(Dmat, Usign, res2: ISCL_MEMORY_ALIGN)
            for (ic=0; ic<Ncols; ic++)
            {
                // Change residual range from [0,2] to [0,1] and weight
                res = sqrtSigmaWt[ipol]*( Dmat[ipol*mblock+ic]
                                        - Usign[ipol*mblock+ic]);
                res2[ic] = res2[ic] + res*res;
            }
        }
        // Compute the variance reduction and put it into objective function
        for (ic=0; ic<Ncols; ic++)
        {
            phi[jmt+ic] = 1.0 - res2[ic];
        }
    } // Loop on moment tensor blocks
    memory_free64f(&Usign);
    memory_free64f(&U);
    memory_free64f(&res2);
    } // end parallel section
    memory_free64f(&sqrtSigmaWt);
    return 0;
}

static int computePadding64f(const int n)
{
    size_t mod, pad;
    int ipad;
    // Set space and make G matrix
    pad = 0;
    mod = ((size_t) n*sizeof(double))%64;
    if (mod != 0)
    {
        pad = (64 - mod)/sizeof(double);
    }
    ipad = (int) pad;
    return ipad;
}
