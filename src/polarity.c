#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
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

int parmt_polarity_gridSearch(const struct parmtData_struct data)
{
    const char *fcnm = "parmt_polarity_gridSearch\0";
    char kt0[8];
    double *deps, *G, *Gw;
    int *observation, *polarity, *waveType, ierr, iloc, iobs, ipol,
        iwav, k, nPolarity;
    G = NULL; 
    // Extract the depths in the grid-search from the first waveform 
    iobs = 0;
    deps = memory_calloc64f(data.nlocs);
    for (iloc=0; iloc<data.nlocs; iloc++)
    {
        k = iobs*data.nlocs + iloc;
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVDP, data.sacGxx[k].header,
                                    &deps[iloc]);
        if (ierr != 0)
        {
            printf("%s: Unable to get event depth\n", fcnm);
            memory_free64f(&deps);
            return -1;
        }
    }
    // Compute the number of polarities
    nPolarity = 0;
    polarity = memory_calloc32i(data.nobs);
    waveType = memory_calloc32i(data.nobs);
    observation = array_set32i(data.nobs, -1, &ierr);
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        ierr = sacio_getIntegerHeader(SAC_INT_UNUSED0,
                                      data.data[iobs].header,
                                      &ipol);
        // Ensure it is defined
        if (ierr != 0)
        {
            ierr = sacio_getCharacterHeader(SAC_CHAR_KT0,
                                            data.data[iobs].header,
                                            kt0);
            if (ierr != 0)
            {
                printf("%s: t0 phase name is not set", fcnm);
                continue;
            }
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
                printf("%s: Could not classify t0 phase %s\n", fcnm, kt0);
                continue;
            }
            if (ipol == 1 || ipol ==-1)
            {
                observation[nPolarity] = iobs;
                waveType[nPolarity] = iwav;
                polarity[nPolarity] = ipol;
                nPolarity = nPolarity + 1;
            }
            else
            {
                printf("%s: Can't determine polarity\n", fcnm);
            }
/*
            if ((iwav == 1 || iwav == 2) && (ipol == 1 || ipol == 2))
            {
                // Loop on the event depths
                parmt_polarity_computeGreensRowFromData(data.data[iobs],
                                                        iwav,
                                                        deps[8],
                                                        parms.ttimesTablesDir,
                                                        parms.ttimesModel,
                                                        &G[LDG]);
            }
*/
        }
    }
    if (nPolarity == 0)
    {
        printf("%s: There are no polarities\n", fcnm);
        memory_free32i(&polarity);
        memory_free32i(&waveType);
        return 0; 
    }
    // Compute the forward modeling matrix
    for (iobs=0; iobs<data.nobs; iobs++)
    {
    //    parmt_polarity_computeGreensRowFromData(data.data[iobs],

    }
    // 
    memory_free64f(&deps);
    memory_free64f(&G);
    memory_free32i(&observation);
    memory_free32i(&polarity);
    memory_free32i(&waveType);
    return 0;
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
    ierr = parmt_polarity_computeGreensRowFromCoordinates(wavetype, icomp,
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
int parmt_polarity_computeGreensRowFromCoordinates(
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
    int i, k;
    const double pi180 = M_PI/180.0;
    const bool lrot = false;
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
    fillBasis(aoiRec, (bazRec - 180.0), lhat, phat, phihat); 
    // Compute geometric factors at receiver 
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
        // Loop on mts terms
        for (i=0; i<6; i++)
        {
            // Extract the (north, east, down) component
            t3 =-up[i]*lhat[0]; // z-down -> z-up
            t2 = up[i]*lhat[1]; // east
            t1 = up[i]*lhat[2]; // north 
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
            G63[i*3+0] = xsign*uz;
            G63[i*3+1] = u1;
            G63[i*3+2] = u2;
        }
    }
    // SH wave
    else
    {
        // Loop on mts terms
        for (i=0; i<6; i++)
        {
            // Extract the (north, east, down) component
            t3 =-usv[i]*phat[0] - ush[i]*phihat[0]; // z-down -> z-up
            t2 = usv[i]*phat[1] + ush[i]*phihat[1]; // east
            t1 = usv[i]*phat[2] + ush[i]*phihat[2]; // north
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
            G63[i*3+0] = xsign*uz;
            G63[i*3+1] = u1;
            G63[i*3+2] = u2;
        }
    }
    // Copy the result - vertical
    if (icomp == 1)
    {
        for (i=0; i<6; i++)
        {
            G[i] = G63[i*3+0];
        }
    }
    // 1 component
    else if (icomp == 2)
    {
        for (i=0; i<6; i++)
        {
            G[i] = G63[i*3+1];
        }
    }
    // 2 component
    else if (icomp == 3)
    {
        for (i=0; i<6; i++)
        {
            G[i] = G63[i*3+2];
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

int polarity_performSearch(const int ldm, const int nmt,
                           const int nobs, const int blockSize,
                           const double *__restrict__ Gxx,
                           const double *__restrict__ Gyy,
                           const double *__restrict__ Gzz,
                           const double *__restrict__ Gxy,
                           const double *__restrict__ Gxz,
                           const double *__restrict__ Gyz, 
                           const double *__restrict__ mts,
                           const double *__restrict__ d,
                           const double *__restrict__ wts,
                           double *__restrict__ phi)
{
    const char *fcnm = "polarity_performSearch\0";
    double *Dmat, *G, *Sigma;
    int i, ic, ierr, pad, mblock;
    const int Mrows = nobs;
    const int Kcols = 6;
    // Compute the blocksize
    pad = computePadding64f(blockSize);
    mblock = blockSize + pad;
    // If the weights aren't defined then set them to 1
    if (wts == NULL)
    {
        Sigma = array_set64f(nobs, 1.0, &ierr);
    }
    else
    {
        Sigma = array_copy64f(nobs, wts, &ierr);
    }
    // Set the data
    Dmat = memory_calloc64f(nobs*mblock);
    for (i=0; i<nobs; i++)
    {
        for (ic=0; ic<mblock; ic++)
        {
            Dmat[i*mblock+ic] = d[i];
        }
    }
    // Assemble the Green's functions matrix
    G = memory_calloc64f(nobs*LDG);
    for (i=0; i<nobs; i++)
    {
        G[LDG*i+0] = Gxx[i];
        G[LDG*i+1] = Gyy[i];
        G[LDG*i+2] = Gzz[i];
        G[LDG*i+3] = Gxy[i];
        G[LDG*i+4] = Gxz[i];
        G[LDG*i+5] = Gyz[i];
    }
    // Perform the polarity search
    ierr = performPolaritySearch64f(nmt, ldm,
                                    nobs,
                                    blockSize, mblock,
                                    Mrows, Kcols,
                                    Dmat, G, Sigma, mts, phi);
    if (ierr != 0)
    {
        printf("%s: Error performing polarity search\n", fcnm);
    }
    // Clean up
    memory_free64f(&G);
    memory_free64f(&Dmat);
    memory_free64f(&Sigma);
    return ierr;
}
//============================================================================//
static int performPolaritySearch64f(const int nmt, const int ldm, 
                                    const int nobs,
                                    const int blockSize, const int mblock,
                                    const int Mrows, const int Kcols,
                                    const double *__restrict__ Dmat,
                                    const double *__restrict__ G,
                                    const double *__restrict__ Sigma,
                                    const double *__restrict__ mts, 
                                    double *__restrict__ phi) 
{
    const char *fcnm = "performPolaritySearch64f\0";
    double *M, *U, *Usign, *res2, *sigmaWt, res, traceSigma;
    int i, ic, idx, imt, ierr, jdx, jmt, kmt, nmtBlocks, Ncols;
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
    __assume_aligned(Sigma, 64);
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
    traceSigma = array_sum64f(nobs, Sigma, &ierr);
    if (ierr != 0){traceSigma = 1.0;}
    sigmaWt = memory_calloc64f(nobs);
#ifdef __INTEL_COMPILER
    __assume_aligned(sigmaWt, 64);
#endif
    for (i=0; i<nobs; i++){sigmaWt[i] = Sigma[i]/traceSigma;}
/*
    #pragma omp parallel \
     private (i, ic, idx, imt, jdx, jmt, M, Ncols, res, res2, U, Usign) \
     shared (G, Dmat, mts, nmtBlocks, phi, sigmaWt) \
     default (none) reduction(+:ierr)
    {
*/
    M = memory_calloc64f(mblock*nobs);
    U = memory_calloc64f(mblock*nobs);
    res2 = memory_calloc64f(mblock);
    Usign = memory_calloc64f(mblock*nobs);
#ifdef __INTEL_COMPILER
    __assume_aligned(M, 64); 
    __assume_aligned(U, 64); 
    __assume_aligned(res2, 64);
    __assume_aligned(Usign, 64);
#endif
/*
    #pragma omp for
*/
    for (kmt=0; kmt<nmtBlocks; kmt++)
    {
        jmt = kmt*blockSize;
        Ncols = MIN(blockSize, nmt - jmt); // Number of columns of M
        // Set the moment tensor matrix 
        for (i=0; i<6; i++)
        {
            #pragma omp simd aligned(M: 64)
            for (ic=0; ic<Ncols; ic++)
            {
                imt = jmt + ic;
                idx = ldm*imt + i;
                jdx = mblock*i;
                M[jdx+ic] = mts[idx];
            }
        }
        // Compute U = GM
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    Mrows, Ncols, Kcols, one, G, LDG,
                    M, mblock, zero, U, mblock);
        // Make the theoretical polarity +1 or -1 to match the data
        array_copysign64f_work(mblock*nobs, U, Usign);
        memset(res2, 0, (size_t) mblock*sizeof(double));
        // Compute the weighted residual
        for (i=0; i<nobs; i++)
        {
            for (ic=0; ic<Ncols; ic++)
            {
                res = Dmat[i*mblock+ic] - Usign[i*mblock+ic];
                res2[ic] = res2[ic] + sigmaWt[i]*(res*res);
            }
        }
        // Compute the variance and put it into objective function
        for (ic=0; ic<Ncols; ic++)
        {
            phi[jmt+ic] = 1.0 - res2[ic];
        }
    }
    memory_free64f(&Usign);
    memory_free64f(&U);
    memory_free64f(&M);
    memory_free64f(&res2);
/*
    } // end parallel section
*/
    memory_free64f(&sigmaWt);
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
