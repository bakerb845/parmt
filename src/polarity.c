#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "parmt_polarity.h"
#include "ttimes.h"
#include "iscl/geodetic/geodetic.h"

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
