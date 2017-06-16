#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

static double dot3(const double *__restrict__ a,
                   const double *__restrict__ b);
static void setM3x3(const int k, double *__restrict__ M);
static void fillBasis(const double i, const double phi,
                      double *__restrict__ gamHat,
                      double *__restrict__ tHat,
                      double *__restrict__ phiHat);
static double computeContraction3x3(const double *__restrict__ a,
                                    const double *__restrict__ M,
                                    const double *__restrict__ b);

int main()
{
    double M[9], gamHat[3], phi[3], theta[3], c, xsum;
    int i, k;
    double toa = 22.0;
    double az = 99.0;
    // Set the basis
    fillBasis(toa, az, gamHat, theta, phi);
    // Check normalization; p1
    for (k=0; k<6; k++)
    {
       setM3x3(k, M);
       xsum = 0.0;
       for (i=0; i<9; i++)
       {
           xsum = xsum + M[i]*M[i];
       }
       xsum = sqrt(xsum)/sqrt(2.0);
       c = computeContraction3x3(gamHat, M, gamHat);
       //printf("%f %f %f\n", c, xsum, sqrt(2.0)*gamHat[2]*gam[2]); 
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
 * @param[in] wavetype  Observation type: a P or S wave.
 *                        =1 -> P wave. \n
 *                        =2 -> S wave.
 * @param[in] icomp     Receiver component of motion: \n
 *                        =1 -> Vertical channel \n
 *                        =2 -> 1 or North channel \n
 *                        =3 -> 2 or East channel
 * @param[in] azSrc     Source to receiver azimuth is measured positive from
 *                      north (degrees).
 * @param[in] toaSrc    Take-off angle (measured positive from x3 where x3 
 *                      points down) (degrees).
 * @param[in] bazRec    Receiver to source back azimuth measured positive
 *                      from north (degrees).
 * @param[in] aoiRec    Angle of incidence (degrees) at receiver.
 * @param[in] cmpaz     Component azimuth (0 north, +90 east).
 * @param[in] cmpinc    Component inclinantion (-90 up, 0 east/north, +90 down).
 *
 * @param[out] G        Row of matrix s.t. G*m produces estimates the polarity
 *                      at the station.  Here m is packed 
 *                      \f$ \{m_{xx}, m_{yy}, m_{zz},
 *                            m_{xy}, m_{xz}, m_{yz} \} \f$.
 *
 * @result 0 indicates success.
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
    double M[9], G63[18], G63ned[18], up[6], ush[6], usv[6],
           gamHat[3], lhat[3], tHat[3], phiHat[3],
           azRec, cosba, cos_cmpaz, cost_rec,
           r11, r12, r13, r21, r22, r23, r31, r32, r33, 
           sinba, sin_cmpaz,
           sint_rec,
           tl, tq, tt, theta, xsign, u1, u2, ue, un, uz;
    int i, ierr, k;
    const double pi180 = M_PI/180.0;
    const bool lrot = true;
    const int P_WAVE = 1; 
    const int S_WAVE = 2; 
    //------------------------------------------------------------------------//
    //
    // Error check
    ierr = 0;
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
    // Fill the basis at the source (Udias pg 78)
    fillBasis(toaSrc, azSrc, gamHat, tHat, phiHat);
    // Compute the dot-products for u_p, u_sv, and u_sh (A&R Eqn 4.97)
    // for all 6 individual moment tensor terms at the source location
    for (k=0; k<6; k++)
    {    
        // Set the moment tensor with the k'th MT term in USE coordinates.
        setM3x3(k, M);
        // Compute contraction
        up[k]  = computeContraction3x3(gamHat, M, gamHat);
        usv[k] = computeContraction3x3(tHat,   M, gamHat);
        ush[k] = computeContraction3x3(phiHat, M, gamHat);
    }
    // Fill the basis at the receiver.  Here we flip the back-azimuth 
    // (i.e. the receiver to source angle) so tHat it points in the direction
    // of the source to the receiver.  This basis should point in the LQT 
    // direction.
    azRec = bazRec + 180.0;
    if (azRec >= 360.0){azRec = azRec - 360.0;}
    fillBasis(aoiRec, azRec, lhat, tHat, phiHat);
    // Compute geometric factors at receiver 
    //printf("%f %f %f %f %f\n", toaSrc, azSrc, bazRec, aoiRec, cmpaz);
    theta = aoiRec*pi180;
    cosba     = cos(bazRec*pi180);
    cost_rec  = cos(aoiRec*pi180); //theta);
    sinba     = sin(bazRec*pi180);
    sint_rec  = sin(aoiRec*pi180); //theta);
    cos_cmpaz = cos(cmpaz*pi180);
    sin_cmpaz = sin(cmpaz*pi180);
    // Set 3D rotation matrix to go from right-handed ray coordinate
    // system to left handed Z, N, E observation system
    r11 = cost_rec;
    r21 = sint_rec;
    r31 = 0.0;
    r12 =-sint_rec*sinba;
    r22 = cost_rec*sinba;
    r32 =         -cosba;
    r13 =-sint_rec*cosba;
    r23 = cost_rec*cosba;
    r33 =          sinba;
    // Flip sign for receivers tHat acquire positive down 
    xsign = 1.0;
    if (fabs(cmpinc - 90.0) < 1.e-4){xsign =-1.0;}
    // Initalize USE 6 x 3 Greens functions 
    memset(G63, 0, 18*sizeof(double));
    // Compute 6 x 3 subforward modeling matrix with row for up (1 channel)
    if (wavetype == P_WAVE)
    {
        // Loop on mt terms
        for (k=0; k<6; k++)
        {
            // Set the vector in Eqn 4.96
            tl = up[k]*lhat[0]; // Radial (L - longitudinal)
            tq = up[k]*lhat[1]; // SV (Q - polarized with L)
            tt = up[k]*lhat[2]; // SH (T - transverse)
            // Rotate into observation coordinates (Z, E, N) where Z is up
            uz = tl*r11 + tq*r21 + tt*r31; // L -> Z
            ue = tl*r12 + tq*r22 + tt*r32; // Q -> E
            un = tl*r13 + tq*r23 + tt*r33; // T -> N
            // Rotate into (1, 2)
            u1 = un*cos_cmpaz + ue*sin_cmpaz;
            u2 =-un*sin_cmpaz + ue*cos_cmpaz;
            // Tabulate the USE MT Green's functions in observation Z12
            G63[0*6+k] = xsign*uz; // Beware the upside down sensor 
            G63[1*6+k] = u1;
            G63[2*6+k] = u2;
        }
    }
    // S wave is superposition of SV and SH wave
    else if (wavetype == S_WAVE)
    {
        // Loop on mt terms
        for (k=0; k<6; k++)
        {
            for (i=0; i<2; i++)
            {
                // Set the vector in Eqn 4.96
                if (i == 0)
                {
                    tl = up[k]*tHat[0]; // Radial (L - longitudinal)
                    tq = up[k]*tHat[1]; // SV (Q - polarized with L)
                    tt = up[k]*tHat[2]; // SH (T - transverse)
                }
                else
                {
                    tl = up[k]*phiHat[0]; // Radial (L - longitudinal)
                    tq = up[k]*phiHat[1]; // SV (Q - polarized with L)
                    tt = up[k]*phiHat[2]; // SH (T - transverse)
                }
                // Rotate into observation coordinates (Z, E, N) where Z is up
                uz = tl*r11 + tq*r21 + tt*r31; // L -> Z
                ue = tl*r12 + tq*r22 + tt*r32; // Q -> E
                un = tl*r13 + tq*r23 + tt*r33; // T -> N
                // Rotate into (1, 2)
                u1 = un*cos_cmpaz + ue*sin_cmpaz;
                u2 =-un*sin_cmpaz + ue*cos_cmpaz;
                // Tabulate the USE MT Green's functions in observation Z12
                G63[0*6+k] += xsign*uz; // Beware the upside down sensor 
                G63[1*6+k] += u1;
                G63[2*6+k] += u2;
            }
        }
    }
    else
    {
         printf("%s: No idea what the wave type %d is\n", fcnm, wavetype);
         ierr = 1;
    }
    // Convert the basis from USE to NED; i.e. instead of inputting: 
    // mrr, mtt, mpp, mrt, mrp, mtp give mxx, myy, mzz, mxy, mxz, myz
    for (i=0; i<3; i++)
    {
        G63ned[i*6+0] = G63[i*6+1]; //  tt -> xx
        G63ned[i*6+1] = G63[i*6+2]; //  pp -> yy
        G63ned[i*6+2] = G63[i*6+0]; //  rr -> zz
        G63ned[i*6+3] =-G63[i*6+5]; // -tp -> xy
        G63ned[i*6+4] = G63[i*6+3]; //  rt -> xz
        G63ned[i*6+5] =-G63[i*6+4]; // -rp -> yz
    }
    // Copy the result - vertical
    if (icomp == 1)
    {
        for (k=0; k<6; k++){G[k] = G63[k];}
    }
    // 1 component
    else if (icomp == 2)
    {
        for (k=0; k<6; k++){G[k] = G63[6+k];}
    }
    // 2 component
    else if (icomp == 3)
    {
        for (k=0; k<6; k++){G[k] = G63[12+k];}
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the radiation patterns for double-couples specified in
 *        terms of strike, dip, and rake.  This is mainly for debugging and
 *        of little use to the general application.  The requisite equations
 *        are detailed in Aki and Richards 4.89-4.91 pgs 108-109 and Udias
 *        pg 78-79.
 *
 * @param[in] strike     Fault strike angle (degrees). This is measured
 *                       positive clockwise from north.
 * @param[in] dip        Fault dip angle (degrees).  This is measured positive
 *                       down from horizontal.
 * @param[in] rake       Fault rake angle (degrees).
 * @param[in] toa        Source take-off angle (degrees).  This is measured
 *                       positive from z+ down.
 * @param[in] phi        Source to receiver azimuth (degrees).  This is measured
 *                       positive clockwise from north.
 *
 * @param[out] frp       Fault radiation pattern for P, SV, and SH waves
 *                       respectively.  This is an array of dimension [3].
 *                       Note, these do not involve the average slip or any
 *                       material properties at the source.
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
void parmt_polarity_computeAki488to491(
    const double strike, const double dip, const double rake,
    const double toa, const double phi,
    double *__restrict__ frp)
{
    double gamHat[3], nu[3],  uSlip[3], phiHat[3], tHat[3];
    double cosDelta, cosDelta_sinLambda,
           cosLambda, cosPhis, delta, gamHat_nu, gamHat_uSlip, lambda, phis,
           sinDelta, sinLambda, sinPhis;
    const double pi180 = M_PI/180.0;
    phis   = pi180*strike;
    lambda = pi180*rake; 
    delta  = pi180*dip;

    cosLambda = cos(lambda);
    cosPhis   = cos(phis);
    cosDelta  = cos(delta);
    sinLambda = sin(lambda);
    sinPhis   = sin(phis);
    sinDelta  = sin(delta); 
    cosDelta_sinLambda = cosDelta*sinLambda; 
    // Slip vector in (North, East, Down)
    uSlip[0] = cosLambda*cosPhis + cosDelta_sinLambda*sinPhis;
    uSlip[1] = cosLambda*sinPhis - cosDelta_sinLambda*cosPhis;
    uSlip[2] =-sinLambda*sinDelta;
    // Fault normal in (North, East, Down)
    nu[0] =-sinDelta*sinPhis;
    nu[1] = sinDelta*cosPhis;
    nu[2] =-cosDelta;
    // P wave direction; SV wave direction; SH wave direction in USE.
    fillBasis(toa, phi, gamHat, tHat, phiHat);
    // Compute the radiation patterns for P, SV, and SH (respectively).
    // Aki and Richards 4.89-4.91 make use of angle addition/subtraction
    // identities to show the compatibility of the two coordinate systems
    // but the dot product requires computation of many fewer trigonometric
    // terms.
    gamHat_nu = dot3(gamHat, nu);
    gamHat_uSlip = dot3(gamHat, uSlip);
    frp[0] = 2.0*gamHat_nu*gamHat_uSlip;
    frp[1] = gamHat_nu*dot3(uSlip, tHat)   + gamHat_uSlip*dot3(nu, gamHat);
    frp[2] = gamHat_nu*dot3(uSlip, phiHat) + gamHat_uSlip*dot3(nu, phiHat);
    return;
}
//============================================================================//
static double dot3(const double *__restrict__ a,
                   const double *__restrict__ b)
{
    double dp;
    dp = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return dp;
}
//============================================================================//
/*!
 * @brief Fills in the basis vectors - Aki and Richards pg 98 or pg 108 Eqn 4.88
 *        or Udias page 78-79.  This is
 *        oriented in \$ \{ \mathbf{r}, \mathbf{\theta}, \mathbf{\phi} \} \$
 *        i.e. the radial then the two transverse directions.
 *
 * @param[in] i           Incidence angle (degrees).  This is measured from
 *                        \f$ x_3 \$.  In ray coordinates this can be taken
 *                        as the take-off angle measured from the downward
 *                        vertical. 
 * @param[in] phi         Azimuth (degrees) measured from \f$ x_1 \f$.
 *                        Aki and Richards take \f$ x_1 \f$ to be in the
 *                        direction of the slip.  In ray coordinates this
 *                        can be taken as the the azimuth measured positive
 *                        from north which means tHat \f$ x_1 \f$ is oriented
 *                        positive north.
 *
 * @param[out] gamHat     Unit vector in \f$ \hat{\textbf{r}} \f$ or 
 *                        (alternatively written as \f$ \hat{\mathbf{\gamma}} \f$).
 *                        In the context of 4.88 this is the P-wave direction.
 * @param[out] tHat       Unit vector in \f$ \hat{\mathbf{\theta}} \f$.
 *                        In the context of 4.88 this is the SV-wave direction.
 * @param[out] phiHat     Unit vector in \f$ \hat{\mathbf{\phi}} \f$.
 *                        In the context of 4.88 this is the SH-wave direction.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
static void fillBasis(const double i, const double phi, 
                      double *__restrict__ gamHat,
                      double *__restrict__ tHat,
                      double *__restrict__ phiHat)
{
    double cosi, sini, cosp, sinp;
    const double pi180 = M_PI/180.0;
    cosi = cos(i*pi180);
    sini = sin(i*pi180);
    cosp = cos(phi*pi180);
    sinp = sin(phi*pi180);
    // P
    gamHat[0] = sini*cosp;
    gamHat[1] = sini*sinp;
    gamHat[2] = cosi;
    // SV 
    tHat[0] = cosi*cosp;
    tHat[1] = cosi*sinp;
    tHat[2] =-sini; 
    // SH
    phiHat[0] =-sinp;
    phiHat[1] = cosp;
    phiHat[2] = 0.0; 
    return;
}
//============================================================================//
/*!
 * @brief Computes the contraction in Aki and Richards Eqn 4.96
 *        for the given basis vectors and moment tensor:
 *        \f$ C = a_i M_{ij} b_j \f$. 
 *
 * @param[in] a    Vector of dimension of [3].
 * @param[in] M    Moment tensor with dimension [3 \times 3] in row major
 *                 order.
 * @param[in] b    Vector of dimension of [3].
 *
 * @result The contraction \f$c = a_i M_{ij} b_j \f$.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
static double computeContraction3x3(const double *__restrict__ a,
                                    const double *__restrict__ M,
                                    const double *__restrict__ b)
{ 
    double c; 
    c = a[0]*(M[0]*b[0] + M[1]*b[1] + M[2]*b[2])
      + a[1]*(M[3]*b[0] + M[4]*b[1] + M[5]*b[2])
      + a[2]*(M[6]*b[0] + M[7]*b[1] + M[8]*b[2]);
    return c;
}
//============================================================================//
/*!
 * @brief Sets the moment tensor matrix where the k'th moment tensor
 *        term is 1 and others are zero.  Here moment tensor is in USE 
 *        coordinates. 
 *
 * @param[in] k     Extraction of the k'th index. \n:
 *                  k=0 extracts \f$ M_{rr} \f$.
 *                  k=1 extracts \f$ M_{\theta \theta} \f$. 
 *                  k=2 extracts \f$ M_{\phi \phi} \f$.
 *                  k=3 extracts \f$ M_{r \theta} \f$.
 *                  k=4 extracts \f$ M_{r \phi} \f$.
 *                  k=5 Extracts \f$ M_{\theta \phi} \f$.
 *
 * @param[out] M    Appropriate moment tensor terms for computing contraction.
 *                  The moment tensor matrix has length 1.  This is an array
 *                  of dimension [9].
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
static void setM3x3(const int k, double *__restrict__ M)
{
    memset(M, 0, 9*sizeof(double));
    // mrr
    if (k == 0)
    {
        M[0] = M_SQRT2;
    }
    // mtt
    else if (k == 1)
    {
        M[4] = M_SQRT2;
    }
    // mpp
    else if (k == 2)
    {
        M[8] = M_SQRT2;
    }
    // mrt and mtr
    else if (k == 3)
    {
        M[1] = 1.0;
        M[3] = 1.0;
    }
    // mrp and mpr
    else if (k == 4)
    {
        M[2] = 1.0;
        M[6] = 1.0;
    }
    // mtp and mpt
    else
    {
       M[5] = 1.0;
       M[7] = 1.0;
    }
    return;
}
