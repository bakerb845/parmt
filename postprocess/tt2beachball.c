#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include "parmt_postProcess.h"
#include "beachball.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "compearth.h"

static void argsort3(const double *__restrict__ x, int *__restrict__ iperm);
static void sortUL(const double *__restrict__ lam,
                   const double *__restrict__ U,
                   double *__restrict__ lamSort,
                   double *__restrict__ Usort);
static void eigenvector2PrincipalAxis(const double eig,
                                      const double *__restrict__ ev, 
                                      double *__restrict__ paxis);
static void vecstrd(const double str, const double dip, double r[3]);
static double dot3(const double *__restrict__ x,
                   const double *__restrict__ y);
static void gemv3(const double *__restrict__ A, const double *__restrict__ x,
                  double *__restrict__ y);
static void crossp(const double *__restrict__ u,
                   const double *__restrict__ v,
                   double *__restrict__ n);
static void strdvec(double r[3], double *str, double *dip);

int postprocess_tnp2beachballPolarity(const int nxp,
                                      const double xc, const double yc,
                                      const double rad,
                                      const double *__restrict__ pAxis,
                                      const double *__restrict__ nAxis,
                                      const double *__restrict__ tAxis,
                                      double *__restrict__ xw1,
                                      double *__restrict__ yw1,
                                      int8_t *__restrict__ pn1)
{
    double ax1[3] __attribute__ ((aligned(64)));
    double ax2[3] __attribute__ ((aligned(64)));
    double ax3[3] __attribute__ ((aligned(64)));
    double rm[9] __attribute__ ((aligned(64)));
    double rmt[9] __attribute__ ((aligned(64)));
    double tm[9] __attribute__ ((aligned(64)));
    double pt[3] __attribute__ ((aligned(64)));
    double ptnew[3] __attribute__ ((aligned(64)));
    double ampMax, az, baz, baz1, bdp, bdp1, bev, bevI, dp,
           dx, dxc, dy, dyc, evI,
           fclvd, fIso, paz, paz1, pdp, pdp1, pev, pevI, plusMin,
           rad2, rpt, rpt2,
           taz, taz1, tdp, tdp1, tev, tevI, theta, x, x0, y, y0;
    int ix, iy, k, nx, ny;
    const double degradi = 180.0/M_PI;
    const double third = 1.0/3.0;
    const double sqrt2 = 1.4142135623730951;
    const double sqrt2i = 1.0/sqrt2;
    taz1 = tAxis[0];
    tdp1 = tAxis[1];
    tevI = tAxis[2];
    baz1 = nAxis[0];
    bdp1 = nAxis[1];
    bevI = nAxis[2];
    paz1 = pAxis[0];
    pdp1 = pAxis[1];
    pevI = pAxis[2];
    // Avoid overwriting memory
    taz = taz1;
    tdp = tdp1;
    baz = baz1;
    bdp = bdp1;
    paz = paz1;
    pdp = pdp1;
    // Step 1: set up axes:
    //  - D axis ('dominant' 1-axis - with largest absolute eigenvalue)
    //  - M axis ('minor'    3-axis)
    //  - B axis (           2-axis)
    evI = (tevI + bevI + pevI)*third;
    tev = tevI - evI;
    bev = bevI - evI;
    pev = pevI - evI;
    vecstrd(paz, pdp, ax1);
    vecstrd(taz, tdp, ax3);
    fclvd = fabs(bev/pev);
    fIso = evI/pev;
    if (fabs(tev) > fabs(pev))
    {
        vecstrd(taz, tdp, ax1);
        vecstrd(paz, pdp, ax3);
        fclvd = fabs(bev/tev);
        fIso = evI/tev;
    }
    //if (fabs(evI) > .03) 
    //{
    //    write(6,901) evI
    //}
    crossp(ax1, ax3, ax2);
    strdvec(ax2, &baz, &bdp);
    cliffsNodes_mtensor(taz, tdp, paz, pdp, tevI, bevI, pevI,
                        tm, rm, rmt, &ampMax);
    nx = nxp;
    ny = nxp;
    memset(xw1, 0, (size_t) (nx*ny)*sizeof(double));
    memset(yw1, 0, (size_t) (nx*ny)*sizeof(double));
    memset(pn1, 0, (size_t) (nx*ny)*sizeof(int8_t));
    dx = 2.0*rad/(double) (nx - 1);
    dy = 2.0*rad/(double) (ny - 1);
    x0 = xc - rad;
    y0 = yc - rad;
    rad2 = rad*rad; // saves a sqrt computation when comparing distance
    // loop on pixel grid
    for (iy=0; iy<ny; iy++)
    {
        for (ix=0; ix<nx; ix++)
        {
            x = x0 + (double) ix*dx;
            y = y0 + (double) iy*dy;
            k = iy*nx + ix;
            pn1[k] = 0;
            xw1[k] = x;
            yw1[k] = y;
            // require (x,y) be in the focal sphere
            dxc = x - xc;
            dyc = y - yc;
            rpt2 = dxc*dxc + dyc*dyc;
            if (rpt2 < rad2)
            {
                rpt = sqrt(rpt2);
                // compute the azimuth
                theta = atan2(y - yc, x - xc);  // measure from x axis
                az = (M_PI_2 - theta)*degradi;  // az is measured from north
                if (az < 0.0){az = az + 360.0;} // convention is [0,360]
                // compute the dip
                dp = (M_PI_2 - 2.0*asin(rpt*sqrt2i))*degradi;
                vecstrd(az, dp, pt);
                gemv3(tm, pt, ptnew);
                plusMin = dot3(pt, ptnew);
                pn1[k] =-1; 
                if (plusMin > 100.0*DBL_EPSILON*ampMax){pn1[k] = 1;}
            } // end check on location
        } // loop on x
    } // loop on y
    return 0;
}
/*!
 * @brief Converts a Tape and Tape 2015 tensor to tension, null, and
 *        plunge eigenvectors in USE coordinates for plotting with
 *        cliffsnodes
 * @param[in] beta    colatitude (radians)
 * @param[in] gamma   longitude (radians)
 * @param[in] kappa   strike angle (radians)
 * @param[in] sigma   slip angle (radians)
 * @param[in] theta   dip angle (radians)
 *
 * @param[out] pAxis  (azimuth, plunge, eigenvalue) of pressure axis [3]
 * @param[out] nAxis  (azimuth, plunge, eigenvalue) of null axis [3]
 * @param[out] tAxis  (azimuth, plunge, eigenvalue) of tension axis [3]
 *
 * @result 0 indicate success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int postprocess_tt2tnp(const double beta, const double gamma,
                       const double kappa, const double sigma,
                       const double theta,
                       double *__restrict__ pAxis,
                       double *__restrict__ nAxis,
                       double *__restrict__ tAxis)
{
    const double R[9] = {
             0.7071067811865476, 0.0,                -0.7071067811865476,
            -0.4082482904638631, 0.8164965809277261, -0.4082482904638631, 
             0.5773502691896258, 0.5773502691896258,  0.5773502691896258};
    const double Yrot[9] = {
             0.7071067811865475, 0.0, 0.7071067811865475,
             0.0, 1.0, 0.0,
            -0.7071067811865475, 0.0, 0.7071067811865475};
    // Convert NWU to USE
    const double Rot[9] = {0, -1, 0, 0, 0, -1, 1, 0, 0};
    //const double RotInv[9] = {0, 0, 1, -1, 0, 0, 0, -1, 0};
    double Zsigma[9], Xtheta[9], Zkappa[9], ZX[9], V[9],
           U[9], Uuse[9], Usort[9], Uw[9],
           lambda[3], v3[3], lamSort[3],
           cosb, cosg, sinb, sing, nKappaDeg, sigmaDeg, thetaDeg;
    const double toDeg = 180.0/M_PI;
    // Compute the eigenvalues from Eqn 7
    sinb = sin(beta);
    sing = sin(gamma);
    cosb = cos(beta);
    cosg = cos(gamma);
    v3[0] = sinb*cosg;
    v3[1] = sinb*sing;
    v3[2] = cosb;
    // Compute R*lam 
    lambda[0] = R[0]*v3[0] + R[3]*v3[1] + R[6]*v3[2];
    lambda[1] = R[1]*v3[0] + R[4]*v3[1] + R[7]*v3[2];
    lambda[2] = R[2]*v3[0] + R[5]*v3[1] + R[8]*v3[2];
    // Eqn 9-10
    sigmaDeg  = sigma*toDeg;
    thetaDeg  = theta*toDeg;
    nKappaDeg =-kappa*toDeg;
    compearth_eulerUtil_rotmat(1, &sigmaDeg, 3, Zsigma);
    compearth_eulerUtil_rotmat(1, &thetaDeg, 1, Xtheta);
    compearth_eulerUtil_rotmat(1, &nKappaDeg, 3, Zkappa);
    //compearth_eulerUtil_rotmatGen(1,-45.0, 2, Yrot);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, Zkappa, 3, Xtheta, 3, 0.0, ZX, 3); 
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, ZX, 3, Zsigma, 3, 0.0, V, 3); 
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, V, 3, Yrot, 3, 0.0, U, 3);
    // M could be computed by computing M = U*Lambda*inv(U) hence
    // U corresponds to the eigenvectors and lambda the scaling.
    // However, the current basis is NWU which must be converted to USE.
    // If performing a proper rotation of M to M_{use} we would compute:
    //   R M inv(R) = R U Lambda inv(U) inv(R) 
    // thus what we must be do is transform the eigenvectors into the
    // apropriate frame.  Intuitively it is the rotation of each 
    // eigenvector (column) of the U matrix which is indeed the matrix
    // matrix multiplication of R*U.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, Rot, 3, U, 3, 0.0, Uuse, 3);
    // If I want to finish this I could compute:
    //   R U Lambda inv(U) inv(R) = R U Lambda U^T R^T
/*
    double lamMat[9], Temp[9], M[9], UtRt[9];
    memset(lamMat, 0, 9*sizeof(double));
    lamMat[0] = lambda[0]; lamMat[4] = lambda[1]; lamMat[8] = lambda[2];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, Uuse, 3, lamMat, 3, 0.0, Temp, 3);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans,
                3, 3, 3, 1.0, U, 3, Rot, 3, 0.0, UtRt, 3); 
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, Temp, 3, UtRt, 3, 0.0, M, 3);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            printf("%f ", M[3*j+i]);
        }
        printf("\n");
     }
*/
    sortUL(lambda, Uuse, lamSort, Usort);
//printf("%f %f %f\n", lamSort[0], lamSort[1], lamSort[2]);
    // Compute pressure, null, and tension axes
    eigenvector2PrincipalAxis(lamSort[0], &Usort[0], pAxis);
    eigenvector2PrincipalAxis(lamSort[1], &Usort[3], nAxis);
    eigenvector2PrincipalAxis(lamSort[2], &Usort[6], tAxis);
    // Explosion - simplification (all positive eigenvalues)
    if (lamSort[0] > 0.0)
    {
        pAxis[2] = 1.0;
        nAxis[2] = 1.0;
        tAxis[2] = 1.0;
    }
    // Implosion - simplification (all negative eigenvalues)
    if (lamSort[2] < 0.0)
    {
        pAxis[0] =-1.0;
        nAxis[1] =-1.0;
        tAxis[2] =-1.0;
    }
    //printf("%f %f %f\n", pAxis[0], pAxis[1], pAxis[2]);
    //printf("%f %f %f\n", nAxis[0], nAxis[1], nAxis[2]);
    //printf("%f %f %f\n", tAxis[0], tAxis[1], tAxis[2]);
    return 0;
}

/*!
 * @brief Convert pressure, null, and tension principal axes to azimuth,
 *        plunge, and length.
 */
static void eigenvector2PrincipalAxis(const double eig,
                                      const double *__restrict__ ev,
                                      double *__restrict__ paxis)
{
    const char *fcnm = "cmopad_Eigenvector2PrincipalAxis\0";
    double v[3], az, plunge;
    double twopi = 2.0*M_PI;
    double pi180i = 180.0/M_PI;
    int ierr;
    //------------------------------------------------------------------------//    //
    // Compute azimuth and plunge angles where ev is USE coordinates
    plunge = asin(-ev[0]);        //arcsin(-z/r)
    az     = atan2(ev[2],-ev[1]); //atan(x/y)
    if (plunge <= 0.0)
    {    
        plunge =-plunge;
        az = az + M_PI;
    }
    if (az < 0.0){az = az + twopi;}   //Shift to [0,360]
    if (az > twopi){az = az - twopi;} //Shift to [0,360]
    // Convert to degrees
    plunge = plunge*pi180i;
    az = az*pi180i;
    // Copy back result
    paxis[0] = az;     //First value is azimuth (degrees)
    paxis[1] = plunge; //Second value is plunge (degrees)
    paxis[2] = eig;    //Final value is eigenvalue (distance)
    return;
}
//===========================================================================//
static void sortUL(const double *__restrict__ lam,
                   const double *__restrict__ U,
                   double *__restrict__ lamSort,
                   double *__restrict__ Usort)
{
    int perm[3], j;
    argsort3(lam, perm);
    lamSort[0] = lam[perm[0]];
    lamSort[1] = lam[perm[1]];
    lamSort[2] = lam[perm[2]];
    for (j=0; j<3; j++)
    {
        Usort[3*j+0] = U[3*perm[j]+0];
        Usort[3*j+1] = U[3*perm[j]+1];
        Usort[3*j+2] = U[3*perm[j]+2];
    }
    return;
}

//============================================================================//
/*!
 * @brief Permutation to sort length 3 array into ascending order
 */
static void argsort3(const double *__restrict__ x, int *__restrict__ iperm)
{
    const char *fcnm = "__cmopad_argsort3\0";
    int i, temp;
    const int a = 0; 
    const int b = 1; 
    const int c = 2; 
    // Copy
    iperm[a] = a; 
    iperm[b] = b; 
    iperm[c] = c; 
    if (x[iperm[a]] > x[iperm[c]])
    {    
        temp = iperm[c];
        iperm[c] = iperm[a];
        iperm[a] = temp;
    }    
    if (x[iperm[a]] > x[iperm[b]])
    {    
        temp = iperm[b];
        iperm[b] = iperm[a];
        iperm[a] = temp;
    }    
    //Now the smallest element is the first one. Just check the 2-nd and 3-rd
    if (x[iperm[b]] > x[iperm[c]])
    {    
        temp = iperm[c];
        iperm[c] = iperm[b];
        iperm[b] = temp;
    }
    return;
}
//============================================================================//
/*!
 * @brief Find components of downward pointing unit vector pole from the 
 *        strike and dip angles.  This is based on Cliff Frolich's vecstrd.
 *
 * @param[in] str    strike angle of pole (degrees)
 * @param[in] dip    dip angle of pole (degrees)
 *
 * @param[out] r     corresponding unit vector pole (3)
 *
 * @author Ben Baker
 *
 * @copyright BSD
 * 
 */
static void vecstrd(const double str, const double dip, double r[3])
{
    double cosdip, cosstr, diprad, sindip, sinstr, strrad;
    const double degrad = M_PI/180.0;
    // convert to radians
    strrad = str*degrad;
    diprad = dip*degrad;
    // compute cosines and sines
    sinstr = sin(strrad);
    cosstr = cos(strrad);
    cosdip = cos(diprad);
    sindip = sin(diprad);
    // compute (x,y,z) from strike and dip for unit radius 
    r[0] = sinstr*cosdip;
    r[1] = cosstr*cosdip;
    r[2] =-sindip;
    return;
}
/*!
 * @brief Utility function for computing dot product of length 3 vector
 */
static double dot3(const double *__restrict__ x,
                   const double *__restrict__ y)
{
    double dot;
    dot = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
    return dot;
}
/*!
 * @brief Utility function for computing y = Ax for 3 x 3 matrix.
 */
static void gemv3(const double *__restrict__ A, const double *__restrict__ x,
                  double *__restrict__ y)
{
    y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
    y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
    y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
    return;
}
/*!
 * @brief Cross product of 2 length 3 vectors u and v.  Returns normal vector n.
 * 
 * @param[in] u    first vector in cross product u x v [3]
 * @param[in] v    second vector in cross product u x v [3]
 * 
 * @param[out] n   normal vector n = u x v [3]
 *
 * @author Ben Baker
 *
 * @copyright BSD
 *
 */
static void crossp(const double *__restrict__ u,
                   const double *__restrict__ v,
                   double *__restrict__ n)
{
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
    return;
}
/*!
 * @brief Finds the strike and dip of downward pointing unit vector.  This
 *        is based on Cliff Frolich's strdvec.
 *
 * @param[in,out] r      on input contains the downward pointing unit
 *                       vector.
 *                       if r[2] is > 0 then the unit vector is pointing
 *                       updward.  in this case then on output r3 will
 *                       be negated.
 *
 * @param[out] str       strike angle (degrees) \f$ s.t. \phi \in [0,360] \f$.
 * @param[out] dip       dip angle (degrees) 
 *
 * @author Ben Baker
 *
 * @copyright BSD
 *
 */
static void strdvec(double r[3], double *str, double *dip)
{
    double rs;
    const double degradi = 180.0/M_PI;
    if (r[2] > 0.0)
    {    
        r[0] =-r[0];
        r[1] =-r[1];
        r[2] =-r[2];
    }        
    *str = atan2(r[0], r[1])*degradi;
    if (*str < 0.){*str = *str + 360.;}
    rs = sqrt(r[0]*r[0] + r[1]*r[1]); 
    *dip=atan2(-r[2], rs)*degradi;
    return;
}
