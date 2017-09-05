#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "iscl/memory/memory.h"

/*!
 * @brief Converts the NED based fundamental faults to Green's functions which
 *        scale the NED based moment tensor terms s.t. \f$ u = G m \f$.
 *        where \f$ m = \{ m_{xx}, m_{yy}, m_{zz}, m_{xy}, m_{xz}, m_{yz} \} \f$
 *        are the NED moment tensor terms and \f$ u \f$ is the estimate for
 *        the icomp'th component with source receiver azimuth, az, and given
 *        channel orientation encapsulated by cmpaz and cmpinc.
 *
 * @param[in] npgrns     number of points in greens functions
 * @param[in] ldg        leading dimension of G (>= 6)
 * @param[in] icomp      =1 for vertical Green's functions.
 *                       =2 for north (channel 2) Green's functions.
 *                       =3 for east (channel 3) Green's functions.
 * @param[in] azin       source to receiver azimuth (degrees)
 * @param[in] bazin      receiver to source azimuth (degrees).  if you 
 *                       do not know then set this to:
 *                         fmod(az + 180.0, 360.0)
 * @param[in] cmpaz      component azimuth (0 north, +90 east)
 * @param[in] cmpinc     component inclinantion (-90 up, 0 east/north, +90 down)
 * @param[in] ZDS        vertical greens fn for 90 degree dip slip [npgrns]
 * @param[in] ZSS        vertical greens fn for vertical strike slip [npgrns]
 * @param[in] ZDD        vertical greens fn for 45 degree dip slip [npgrns]
 * @param[in] ZEX        vertical greens fn for an explosion [npgrns]
 * @param[in] RDS        radial greens fn for 90 degree dip slip [npgrns]
 * @param[in] RSS        radial greens fn for vertical strike slip [npgrns]
 * @param[in] RDD        radial greens fn for 45 degree dip slip [npgrns]
 * @param[in] REX        radial greens fn for an explosion [npgrns]
 * @param[in] TDS        transverse greens fn for 90 degree dip slip [npgrns]
 * @param[in] TSS        transverse greens fn for vertical strike slip [npgrns]
 *
 * @param[out] G         row major Green's functions matrix [npgrns x ldg].
 *                       the columns are ordered: 
 *                        \f$ \{ G_{xx}, G_{yy}, G_{zz},
 *                               G_{xy}, G_{xz}, G_{yz} \}  \f$
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_utils_ff2mtGreensMatrix64f(const int npgrns, const int ldg,
                                     const int icomp,
                                     const double azin, const double bazin,
                                     const double cmpaz, const double cmpinc,
                                     const double *__restrict__ ZDS,
                                     const double *__restrict__ ZSS,
                                     const double *__restrict__ ZDD,
                                     const double *__restrict__ ZEX,
                                     const double *__restrict__ RDS,
                                     const double *__restrict__ RSS,
                                     const double *__restrict__ RDD,
                                     const double *__restrict__ REX,
                                     const double *__restrict__ TDS,
                                     const double *__restrict__ TSS,
                                     double *__restrict__ G)
{
    double *Gxx, *Gyy, *Gzz, *Gxy, *Gxz, *Gyz;
    int ierr;
    ierr = 0;
    if (npgrns < 1 || ldg < 6 || G == NULL)
    {
        if (npgrns < 1)
        {
            fprintf(stderr, "%s: No points in greens fns\n", __func__);
        }
        if (ldg < 6){fprintf(stderr, "%s: ldg must be >= 6\n", __func__);}
        if (G == NULL){fprintf(stderr, "%s: Error G is NULL\n", __func__);}
        return -1;
    }
    Gxx = memory_calloc64f(npgrns);
    Gyy = memory_calloc64f(npgrns);
    Gzz = memory_calloc64f(npgrns);
    Gxy = memory_calloc64f(npgrns);
    Gxz = memory_calloc64f(npgrns);
    Gyz = memory_calloc64f(npgrns);
    ierr = parmt_utils_ff2mtGreens64f(npgrns, icomp, azin, bazin,
                                      cmpaz, cmpinc,
                                      ZDS, ZSS, ZDD, ZEX,
                                      RDS, RSS, RDD, REX,
                                      TDS, TSS,
                                      Gxx, Gyy, Gzz,
                                      Gxy, Gxz, Gyz);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to compute greens functions columns\n",
                __func__);
    } 
    else
    {
        cblas_dcopy(npgrns, Gxx, 1, &G[0], ldg);
        cblas_dcopy(npgrns, Gyy, 1, &G[1], ldg);
        cblas_dcopy(npgrns, Gzz, 1, &G[2], ldg);
        cblas_dcopy(npgrns, Gxy, 1, &G[3], ldg);
        cblas_dcopy(npgrns, Gxz, 1, &G[4], ldg);
        cblas_dcopy(npgrns, Gyz, 1, &G[5], ldg);
    }
    memory_free64f(&Gxx);
    memory_free64f(&Gyy);
    memory_free64f(&Gzz);
    memory_free64f(&Gxy);
    memory_free64f(&Gxz);
    memory_free64f(&Gyz);
    return ierr;
}
//============================================================================//
/*!
 * @brief Converts the NED based fundamental faults to Green's functions which
 *        scale the NED based moment tensor terms s.t. \f$ u = G m \f$.
 *        where \f$ m = \{ m_{xx}, m_{yy}, m_{zz}, m_{xy}, m_{xz}, m_{yz} \} \f$
 *        are the NED moment tensor terms and \f$ u \f$ is the estimate for
 *        the icomp'th component with source receiver azimuth, az, and given
 *        channel orientation encapsulated by cmpaz and cmpinc. 
 *
 * @param[in] npgrns     number of points in greens functions
 * @param[in] icomp      =1 for vertical Green's functions.
 *                       =2 for north (channel 2) Green's functions.
 *                       =3 for east (channel 3) Green's functions.
 * @param[in] azin       source to receiver azimuth (degrees)
 * @param[in] bazin      receiver to source azimuth (degrees).  if you 
 *                       do not know then set this to:
 *                         fmod(az + 180.0, 360.0)
 * @param[in] cmpaz      component azimuth (0 north, +90 east)
 * @param[in] cmpinc     component inclinantion (-90 up, 0 east/north, +90 down)
 * @param[in] ZDS        vertical greens fn for 90 degree dip slip [npgrns]
 * @param[in] ZSS        vertical greens fn for vertical strike slip [npgrns]
 * @param[in] ZDD        vertical greens fn for 45 degree dip slip [npgrns]
 * @param[in] ZEX        vertical greens fn for an explosion [npgrns]
 * @param[in] RDS        radial greens fn for 90 degree dip slip [npgrns]
 * @param[in] RSS        radial greens fn for vertical strike slip [npgrns]
 * @param[in] RDD        radial greens fn for 45 degree dip slip [npgrns]
 * @param[in] REX        radial greens fn for an explosion [npgrns]
 * @param[in] TDS        transverse greens fn for 90 degree dip slip [npgrns]
 * @param[in] TSS        transverse greens fn for vertical strike slip [npgrns]
 *
 * @param[out] Gxx       greens function scaling mxx moment tensor term [npgrns]
 * @param[out] Gyy       greens function scaling myy moment tensor term [npgrns]
 * @param[out] Gzz       greens function scaling mzz moment tensor term [npgrns]
 * @param[out] Gxy       greens function scaling mxy moment tensor term [npgrns]
 * @param[out] Gxz       greens function scaling mxz moment tensor term [npgrns]
 * @param[out] Gyz       greens function scaling myz moment tensor term [npgrns]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_utils_ff2mtGreens64f(const int npgrns, const int icomp,
                               const double azin, const double bazin,
                               const double cmpaz, const double cmpinc,
                               const double *__restrict__ ZDS,
                               const double *__restrict__ ZSS,
                               const double *__restrict__ ZDD,
                               const double *__restrict__ ZEX,
                               const double *__restrict__ RDS,
                               const double *__restrict__ RSS,
                               const double *__restrict__ RDD,
                               const double *__restrict__ REX,
                               const double *__restrict__ TDS,
                               const double *__restrict__ TSS,
                               double *__restrict__ Gxx,
                               double *__restrict__ Gyy,
                               double *__restrict__ Gzz,
                               double *__restrict__ Gxy,
                               double *__restrict__ Gxz,
                               double *__restrict__ Gyz)
{
    double az, baz, c2a, ca, caz, cos_caz, cost,
           gxx_e, gxy_e, gxz_e, gyy_e, gyz_e, gzz_e,
           gxx_n, gxy_n, gxz_n, gyy_n, gyz_n, gzz_n,
           gxx_r, gxy_r, gxz_r, gyy_r, gyz_r, gzz_r,
           gxx_t, gxy_t, gxz_t, gyy_t, gyz_t, gzz_t,
           gxx_z, gxy_z, gxz_z, gyy_z, gyz_z, gzz_z,
           rex, rdd, rds, rss, s2a, sa, sin_caz, sint,
           tds, tss, xscal,
           zex, zdd, zds, zss;
    int i;
    const double pi180 = M_PI/180.0;
    const double half = 0.5; 
    const double sixth = 1.0/6.0;
    const double third = 1.0/3.0;
    //------------------------------------------------------------------------//
    if (icomp < 1 || icomp > 3)
    {
        fprintf(stderr, "%s: Invalid component: %d\n", __func__, icomp);
        return -1;
    }
    if (icomp == 1)
    {
        if (ZDS == NULL || ZSS == NULL || ZDD == NULL || ZEX == NULL)
        {
            if (ZDS == NULL){fprintf(stderr, "%s: zds is NULL\n", __func__);}
            if (ZSS == NULL){fprintf(stderr, "%s: zss is NULL\n", __func__);}
            if (ZDD == NULL){fprintf(stderr, "%s: zdd is NULL\n", __func__);}
            if (ZEX == NULL){fprintf(stderr, "%s: zex is NULL\n", __func__);}
            return -1;
        }
    }
    else
    {
        if (RDS == NULL || RSS == NULL || RDD == NULL || REX == NULL ||
            TDS == NULL || TSS == NULL)
        {
            if (RDS == NULL){fprintf(stderr, "%s: rds is NULL\n", __func__);}
            if (RSS == NULL){fprintf(stderr, "%s: rss is NULL\n", __func__);}
            if (RDD == NULL){fprintf(stderr, "%s: rdd is NULL\n", __func__);}
            if (REX == NULL){fprintf(stderr, "%s: rex is NULL\n", __func__);}
            if (TDS == NULL){fprintf(stderr, "%s: tds is NULL\n", __func__);}
            if (TSS == NULL){fprintf(stderr, "%s: tss is NULL\n", __func__);}
            return -1;
        }
    } 
    // Get the backazimuth
    az = azin;
    baz = bazin; //fmod(az + 180.0, 360.0);
    // Convert to radians
    az = az*pi180;
    baz = baz*pi180;
    caz = cmpaz*pi180;
    // Compute the geometric factors
    sa = sin(az);
    ca = cos(az);
    s2a = sin(2.*az);
    c2a = cos(2.*az);
    sint = sin(baz);
    cost = cos(baz);
    cos_caz = cos(caz);
    sin_caz = sin(caz);
    // the cmpaz will come from the metadata so fix that
    if (icomp == 3)
    {
       cos_caz = cos(caz - M_PI/2.0);
       sin_caz = sin(caz - M_PI/2.0);
    }
    // Set the columns (Minson et. al. 2008)
    xscal = 1.0;
    if (fabs(cmpinc - 90.0) < 1.e-4){xscal =-1.0;}
    if (icomp == 1)
    {
        #pragma omp simd
        for (i=0; i<npgrns; i++)
        {
            zss = ZSS[i];
            zdd = ZDD[i];
            zds = ZDS[i];
            zex = ZEX[i];
            // compute the greens function corresponding to the mt
            gxx_z = half*(zss*c2a) - sixth*zdd + third*zex;
            gyy_z =-half*(zss*c2a) - sixth*zdd + third*zex;
            gzz_z = third*(zdd + zex);
            gxy_z = zss*s2a;
            gxz_z = zds*ca;
            gyz_z = zds*sa;
            // copy and include polarity so it matches the observation
            Gxx[i] = xscal*gxx_z;
            Gyy[i] = xscal*gyy_z;
            Gzz[i] = xscal*gzz_z;
            Gxy[i] = xscal*gxy_z;
            Gxz[i] = xscal*gxz_z;
            Gyz[i] = xscal*gyz_z;
        } // Loop on points
    }
    else
    {
        // North or 2 component
        if (icomp == 2)
        {
            // Loop on rows (data points)
            //#pragma omp simd
            for (i=0; i<npgrns; i++)
            {
                rss = RSS[i];
                rdd = RDD[i];
                rds = RDS[i];
                rex = REX[i];
                tds = TDS[i];
                tss = TSS[i];

                gxx_r = half*(rss*c2a) - sixth*rdd + third*rex;
                gyy_r =-half*(rss*c2a) - sixth*rdd + third*rex;
                gzz_r = third*(rdd + rex);
                gxy_r = rss*s2a;
                gxz_r = rds*ca;
                gyz_r = rds*sa;
 
                gxx_t = half*tss*s2a;
                gyy_t =-half*tss*s2a;
                gzz_t = 0.0;
                gxy_t =-tss*c2a;
                gxz_t = tds*sa; 
                gyz_t =-tds*ca;

                // Rotate from (r, t) -> (n, e)
                //   n = cos(baz-180)*r - sin(baz-180)*t
                //     =-cos(baz)*r + sin(baz)*t
                //   e = sin(baz-180)*r + cos(baz-180)*t
                //     =-sin(baz)*r - cos(baz)*t
                gxx_n =-cost*gxx_r + sint*gxx_t;
                gxx_e =-sint*gxx_r - cost*gxx_t;
                gyy_n =-cost*gyy_r + sint*gyy_t;
                gyy_e =-sint*gyy_r - cost*gyy_t;
                gzz_n =-cost*gzz_r + sint*gzz_t;
                gzz_e =-sint*gzz_r - cost*gzz_t;
                gxy_n =-cost*gxy_r + sint*gxy_t;
                gxy_e =-sint*gxy_r - cost*gxy_t;
                gxz_n =-cost*gxz_r + sint*gxz_t;
                gxz_e =-sint*gxz_r - cost*gxz_t;
                gyz_n =-cost*gyz_r + sint*gyz_t;
                gyz_e =-sint*gyz_r - cost*gyz_t;
                // Rotate from (n, e) -> (1, 2) around az
                //   1 = cos(comp_az)*n + sin(comp_az)*e
                Gxx[i] = cos_caz*gxx_n + sin_caz*gxx_e;
                Gyy[i] = cos_caz*gyy_n + sin_caz*gyy_e;
                Gzz[i] = cos_caz*gzz_n + sin_caz*gzz_e;
                Gxy[i] = cos_caz*gxy_n + sin_caz*gxy_e;
                Gxz[i] = cos_caz*gxz_n + sin_caz*gxz_e;
                Gyz[i] = cos_caz*gyz_n + sin_caz*gyz_e;

            } // Loop on points
        }
        // East or 3 component
        else
        {
            // Loop on rows (data points) 
            #pragma omp simd
            for (i=0; i<npgrns; i++)
            {
                rss = RSS[i];
                rdd = RDD[i];
                rds = RDS[i];
                rex = REX[i];
                tds = TDS[i];
                tss = TSS[i];

                gxx_r = half*(rss*c2a) - sixth*rdd + third*rex;
                gyy_r =-half*(rss*c2a) - sixth*rdd + third*rex;
                gzz_r = third*(rdd + rex);
                gxy_r = rss*s2a;
                gxz_r = rds*ca;
                gyz_r = rds*sa;

                gxx_t = half*tss*s2a;
                gyy_t =-half*tss*s2a;
                gzz_t = 0.0;
                gxy_t =-tss*c2a;
                gxz_t = tds*sa;
                gyz_t =-tds*ca;

                // Rotate from (r, t) -> (n, e)
                //   n = cos(baz-180)*r - sin(baz-180)*t
                //     =-cos(baz)*r + sin(baz)*t
                //   e = sin(baz-180)*r + cos(baz-180)*t
                //     =-sin(baz)*r - cos(baz)*t
                gxx_n =-cost*gxx_r + sint*gxx_t;
                gxx_e =-sint*gxx_r - cost*gxx_t;
                gyy_n =-cost*gyy_r + sint*gyy_t;
                gyy_e =-sint*gyy_r - cost*gyy_t;
                gzz_n =-cost*gzz_r + sint*gzz_t;
                gzz_e =-sint*gzz_r - cost*gzz_t;
                gxy_n =-cost*gxy_r + sint*gxy_t;
                gxy_e =-sint*gxy_r - cost*gxy_t;
                gxz_n =-cost*gxz_r + sint*gxz_t;
                gxz_e =-sint*gxz_r - cost*gxz_t;
                gyz_n =-cost*gyz_r + sint*gyz_t;
                gyz_e =-sint*gyz_r - cost*gyz_t;
                // Rotate from (n, e) -> (1, 2) around az
                //   2 =-sin(comp_az)*n + cos(comp_az)*e
                Gxx[i] =-sin_caz*gxx_n + cos_caz*gxx_e;
                Gyy[i] =-sin_caz*gyy_n + cos_caz*gyy_e;
                Gzz[i] =-sin_caz*gzz_n + cos_caz*gzz_e;
                Gxy[i] =-sin_caz*gxy_n + cos_caz*gxy_e;
                Gxz[i] =-sin_caz*gxz_n + cos_caz*gxz_e;
                Gyz[i] =-sin_caz*gyz_n + cos_caz*gyz_e;
            } // Loop on points
        } // End check on 1 or 2 component
    } // End check on component
    return 0;
}
