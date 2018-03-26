#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

#ifndef M_PI 
#define M_PI   3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

static int diffConstantArray(const int nx,
                             const double xmin, const double xmax,
                             const double *__restrict__ x,
                             double *__restrict__ dx);

/*!
 * @brief Computes the non-uniform grid spacing in beta.  This supposes that
 *        the discretization grid in u is uniform and that each u represents
 *        a cell center.  For the purposes of the quadrature it is assumed that
 *        the objective function is constant over the cell.
 *
 * @param[in] nb       number of betas 
 * @param[in] betas    lune colatitudes (radians) at cell `centers' [nb]
 *
 * @param[out] db      beta cell-spacings (radians) [nb]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int postprocess_computeBetaCellSpacing(const int nb,
                                       const double *__restrict__ betas,
                                       double *__restrict__ db)
{
    double *diff, *u, b1, b2, du, du2, u1, u2;
    int i, ierr; 
    const double threePi4 = 0.75*M_PI;
    const double uMin = 0.0;
    const double uMax = threePi4;
    const double betaMin = 0.0;
    const double betaMax = M_PI;
    // Base case - set to entire interval
    if (nb == 1)
    {
        db[0] = (betaMax - betaMin);
        return 0;
    }
    // Convert betas to u
    u = memory_calloc64f(nb);
    compearth_beta2u(nb, betas, u);
    // Compute the average grid spacing in u
    diff = array_diff64f(nb, u, &ierr); 
    du = array_sum64f(nb-1, diff, &ierr)/(double) (nb - 1);
    du2 = 0.5*du;
    // Turn the cellular betas into cell spacings
    for (i=0; i<nb; i++)
    {
        u1 = u[i] - du2;
        u2 = u[i] + du2; 
        if (i == 0){u1 = uMin;}
        if (i == nb - 1){u2 = uMax;}
        u1 = fmax(uMin, u1);
        u2 = fmin(uMax, u2);
        compearth_u2beta(1, &u1, &b1); //compearth_u2beta(1, 50, 2, &u1, 1.e-8, &b1);
        compearth_u2beta(1, &u2, &b2); //compearth_u2beta(1, 50, 2, &u2, 1.e-8, &b2);
        db[i] = b2 - b1;
        //printf("%f %f %f\n", db[i], u1, u2);
    }
    memory_free64f(&diff);
    memory_free64f(&u);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the non-uniform grid spacing in gamma.  This supposes that
 *        the discretization grid in v is uniform and that each v represents
 *        a cell center.  For the purposes of the quadrature it is assumed that
 *        the objective function is constant over the cell.
 *
 * @param[in] ng       number of gammas
 * @param[in] gammas   lune longitudes (radians) at cell `centers' [ng]
 *
 * @param[out] dg      gamma cell-spacings (radians) [ng]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int postprocess_computeGammaCellSpacing(const int ng,
                                        const double *__restrict__ gammas,
                                        double *__restrict__ dg)
{
    double *diff, *v, g1, g2, dv, dv2, v1, v2;
    int i, ierr; 
    const double vMin =-1.0/3.0;
    const double vMax = 1.0/3.0;
    const double gammaMin =-M_PI/6.0;
    const double gammaMax = M_PI/6.0;
    // Base case - set to entire interval
    if (ng == 1)
    {
        dg[0] = (gammaMax - gammaMin);
        return 0;
    }
    // Convert betas to u
    v = memory_calloc64f(ng);
    compearth_gamma2v(ng, gammas, v);
    // Compute the average grid spacing in v
    diff = array_diff64f(ng, v, &ierr);
    dv = array_sum64f(ng-1, diff, &ierr)/(double) (ng - 1);
    dv2 = 0.5*dv;
    // Turn the cellular gammas into cell spacings
    for (i=0; i<ng; i++)
    {
        v1 = v[i] - dv2;
        v2 = v[i] + dv2; 
        if (i == 0){v1 = vMin;}
        if (i == ng - 1){v2 = vMax;}
        v1 = fmax(vMin, v1);
        v2 = fmin(vMax, v2);
        compearth_v2gamma(1, &v1, &g1);
        compearth_v2gamma(1, &v2, &g2);
        dg[i] = g2 - g1; 
        //printf("%f %f %f\n", dg[i], v1, v2);
    }
    memory_free64f(&diff);
    memory_free64f(&v);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the non-uniform grid spacing in theta.  This supposes that
 *        the discretization grid in h is uniform and that each h represents
 *        a cell center.  For the purposes of the quadrature it is assumed that
 *        the objective function is constant over the cell.
 *
 * @param[in] nt       number of thetas 
 * @param[in] thetas   dip angles (radians) at cell `centers' [nt]
 *
 * @param[out] dt      theta cell-spacings (radians) [nt]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int postprocess_computeThetaCellSpacing(const int nt,
                                        const double *__restrict__ thetas,
                                        double *__restrict__ dt)
{
    double *diff, *h, h1, h2, dh, dh2, t1, t2;
    int i, ierr;
    const double hMin = 1.0;
    const double hMax = 0.0;
    const double thetaMin = 0.0;
    const double thetaMax = M_PI_2;
    // Base case - set to entire interval
    if (nt == 1)
    {
        dt[0] = (thetaMax - thetaMin);
        return 0;
    }
    // Convert betas to u
    h = memory_calloc64f(nt);
    compearth_theta2h(nt, thetas, h);
    // Compute the average grid spacing in v
    diff = array_diff64f(nt, h, &ierr);
    dh = array_sum64f(nt-1, diff, &ierr)/(double) (nt - 1);
    dh2 = 0.5*dh;
    // Turn the cellular thetas into cell spacings
    for (i=0; i<nt; i++)
    {
        h1 = h[i] - dh2;
        h2 = h[i] + dh2;
        if (i == 0){h1 = hMin;}
        if (i == nt - 1){h2 = hMax;}
        h1 = fmin(hMin, h1);
        h2 = fmax(hMax, h2);
        compearth_h2theta(1, &h1, &t1);
        compearth_h2theta(1, &h2, &t2);
        dt[i] = t2 - t1;
        //printf("%f %f %f\n", dt[i], h1, h2);
    }
    memory_free64f(&diff);
    memory_free64f(&h);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the cell spacing in strike angles.  For the purposes of the
 *        quadrature the objective function is assumed constant over the cell.
 *
 * @param[in] nk      number of strikes 
 * @param[in] kappas  kappa angles (radians) [ns]
 *
 * @param[out] dk     strike cell spacings (sigmas) [ns]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int postprocess_computeKappaCellSpacing(const int nk, 
                                        const double *__restrict__ kappas,
                                        double *__restrict__ dk) 
{
    int ierr;
    const double kappaMin = 0.0;
    const double kappaMax = 2.0*M_PI;
    ierr = diffConstantArray(nk, kappaMin, kappaMax, kappas, dk);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the cell spacing in slip angles.  For the purposes of the
 *        quadrature the objective function is assumed constant over the cell.
 *
 * @param[in] ns      number of slips 
 * @param[in] sigmas  slip angles (radians) [ns]
 *
 * @param[out] ds     slip cell spacings (sigmas) [ns]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int postprocess_computeSigmaCellSpacing(const int ns,
                                        const double *__restrict__ sigmas,
                                        double *__restrict__ ds)
{
    int ierr;
    const double sigmaMin =-M_PI_2;
    const double sigmaMax = M_PI_2;
    ierr = diffConstantArray(ns, sigmaMin, sigmaMax, sigmas, ds);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the cell spacing in scalar moment.  For the purposes of
 *        the quadrature the objective function is assumed constant over the
 *        cell.
 *
 * @param[in] nm    number of scalar moments
 * @param[in] M0s   scalar moments (Newton-meters) [nm]
 *
 * @param[out] dm   scalar moment cell spacings (Newton-meters) [nm] 
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int postprocess_computeM0CellSpacing(const int nm,
                                     const double *__restrict__ M0s,
                                     double *__restrict__ dm)
{
    int ierr;
    const double m0Min = 0.0;
    const double m0Max = 1.0;
    ierr = diffConstantArray(nm, m0Min, m0Max, M0s, dm);
    return ierr;
}
//============================================================================//
/*!
 * @brief Private utility for assigning constant cell spacings to the
 *        uniformly spaced arrays
 */
static int diffConstantArray(const int nx,
                             const double xmin, const double xmax,
                             const double *__restrict__ x,
                             double *__restrict__ dx)
{
    double *diff, d, d2;
    int ierr;
    if (nx == 1)
    {
        dx[0] = xmax - xmin;
        return 0;
    }
    diff = array_diff64f(nx, x, &ierr);
    d = array_sum64f(nx-1, diff, &ierr)/(double) (nx - 1);
    ierr = array_set64f_work(nx, d, dx);
    memory_free64f(&diff);
    return ierr;
}
