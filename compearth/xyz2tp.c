#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/*!
 * @brief This function takes a matrix of 3-vectors in xyz and
 *        returns a vector of theta and phi values.
 *
 * @param[in] n      number of ordinates
 * @param[in] x      x coordinates [n]
 * @param[in] y      y coordinates [n]
 * @param[in] z      z coordiantes [n]
 *
 * @param[out] th    polar angles in radians [n]
 * @param[out] ph    azimuthal angles in radians [n] 
 * @param[out] rho   radii in units of (x, y, z) [n]
 *
 * @author Carl Tape and converted to C by Ben Baker
 *
 */
void compearth_xyz2tp(const int n,
                      const double *__restrict__ x,
                      const double *__restrict__ y,
                      const double *__restrict__ z,
                      double *__restrict__ ph,
                      double *__restrict__ th,
                      double *__restrict__ rho)
{
    const double pi2 = M_PI/2.0;
    int i;
    // compute the colatitude, longitude, and radius from the
    // x,y,z coordinates
    compearth_matlab_cart2sph(n, x, y, z,
                              th, ph, rho);
    // convert to latitude
    for (i=0; i<n; i++)
    {
        th[i] = M_PI_2 - th[i];
    }
    return;
}
