#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct luneImage_struct
{
    double *xWesternLine;  /*!< x locations defining the eastern line
                                [nWesternLine] */
    double *yWesternLine;  /*!< y locations defining the western line
                                [nWesternLine] */
    double *xEasternLine;  /*!< x locations defining the eastern line
                                [nEasternLine] */
    double *yEasternLine;  /*!< y locations defining the eastern line
                                [nEasternLine] */
    double *uPts;          /*!< u locations in gridsearch */
    double *vPts;          /*!< v locations in gridsearch */ 
    double *xPts;          /*!< x points in lune image */
    double *yPts;          /*!< y points in lune image */ 
    int nWesternLine;      /*!< Number of points in western line */
    int nEasternLine;      /*!< Number of points in eastern line */
    int npts;              /*!< Number of points (pixels) in lune image */
    char pad[4];
};

void projections_hammer_latLonToXY(const double lambda, const double phi,
                                   double *x, double *y);
void projections_hammer_xyToLatLon(const double x, const double y,
                                   double *lambda, double *phi);
static int findNextNorthPointForFixedLongitude(
    const double lam0, const double phi0,
    const double dy, double *phi1);
static int findNextEastPointForFixedLatitude(
    const double lam0, const double phi0,
    const double dx, double *lam1);
#pragma omp declare simd
static double transform(const double a, const double b,
                        const double c, const double d, const double x);

int main()
{
    double dPhi, dLam, x, y, lambda, lambda1, phi, lamMin, lamMax, phiMin, phiMax;
    double *x0s, *y0s, dx, dy, dyref, phi1, u, umin, umax, v, vmin, vmax, xmin, xmax, xwidth, y0;
    double xpixel, ypixel;
    int ierr, ix, iy, nlat, nlon;

    int nyPixel = 122; // recommend this is even
    int nxPixel = nyPixel/3.0; // aspect ratio is 1:3
    const double xlow  =-5.22104768880206e-01; //-sqrt(2.0);
    const double xhigh = 5.22104768880206e-01; //sqrt(2.0);
    const double ylow  =-sqrt(2.0);
    const double yhigh = sqrt(2.0);

    phiMin =-M_PI_2 + 1.e-14; // numerical stability problem
    phiMax = M_PI_2 - 1.e-14; // numerical stability problem
    lamMin =-M_PI/6.0;
    lamMax = M_PI/6.0;
/*
projections_hammer_latLonToXY(0, phiMin, &x, &y);
printf("%16.14e %16.14e\n", x, y);
projections_hammer_latLonToXY(0, phiMax, &x, &y);
printf("%16.14e %16.14e\n", x, y);
getchar();
*/
    nlon = nxPixel;
    nlat = nyPixel;
    dPhi = (phiMax - phiMin)/(double) (nlat - 1);
    dLam = (lamMax - lamMin)/(double) (nlon - 1);
    x0s = (double *) aligned_alloc(64, (size_t) (nlat+1)*sizeof(double));
    y0s = (double *) aligned_alloc(64, (size_t) (nlat+1)*sizeof(double));
/*
    for (iy=0; iy<nlat; iy++)
    {
        for (ix=0; ix<nlon; ix++)
        {
            lambda = lamMin + (double) ix*dLam;
            phi = phiMin + (double) iy*dPhi;
            projections_hammer_latLonToXY(lambda, phi, &x, &y);
            //projections_hammer_xyToLatLon(x, y, &l, &p);
//printf("%f %f %f\n", x, y, sqrt(x*x + y*y));
        }
//printf("\n");
if (iy < nlat - 1){findNextNorthPointForFixedLongitude(lamMin, phi, 0.1, &phi1);}
    } 
*/
    // evaluate the latitudes on the left side
    dy = (phiMax - phiMin)/(double) (nlat - 1);
    //dy = M_PI/(double) (nlat - 1); //(yylow
    //dy = 1.0/(double) (nlat - 1);
    //dy = 1.0/(double) (nyPixel - 1);
//printf("%f\n", dy);
    dy = transform(0, (double) (nyPixel - 1), ylow, yhigh, 1.0) //phiMin, phiMax, 1.0)
       - transform(0, (double) (nyPixel - 1), ylow, yhigh, 0.0);//phiMin, phiMax, 0.0);
//printf("%f\n", dy);
//printf("%f\n", dy);
//return 0;
    phi = phiMin;
    projections_hammer_latLonToXY(lamMin, phiMin, &x0s[0], &y0s[0]);
//printf("%f\n", phi);
    for (iy=0; iy<nlat; iy++)
    {
        ierr = findNextNorthPointForFixedLongitude(lamMin, phi, dy, &phi1);
        // out of space - draw the top point
        if (ierr == 1)
        {
            projections_hammer_latLonToXY(lamMin, phiMax, &x, &y);
            iy = iy + 1;
            x0s[iy+1] = x;
            y0s[iy+1] = y; 
            nlat = iy + 2;
            break;
        }
        projections_hammer_latLonToXY(lamMin, phi1, &x, &y);
        if (y > yhigh)
        {
            projections_hammer_latLonToXY(lamMin, phiMax, &x, &y);
            x0s[iy+1] = x;
            y0s[iy+1] = y;
            nlat = iy + 2;
            break;
        }
        x0s[iy+1] = x;
        y0s[iy+1] = y;
        dyref = y0s[iy+1] - y0s[iy];
//printf("%f\n", dyref);
        phi = phi1;
    }
//getchar();
    // draw the other half 
    for (iy=0; iy<nlat; iy++)
    {
//printf("%f %f\n", x0s[iy], y0s[iy]);
        //printf("%f %f\n", -x0s[iy], y0s[iy]);
        //projections_hammer_xyToLatLon(-x0s[iy], y0s[iy],
        //                              &lambda, &phi);
    }
    //xwidth = transform(-M_PI, M_PI, 0, (double) (nxPixel - 1), M_PI);
    // Evaluate the longitudes for each latitude
    dx = 1.0; // pixel width is unity
    for (iy=0; iy<nlat; iy++)
    {
        y = y0s[iy];
        // transform x0 into the pixel domain 
        //xmin = transform(-M_PI/6.0, M_PI/6.0, 0, (double) (nxPixel - 1), x0s[iy]);
        //xmax = transform(-M_PI/6.0, M_PI/6.0, 0, (double) (nxPixel - 1),-x0s[iy]);
        xmin = transform(xlow, xhigh, 0, (double) (nxPixel - 1), x0s[iy]);
        xmax = transform(xlow, xhigh, 0, (double) (nxPixel - 1),-x0s[iy]);
        xmin = (double) ((int) xmin) + dx;
        xmax = (double) ((int) xmax) - dx;
        nlon = (int) ((xmax - xmin)/dx + 0.5) + 1;
//printf("%f\n", x);
//getchar();
        for (ix=0; ix<nlon; ix++)
        {
            //if (x > xmax){break;}
            //xpixel = transform(-M_PI, M_PI, 0, (double) (nxPixel - 1), x);
            xpixel = xmin + (double) ix*dx;
            //ypixel = transform(-M_PI_2, M_PI_2, 0, (double) (nyPixel - 1), y);
            ypixel = transform(ylow, yhigh, 0, (double) (nyPixel - 1), y);
            xpixel = (double) ((int) xpixel);
            //ypixel = (double) ((int) ypixel);
printf("%f %f %f\n", xpixel, ypixel, y0s[iy]);
            //x = transform(0, (double) (nxPixel - 1),-M_PI/6.0, M_PI/6.0, xpixel);
            //y = transform(0, (double) (nyPixel - 1),-M_PI_2, M_PI_2, ypixel); 
            x = transform(0, (double) (nxPixel - 1), xlow, xhigh, xpixel);
            y = transform(0, (double) (nyPixel - 1), ylow, yhigh, ypixel);
            //if (x >-x0s[iy]){break;}
            projections_hammer_xyToLatLon(x, y,
                                          &lambda, &phi);
//if (lambda > 1.0){printf("%f %f %f %f\n", x, y, lambda, phi); getchar();}
            // (u, v) space
            u = 0.75*sin(M_PI_2 - phi)
              - 0.5*sin(2.0*(M_PI_2 - phi))
              + 0.0625*sin(4.0*(M_PI_2 - phi));
            v = 1.0/3.0*sin(3.0*lambda);
//printf("%f %f\n", lambda, phi);
//printf("%f %f\n", u, v);
       
 
        }
    }
    // evaluate the longitudes
return 0;
    dx = (lamMax - lamMin)/(double) (nlon - 1);
lambda = lamMin;
lambda = lamMin;
phi = 0.5;
dx = 0.1;
printf("\n");
    for (ix=0; ix<101; ix++)
    {
        ierr = findNextEastPointForFixedLatitude(lambda, phi, dx, &lambda1);
 if (lambda > M_PI){break;}
        if (ierr == 1){break;}
        projections_hammer_latLonToXY(lambda1, phi, &x, &y);
printf("%f %f\n", x, y); //lambda1, phi);
        lambda = lambda1;
   }

    umin =-M_PI_2;
    umax = M_PI_2;
    vmin =-M_PI/6.0;
    vmax = M_PI/6.0;
    return 0;
}

/*!
 * @brief Computes the next northern point such that y_1 = y_0 + dy maps to 
 *        dy = y(phi_1) - y(phi_0) with unknown phi_1.  In effect this is 
 *        solving:
 *
 *        dy = f(\phi_1) - f(\phi_0)
 *           = \frac{\sqrt{2} \sin \phi_1}
 *                  {\sqrt{1 + \cos \phi_1 \cos \frac{\lambda}{2}}
 *           - \frac{\sqrt{2} \sin \phi_0}
 *                  {\sqrt{1 + \cos \phi_0 \cos \frac{\lambda}{2}}}
 *
 *        Calling the second term on the right hand side c and rearraing
 *        has that 
 *
 *        dy + f = \frac{\sqrt{2} \sin \phi_1}
 *                      {\sqrt{1 + \cos \phi_1 \cos \frac{\lambda}{2}}}
 *
 *        or, after squaring, 
 *
 *        (dy + f)^2 = \frac{2 (1 - \cos^2 \phi_1)}
 *                          {1 + \cos \phi_1 \cos \frac{\lambda}{2}}
 *
 *        Rearranging in the form of a quadratic for \cos \phi_1 has that
 *
 *        0 = \cos^2 \phi_1
 *          + \frac{(dy + f)^2}{2} \cos \frac{\lambda}{2} \cos \phi_1
 *          + \frac{(dy + f)^2 - 2}{2} 
 *
 *        After solving for the roots \pm x = \cos(\phi_1) all that 
 *        remains is to compute \phi_1 = acos(x) where it would appear that
 *        either the positive or negative root cah be chosen. 
 *
 * @param[in] lam0    fixed longitude (radians)
 * @param[in] phi0    latitude (-pi/2, pi,2)
 *
 * @param[in] dy      desired grid spacing (radians)
 *
 * @param[out] phi1   new latitude so that y(phi_1) - y(phi_0) = dy
 *
 * @result -1 indicates that phi0 + dy is greater than pi/2.
 *          0 indicates succcess.
 *          1 indicates a failure to find a suitable phi1
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2 
 */ 
static int findNextNorthPointForFixedLongitude(
    const double lam0,  const double phi0,
    const double dy, double *phi1)
{
    double b, c, cosHalfLam0, det, dyf, dyf2, f, sqrtDet, x, y0, y1;
    const double sqrt2 = 1.4142135623730951;
    //------------------------------------------------------------------------//
    //
    // precompute some values to save a few cycles
    *phi1 = M_PI_2; // choose as the upper limit the most northern point
    cosHalfLam0 = cos(0.5*lam0);
    f = sqrt2*sin(phi0)/(sqrt(1.0 + cos(phi0)*cos(0.5*lam0)));
    dyf = dy + f;
    dyf2 = dyf*dyf;
    // set the b and c terms in the quadratic equation where a = 1
    b = 0.5*dyf2*cosHalfLam0;
    c = 0.5*(dyf2 - 2.0);
    det = b*b - 4.0*c;
    // this means that \phi_1 + \phi_2 > pi/2 i.e. we're out of bounds 
    // return with the maximum pi/2 as defind above
    if (fabs(det) < 1.e-10){det = 0.0;}
    if (det < 0.0){return 1;}
    sqrtDet = sqrt(det);
    *phi1 = acos(0.5*(-b + sqrtDet)); // choose the positive root
    // resolve the sign ambiguity generated when we squared phi_0 + f
    // so that phi_1 is > phi_0.  i feel tighter bounds could be made. 
    if (phi0 + dy < 0.0){*phi1 =-*phi1;}
    // verification step
    y0 = f;//projections_hammer_latLonToXY(lam0,  phi0, &x, &y0);
    projections_hammer_latLonToXY(lam0, *phi1, &x, &y1);
    if (fabs(y1 - y0 - dy) > 1.e-10)
    {
        // try to resolve the sign ambiguity when crossing the equator
        *phi1 =-*phi1;
        projections_hammer_latLonToXY(lam0, *phi1, &x, &y1);
        if (fabs(y1 - y0 - dy) > 1.e-10)
        {
            printf("warning in computing phi_1: %f - %f = %f /= dy=%f\n",
                   y1, y0, y1 - y0, dy);
            return -1;
        }
    }
    return 0;
}
//============================================================================//

/*!
 * @brief Computes the next eastern point such that x_1 = x_0 + dx maps to
 *        dx = x(\lambda_1) - x(\lambda_0) with unknown lambda_1.  This amounts
 *        to solving:
 *
 *        dx = \frac{2 \sqrt{2} \cos \phi \sin \frac{\lambda_1}{2}}
 *                  {\sqrt{1 + \cos \phi \cos \frac{\lambda_1}{2} }} 
 *           - \frac{2 \sqrt{2} \cos \phi \sin \frac{\lambda_0}{2}}
 *                  {\sqrt{1 + \cos \phi \cos \frac{\lambda_0}{2} }}
 *
 *        or
 *
 *        dx + f = \frac{2 \sqrt{2} \cos \phi \sin \frac{\lambda_1}{2}}
 *                      {\sqrt{1 + \cos \phi \cos \frac{\lambda_1}{2} }}
 *
 *        After squaring this becomes
 *
 *        (dx + f)^2 = \frac{8 \cos^2 \phi (1 - \cos^2 \frac{\lambda_1}{2})}
 *                          {1 + \cos \phi \cos \frac{\lambda_1}}{2} }
 *
 *        Using cos^2 + sin^2 = 1 and rearranging as a quadratic equation has
 *        that
 *
 *        0 = \cos^2 \frac{\lambda_1}{2}
 *          + \frac{(dx + f)^2}{8} \cos \phi \cos \frac{\lambda_1}{2}
 *          + \frac{(dx + f)^2 - 8 \cos^2 \phi}{8}
 *
 *        
 *
 */ 
static int findNextEastPointForFixedLatitude(
    const double lam0, const double phi0,
    const double dx, double *lam1)
{
    double b, c, cosPhi, cosPhi2, den1, den2, det, dxf, dxf2, f,
           halfLam, sqrtDet, x0, x1, y;
    const double twoSqrt2 = 2.8284271247461903;
    const double eigth = 0.125;
    //------------------------------------------------------------------------//
    *lam1 = M_PI; // define the maximum point
    cosPhi = cos(phi0);
    cosPhi2 = cosPhi*cosPhi;
    den1 = 1.0/(8.0*cosPhi);
    den2 = 1.0/(8.0*cosPhi2);
    halfLam = 0.5*lam0; 
    f = twoSqrt2*cosPhi*sin(halfLam)/sqrt(1.0 + cosPhi*cos(halfLam));
    dxf = dx + f;
    dxf2 = dxf*dxf;
    // create the quadratic
    b = den1*dxf2;
    c = den2*dxf2 - 1.0;
    det = b*b - 4.0*c;
    // this corresponds to looping around the world
    if (fabs(det) < 1.e-10){det = 0.0;}
    if (det < 0.0){return 1;}
    sqrtDet = sqrt(det);
    *lam1 = 2.0*acos(0.5*(-b+sqrtDet)); // choose larger root
    if (lam0 + dx < 0.0){*lam1 =-*lam1;}
    // verification step
    x0 = f; //projections_hammer_latLonToXY( lam0, phi0, &x0, &y);
    projections_hammer_latLonToXY(*lam1, phi0, &x1, &y);
    if (fabs(x1 - f - dx) > 1.e-10)
    {
        *lam1 =-*lam1;
        projections_hammer_latLonToXY(*lam1, phi0, &x1, &y);
        // cycling around
        if (fabs(x1 - x0 - dx) > 1.e-10)
        {
            return 1;
            //printf("warning in computing lam_1: %f - %f = %f /= dy=%f\n ",
            //       x1, x0, x1 - x0, dx);
            //return -1; 
        }
    }
    return 0; 
}

int lune_findLuneInterpolationPoints(const double umin, const double umax,
                                     const double vmin, const double vmax,
                                     const int nxPixel, const int nyPixel)
{
    const char *fcnm = "lune_findLuneInterpolationPoints\0";
    double dx, dy, lam0, phi0;
    int ix, iy;
    const double sqrt2 = 1.4142135623730951;
    if (nxPixel < 2 || nyPixel < 2)
    {
        if (nxPixel < 2){printf("%s: Insufficient number of x pixels\n", fcnm);}
        if (nyPixel < 2){printf("%s: Insufficient number of y pixels\n", fcnm);}
        return -1;
    }
    // Start at the south-western most point and work up the `west'
    projections_hammer_latLonToXY(vmin, umin, &lam0, &phi0);
    // iterate until the northern most point is reached (but not exceeded)
    dx = 1.0/(double) (nxPixel - 1);
    dy = 1.0/(double) (nyPixel - 1);
    if (umax <= umin || vmax <= vmin)
    {

    }
    // map from cartesian onto unit sphere
    for (iy=0; iy<nyPixel; iy++)
    {
        for (ix=0; ix<nxPixel; ix++)
        {
/*
            x = (double) ix*dx;
            y = (double) iy*dy; 
            // Get the coordinate 
            r = sqrt(x*x + y*y + 1.0);
            theta = arccos(1.0/r);
            phi = atan2(y, x); 
*/
        }
    }

    return 0;
}

/*!
 * @brief Conversion of x and y to longitude and latitude in the Hammer
 *        (equal-area) projection
 *
 * @param[in] x        x location
 * @param[in] y        y location
 *
 * @param[out] lambda  longitude (radians)
 * @param[out] phi     latitude (radians)
 *
 * @url https://en.wikipedia.org/wiki/Hammer_projection
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
void projections_hammer_xyToLatLon(const double x, const double y,
                                   double *lambda, double *phi)
{
    double z;
    // Compute intermediate variable
    z = sqrt(1.0 - 0.0625*(x*x) - 0.25*(y*y));
    // Compute longitude
    *lambda = 2.0*atan2(z*x, 4.0*(z*z) - 2.0);
    *phi = asin(z*y);
    return;
}

/*!
 * @brief Conversion of longitude and latitude to x and y in the
 *        Hammer (equal-area) projection.
 *
 * @param[in] lambda   longitude (radians).  \f$ \lambda \in [-\pi, \pi ] \f$.
 * @param[in] phi      latitude (radians). \f$ \phi \in [0, \pi] \f$.
 *
 * @param[out] x       x corresponding to longitude.
 *                     \f$ x \in [-2 \sqrt{2}, 2 \sqrt{2}] \f$.
 * @param[out] y       y corresponding to latitude.
 *                     \f$ y \in \left [-\sqrt{2], \sqrt{2} ] \right \f$.
 * 
 * @url https://en.wikipedia.org/wiki/Hammer_projection
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
void projections_hammer_latLonToXY(const double lambda, const double phi,
                                   double *x, double *y)
{
    double cosLam2, cosPhi, den, halfLam, sinLam2, sinPhi;
    const double sqrt2 = 1.4142135623730951;
    const double twoSqrt2 = 2.8284271247461903;
    halfLam = 0.5*lambda;
    cosPhi = cos(phi);
    sinPhi = sin(phi);
    sinLam2 = sin(halfLam);
    cosLam2 = cos(halfLam);
    den = 1.0/sqrt(1.0 + cosPhi*cosLam2);
    *x = twoSqrt2*cosPhi*sinLam2*den;
    *y = sqrt2*sinPhi*den;
    return;
}
//============================================================================//
/*!
 * @brief Interpolate function values in rectilinear (v, u) space
 *
 * @param[in] nv     number of v coordinates (> 1) 
 * @param[in] v      monotonic increasing v points which
 *                   correspond to longitudes [nv].  Note that
 *                   \f$ v \in \left [-\frac{1}{3}, \frac{1}{3} \right ] \f$.
 * @param[in] nu     number of u coordinates (> 1)
 * @param[in] u      monotonic increasing u points which
 *                   correspond to colatitudes [nu].  Note that 
 *                   \f$ u \in \left [ 0, \frac{3 \pi}{4} \right ] \f$.
 * @param[in] f      values at f(v,u).  
 * @param[in] nInt   number of interpolation points
 * @param[in] vInt   colatitude-like interpolation points [nInt]
 * @param[in] uInt   longitude-like interpolation points [nInt]
 * @param[out] fInt  interpolated values of f at f(vInt, uInt).
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int lune_interpolateInVUSpace(const int nv, const double *__restrict__ v,
                              const int nu, double *__restrict__ u,
                              const double *__restrict__ f,
                              const int nInt,
                              const double *__restrict__ vInt,
                              const double *__restrict__ uInt,
                              double *__restrict__ fInt)
{
    // nothing to do
    if (nInt < 1){return 0;}
    // input errors
    if (nu < 2 || nv < 2 || u == NULL || v == NULL || 
        vInt == NULL || uInt == NULL || fInt == NULL)
    {
        return -1;
    }
    // verify u and v are sorted

    return 0;
}
//============================================================================//
/*!
 * @brief Convert from \f$ x \in [a,b] \f$ to \f$ \xi \in [c, d] \f$.
 *
 * @param[in] a    lower value in from interval [a,b] s.t. \f$ a \le x \le b \f$
 * @param[in] b    upper bound in from interval [a,b] s.t. \f$ a \le x \le b \f$
 * @param[in] c    lower bound in to interval [c,d] s.t. \f$ c \le \xi \le d \f$
 * @param[in] d    upper bound in to interval [c,d] s.t. \f$ c \le \xi \le d \f$
 *
 * @result corresponding transformed variable of x which now resides in the
 *         new interval \f$ \xi \in [c, d]\f$.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 */
#pragma omp declare simd
static double transform(const double a, const double b,
                        const double c, const double d, const double x)
{
    double c1, c2, det, xi;
    det = 1.0/(b - a);
    c1 = det*(b*c - a*d);
    c2 = det*(d - c);
    xi = c1 + x*c2;
    return xi;
}
