#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "parmt_postProcess.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

static double *diff(const int n, const double *__restrict__ x, int *ierr);
static int computePadding8i(const int n);

int marginal_getOptimum(const int nloc, const int nm, const int nb,
                        const int ng, const int nk, const int ns, const int nt,
                        const double *__restrict__ phi,
                        int *jloc, int *jm, int *jb, int *jg,
                        int *jk, int *js, int *jt)
{
    int ib, ig, ik, iloc, imax, im, imt, indx, is, it, nmt;
    enum isclError_enum ierr;
    nmt = nm*nb*ng*nk*ns*nt;
    imax = array_argmax64f(nmt*nloc, phi, &ierr);
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            for (ib=0; ib<nb; ib++)
            {
                for (ig=0; ig<ng; ig++)
                {
                    for (ik=0; ik<nk; ik++)
                    {
                        for (is=0; is<ns; is++)
                        {
                            for (it=0; it<nt; it++)
                            {
                                imt = im*nb*ng*nk*ns*nt
                                    + ib*ng*nk*ns*nt
                                    + ig*nk*ns*nt
                                    + ik*ns*nt
                                    + is*nt
                                    + it;
                                indx = iloc*nmt + imt;
                                if (indx == imax)
                                {
                                    *jloc = iloc;
                                    *jm = im;
                                    *jb = ib;
                                    *jg = ig;
                                    *jk = ik;
                                    *js = is;
                                    *jt = it;
                                }
                            }
                        }
                    }
                }
            }
         }
     }
     return 0;
}
/*!
 * @brief Computes the marginal PDF at lune points using Eqn 48 of
 *        Tape and Tape 2015.
 *
 * @param[in] nloc       number of locations in grid-search.  
 * @param[in] nm         number of moment magnitudes in grid search
 * @param[in] nb         number of colatitudes
 * @param[in] betas      colatitudes (radians) of lune points [nb]
 * @param[in] ng         number of longitudes
 * @param[in] gammas     longitudes (radians) of lune points [ng]
 * @param[in] nk         number of strike angles
 * @param[in] kappas     fault strike angles (radians) [nk]
 * @param[in] ns         number of slip angles
 * @param[in] sigmas     fault slip angles (radians) [ns]
 * @param[in] nh         number of dip angles
 * @param[in] h          related to dip angles (radians) [nh]
 * @param[in] phi        objective function computed throughout grid-search
 *                       space to marginalize [nloc*nm*nb*ng*nk*ns*nh]
 * @param[out] luneMPDF  marginalized lune [nb x ng].  the ib, ig'th 
 *                       colatitude and longitude is given by ib*ng + ig.
 *   
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 * @bug There is inadequate handling of moment tensor scaling and locs
  *     is a misnomer.  Technically this can only loop in 1D and I would
 *      require the cell volumes for position to loop on locations.
 */
int marginal_computeLuneMPDF(const int nloc,
                             const int nm,
                             const int nb, const double *__restrict__ betas,
                             const int ng, const double *__restrict__ gammas,
                             const int nk, const double *__restrict__ kappas,
                             const int ns, const double *__restrict__ sigmas,
                             const int nh, const double *__restrict__ h,
                             const double *__restrict__ phi,
                             double *__restrict__ luneMPDF)
{
    const char *fcnm = "marginal_computeLuneMPDF\0";
    int ib, ierr, ig, igIb, ih, ik, iloc, im, imt, indx, is, nmt;
    double *cos3g, *dh, *dk, *ds, *twoSinB3;
    double sinb, sinb3, dkdsdh, dV;
    nmt = nm*nb*ng*nk*ns*nh;
    memset(luneMPDF, 0, (size_t) (ng*nb)*sizeof(double));
    if (nk < 2 || ns < 2 || nh < 2)
    {
        if (nk < 2){printf("%s: Error nk must be > 1\n", fcnm);}
        if (ns < 2){printf("%s: Error ns must be > 1\n", fcnm);}
        if (nh < 2){printf("%s: Error nh must be > 1\n", fcnm);}
    }
    // Differentiate weights for Riemann quadrature
    dk = diff(nk, kappas, &ierr);
    ds = diff(ns, sigmas, &ierr);
    dh = diff(nh, h, &ierr);
    // theta goes from 0 to pi/2 as h goes from 1 to 0 hence it's 
    // orientation is reversed from our ini file (recall that it asks
    // for theta and not h) and hence -dh is not conducive to quadrature.
    cblas_dscal(nh-1, -1.0, dh, 1);
    // Precompute the trigonemetric terms
    twoSinB3 = memory_calloc64f(nb);
    for (ib=0; ib<nb; ib++)
    {
        twoSinB3[ib] = 2.0*pow(sin(betas[ib]), 3);
    }
    cos3g = memory_calloc64f(ng);
    for (ig=0; ig<ng; ig++)
    {
        cos3g[ig] = cos(3.0*gammas[ig]);
    }
    // Marginalize the lune 
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            for (ib=0; ib<nb; ib++)
            {
                for (ig=0; ig<ng; ig++)
                {
                    for (ik=0; ik<nk-1; ik++)
                    {
                        for (is=0; is<ns-1; is++)
                        {
                            for (ih=0; ih<nh-1; ih++)
                            {
                                imt = im*nb*ng*nk*ns*nh
                                    + ib*ng*nk*ns*nh
                                    + ig*nk*ns*nh
                                    + ik*ns*nh
                                    + is*nh
                                    + ih;
                                igIb = ib*ng + ig;
                                indx = iloc*nmt + imt;
                                //sinb = sin(betas[ib]);
                                //sinb3 = (sinb*sinb)*sinb;
                                //twoSinb3 = 2.0*sinb3;
                                //cos3g = cos(3.0*gammas[ig]);
                                dkdsdh = (dk[ik]*ds[is])*dh[ih];
                                dV = twoSinB3[ib]*cos3g[ig]*dkdsdh;
                                luneMPDF[igIb] = luneMPDF[igIb] + phi[indx]*dV;
                            }
                        }
                    }
                }
            }
        }
    }
    memory_free64f(&dh);
    memory_free64f(&dk);
    memory_free64f(&ds);
    memory_free64f(&twoSinB3);
    memory_free64f(&cos3g);
    return 0;
}

int marginal_computeLuneUVMPDF(const int nloc,
                               const int nm, 
                               const int nu, const double *__restrict__ us,
                               const int nv, const double *__restrict__ vs,
                               const int nk, const double *__restrict__ kappas,
                               const int ns, const double *__restrict__ sigmas,
                               const int nh, const double *__restrict__ hs,
                               const double *__restrict__ phi,
                               double *__restrict__ luneMPDF)
{
    const char *fcnm = "marginal_computeLuneMPDF\0";
    int ierr, ih, ik, iloc, im, imt, indx, is, iu, iuIv, iv, nmt;
    double *dh, *dk, *ds, *du, *dv;
    double dudv, dkdsdh, dV; 
    const double sqrt8 = sqrt(8.0);
    nmt = nm*nu*nv*nk*ns*nh;
    memset(luneMPDF, 0, (size_t) (nu*nv)*sizeof(double));
    if (nk < 2 || ns < 2 || nh < 2)
    {   
        if (nk < 2){printf("%s: Error nk must be > 1\n", fcnm);}
        if (ns < 2){printf("%s: Error ns must be > 1\n", fcnm);}
        if (nh < 2){printf("%s: Error nh must be > 1\n", fcnm);}
    }   
    // Differentiate weights for Riemann quadrature
    dk = diff(nk, kappas, &ierr);
    ds = diff(ns, sigmas, &ierr);
    dh = diff(nh, hs,     &ierr);
/*
    du = diff(nu, us,     &ierr);
    dv = diff(nv, vs,     &ierr);
for (int i=0; i<nu; i++)
{
  printf("du: %f: \n", du[i]);
}
for (int i=0; i<nv; i++)
{
 printf("dv: %f: \n", dv[i]);
}
*/
    // theta goes from 0 to pi/2 as h goes from 1 to 0 hence it's 
    // orientation is reversed from our ini file (recall that it asks
    // for theta and not h) and hence -dh is not conducive to quadrature.
    cblas_dscal(nh-1, -1.0, dh, 1); 
    // Marginalize the lune 
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            for (iu=0; iu<nu; iu++)
            {
                for (iv=0; iv<nv; iv++)
                {   
                    for (ik=0; ik<nk-1; ik++)
                    {   
                        for (is=0; is<ns-1; is++)
                        {   
                            for (ih=0; ih<nh-1; ih++)
                            {
                                imt = im*nu*nv*nk*ns*nh
                                    + iu*nv*nk*ns*nh
                                    + iv*nk*ns*nh
                                    + ik*ns*nh
                                    + is*nh
                                    + ih;
                                iuIv = iv*nu + iu;
                                indx = iloc*nmt + imt;
                                //sinb = sin(betas[ib]);
                                //sinb3 = (sinb*sinb)*sinb;
                                //twoSinb3 = 2.0*sinb3;
                                //cos3g = cos(3.0*gammas[ig]);
                                dudv = sqrt8; //du[iu]*dv[iv];
                                dkdsdh = (dk[ik]*ds[is])*dh[ih];
                                dV = dudv*dkdsdh;
                                luneMPDF[iuIv] = luneMPDF[iuIv] + phi[indx]*dV;
                            }
                        }
                    }
                }
            }
        }
    }
/*
    memory_free64f(&du);
    memory_free64f(&dv);
*/
    memory_free64f(&dk);
    memory_free64f(&ds);
    memory_free64f(&dh);
    return 0;
}
//============================================================================//
int marginal_computeDepthMPDF(const int nloc,
                              const int nm, const double *__restrict__ M0s,
                              const int nb, const double *__restrict__ betas,
                              const int ng, const double *__restrict__ gammas,
                              const int nk, const double *__restrict__ kappas,
                              const int ns, const double *__restrict__ sigmas,
                              const int nt, const double *__restrict__ thetas,
                              const double *__restrict__ phi,
                              double *__restrict__ depMagMPDF,
                              double *__restrict__ depMPDF)
{
    const char *fcnm = "marginal_computeDepthMPDF\0";
    double *cos3g, *db,  *dg, *dk, *dm, *ds, *dt, *twoSinB4, *sint,
           arg, dbdg, dkdsdt, dV5, geom, sum;
    int ib, ierr, ig, ik, iloc, im, imt, indx, is, it, jloc;
    ierr = 0;
    memset(depMPDF, 0, (size_t) nloc*sizeof(double));
    memset(depMagMPDF, 0, (size_t) (nloc*nm)*sizeof(double));
    // Cell spacings for for Riemann quadrature
    dm = memory_calloc64f(nm);
    db = memory_calloc64f(nb);
    dg = memory_calloc64f(ng);
    dk = memory_calloc64f(nk);
    ds = memory_calloc64f(ns);
    dt = memory_calloc64f(nt);
    ierr += postprocess_computeBetaCellSpacing(nb,  betas, db);
    ierr += postprocess_computeGammaCellSpacing(ng, gammas, dg);
    ierr += postprocess_computeKappaCellSpacing(nk, kappas, dk);
    ierr += postprocess_computeSigmaCellSpacing(ns, sigmas, ds);
    ierr += postprocess_computeThetaCellSpacing(nt, thetas, dt);
    ierr += postprocess_computeM0CellSpacing(nm, M0s, dm);
printf("%s: for now i'm setting dm = 1 and dl = 1\n", fcnm);
array_set64f_work(nm, 1.0, dm);
    if (ierr != 0)
    {
        printf("%s: Failed to compute cell spacing\n", fcnm);
        return -1;
    }
    if (array_min64f(nm, dm, &ierr) <= 0.0 ||
        array_min64f(nb, db, &ierr) <= 0.0 ||
        array_min64f(ng, dg, &ierr) <= 0.0 ||
        array_min64f(nk, dk, &ierr) <= 0.0 ||
        array_min64f(ns, ds, &ierr) <= 0.0 ||
        array_min64f(nt, dt, &ierr) <= 0.0)
    {
        printf("%s: Negative jacobian\n", fcnm);
        return -1;
    }
    // Compute the geometric factors
    twoSinB4 = memory_calloc64f(nb);
    for (ib=0; ib<nb; ib++)
    {
        twoSinB4[ib] = 2.0*pow(sin(betas[ib]), 4);
    }
    sint = memory_calloc64f(nt);
    for (it=0; it<nt; it++)
    {
        sint[it] = sin(thetas[it]);
    }
    cos3g = memory_calloc64f(ng);
    for (ig=0; ig<ng; ig++)
    {
        cos3g[ig] = cos(3.0*gammas[ig]);
    }
    if (array_min64f(nb, twoSinB4, &ierr) <= 0.0 ||
        array_min64f(nt, sint, &ierr) <= 0.0 ||
        array_min64f(ng, cos3g, &ierr) <= 0.0)
    {
        printf("%s: Warning negative jacobian from geometric factors\n", fcnm);
    }
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        { 
            for (ib=0; ib<nb; ib++)
            { 
                for (ig=0; ig<ng; ig++)
                { 
                    for (ik=0; ik<nk; ik++)
                    {
                        for (is=0; is<ns; is++)
                        {
                            for (it=0; it<nt; it++)
                            {
                                imt = iloc*nm*nb*ng*nk*ns*nt
                                    +      im*nb*ng*nk*ns*nt
                                    +         ib*ng*nk*ns*nt
                                    +            ig*nk*ns*nt
                                    +               ik*ns*nt
                                    +                  is*nt
                                    +                     it;
                                dkdsdt = (dk[ik]*ds[is])*dt[it];
                                dbdg = db[ib]*dg[ig];
                                geom = (twoSinB4[ib]*cos3g[ig])*sint[it];
                                dV5 = geom*(dbdg*dkdsdt);
                                //if (dV5 < 0.0){printf("error\n");}
                                arg = phi[imt]*dV5;
                                jloc = im*nloc + iloc;
                                depMagMPDF[jloc] = depMagMPDF[jloc] + arg;
                            }
                        }
                    }
                }
            }
        }
    }
    // Integrate out the depMagMPDF
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            jloc = im*nloc + iloc;
            depMPDF[iloc] = depMPDF[iloc] +  depMagMPDF[jloc]*dm[im];
        }
    }
    memory_free64f(&cos3g);
    memory_free64f(&sint);
    memory_free64f(&twoSinB4);
    memory_free64f(&db);
    memory_free64f(&dg);
    memory_free64f(&dk);
    memory_free64f(&ds);
    memory_free64f(&dt);
    memory_free64f(&dm);
    return ierr;
}
//============================================================================//
double marginal_computeNormalization(
    const int nloc, const double *__restrict__ deps,
    const int nm, const double *__restrict__ M0s,
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas,
    const double *__restrict__ phi,
    int *ierr)
{
    const char *fcnm = "marginal_computeNormalization\0";
    double *cos3g, *db,  *dg, *dk, *dl, *dm, *ds, *dt, *twoSinB4, *sint,
           dbdg, dkdsdt, dV2, dV5, dV7, geom, sum;
    int ib, ig, ik, iloc, im, imt, is, it, nmt;
    *ierr = 0;
    sum = 0.0;
    nmt = nm*nb*ng*nk*ns*nt;
    // Cell spacings for for Riemann quadrature
    dl = memory_calloc64f(nloc);
    dm = memory_calloc64f(nm);
    db = memory_calloc64f(nb);
    dg = memory_calloc64f(ng);
    dk = memory_calloc64f(nk);
    ds = memory_calloc64f(ns);
    dt = memory_calloc64f(nt);
    *ierr += postprocess_computeBetaCellSpacing(nb,  betas, db);
    *ierr += postprocess_computeGammaCellSpacing(ng, gammas, dg);
    *ierr += postprocess_computeKappaCellSpacing(nk, kappas, dk);
    *ierr += postprocess_computeSigmaCellSpacing(ns, sigmas, ds);
    *ierr += postprocess_computeThetaCellSpacing(nt, thetas, dt);
    *ierr += postprocess_computeM0CellSpacing(nm, M0s, dm);
    *ierr += postprocess_computeM0CellSpacing(nloc, deps, dl); 
printf("%s: for now i'm setting dm = 1 and dl = 1\n", fcnm);
array_set64f_work(nm, 1.0, dm);
array_set64f_work(nloc, 1.0, dl);
    if (*ierr != 0)
    {
        printf("%s; Error computing normalization\n", fcnm);
        return sum;
    }
    if (array_min64f(nloc, dl, ierr) <= 0.0 ||
        array_min64f(nm, dm, ierr) <= 0.0 ||
        array_min64f(nb, db, ierr) <= 0.0 || 
        array_min64f(ng, dg, ierr) <= 0.0 ||
        array_min64f(nk, dk, ierr) <= 0.0 ||
        array_min64f(ns, ds, ierr) <= 0.0 ||
        array_min64f(nt, dt, ierr) <= 0.0)
    {
        printf("%s: Negative jacobian\n", fcnm);
        *ierr = 1;
        return sum;
    } 
    // Compute the geometric factors
    twoSinB4 = memory_calloc64f(nb);
    for (ib=0; ib<nb; ib++)
    {
        twoSinB4[ib] = 2.0*pow(sin(betas[ib]), 4);
    }
    sint = memory_calloc64f(nt);
    for (it=0; it<nt; it++)
    {
        sint[it] = sin(thetas[it]);
    }
    cos3g = memory_calloc64f(ng);
    for (ig=0; ig<ng; ig++)
    {
        cos3g[ig] = cos(3.0*gammas[ig]);
    }
    sum = 0.0;
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            for (ib=0; ib<nb; ib++)
            {
                for (ig=0; ig<ng; ig++)
                {
                    for (ik=0; ik<nk; ik++)
                    {
                        for (is=0; is<ns; is++)
                        {
                            for (it=0; it<nt; it++)
                            {
                                imt = iloc*nm*nb*ng*nk*ns*nt
                                    +      im*nb*ng*nk*ns*nt
                                    +         ib*ng*nk*ns*nt
                                    +            ig*nk*ns*nt
                                    +               ik*ns*nt
                                    +                  is*nt
                                    +                     it;
                                dkdsdt = (dk[ik]*ds[is])*dt[it];
                                dbdg = db[ib]*dg[ig];
                                geom = (twoSinB4[ib]*cos3g[ig])*sint[it];
                                dV5 = geom*(dbdg*dkdsdt);
                                dV2 = dl[iloc]*dm[im];   
                                dV7 = dV5*dV2;
                                sum = sum + phi[imt]*dV7;
                            }
                        }
                    }
                }
            }
        }
    }
    memory_free64f(&cos3g);
    memory_free64f(&sint);
    memory_free64f(&twoSinB4);
    memory_free64f(&db);
    memory_free64f(&dg);
    memory_free64f(&dk);
    memory_free64f(&ds);
    memory_free64f(&dt);
    memory_free64f(&dl);
    memory_free64f(&dm);
    return sum;
}
//============================================================================//
int marginal_computeMarginalBeachball(
    const int nloc,
    const int nm, 
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas,
    const double *__restrict__ phi
    )
{
    const char *fcnm = "marginal_computeMarginalBeachball\0";
    const double M0loc = 1.0/sqrt(2.0);
    double pAxis[3], nAxis[3], tAxis[3];
    double *bb, *bbAvg, *bbWt, *cos3g, *db,  *dg, *dk, *ds, *dt,
           *twoSinB4, *sint, *xw1, *yw1, dbdg, dkdsdt, dV5, geom, xscal;
    int ib, ig, ik, ierr, iloc, im, imt, indx, is, it, j, jmt, jndx,
        ldi, nmt, nmtBB;
    int8_t *bmap;
    size_t nwork;
    const int nxp = 51;
    const int nyp = nxp;
    const double xc = 1.5;
    const double yc = xc;
    const double rad = 1.0;
    nmt = nm*nb*ng*nk*ns*nt;
    nmtBB = nb*ng*nk*ns*nt;
    ldi = nxp*nyp + computePadding8i(nxp*nyp);
printf("%d %d\n", ldi, nxp*nyp);
    nwork = (size_t)(nmtBB)*(size_t) (ldi);
    if (nwork > INT_MAX)
    {
        printf("%s: Insufficient space for bb\n", fcnm);
        return -1;
    }
    xw1 = memory_calloc64f(nxp*nyp);
    yw1 = memory_calloc64f(nxp*nyp);
    bmap = (int8_t *) calloc(nwork, sizeof(int8_t));
    bbAvg = memory_calloc64f(nxp*nyp);
    bbWt = memory_calloc64f(nxp*nyp);
    bb = memory_calloc64f(nxp*nyp);
    printf("%s: Drawing...\n", fcnm);
    for (ib=0; ib<nb; ib++)
    {
        printf("%s: Drawing beta: %d\n", fcnm, ib+1);
        for (ig=0; ig<ng; ig++)
        {
            for (ik=0; ik<nk; ik++)
            {
                for (is=0; is<ns; is++)
                {
                    for (it=0; it<nt; it++)
                    {
                        postprocess_tt2tnp(betas[ib], gammas[ig], kappas[ik],
                                           sigmas[is], thetas[it],
                                           pAxis, nAxis, tAxis);
                        jmt = ib*ng*nk*ns*nt
                            + ig*nk*ns*nt
                            + ik*ns*nt
                            + is*nt
                            + it;
                        jndx = jmt*ldi;
                        postprocess_tnp2beachballPolarity(nxp, xc, yc, rad,
                                                       pAxis, nAxis, tAxis,
                                                       xw1, yw1, &bmap[jndx]);
                    }
                }
            }
        }
    }
    // Differentiate weights for Riemann quadrature
    db = diff(nb, betas,  &ierr);
    dg = diff(ng, gammas, &ierr);
    dk = diff(nk, kappas, &ierr);
    ds = diff(ns, sigmas, &ierr);
    dt = diff(nt, thetas, &ierr);
    // Compute the geometric factors
    twoSinB4 = memory_calloc64f(nb);
    for (ib=0; ib<nb; ib++)
    {
        twoSinB4[ib] = 2.0*pow(sin(betas[ib]), 4);
    }
    sint = memory_calloc64f(nt);
    for (it=0; it<nt; it++)
    {
        sint[it] = sin(thetas[it]);
    }
    cos3g = memory_calloc64f(ng);
    for (ig=0; ig<ng; ig++)
    {
        cos3g[ig] = cos(3.0*gammas[ig]);
    }
printf("integrating\n");
    // Stack the result
    for (iloc=0; iloc<nloc; iloc++)
    {
        for (im=0; im<nm; im++)
        {
            for (ib=0; ib<nb; ib++)
            {
                for (ig=0; ig<ng; ig++)
                {
                    for (ik=0; ik<nk; ik++)
                    {
                        for (is=0; is<ns; is++)
                        {
                            for (it=0; it<nt; it++)
                            {
                                imt = im*nb*ng*nk*ns*nt
                                    + ib*ng*nk*ns*nt
                                    + ig*nk*ns*nt
                                    + ik*ns*nt
                                    + is*nt
                                    + it; 
                                indx = iloc*nmt + imt;
                                dkdsdt = (dk[ik]*ds[is])*dt[it];
                                dbdg = db[ib]*dg[ig];
                                geom = (twoSinB4[ib]*cos3g[ig])*sint[it];
                                dV5 = geom*(dbdg*dkdsdt); 
//if (dV5 < 0.0){printf("error %f %f %f %f %f %f %f\n", dV5, dbdg, dkdsdt, geom, twoSinB4[ib], cos3g[ig], sint[it]);}
                                xscal = (phi[indx] - 1);//*dV5;
                                jmt = ib*ng*nk*ns*nt
                                    + ig*nk*ns*nt
                                    + ik*ns*nt
                                    + is*nt
                                    + it;
                                jndx = jmt*ldi;
                                for (j=0; j<nxp*nxp; j++)
                                {
                                    bb[j] = bb[j] + xscal*(double) bmap[jndx+j];
//                                    bbAvg[j] = bbAvg[j] + (double) bmap[jndx+j];
//                                    bbWt[j] = bbWt[j] + (double) bmap[jndx+j]*dV5;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
printf("dumping result\n");
int ix, iy;
    FILE *fwork;
    fwork = fopen("depmag/beachball.txt", "w");
    for (iy=0; iy<nyp; iy++)
    {
        for (ix=0; ix<nxp; ix++)
        {
            int k = iy*nxp + ix;
            fprintf(fwork, "%e %e %e %e %e %d\n", xw1[k], yw1[k], bb[k], bbAvg[k], bbWt[k], bmap[nmtBB/2*ldi+k]);
        }
        fprintf(fwork, "\n");
    }
    fclose(fwork);
    memory_free64f(&bbAvg);
    memory_free64f(&bbWt);
    memory_free64f(&xw1);
    memory_free64f(&yw1);
    memory_free64f(&cos3g);
    memory_free64f(&sint);
    memory_free64f(&twoSinB4);
    memory_free64f(&db);
    memory_free64f(&dg);
    memory_free64f(&dk);
    memory_free64f(&ds);
    memory_free64f(&dt);
    free(bmap);
    return 0;
}

static double *diff(const int n, const double *__restrict__ x, int *ierr)
{
    double *d;
    int i;
    d = array_set64f(n, 1.0, ierr);
    for (i=0; i<n; i++)
    {
        if (i < n - 1)
        {
            d[i] = x[i+1] - x[i];
        }
        else
        {
            d[i] = x[i] - x[i-1];
        }
    }
    return d;
}

/*! Related to beta */
/*
static double g11(void)
{
    return 1;
}
*/

/*! Related to gamma */
/*
static double g22(const double sinBeta)
{
    return sinBeta*sinBeta;
}
*/

/*! Related to kappa */
/*
static double g33( )
{
    const double twoSqrt3 = 3.4641016151377544;
    coss2 = coss*coss;
    g33 = 0.5*((4.0 + (-2.0 + 3.0*(1.0 - cos2t))*coss2)*cos2g
        + twoSqrt3*sin2g*sins*sin2t)*sin2b;
}
*/

/*! Related to sigma */
/*
static void g44( )
{
    for (ib=0; ib<nb; ib++)
    {
        for (ig=0; ig<ng; ig++)
        {
            g33 = (2.0 - cos(2.0*gamma))*sinb2;
            for (ik=0; ik<nk; ik++)
            { 
                for (is=0; is<ns; is++)
                {
                    for (it=0; it<nt; it++)
                    {

                    }
                }
            }
        }
    }
    return;
}
*/

/*!
 * @brief Computes the determinant of the upper 3 x 3 matrix 
 */
/*
static void computeDet123(const int nmt,
                          const double *__restrict__ g11,
                          const double *__restrict__ g22,
                          const double *__restrict__ g33,
                          double *__restrict__ det)
{
    int i;
    for (i=0; i<nmt; i++)
    {
        det[i] = (g11[i]*g22[i])*g33;
    }
    return 
}

static void computeDet124(const int nmt,
                          const double *__restrict__ g11,
                          const double *__restrict__ g22,
                          const double *__restrict__ g44,
                          double *__restrict__ det)
{
    computeDet123(nmt, g11, g22, g44, det);
    return;
}

static void computeDet125(const int nmt,
                          const double *__restrict__ g11,
                          const double *__restrict__ g22,
                          const double *__restrict__ g55,
                          double *__restrict__ det)
{
    computeDet123(nmt, g11, g22, g55, det);
    return;
}
*/



int marginal_write1DToGnuplot(const char *fname,
                              const int nx, const double *__restrict__ xlocs,
                              const double *__restrict__ vals)
{
    FILE *fout;
    int ix;
    fout = fopen(fname, "w");
    for (ix=0; ix<nx; ix++)
    {
        fprintf(fout, "%.8e %.10e\n", xlocs[ix], vals[ix]);
    }
    fclose(fout);
    return 0;
}

int marginal_write2DToGnuplot(const char *fname,
                              const int nx, const double *__restrict__ xlocs,
                              const int ny, const double *__restrict__ ylocs,
                              const double *__restrict__ vals)
{
    FILE *fout;
    int ix, ixy, iy;
    fout = fopen(fname, "w");
    for (iy=0; iy<ny; iy++)
    {
        for (ix=0; ix<nx; ix++)
        {
            ixy = iy*nx + ix;
            fprintf(fout, "%.8e %.8e %.10e\n", xlocs[ix], ylocs[iy], vals[ixy]);
        }
        fprintf(fout, "\n");
    } 
    fclose(fout);
    return 0;
}
//============================================================================//
static int computePadding8i(const int n)
{
    size_t mod, pad;
    int ipad;
    // Set space and make G matrix
    pad = 0;
    mod = ((size_t) n*sizeof(int8_t))%64;
    if (mod != 0)
    {
        pad = (64 - mod)/sizeof(int8_t);
    }
    ipad = (int) pad;
    return ipad;
}
