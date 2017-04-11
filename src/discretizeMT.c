#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mpi.h>
#include "compearth.h"
#include "parmt_mtsearch.h"
#ifdef PARMT_USE_INTEL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

/*
int main_temp()
{
    double *betas, *gammas, *kappas, *M0s, *mts, *sigmas, *thetas;
    int i, ierr, nmt;
    const int ng = 29;
    const int nb = 29;
    const int nk = 6;
    const int ns = 6;
    const int nt = 6;
    const int nm = 1;
    const int ldm = 8;
    const double oneDeg = M_PI/180.0;
    const double betaMin = 0.0 + oneDeg;
    const double betaMax = M_PI - oneDeg;
    const double gammaMin =-M_PI/6.0 + oneDeg;
    const double gammaMax = M_PI/6.0 - oneDeg;
    const double kappaMin = 0.0 + oneDeg;
    const double kappaMax = 2.0*M_PI - oneDeg;
    const double thetaMin = 0.0 + oneDeg;
    const double thetaMax = 0.5*M_PI - oneDeg;
    const double sigmaMin =-0.5*M_PI + oneDeg;
    const double sigmaMax = 0.5*M_PI - oneDeg;
    nmt = nm*nb*ng*nk*ns*nt; 
    mts = memory_calloc64f(nmt*ldm); //(double *)calloc((size_t) (nmt*ldm), sizeof(double)); 
    M0s = memory_calloc64f(nm); //(double *)calloc((size_t) nm, sizeof(double)); 
    betas = memory_calloc64f(nb); //(double *)calloc((size_t) nb, sizeof(double));
    gammas = memory_calloc64f(ng); //(double *)calloc((size_t) ng, sizeof(double));
    kappas = memory_calloc64f(nk); //(double *)calloc((size_t) nk, sizeof(double));
    sigmas = memory_calloc64f(ns); //(double *)calloc((size_t) ns, sizeof(double));
    thetas = memory_calloc64f(nt); //(double *)calloc((size_t) nt, sizeof(double));

    ierr = array_linspace64f_work(1.0/sqrt(2.0), 1.0/sqrt(2.0), nm, M0s);
    ierr = array_linspace64f_work(betaMin+0.001,  betaMax-0.001,  nb, betas); 
    ierr = array_linspace64f_work(gammaMin, gammaMax, ng, gammas);
    ierr = array_linspace64f_work(kappaMin, kappaMax, nk, kappas);
    ierr = array_linspace64f_work(thetaMin, thetaMax, nt, thetas);
    ierr = array_linspace64f_work(sigmaMin, sigmaMax, ns, sigmas);
    // avoid warnings on theta in tt2cmt
    ierr = gridSearch_discretizeMT64f(ng, gammas, nb, betas,  nm, M0s,
                                      nk, kappas, nt, thetas, ns, sigmas,
                                      ldm, nmt, mts);

double Mdc[6], lam[3], U[6], theta;
    compearth_tt2cmt(0.0, 0.0, 1./sqrt(2.),
                     0.0, 0.1, 0.0, Mdc, lam, U);

    FILE *ofl = fopen("mts.txt", "w");
    for (i=0; i<nmt; i++) 
    {
        compearth_angleMT(1, Mdc, &mts[ldm*i], &theta);

        fprintf(ofl, "%f %f %f %f %f %f %f\n", //i+1,
                mts[ldm*i+0], mts[ldm*i+1], mts[ldm*i+2], 
                mts[ldm*i+3], mts[ldm*i+4], mts[ldm*i+5], theta);
    }
    fclose(ofl);

    free(mts);
    free(M0s);
    free(betas);
    free(gammas);
    free(kappas);
    free(sigmas);
    free(thetas); 
    return 0;
}
*/
//============================================================================//
/*!
 * @brief Computes the moment tensors in the grid search. 
 *        The grid search ordering is:
 *         Loop on magnitudes
 *          Loop on colatitudes
 *           Loop on longitudes
 *            Loop on strikes
 *             Loop on rakes
 *              Loop on dips 
 *
 * @param[in] comm    MPI communicator on which moment tensor will be split
 * @param[in] ng      number of longitudes
 * @param[in] gammas  longitudes (radians) on lune s.t. 
 *                    \f$ \gamma \in [-\pi/6, \pi/6] \f$ [ng]
 * @param[in] nb      number of colatitudes
 * @param[in] betas   colatitudes (radians) on lune s.t. 
 *                    \f$ \beta \in [0, \pi] \f$ [nb]
 * @param[in] nm      number of scalar moments
 * @param[in] M0s     scalar moments in Newton-meters [nm]. 
 *                    for a `unit' moment tensor choose M0 = 1/sqrt(2)
 * @param[in] nk      number of strike angles
 * @param[in] kappas  strike angles (radians) s.t.
 *                    \f$ \kappa \in [0, 2\pi] \f$ [nk]
 * @param[in] nt      number of dip angles
 * @param[in] thetas  dips angles (radians) s.t.
 *                    \f$ \theta \in [0, \pi/2] \f$ [nt]
 * @param[in] ns      number of rake angles
 * @param[in] sigmas  rake angles (radians) s.t. 
 *                    \f$ \sigma \in [-\pi/2, \pi/2] \f$ [ns]
 * @param[in] ldm     leading dimension of mts.  must be >= 6.
 * @param[in] nmt     number of moment tensors (=nm*nb*ng*nk*ns*nt)
 *
 * @param[out] mts    moment tensors (Newton-meters) in north,east,down s.t.
 *                    \f$ \textbf{m}
 *                      = \f$ \{m_{xx}, m_{yy}, m_{zz}, 
 *                              m_{xy}, m_{xz}, m_{yz}\} \f$.
 *                    This has dimensions [nm*nb*ng*nk*ns*nt x ldm].  Observe
 *                    the ordering as defind above.
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI licensed under Apached 2
 *
 */
int parmt_discretizeMT64f_MPI(const MPI_Comm comm,
                              const int ng,
                              const double *__restrict__ gammas,
                              const int nb,
                              const double *__restrict__ betas,
                              const int nm,
                              const double *__restrict__ M0s,
                              const int nk,
                              const double *__restrict__ kappas,
                              const int nt,
                              const double *__restrict__ thetas,
                              const int ns,
                              const double *__restrict__ sigmas,
                              const int ldm,
                              struct localMT_struct *mts)
{
    const char *fcnm = "parmt_discretizeMT64f_MPI\0";
    double lam[3], Muse[6], U[9] __attribute__ ((aligned (64)));
    double *bloc, *gloc, *kloc, *mloc, *sloc, *tloc, *mtWork;
    int64_t nmt64;
    int *nmtProc, *offset;
    int dmt, i, ib, ib0, ierr, ierr1, ig, ig0, ik, ik0, im, im0, is, is0,
        it, it0, imt, imtAll, imtloc, imt1, imt2, jmt, myid, nbloc, ngloc,
        nkloc, nmloc, nmtall, nmt, nmtWork, nprocs, nsloc, ntloc;
    bool *lbloc, *lgloc, *lkloc, *lmloc, *lsloc, *ltloc;
    const int master = 0;
    const double pi180i = 180.0/M_PI;
    //------------------------------------------------------------------------//
    ierr = 0;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    memset(mts, 0, sizeof(struct localMT_struct));
    mts->comm = comm;
    // Check for integer overflow
    nmt64 = (int64_t) ng * (int64_t) nb * (int64_t) nm
           *(int64_t) nk * (int64_t) nt * (int64_t) ns; 
    if (nmt64 > INT_MAX)
    {   
        printf("%s: Error integer overflow - gridsearch too large\n", fcnm);
        return -1; 
    }   
    nmt = (int) nmt64;
    // Verify the inputs
    if (ldm < 6 || ng < 1 || nb < 1 || nm < 1 || nk < 1 || nt < 1 || ns < 1 ||
        gammas == NULL || betas == NULL || M0s == NULL || kappas == NULL ||
        thetas == NULL || sigmas == NULL)
    {
        if (ldm < 6){printf("%s: Invalid leading dimension\n", fcnm);}
        if (ng < 1){printf("%s: ng must be positive\n", fcnm);}
        if (nb < 1){printf("%s: nb must be positive\n", fcnm);}
        if (nm < 1){printf("%s: nm must be positive\n", fcnm);}
        if (nk < 1){printf("%s: nk must be positive\n", fcnm);}
        if (nt < 1){printf("%s: nt must be positive\n", fcnm);}
        if (ns < 1){printf("%s: ns must be positive\n", fcnm);}
        if (gammas == NULL){printf("%s: gammas is NULL\n", fcnm);}
        if (betas == NULL){printf("%s: betas is NULL\n", fcnm);}
        if (M0s == NULL){printf("%s: M0s is NULL\n", fcnm);}
        if (kappas == NULL){printf("%s: kappas is NULL\n", fcnm);}
        if (thetas == NULL){printf("%s: thetas is NULL\n", fcnm);}
        if (sigmas == NULL){printf("%s: sigmas is NULL\n", fcnm);}
        return -1;
    }
    // divide the moment tensor grid
    dmt = MAX(nmt/nprocs, 1);
    imt1 = myid*dmt;
    imt2 = (myid + 1)*dmt;
    if (myid == nprocs - 1){imt2 = nmt;}
    mts->nmt = imt2 - imt1;
    mts->ldm = ldm;
    MPI_Allreduce(&mts->nmt, &nmtall, 1, MPI_INTEGER, MPI_SUM, mts->comm);
    if (nmtall != nmt)
    {
        if (myid == master)
        {
            printf("%s: Failed to partition domain %d %d\n", fcnm, nmtall, nmt);
        }
        return -1;
    }
    // create the requisite information for an MPI_GATHER
    mts->nmtProc = memory_calloc32i(nprocs);
    mts->offset = memory_calloc32i(nprocs);
    nmtProc = memory_calloc32i(nprocs);
    offset = memory_calloc32i(nprocs);
    nmtProc[myid] = mts->nmt; 
    offset[myid] = imt1;
    MPI_Allreduce(nmtProc, mts->nmtProc, nprocs, MPI_INTEGER,
                  MPI_SUM, mts->comm);
    MPI_Allreduce(offset, mts->offset,  nprocs, MPI_INTEGER,
                  MPI_SUM, mts->comm); 
    mts->commSize = nprocs;
    mts->nmtAll = nmtall;
    mts->myid = myid;
    memory_free32i(&nmtProc);
    memory_free32i(&offset);
    // finall do the division on each process
    mts->l2g = memory_calloc32i(MAX(1, mts->nmt));
    lmloc = memory_calloc8l(nm);
    lbloc = memory_calloc8l(nb);
    lgloc = memory_calloc8l(ng);
    lsloc = memory_calloc8l(ns);
    lkloc = memory_calloc8l(nk);
    ltloc = memory_calloc8l(nt);
    im0 = INT_MAX;
    ib0 = INT_MAX;
    ig0 = INT_MAX;
    is0 = INT_MAX;
    ik0 = INT_MAX;
    it0 = INT_MAX;
    // get my local moment tensor terms 
    imtloc = 0;
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
                            // index
                            imt = im*nb*ng*nk*ns*nt
                                + ib*ng*nk*ns*nt
                                + ig*nk*ns*nt
                                + ik*ns*nt
                                + is*nt
                                + it; 
                            // verify it is in bounds
                            if (imt >= imt1 && imt < imt2)
                            {
                                im0 = MIN(im0, im);
                                ib0 = MIN(ib0, ib);
                                ig0 = MIN(ig0, ig);
                                is0 = MIN(is0, is);
                                ik0 = MIN(ik0, ik);
                                it0 = MIN(it0, it);
                                lmloc[im] = true;
                                lbloc[ib] = true;
                                lgloc[ig] = true;
                                lsloc[is] = true;
                                lkloc[ik] = true;
                                ltloc[it] = true;
                                // update the pointers
                                mts->l2g[imtloc] = imt;
                                imtloc = imtloc + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_Allreduce(&imtloc, &imtAll, 1, MPI_INT, MPI_SUM, mts->comm);
    if (imtAll != mts->nmtAll)
    {
        printf("%s: Error imtAll != nmtAll\n", fcnm);
        return -1;
    }
    // Compute a superset of my grid
    // get the local elements of the moment tensor grid search
    mloc = array_mask64f(nm, M0s,    lmloc, &nmloc, &ierr);
    bloc = array_mask64f(nb, betas,  lbloc, &nbloc, &ierr);
    gloc = array_mask64f(ng, gammas, lgloc, &ngloc, &ierr);
    kloc = array_mask64f(nk, kappas, lkloc, &nkloc, &ierr);
    sloc = array_mask64f(ns, sigmas, lsloc, &nsloc, &ierr);
    tloc = array_mask64f(nt, thetas, ltloc, &ntloc, &ierr);
    nmtWork = MAX(1, nmloc)*MAX(1, nbloc)*MAX(1, ngloc)
             *MAX(1, nkloc)*MAX(1, nsloc)*MAX(1, ntloc);
    mtWork = memory_calloc64f(mts->ldm*MAX(1, nmtWork));
    // discretize the moment tensor space
    ierr = parmt_discretizeMT64f(ngloc, gloc,
                                 nbloc, bloc,
                                 nmloc, mloc,
                                 nkloc, kloc,
                                 ntloc, tloc,
                                 nsloc, sloc,
                                 mts->ldm, nmtWork, mtWork);
    if (ierr != 0)
    {
        printf("%s: Error discretizing MT on process %d\n", fcnm, myid);
        return -1;
    }
    mts->mts = memory_calloc64f(mts->ldm*MAX(1, mts->nmt));
    // Extract my moment tensors from superset 
    imtloc = 0;
    jmt =-1;
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
                            // index
                            imt = im*nb*ng*nk*ns*nt
                                + ib*ng*nk*ns*nt
                                + ig*nk*ns*nt
                                + ik*ns*nt
                                + is*nt
                                + it;
                            jmt = (im - im0)*nbloc*ngloc*nkloc*nsloc*ntloc
                                + (ib - ib0)*ngloc*nkloc*nsloc*ntloc
                                + (ig - ig0)*nkloc*nsloc*ntloc
                                + (ik - ik0)*nsloc*ntloc
                                + (is - is0)*ntloc
                                + (it - it0);
                            // verify it is in bounds
                            if (imt >= imt1 && imt < imt2)
                            {
                                if (jmt < 0 || jmt >= nmtWork)
                                {
                                    printf("%s: Sizing error %d %d %d\n",
                                           fcnm, mts->myid, jmt, nmtWork);
                                    ierr = 1;
                                }
                                for (i=0; i<6; i++)
                                {
                                    mts->mts[mts->ldm*imtloc+i]
                                         = mtWork[mts->ldm*jmt+i];
                                }
                                imtloc = imtloc + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // free memory
    memory_free8l(&lmloc);
    memory_free8l(&lbloc);
    memory_free8l(&lgloc);
    memory_free8l(&lkloc);
    memory_free8l(&lsloc);
    memory_free8l(&ltloc);
    memory_free64f(&mloc);
    memory_free64f(&bloc);
    memory_free64f(&gloc);
    memory_free64f(&kloc);
    memory_free64f(&sloc);
    memory_free64f(&tloc);
    memory_free64f(&mtWork);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the moment tensors in the grid search. 
 *        The grid search ordering is:
 *         Loop on magnitudes
 *          Loop on colatitudes
 *           Loop on longitudes
 *            Loop on strikes
 *             Loop on rakes
 *              Loop on dips 
 *
 * @param[in] ng      number of longitudes
 * @param[in] gammas  longitudes (radians) on lune s.t. 
 *                    \f$ \gamma \in [-\pi/6, \pi/6] \f$ [ng]
 * @param[in] nb      number of colatitudes
 * @param[in] betas   colatitudes (radians) on lune s.t. 
 *                    \f$ \beta \in [0, \pi] \f$ [nb]
 * @param[in] nm      number of scalar moments
 * @param[in] M0s     scalar moments in Newton-meters [nm]. 
 *                    for a `unit' moment tensor choose M0 = 1/sqrt(2)
 * @param[in] nk      number of strike angles
 * @param[in] kappas  strike angles (radians) s.t.
 *                    \f$ \kappa \in [0, 2\pi] \f$ [nk]
 * @param[in] nt      number of dip angles
 * @param[in] thetas  dips angles (radians) s.t.
 *                    \f$ \theta \in [0, \pi/2] \f$ [nt]
 * @param[in] ns      number of rake angles
 * @param[in] sigmas  rake angles (radians) s.t. 
 *                    \f$ \sigma \in [-\pi/2, \pi/2] \f$ [ns]
 * @param[in] ldm     leading dimension of mts.  must be >= 6.
 * @param[in] nmt     number of moment tensors (=nm*nb*ng*nk*ns*nt)
 *
 * @param[out] mts    moment tensors (Newton-meters) in north,east,down s.t.
 *                    \f$ \textbf{m}
 *                      = \f$ \{m_{xx}, m_{yy}, m_{zz}, 
 *                              m_{xy}, m_{xz}, m_{yz}\} \f$.
 *                    This has dimensions [nm*nb*ng*nk*ns*nt x ldm].  Observe
 *                    the ordering as defind above.
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI licensed under Apached 2
 *
 */
int parmt_discretizeMT64f(const int ng,
                          const double *__restrict__ gammas,
                          const int nb,
                          const double *__restrict__ betas,
                          const int nm,
                          const double *__restrict__ M0s,
                          const int nk,
                          const double *__restrict__ kappas,
                          const int nt,
                          const double *__restrict__ thetas,
                          const int ns,
                          const double *__restrict__ sigmas,
                          const int ldm,
                          const int nmt, double *__restrict__ mts)
{
    const char *fcnm = "parmt_discretizeMT64f\0";
    double lam[3], Muse[6], U[9] __attribute__ ((aligned (64)));
    double deltaDeg, gammaDeg, kappaDeg, M0, sigmaDeg, thetaDeg;
    int i, ierr, ierr1, ib, ig, im, imt, ik, is, it;
    const double pi180i = 180.0/M_PI;
    const double betaMin = 0.0;
    const double betaMax = M_PI;
    const double gammaMin =-M_PI/6.0;
    const double gammaMax = M_PI/6.0;
    const double kappaMin = 0.0; 
    const double kappaMax = 2.0*M_PI;
    const double thetaMin = 0.0;
    const double thetaMax = M_PI_2;
    const double sigmaMin =-M_PI_2;
    const double sigmaMax = M_PI_2;
    //------------------------------------------------------------------------//
    ierr = 0;
    // Verify the inputs
    if (ldm < 6 || ng < 1 || nb < 1 || nm < 1 || nk < 1 || nt < 1 || ns < 1 ||
        nmt != ng*nb*nm*nk*nt*ns ||
        gammas == NULL || betas == NULL || M0s == NULL || kappas == NULL ||
        thetas == NULL || sigmas == NULL)
    {
        if (ldm < 6){printf("%s: Invalid leading dimension\n", fcnm);}
        if (ng < 1){printf("%s: ng must be positive\n", fcnm);}
        if (nb < 1){printf("%s: nb must be positive\n", fcnm);}
        if (nm < 1){printf("%s: nm must be positive\n", fcnm);}
        if (nk < 1){printf("%s: nk must be positive\n", fcnm);}
        if (nt < 1){printf("%s: nt must be positive\n", fcnm);}
        if (ns < 1){printf("%s: ns must be positive\n", fcnm);}
        if (nmt != ng*nb*nm*nk*nt*ns)
        {
            printf("%s: nmt != ng*nb*nm*nk*nt*ns %d %d\n",
                   fcnm, nmt, ng*nb*nm*nk*nt*ns);
        }
        if (gammas == NULL){printf("%s: gammas is NULL\n", fcnm);}
        if (betas == NULL){printf("%s: betas is NULL\n", fcnm);}
        if (M0s == NULL){printf("%s: M0s is NULL\n", fcnm);}
        if (kappas == NULL){printf("%s: kappas is NULL\n", fcnm);}
        if (thetas == NULL){printf("%s: thetas is NULL\n", fcnm);}
        if (sigmas == NULL){printf("%s: sigmas is NULL\n", fcnm);}
        return -1;
    }
    for (i=0; i<ng; i++)
    {
        if (gammas[i] < gammaMin || gammas[i] > gammaMax)
        {
            printf("%s: gammas %f out or range [%f,%f]\n",
                   fcnm, gammas[i], gammaMin, gammaMax); 
            return -1;
        }
    }
    for (i=0; i<nb; i++)
    {   
        if (betas[i] < betaMin || betas[i] > betaMax)
        {
            printf("%s: betas %f out or range [%f,%f]\n",
                   fcnm, betas[i], betaMin, betaMax);
            return -1;
        }
    }
    for (i=0; i<nk; i++)
    {   
        if (kappas[i] < kappaMin || kappas[i] > kappaMax)
        {
            printf("%s: kappas %f out or range [%f,%f]\n",
                   fcnm, kappas[i], kappaMin, kappaMax);
            return -1;
        }
    }
    for (i=0; i<nt; i++)
    {   
        if (thetas[i] < thetaMin || thetas[i] > thetaMax)
        {
            printf("%s: thetas %f out or range [%f,%f]\n",
                   fcnm, thetas[i], thetaMin, thetaMax);
            return -1;
        }
    }
    for (i=0; i<ns; i++)
    {   
        if (sigmas[i] < sigmaMin || sigmas[i] > sigmaMax)
        {
            printf("%s: sigmas %f out or range [%f,%f]\n",
                   fcnm, sigmas[i], sigmaMin, sigmaMax);
            return -1;
        }
    }
    // Loop on scalar moments
    #pragma omp parallel for collapse(6) \
     firstprivate (lam, Muse, U) \
     private (deltaDeg, gammaDeg, kappaDeg, M0, \
              sigmaDeg, thetaDeg, ierr1, im, imt, ib, ig, ik, is, it) \
     shared (fcnm, M0s, betas, gammas, kappas, sigmas, mts, thetas) \
     reduction (max:ierr) \
     default (none)
    for (im=0; im<nm; im++)
    {
        // Loop on colatitudes
        for (ib=0; ib<nb; ib++)
        {
            // Loop on longitudes
            for (ig=0; ig<ng; ig++)
            {
                // Loop on strike
                for (ik=0; ik<nk; ik++)
                {
                    // Loop on rake
                    for (is=0; is<ns; is++)
                    {
                        // Loop on dip
                        for (it=0; it<nt; it++)
                        {
                            // tape**2 space term
                            M0 = M0s[im];
                            deltaDeg = (M_PI/2.0 - betas[ib])*pi180i;
                            gammaDeg = gammas[ig]*pi180i;
                            kappaDeg = kappas[ik]*pi180i;
                            sigmaDeg = sigmas[is]*pi180i;
                            thetaDeg = thetas[it]*pi180i;
                            // index
                            imt = im*nb*ng*nk*ns*nt
                                + ib*ng*nk*ns*nt
                                + ig*nk*ns*nt
                                + ik*ns*nt
                                + is*nt
                                + it;
                            // Compute the corresponding moment tensor
                            ierr1 = compearth_tt2cmt(gammaDeg, deltaDeg, M0,
                                                     kappaDeg, thetaDeg,
                                                     sigmaDeg,
                                                     Muse, lam, U);
                            if (ierr1 != 0)
                            {
                                printf("%s: Error calling tt2cmt %d\n",
                                       fcnm, ierr1);
                                ierr = ierr + 1;
                                continue;
                            }
                            ierr1 = compearth_convertMT(USE, NED, Muse,
                                                        &mts[imt*ldm]);
                            if (ierr1 != 0)
                            {
                                printf("%s: Error changing coords\n", fcnm);
                                ierr = ierr + 1;
                                continue;
                            }
                        } // loop on dip
                    } // loop on rake
                } // loop on strike
            } // loop on longitudes
        } // loop on latitudes
    } // loop on scalar moments 
    return 0;
}
