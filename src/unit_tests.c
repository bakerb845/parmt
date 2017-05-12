#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <float.h>
#include <math.h>
#include "sacio.h"
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "parmt_mtsearch.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/random/random.h"
#include "iscl/signal/convolve.h"
#include "iscl/time/time.h"

int testBlockSize(void);

int main()
{
    char *ffList[10] = {"rotation_data/B00103ZDS.sac\0",
                        "rotation_data/B00106ZSS.sac\0",
                        "rotation_data/B00101ZDD.sac\0",
                        "rotation_data/B00109ZEX.sac\0",
                        "rotation_data/B00104RDS.sac\0",
                        "rotation_data/B00107RSS.sac\0",
                        "rotation_data/B00102RDD.sac\0",
                        "rotation_data/B00110REX.sac\0",
                        "rotation_data/B00105TDS.sac\0",
                        "rotation_data/B00108TSS.sac\0"};
    struct parmtNoiseBasis_struct basis;
    struct sacData_struct sacZDS, sacZSS, sacZDD, sacZEX,
                          sacRDS, sacRSS, sacRDD, sacREX,
                          sacTDS, sacTSS, sac;
    double *CeInv, *CeInvDiag, *d, *dn, *G, *R, *X,
           *mts, *betas, *gammas, *kappas, *M0s,
           *phi, *sigmas, *thetas, *xrand, *var, *eye;
    double baz, betaMax, betaMin, gammaMax, gammaMin, kappaMin, kappaMax,
           m0Min, m0Max, sigmaMax, sigmaMin, thetaMax, thetaMin;
    double cmpaz, cmpinc, sigma, xnorm;
    int *lags, i, icomp, ierr, imt, imtopt, j, ldc, ldz, ng, nb, nk,
        nlags, ns, nt, nm, nmt, npgrns, npts, rank;
    const double oneDeg = M_PI/180.0;  // one degree for padding
    const double az = 25.0;
    const int ldm = 8;
    const int blockSize = 32;
    iscl_init();
    // Discretize the moment tensor space
    ng = 11; // longitudes
    nb = 11; // colatitudes
    nk = 11;  // strike angles
    ns = 11;  // rake
    nt = 11;  // dip
    nm = 4;  // magnitudes
    betaMin = 0.0 + oneDeg;
    betaMax = M_PI - oneDeg;
    gammaMin =-M_PI/6.0 + oneDeg;
    gammaMax = M_PI/6.0 - oneDeg;
    kappaMin = 0.0 + oneDeg;
    kappaMax = 2.0*M_PI - oneDeg;
    thetaMin = 0.0 + oneDeg;
    thetaMax = 0.5*M_PI - oneDeg;
    sigmaMin =-0.5*M_PI + oneDeg;
    sigmaMax = 0.5*M_PI - oneDeg;
    m0Min = 1.0/sqrt(2.0);
    m0Max = 2.0/sqrt(2.0);
    //m0Max = 1.e8/sqrt(2.0) + (double) (nm - 1)*1.e8/sqrt(2.0);
    //omp_set_num_threads(nthreads);
    mkl_set_num_threads(1);
    // Discretize the MT space
    M0s    = array_linspace64f(m0Min,    m0Max,    nm, &ierr);
    betas  = array_linspace64f(betaMin,  betaMax,  nb, &ierr);
    gammas = array_linspace64f(gammaMin, gammaMax, ng, &ierr);
    kappas = array_linspace64f(kappaMin, kappaMax, nk, &ierr);
    thetas = array_linspace64f(thetaMin, thetaMax, nt, &ierr);
    sigmas = array_linspace64f(sigmaMin, sigmaMax, ns, &ierr);
    nmt = ng*nb*nk*ns*nt*nm; 
    mts = memory_calloc64f(nmt*ldm);
    ierr = parmt_discretizeMT64f(ng, gammas,
                                 nb, betas,
                                 nm, M0s,
                                 nk, kappas,
                                 nt, thetas,
                                 ns, sigmas,
                                 ldm, nmt, mts);
    if (ierr != 0)
    {
        printf("Failed to discretize MTs\n");
        return EXIT_FAILURE;
    }
    memory_free64f(&M0s);
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&thetas);
    memory_free64f(&sigmas);
    // Read the green's functions generated by CPS
    memset(&sacZDS, 0, sizeof(struct sacData_struct));
    memset(&sacZSS, 0, sizeof(struct sacData_struct));
    memset(&sacZDD, 0, sizeof(struct sacData_struct));
    memset(&sacZEX, 0, sizeof(struct sacData_struct));
    memset(&sacRDS, 0, sizeof(struct sacData_struct));
    memset(&sacRSS, 0, sizeof(struct sacData_struct));
    memset(&sacRDD, 0, sizeof(struct sacData_struct));
    memset(&sacREX, 0, sizeof(struct sacData_struct));
    memset(&sacTDS, 0, sizeof(struct sacData_struct));
    memset(&sacTSS, 0, sizeof(struct sacData_struct));
    ierr  = sacio_readTimeSeriesFile(ffList[0], &sacZDS);
    ierr += sacio_readTimeSeriesFile(ffList[1], &sacZSS);
    ierr += sacio_readTimeSeriesFile(ffList[2], &sacZDD);
    ierr += sacio_readTimeSeriesFile(ffList[3], &sacZEX);
    ierr += sacio_readTimeSeriesFile(ffList[4], &sacRDS);
    ierr += sacio_readTimeSeriesFile(ffList[5], &sacRSS);
    ierr += sacio_readTimeSeriesFile(ffList[6], &sacRDD);
    ierr += sacio_readTimeSeriesFile(ffList[7], &sacREX);
    ierr += sacio_readTimeSeriesFile(ffList[8], &sacTDS);
    ierr += sacio_readTimeSeriesFile(ffList[9], &sacTSS);
    if (ierr != 0)
    {
        printf("Failed to load fundmaental faults\n");
        return EXIT_FAILURE;
    }
    // Set the green's functions
    npgrns = MIN(64, sacZDS.npts);
    npts = npgrns;
    icomp = 1;     // vertical
    cmpaz = 0.0;   // channel azimuth
    cmpinc =-90.0; // channel orientation (-90+ up 0 is horizontal)
    G = memory_calloc64f(6*npgrns);
    printf("Setting Green's functions...\n");
    baz = fmod(az + 180, 360.0);
    ierr = parmt_utils_ff2mtGreens64f(npgrns, icomp,
                                      az, baz,
                                      cmpaz, cmpinc,
                                      sacZDS.data, sacZSS.data,
                                      sacZDD.data, sacZEX.data,
                                      sacRDS.data, sacRSS.data,
                                      sacRDD.data, sacREX.data,
                                      sacTDS.data, sacTSS.data,
                                      &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                                      &G[3*npgrns], &G[4*npgrns], &G[5*npgrns]);
    if (ierr != 0)
    {
        printf("Failed to set Green's functions\n");
        return EXIT_FAILURE;
    }
    // Free intermediate memory
    sacio_free(&sacZDS);
    sacio_free(&sacZSS);
    sacio_free(&sacZDD);
    sacio_free(&sacZEX);
    sacio_free(&sacRDS);
    sacio_free(&sacRSS);
    sacio_free(&sacRDD);
    sacio_free(&sacREX);
    sacio_free(&sacTDS);
    sacio_free(&sacTSS);
    // make a noise matrix 
    ldc = npts;
    CeInv = memory_calloc64f(npts*ldc);
    R = memory_calloc64f(npts*ldc);
    X = memory_calloc64f((2*npts-1)*ldc);
    //xrand = random_randn64f(npts, 0.0, 0.01, &ierr); 
    xrand = random_rand64f(npts, &ierr);
    xrand[0] = 0.0;
    xrand[npts-1] = 0.0;
    xnorm = cblas_dnrm2(npts, xrand, 1);
    cblas_dscal(npts, 1.0/xnorm, xrand, 1);
    ierr = convolve_corrmtx64f_work(npts, xrand,
                                    npts-1, ldc, X,
                                    CORRMTX_AUTOCORRELATION,
                                    true, ldc, R);
FILE *fmat = fopen("xcmat.txt", "w");
for (i=0; i<npts; i++)
{
  for (j=0; j<npts; j++)
  {
   fprintf(fmat, "%d %d %e\n", i+1, npts-j, R[i*ldc+j]);
  }
fprintf(fmat,"\n");
}
fclose(fmat);
    memory_free64f(&X);
    cblas_dcopy(npts*ldc, R, 1, CeInv, 1);
    // Make a noise basis
    ierr = parmt_utils_getNoiseBasisFromCtC64f(npts, ldc, R, &basis);
    if (ierr != 0)
    {
        printf("Failed to compute eigenbasis\n");
        return EXIT_FAILURE;
    }
    // Compute inverse of R
    ierr = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', npts, CeInv, ldc);
    if (ierr != 0)
    {
        printf("Failed to compute chol(CeInv)\n");
        return EXIT_FAILURE;
    } 
    ierr = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', npts, CeInv, ldc);
    if (ierr != 0)
    {
        printf("Failed to compute CeInv %d\n", ierr);
        return EXIT_FAILURE;
    }
    // The KL basis uses sqrt(eig) but for consistency with inv(Ce) need
    // just the regular old eigenvalues
    for (i=0; i<basis.nvals; i++)
    {
        basis.sqrtEvals[i] = pow(basis.sqrtEvals[i], 2);
    }
    // Copy upper to lower
    for (i=0; i<npts; i++)
    {
        for (j=0; j<i; j++)
        {
            CeInv[i*ldc+j] = CeInv[j*ldc+i];
        }
    }
    CeInvDiag = memory_calloc64f(npts);
    for (i=0; i<npts; i++){CeInvDiag[i] = CeInv[ldc*i+i];}
    // Make a synthetic
    imtopt = MIN(nmt-1, MAX(0, (3*nmt)/8));
    d = memory_calloc64f(npts);
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                npgrns, 6, 1.0, G, npgrns,
                &mts[ldm*imtopt], 1, 0.0, d, 1);
    // make a noisy synthetic so Ei doesn't blow up
double relError = 0.2;
rank = parmt_utils_getEffectiveRank(relError, basis);
printf("%d\n", rank);
//    rank = basis.npts/4;
//printf("%d\n", rank);
//getchar();
printf("Effective rank used: %d\n", rank);
    sigma = 1.0;
    double snr = 15.0;
    dn = memory_calloc64f(npts);
    parmt_utils_makeKLNoise(rank, npts, basis, dn);
    double xsum = array_sum64f(npts, dn, &ierr);
    for (i=0; i<npts; i++){dn[i] = dn[i] - xsum/(double) npts;} 
    double signalEnergy = cblas_dnrm2(npts, d, 1); // signal variance
    double noiseEnergy = cblas_dnrm2(npts, dn, 1); // noise variance
printf("%f %f\n", signalEnergy, noiseEnergy);
    double xvar = 0.0;
    double xfact = (1.0/noiseEnergy)*signalEnergy/pow(10.0, 0.05*snr);
    for (i=0; i<npts; i++)
    {
        dn[i] = xfact*dn[i] + d[i];
    }
    sigma = 0.0;
    for (i=0; i<npts; i++)
    {
        sigma = sigma + pow(dn[i] - d[i], 2); // noise variance
    }
    sigma = sqrt(sigma);
    //sigma = pow(cblas_dnrm2(npts, dn, 1), 2);
printf("variance: %f\n", sigma);
//sigma = 0.002;
printf("%e %e %e %e\n", xfact, sigma, d[14], dn[14]); 
    //------------------------------------------------------------------------//
    //                        unlagged conventional                           //
    //------------------------------------------------------------------------//
    // run the inversion for a diagonal matrix
    phi = memory_calloc64f(nmt);
    time_tic();
    nlags = 0;
    ierr = parmt_mtSearch64f(ldm, nlags, nmt,
                             npts, ldc,
                             true,
                             &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                             &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                             CeInvDiag, mts, d, phi);
    if (ierr != 0)
    {
        printf("Failed calling mtSearch64f 1\n");
        return EXIT_FAILURE;
    }
    imt = array_argmin64f(nmt, phi);
    if (imt != imtopt)
    {
        printf("Failed to recover optimal index in test 1\n");
        return EXIT_FAILURE;
    }
    printf("Unlagged diagonal matrix time: %f (s)\n", time_toc());
    time_tic();
 printf("est + true optimal indices: %d %d %e\n", imt, imtopt, phi[imt]);

    ierr = parmt_mtSearch64f(ldm, nlags, nmt,
                             npts, ldc,
                             false,
                             &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                             &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                             CeInv, mts, dn, phi);
    if (ierr != 0)
    {   
        printf("Failed calling mtSearch64f 2\n");
        return EXIT_FAILURE;
    }
    imt = array_argmin64f(nmt, phi);
    if (imt != imtopt)
    {   
        printf("Failed to recover optimal index in test 2 %d %d\n", imt, imtopt);
        //return EXIT_FAILURE;
    }
imt = array_argmin64f(nmt, phi);
    printf("est + true optimal indices: %d %d %e\n", imt, imtopt, phi[imt]);
    printf("Unlagged time: %f (s)\n", time_toc());
    lags = memory_calloc32i(nmt);
    var = memory_calloc64f(npts);
    nlags = 0;
    eye = array_set64f(npts, 1.0, &ierr);
    time_tic();
    ierr = parmt_mtSearchL164f(ldm, nmt, npts, blockSize,
                               nlags, false,
                               &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                               &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                               eye, mts, d, phi, var, lags);
    imt = array_argmax64f(nmt, phi);
    printf("L1 est + true optimal indices: %d %d %e %e\n",
           imt, imtopt, phi[imt], phi[imtopt]);
    printf("Unlagged general matrix time: %f (s)\n", time_toc());
    nlags = 4;
omp_set_num_threads(1);
    ierr = parmt_mtSearchL164f(ldm, nmt, npts, blockSize,
                               nlags, true,
                               &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                               &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                               eye, mts, d, phi, var, lags);
    imt = array_argmax64f(nmt, phi);
    printf("L1 lagged est + true optimal indices: %d %d %d %e %e\n",
           imt, imtopt, lags[imt], phi[imt], phi[imtopt]);
    printf("Lagged time: %f (s)\n", time_toc());
    // repeat with lags 
    ierr = testBlockSize();
/*
    nlags = 10;
    eye = array_set64f(npts, 1.0, &ierr);
    ierr = parmt_mtSearchL164f(ldm, nmt, npts, blockSize,
                               nlags, true,
                               &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                               &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                               eye, mts, d, phi, var, lags);
    imt = array_argmax64f(nmt, phi);
    printf("L1 est + true optimal indices: %d %d %e %d\n", imt, imtopt, phi[imt], lags[imt]);
    printf("Unlagged general matrix time: %f (s)\n", time_toc());
*/
/*
    //------------------------------------------------------------------------//
    //                             unlagged kl                                //
    //------------------------------------------------------------------------//
    time_tic();
    ldz = basis.lde;
    ierr = parmt_mtSearchKL64f(nlags, nmt,
                               npts, ldz,
                               rank,
                               sigma,
                               &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                               &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                               basis.sqrtEvals,
                               basis.evects,
                               mts, dn, phi);
    if (ierr != 0)
    {
        printf("Failed calling mtSearchKL64f unlagged\n");
        return EXIT_FAILURE;
    }
    imt = array_argmax64f(nmt, phi);
*/
/*
    if (imt != imtopt)
    {
        printf("Failed to locate unlagged kl optimum\n");
        return EXIT_FAILURE;
    }
*/
/*
    printf("est + true optimal indices: %d %d %e %e %e\n", imt, imtopt,
           phi[imt], phi[imtopt], phi[array_argmin64f(nmt, phi)]);
    printf("%e\n", array_sum64f(nmt, phi, &ierr));
    printf("Unlagged general kl time %f (s)\n", time_toc());
*/
//TEST:;
/*
double *est = memory_calloc64f(npts);
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                npgrns, 6, 1.0, G, npgrns,
                &mts[ldm*imt], 1, 0.0, est, 1);
FILE *fn = fopen("noise.txt","w");
for (i=0; i<npts; i++)
{
 fprintf(fn, "%d %e %e %e %e\n", i, d[i], dn[i], est[i], d[i] - dn[i]); 
}
fclose(fn);
*/
    // Free memory
    parmt_utils_freeNoiseBasis(&basis);
    memory_free64f(&G);
    memory_free64f(&CeInv);
    memory_free64f(&CeInvDiag);
    memory_free64f(&mts);
    memory_free64f(&d);
    memory_free64f(&phi);
    memory_free64f(&var);
    memory_free64f(&eye);
    iscl_finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
int testPolarityMTSearch(void)
{
    const char *fcnm = "testPolarityMTSearch\0";
    const int blockSizes[4] = {1, 2, 32, 64};
    int ierr, nb, ng, nk, nm, ns, nt;
    double *betas, *gammas, *kappas, *mts, *M0s, *thetas, *sigmas;
    double betaMin, betaMax, gammaMin, gammaMax, kappaMin, kappaMax,
           m0Min, m0Max, sigmaMin, sigmaMax, thetaMin, thetaMax;
    const double oneDeg = M_PI/180.0;  // one degree for padding
    // Discretize the moment tensor space
    ng = 5; // longitudes
    nb = 6; // colatitudes
    nk = 7; // strike angles
    ns = 8; // rake
    nt = 9; // dip
    nm = 1; // magnitudes
    betaMin = 0.0 + oneDeg;
    betaMax = M_PI - oneDeg;
    gammaMin =-M_PI/6.0 + oneDeg;
    gammaMax = M_PI/6.0 - oneDeg;
    kappaMin = 0.0 + oneDeg;
    kappaMax = 2.0*M_PI - oneDeg;
    thetaMin = 0.0 + oneDeg;
    thetaMax = 0.5*M_PI - oneDeg;
    sigmaMin =-0.5*M_PI + oneDeg;
    sigmaMax = 0.5*M_PI - oneDeg;
    m0Min = 1.0/sqrt(2.0);
    m0Max = 2.0/sqrt(2.0);
    ierr = parmt_discretizeCells64f(nb, betaMin, betaMax,
                                    ng, gammaMin, gammaMax,
                                    nk, kappaMin, kappaMax,
                                    ns, sigmaMin, sigmaMax,
                                    nt, thetaMin, thetaMax,
                                    nm, m0Min, m0Max,
                                    &betas, &gammas, &kappas,
                                    &sigmas, &thetas, &M0s); 
    if (ierr != 0)
    {
        printf("%s: Error making cell discretization\n", fcnm);
        return EXIT_FAILURE;
    }
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&sigmas);
    memory_free64f(&thetas);
    memory_free64f(&M0s);
    // Investigate different block sizes
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Tests the effect of blocksize in gemm on performance
 */ 
int testBlockSize(void)
{
    const char *fcnm = "testBlockSize\0";
    FILE *fout;
    char fname[PATH_MAX];
    double *G, *d, *eye, *phi, *phiOri, *var, time, dnorm, xmin;
    int *lags = NULL;
    int i, ierr, imt, nmt, nt0, npmax, npts, npgrns, pSize;
    const int nblocks = 10;
    const int nps = 6;
    int pSizes[6] = {32, 64, 128, 256, 512, 1024};
    int blockSizes[10] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512};
    const int nlags = 0;
    const int ldm = 6;
    const bool lwantLags = false;
    int imtopt, ng, nb, nk, ns, nt, nm;
    double *betas, *gammas, *kappas, *mts, *M0s, *thetas, *sigmas;
    double betaMin, betaMax, gammaMin, gammaMax, kappaMin, kappaMax,
           m0Min, m0Max, sigmaMin, sigmaMax, thetaMin, thetaMax;
    const double oneDeg = M_PI/180.0;  // one degree for padding
    // disretize the moment tensor
    // Discretize the moment tensor space
    ng = 5; // longitudes
    nb = 5; // colatitudes
    nk = 7; // strike angles
    ns = 7; // rake
    nt = 7; // dip
    nm = 1; // magnitudes
    betaMin = 0.0 + oneDeg;
    betaMax = M_PI - oneDeg;
    gammaMin =-M_PI/6.0 + oneDeg;
    gammaMax = M_PI/6.0 - oneDeg;
    kappaMin = 0.0 + oneDeg;
    kappaMax = 2.0*M_PI - oneDeg;
    thetaMin = 0.0 + oneDeg;
    thetaMax = 0.5*M_PI - oneDeg;
    sigmaMin =-0.5*M_PI + oneDeg;
    sigmaMax = 0.5*M_PI - oneDeg;
    m0Min = 1.0/sqrt(2.0);
    m0Max = 2.0/sqrt(2.0);
    // Discretize the MT space
    M0s    = array_linspace64f(m0Min,    m0Max,    nm, &ierr);
    betas  = array_linspace64f(betaMin,  betaMax,  nb, &ierr);
    gammas = array_linspace64f(gammaMin, gammaMax, ng, &ierr);
    kappas = array_linspace64f(kappaMin, kappaMax, nk, &ierr);
    thetas = array_linspace64f(thetaMin, thetaMax, nt, &ierr);
    sigmas = array_linspace64f(sigmaMin, sigmaMax, ns, &ierr);
    nmt = ng*nb*nk*ns*nt*nm; 
    imtopt = nmt/2;
    mts = memory_calloc64f(nmt*ldm);
    ierr = parmt_discretizeMT64f(ng, gammas,
                                 nb, betas,
                                 nm, M0s,
                                 nk, kappas,
                                 nt, thetas,
                                 ns, sigmas,
                                 ldm, nmt, mts);
    if (ierr != 0)
    {   
        printf("Failed to discretize MTs\n");
        return EXIT_FAILURE;
    }   
    memory_free64f(&M0s);
    memory_free64f(&betas);
    memory_free64f(&gammas);
    memory_free64f(&kappas);
    memory_free64f(&thetas);
    memory_free64f(&sigmas);
    // make some random data for testing
    npmax = array_max32i(nps, pSizes); 
    G = random_rand64f(6*npmax, &ierr);
    // set space
    phiOri = memory_calloc64f(nmt);
    phi = memory_calloc64f(nmt);
    var = memory_calloc64f(npmax);
    d = memory_calloc64f(npmax);
    eye = array_set64f(npmax, 1.0, &ierr);
    // fix the number of threads
    nt0 = omp_get_num_threads();
    omp_set_num_threads(1);
    // loop on point sizes
    for (pSize=0; pSize<nps; pSize++)
    {
        // make synthetic data
        memset(fname, 0, PATH_MAX*sizeof(char));
        sprintf(fname, "blockSize_%d.txt", pSizes[pSize]);
        fout = fopen(fname, "w");
        npts = pSizes[pSize];
        npgrns = npts;
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    npgrns, 6, 1.0, G, npgrns, &mts[ldm*imtopt], 1, 0.0, d, 1);
        for (i=0; i<nblocks; i++)
        {
            time_tic();
            ierr = parmt_mtSearchL164f(ldm, nmt, npts, blockSizes[i],
                                       nlags, lwantLags, 
                                       &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                                       &G[3*npgrns], &G[4*npgrns], &G[5*npgrns],
                                       eye, mts, d, phi, var, lags);
            time = time_toc();
            fprintf(fout, "%d %d %e %e %d\n",
                    npts, blockSizes[i], 1.e6*time/(double) nmt, time/(double) nmt, nmt);
            if (i == 0)
            {
                ierr = array_copy64f_work(nmt, phi, phiOri);
            }
            dnorm = array_normDiff64f(nmt, phi, phiOri, ONE_NORM, 1.0, &ierr);
            xmin = array_min64f(nmt, phi);
            if (dnorm > DBL_EPSILON*(double) nmt)
            {
                printf("Error computing blocksize\n");
                ierr = 1; 
                goto ERROR;
            }
            imt = array_argmax64f(nmt, phi);
            if (imt != imtopt)
            {
                printf("Failed to locate optimum\n");
                ierr = 1;
                goto ERROR;
            }
        }
        fclose(fout);
    }
ERROR:;
    // restore number of threads
    omp_set_num_threads(nt0);
    // clear memory
    memory_free64f(&d);
    memory_free64f(&phiOri);
    memory_free64f(&phi);
    memory_free64f(&eye);
    memory_free64f(&var);
    memory_free64f(&mts);
    memory_free64f(&G);
    return 0;
}
