#ifndef PARMT_MTSEARCH_H__
#define PARMT_MTSEARCH_H__ 1
#include <stdbool.h>
#include <mpi.h>
#include "parmt_config.h"
#include "parmt_utils.h"

struct localMT_struct
{
    double *mts;   /*!< Local moment tensors [ldm*nmt] */
    //int *l2g;      /*!< Maps local to global moment tensor indices */
    int *nmtProc;  /*!< Number of moment tensors on each process in
                        communicator [commSize] */
    int *offset;   /*!< Process' global start index offset [commSize]  */
    MPI_Comm comm; /*!< MPI communicator */
    int nmtAll;    /*!< Total number of moment tensors */
    int nmt;       /*!< Number of local moment tensors */
    int ldm;       /*!< Leading dimension of mts */
    int commSize;  /*!< Communicator size */
    int myid;      /*!< my process ID on communicator */
};

#ifdef __cplusplus
extern "C"
{
#endif

//============================================================================//
//                                MPI communication                           //
//============================================================================//
int parmt_broadcast_parmtPolarityParms(struct parmtPolarityParms_struct *parms,
                                        const int root, const MPI_Comm comm);
int parmt_broadcast_mtSearchParms(struct parmtMtSearchParms_struct *parms,
                                  const int root, const MPI_Comm comm); 
int parmt_broadcast_data(struct parmtData_struct *data,
                         const int root, const MPI_Comm comm);
int parmt_broadcast_generalParms(struct parmtGeneralParms_struct *parms,
                                 const int root, const MPI_Comm comm);
int parmt_splitComm(const MPI_Comm globalComm,
                    const int npInObsGroups,
                    const int npInLocGroups,
                    const int npInMTGroups,
                    bool *linObsComm, MPI_Comm *obsComm,
                    bool *linLocComm, MPI_Comm *locComm,
                    bool *linMTComm,  MPI_Comm *mtComm);

//============================================================================//
//                             compute synthetics                             //
//============================================================================//
int parmt_computeSynthetic(const int iobs,
                           const int ilocOpt,
                           const int imtOpt,
                           const int ldm,
                           const double *__restrict__ mts,
                           struct parmtData_struct *data);

//============================================================================//
//                        discretize moment tensor space                      //
//============================================================================//
/*! Cell centers in grid-search */
int parmt_discretizeCells64f(
    const int nb, const double betaLower, const double betaUpper,
    const int ng, const double gammaLower, const double gammaUpper,
    const int nk, const double kappaLower, const double kappaUpper,
    const int ns, const double sigmaLower, const double sigmaUpper,
    const int nt, const double thetaLower, const double thetaUpper,
    const int nm, const double m0Lower, const double m0Upper,
    double **betas,  double **gammas, double **kappas,
    double **sigmas, double **thetas, double **M0s);
/*! Discretize the moment tensor space with MPI */
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
                              struct localMT_struct *mts);
/*! Discretize the moment tensor space */
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
                          const int nmt, double *__restrict__ mts);
//============================================================================//
//                            loop on observations                            //
//============================================================================//
int parmt_obsSearch64f(const MPI_Comm globalComm,
                       const MPI_Comm obsComm, const MPI_Comm locComm,
                       const bool linObsComm, const bool linLocComm,
                       const struct parmtGeneralParms_struct parms,
                       struct parmtData_struct data,
                       struct localMT_struct mtloc,
                       double *__restrict__ phi);
//============================================================================//
//                            location grid search                            //
//============================================================================//
int parmt_locSearchL164f(const MPI_Comm locComm,
                         const int iobs, const int blockSize,
                         const int nlags, const bool lwantLags,
                         struct localMT_struct mtloc,
                         struct parmtData_struct *data,
                         const double *__restrict__ CeInv,
                         double *__restrict__ phi,
                         double *__restrict__ var,
                         int *__restrict__ lags);

int parmt_locSearchXC64f(const MPI_Comm locComm,
                         const int iobs, const int nblock,
                         const int nlags, const bool lwantLags,
                         struct localMT_struct mtloc,
                         struct parmtData_struct *data,
                         double *__restrict__ phi, int *__restrict__ lags);
//============================================================================//
//                             objective function                             //
//============================================================================//
/*! log(sum(exp(-x))) */
double mtsearch_logSumExp64f(const int n, const double *__restrict__ x);
/*! Cross-correlation based objective function */
int parmt_computeStackedCrossCorrelation_MPI(
    const MPI_Comm comm,
    const int npgrns,
    const double *__restrict__ Gxx,
    const double *__restrict__ Gyy,
    const double *__restrict__ Gzz,
    const double *__restrict__ Gxy,
    const double *__restrict__ Gxz,
    const double *__restrict__ Gyz,
    const int ldm, const int nmt, const double *__restrict__ mt, 
    const int npts, const double *__restrict__ data,
    const int lxc, double *__restrict__ xc);
int parmt_computeStackedCrossCorrelation(
    const int npgrns,
    const double *__restrict__ Gxx,
    const double *__restrict__ Gyy,
    const double *__restrict__ Gzz,
    const double *__restrict__ Gxy,
    const double *__restrict__ Gxz,
    const double *__restrict__ Gyz,
    const int ldm, const int nmt, const double *__restrict__ mt,
    const int npts, const double *__restrict__ data,
    const int lxc, double *__restrict__ xc);

/*! Moment tensor search driver */
int parmt_mtSearch64f(const int ldm, const int nlags, const int nmt, 
                      const int npts, const int ldc, 
                      const bool ldiag,
                      const double *__restrict__ Gxx, 
                      const double *__restrict__ Gyy, 
                      const double *__restrict__ Gzz, 
                      const double *__restrict__ Gxy, 
                      const double *__restrict__ Gxz, 
                      const double *__restrict__ Gyz, 
                      const double *__restrict__ CeInv,
                      const double *__restrict__ mts, 
                      const double *__restrict__ d,
                      double *__restrict__ phi);
int parmt_mtSearch32f(const int ldm, const int nlags, const int nmt, 
                      const int npts, const int ldc, 
                      const bool ldiag,
                      const float *__restrict__ Gxx, 
                      const float *__restrict__ Gyy, 
                      const float *__restrict__ Gzz, 
                      const float *__restrict__ Gxy, 
                      const float *__restrict__ Gxz, 
                      const float *__restrict__ Gyz, 
                      const float *__restrict__ CeInv,
                      const float *__restrict__ mts, 
                      const float *__restrict__ d,
                      float *__restrict__ phi);
int parmt_mtSearchKL64f(const int ldm, const int nlags, const int nmt, 
                        const int npts, const int ldz, 
                        const int rank,
                        const double sigma,
                        const double *__restrict__ Gxx, 
                        const double *__restrict__ Gyy, 
                        const double *__restrict__ Gzz, 
                        const double *__restrict__ Gxy, 
                        const double *__restrict__ Gxz, 
                        const double *__restrict__ Gyz, 
                        const double *__restrict__ sqrtEig,
                        const double *__restrict__ Z,
                        const double *__restrict__ mts, 
                        const double *__restrict__ d,
                        double *__restrict__ phi);
/*! Cross-correlation */
int parmt_mtSearchXC64f(const int ldm, const int nmt,
                        const int npts, const int nblock,
                        const int maxlag, const bool lwantLag,
                        const double *__restrict__ Gxx,
                        const double *__restrict__ Gyy,
                        const double *__restrict__ Gzz,
                        const double *__restrict__ Gxy,
                        const double *__restrict__ Gxz,
                        const double *__restrict__ Gyz,
                        const double *__restrict__ mts,
                        const double *__restrict__ d,
                        double *__restrict__ phi,
                        int *__restrict__ lags);
/*! L1 distance */
int parmt_mtSearchL164f(const int ldm, const int nmt,
                        const int npts, const int blockSize,
                        const int nlags, const bool lwantLags,
                        const double *__restrict__ Gxx, 
                        const double *__restrict__ Gyy, 
                        const double *__restrict__ Gzz, 
                        const double *__restrict__ Gxy, 
                        const double *__restrict__ Gxz, 
                        const double *__restrict__ Gyz, 
                        const double *__restrict__ CeInv,
                        const double *__restrict__ mts, 
                        const double *__restrict__ d,
                        double *__restrict__ phi,
                        double *__restrict__ var,
                        int *__restrict__ lags);
/*! Mahalanobis distance */
int parmt_mtSearch_mahalanobis64f(const int ldm, const int nmt, const int npts,
                                  const int ldc, const bool ldiag,
                                  const double *__restrict__ G,
                                  const double *__restrict__ CeInv,
                                  const double *__restrict__ mts, 
                                  const double *__restrict__ d,
                                  double *__restrict__ phi);
int parmt_mtSearch_mahalanobis32f(const int ldm, const int nmt, const int npts,
                                  const int ldc, const bool ldiag,
                                  const float *__restrict__ G,
                                  const float *__restrict__ CeInv,
                                  const float *__restrict__ mts,
                                  const float *__restrict__ d,
                                  float *__restrict__ phi);
/*! Mahalanobis distance with time lags */
int parmt_mtSearch_shiftedMahalanobis64f(const int nlags, const int nmt,
                                         const int npts, const int ldc,
                                         const bool ldiag,
                                         const double *__restrict__ G,
                                         const double *__restrict__ CeInv,
                                         const double *__restrict__ mts,
                                         const double *__restrict__ d,
                                         double *__restrict__ phi);
int parmt_mtSearch_shiftedMahalanobis32f(const int nlags, const int nmt,
                                         const int npts, const int ldc,
                                         const bool ldiag,
                                         const float *__restrict__ G,
                                         const float *__restrict__ CeInv,
                                         const float *__restrict__ mts,
                                         const float *__restrict__ d,
                                         float *__restrict__ phi);
/*! Karhunen Loeve transform */
int parmt_mtSearch_kl32f(const int nlags, const int nmt,
                         const int npts, const int ldz,
                         const int rank,
                         const float sigma,
                         const float *__restrict__ Gxx,
                         const float *__restrict__ Gyy, 
                         const float *__restrict__ Gzz, 
                         const float *__restrict__ Gxy, 
                         const float *__restrict__ Gxz, 
                         const float *__restrict__ Gyz, 
                         const float *__restrict__ sqrtEig,
                         const float *__restrict__ Z,
                         const float *__restrict__ mts, 
                         const float *__restrict__ d,
                         float *__restrict__ phi);

#ifdef __cplusplus
}
#endif
#endif
