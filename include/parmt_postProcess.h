#ifndef PARMT_POSTPROCESS_H__
#define PARMT_POSTPROCESS_H__ 1
#include <stdint.h>
#include "parmt_config.h"

#ifdef __cplusplus
extern "C"
{
#endif

// for computing cell spacings
int postprocess_computeBetaCellSpacing(const int nb, 
                                       const double *__restrict__ betas,
                                       double *__restrict__ db);
int postprocess_computeGammaCellSpacing(const int ng,
                                        const double *__restrict__ gammas,
                                        double *__restrict__ dg);
int postprocess_computeThetaCellSpacing(const int nt,
                                        const double *__restrict__ thetas,
                                        double *__restrict__ dt);
int postprocess_computeKappaCellSpacing(const int nk,
                                        const double *__restrict__ kappas,
                                        double *__restrict__ dk);
int postprocess_computeSigmaCellSpacing(const int ns,
                                        const double *__restrict__ sigmas,
                                        double *__restrict__ ds);
int postprocess_computeM0CellSpacing(const int nm,
                                     const double *__restrict__ M0s,
                                     double *__restrict__ dm);

int marginal_computeLuneMPDF(const int nloc,
                             const int nm,
                             const int nb, const double *__restrict__ betas,
                             const int ng, const double *__restrict__ gammas,
                             const int nk, const double *__restrict__ kappas,
                             const int ns, const double *__restrict__ sigmas,
                             const int nh, const double *__restrict__ h,
                             const double *__restrict__ phi,
                             double *__restrict__ luneMPDF);
int marginal_computeLuneUVMPDF(const int nloc,
                               const int nm, 
                               const int nu, const double *__restrict__ us, 
                               const int nv, const double *__restrict__ vs, 
                               const int nk, const double *__restrict__ kappas,
                               const int ns, const double *__restrict__ sigmas,
                               const int nh, const double *__restrict__ hs, 
                               const double *__restrict__ phi,
                               double *__restrict__ luneMPDF);
int marginal_computeDepthMPDF(const int nloc,
                              const int nm, const double *__restrict__ M0s,
                              const int nb, const double *__restrict__ betas,
                              const int ng, const double *__restrict__ gammas,
                              const int nk, const double *__restrict__ kappas,
                              const int ns, const double *__restrict__ sigmas,
                              const int nt, const double *__restrict__ thetas,
                              const double *__restrict__ phi,
                              double *__restrict__ depMagMPDF,
                              double *__restrict__ depMPDF);
int marginal_getOptimum(const int nloc, const int nm, const int nb, 
                        const int ng, const int nk, const int ns, const int nt, 
                        const double *__restrict__ phi,
                        int *jloc, int *jm, int *jb, int *jg,
                        int *jk, int *js, int *jt);
int marginal_computeMarginalBeachball(
    const int nloc,
    const int nm, 
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas,
    const double *__restrict__ phi 
    );
int postprocess_tt2tnp(const double beta, const double gamma,
                       const double kappa, const double sigma,
                       const double theta,
                       double *__restrict__ pAxis,
                       double *__restrict__ nAxis,
                       double *__restrict__ tAxis);
int postprocess_tnp2beachballPolarity(const int nxp,
                                      const double xc, const double yc, 
                                      const double rad,
                                      const double *__restrict__ pAxis,
                                      const double *__restrict__ nAxis,
                                      const double *__restrict__ tAxis,
                                      double *__restrict__ xw1,
                                      double *__restrict__ yw1,
                                      int8_t *__restrict__ pn1);
double marginal_computeNormalization(
    const int nloc, const double *__restrict__ deps,
    const int nm, const double *__restrict__ M0s,
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas,
    const double *__restrict__ phi,
    int *ierr);

int marginal_write1DToGnuplot(const char *fname,
                              const int nx, const double *__restrict__ xlocs,
                              const double *__restrict__ vals);
int marginal_write2DToGnuplot(const char *fname,
                              const int nx, const double *__restrict__ xlocs,
                              const int ny, const double *__restrict__ ylocs,
                              const double *__restrict__ vals);

#ifdef __cplusplus
}
#endif

#endif
