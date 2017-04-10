#ifndef PARMT_POLARITY_H__
#define PARMT_POLARITY_H__ 1
#include "sacio.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! High level driver which extracts requisite data from the SAC header
    for computing the take-off angle and incidence angle
    from ttimes for a 1D layered earth model for computation of 
    row of the greens function polarity matrix */
int parmt_polarity_computeGreensRowFromData(const struct sacData_struct data,
                                            const int wavetype,
                                            const double evdp,
                                            const char *dirnm,
                                            const char *model,
                                            double *__restrict__ G);
/*! Boot straps the take-off angle and incidence angle from ttimes for
    1D layered earth models at global scales for computation of a 
    row of the greens function polarity matrix */
int parmt_polarity_computeGreensRowFromCoordinates(
    const int wavetype, const int icomp,
    const double evla, const double evlo, const double evdp,
    const double stla, const double stlo,
    const double cmpinc, const double cmpaz,
    const char *dirnm, const char *model, 
    double *__restrict__ G);
/*! Computes a row of the Green's functions polarity matrix */
int parmt_polarity_computeGreensMatrixRow(const int wavetype,
                                          const int icomp,
                                          const double azSrc,
                                          const double toaSrc,
                                          const double bazRec,
                                          const double aoiRec,
                                          const double cmpinc,
                                          const double compaz,
                                          double *__restrict__ G);


#ifdef __cplusplus
}
#endif
#endif
