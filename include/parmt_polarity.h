#ifndef PARMT_POLARITY_H__
#define PARMT_POLARITY_H__ 1
#include "sacio.h"
#include "parmt_mtsearch.h"

struct polarityData_struct
{
    double *Gxx;      /*!< Polarity Green's functions scaling mxx term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           mxx is in NED system. */
    double *Gyy;      /*!< Polarity Green's functions scaling myy term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           myy is in NED system. */
    double *Gzz;      /*!< Polarity Green's functions scaling mzz term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           mzz is in NED system. */
    double *Gxy;      /*!< Polarity Green's functions scaling mxy term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           mxy is in NED system. */
    double *Gxz;      /*!< Polarity Green's functions scaling mxz term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           mxz is in NED system. */
    double *Gyz;      /*!< Polarity Green's functions scaling myz term.
                           [nlocs x nPolarity] with leading dimension nPolarity.
                           myz is in NED system. */
    double *polarity; /*!< Polarity (+/-) for ip'th polarity [nPolarity] */
    double *wts;      /*!< Data weights (e.g. pick quality) for each
                           observation (0,1]) */
    int *obsMap;      /*!< Maps from local polarity observation to original
                           waveform observation [nPolarity] */
    int nPolarity;    /*!< Number of polarities */
    int nobs;         /*!< Number of observations (>= nPolarity) */
    int nlocs;        /*!< Number of locations */
};

#ifdef __cplusplus
extern "C"
{
#endif


int polarity_performLocationSearch64f(const MPI_Comm locComm,
                                      const int blockSize,
                                      struct polarityData_struct polarityData, 
                                      struct localMT_struct mtloc,
                                      double *__restrict__ phi);
/*! Highest level driver which computes Green's functions matrix for
    estimation */
int parmt_polarity_computeTTimesGreens(
    const MPI_Comm globalComm,
    const struct parmtPolarityParms_struct parms,
    const struct parmtData_struct data,
    struct polarityData_struct *polarityData);
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
int parmt_polarity_computeGreensRowFromTtimes(
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
