#ifndef PARMT_UTILS_H__
#define PARMT_UTILS_H__ 1
#include "parmt_config.h"
#include <sacio.h>
#include <hdf5.h>
#include <limits.h>

#define SAC_WEIGHT_HDR SAC_FLOAT_USER0
#define SAC_MAXLAG_HDR SAC_FLOAT_USER1
#define SAC_POLWGT_HDR SAC_FLOAT_USER2

#define PARMT_ID 1
#define POLARMT_ID 2
#define POSTMT_ID 3

struct parmtNoiseBasis_struct
{
    double *sqrtEvals; /*!< sqrt(Eigenvalues) in descending order [nvals] */
    double *evects;    /*!< Column major corresponding eigenvectors
                            [lde*nvals] */
    int lde;           /*!< Leading dimension of eigenvectors */
    int nvals;         /*!< Number of eigenvalues */
    int npts;          /*!< Number of points in noise signal */
};

struct parmtGeneralParms_struct
{
    char projnm[256];          /*!< Project name */
    char programName[256];     /*!< Program name */
    char dataFile[PATH_MAX];   /*!< Input data/greens functions file */
    char resultsFileSuffix[PATH_MAX]; /*!< Results file suffix */
    char polarityFileSuffix[PATH_MAX]; /*!< Polarity file suffix */
    char resultsDir[PATH_MAX];        /*!< Results file directory */
    char parmtArchive[PATH_MAX];    /*!< Parmt archive name */
    char polarmtArchive[PATH_MAX];  /*!< Polarity archive name */
    char postmtFile[PATH_MAX];      /*!< Input into postmt */
    double defaultMaxLagTime;  /*!< Maximum max lag time (seconds) that
                                    can be applied to any observation */
    int blockSize;             /*!< Block size in matrix-matrix GM=U
                                    multliplications */
    int objFnType;             /*!< If 1 then this is L1 norm.
                                    If 2 then this is a cross-correlation */
    int programID; 
    bool lwantLags;            /*!< If true then the program is to optimize
                                    at different synthetic/data lags */
    char pad[3];
};

struct parmtPolarityParms_struct
{
    char ttimesTablesDir[PATH_MAX]; /*!< ttimes travel-times table directory*/
    char ttimesModel[64];           /*!< ttimes model (default is ak135) */
    bool lcomputePolarity;          /*!< if true compute polarity info */
    char pad[3];
};

struct parmtData_struct
{
    struct sacData_struct *data;   /*!< Observed waveform [nobs] */
    struct sacData_struct *sacGxx; /*!< Green's functions scaling mxx
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *sacGyy; /*!< Green's functions scaling myy
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *sacGzz; /*!< Green's functions scaling mzz
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *sacGxy; /*!< Green's functions scaling mxy
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *sacGxz; /*!< Green's functions scaling mxz
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *sacGyz; /*!< Green's functions scaling myz
                                        moment tensor term [nobs x nlocs] */
    struct sacData_struct *est;    /*!< Corresponding best fitting estimate */
    int nlocs;     /*!< Number of locations in grid-search */
    int nobs;      /*!< Number of waveforms in grid-search */
};

struct parmtMtSearchParms_struct
{
    double betaLower;  /*!< Lower colatitude \f$ \beta \in [0, \pi] \f$ */
    double betaUpper;  /*!< Upper colatitude \f$ \beta \in [0, \pi] \f$ */
    double gammaLower; /*!< Lower longitude
                            \f$ \gamma \in [-\frac{\pi}{6},\frac{\pi}{6}] \f$ */
    double gammaUpper; /*!< Upper longitude
                            \f$ \gamma \in [-\frac{\pi}{6},\frac{\pi}{6}] \f$ */
    double kappaLower; /*!< Lower strike angle \f$ \kappa \in [0, 2\pi] \f$ */
    double kappaUpper; /*!< Upper strike nagle \f$ \kappa \in [0, 2\pi] \f$ */
    double sigmaLower; /*!< Lower rake angle
                            \f$ \kappa \in [-\frac{\pi}{2},\frac{\pi}{2}] \f$ */
    double sigmaUpper; /*!< Uppwer rake angle
                            \f$ \kappa \in [-\frac{\pi}{2},\frac{\pi}{2}] \f$ */
    double thetaLower; /*!< Lower dip angle
                            \f$ \theta \in [0, \frac{\pi}{2} ] \f$ */
    double thetaUpper; /*!< Upper dip angle
                            \f$ \theta \in [0, \frac{\pi}{2} ] \f$ */
    double m0Lower;    /*!< Lower scalar moment (N-m) */
    double m0Upper;    /*!< Upper scalar moment (N-m) */
    int nb; /*!< Number of colatitudes */
    int ng; /*!< Number of longitudes */ 
    int nk; /*!< Number of strike angles */
    int ns; /*!< Number of rake angles */
    int nt; /*!< Number of dip angles */
    int nm; /*!< Number of scalar moments */
    bool luseLog; /*!< If true then use a log scale when discretizing MT. */
    char pad[3];
};

struct parmtDatas_struct
{
    struct parmtData_struct *data;
    int nobs;
};

#ifdef __cplusplus
extern "C"
{
#endif

int parmt_utils_readGeneralParms(const char *iniFile,
                                 struct parmtGeneralParms_struct *parms);
int parmt_utils_readPolarityParms(const char *iniFile,
                                  struct parmtPolarityParms_struct *parms);
int parmt_utils_readMtSearch(const char *iniFile,
                             struct parmtMtSearchParms_struct *parms);
int parmt_utils_getEffectiveRank(const double pct,
                                 const struct parmtNoiseBasis_struct basis);
int parmt_utils_makeKLNoise(const int rank,
                            const int npts,
                            const struct parmtNoiseBasis_struct basis,
                            double *__restrict__ xsig);
int parmt_utils_freeNoiseBasis(struct parmtNoiseBasis_struct *basis);
int parmt_utils_getNoiseBasis64f(const int npts,
                                 const double *__restrict__ data,
                                 struct parmtNoiseBasis_struct *basis);
int parmt_utils_getNoiseBasisFromCtC64f(const int npts, const int ldctc,
                                        const double *__restrict__ CtC,
                                        struct parmtNoiseBasis_struct *basis);

int parmt_utils_ff2mtGreensMatrix64f(const int npgrns, const int ldg,
                                     const int icomp,
                                     const double azin, const double bazin,
                                     const double cmpaz, const double cmpinc,
                                     const double *__restrict__ ZDS,
                                     const double *__restrict__ ZSS,
                                     const double *__restrict__ ZDD,
                                     const double *__restrict__ ZEX,
                                     const double *__restrict__ RDS,
                                     const double *__restrict__ RSS,
                                     const double *__restrict__ RDD,
                                     const double *__restrict__ REX,
                                     const double *__restrict__ TDS,
                                     const double *__restrict__ TSS,
                                     double *__restrict__ G);
int parmt_utils_ff2mtGreens64f(const int npgrns, const int icomp,
                               const double azin, const double bazin,
                               const double cmpaz, const double cmpinc,
                               const double *__restrict__ ZDS,
                               const double *__restrict__ ZSS,
                               const double *__restrict__ ZDD,
                               const double *__restrict__ ZEX,
                               const double *__restrict__ RDS,
                               const double *__restrict__ RSS,
                               const double *__restrict__ RDD,
                               const double *__restrict__ REX,
                               const double *__restrict__ TDS,
                               const double *__restrict__ TSS,
                               double *__restrict__ Gxx,
                               double *__restrict__ Gyy,
                               double *__restrict__ Gzz,
                               double *__restrict__ Gxy,
                               double *__restrict__ Gxz,
                               double *__restrict__ Gyz);
//-------------------------------synthetics-----------------------------------//
double *parmt_utils_sacGrns2GrnsMat(const struct sacData_struct sacGxx,
                                    const struct sacData_struct sacGyy,
                                    const struct sacData_struct sacGzz,
                                    const struct sacData_struct sacGxy,
                                    const struct sacData_struct sacGxz,
                                    const struct sacData_struct sacGyz,
                                    int *mrows, int *ierr);
int parmt_utils_sacGrnsToEst(const struct sacData_struct data,
                             const struct sacData_struct sacGxx,
                             const struct sacData_struct sacGyy,
                             const struct sacData_struct sacGzz,
                             const struct sacData_struct sacGxy,
                             const struct sacData_struct sacGxz,
                             const struct sacData_struct sacGyz,
                             const double *__restrict__ mt, 
                             struct sacData_struct *est);
//-------------------------------file io--------------------------------------//
int parmt_io_createObjfnArchive64f(
    //const char *archiveDir, const char *projnm, const char *suffix,
    const char *programName, const char *flname,
    const int nobs,
    const int nlocs,
    const double *__restrict__ deps,
    const int nm, const double *__restrict__ M0s,
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas);
int parmt_io_writeVarianceForWaveform64f(
    //const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const char *flname,
    const int waveformNumber,
    const int npts, const double *__restrict__ var);
int parmt_io_writeObjectiveFunction64f(
    //const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const char *flname,
    const int nmt, const double *__restrict__ phi);
int parmt_io_readObjfnArchive64f(
    //const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const char *flname,
    char programName[256],
    int *nlocs, double **deps,
    int *nm, double **M0s,
    int *nb, double **betas,
    int *ng, double **gammas,
    int *nk, double **kappas,
    int *ns, double **sigmas,
    int *nt, double **thetas,
    int *nmt, double **phi);
int parmt_utils_setDataOnG(const int iobs, const int iloc, const int npmax,
                           const struct parmtData_struct data, 
                           double *__restrict__ G);
int utils_dataArchive_readAllWaveforms(const char *dataFile,
                                       struct parmtData_struct *data);
int utils_dataArchive_getLocationID(const hid_t h5fl,
                                    const double evla,
                                    const double evlo,
                                    const double evdp);
int utils_dataArchive_getObservationID(const hid_t h5fl,
                                       const struct sacData_struct obs);
int utils_dataArchive_addObservation(const hid_t h5fl,
                                     const struct sacData_struct obs);
int utils_dataArchive_setFileName(const char *dirnm, const char *projnm,
                                 char fname[PATH_MAX]);
int utils_dataArchive_initialize(const char *fname,
                                 const int nlocs,
                                 const double *__restrict__ evlas,
                                 const double *__restrict__ evlos,
                                 const double *__restrict__ evdps);
int utils_dataArchive_getNumberOfObservations(const hid_t h5fl);
int utils_dataArchive_getObservationID(const hid_t h5fl,
                                       const struct sacData_struct obs);
int utils_dataArchive_addGreensFunctions(const hid_t h5fl,
                                         const int waveformID,
                                         const int locationID,
                                         const struct sacData_struct sacGxx,
                                         const struct sacData_struct sacGyy,
                                         const struct sacData_struct sacGzz,
                                         const struct sacData_struct sacGxy,
                                         const struct sacData_struct sacGxz,
                                         const struct sacData_struct sacGyz);
int utils_dataArchive_loadGreensFunctions(const hid_t h5fl,
                                          const int waveformID,
                                          const int locationID,
                                          struct sacData_struct *sacGxx,
                                          struct sacData_struct *sacGyy,
                                          struct sacData_struct *sacGzz,
                                          struct sacData_struct *sacGxy,
                                          struct sacData_struct *sacGxz,
                                          struct sacData_struct *sacGyz);
//--------------------------------special header variables--------------------//
/* Gets and sets observation weight */
double parmt_utils_getWeight(const struct sacData_struct obs,
                             const double defaultWeight,
                             bool *ldefault);
int parmt_utils_setPolarityWeight(const double weight,
                                  struct sacData_struct *obs);
int parmt_utils_setWeight(const double weight,
                          struct sacData_struct *obs);
double parmt_utils_getPolarityWeight(const struct sacData_struct obs,
                                     const double defaultWeight,
                                     bool *ldefault);
/* Gets and sets observation max lag time */
int parmt_utils_setLagTime(const double maxLagTime,
                           struct sacData_struct *obs);
double parmt_utils_getLagTime(const struct sacData_struct obs,
                              const double defaultMaxLag,
                              bool *ldefault);

#ifdef __cplusplus
}
#endif
#endif
