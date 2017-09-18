#ifndef PREPMT_GREENS_H__
#define PREPMT_GREENS_H__ 1
#include "prepmt/prepmt_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

int prepmt_greens_writeArchive(const char *archiveName, 
                               const int nwaves, const int ndepths,
                               const double evla, const double evlo,
                               const double *__restrict__ depths,
                               const struct sacData_struct *sac,
                               const struct sacData_struct *sacGrns);
int prepmt_greens_ffGreensToGreens(const int nobs,
                                   const struct sacData_struct *obs,
                                   const int ndepth,
                                   const int ntstar,
                                   const struct sacData_struct *ffGrns,
                                   struct sacData_struct *grns);

int prepmt_greens_getHudson96GreensFunctionIndex(
    const enum prepmtGreens_enum GMT_TERM,
    const int nobs, const int ntstar, const int ndepth, 
    const int iobs, const int itstar, const int idepth);

int prepmt_greens_getHudson96GreensFunctionsIndices(
    const int nobs, const int ntstar, const int ndepth,
    const int iobs, const int itstar, const int idepth,
    int indices[6]);

int prepmt_greens_processHudson96Greens(
    const int nobs, const int ntstar, const int ndepth,
    const struct prepmtCommands_struct cmds,
    struct sacData_struct *grns);
int prepmt_greens_cutHudson96FromData(const int nobs,
                                      const struct sacData_struct *data,
                                      const int ndepth, const int ntstar,
                                      struct sacData_struct *grns);
int prepmt_greens_repickGreensWithSTALTA(
    const double sta, const double lta, const double threshPct,
    struct sacData_struct *grns);
int prepmt_greens_xcAlignGreensToData(const struct sacData_struct data,
                                      const bool luseEnvelope, const bool lnorm,
                                      const double maxTimeLag,
                                      struct sacData_struct *grns);
int prepmt_greens_xcAlignGreensToData_work(
    const int npts, double *__restrict__ data,
    const int npgrns,
    const bool luseEnvelope, const bool norm,
    const int maxShift,
    double *__restrict__ Gxx,
    double *__restrict__ Gyy,
    double *__restrict__ Gzz,
    double *__restrict__ Gxy,
    double *__restrict__ Gxz,
    double *__restrict__ Gyz,
    const int lxc, double *__restrict__ xc);

double prepmt_greens_scoreXCAlignment(const int npts,
                                      const bool luseEnvelope,
                                      const bool lnorm,
                                      const double *__restrict__ data,
                                      const double *__restrict__ Gxx,
                                      const double *__restrict__ Gyy,
                                      const double *__restrict__ Gzz,
                                      const double *__restrict__ Gxy,
                                      const double *__restrict__ Gxz,
                                      const double *__restrict__ Gyz,
                                      int *ierr);

#ifdef __cplusplus
}
#endif

#endif
