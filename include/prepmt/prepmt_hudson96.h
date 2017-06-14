#ifndef PREPMT_HUDSON96_H__
#define PREPMT_HUDSON96_H__
#include <cps.h>
#include <sacio.h>
#include "prepmt/prepmt_config.h"

#ifdef __cplusplus
extern "C"
{
#endif

int prepmt_hudson96_readHudson96Parameters(const char *iniFile,
                                           const char *section,
                                           struct hudson96_parms_struct *parms);

struct sacData_struct *prepmt_hudson96_computeGreensFF(
    const struct hudson96_parms_struct hudson96Parms,
    const struct hpulse96_parms_struct hpulse96Parms,
    //const bool luseCrust1, const char *crust1Dir,
    //const bool luseSrcModel, const char *srcModelName,
    const struct vmodel_struct telmod,
    const struct vmodel_struct srcmod,
    const struct vmodel_struct *recmod,
    const int ntstar, const double *__restrict__ tstars,
    const int ndepth, const double *__restrict__ depths,
    const int nobs, const struct sacData_struct *obs, int *ierr);

int prepmt_hudson96_observationDepthTstarToIndex(
    const int ndepth, const int ntstar,
    const int iobs, const int idep, const int it);


int hudson96_getModels(const int nobs, const struct sacData_struct *obs,
                       const bool luseCrust1, const char *crust1Dir,
                       const bool luseSrcModel, const char *srcModelName,
                       struct vmodel_struct *telmod,
                       struct vmodel_struct *srcmod, 
                       struct vmodel_struct *recmod);

int prepmt_hudson96_grd2ijk(const int igrd,
                            const int n1, const int n2, const int n3, 
                            int *i, int *j, int *k);

#ifdef __cplusplus
}
#endif

#endif
