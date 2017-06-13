#ifndef PARMT_PREPDATA_H__
#define PARMT_PREPDATA_H__ 1
#include <stdbool.h>
#include "sacio.h"
#include "prepmt/prepmt_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Utility function to process data */
int prepmt_prepData_process(const struct prepmtCommands_struct cmds,
                            const int nobs, struct sacData_struct *data);
/* Utility for reading the default sampling period and cut start/end time */
int prepmt_prepData_getDefaultDTAndWindowFromIniFile(const char *iniFile,
                                                     const char *section,
                                                     double *targetDt,
                                                     double *cutStart,
                                                     double *cutEnd);
/* Utility function for reading the data and pole-zero file names. */
int prepmt_prepData_readDataListFromIniFile(const char *iniFile,
                                            const char *prefix,
                                            int *nfiles,
                                            char ***sacFiles,
                                            char ***sacpzFiles);
/* Utility function for setting the event information */
int prepmt_prepData_setEventInformation(const double evla,
                                        const double evlo,
                                        const double evdp,
                                        const double evtime,
                                        const int nobs,
                                        struct sacData_struct *data);
/* Utility function for verifying event / station distances */
int prepmt_prepData_verifyTeleseismicDistance(const double *dminIn,
                                              const double *dmaxIn,
                                              const struct sacData_struct data);
int prepmt_prepData_verifyRegionalDistance(const double *dminIn,
                                           const double *dmaxIn,
                                           const struct sacData_struct data);
/* Utilities for computing theoretical travel times. */
int prepmt_prepData_setTheoreticalSurfaceWaveArrivalTime(
    const double vel, const bool lr,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    const int nobs, struct sacData_struct *data);
int prepmt_prepData_setPrimaryPorSPickFromTheoreticalTime(
    const char *dirnm, const char *model,
    const bool ldop,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    const int nobs, struct sacData_struct *data);
int prepmt_prepData_computeTheoreticalPorSPickTimes(
    const char *dirnm, const char *model,
    const bool ldop,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes);
int prepmt_prepData_computeTheoreticalPPickTimes(
    const char *dirnm, const char *model,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes);
int prepmt_prepData_computeTheoreticalSPickTimes(
    const char *dirnm, const char *model,
    const int nobs, const struct sacData_struct *data,
    double *__restrict__ ptimes);
#ifdef __cplusplus
}
#endif
#endif
