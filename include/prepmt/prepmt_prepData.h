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
int prepmt_prepData_verifyTeleseismicDistance(const bool lisP,
                                              const double *dminIn,
                                              const double *dmaxIn,
                                              const struct sacData_struct data);
int prepmt_prepData_verifyRegionalDistance(const double *dminIn,
                                           const double *dmaxIn,
                                           const struct sacData_struct data);
/* Archive utility */
int prepmt_prepData_archiveWaveforms(const char *archiveFile,
                                     const int nobs, 
                                     const struct sacData_struct *data);
struct sacData_struct *prepmt_prepData_readArchivedWaveforms(
    const char *archiveFile, int *nobs, int *ierr);
/* Reads the method to set the pick information */
int prepmt_prepData_readPickModel(const char *iniFile,
                                  const char *prefix,
                                  bool *lsetNewPicks,
                                  bool *lusePickFile,
                                  char pickFile[PATH_MAX],
                                  char ttimesTableDir[PATH_MAX],
                                  char ttimesModel[128]);
/* Reads intermediate file options */
int prepmt_prepData_intermediateFileOptions(const char *iniFile,
                                            const char *section,
                                            bool *lwrtIntFiles,
                                            char wfDir[PATH_MAX],
                                            char wfSuffix[128],
                                            char archiveFile[PATH_MAX]);
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
