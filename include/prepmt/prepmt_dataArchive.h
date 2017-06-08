#ifndef PREPMT_DATAARCHIVE_H__
#define PREPMT_DATAARCHIVE_H__
#include <hdf5.h>
#include "sacio.h"

#ifdef __cplusplus
extern "C"
{
#endif

int prepmt_dataArchive_closeArchive(const hid_t h5fl);
int prepmt_dataArchive_createArchive(
    const char *fname, //const char *dirnm, const char *projnm,
    const int ndeps,
    const double evla,
    const double evlo,
    const double *evdps);
hid_t prepmt_dataArchive_openArchive(
    const char *fname, int *ierr); //const char *dirnm, const char *projnm, int *ierr);
int prepmt_dataArchive_addObservation(
    const hid_t h5fl, const struct sacData_struct obs);
int prepmt_dataArchive_addGreensFunctions(const hid_t h5fl,
                                          const struct sacData_struct sac,
                                          const struct sacData_struct sacGxx,
                                          const struct sacData_struct sacGyy,
                                          const struct sacData_struct sacGzz,
                                          const struct sacData_struct sacGxy,
                                          const struct sacData_struct sacGxz,
                                          const struct sacData_struct sacGyz);


#ifdef __cplusplus
}
#endif

#endif
