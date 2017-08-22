#ifndef PREPMT_PICKFILE_H__
#define PREPMT_PICKFILE_H__ 1
#include <sacio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*! Reads NonLinLoc phase arrivals onto SAC data */
int prepmt_pickFile_nonLinLoc2sac(const char *pickFile,
                                  const enum sacHeader_enum pickVar,
                                  const enum sacHeader_enum pickNameVar,
                                  const int nobs, struct sacData_struct *data);
/* SAC data to line in NonLinLoc phase arrival file */
int prepmt_pickFile_sacToNonLinLocLine(const struct sacData_struct data,
                                       const enum sacHeader_enum pickVar,
                                       const enum sacHeader_enum pickNameVar,
                                       const enum sacHeader_enum errMagVar,
                                       const enum sacHeader_enum codaDurVar,
                                       const enum sacHeader_enum ampVar,
                                       const enum sacHeader_enum periodVar,
                                       const enum sacHeader_enum weightVar,
                                       char line[256]);

#ifdef __cplusplus
}
#endif
#endif
