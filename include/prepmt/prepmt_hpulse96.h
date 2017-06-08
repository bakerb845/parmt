#ifndef PREPMT_HPULSE96_H__
#define PREPMT_HPULSE96_H__
#include <cps.h>

#ifdef __cplusplus
extern "C"
{
#endif

int prepmt_hpulse96_readHpulse96Parameters(const char *iniFile,
                                           const char *section,
                                           struct hpulse96_parms_struct *parms);

#ifdef __cplusplus
}
#endif

#endif
