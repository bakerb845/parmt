#ifndef PARMT_INVERT_H__
#define PARMT_INVERT_H__ 1
#include "parmt_utils.h"

#ifdef __cplusplus
extern "C"
{
#endif

int parmt_invertMT64f(const int nobs,
                      const bool lcDev,
                      const int *dataPtr,
                      const double *__restrict__ wts,
                      const double *__restrict__ gxxAll,
                      const double *__restrict__ gyyAll,
                      const double *__restrict__ gzzAll,
                      const double *__restrict__ gxyAll,
                      const double *__restrict__ gxzAll,
                      const double *__restrict__ gyzAll,
                      const double *__restrict__ obsAll,
                      double *phi, double *__restrict__ m);

#ifdef __cplusplus
}
#endif
#endif
