#ifndef PREPMT_EVENT_H__
#define PREPMT_EVENT_H__ 1
#include "prepmt/prepmt_struct.h"
#include "compearth.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C"
{
#endif
int prepmt_event_broadcast(struct prepmtEventParms_struct *event,
                           const int root, const MPI_Comm comm);
int prepmt_event_setDefaults(struct prepmtEventParms_struct *event);
int prepmt_event_initializeFromIniFile(
    const char *fname,
    struct prepmtEventParms_struct *event);

#ifdef __cplusplus
}
#endif
#endif
