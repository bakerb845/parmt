#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <iniparser.h>
#include <mpi.h>
#include "prepmt/prepmt_event.h"
#include "iscl/log/log.h"
#include "iscl/os/os.h"
#include "iscl/time/time.h"

#define M11 1.0        /*! Default to explosion - mxx or mtt */
#define M22 1.0        /*! Default to explosion - myy or mpp */
#define M33 1.0        /*! Default to explosion - mzz or mrr */
#define M12 0.0        /*! Default to explosion - */
#define M13 0.0        /*! Default to explosion - */
#define M23 0.0        /*! Default to explosion - */
#define EVLA 41.30800  /*! Default to North Korea test site latitude */
#define EVLO 129.0760  /*! Default to North Korea test site longitude */
#define EVDP 1.0       /*! Default to shallow depth */
#define EVTIME 0.0     /*! Default origin time is unknown */

/*!
 * @brief Broadcasts the event information from root to all processes on 
 *        communicator.
 *
 * @param[in,out] event   On input the root has the event information.
 *                        On output all processes have the event information.
 *
 * @param[in] root        Process rank that is broadcasting on communicator.
 * @param[in] comm        MPI communicator.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_event_broadcast(struct prepmtEventParms_struct *event,
                           const int root, const MPI_Comm comm)
{
    int ibasis, myid;
    MPI_Comm_rank(comm, &myid);
    MPI_Bcast(&event->m11,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->m22,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->m33,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->m12,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->m13,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->m23,       1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->latitude,  1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->longitude, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->depth,     1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&event->time,      1, MPI_DOUBLE, root, comm);
    if (myid == root){ibasis = (int) event->basis;}
    MPI_Bcast(&ibasis, 1, MPI_INT, root, comm);
    if (myid != root){event->basis = (enum compearthCoordSystem_enum) ibasis;}
    return 0;
}
/*!
 * @brief Sets the default event information for North Korea.
 *
 * @param[out] event   Default event is an explosion in North Korea.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_event_setDefaults(struct prepmtEventParms_struct *event)
{
    memset(event, 0, sizeof(struct prepmtEventParms_struct));
    event->m11 = M11;
    event->m22 = M22;
    event->m33 = M33;
    event->m12 = M12;
    event->m13 = M13;
    event->m23 = M23;
    event->basis = NED;
    event->latitude = EVLA;
    event->longitude = EVLO;
    event->depth = EVDP;
    event->time = EVTIME; 
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the event information from the ini file.
 *
 * @param[in] fname    Name of ini file.
 *
 * @param[out] event   Event information (lat, lon, depth, origin time, and
 *                     moment tensor).
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_event_initializeFromIniFile(
    const char *fname,
    struct prepmtEventParms_struct *event)
{
    const char *fcnm = "prepmt_event_initializeFromIniFile\0";
    const char *s;
    double second;
    int dom, month, nzhour, nzmin, nzsec, nzmusec, nzyear;
    bool lread;
    dictionary *ini;
    //------------------------------------------------------------------------//
    //
    // Set the defaults
    prepmt_event_setDefaults(event);
    // Verify the file exists
    if (!os_path_isfile(fname))
    {
        log_errorF("%s: Error ini file %s doesn't exist\n", fcnm, fname);
        return -1;
    }
    ini = iniparser_load(fname);
    // Read the moment tensor terms
    lread = false;
    s = iniparser_getstring(ini, "eventInfo:basis\0", NULL);
    if (s != NULL)
    {
        if (strcasecmp(s, "USE\0") == 0)
        {
            lread = true;
            event->basis = USE;
        }
        if (strcasecmp(s, "NED\0") == 0)
        {
            lread = true;
            event->basis = NED;
        }
    }
    if (lread && event->basis == USE)
    {
        event->m11 = iniparser_getdouble(ini, "eventInfo:mrr\0", M11);
        event->m22 = iniparser_getdouble(ini, "eventInfo:mtt\0", M22);
        event->m33 = iniparser_getdouble(ini, "eventInfo:mpp\0", M33);
        event->m12 = iniparser_getdouble(ini, "eventInfo:mrt\0", M12);
        event->m13 = iniparser_getdouble(ini, "eventInfo:mrp\0", M13);
        event->m23 = iniparser_getdouble(ini, "eventInfo:mtp\0", M23);
    }
    if (lread && event->basis == NED)
    {
        event->m11 = iniparser_getdouble(ini, "eventInfo:mxx\0", M11);
        event->m22 = iniparser_getdouble(ini, "eventInfo:myy\0", M22);
        event->m33 = iniparser_getdouble(ini, "eventInfo:mzz\0", M33);
        event->m12 = iniparser_getdouble(ini, "eventInfo:mxy\0", M12);
        event->m13 = iniparser_getdouble(ini, "eventInfo:mxz\0", M13);
        event->m23 = iniparser_getdouble(ini, "eventInfo:myz\0", M23);
    }
    event->latitude  = iniparser_getdouble(ini, "eventInfo:latitude\0",  EVLA);
    event->longitude = iniparser_getdouble(ini, "eventInfo:longitude\0", EVLO);
    event->depth     = iniparser_getdouble(ini, "eventInfo:depth\0",     EVDP);
    s = iniparser_getstring(ini, "eventInfo:time\0", NULL);
    if (s != NULL)
    {
        sscanf(s, "%d-%d-%d:%d:%d:%lf",
               &nzyear, &month, &dom, &nzhour, &nzmin, &second);
        nzsec = (int) second;
        nzmusec = (int) ((second - (double) nzsec)*1.e6);
        event->time = time_calendar2epoch2(nzyear, month, dom, nzhour,
                                           nzmin, nzsec, nzmusec);
    }
    iniparser_freedict(ini);
    return 0;
}
