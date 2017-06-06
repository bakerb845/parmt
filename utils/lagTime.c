#include <stdio.h>
#include <stdlib.h>
#include "parmt_utils.h"
#include "sacio.h"
#include "iscl/log/log.h"

/*!
 * @brief Sets the lag-time for this observation.
 *
 * @param[in] maxLagTime   Maximum lag time (seconds).  If negative then
 *                         this variable will be ignored.
 *
 * @param[out] obs      Holds the data lag-time in the header. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_setLagTime(const double maxLagTime,
                           struct sacData_struct *obs)
{
    const char *fcnm = "parmt_utils_setLagTime\0";
    sacio_setFloatHeader(SAC_MAXLAG_HDR, maxLagTime, &obs->header);
    return 0;
}
//============================================================================//
/*!
 * @brief Gets the lag-time for this observation.
 *
 * @param[in] obs            SAC observation.
 * @param[in] defaultMaxLag  If the max lag-time has not been set then this
 *                           will be the default. 
 *
 * @param[out] ldefault      If true then the lag-time could not be read
 *                           from the header and the default value is set.
 *
 * @result The data max lag-time (seconds) corresponding to the given
 *         data observation.
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
double parmt_utils_getLagTime(const struct sacData_struct obs,
                              const double defaultMaxLag,
                              bool *ldefault)
{
    const char *fcnm = "parmt_utils_getLagTime\0";
    double lagtime;
    int ierr;
    *ldefault = true;
    lagtime = defaultMaxLag;
    ierr = sacio_getFloatHeader(SAC_MAXLAG_HDR, obs.header, &lagtime);
    if (ierr != 0 || lagtime < 0.0)
    {
        *ldefault = true;
        lagtime = defaultMaxLag;
    }
    return lagtime;
}
