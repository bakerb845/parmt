#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "parmt_utils.h"
#include "sacio.h"
#include "iscl/log/log.h"

/*!
 * @brief Sets the max lag time for the given network, station, channel,
 *        location in the H5 data archive.
 *
 * @param[in] h5fl        HDF5 data file handle.  The file must be opened with
 *                        H5F_ACC_RDWR.
 * @param[in] network     Network name.
 * @param[in] station     Station name.
 * @param[in] channel     Channel name.
 * @param[in] location    Location name.
 * @param[in] maxLagTime  Max lag time (seconds).  If negative then this
 *                        will be ignored.
 * 
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int parmt_utils_setLagTimeInH5(const hid_t h5fl,
                               const char *network, const char *station,
                               const char *channel, const char *location,
                               const double maxLagTime)
{
    const char *fcnm = "parmt_utils_setLagTimeInH5\0";
    const char *sacFileName = "Observation";
    char varname[256], knetwk[8], kstnm[8], kcmpnm[8], khole[8];
    bool lupd;
    struct sacData_struct sac;
    hid_t obsGroup;
    int ierr, ifound, iobs, nobs;

    // Get number of objects in group
    nobs = utils_dataArchive_getNumberOfObservations(h5fl);
    if (nobs <= 0)
    {
        log_errorF("%s: No observations\n", fcnm);
        return -1;
    }
    lupd = false;
    memset(&sac, 0, sizeof(struct sacData_struct));
    // Loop on observations
    for (iobs=0; iobs<nobs; iobs++)
    {
        // Open the group
        sprintf(varname, "/Observations/Observation_%d", iobs);
        obsGroup = H5Gopen(h5fl, "Observations", H5P_DEFAULT);
        // Load the SAC file
        ierr = sacioh5_readTimeSeries2("Observation\0", obsGroup, &sac);
        if (ierr != 0)
        {
            log_errorF("%s: Error loading time series for obs %d\n",
                       fcnm, iobs+1);
        }
        // Check if this station matches
        ifound = 0;
        sacio_getCharacterHeader(SAC_CHAR_KNETWK, sac.header, knetwk);
        if (strcasecmp(network, knetwk) == 0){ifound = ifound + 1;}
        sacio_getCharacterHeader(SAC_CHAR_KSTNM,  sac.header, kstnm);
        if (strcasecmp(station, kstnm) == 0){ifound = ifound + 1;}
        sacio_getCharacterHeader(SAC_CHAR_KCMPNM, sac.header, kcmpnm);
        if (strcasecmp(channel, kcmpnm) == 0){ifound = ifound + 1;}
        sacio_getCharacterHeader(SAC_CHAR_KHOLE, sac.header, khole);  
        if (strcasecmp(location, khole) == 0){ifound = ifound + 1;}
        if (ifound == 4)
        {
            ierr = parmt_utils_setLagTime(maxLagTime, &sac);
            if (ierr != 0)
            {
                printf("%s: Failed to set lag time\n", fcnm);
            }
            else
            {
                sacioh5_writeTimeSeries2("Observation\0", obsGroup, sac);
            }
            lupd = true;
        }
        sacio_free(&sac);
        H5Gclose(obsGroup);
        if (ierr != 0){return -1;}
        if (lupd){break;}
    }
    if (!lupd)
    {
        log_errorF("%s: Couldn't find %s.%s.%s.%s in archive\n",
                   network, station, channel, location);
        return -1;
    }
    return 0;
}
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
