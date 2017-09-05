#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "parmt_utils.h"
#include "sacio.h"

#define DATA_WEIGHT_JOB 1
#define POLARITY_WEIGHT_JOB 2

/*!
 * @brief Utility for setting data weight in H5 file.
 *
 * @param[in] h5fl        HDF5 data file handle.  The file must be opened with
 *                        H5F_ACC_RDWR.
 * @param[in] network     Network name.
 * @param[in] station     Station name.
 * @param[in] channel     Channel name.
 * @param[in] location    Location name.
 * @param[in] weight      Data weight to be set in SAC header.
 *                        This cannot be negative.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int parmt_utils_setWeightInH5(const hid_t h5fl,
                              const char *network, const char *station,
                              const char *channel, const char *location,
                              const double weight)
{
    const int job = DATA_WEIGHT_JOB;
    int ierr;
    ierr = parmt_utils_setAnyWeightInH5(h5fl, job,
                                        network, station, channel, location,
                                        weight);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting data weight\n", __func__);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Utility for setting polarity weight in H5 file.
 *
 * @param[in] h5fl        HDF5 data file handle.  The file must be opened with
 *                        H5F_ACC_RDWR.
 * @param[in] network     Network name.
 * @param[in] station     Station name.
 * @param[in] channel     Channel name.
 * @param[in] location    Location name.
 * @param[in] weight      Polarity weight to be set in SAC header.
 *                        This cannot be negative.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int parmt_utils_setPolarityWeightInH5(const hid_t h5fl,
                                      const char *network, const char *station,
                                      const char *channel, const char *location,
                                      const double weight)
{
    const int job = POLARITY_WEIGHT_JOB;
    int ierr;
    ierr = parmt_utils_setAnyWeightInH5(h5fl, job,
                                        network, station, channel, location,
                                        weight);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error setting data weight\n", __func__);
    }   
    return ierr;
}
//============================================================================//
/*!
 * @brief Utility for setting data or polarity weight in H5 file.  It is
 *        recommended one use setDataWeightInH5 or setPolarityWeightInH5 
 *        instead of this function.
 *
 * @param[in] h5fl        HDF5 data file handle.  The file must be opened with
 *                        H5F_ACC_RDWR.
 * @param[in] job         If job == 1 then this sets the data weight.
 *                        If job == 2 the nthis sets the polarity weight.
 * @param[in] network     Network name.
 * @param[in] station     Station name.
 * @param[in] channel     Channel name.
 * @param[in] location    Location name.
 * @param[in] weight      Data weight to be set in SAC header.
 *                        This cannot be negative.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int parmt_utils_setAnyWeightInH5(const hid_t h5fl, const int job,
                                 const char *network, const char *station,
                                 const char *channel, const char *location,
                                 const double weight)
{
    const char *sacFileName = "Observation";
    char varname[256], knetwk[8], kstnm[8], kcmpnm[8], khole[8];
    bool lupd;
    struct sacData_struct sac;
    hid_t obsGroup;
    int ierr, ifound, iobs, nobs;
    if (weight < 0.0)
    {
        fprintf(stderr, "%s: Error weight cannot be negative\n", __func__);
        return -1;
    }
    if (job < 1 || job > 2)
    {
        fprintf(stderr, "%s: job must be 1 or 2\n", __func__);
        return -1;
    }
    // Get number of objects in group
    nobs = utils_dataArchive_getNumberOfObservations(h5fl);
    if (nobs <= 0)
    {   
        fprintf(stderr, "%s: No observations\n", __func__);
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
            fprintf(stderr, "%s: Error loading time series for obs %d\n",
                    __func__, iobs+1);
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
            if (job == 1)
            {
                ierr = parmt_utils_setWeight(weight, &sac);
            }
            else
            {
                ierr = parmt_utils_setPolarityWeight(weight, &sac);
            }
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Failed to set weight for job %d\n",
                        __func__, job);
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
        fprintf(stderr, "%s: Couldn't find %s.%s.%s.%s in archive\n",
                __func__, network, station, channel, location);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the weight for this observation.
 *
 * @param[in] weight    Data weight (cannot be negative).
 *
 * @param[out] obs      Holds the data weight in the header. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_setWeight(const double weight,
                          struct sacData_struct *obs)
{
    sacio_setFloatHeader(SAC_WEIGHT_HDR, -12345.0, &obs->header);
    if (weight < 0.0)
    {
        fprintf(stderr, "%s: Error the data weight %f cannot be negative\n",
                __func__, weight);
        return -1;
    }
    sacio_setFloatHeader(SAC_WEIGHT_HDR, weight, &obs->header);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the weight for this observation.
 *
 * @param[in] weight    Polarity weight (cannot be negative).
 *
 * @param[out] obs      Holds the polarity weight in the header. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_setPolarityWeight(const double weight,
                                  struct sacData_struct *obs)
{
    sacio_setFloatHeader(SAC_POLWGT_HDR, -12345.0, &obs->header);
    if (weight < 0.0)
    {
        fprintf(stderr, "%s: Error the polarity weight %f cannot be negative\n",
                __func__, weight);
        return -1;
    }
    sacio_setFloatHeader(SAC_POLWGT_HDR, weight, &obs->header);
    return 0;
}
//============================================================================//
/*!
 * @brief Gets the weight for this observation.
 *
 * @param[in] obs            SAC observation.
 * @param[in] defaultWeight  If the weight has not been set then this
 *                           will be the default weight (probably 1).
 *
 * @param[out] ldefault      If true then the weight could not be read
 *                           from the header and the default value is set.
 *
 * @result The data weight corresponding to the given data observation.
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
double parmt_utils_getWeight(const struct sacData_struct obs,
                             const double defaultWeight,
                             bool *ldefault)
{
    double weight;
    int ierr;
    if (defaultWeight <= 0.0)
    {
        fprintf(stdout, "%s: It's strange your default weight isn't positive\n",
                __func__);
    }
    *ldefault = true;
    weight = defaultWeight;
    ierr = sacio_getFloatHeader(SAC_WEIGHT_HDR, obs.header, &weight);
    if (ierr == 0)
    {
        *ldefault = false;
        if (weight < 0.0)
        {
            fprintf(stderr, "%s: Error invalid weight=%f; defaulting to %f\n",
                    __func__, weight, defaultWeight);
            weight = defaultWeight;
            *ldefault = true;
        }
    }
    else
    {
        weight = defaultWeight;
    }
    return weight;
}
//============================================================================//
/*!
 * @brief Gets the polarity weight for this observation.
 *
 * @param[in] obs            SAC observation.
 * @param[in] defaultWeight  If the weight has not been set then this
 *                           will be the default weight (probably 1).
 *
 * @param[out] ldefault      If true then the weight could not be read
 *                           from the header and the default value is set.
 *
 * @result The polarity weight corresponding to the given data observation.
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
double parmt_utils_getPolarityWeight(const struct sacData_struct obs,
                                     const double defaultWeight,
                                     bool *ldefault)
{
    double weight;
    int ierr;
    if (defaultWeight <= 0.0)
    {
        fprintf(stdout, "%s: It's strange your default weight isn't positive\n",
                __func__);
    }
    *ldefault = true;
    weight = defaultWeight;
    ierr = sacio_getFloatHeader(SAC_POLWGT_HDR, obs.header, &weight);
    if (ierr == 0)
    {
        *ldefault = false;
        if (weight < 0.0)
        {
            fprintf(stderr, "%s: Error invalid weight=%f; defaulting to %f\n",
                     __func__, weight, defaultWeight);
            weight = defaultWeight;
            *ldefault = true;
        }
    }
    else
    {
        weight = defaultWeight;
    }
    return weight;
}
