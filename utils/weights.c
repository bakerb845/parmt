#include <stdio.h>
#include <stdlib.h>
#include "parmt_utils.h"
#include "sacio.h"
#include "iscl/log/log.h"

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
    const char *fcnm = "parmt_utils_setWeight\0";
    sacio_setFloatHeader(SAC_WEIGHT_HDR, -12345.0, &obs->header);
    if (weight < 0.0)
    {
        log_errorF("%s: Error the data weight %f cannot be negative\n",
                   fcnm, weight);
        return -1;
    }
    sacio_setFloatHeader(SAC_WEIGHT_HDR, weight, &obs->header);
    return 0;
}
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
    const char *fcnm = "parmt_utils_setPolarityWeight\0";
    sacio_setFloatHeader(SAC_POLWGT_HDR, -12345.0, &obs->header);
    if (weight < 0.0)
    {
        log_errorF("%s: Error the polarity weight %f cannot be negative\n",
                   fcnm, weight);
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
    const char *fcnm = "parmt_utils_getWeight\0";
    double weight;
    int ierr;
    if (defaultWeight <= 0.0)
    {
        log_warnF("%s: It's strange that your default weight isn't positive\n",
                  fcnm);
    }
    *ldefault = true;
    weight = defaultWeight;
    ierr = sacio_getFloatHeader(SAC_WEIGHT_HDR, obs.header, &weight);
    if (ierr == 0)
    {
        *ldefault = false;
        if (weight < 0.0)
        {
            log_errorF("%s: Error weight %f is invalid - defaulting to %f\n",
                       fcnm, weight, defaultWeight);
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
    const char *fcnm = "parmt_utils_getWeight\0";
    double weight;
    int ierr;
    if (defaultWeight <= 0.0)
    {
        log_warnF("%s: It's strange that your default weight isn't positive\n",
                  fcnm);
    }
    *ldefault = true;
    weight = defaultWeight;
    ierr = sacio_getFloatHeader(SAC_POLWGT_HDR, obs.header, &weight);
    if (ierr == 0)
    {
        *ldefault = false;
        if (weight < 0.0)
        {
            log_errorF("%s: Error weight %f is invalid - defaulting to %f\n",
                       fcnm, weight, defaultWeight);
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
