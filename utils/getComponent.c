#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parmt_utils.h"

/*!
 * @brief Classifies the component in the left-handed coordinate system
 *        ZNE or Z12 or Z23.
 *
 * @param[in] kcmpnm    Channel name.  What is important here is the third
 *                      letter.
 * @param[in] cmpinc    Component incidence angle from vertical.  This is
 *                      measured in degrees.  0 degrees is vertical up,
 *                      180 degrees is vertical down, and 90 degrees is
 *                      horizontal.  This is the SAC convention.  To convert
 *                      from the SEED convention use cmpinc + 90.
 * @param[out] icomp    Component where:
 * @param[out] icomp    1 is vertical.
 * @param[out] icomp    2 is north.
 * @param[out] icomp    3 is east.
 *
 * @result 0 indicates success.
 *
 */
int parmt_utils_getComponent(const char *kcmpnm,
                             const double cmpinc,
                             int *icomp)
{
    int ierr;
    ierr = 0;
    *icomp = 0;
    if (fabs(cmpinc - 0.0) < 1.e-4 || fabs(cmpinc - 180.0) < 1.e-4)
    {
        if (kcmpnm[2] == 'Z' || kcmpnm[2] == 'z' || kcmpnm[2] == '1')
        {
            *icomp = 1;
        }
        else
        {
            fprintf(stderr, "%s: cmpinc is 0 but channel=%s isn't vertical",
                    __func__, kcmpnm);
            ierr = 1;
        }
    }
    else if (fabs(cmpinc - 90.0) < 1.e-4)
    {
        *icomp = 2;
        if (kcmpnm[2] == 'N' || kcmpnm[2] == 'n' ||
            kcmpnm[2] == '1' || kcmpnm[2] == '2')
        {
            *icomp = 2;
        }
        else if (kcmpnm[2] == 'N' || kcmpnm[2] == 'n' ||
                 kcmpnm[2] == '2' || kcmpnm[2] == '3')
        {
            *icomp = 3;
        }
        else
        {
            fprintf(stderr, "%s: cmpinc is 90 but channel=%s is weird",
                    __func__, kcmpnm);
            ierr = 1;
        }
    }
    else
    {
        fprintf(stderr, "%s: Can't classify component %s with cmpinc=%f\n",
                __func__, kcmpnm, cmpinc);
        ierr = 1;
    }
    if (ierr != 0)
    {
        fprintf(stderr, "%s: WARNING: icomp was set to vertical\n", __func__);
        *icomp = 1;
    }
    return ierr;
}
