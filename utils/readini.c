#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <iniparser.h>
#include "parmt_utils.h"
#include "ttimes_config.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "compearth.h"

#define DEFAULT_PROJECT_NAME "parmt"

/*!
 * @brief Reads the general parmt parameters from the [general] section.
 *
 * @param[in] iniFile     Initialization file to read.
 *
 * @param[out] parms      The parmt general parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_readGeneralParms(const char *iniFile,
                                 struct parmtGeneralParms_struct *parms)
{
    const char *s;
    char *dirName;
    int ierr;
    dictionary *ini;
    ierr = 0;
    memset(parms, 0, sizeof(struct parmtGeneralParms_struct));
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file %s does not exist\n",
                __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    // project name
    s = iniparser_getstring(ini, "general:projectName", DEFAULT_PROJECT_NAME);
    if (s == NULL)
    {
        ierr = 1;
    }
    else
    {
        if (strlen(s) == 0){ierr = 1;}
    }
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error no project name\n", __func__);
        goto END;
    }
    strcpy(parms->projnm, s);
    // data file
    s = iniparser_getstring(ini, "general:dataFile\0", NULL);
    if (s == NULL)
    {
        sprintf(parms->dataFile, "%s.h5", parms->projnm);
    }
    else
    {
        strcpy(parms->dataFile, s);
    }
    if (!os_path_isfile(parms->dataFile))
    {
        fprintf(stderr, "%s: Error data file does %s not exist\n",
                __func__, parms->dataFile);
        ierr = 1;
        goto END;
    }
    s = iniparser_getstring(ini, "general:parmtArchive\0", NULL);
    if (s == NULL)
    {
        //printf("no default archive file yet - not done yet\n");
        strcpy(parms->parmtArchive, "parmtOutput.h5");
    } 
    else
    {
        strcpy(parms->parmtArchive, s);
    }
    dirName = os_dirname(parms->parmtArchive, &ierr);
    if (!os_path_isdir(dirName))
    {   
        ierr = os_makedirs(dirName);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to make directory %s\n",
                     __func__, dirName);
            goto END;
        }
    }
    memory_free8c(&dirName);

    s = iniparser_getstring(ini, "general:polarmtArchive\0", NULL);
    if (s == NULL)
    {
        printf("no polarmt file not done yet\n");
        strcpy(parms->polarmtArchive, "polarmtOutput.h5");
    }
    else
    {
        strcpy(parms->polarmtArchive, s);
    }
    dirName = os_dirname(parms->polarmtArchive, &ierr);
    if (!os_path_isdir(dirName))
    {
        ierr = os_makedirs(dirName);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to make directory %s\n",
                    __func__, dirName);
            goto END;
        }
    }
    memory_free8c(&dirName);

    s = iniparser_getstring(ini, "general:postmtFile\0", NULL);
    strcpy(parms->postmtFile, s);
/*
    // output file directory
    s = iniparser_getstring(ini, "general:resultsFileDirectory\0", "./");
    strcpy(parms->resultsDir, s);
    // output file suffix
    s = iniparser_getstring(ini, "general:resultsFileSuffix\0", NULL);
    if (s != NULL)
    {
        if (strlen(s) > 0){strcpy(parms->resultsFileSuffix, s);}
    }
    // polarity file suffix
    s = iniparser_getstring(ini, "general:polarityFileSuffix\0", NULL);
    if (s == NULL)
    {
        sprintf(parms->polarityFileSuffix, "%s_polarity",
                parms->resultsFileSuffix);
    }
    else
    {
        if (strlen(s) > 0)
        {
            strcpy(parms->polarityFileSuffix, s);
        }
        else
        {
            sprintf(parms->polarityFileSuffix, "%s_polarity",
                    parms->resultsFileSuffix);
        }
    }
*/
    // read the block size
    parms->blockSize = iniparser_getint(ini, "general:blockSize\0", 32);
    if (parms->blockSize < 1)
    {
        fprintf(stdout, "%s: blockSize = %d is invalid; setting to 32\n",
                __func__, parms->blockSize);
        parms->blockSize = 32;
    }
    // determine if i want to use lags or not
    parms->lwantLags = iniparser_getboolean(ini, "general:lwantLags\0", 0);
    if (parms->lwantLags)
    {
        parms->defaultMaxLagTime
           = iniparser_getdouble(ini, "general:defaultMaxLagTime", 2.0);
        if (parms->defaultMaxLagTime < 0.0)
        {
            fprintf(stdout, "%s: Invalid lag time %f - setting to 2\n",
                    __func__, parms->defaultMaxLagTime);
        }
    }
    // determine the objective function
    parms->objFnType = iniparser_getint(ini, "general:objFnType\0", 1);
    if (parms->objFnType < 1 || parms->objFnType > 3)
    {
        fprintf(stderr, "%s: Invalid obj fn type %d\n",
                 __func__, parms->objFnType);
        return -1;
    }
    // determine if the l1 lag should be rescaled
    parms->lrescale = iniparser_getboolean(ini, "general:lrescale\0", 0);
END:;
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Reads the parameters necessary for computing the polarities in
 *        a global 1D earth. 
 *
 * @param[in] iniFile    Initialization file to read.
 *
 * @param[out] parms     Polarity grid search parameters. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_readPolarityParms(const char *iniFile,
                                  struct parmtPolarityParms_struct *parms)
{
    const char *s;
    int ierr;
    dictionary *ini;
    ierr = 0;
    memset(parms, 0, sizeof(struct parmtPolarityParms_struct));
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file %s does not exist\n",
                __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    // do i even want to do this analysis?
    if (!parms->lcomputePolarity)
    {
        strcpy(parms->ttimesTablesDir, TTIMES_DEFAULT_TABLE_DIRECTORY);
        strcpy(parms->ttimesModel, TTIMES_DEFAULT_MODEL);
        goto END;
    }
    // ttimes tables directory
    s = iniparser_getstring(ini, "ttimes:ttimesTableDir\0",
                            TTIMES_DEFAULT_TABLE_DIRECTORY);
    if (!os_path_isdir(s))
    {   
        fprintf(stderr, "%s: Error ttimes tables directory doesn't exist\n",
                __func__);
    }   
    else
    {   
        strcpy(parms->ttimesTablesDir, s); 
    }   
    // ttimes model
    s = iniparser_getstring(ini, "ttimes:ttimesModel\0", TTIMES_DEFAULT_MODEL);
    if (s != NULL)
    {   
        strcpy(parms->ttimesModel, s); 
    }   
    else
    {   
        fprintf(stderr, "%s: Error ttimesModel not specified\n", __func__);
    }

END:;
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Reads the moment tensor grid search parameters.
 *
 * @param[in] iniFile    Initialization file to read.
 *
 * @param[out] parms     MT grid search parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int parmt_utils_readMtSearch(const char *iniFile, 
                             struct parmtMtSearchParms_struct *parms)
{
    double betaLower, betaUpper, mw;
    int ierr;
    dictionary *ini;
    bool isMw, lrescale;
    const double pi180 = M_PI/180.0;
    ierr = 0;
    memset(parms, 0, sizeof(struct parmtMtSearchParms_struct));
    if (!os_path_isfile(iniFile))
    {   
        fprintf(stderr, "%s: Error ini file %s does not exist\n",
                __func__, iniFile);
        return -1; 
    }   
    ini = iniparser_load(iniFile);
    // Colatitudes
    parms->betaLower =-90.0;
    parms->betaUpper = 90.0;
    parms->nb = iniparser_getint(ini, "mtsearch:nlat", -1);
    if (parms->nb ==-1)
    {
        fprintf(stdout, "%s: Setting five latitudes in grid search\n",
                __func__);
        parms->nb = 5;
    }
    else
    {
        // Try to read the low 
        parms->betaLower = iniparser_getdouble(ini, "mtsearch:latLower",-90.0);
        parms->betaUpper = iniparser_getdouble(ini, "mtsearch:latUpper", 90.0);
        if (parms->betaLower <-90.0)
        {
            fprintf(stdout, "%s: Overriding mtlat_lower to -90 degrees\n",
                     __func__);
            parms->betaLower =-90.0;
        }
        if (parms->betaUpper > 90.0)
        {
            fprintf(stdout, "%s: Overriding mtlat_upper to 90 degrees\n",
                    __func__);
            parms->betaUpper = 90.0;
        }
    }
    // Reverse so that 180 is max and 0 is min
    betaUpper = parms->betaUpper;
    betaLower = parms->betaLower;
    parms->betaLower = (90.0 - betaUpper)*pi180; // lat -> colatitude
    parms->betaUpper = (90.0 - betaLower)*pi180; // lat -> colatitude
    if (parms->nb == 1){parms->betaUpper = parms->betaLower;}
    // Longitudes
    parms->gammaLower =-30.0;
    parms->gammaUpper = 30.0;
    parms->ng = iniparser_getint(ini, "mtsearch:nlon", -1);
    if (parms->ng ==-1)
    {
        fprintf(stdout, "%s: Setting five longitudes in grid search\n",
                __func__);
        parms->ng = 5;
    }
    else
    {
        // Try to read the low 
        parms->gammaLower = iniparser_getdouble(ini,
                                                "mtsearch:lonLower",-30.0);
        parms->gammaUpper = iniparser_getdouble(ini,
                                                "mtsearch:lonUpper", 30.0);
        if (parms->gammaLower <-30.0)
        {
            fprintf(stdout, "%s: Overriding mtlon_lower to -30 degrees\n",
                    __func__);
            parms->gammaLower =-30.0;
        }
        if (parms->gammaUpper > 30.0)
        {
            fprintf(stdout, "%s: Overriding mtlon_upper to 30 degrees\n",
                    __func__);
            parms->gammaUpper = 30.0;
        }
    }
    parms->gammaLower = parms->gammaLower*pi180;
    parms->gammaUpper = parms->gammaUpper*pi180;
    if (parms->ng == 1){parms->gammaUpper = parms->gammaLower;}
    // Strike angles
    parms->kappaLower = 0.0;
    parms->kappaUpper = 360.0;
    parms->nk = iniparser_getint(ini, "mtsearch:nstrike", -1);
    if (parms->nk ==-1)
    {
        fprintf(stdout, "%s: Setting five strike angles in grid search\n",
                __func__);
        parms->nk = 5;
    }
    else
    {
        // Try to read the low 
        parms->kappaLower = iniparser_getdouble(ini,
                                                "mtsearch:strikeLower", 0.0);
        parms->kappaUpper = iniparser_getdouble(ini,
                                                "mtsearch:strikeUpper", 360.0);
        if (parms->kappaLower < 0.0)
        {
            fprintf(stdout, "%s: Overriding strike_lower to 0 degrees\n",
                    __func__);
            parms->kappaLower = 0.0;
        }
        if (parms->kappaUpper > 360.0)
        {
            fprintf(stdout, "%s: Overriding strike_upper to 360 degrees\n",
                    __func__);
            parms->kappaUpper = 360.0;
        }
    }
    parms->kappaLower = parms->kappaLower*pi180;
    parms->kappaUpper = parms->kappaUpper*pi180;
    if (parms->nk == 1){parms->kappaUpper = parms->kappaLower;}
    // Dip angles
    parms->thetaLower =  0.0;
    parms->thetaUpper = 90.0;
    parms->nt = iniparser_getint(ini, "mtsearch:ndip", -1);
    if (parms->nt ==-1)
    {
        fprintf(stdout, "%s: Setting five dip angles in grid search\n",
                __func__);
        parms->nt = 5;
    }
    else
    {
        // Try to read the low 
        parms->thetaLower = iniparser_getdouble(ini,
                                                "mtsearch:thetaLower", 0.0);
        parms->thetaUpper = iniparser_getdouble(ini,
                                                "mtsearch:thetaUpper", 90.0);
        if (parms->thetaLower < 0.0)
        {
            fprintf(stdout, "%s: Overriding theta_lower to 0 degrees\n",
                    __func__);
            parms->thetaLower = 0.0;
        }
        if (parms->thetaUpper > 90.0)
        {
            fprintf(stdout, "%s: Overriding theta_upper to 90 degrees\n",
                    __func__);
            parms->thetaUpper = 90.0;
        }
    }
    parms->thetaLower = parms->thetaLower*pi180;
    parms->thetaUpper = parms->thetaUpper*pi180;
    if (parms->nt == 1){parms->thetaUpper = parms->thetaLower;}
    // Rake angles
    parms->sigmaLower =-90.0;
    parms->sigmaUpper = 90.0;
    parms->ns = iniparser_getint(ini, "mtsearch:nrake", -1);
    if (parms->ns ==-1)
    {
        fprintf(stdout, "%s: Setting five rake angles in grid search\n",
                __func__);
        parms->ns = 5;
    }
    else
    {
        // Try to read the low 
        parms->sigmaLower = iniparser_getdouble(ini,
                                                "mtsearch:rakeLower",-90.0);
        parms->sigmaUpper = iniparser_getdouble(ini,
                                                "mtsearch:rakeUpper", 90.0);
        if (parms->sigmaLower < -90.0)
        {
            fprintf(stdout, "%s: Overriding rake_lower to -90 degrees\n",
                    __func__);
            parms->sigmaLower =-90.0;
        }
        if (parms->sigmaUpper > 90.0)
        {
            fprintf(stdout, "%s: Overriding strike_upper to 90 degrees\n",
                    __func__);
            parms->sigmaUpper = 90.0;
        }
    }
    parms->sigmaLower = parms->sigmaLower*pi180;
    parms->sigmaUpper = parms->sigmaUpper*pi180;
    if (parms->ns == 1){parms->sigmaUpper = parms->sigmaLower;}
    // Scalar moment
    isMw = true;
    parms->m0Lower = 5.0;
    parms->m0Upper = 6.0;
    parms->nm = iniparser_getint(ini, "mtsearch:nm", -1);
    if (parms->nm < 1)
    {
        fprintf(stdout, "%s: Setting one scalar moment to mag 5.5\n", __func__);
        parms->nm = 1;
    }
    else
    {
        isMw = iniparser_getboolean(ini, "parmt:isMw", 1);
        parms->m0Lower = iniparser_getdouble(ini,
                                             "mtsearch:m0Lower", 5.0);
        parms->m0Upper = iniparser_getdouble(ini,
                                             "mtsearch:m0Upper", 6.0);
    }
    if (parms->nm == 1){parms->m0Upper = parms->m0Lower;}
    if (isMw)
    {
        //parms->m0Lower = pow(10.0, 1.5*parms->m0Lower + 9.1);
        //parms->m0Upper = pow(10.0, 1.5*parms->m0Upper + 9.1);
        mw = parms->m0Lower; 
        compearth_mw2m0(1, CE_KANAMORI_1978, &mw, &parms->m0Lower);
        mw = parms->m0Upper;
        compearth_mw2m0(1, CE_KANAMORI_1978, &mw, &parms->m0Upper);
    }
    parms->luseLog = iniparser_getboolean(ini, "mtsearch:luseLog", false);

    // determine if i want to rescale the synthetics to the observations or ont
    lrescale = iniparser_getboolean(ini, "general:lrescale\0", 0);
    if (lrescale)
    {   
        if (parms->nm > 1)
        {
            fprintf(stdout, "%s: lrescale invalidates magnitude!\n", __func__);
            fprintf(stdout, "%s: nm should be set to 1!\n", __func__);
        }
    }

    iniparser_freedict(ini);
    return ierr;
}
