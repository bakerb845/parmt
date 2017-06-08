#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "prepmt/prepmt_hpulse96.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

/*!
 * @brief Reads the ini file for the Computer Programs in Seismology hpulse96
 *        modeling variables.
 *
 * @param[in] iniFile     Name of ini file.
 * @param[in] section     Section of ini file with hpulse parameters.
 *
 * @param[out] parms      The hpulse96 parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_hpulse96_readHpulse96Parameters(const char *iniFile,
                                           const char *section,
                                           struct hpulse96_parms_struct *parms)
{
    const char *fcnm = "prepmt_hpulse96_readHpulse96Parameters\0";
    const char *s;
    char vname[256];
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHpulse96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        log_errorF("%s: ini file: %s does not exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        log_errorF("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:ntau", section);
    parms->ntau = iniparser_getint(ini, vname, -1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:alp", section);
    parms->alp = iniparser_getdouble(ini, vname, -1.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:ipt", section);
    parms->ipt = iniparser_getint(ini, vname, -1);
    if (parms->ipt == 4)
    {
        memset(vname, 0, 256*sizeof(char));
        sprintf(vname, "%s:rfile", section);
        s = iniparser_getstring(ini, vname, NULL);
        if (s != NULL)
        {
            strcpy(parms->rfile, s);
            if (!os_path_isfile(parms->rfile))
            {
                log_errorF("%s: Response file does not exist!\n", fcnm);
                return -1;
            }
        }
        else
        {
            log_errorF("%s: Response file not specified!\n", fcnm);
            return -1;
        }
    }
    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:idva", section);
    parms->idva = iniparser_getint(ini, vname, 1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:iodva", section);
    parms->iodva = iniparser_getint(ini, vname, -1);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:xmult", section);
    parms->xmult = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 256*sizeof(char));
    sprintf(vname, "%s:dozero", section);
    parms->dozero = iniparser_getboolean(ini, vname, 0);

    if (parms->ipt >= 0 && parms->ipt <= 1 && parms->ntau <= 0)
    {
        parms->ntau = 1;
    }
    if (parms->ipt == 2 && parms->alp <= 0.0){parms->alp = 1.0;}
    if (parms->ipt == 2 && parms->alp < 0.0)
    {
        log_errorF("%s: No alpha for Ohnaka pulse\n", fcnm);
        return -1;
    }
    if (parms->ipt < 0)
    {
        log_errorF("%s: No pulse shape defined\n", fcnm);
        return -1;
    }
    if (parms->iodva < 0){parms->iodva = parms->idva;}
    iniparser_freedict(ini);
    return 0;
}
