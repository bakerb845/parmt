#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "parmt_postProcess.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "compearth.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"


static void getBaseAndExp(const double val, double *base, int *exp);
static void setFillColor(const int i, const int iopt, char color[32]);

/*!
 * @brief Writes the global station distribution of stations, the
 *        station names, the epicenter or moment tensor, and, if
 *        desired, some indication of station polarity.
 */
int postmt_gmtHelper_writeGlobalMap(
    struct globalMapOpts_struct globalMap,
    const int nobs, const struct sacData_struct *data)
{
    const char *fcnm = "postmt_gmtHelper_writeGlobalMap\0";
    FILE *ofl;
    char *dirName, line[256], cpick[8];
    int *pol, i, ierr, l, nw, nwd, nwu;
    size_t lenos;
    const char *forwardSlash = "/";
    const int nTimeVars = 11;
    const enum sacHeader_enum timeVarNames[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};

    dirName = os_dirname(globalMap.outputScript);
    if (!os_path_isdir(dirName))
    {
        ierr = os_makedirs(dirName); 
        if (ierr != 0)
        {
            printf("%s: Failed to make output directory: %s\n", fcnm, dirName);
            return -1;
        }
    }
    if (dirName != NULL){free(dirName);}

    ofl = fopen(globalMap.outputScript, "w"); 
    fprintf(ofl, "#!/bin/bash\n");
    fprintf(ofl, "outps=%s\n", globalMap.psFile);
    fprintf(ofl, "olat=%f\n", globalMap.evla);
    fprintf(ofl, "olon=%f\n", globalMap.evlo);
    fprintf(ofl, "J=-JH${olon}%s6i\n", forwardSlash);
    fprintf(ofl, "R=-Rg\n");
    fprintf(ofl, "gmt pscoast $J $R -B0g30 -Di -Ggray -P -K > ${outps}\n");
    fprintf(ofl, "ts=0.15i\n");
  
    fprintf(ofl, "# Draw great circle arcs between source and receivers\n");
    fprintf(ofl, "gmt psxy $J $R -W1p -O -K << EOF >> ${outps}\n");
    for (i=0; i<nobs; i++)
    {   
        fprintf(ofl, "%8.3f %8.3f\n", data[i].header.stlo, data[i].header.stla);
        fprintf(ofl, "%8.3f %8.3f\n", globalMap.evlo, globalMap.evla);
        if (i < nobs - 1){fprintf(ofl, "\n");}
    }   
    fprintf(ofl, "EOF\n");

    fprintf(ofl, "# Plot the moment tensor\n");
    postmt_gmtHelper_makePsmecaLine(globalMap.basis, globalMap.mts,
                                    globalMap.evla, globalMap.evlo,
                                    globalMap.evdp, "Optimum\0", line);
    if (globalMap.lwantMT)
    {
        fprintf(ofl, "gmt psmeca $J $R -Sm0.15i -M -N -Gblue -W0.5p,black -O -K << EOF >> ${outps}\n");
        fprintf(ofl, "%s", line);
        fprintf(ofl, "EOF\n");
    }
    else
    {
        fprintf(ofl, "# Plot the epicenter\n");
        fprintf(ofl, "gmt psxy $J $R -Sa${ts} -Gblue -Wblack -O -K << EOF >> ${outps}\n");
        fprintf(ofl, "%8.3f %8.3f\n", globalMap.evlo, globalMap.evla);
        fprintf(ofl, "EOF\n");
    }

    if (!globalMap.lwantPolarity)
    {
        fprintf(ofl, "# Plot the stations\n");
        fprintf(ofl, "gmt psxy $J $R -St${ts} -Gred -Wblack -O -K << EOF >> ${outps}\n");
        for (i=0; i<nobs; i++)
        {
            fprintf(ofl, "%8.3f %8.3f %s\n", data[i].header.stlo,
                    data[i].header.stla, data[i].header.kstnm);
        }
        fprintf(ofl, "EOF\n");
    }
    else
    {
        // Get the polarities
        nw = 0;
        nwd = 0;
        nwu = 0;
        pol = (int *) calloc((size_t) nobs, sizeof(int));
        for (i=0; i<nobs; i++)
        {
            memset(cpick, 0, 8*sizeof(char));
            for (l=0; l<nTimeVars; l++)
            {
                ierr = sacio_getCharacterHeader(timeVarNames[l],
                                                data[i].header, cpick);
                if (ierr == 0){break;}
            }
            lenos = strlen(cpick);
            if (lenos > 0)
            {
                if (cpick[lenos-1] == '+')
                {
                    nwu = nwu + 1;
                    pol[i] = 1;
                }
                else if (cpick[lenos-1] == '-')
                {
                    nwd = nwd + 1;
                    pol[i] =-1;
                }
                else
                {
                    nw = nw + 1;
                }
            }
            else
            {
                nw = nw + 1;
            }
        }
        // Write the unknowns
        if (nw > 0)
        {
            fprintf(ofl, "# Plot the indeterminant stations\n");
            fprintf(ofl, "gmt psxy $J $R -Ss${ts} -Gred -Wblack -O -K << EOF >> ${outps}\n");
            for (i=0; i<nobs; i++)
            {
                if (pol[i] == 0)
                {
                    fprintf(ofl, "%8.3f %8.3f %s\n", data[i].header.stlo,
                            data[i].header.stla, data[i].header.kstnm);
                }
            }   
            fprintf(ofl, "EOF\n");
        }
        // Write the ups
        if (nwu > 0)
        {
            fprintf(ofl, "# Plot the upward stations\n");
            fprintf(ofl, "gmt psxy $J $R -St${ts} -Gred -Wblack -O -K << EOF >> ${outps}\n");
            for (i=0; i<nobs; i++)
            {
                if (pol[i] == 1)
                {
                    fprintf(ofl, "%8.3f %8.3f %s\n", data[i].header.stlo,
                            data[i].header.stla, data[i].header.kstnm);
                }
            }
            fprintf(ofl, "EOF\n");
        }
        // Write the downs
        if (nwd > 0)
        {
            fprintf(ofl, "# Plot the down stations\n");
            fprintf(ofl, "gmt psxy $J $R -Si${ts} -Gred -Wblack -O -K << EOF >> ${outps}\n");
            for (i=0; i<nobs; i++)
            {
                if (pol[i] ==-1)
                {
                    fprintf(ofl, "%8.3f %8.3f %s\n", data[i].header.stlo,
                            data[i].header.stla, data[i].header.kstnm);
                }
            }
            fprintf(ofl, "EOF\n");
        }
        free(pol);
    }

    fprintf(ofl, "# Plot the station names\n");
    fprintf(ofl, "gmt pstext $J $R -F+f8p+a0+jCT -D0%s-4p -O << EOF >> ${outps}\n",
            forwardSlash);
    for (i=0; i<nobs; i++)
    {
        fprintf(ofl, "%8.3f %8.3f %s\n", data[i].header.stlo,
                data[i].header.stla, data[i].header.kstnm);
    }
    fprintf(ofl, "EOF\n");
    fprintf(ofl, "gmt psconvert -A -Tj ${outps}\n");
    fprintf(ofl, "rm ${outps}\n");
    fclose(ofl);
    chmod(globalMap.outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeThetaBoxes(const bool lappend, const bool lclose,
                                     const char *outputScript,
                                     const char *psFile,
                                     const int nt,
                                     const int ithetaOpt,
                                     const double *__restrict__ thetas,
                                     const double *__restrict__ thetaHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumTheta, *h, hAvg, thetaAvg, tmax;
    int i, ierr;
    const bool lwritePrior = true;
    tmax = array_max64f(nt, thetaHist);
    // Compute h
    h = memory_calloc64f(nt);
    compearth_theta2h(nt, thetas, h);
    thetaAvg = -30.0;
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R0/90/0/%.2f -Bg10a10:\"Dips (deg)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, tmax*1.1, tmax*0.2, app, more);
    setFillColor(0, ithetaOpt, color);
    if (nt > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%f %f\n", 0.0, 0.0);
    fprintf(ofl, "%f %f\n", 0.0, thetaHist[0]);
    for (i=0; i<nt-1; i++)
    {
        hAvg = 0.5*(h[i] + h[i+1]);
        compearth_h2theta(1, &hAvg, &thetaAvg);
        fprintf(ofl, "%f %f\n", thetaAvg*180.0/M_PI, thetaHist[i]);
        fprintf(ofl, "%f %f\n", thetaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, ithetaOpt, color);
        if (i < nt - 2 || lwritePrior)
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%f %f\n", thetaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "%f %f\n", thetaAvg*180.0/M_PI, thetaHist[i+1]);
    }
    fprintf(ofl, "%f %f\n", 90.0, thetaHist[nt-1]);
    fprintf(ofl, "%f %f\n", 90.0, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {   
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n", 00.0, 1.0/(double) nt);
        fprintf(ofl, "%f %f\n", 90.0, 1.0/(double) nt);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R0/90/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n");
    cumTheta = array_cumsum64f(nt, thetaHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<nt; i++)
    {   
        fprintf(ofl, "%f %f\n", thetas[i]*180/M_PI, cumTheta[i]);
    }   
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&cumTheta);
    memory_free64f(&h);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeSigmaBoxes(const bool lappend, const bool lclose,
                                     const char *outputScript,
                                     const char *psFile,
                                     const int ns,
                                     const int isigmaOpt,
                                     const double *__restrict__ sigmas,
                                     const double *__restrict__ sigmaHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumSigma, ds, smax;
    int i, ierr;
    const bool lwritePrior = true;
    smax = array_max64f(ns, sigmaHist);
    // Take averages
    ds = 0.0;
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
           "gmt psbasemap -JX5i/1i %s -R-90/90/0/%.2f -Bg15a15:\"Slips (deg)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, smax*1.1, smax*0.2, app, more);
    setFillColor(0, isigmaOpt, color);
    if (ns > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
       ds = (sigmas[1] - sigmas[0])*180.0/M_PI;
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%.2f %f\n", -90.0, 0.0);
    fprintf(ofl, "%.2f %f\n", -90.0, sigmaHist[0]);
    for (i=0; i<ns-1; i++)
    {   
        fprintf(ofl, "%.2f %f\n", -90.0 + (double) (i+1)*ds, sigmaHist[i]);
        fprintf(ofl, "%.2f %f\n", -90.0 + (double) (i+1)*ds, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, isigmaOpt, color);
        if (i < ns - 2 || lwritePrior)
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%.2f %f\n", -90.0 + (double) (i+1)*ds, 0.0);
        fprintf(ofl, "%.2f %f\n", -90.0 + (double) (i+1)*ds, sigmaHist[i+1]);
    }
    fprintf(ofl, "%.2f %f\n", 90.0, sigmaHist[ns-1]);
    fprintf(ofl, "%.2f %f\n", 90.0, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {   
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n",-90.0, 1.0/(double) ns);
        fprintf(ofl, "%f %f\n", 90.0, 1.0/(double) ns);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R-90/90/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n");
    cumSigma = array_cumsum64f(ns, sigmaHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<ns; i++)
    {
        fprintf(ofl, "%f %f\n", sigmas[i]*180/M_PI, cumSigma[i]);
    }   
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&cumSigma);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}

//============================================================================//
int postmt_gmtHelper_writeKappaBoxes(const bool lappend, const bool lclose,
                                     const char *outputScript,
                                     const char *psFile,
                                     const int nk,
                                     const int kappaOpt,
                                     const double *__restrict__ kappas,
                                     const double *__restrict__ kappaHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumKappa, dk, kmax;
    int i, ierr;
    const bool lwritePrior = true;
    kmax = array_max64f(nk, kappaHist);
    // Take averages
    dk = 0.0;
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R0/360/0/%.2f -Bg30a30:\"Strike (deg)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, kmax*1.1, kmax*0.2, app, more);
    setFillColor(0, kappaOpt, color);
    if (nk > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
       dk = (kappas[1] - kappas[0])*180.0/M_PI;
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%.2f %f\n", 0.0, 0.0);
    fprintf(ofl, "%.2f %f\n", 0.0, kappaHist[0]);
    for (i=0; i<nk-1; i++)
    {
        fprintf(ofl, "%.2f %f\n", 0.0 + (double) (i+1)*dk, kappaHist[i]);
        fprintf(ofl, "%.2f %f\n", 0.0 + (double) (i+1)*dk, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, kappaOpt, color);
        if (i < nk - 2 || lwritePrior)
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%.2f %f\n", 0.0 + (double) (i+1)*dk, 0.0);
        fprintf(ofl, "%.2f %f\n", 0.0 + (double) (i+1)*dk, kappaHist[i+1]);
    }   
    fprintf(ofl, "%.2f %f\n", 360.0, kappaHist[nk-1]);
    fprintf(ofl, "%.2f %f\n", 360.0, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {   
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n", 0.0,   1.0/(double) nk);
        fprintf(ofl, "%f %f\n", 360.0, 1.0/(double) nk);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R0/360/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n");
    cumKappa = array_cumsum64f(nk, kappaHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<nk; i++)
    {   
        fprintf(ofl, "%f %f\n", kappas[i]*180.0/M_PI, cumKappa[i]);
    }   
    fprintf(ofl, "EOF\n");
    if (lclose)
    {
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&cumKappa);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeDepthBoxes(const bool lappend, const bool lclose,
                                     const char *outputScript,
                                     const char *psFile,
                                     const int nd,
                                     const int idepOpt,
                                     const double *__restrict__ deps,
                                     const double *__restrict__ depHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumDep, dd, dmax, depMin, depMax;
    int i, ierr;
    const bool lwritePrior = true;
    // Compute moment magnitudes 
    dmax = array_max64f(nd, depHist);
    dd = 0.1;
    if (nd > 1){dd = deps[1] - deps[0];};
    depMin = fmax(0.0, deps[0]    - dd/2.0);
    depMax = fmax(0.0, deps[nd-1] + dd/2.0);
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R%.1f/%.1f/0/%.2f -Bg%fa%f:\"Depths (km)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, depMin, depMax, dmax*1.1, dd, (nd-1)*dd/5.0, dmax*0.2, app, more);
    setFillColor(0, idepOpt, color);
    if (nd > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%.2f %f\n", depMin, 0.0);
    fprintf(ofl, "%.2f %f\n", depMin, depHist[0]);
    for (i=0; i<nd-1; i++)
    {
        fprintf(ofl, "%.2f %f\n", depMin + (double) (i+1)*dd, depHist[i]);
        fprintf(ofl, "%.2f %f\n", depMin + (double) (i+1)*dd, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, idepOpt, color);
        if (i < nd - 2 || lwritePrior)
        {   
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {   
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%.2f %f\n", depMin + (double) (i+1)*dd, 0.0);
        fprintf(ofl, "%.2f %f\n", depMin + (double) (i+1)*dd, depHist[i+1]);
    }
    fprintf(ofl, "%.2f %f\n", depMax, depHist[nd-1]);
    fprintf(ofl, "%.2f %f\n", depMax, 0.0);
    fprintf(ofl, "EOF\n"); 
    // Write the prior distribution
    if (lwritePrior) 
    {   
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n", depMin, 1.0/(double) nd);
        fprintf(ofl, "%f %f\n", depMax, 1.0/(double) nd);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R%f/%f/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n", depMin, depMax);
    cumDep = array_cumsum64f(nd, depHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<nd; i++)
    {
        fprintf(ofl, "%f %f\n", deps[i], cumDep[i]);
    }
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&cumDep);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeMagnitudeBoxes(const bool lappend, const bool lclose,
                                         const char *outputScript,
                                         const char *psFile,
                                         const int nm,
                                         const int magOpt,
                                         const double *__restrict__ M0s,
                                         const double *__restrict__ magHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumMag, *Mw, dm, mmax, mwMin, mwMax;
    int i, ierr;
    const bool lwritePrior = true;
    // Compute moment magnitudes 
    Mw = memory_calloc64f(nm);
    mmax = array_max64f(nm, magHist);
    compearth_m02mw(nm, KANAMORI_1978, M0s, Mw);
    dm = 0.1;
    if (nm > 1){dm = Mw[1] - Mw[0];};
    mwMin = Mw[0] - dm/2.0;
    mwMax = Mw[nm-1] + dm/2.0;
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R%.2f/%.2f/0/%.2f -Bg%fa%f:\"Magnitude (Mw)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, mwMin, mwMax, mmax*1.1, dm, (nm-1)*dm/5.0, mmax*0.2, app, more);
    setFillColor(0, magOpt, color);
    if (nm > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }   
    fprintf(ofl, "%.2f %f\n", mwMin, 0.0);
    fprintf(ofl, "%.2f %f\n", mwMin, magHist[0]);
    for (i=0; i<nm-1; i++)
    {
        fprintf(ofl, "%.2f %f\n", mwMin + (double) (i+1)*dm, magHist[i]);
        fprintf(ofl, "%.2f %f\n", mwMin + (double) (i+1)*dm, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, magOpt, color);
        if (i < nm - 2 || lwritePrior)
        {   
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%.2f %f\n", mwMin + (double) (i+1)*dm, 0.0);
        fprintf(ofl, "%.2f %f\n", mwMin + (double) (i+1)*dm, magHist[i+1]);
    }
    fprintf(ofl, "%.2f %f\n", mwMax, magHist[nm-1]);
    fprintf(ofl, "%.2f %f\n", mwMax, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {   
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n", mwMin, 1.0/(double) nm);
        fprintf(ofl, "%f %f\n", mwMax, 1.0/(double) nm);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R%f/%f/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n", mwMin, mwMax);
    cumMag = array_cumsum64f(nm, magHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<nm; i++)
    {
        fprintf(ofl, "%f %f\n", Mw[i], cumMag[i]);
    }
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    fclose(ofl);
    memory_free64f(&cumMag);
    memory_free64f(&Mw);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeGammaBoxes(const bool lappend, const bool lclose,
                                     const char *outputScript,
                                     const char *psFile,
                                     const int ng,
                                     const int igammaOpt,
                                     const double *__restrict__ gammas,
                                     const double *__restrict__ gammaHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *cumGamma, *v, gmax, gammaAvg, vAvg;
    int i, ierr;
    const bool lwritePrior = true;
    gmax = array_max64f(ng, gammaHist);
    // Compute v
    v = memory_calloc64f(ng);
    compearth_gamma2v(ng, gammas, v);
    // Take averages
    gammaAvg = -30.0;
    memset(app, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R-30/30/0/%.2f -Bg5a5:\"Longitude (deg)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, gmax*1.1, gmax*0.2, app, more);
    setFillColor(0, igammaOpt, color);
    if (ng > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%f %f\n", -30.0, 0.0);
    fprintf(ofl, "%f %f\n", -30.0, gammaHist[0]);
    for (i=0; i<ng-1; i++)
    {
        vAvg = 0.5*(v[i] + v[i+1]);
        compearth_v2gamma(1, &vAvg, &gammaAvg);
        fprintf(ofl, "%f %f\n", gammaAvg*180.0/M_PI, gammaHist[i]);
        fprintf(ofl, "%f %f\n", gammaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, igammaOpt, color);
        if (i < ng - 2 || lwritePrior)
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%f %f\n", gammaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "%f %f\n", gammaAvg*180.0/M_PI, gammaHist[i+1]);
    }
    fprintf(ofl, "%f %f\n", 30.0, gammaHist[ng-1]);
    fprintf(ofl, "%f %f\n", 30.0, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n",-30.0, 1.0/(double) ng);
        fprintf(ofl, "%f %f\n", 30.0, 1.0/(double) ng);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF
    fprintf(ofl, "gmt psbasemap -R-30/30/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n");
    cumGamma = array_cumsum64f(ng, gammaHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n");
    }
    for (i=0; i<ng; i++)
    {   
        fprintf(ofl, "%f %f\n", gammas[i]*180.0/M_PI, cumGamma[i]);
    }   
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&cumGamma);
    memory_free64f(&v);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_writeBetaBoxes(const bool lappend, const bool lclose,
                                    const char *outputScript,
                                    const char *psFile,
                                    const int nb,
                                    const int ibetaOpt,
                                    const double *__restrict__ betas,
                                    const double *__restrict__ betaHist)
{
    FILE *ofl;
    char color[32], app[8], more[8], shift[8];
    double *betaCum, *u, bmax, betaAvg, uAvg;
    int i, ierr;
    const bool lwritePrior = true;
    bmax = array_max64f(nb, betaHist);
    // Compute u
    u = memory_calloc64f(nb);
    compearth_beta2u(nb, betas, u);
    // Take averages
    betaAvg = 0.0;
    memset(app, 0, 8*sizeof(char));
    memset(more, 0, 8*sizeof(char));
    memset(shift, 0, 8*sizeof(char));
    if (!lappend)
    {
        ofl = fopen(outputScript, "w");
        fprintf(ofl, "#!/bin/bash\n");
        fprintf(ofl, "gmt gmtset FONT_LABEL 12p\n");
        fprintf(ofl, "gmt gmtset MAP_LABEL_OFFSET 0.1c\n");
        fprintf(ofl, "psfl=%s\n", psFile);
    }
    else
    {
        ofl = fopen(outputScript, "a");
        strcpy(app, "-O\0");
        strcpy(more, ">\0");
        strcpy(shift, "-Y4.0\0");
    }
    fprintf(ofl,
            "gmt psbasemap -JX5i/1i %s -R-90/90/0/%.3f -Bg15a15:\"Latitude (deg)\":/a%.2f:\"Likelihood\":WSn -P %s -K >%s ${psfl}\n",
            shift, bmax*1.1, bmax*0.2, app, more);
    setFillColor(0, ibetaOpt, color);
    if (nb > 1)
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
    }
    else
    {
       fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
    }
    fprintf(ofl, "%f %f\n", 90.0 - 0.0, 0.0);
    fprintf(ofl, "%f %f\n", 90.0 - 0.0, betaHist[0]); 
    for (i=0; i<nb-1; i++)
    {
        uAvg = 0.5*(u[i] + u[i+1]);
        compearth_u2beta(1, 20, 2, &uAvg, 1.e-7, &betaAvg);
        fprintf(ofl, "%f %f\n", 90.0 - betaAvg*180.0/M_PI, betaHist[i]);
        fprintf(ofl, "%f %f\n", 90.0 - betaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "EOF\n");
        setFillColor(i+1, ibetaOpt, color);
        if (i < nb - 2 || lwritePrior)
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O -K << EOF >> ${psfl}\n", color);
        }
        else
        {
            fprintf(ofl, "gmt psxy -R -J -Wblack %s -O << EOF >> ${psfl}\n", color);
        }
        fprintf(ofl, "%f %f\n", 90.0 - betaAvg*180.0/M_PI, 0.0);
        fprintf(ofl, "%f %f\n", 90.0 - betaAvg*180.0/M_PI, betaHist[i+1]);
    }
    fprintf(ofl, "%f %f\n", 90.0 - M_PI*180.0/M_PI, betaHist[nb-1]);
    fprintf(ofl, "%f %f\n", 90.0 - M_PI*180.0/M_PI, 0.0);
    fprintf(ofl, "EOF\n");
    // Write the prior distribution
    if (lwritePrior)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,black -O -K << EOF >> ${psfl}\n");
        fprintf(ofl, "%f %f\n",-90.0, 1.0/(double) nb);
        fprintf(ofl, "%f %f\n", 90.0, 1.0/(double) nb);
        fprintf(ofl, "EOF\n");
    }
    // Write the CDF 
    fprintf(ofl, "gmt psbasemap -R-90/90/0/1.05 -J -Bp0.2/a0.2:\"CDF\":E -O -K >> ${psfl}\n");
    betaCum = array_cumsum64f(nb, betaHist, &ierr);
    if (lclose)
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O << EOF >> ${psfl}\n");
    }
    else
    {
        fprintf(ofl, "gmt psxy -R -J -W1,blue -O -K << EOF >> ${psfl}\n"); 
    }
    for (i=0; i<nb; i++)
    {
        // the cdf is written backwards because i flip the axis to make it
        // increase left to right; this is opposite of what colatitude wants
        // to do
        fprintf(ofl, "%f %f\n", 90.0 - betas[i]*180.0/M_PI, betaCum[nb-1-i]);
    }   
    fprintf(ofl, "EOF\n");
    if (lclose)
    {    
        //fprintf(ofl, "gmt psconvert -A -Tj ${psfl}\n");
        //fprintf(ofl, "rm ${psfl}\n");
    }
    memory_free64f(&u);
    memory_free64f(&betaCum);
    fclose(ofl);
    chmod(outputScript, 0755);
    return 0;
}
//============================================================================//
int postmt_gmtHelper_makeRegularHistograms(
    const int nlocs, const int nm,
    const int nb, const int ng, const int nk,
    const int ns, const int nt,
    const int nmt, const double *__restrict__ phi,
    double **__restrict__ locHist,
    double **__restrict__ magHist,
    double **__restrict__ betaHist,
    double **__restrict__ gammaHist,
    double **__restrict__ kappaHist,
    double **__restrict__ sigmaHist,
    double **__restrict__ thetaHist)
{
    double *betaH, *gammaH, *kappaH, *locH, *magH, *sigmaH, *thetaH, xsum;
    int ib, ierr, ig, ik, il, im, imt, is, it; 
    locH = memory_calloc64f(nlocs);
    magH = memory_calloc64f(nm);
    betaH = memory_calloc64f(nb);
    gammaH = memory_calloc64f(ng);
    kappaH = memory_calloc64f(nk);
    sigmaH = memory_calloc64f(ns);
    thetaH = memory_calloc64f(nt);
    for (il=0; il<nlocs; il++)
    {
        for (im=0; im<nm; im++)
        {
            for (ib=0; ib<nb; ib++)
            {
                for (ig=0; ig<ng; ig++)
                {
                    for (ik=0; ik<nk; ik++)
                    {
                        for (is=0; is<ns; is++)
                        {
                            for (it=0; it<nt; it++)
                            {
                                imt = il*nm*nb*ng*nk*ns*nt
                                    + im*nb*ng*nk*ns*nt
                                    + ib*ng*nk*ns*nt
                                    + ig*nk*ns*nt
                                    + ik*ns*nt
                                    + is*nt
                                    + it;
                                locH[il]   = locH[il]   + phi[imt];
                                magH[im]   = magH[im]   + phi[imt];
                                betaH[ib]  = betaH[ib]  + phi[imt];
                                gammaH[ig] = gammaH[ig] + phi[imt];
                                kappaH[ik] = kappaH[ik] + phi[imt];
                                sigmaH[is] = sigmaH[is] + phi[imt];
                                thetaH[it] = thetaH[it] + phi[imt];
                            }
                        }
                    }
                }
            }
        }
    }
    xsum = array_sum64f(nmt, phi, &ierr);
    cblas_dscal(nlocs, 1.0/xsum, locH, 1);
    cblas_dscal(nm, 1.0/xsum, magH, 1);
    cblas_dscal(nb, 1.0/xsum, betaH, 1);
    cblas_dscal(ng, 1.0/xsum, gammaH, 1);
    cblas_dscal(nk, 1.0/xsum, kappaH, 1);
    cblas_dscal(ns, 1.0/xsum, sigmaH, 1);
    cblas_dscal(nt, 1.0/xsum, thetaH, 1);
    *locHist = locH;
    *magHist = magH;
    *betaHist = betaH;
    *gammaHist = gammaH;
    *kappaHist = kappaH;
    *sigmaHist = sigmaH;
    *thetaHist = thetaH;
    return 0;
}
//============================================================================//
static void getBaseAndExp(const double val, double *base, int *exp)
{
    char cval[64], cexp[64], temp[64];
    int i, k, l;
    bool lp1;
    memset(cval, 0, 64*sizeof(char));
    memset(cexp, 0, 64*sizeof(char));
    memset(temp, 0, 64*sizeof(char));
    sprintf(temp, "%04.2e", val);
    lp1 = true;
    k = 0;
    l = 0;
    for (i=0; i<strlen(temp); i++)
    {
        if (temp[i] == 'e')
        {
            lp1 = false;
            continue;
        }
        else
        {
            if (lp1)
            {
                cval[k] = temp[i];
                k = k + 1;
            }
            else
            {
                cexp[l] = temp[i];
                l = l + 1;
            }
        }   
    }
    *exp = atoi(cexp);
    *base = atof(cval);
    return;
}

static void getBaseAndExpMT(const double *mtIn, double *mtOut, int *expOut)
{
    double xfact;
    int expWork, i;
    getBaseAndExp(mtIn[0], &mtOut[0], expOut);
    for (i=1; i<6; i++)
    {
        getBaseAndExp(mtIn[i], &mtOut[i], &expWork);
        if (expWork > *expOut){*expOut = expWork;}
    }
    // Rescale
    for (i=0; i<6; i++)
    {
        getBaseAndExp(mtIn[i], &mtOut[i], &expWork);
        xfact = pow(10.0, expWork - *expOut);
        mtOut[i] = mtOut[i]*xfact;
    }
    return;
}

int postmt_gmtHelper_makePsmecaLine(const enum compearthCoordSystem_enum basis,
                                    const double *mt,
                                    const double evla, const double evlo,
                                    const double evdp, const char *evid,
                                    char line[128])
{
    const char *fcnm = "postmt_gmtHelper_makePsmecaLine\0";
    double mtUSE[6], mtGMT[6];
    int exp, ierr;
    memset(line, 0, 128*sizeof(char));
    ierr = compearth_convertMT(basis, USE, mt, mtUSE);
    if (ierr != 0)
    {
        printf("%s: Error switching basis\n", fcnm);
        return -1;
    }
    getBaseAndExpMT(mtUSE, mtGMT, &exp);
    sprintf(line, "%f %f %f %f %f %f %f %f %f %d %s\n",
            evlo, evla, evdp,
            mtGMT[0], mtGMT[1], mtGMT[2], mtGMT[3], mtGMT[4], mtGMT[5],
            exp, evid);
    return 0;
}

static void setFillColor(const int i, const int iopt, char color[32])
{
    memset(color, 0, 32*sizeof(char));
    if (i == iopt)
    {
        strcpy(color, "-Gyellow\0");
    }
    else
    {
        strcpy(color, "-Gred\0");
    }
    return;
}
//============================================================================//
/*
    sprintf("%7.2f %5.2f %f \n"< 
            evlo, evla, evdp,
            mrr, mtt, mpp,  
            
Columns: lon lat depth mrr mtt mpp mrt mrp mtp iexp name
-176.96 -29.25 48 7.68 0.09 -7.77 1.39 4.52 -3.26 26 X Y 010176A  
*/
