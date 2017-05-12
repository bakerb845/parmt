#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mkl_cblas.h>
#include "sacio.h"
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"

int main()
{
    char *ffList[10] = {"rotation_data/B00103ZDS.sac\0",
                        "rotation_data/B00106ZSS.sac\0",
                        "rotation_data/B00101ZDD.sac\0",
                        "rotation_data/B00109ZEX.sac\0",
                        "rotation_data/B00104RDS.sac\0",
                        "rotation_data/B00107RSS.sac\0",
                        "rotation_data/B00102RDD.sac\0",
                        "rotation_data/B00110REX.sac\0",
                        "rotation_data/B00105TDS.sac\0",
                        "rotation_data/B00108TSS.sac\0"};
    char *vList[6] = {"rotation_data/xx101.hhz.sac\0",
                      "rotation_data/yy101.hhz.sac\0",
                      "rotation_data/zz101.hhz.sac\0",
                      "rotation_data/xy101.hhz.sac\0",
                      "rotation_data/xz101.hhz.sac\0",
                      "rotation_data/yz101.hhz.sac\0"};
    char *nList[6] = {"rotation_data/xx101.hhn.sac\0",
                      "rotation_data/yy101.hhn.sac\0",
                      "rotation_data/zz101.hhn.sac\0",
                      "rotation_data/xy101.hhn.sac\0",
                      "rotation_data/xz101.hhn.sac\0",
                      "rotation_data/yz101.hhn.sac\0"};
    char *eList[6] = {"rotation_data/xx101.hhe.sac\0",
                      "rotation_data/yy101.hhe.sac\0",
                      "rotation_data/zz101.hhe.sac\0",
                      "rotation_data/xy101.hhe.sac\0",
                      "rotation_data/xz101.hhe.sac\0",
                      "rotation_data/yz101.hhe.sac\0"};
    struct sacData_struct sacZDS, sacZSS, sacZDD, sacZEX,
                          sacRDS, sacRSS, sacRDD, sacREX,
                          sacTDS, sacTSS, sac;
    double *G, baz, cmpaz, cmpinc, xnorm;
    int ierr, imt, indx, npgrns;
    const double az = 25.0;
    const double M0grns = 1.e20; //3.5481338923357603e+23;
    ierr = 0;
    // Read the green's functions generated by CPS
    memset(&sacZDS, 0, sizeof(struct sacData_struct));
    memset(&sacZSS, 0, sizeof(struct sacData_struct));
    memset(&sacZDD, 0, sizeof(struct sacData_struct));
    memset(&sacZEX, 0, sizeof(struct sacData_struct));
    memset(&sacRDS, 0, sizeof(struct sacData_struct));
    memset(&sacRSS, 0, sizeof(struct sacData_struct));
    memset(&sacRDD, 0, sizeof(struct sacData_struct));
    memset(&sacREX, 0, sizeof(struct sacData_struct));
    memset(&sacTDS, 0, sizeof(struct sacData_struct));
    memset(&sacTSS, 0, sizeof(struct sacData_struct));
    memset(&sac,    0, sizeof(struct sacData_struct));
    ierr += sacio_readTimeSeriesFile(ffList[0], &sacZDS);
    ierr += sacio_readTimeSeriesFile(ffList[1], &sacZSS);
    ierr += sacio_readTimeSeriesFile(ffList[2], &sacZDD);
    ierr += sacio_readTimeSeriesFile(ffList[3], &sacZEX);
    ierr += sacio_readTimeSeriesFile(ffList[4], &sacRDS);
    ierr += sacio_readTimeSeriesFile(ffList[5], &sacRSS);
    ierr += sacio_readTimeSeriesFile(ffList[6], &sacRDD);
    ierr += sacio_readTimeSeriesFile(ffList[7], &sacREX);
    ierr += sacio_readTimeSeriesFile(ffList[8], &sacTDS);
    ierr += sacio_readTimeSeriesFile(ffList[9], &sacTSS);
    if (ierr != 0)
    {
        printf("Failed to load fundmaental faults\n");
        return EXIT_FAILURE;
    }
    // Set space
    npgrns = sacZDS.npts;
    G = memory_calloc64f(6*npgrns);
    // Test the vertical
    cmpaz = 0.0;
    cmpinc =-90.0;
    baz = fmod(az + 180.0, 360.0); 
    ierr = parmt_utils_ff2mtGreens64f(npgrns, 1, az, baz,
                                      cmpaz, cmpinc,
                                      sacZDS.data, sacZSS.data,
                                      sacZDD.data, sacZEX.data,
                                      sacRDS.data, sacRSS.data,
                                      sacRDD.data, sacREX.data,
                                      sacTDS.data, sacTSS.data,
                                      &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                                      &G[3*npgrns], &G[4*npgrns], &G[5*npgrns]);
    if (ierr != 0)
    {
        printf("Failed to called vertical\n");
        return EXIT_FAILURE;
    }
    for (imt=0; imt<6; imt++)
    {
        ierr = sacio_readTimeSeriesFile(vList[imt], &sac);
        if (ierr != 0)
        {
            printf("Failed to read %s\n", vList[imt]);
            return EXIT_FAILURE;
        }
        if (sac.npts != npgrns)
        {
            printf("Veritcal size mismatch\n");
            return EXIT_FAILURE;
        }
        cblas_dscal(sac.npts, M0grns, sac.data, 1);
        indx = imt*npgrns;
        xnorm = array_normDiff64f(npgrns, sac.data, &G[indx],
                                  ONE_NORM, 1.0, &ierr);
        if (fabs(xnorm) > 1.e-10)
        {
            printf("Failed vertical test\n");
            return EXIT_FAILURE;
        }
        sacio_free(&sac);
    }
    // Test the north component
    cmpaz = 0.0;
    cmpinc = 90.0;
    baz = fmod(az + 180.0, 360.0);
    ierr = parmt_utils_ff2mtGreens64f(npgrns, 2, az, baz,
                                      cmpaz, cmpinc,
                                      sacZDS.data, sacZSS.data,
                                      sacZDD.data, sacZEX.data,
                                      sacRDS.data, sacRSS.data,
                                      sacRDD.data, sacREX.data,
                                      sacTDS.data, sacTSS.data,
                                      &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                                      &G[3*npgrns], &G[4*npgrns], &G[5*npgrns]);
    if (ierr != 0)
    {
        printf("Failed to called north\n");
        return EXIT_FAILURE;
    }
    for (imt=0; imt<6; imt++)
    {
        ierr = sacio_readTimeSeriesFile(nList[imt], &sac);
        if (ierr != 0)
        {
            printf("Failed to read %s\n", nList[imt]);
            return EXIT_FAILURE;
        }
        if (sac.npts != npgrns)
        {
            printf("north size mismatch\n");
            return EXIT_FAILURE;
        }
        cblas_dscal(sac.npts, M0grns, sac.data, 1); 
        indx = imt*npgrns;
        xnorm = array_normDiff64f(npgrns, sac.data, &G[indx],
                                  ONE_NORM, 1.0, &ierr);
        if (fabs(xnorm) > 1.e-10)
        {
            printf("Failed north test\n");
            return EXIT_FAILURE;
        }
        sacio_free(&sac);
    }
    // Test the east component
    cmpaz = 90.0;
    cmpinc = 0.0;
    baz = fmod(az + 180.0, 360.0);
    ierr = parmt_utils_ff2mtGreens64f(npgrns, 3, az, baz,
                                      cmpaz, cmpinc, 
                                      sacZDS.data, sacZSS.data,
                                      sacZDD.data, sacZEX.data,
                                      sacRDS.data, sacRSS.data,
                                      sacRDD.data, sacREX.data,
                                      sacTDS.data, sacTSS.data,
                                      &G[0*npgrns], &G[1*npgrns], &G[2*npgrns],
                                      &G[3*npgrns], &G[4*npgrns], &G[5*npgrns]);
    if (ierr != 0)
    {
        printf("Failed to called east\n");
        return EXIT_FAILURE;
    }
    for (imt=0; imt<6; imt++)
    {
        ierr = sacio_readTimeSeriesFile(eList[imt], &sac);
        if (ierr != 0)
        {
            printf("Failed to read %s\n", eList[imt]);
            return EXIT_FAILURE;
        }
        if (sac.npts != npgrns)
        {
            printf("east size mismatch\n");
            return EXIT_FAILURE;
        }
        cblas_dscal(sac.npts, M0grns, sac.data, 1);
        indx = imt*npgrns;
        xnorm = array_normDiff64f(npgrns, sac.data, &G[indx],
                                  ONE_NORM, 1.0, &ierr);
        if (fabs(xnorm) > 1.e-10)
        {
            printf("Failed east test\n");
            return EXIT_FAILURE;
        }
/*
if (imt == 1)
{
FILE *test = fopen("test.txt","w");
for (i=0; i<npgrns; i++)
{
fprintf(test, "%d %e %e\n", i, G[indx+i], sac.data[i]);
}
fclose(test);
return 0;
}
*/
        sacio_free(&sac);
    }
    // Free space
    memory_free64f(&G);
    sacio_free(&sacZDS);
    sacio_free(&sacZSS);
    sacio_free(&sacZDD);
    sacio_free(&sacZEX);
    sacio_free(&sacRDS);
    sacio_free(&sacRSS);
    sacio_free(&sacRDD);
    sacio_free(&sacREX);
    sacio_free(&sacTDS);
    sacio_free(&sacTSS);
    return EXIT_SUCCESS;
}
