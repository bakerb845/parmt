#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>
#include <iniparser.h>
#include "prepmt/prepmt_hpulse96.h"
#include "prepmt/prepmt_event.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "sacio.h"
#include "cps.h"
#include "cps_mpi.h"
#include "cps_utils.h"
#include "cps_defaults.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/string/string.h"

#define PROGRAM_NAME "hspec96"

static void grns2sac(const struct hwave_greens_struct grns,
                     struct sacData_struct sac[10]);
struct precompute_struct
{
    char sourceModel[PATH_MAX]; /*!< Source model */
    char crustDir[PATH_MAX];  /*!< Crust1.0 directory */
    char hspecDistanceFile[PATH_MAX];  /*!< Table of distances for hspec */
    char hspecArchiveFile[PATH_MAX];  /*!< Name of archive file with ff green's
                                           functions. */
    double *depths; /*!< Depths (km) at which to compute greens fns */
    double *dists;  /*!< Distances (degrees) at which to compute
                         surface wave greens fns */
    double slat;    /*!< Source latitude (degrees) */
    double slon;    /*!< Source longitude (degrees) */
    double dt;   /*!< Sampling period */
    int npts;    /*!< Number of samples */
    int ndepths; /*!< Number of depths at which to compute greens fns */
    int ndists;  /*!< Number of distances to compute surface waves */
    bool luseSourceMod; /*!< If true then use the local source model */ 
    bool luseCrust; /*!< If true then use crust1.0 at teleseismic distances */
    bool luseDistanceTable; /*!< If true then use a distance table.
                                 Otherwise linearly interpolate distances */
};
int prepmt_hspec96_readParameters(const char *iniFile,
                                  const char *section,
                                  struct precompute_struct *precompute);
int prepmt_hspec96_initializeArchive(
    const char *archiveName, const char *model,
    const int ndepths, const double *__restrict__ depths,
    const int ngcarcs,  const double *__restrict__ gcarcs);
int prepmt_hspec96_readHspec96Parameters(
    const char *iniFile,
    struct hprep96_parms_struct *hprepParms,
    struct hspec96_parms_struct *hspecParms,
    struct hpulse96_parms_struct *hpulseParms);
static void printUsage(void);
static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX]);

/*!
 * @brief Makes the regional fullwaveform fundamental fault Green's 
 *        functions archive.
 *
 */
int main(int argc, char **argv)
{
    char fname[PATH_MAX], iniFile[PATH_MAX], groupName[512], *modelName;
    struct prepmtEventParms_struct event;
    struct sacData_struct *sac;
    struct precompute_struct precompute;
    struct hwave_greens_struct *grns;
    struct hpulse96_data_struct *zresp;
    struct hprep96_parms_struct hprepParms;
    struct hspec96_parms_struct hspecParms;
    struct hpulse96_parms_struct hpulseParms;
    struct vmodel_struct recmod, vmodel;
    double *depthr, *depths, *r, *tshift, *vred, dt;
    hid_t groupID, h5fl;
    int idep, idist, ierr, k, l, myid, ndeps, ndists, nprocs, npts, nwaves;
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    const char *cfaults[10] = {"ZDS\0", "ZSS\0", "ZDD", "ZEX\0",
                               "RDS\0", "RSS\0", "RDD", "REX\0",
                               "TDS\0", "TSS\0"};
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
    const int master = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    sac = NULL;
    grns = NULL;
    zresp = NULL;
    modelName = NULL;
    r = NULL;
    depthr = NULL;
    depths = NULL;
    vred = NULL;
    tshift = NULL; 
    memset(&hprepParms, 0, sizeof(struct hprep96_parms_struct));
    memset(&hspecParms, 0, sizeof(struct hspec96_parms_struct));
    memset(&hpulseParms, 0, sizeof(struct hpulse96_parms_struct));
    memset(&vmodel, 0, sizeof(struct vmodel_struct));
    if (myid == master)
    {
        // Parse the input arguments
        ierr = parseArguments(argc, argv, iniFile);
        if (ierr != 0)
        {
            if (ierr !=-1)
            {
                printf("%s: Failed to parse input commands\n", PROGRAM_NAME);
            }
            goto FINISH_EARLY;
        } 
        // Read the CPS variables
        ierr = prepmt_hspec96_readHspec96Parameters(iniFile,
                                                    &hprepParms, &hspecParms,
                                                    &hpulseParms);
        if (ierr != 0)
        {
            printf("Error reading hprep parameters\n");
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        // Read the additional contrived forward modeling variables
        ierr = prepmt_hspec96_readParameters(iniFile,
                                             "precompute", &precompute);
        if (ierr != 0)
        {
            printf("Error reading precompute parameters\n");
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        ndists = precompute.ndists;
        ndeps = precompute.ndepths;
        // expand to make a grid
        r = (double *) calloc((size_t) (ndists*ndeps), sizeof(double));
        depths = (double *) calloc((size_t) (ndists*ndeps), sizeof(double));
        for (idep=0; idep<ndeps; idep++)
        {
            for (idist=0; idist<ndists; idist++)
            {
                k = idist*ndeps + idep;
                r[k] = precompute.dists[idist]*111.195;
                depths[k] = precompute.depths[idep];
                printf("%f %f\n", r[k]/111.195, depths[k]);
            }
        }
        npts = precompute.npts;
        dt = precompute.dt;
//ndists = 31;
//ndeps = 26;
//npts = 1024; //2048;
//dt = 1.0; //0.5;
        printf("%d %f %s\n", npts, dt, hprepParms.mfile);
        // load the model into memory
        if (precompute.luseSourceMod)
        {
            printf("%s: Loading source model %s...\n",
                   PROGRAM_NAME, precompute.sourceModel);
            strcpy(hprepParms.mfile, precompute.sourceModel);
            ierr = cps_getmod(hprepParms.mfile, &vmodel);
            if (ierr != 0)
            {
                printf("Failed to load source model\n");
                goto FINISH_EARLY;
            }
        }
        else if (precompute.luseCrust)
        {
            ierr = prepmt_event_initializeFromIniFile(iniFile, &event);
            if (ierr != 0)
            {
                printf("%s: Error getting event info\n", PROGRAM_NAME);
                goto FINISH_EARLY;
            } 
            printf("%s: Loading CPS model (lat,lon)=(%f,%f)\n", PROGRAM_NAME,
                   event.latitude, event.longitude);
            ierr = cps_crust1_getCrust1ForHerrmann(precompute.crustDir, false,
                                              event.latitude, event.longitude,
                                              &event.latitude, &event.longitude,
                                              1, &vmodel, &recmod);
            cps_utils_freeVmodelStruct(&recmod);
        }
        else
        {
            printf("%s: Model is not set\n", PROGRAM_NAME);
            ierr = 1;
            goto FINISH_EARLY;
        }
        // Initialize the archive - remove NULL terminator
        modelName = string_strip(NULL, vmodel.title);
        printf("%s: Initializing archive file %s with model %s\n",
               PROGRAM_NAME, precompute.hspecArchiveFile, modelName); 
        ierr = prepmt_hspec96_initializeArchive(precompute.hspecArchiveFile,
                                                modelName,
                                                ndeps, depths,
                                                ndists, precompute.dists);
        if (ierr != 0)
        {
            printf("%s: Failed to initialize archive\n", PROGRAM_NAME);
            goto FINISH_EARLY;
        }
        free(precompute.depths);
        free(precompute.dists);
    }
FINISH_EARLY:;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ierr, 1, MPI_INTEGER, master, MPI_COMM_WORLD);
    if (ierr != 0){goto FINISH;}
    cps_broadcast_vmodelStruct(MPI_COMM_WORLD, master, &vmodel);
    cps_broadcast_hprep96ParmsStruct(MPI_COMM_WORLD, master, &hprepParms);
    cps_broadcast_hspec96ParmsStruct(MPI_COMM_WORLD, master, &hspecParms);
    cps_broadcast_hpulse96ParmsStruct(MPI_COMM_WORLD, master, &hpulseParms);
    MPI_Bcast(&ndists, 1, MPI_INTEGER, master, MPI_COMM_WORLD);
    MPI_Bcast(&ndeps,  1, MPI_INTEGER, master, MPI_COMM_WORLD);
    MPI_Bcast(&npts, 1, MPI_INTEGER, master, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
    // initialize receiver information
    nwaves = ndists*ndeps;
    depthr = (double *) calloc((size_t) nwaves, sizeof(double));
    if (myid != master)
    {
        depths = (double *) calloc((size_t) nwaves, sizeof(double));
        r = (double *) calloc((size_t) nwaves, sizeof(double));
    }
    tshift = (double *) calloc((size_t) nwaves, sizeof(double));
    vred = (double *) calloc((size_t) nwaves, sizeof(double));
    MPI_Bcast(depths, nwaves, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(r, nwaves, MPI_DOUBLE, master, MPI_COMM_WORLD);
    // set the receiver distances
//r[0] = 55.597; //111.195;
//r[1] = 111.195;
//r[2] = 166.792;
//depths[0] = 7.0;
//for (idep=0; idep<ndeps; idep++)
//{
//    for (idist=0;  idist<ndists; idist++)
//    {
//        k = idist*ndeps + idep;
//        r[k] = 55.597*(double) (idist + 1); 
//        depths[k] = (double) idep;
//    }
//}
/*
printf("premature end\n");
MPI_Finalize();
return 0;
*/
    // call hspec
    zresp = hspec96_mpi_interface(MPI_COMM_WORLD, nwaves, npts, dt,
                                  r, tshift, vred, depthr, depths,
                                  hprepParms, hspecParms,
                                  vmodel, &ierr);
    MPI_Barrier(MPI_COMM_WORLD);
    if (ierr != 0)
    {
        printf("Error calling hspec\n");
        goto FINISH;
    }
    // convolve STF and convert back to time domain
    grns = hpulse96_mpi_interface(MPI_COMM_WORLD, nwaves, hpulseParms,
                                  zresp, &ierr); 
    if (ierr != 0)
    {
        printf("Error computing greens fns\n");
        goto FINISH;
    }
    // dump the results
    if (myid == master)
    {
        printf("Writing results...\n");
        h5fl = H5Fopen(precompute.hspecArchiveFile, H5F_ACC_RDWR, H5P_DEFAULT);
        sac = (struct sacData_struct *)
              calloc(10, sizeof(struct sacData_struct));
        k = 0;
        //for (idist=0; idist<ndists; idist++)
        for (idep=0; idep<ndeps; idep++)
        {
            memset(fname, 0, PATH_MAX*sizeof(char));
            sprintf(fname, "dump/depth_%d", idep);
            os_makedirs(fname);

            memset(groupName, 0, 512*sizeof(char));
            sprintf(groupName, "/%s/Depth_%d", modelName, idep);
            groupID = H5Gopen(h5fl, groupName, H5P_DEFAULT);
            //for (idep=0; idep<ndeps; idep++)
            for (idist=0; idist<ndists; idist++)
            {
                k = idist*ndeps + idep;
                grns2sac(grns[k], sac);
                memset(fname, 0, PATH_MAX*sizeof(char));
                sprintf(fname, "dump/depth_%d/%04d01ZDD.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[0]);
                sprintf(fname, "dump/depth_%d/%04d02RDD.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[1]);
                sprintf(fname, "dump/depth_%d/%04d03ZDS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[2]);
                sprintf(fname, "dump/depth_%d/%04d04RDS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[3]);
                sprintf(fname, "dump/depth_%d/%04d05TDS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[4]);
                sprintf(fname, "dump/depth_%d/%04d06ZSS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[5]);
                sprintf(fname, "dump/depth_%d/%04d07RSS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[6]);
                sprintf(fname, "dump/depth_%d/%04d08TSS.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[7]);
                sprintf(fname, "dump/depth_%d/%04d09ZEX.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[8]);
                sprintf(fname, "dump/depth_%d/%04d10REX.SAC", idep, idist+1);
                sacio_writeTimeSeriesFile(fname, sac[9]);
                // Dump to the H5 archive and release memory
                for (l=0; l<10; l++)
                {
                    // Rescale to unit magnitude
                    cblas_dscal(sac[l].npts, xscal, sac[l].data, 1);
                    sprintf(fname, "Greens_%s_%d", cfaults[l], idist);
                    sacioh5_writeTimeSeries2(fname, groupID, sac[l]);
                    sacio_free(&sac[l]);
                }
//break;
            }
            H5Gclose(groupID);
//break;
        }
        if (sac != NULL){free(sac);}
        H5Fclose(h5fl);
        printf("%s: Finishing program\n", PROGRAM_NAME);
    }
    MPI_Barrier(MPI_COMM_WORLD);
FINISH:;
    if (grns != NULL)
    {
        for (k=0; k<nwaves; k++)
        {
            cps_utils_freeHwaveGreensStruct(&grns[k]);
        }
        free(grns);
        grns = NULL;
    }
    if (zresp != NULL)
    {
        for (k=0; k<nwaves; k++)
        {
            cps_utils_freeHpulse96DataStruct(&zresp[k]);
        }
        free(zresp);
        zresp = NULL;
    }
    // release memory
    memory_free8c(&modelName);
    free(depthr);
    free(depths);
    free(r);
    free(tshift);
    free(vred); 
    cps_utils_freeHpulse96ParmsStruct(&hpulseParms);
    cps_utils_freeVmodelStruct(&vmodel);
    MPI_Finalize();
    return 0;
}
//============================================================================//
int prepmt_hspec96_initializeArchive(
    const char *archiveName, const char *model,
    const int ndepths, const double *__restrict__ depths,
    const int ngcarcs,  const double *__restrict__ gcarcs)
{
    char modelName[512], depthName[512]; //, varName[512];
    hid_t h5fl, dataSet, dataSpace, depthGroup, modelGroup;
    hsize_t dims[1];
    int idep;
    //------------------------------------------------------------------------//
    h5fl = H5Fcreate(archiveName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    memset(modelName, 0, 512*sizeof(char));
    strcpy(modelName, "/\0");
    strcat(modelName, model);
    modelGroup = H5Gcreate2(h5fl, modelName, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
    // write the depths and distances
    dims[0] = (hsize_t) ndepths;
    dataSpace = H5Screate_simple(1, dims, NULL);
    dataSet = H5Dcreate2(modelGroup, "Depths\0", H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
             H5S_ALL, H5P_DEFAULT, depths);
    H5Dclose(dataSet);
    H5Sclose(dataSpace);
    // Write the distances
    dims[0] = (hsize_t) ngcarcs;
    dataSpace = H5Screate_simple(1, dims, NULL);
    dataSet = H5Dcreate2(modelGroup, "Distances\0", H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
             H5S_ALL, H5P_DEFAULT, gcarcs);
    H5Dclose(dataSet);
    H5Sclose(dataSpace);
    // Create the depths groups
    for (idep=0; idep<ndepths; idep++)
    {
        memset(depthName, 0, 512*sizeof(char));
        sprintf(depthName, "%s/Depth_%d", modelName, idep);
        depthGroup = H5Gcreate2(h5fl, depthName, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(depthGroup);
    }
    H5Gclose(modelGroup);
    H5Fclose(h5fl);
    return 0;
}
//============================================================================//
int prepmt_hspec96_readParameters(const char *iniFile,
                                  const char *section,
                                  struct precompute_struct *precompute)
{
    const char *fcnm = "prepmt_hspec96_readParameters\0";
    FILE *dfl;
    const char *s;
    char *dirName, cline[256], vname[128];
    double depth0, depth1, dist0, dist1;
    int i, ierr;
    dictionary *ini;
    memset(precompute, 0, sizeof(struct precompute_struct));
    if (!os_path_isfile(iniFile))
    {
        printf("%s: ini file: %s does not exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        printf("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:ndepths", section);
    precompute->ndepths = iniparser_getint(ini, vname, 1);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:depthMin", section);
    depth0 = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:depthMax", section);
    depth1 = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:ndists", section);
    precompute->ndists = iniparser_getint(ini, vname, 1);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:minDistance", section);
    dist0 = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:maxDistance", section);
    dist1 = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:dtHspec", section);
    precompute->dt = iniparser_getdouble(ini, vname, 1.0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:nptsHspec", section);
    precompute->npts = iniparser_getint(ini, vname, 1024);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:luseCrust1", section);
    precompute->luseCrust = iniparser_getboolean(ini, vname, 0);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:luseSourceModel", section);
    precompute->luseSourceMod = iniparser_getboolean(ini, vname, 1);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:sourceModel", section);
    s = iniparser_getstring(ini, vname, NULL);
    if (s != NULL)
    {
        strcpy(precompute->sourceModel, s);
        if (!os_path_isfile(precompute->sourceModel) &&
            precompute->luseSourceMod)
        {
            printf("%s: Source model doesn't %s exist\n", fcnm, s);
            precompute->luseSourceMod = false;
        }
    }
    else
    {
        precompute->luseSourceMod = false;
    }

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:hspecArchiveFile", section);
    s = iniparser_getstring(ini, vname, "hspecFFgreens.h5");
    dirName = os_dirname(s);
    if (!os_path_isdir(dirName))
    {
        ierr = os_makedirs(dirName);
        if (ierr != 0)
        {
            printf("%s: Failed to make output directory: %s\n", fcnm, dirName);
            return -1;
        }
    }
    free(dirName);
    strcpy(precompute->hspecArchiveFile, s);

    //memset(vname, 0, 128*sizeof(char));
    //sprintf(vname, "%s:modelIdentifier", section);
    //s = iniparser_getstring(ini, vname, "UNKNOWN");
    //strcpy(precompute->modelIdentifer, s);

    memset(vname, 0, 128*sizeof(char));
    sprintf(vname, "%s:crustDir", section);
    s = iniparser_getstring(ini, vname, CPS_DEFAULT_CRUST1_DIRECTORY);
    if (s != NULL)
    {
        strcpy(precompute->crustDir, s);
        if (!os_path_isdir(precompute->crustDir) && precompute->luseCrust)
        {
            printf("%s: crust1.0 directory %s doesn't exist\n", fcnm, s);
            precompute->luseCrust = false;
        }
    }
    else
    {
        precompute->luseCrust = false;
    }

    memset(vname, 0, 128*sizeof(char)); 
    sprintf(vname, "%s:luseDistanceTable", section);
    precompute->luseDistanceTable = iniparser_getboolean(ini, vname, false);
    if (precompute->luseDistanceTable)
    {
        memset(vname, 0, 128*sizeof(char));
        sprintf(vname, "%s:hspecDistanceFile", section);
        s = iniparser_getstring(ini, vname, NULL);
        if (!os_path_isfile(s))
        {
            printf("%s: files %s doesn't exist\n", fcnm, s);
            return -1;
        }
        strcpy(precompute->hspecDistanceFile, s);
        dfl = fopen(s, "r");
        memset(cline, 0, 256*sizeof(char));
        fgets(cline, 256, dfl);
        memset(cline, 0, 256*sizeof(char));
        fgets(cline, 256, dfl);
        precompute->ndists = atoi(cline);
        if (precompute->ndists < 1)
        {
            printf("%s: Invalid number of distances\n", fcnm);
            return -1;
        }
        precompute->dists
            = (double *) calloc((size_t) precompute->ndists, sizeof(double));
        for (i=0; i<precompute->ndists; i++)
        {
            memset(cline, 0, 256*sizeof(char));
            if (fgets(cline, 256, dfl) == NULL)
            {
                printf("%s: premature end of distance file\n", fcnm);
                return -1;
            }
            precompute->dists[i] = atof(cline);
        }
        fclose(dfl);
    }
    else
    {
        precompute->dists
            = (double *) calloc((size_t) precompute->ndists, sizeof(double));
        ierr = array_linspace64f_work(dist0, dist1,
                                      precompute->ndists, precompute->dists);
    }
    // depths
    precompute->depths
        = (double *) calloc((size_t) precompute->ndepths, sizeof(double));
    ierr = array_linspace64f_work(depth0, depth1,
                                  precompute->ndepths, precompute->depths);
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
int prepmt_hspec96_readHspec96Parameters(
    const char *iniFile,
    struct hprep96_parms_struct *hprepParms,
    struct hspec96_parms_struct *hspecParms,
    struct hpulse96_parms_struct *hpulseParms)
{
    const char *fcnm = "prepmt_hspec96_readHspec96Parameters\0";
    int ierr;
    cps_setHprep96Defaults(hprepParms);
    cps_setHspec96Defaults(hspecParms);
    ierr = prepmt_hpulse96_readHpulse96Parameters(iniFile,
                                                  "hpulse96\0", hpulseParms);
    if (ierr != 0)
    {
        printf("%s: Error reading hpulse parameters\n", fcnm);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Parses the command line arguments for the ini file.
 *
 * @param[in] argc      Number of arguments input to the program.
 * @param[in] argv      Command line arguments.
 *
 * @param[out] iniFile  If the result is 0 then this is the ini file.
 *
 * @result 0 indicates success. \n
 *        -1 indicates the user inquired about the usage. \n
 *        -2 indicates the command line argument is invalid.
 *
 * @author Ben Baker, ISTI
 *
 */
static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX])
{
    bool linFile;
    linFile = false;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    while (true)
    {   
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"ini_file", required_argument, 0, 'i'},
            {"section", required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            linFile = true;
        }
        else if (c == 'h' || c == '?')
        {
            printUsage();
            return -2;
        }
        else
        {
            printf("%s: Unknown options: %s\n",
                   PROGRAM_NAME, argv[optionIndex]);
        }
    }
    if (!linFile)
    {
        printf("%s: Error must specify ini file\n\n", PROGRAM_NAME);
        printUsage();
        return -1;
    }
    else
    {
        if (!os_path_isfile(iniFile))
        {
            printf("%s: Error - ini file: %s does not exist\n",
                   PROGRAM_NAME, iniFile);
            return EXIT_FAILURE;
        }
    }
    return 0;
}
//============================================================================//
static void printUsage(void)
{
    printf("Usage:\n   %s -i iniFile\n", PROGRAM_NAME); 
    printf("or\n");
    printf("    mpirun -np nProcessors %s -i iniFile\n", PROGRAM_NAME);
    printf("Required arguments:\n");
    printf("    -i ini_file specifies the initialization file\n");
    printf("Optional arguments:\n");
    printf("    -h displays this message\n");
    printf("    -np Number of computer processors to use\n");
    return;
}
//============================================================================//
static void grns2sac(const struct hwave_greens_struct grns,
                     struct sacData_struct sac[10])
{
    const char *kcmpnm[10] = {"ZDD\0", "RDD\0", "ZDS\0", "RDS\0", "TDS\0",
                              "ZSS\0", "RSS\0", "TSS\0", "ZEX\0", "REX\0"};
    double cmpinc, cmpaz;
    int i, k;
    bool lsh;
    for (k=0; k<10; k++)
    {
        memset(&sac[k], 0, sizeof(struct sacData_struct));
        sacio_setDefaultHeader(&sac[k].header);
        lsh = false;
        sac[k].header.npts = grns.npts;
        //sac[k].header.b = grns.b;
        sac[k].npts = grns.npts;
        sac[k].header.b = grns.t0;
        sac[k].header.e = grns.t0 + (double) (grns.npts - 1)*grns.dt;
        sac[k].header.delta = grns.dt;
        sac[k].header.dist = grns.dist;
        sac[k].header.evdp = grns.evdep;
        sac[k].header.stel = grns.stelel;
        sac[k].header.gcarc = grns.dist/111.195;
        sac[k].header.o =-grns.t0;
        strcpy(sac[k].header.kevnm, "SYNTHETIC");
        strcpy(sac[k].header.ko, "O\0");
        sac[k].header.a = grns.timep;
        strcpy(sac[k].header.ka, "P\0");
        sac[k].header.az = 0.0;
        sac[k].header.baz = 180.0;
        sac[k].header.lcalda = 0;
        sac[k].header.norid = 0;
        sac[k].header.nevid = 0;
        sac[k].header.nzyear = 1970;
        sac[k].header.nzjday = 1;
        sac[k].header.nzhour = 0;
        sac[k].header.nzmin = 0;
        sac[k].header.nzsec = 0;
        sac[k].header.nzmsec = 0;
        sac[k].header.lhaveHeader = true;
        strcpy(sac[k].header.kcmpnm, kcmpnm[k]); 
        sac[k].data = (double *) calloc((size_t) grns.npts, sizeof(double));
        cmpaz = 0.0;
        cmpinc = 0.0;
        if (k == 0)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.zdd[i];}
        }
        else if (k == 1)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.rdd[i];}
            cmpinc = 90.0;
        }
        else if (k == 2)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.zds[i];}
        }
        else if (k == 3)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.rds[i];}
            cmpinc = 90.0;
        }
        else if (k == 4)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.tds[i];}
            lsh = true;
            cmpaz = 90.0;
            cmpinc = 90.0;
        }
        else if (k == 5)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.zss[i];}
        }
        else if (k == 6)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.rss[i];}
            cmpinc = 90.0;
        }
        else if (k == 7)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.tss[i];}
            lsh = true;
            cmpaz = 90.0;
            cmpinc = 90.0;
        }
        else if (k == 8)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.zex[i];}
        }
        else if (k == 9)
        {
            for (i=0; i<grns.npts; i++){sac[k].data[i] = grns.rex[i];}
            cmpinc = 90.0;
        }
        if (lsh)
        {
            sac[k].header.t1 = grns.timesh;
            strcpy(sac[k].header.kt0, "SH\0");
        }
        else
        {
            sac[k].header.t1 = grns.timesv;
            strcpy(sac[k].header.kt0, "SV\0");
        }
        sac[k].header.cmpaz = cmpaz;
        sac[k].header.cmpinc = cmpinc;
    }
    return;
}
