#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "parmt_utils.h"
#include "sacio.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

#define MAXOBS 150000
#define LATITUDES "Latitudes"
#define LONGITUDES "Longitudes"
#define DEPTHS "Depths"

struct h5_dataStruct
{
    char network[64];       /*!< network name */
    char station[64];       /*!< Station name */
    char channel[64];       /*!< Channel name */
    char locationCode[64];  /*!< Location code */
    char phaseID[64];       /*!< Phase identifier - e.g. P, Rayleigh, S */ 
    hvl_t data;             /*!< Data */
    double dt;              /*!< Sampling period (seconds) */
    double weight;          /*!< Relative weight [0,1] */
    double latitude;        /*!< Station latitude (degrees) */
    double longitude;       /*!< Station longitude (degrees) */
    double elevation;       /*!< Station elevation above sea-level (meters) */
    double pick;            /*!< Seconds relative to trace start to denote
                                 a body wave phase arrival. */
    int npts;               /*!< Number of data points */
    int polarity;           /*!< First motition polarity:
                                 +1 is up.
                                  0 is unknown.
                                 -1 is down. */
    int maxLags;            /*!< Max number of time lags in optimization */
};

struct h5_greensStruct
{
    char network[64];      /*!< Network name */
    char station[64];      /*!< Station name */
    char channel[64];      /*!< Channel name */
    char locationCode[64]; /*!< Location code */
    hvl_t gxx;             /*!< Green's function scaled by NED mxx moment
                                tensor term.  Units are equivalent to data. */
    hvl_t gyy;             /*!< Green's function scaled by NED myy moment
                                tensor term.  Units are equivalent to data. */
    hvl_t gzz;             /*!< Green's funtcion scaled by NED mzz moment
                                tensor term.  Units are equivalent to data. */
    hvl_t gxy;             /*!< Green's function scaled by NED mxy moment
                                tensor term.  Units are equivalent to data. */
    hvl_t gxz;             /*!< Green's function scaled by NED mxz moment
                                tensor term.  Units are equivalent to data. */
    hvl_t gyz;             /*!< Green's function scaled by NED myz moment
                                tensor term */
    double dt;             /*!< Sampling period (seconds) */
    double pick;           /*!< If this is a body wave Green's function then
                                this is the seconds relative to the trace start
                                time to denote a body wave phase arrival. This
                                should match the data. */
    int npts;              /*!< Number of points in Green's functions - this
                                should match the data */
    int locationID;        /*!< Location ID */
};

int utils_dataArchive_setFileName(const char *dirnm, const char *projnm,
                                  char fname[PATH_MAX])
{
    const char *fcnm = "utils_dataArchive_setFileName\0";
    size_t lenos; 
    memset(fname, 0, PATH_MAX*sizeof(char));
    if (projnm == NULL)
    {
        log_errorF("%s: Error project name can't be NULL\n", fcnm);
        return -1;
    }
    if (strlen(projnm) == 0)
    {
        log_errorF("%s: Error project name undefined\n", fcnm);
        return -1;
    }
    if (dirnm == NULL)
    {
        strcpy(fname, "./\0");
    }
    else
    {
        lenos = strlen(fname);
        if (lenos == 0)
        {
            strcpy(fname, "./\0");
        }
        else
        {
            strcpy(fname, dirnm);
            if (fname[lenos-1] != '/'){strcat(fname, "/\0");}
        }
    }
    strcat(fname, projnm);
    strcat(fname, ".h5\0");
    return 0;
}
//============================================================================//
static int utils_archive_createDataStructure(const hid_t groupID)
{
    const char *fcnm = "utils_archive_createDataStructure\0";
    hid_t dataType, string64Type, vlenCData, vlenDData;
    herr_t status;
    //------------------------------------------------------------------------//
    //
    // define the variable length types
    status = 0;
    string64Type = H5Tcopy(H5T_C_S1);  
    H5Tset_size(string64Type, 64);
    vlenCData = H5Tvlen_create(string64Type);
    vlenDData = H5Tvlen_create(H5T_NATIVE_DOUBLE);
    // create the type
    dataType = H5Tcreate(H5T_COMPOUND,
                         sizeof(struct h5_dataStruct));
    // insert the components of the structure
    status += H5Tinsert(dataType, "Network\0",
                        HOFFSET(struct h5_dataStruct, network),
                        vlenCData); 
    status += H5Tinsert(dataType, "Station\0",
                        HOFFSET(struct h5_dataStruct, station),
                        vlenCData);
    status += H5Tinsert(dataType, "Channel\0",
                        HOFFSET(struct h5_dataStruct, channel),
                        vlenCData);
    status += H5Tinsert(dataType, "LocationCode\0",
                        HOFFSET(struct h5_dataStruct, locationCode),
                        vlenCData);
    status += H5Tinsert(dataType, "PhaseIdentifer\0",
                        HOFFSET(struct h5_dataStruct, phaseID),
                        vlenCData);
    status += H5Tinsert(dataType, "Data\0",
                        HOFFSET(struct h5_dataStruct, data),
                        vlenDData);
    status += H5Tinsert(dataType, "SamplingPeriod\0",
                        HOFFSET(struct h5_dataStruct, dt),
                        H5T_NATIVE_DOUBLE);
    status += H5Tinsert(dataType, "DataWeight\0",
                        HOFFSET(struct h5_dataStruct, weight),
                        H5T_NATIVE_DOUBLE);
    status += H5Tinsert(dataType, "StationLatitude\0",
                        HOFFSET(struct h5_dataStruct, latitude),
                        H5T_NATIVE_DOUBLE);
    status += H5Tinsert(dataType, "StationLongitude\0",
                        HOFFSET(struct h5_dataStruct, longitude),
                        H5T_NATIVE_DOUBLE);
    status += H5Tinsert(dataType, "StationElevation\0",
                        HOFFSET(struct h5_dataStruct, elevation),
                        H5T_NATIVE_DOUBLE); 
    status += H5Tinsert(dataType, "PhaseArrivalTime\0",
                        HOFFSET(struct h5_dataStruct, pick),
                        H5T_NATIVE_DOUBLE);
    status += H5Tinsert(dataType, "NumberOfPoints\0",
                        HOFFSET(struct h5_dataStruct, npts),
                        H5T_NATIVE_INT);
    status += H5Tinsert(dataType, "Polarity\0",
                        HOFFSET(struct h5_dataStruct, polarity),
                        H5T_NATIVE_INT);
    status += H5Tinsert(dataType, "MaxNumberOfTimeLags\0",
                        HOFFSET(struct h5_dataStruct, maxLags),
                        H5T_NATIVE_INT);
    if (status != 0)
    {
        log_errorF("%s: Failed to pack type\n", fcnm);
        return -1;
    }
    // Commit it
    status = H5Tcommit2(groupID, "dataStructure\0", dataType,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (status != 0)
    {
        log_errorF("%s: Failed to commit structure\n", fcnm);
        return -1;
    } 
    // Release memory
    status += H5Tclose(dataType); 
    status += H5Tclose(vlenCData);
    status += H5Tclose(vlenDData);
    status += H5Tclose(string64Type); 
    if (status != 0)
    {
        log_errorF("%s: Failed to close types\n", fcnm);
        return -1;
    }
    return status;
}
//============================================================================//
int utils_dataArchive_getLocationID(const hid_t h5fl,
                                    const double evla,
                                    const double evlo,
                                    const double evdp)
{
    double *data, *evdps, *evlas, *evlos;
    hsize_t dims[1];
    int i, iloc, nlocs, rank;
    const char *names[3] = {LATITUDES, LONGITUDES, DEPTHS};
    hid_t dataSet, dataSpace, groupID, memSpace;
    iloc =-1;
    nlocs = 0;
    evdps = NULL;
    evlas = NULL;
    evlos = NULL;
    groupID = H5Gopen(h5fl, "/Locations\0", H5P_DEFAULT);
    for (i=0; i<3; i++)
    {
        dataSet = H5Dopen(groupID, names[i], H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        rank = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
        if (rank != 1)
        {
            printf("Invalid rank\n");
        }
        nlocs = (int) dims[0];
        if (i == 0)
        {
            evlas = memory_calloc64f(nlocs);
            evlos = memory_calloc64f(nlocs);
            evdps = memory_calloc64f(nlocs);
        }
        if (i == 0)
        {
            data = evlas;
        }
        else if (i == 1)
        {
            data = evlos;
        }
        else
        {
            data = evdps;
        }
        memSpace = H5Screate_simple(rank, dims, NULL);
        H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
                H5P_DEFAULT, data);
        H5Sclose(memSpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        data = NULL;
    }
    H5Gclose(groupID);
    // find the data
    for (i=0; i<nlocs; i++)
    {
        if (fabs(evdps[i] - evdp) < 1.e-10)
        {
            if (fabs(evlas[i] - evla) < 1.e-10)
            {
                if (fabs(evlos[i] - evlo) < 1.e-10 ||
                    fabs(evlos[i] + 360.0 - evlo) < 1.e-10 ||
                    fabs(evlos[i] - 366.0 - evlo) < 1.e-10)
                {
                    iloc = i;
                    break;
                }
            }
        }
    }
    memory_free64f(&evdps);
    memory_free64f(&evlas);
    memory_free64f(&evlos);
    return iloc;
}
//============================================================================//
int utils_dataArchive_readAllWaveforms(const char *dataFile,
                                       struct parmtData_struct *data)
{
    const char *fcnm = "utils_dataArchive_readAllWaveforms\0";
    char varname[512];
    hid_t dataSet, dataSpace, groupID, h5fl;
    hsize_t dims[1];
    int ierr, iloc, iobs, k, ndims, nlocs, nobs, nread;
    // Check the file exists and open it
    ierr = 0;
    memset(data, 0, sizeof(struct parmtData_struct));
    if (!os_path_isfile(dataFile))
    {
        log_errorF("%s: Error data file %s does not exist\n", fcnm, dataFile);
        return -1;
    }
    h5fl = H5Fopen(dataFile, H5F_ACC_RDONLY, H5P_DEFAULT);
    // Get the number of locations
    dataSet = H5Dopen(h5fl, "/Locations/Depths\0", H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    ndims = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    if (ndims != 1)
    {
        log_errorF("%s: Error only 1D arrays supported\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    nlocs = (int) dims[0];
    if (nlocs < 1)
    {
        log_errorF("%s: Error no locations in grid-search\n", fcnm);
        ierr = 1;
        goto ERROR;
    } 
    // Get the number of observations
    nobs = utils_dataArchive_getNumberOfObservations(h5fl);
    if (nobs < 1)
    {
        printf("%s: Error no observations\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    nread = nobs*nlocs;
    if (nread < 1)
    {
        log_errorF("%s: Error nothing to read\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    data->nlocs = nlocs;
    data->nobs = nobs;
    data->data = (struct sacData_struct *)
                 calloc((size_t) nobs, sizeof(struct sacData_struct));
    data->sacGxx = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    data->sacGyy = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    data->sacGzz = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    data->sacGxy = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    data->sacGxz = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    data->sacGyz = (struct sacData_struct *)
                   calloc((size_t) nread, sizeof(struct sacData_struct));
    // read the greens functions
    k = 0;
    for (iobs=0; iobs<nobs; iobs++)
    {
        // read the data
        memset(varname, 0, 512*sizeof(char));
        sprintf(varname, "Observations/Observation_%d", iobs); 
        groupID = H5Gopen2(h5fl, varname, H5P_DEFAULT);
        ierr = sacioh5_readTimeSeries2("Observation", groupID,
                                       &data->data[iobs]);
        if (ierr != 0)
        {
            printf("%s: Error reading data\n", fcnm);
            goto ERROR;
        }
        H5Gclose(groupID);
        for (iloc=0; iloc<nlocs; iloc++)
        {
            k = iobs*nlocs + iloc;
            ierr = utils_dataArchive_loadGreensFunctions(h5fl,
                                                         iobs, iloc,
                                                         &data->sacGxx[k],
                                                         &data->sacGyy[k],
                                                         &data->sacGzz[k],
                                                         &data->sacGxy[k],
                                                         &data->sacGxz[k],
                                                         &data->sacGyz[k]);
            if (ierr != 0)
            {
                log_errorF("%s: Error loading obs/location %d %d\n",
                           fcnm, iobs, iloc);
                goto ERROR;
            }
            //sacio_freeData(struct sacData_struct *sac);
        }
    }
ERROR:;
    H5Fclose(h5fl);
    return ierr;
}
//============================================================================//
/*!
 * @brief Defines the initial H5 archive which will contain the data,
 *        Green's functions, and parameters
 *
 */
int utils_dataArchive_initialize(const char *dirnm, const char *projnm,
                                 const int nlocs,
                                 const double *__restrict__ evlas,
                                 const double *__restrict__ evlos,
                                 const double *__restrict__ evdps)
{
    const char *fcnm = "utils_dataArchive_initialize\0";
    char fname[PATH_MAX], varname[256];
    int i, ierr;
    const hsize_t dims[1] = {nlocs};
    hid_t attribute, dataSet, dataSpace, fid, groupID, locID;
    // Check that there is an output archive
    if (!os_path_isdir(dirnm))
    {
        ierr = os_makedirs(dirnm);
        if (ierr != 0)
        {
            log_errorF("%s: Error making output directory: %s\n", fcnm, dirnm);
            return -1;
        }
    }
    // Make the filename
    ierr = utils_dataArchive_setFileName(dirnm, projnm, fname);
    if (ierr != 0)
    {
        log_errorF("%s: Error setting filename\n", fcnm);
        return -1;
    }
    if (os_path_isfile(fname))
    {
        log_warnF("%s: Warning overwriting file %s\n", fcnm, fname);
    }
    // Open the H5 file
    fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Create the heirarchy
    //groupID = H5Gcreate2(fid, "/Structures\0",
    //                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //H5Gclose(groupID);
    groupID = H5Gcreate2(fid, "/Observations\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(groupID);
    groupID = H5Gcreate2(fid, "/Locations\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataSpace = H5Screate_simple(1, dims, NULL);
    dataSet = H5Dcreate2(groupID, DEPTHS, H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, evdps);
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(groupID, LATITUDES, H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, evlas);
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(groupID, LONGITUDES, H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, evlos);
    H5Dclose(dataSet);

    H5Sclose(dataSpace);
/*
    for (i=0; i<nlocs; i++)
    {
         memset(varname, 0, 256*sizeof(char));
         sprintf(varname, "Location_%d", i);
         locID = H5Gcreate2(groupID, varname,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         dataSpace = H5Screate_simple(1, dims, NULL);
         attribute = H5Acreate2(locID, "Depth\0", H5T_NATIVE_DOUBLE,
                                dataSpace, H5P_DEFAULT, H5P_DEFAULT);
         H5Awrite(attribute, H5T_NATIVE_DOUBLE, &evdps[i]); 
         H5Aclose(attribute);

         attribute = H5Acreate2(locID, "Latitude\0", H5T_NATIVE_DOUBLE,
                                dataSpace, H5P_DEFAULT, H5P_DEFAULT);
         H5Awrite(attribute, H5T_NATIVE_DOUBLE, &evlas[i]); 
         H5Aclose(attribute);

         attribute = H5Acreate2(locID, "Longitude\0", H5T_NATIVE_DOUBLE,
                                dataSpace, H5P_DEFAULT, H5P_DEFAULT);
         H5Awrite(attribute, H5T_NATIVE_DOUBLE, &evlos[i]);
         H5Aclose(attribute);
        
         H5Sclose(dataSpace);
         H5Gclose(locID);
    }
*/
    H5Gclose(groupID);
    // Close the file
    H5Fclose(fid);
    return 0;
}

int utils_dataArchive_getNumberOfObservations(const hid_t h5fl)
{
    const char *fcnm = "utils_dataArchive_getNumberOfObservations\0";
    char varname[256];
    int nobs;
    int iobs;
    nobs =-1;
    for (iobs=0; iobs<MAXOBS; iobs++)
    {
        memset(varname, 0, 256);
        sprintf(varname, "/Observations/Observation_%d", iobs);
        if (H5Lexists(h5fl, varname, H5P_DEFAULT) <= 0){goto END;}
        nobs = iobs + 1;
    }
    log_errorF("%s: Exceeded max number of observations %d\n", fcnm, MAXOBS);
END:;
    return nobs;
}
//============================================================================//
/*!
 * @brief Given the SAC data structure this returns the corresponding 
 *        observation ID in the H5 archive
 *
 * @param[in] h5fl    HDF5 file handle
 * @param[in] obs     observation to match
 *
 * @result if -1 there are no observations.
 *         if the return value is equal to the number of observations
 *         then this is a new observation. 
 *         otherwise [0,nobs-1] is the observation ID.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int utils_dataArchive_getObservationID(const hid_t h5fl,
                                       const struct sacData_struct obs)
{
    const char *fcnm = "utils_dataArchive_getObservationID\0";
    struct sacData_struct sac;
    char varname[256];
    int ierr, iobs, nobs;
    hid_t groupID, obsID;
    //------------------------------------------------------------------------//
    nobs =-1;
    // Count the number of observations
    nobs = utils_dataArchive_getNumberOfObservations(h5fl) - 1;
    if (nobs < 0)
    {
        log_errorF("%s: No observations\n", fcnm);
        return -1;
    }
    // Loop through the observations and match 
    memset(&sac, 0, sizeof(struct sacData_struct));
    groupID = H5Gopen2(h5fl, "/Observations\0", H5P_DEFAULT);
    for (iobs=0; iobs<nobs; iobs++)
    {
        memset(varname, 0, 256*sizeof(char));
        sprintf(varname, "Observation_%d", iobs);
        // doesn't exist
        if (H5Lexists(h5fl, varname, H5P_DEFAULT) <= 0){break;}
        // open it
        obsID = H5Gopen2(groupID, varname, H5P_DEFAULT);
        ierr = sacioh5_readTimeSeries2("Observation\0", obsID, &sac);
        if (ierr != 0)
        {
            log_errorF("%s: Error loading observation %d\n", fcnm, iobs);
            return -1;
        }
        // check it
        if (strcasecmp(obs.header.knetwk, sac.header.knetwk) == 0 &&
            strcasecmp(obs.header.kstnm,  sac.header.kstnm)  == 0 &&
            strcasecmp(obs.header.kcmpnm, sac.header.kcmpnm) == 0 &&
            strcasecmp(obs.header.khole,  sac.header.khole)  == 0)
        {
            nobs = iobs;
            break;
        }
        sacio_free(&sac);
        H5Gclose(obsID);
    }
    H5Gclose(groupID);
    return nobs;
}

int utils_dataArchive_addObservation(const hid_t h5fl,
                                     const struct sacData_struct obs)
{
    char varname[256];
    hid_t groupID;
    int ierr, iobs, nobs;
    // Get the number of observations
    nobs = utils_dataArchive_getNumberOfObservations(h5fl);
    // no observations
    if (nobs < 0)
    {
        iobs = 0;
    }
    else
    {
        iobs = nobs;
    }
    memset(varname, 0, 256*sizeof(char));
    sprintf(varname, "/Observations/Observation_%d", iobs);
    groupID = H5Gcreate2(h5fl, varname,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ierr = sacioh5_writeTimeSeries2("Observation", groupID, obs);
    if (ierr != 0){printf("Failed to write observation\n");}
    H5Gclose(groupID);
    return ierr;
}

int utils_dataArchive_loadGreensFunctions(const hid_t h5fl,
                                          const int waveformID,
                                          const int locationID,
                                          struct sacData_struct *sacGxx,
                                          struct sacData_struct *sacGyy,
                                          struct sacData_struct *sacGzz,
                                          struct sacData_struct *sacGxy,
                                          struct sacData_struct *sacGxz,
                                          struct sacData_struct *sacGyz)
{
    const char *fcnm = "utils_dataArchive_loadGreensFunctions\0";
    char varname[256];
    int ierr;
    hid_t groupID;
    ierr = 0;
    memset(sacGxx, 0, sizeof(struct sacData_struct));
    memset(sacGyy, 0, sizeof(struct sacData_struct));
    memset(sacGzz, 0, sizeof(struct sacData_struct));
    memset(sacGxy, 0, sizeof(struct sacData_struct));
    memset(sacGxz, 0, sizeof(struct sacData_struct));
    memset(sacGyz, 0, sizeof(struct sacData_struct));
    sprintf(varname, "/Observations/Observation_%d", waveformID);
    if (H5Lexists(h5fl, varname, H5P_DEFAULT) <= 0)
    {
        printf("%s: Observation does not exist\n", fcnm);
        return -1;
    }
    groupID = H5Gopen2(h5fl, varname, H5P_DEFAULT);
    sprintf(varname, "Gxx_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGxx);
    if (ierr != 0)
    {
        printf("%s: Error reading Gxx\n", fcnm);
        goto ERROR;
    }
    sprintf(varname, "Gyy_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGyy);
    if (ierr != 0)
    {
        printf("%s: Error reading Gyy\n", fcnm);
        goto ERROR;
    }
    sprintf(varname, "Gzz_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGzz);
    if (ierr != 0)
    {
        printf("%s: Error reading Gzz\n", fcnm);
        goto ERROR;
    }
    sprintf(varname, "Gxy_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGxy);
    if (ierr != 0)
    {
        printf("%s: Error reading Gxy\n", fcnm);
        goto ERROR;
    }
    sprintf(varname, "Gxz_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGxz);
    if (ierr != 0)
    {
        printf("%s: Error reading Gxz\n", fcnm);
        goto ERROR;
    }
    sprintf(varname, "Gyz_%d", locationID);
    ierr = sacioh5_readTimeSeries2(varname, groupID, sacGyz);
    if (ierr != 0)
    {
        printf("%s: Error reading Gyz\n", fcnm);
        goto ERROR;
    }

ERROR:;
    H5Gclose(groupID);
    return ierr;
}

int utils_dataArchive_addGreensFunctions(const hid_t h5fl,
                                         const int waveformID,
                                         const int locationID,
                                         const struct sacData_struct sacGxx,
                                         const struct sacData_struct sacGyy,
                                         const struct sacData_struct sacGzz,
                                         const struct sacData_struct sacGxy,
                                         const struct sacData_struct sacGxz,
                                         const struct sacData_struct sacGyz)
{
    const char *fcnm = "utils_dataArchive_addGreensFunctions\0";
    char varname[256];
    int ierr;
    hid_t groupID;
    sprintf(varname, "/Observations/Observation_%d", waveformID);
    if (H5Lexists(h5fl, varname, H5P_DEFAULT) <= 0)
    {
        printf("%s: Observation does not exist\n", fcnm);
        return -1;
    }
    groupID = H5Gopen2(h5fl, varname, H5P_DEFAULT);
    sprintf(varname, "Gxx_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGxx);
    if (ierr != 0)
    {
        printf("%s: Error writing Gxx\n", fcnm);
        goto ERROR;
    }

    sprintf(varname, "Gyy_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGyy);
    if (ierr != 0)
    {   
        printf("%s: Error writing Gyy\n", fcnm);
        goto ERROR;
    }

    sprintf(varname, "Gzz_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGzz);
    if (ierr != 0)
    {
        printf("%s: Error writing Gzz\n", fcnm);
        goto ERROR;
    }

    sprintf(varname, "Gxy_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGxy);
    if (ierr != 0)
    {
        printf("%s: Error writing Gxy\n", fcnm);
        goto ERROR;
    }

    sprintf(varname, "Gxz_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGxz);
    if (ierr != 0)
    {
        printf("%s: Error writing Gxz\n", fcnm);
        goto ERROR;
    }

    sprintf(varname, "Gyz_%d", locationID);
    ierr = sacioh5_writeTimeSeries2(varname, groupID, sacGyz);
    if (ierr != 0)
    {
        printf("%s: Error writing Gyz\n", fcnm);
        goto ERROR;
    }

ERROR:;
    H5Gclose(groupID); 
    return ierr;
}


