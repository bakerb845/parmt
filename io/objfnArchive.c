#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <hdf5.h>
#include "iscl/os/os.h"

#define OPEN_FILE 1
#define CREATE_FILE 2
#define COLATITUDES "colatitudes"
#define LONGITUDES "longitudes"
#define SCALAR_MOMENTS "scalarMoments"
#define STRIKES "strikes"
#define SLIPS "slips"
#define DIPS "dips"
#define DEPTHS "depths"

static herr_t writeDoubleArrayWithUnits(
    const hid_t groupID, const char *arrayName, const char *units,
    const int n, const double *__restrict__ x);
static herr_t readDoubleArray(
    const hid_t groupID, const char *arrayName, const int nwork,
    int *n, double *__restrict__ x);

/*!
 * @brief Sets the file name for the output archive
 *
 * @param[in] job          if job == 1 then this file already exists.
 *                         if job == 2 then this file need not exist.
 * @param[in] archiveDir   directory in which archive will reside
 * @param[in] projnm       name of project
 * @param[in] suffix       project suffix
 *
 * @param[out] flname      name of h5 archive file 
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_io_createObjfnArchiveName(const int job,
                                    const char *archiveDir,
                                    const char *projnm,
                                    const char *suffix,
                                    char flname[PATH_MAX])
{
    const char *fcnm = "parmt_io_createObjfnArchiveName\0";
    int ierr;
    size_t lenos;
    memset(flname, 0, PATH_MAX*sizeof(char));
    if (archiveDir == NULL)
    {
        strcpy(flname, "./\0");
    }
    else
    {
        lenos = strlen(archiveDir);
        if (lenos > 0)
        {
            strcpy(flname, archiveDir);
            if (flname[lenos-1] != '/'){flname[lenos] = '/';}
        }
        else
        {
            strcpy(flname, "./\0");
        }
    }
    if (!os_path_isdir(archiveDir))
    {
        if (job == OPEN_FILE)
        {
            printf("%s: Error directory %s does not exist\n", fcnm, flname);
            return -1;
        }
        else
        {
            ierr = os_makedirs(archiveDir);
            if (ierr != 0)
            {
                printf("%s: Failed to create archive dir %s\n", fcnm, flname);
                return -1;
            }
        }
    }
    if (projnm == NULL)
    {
        printf("%s: Error project name cannot be NULL\n", fcnm);
        return -1;
    }
    strcat(flname, projnm); 
    if (suffix == NULL)
    {
        strcat(flname, ".h5\0");
    }
    else
    {
        if (strlen(suffix) > 0)
        {
            strcat(flname, "_\0");
            strcat(flname, suffix);
        }
        strcat(flname, ".h5\0");
    }
    if (job == OPEN_FILE && !os_path_isfile(flname))
    {
        printf("%s: Error archive %s does not exist\n", fcnm, flname);
        return -1; 
    }
    return 0;
}
//============================================================================//
int parmt_io_writeVarianceForWaveform64f(
    const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const int waveformNumber,
    const int npts, const double *__restrict__ var)
{
    const char *fcnm = "parmt_io_writeVarianceForWaveform64f\0";
    char flname[PATH_MAX], arrayName[512];
    int ierr;
    hid_t groupID, h5fl;
    ierr = parmt_io_createObjfnArchiveName(OPEN_FILE,
                                           resultsDir, projnm,
                                           resultsSuffix, flname);
    if (ierr != 0)
    {   
        printf("%s: Failed to open archive\n", fcnm);
        return -1; 
    }
    h5fl = H5Fopen(flname, H5F_ACC_RDWR, H5P_DEFAULT);
    groupID = H5Gopen2(h5fl, "/Variances\0", H5P_DEFAULT);
    memset(arrayName, 0, 512*sizeof(char));
    sprintf(arrayName, "waveformVariance_%d", waveformNumber);
    writeDoubleArrayWithUnits(groupID, arrayName, NULL,  npts, var);
    H5Gclose(groupID);
    H5Fclose(h5fl);
    return 0;
}
//============================================================================//
/*!
 * @brief Writes the the objective function to the initialized archive file 
 *
 * @param[in] resultsDir      directory in which archive will reside
 * @param[in] projnm          name of project
 * @param[in] resultsSuffix   project suffix
 * @param[in] nmt             number of moment tensors
 * @param[in] phi             objective function at all moment tensor points
 *                            [nmt]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_io_writeObjectiveFunction64f(
    const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const int nmt, const double *__restrict__ phi)
{
    const char *fcnm = "parmt_io_writeObjectiveFunction64f\0";
    char flname[PATH_MAX];
    int ierr;
    hid_t dataSet, dataSpace, groupID, h5fl;
    herr_t status;
    ierr = parmt_io_createObjfnArchiveName(OPEN_FILE,
                                           resultsDir, projnm,
                                           resultsSuffix, flname);
    if (ierr != 0)
    {
        printf("%s: Failed to open objective function archive\n", fcnm);
        return -1; 
    }
    h5fl = H5Fopen(flname, H5F_ACC_RDWR, H5P_DEFAULT);
    groupID = H5Gopen2(h5fl, "/ObjectiveFunction\0", H5P_DEFAULT);
    dataSet = H5Dopen2(groupID, "waveformFitFunction\0", H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, dataSpace,
                      H5P_DEFAULT, phi);
    if (status != 0)
    {
        printf("%s: Error writing phi\n", fcnm);
    } 
    status += H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    status += H5Gclose(groupID);
    status += H5Fclose(h5fl);
    if (status != 0)
    {
        printf("%s: Error writing objective function\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the objective function archive file
 *
 * @bug There is a deficiency and locations currently means depths 
 *
 * @param[in] resultsDir      directory in which archive will reside
 * @param[in] projnm          name of project
 * @param[in] resultsSuffix   project suffix
 * @param[in] nobs            number of observations
 * @param[in] nlocs           number of locations
 * @param[in] deps            depths (km) [nlocs]
 * @param[in] nm              number of scalar moments
 * @param[in] M0s             scalar moments (Nm) [nm]
 * @param[in] nb              number of betas
 * @param[in] betas           lune colatitudes (radians) [nb]
 * @param[in] ng              number of gammas
 * @param[in] gammas          lune longitudes (radians) [ng]
 * @param[in] nk              number of kappas
 * @param[in] kappas          strike angles (radians) [nk]
 * @param[in] ns              number of sigmas 
 * @param[in] sigmas          slip angles (radians) [ns] 
 * @param[in] nt              number of thetas
 * @param[in] thetas          dip angles (radians) [nt]
 *
 * @result 0 indicates success success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_io_createObjfnArchive64f(
    const char *resultsDir, const char *projnm, const char *resultsSuffix,
    const int nobs,
    const int nlocs,
    const double *__restrict__ deps,
    const int nm, const double *__restrict__ M0s,
    const int nb, const double *__restrict__ betas,
    const int ng, const double *__restrict__ gammas,
    const int nk, const double *__restrict__ kappas,
    const int ns, const double *__restrict__ sigmas,
    const int nt, const double *__restrict__ thetas)
{
    const char *fcnm = "parmt_io_createObjfnArchive\0";
    char flname[PATH_MAX];
    int ierr;
    hid_t dataSpace, dataSet, groupID, h5fl, plist;
    herr_t status;
    hsize_t dims[7];
    const hsize_t rank = 7;
    const double fillValue =-1;
    status = 0;
    ierr = parmt_io_createObjfnArchiveName(CREATE_FILE,
                                           resultsDir, projnm,
                                           resultsSuffix, flname);
    if (ierr != 0)
    {
        printf("%s: Failed to make objective function archive\n", fcnm);
        return -1;
    }
    if (os_path_isfile(flname))
    {
        printf("%s: Warning - overwriting %s\n", fcnm, flname);
    }
    h5fl = H5Fcreate(flname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Create a hyperslab where we recall the order of the inverse problem is:
    // location, magnitude, beta, gamma, kappa, sigma, theta
    groupID = H5Gcreate2(h5fl, "/ObjectiveFunction\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = (int) nlocs;
    dims[1] = (int) nm;
    dims[2] = (int) nb; 
    dims[3] = (int) ng;
    dims[4] = (int) nk;
    dims[5] = (int) ns;
    dims[6] = (int) nt;
    dataSpace = H5Screate_simple(rank, dims, NULL);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_fill_value(plist, H5T_NATIVE_DOUBLE, &fillValue);
    dataSet = H5Dcreate2(groupID, "waveformFitFunction\0", H5T_NATIVE_DOUBLE,
                         dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status += H5Pclose(plist);
    status += H5Dclose(dataSet);  
    status += H5Sclose(dataSpace);
    status += H5Gclose(groupID);
    if (ierr != 0)
    {
        printf("%s: Error creating dataset\n", fcnm);
        return -1;
    }
    // Create the search space 
    groupID = H5Gcreate2(h5fl, "/SearchSpace\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = writeDoubleArrayWithUnits(groupID, DEPTHS,
                                       "kilometers\0", nlocs, deps);
    status = writeDoubleArrayWithUnits(groupID, COLATITUDES,
                                       "radians\0", nb, betas);
    status = writeDoubleArrayWithUnits(groupID, LONGITUDES,
                                       "radians\0", ng, gammas);
    status = writeDoubleArrayWithUnits(groupID, SCALAR_MOMENTS,
                                       "newtonMeters\0", nm, M0s);
    status = writeDoubleArrayWithUnits(groupID, STRIKES,
                                       "radians\0", nk, kappas);
    status = writeDoubleArrayWithUnits(groupID, SLIPS,
                                       "radians\0", ns, sigmas);
    status = writeDoubleArrayWithUnits(groupID, DIPS,
                                       "radians\0", nt, thetas);
    H5Gclose(groupID);
    // Make a group for variances
    groupID = H5Gcreate2(h5fl, "/Variances\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Gclose(groupID);
    status = H5Fclose(h5fl);
    if (status != 0)
    {
        printf("%s: Error closing file\n", fcnm);
        return -1;
    } 
    return 0;
}
//============================================================================//
/*!
 * @brief Loads the objective function and grid search details
 *
 */
int parmt_io_readObjfnArchive64f(
    const char *resultsDir, const char *projnm, const char *resultsSuffix,
    int *nlocs, double **deps,
    int *nm, double **M0s,
    int *nb, double **betas,
    int *ng, double **gammas,
    int *nk, double **kappas,
    int *ns, double **sigmas,
    int *nt, double **thetas,
    int *nmt, double **phi  )
{
    const char *fcnm = "parmt_io_readObjfnArchive64f\0";
    char flname[PATH_MAX];
    hid_t groupID, h5fl;
    herr_t status;
    int ierr;
    //------------------------------------------------------------------------//
    //
    // set outputs 
    *nm = 0;
    *nb = 0;
    *ng = 0;
    *nk = 0;
    *ns = 0;
    *nt = 0;
    *nmt = 0;

    *M0s = NULL;
    *betas = NULL;
    *gammas = NULL;
    *kappas = NULL;
    *sigmas = NULL;
    *thetas = NULL; 
    *phi = NULL;
    ierr = parmt_io_createObjfnArchiveName(OPEN_FILE,
                                           resultsDir, projnm,
                                           resultsSuffix, flname);
    if (ierr != 0)
    {
        printf("%s: Failed to open objective function archive\n", fcnm);
        return -1;
    }
    h5fl = H5Fopen(flname, H5F_ACC_RDONLY, H5P_DEFAULT);
    groupID = H5Gopen2(h5fl, "/ObjectiveFunction\0", H5P_DEFAULT);

    status = H5Gclose(groupID);

    status = H5Fclose(h5fl); 
    return 0;
}
//============================================================================//
/*!
 * @brief Utility function for writing a double array with corresponding units
 *
 */
static herr_t writeDoubleArrayWithUnits(
    const hid_t groupID, const char *arrayName, const char *units,
    const int n, const double *__restrict__ x)
{
    const char *fcnm = "writeDoubleArrayWithUnits\0";
    hid_t attrID, attrSpace, attrType, dataSet, dataSpace;
    herr_t status;
    hsize_t dims[1] = {(hsize_t) n};
    if (n < 1 || x == NULL)
    {
        printf("%s: Nothing to write\n", fcnm);
        return 0;
    }
    dataSpace = H5Screate_simple(1, dims, NULL);
    dataSet = H5Dcreate2(groupID, arrayName, H5T_NATIVE_DOUBLE,
                         dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
                       H5S_ALL, H5P_DEFAULT, x);
    if (units != NULL)
    {
        attrSpace = H5Screate(H5S_SCALAR);
        attrType = H5Tcopy(H5T_C_S1);
        H5Tset_size(attrType, strlen(units));
        H5Tset_strpad(attrType, H5T_STR_NULLTERM);
        attrID = H5Acreate2(dataSet, "units\0", attrType, attrSpace,
                            H5P_DEFAULT, H5P_DEFAULT);
        status += H5Awrite(attrID, attrType, units);
        status += H5Sclose(attrSpace);
        status += H5Tclose(attrType);
        status += H5Aclose(attrID);
    }
    status += H5Dclose(dataSet);
    status += H5Sclose(dataSpace);
    return status;
}
//============================================================================//
static herr_t readDoubleArray(
    const hid_t groupID, const char *arrayName, const int nwork,
    int *n, double *__restrict__ x)
{
    const char *fcnm = "readDoubleArray\0";
    int rank;
    hid_t dataSet, dataSpace, memSpace;
    hsize_t dims[1];
    herr_t status;
    status = 0;
    if (!H5Lexists(groupID, arrayName, H5P_DEFAULT))
    {
        printf("%s: Error array %s does not exist\n", fcnm, arrayName);
        return -1; 
    }
    // Open the dataset and get sizes
    dataSet = H5Dopen(groupID, arrayName, H5P_DEFAULT); 
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank != 1)
    {
        printf("%s: Error can only handle rank 1 arrays\n", fcnm);
    }
    status = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    *n = (int) dims[0];
    // Workspace query
    if (nwork ==-1)
    {
        status = H5Sclose(dataSpace);
        status = H5Dclose(dataSet);
        return 0;
    }
    // Verify sufficient space exists to hold output
    if (x == NULL || nwork < *n)
    {
        if (x == NULL){printf("%s: Error x is NULL\n", fcnm);}
        if (nwork < *n){printf("%s: Insufficient space in x\n", fcnm);}
    }
    memSpace = H5Screate_simple(rank, dims, NULL);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
                     H5P_DEFAULT, x);
    if (status != 0)
    {
        printf("%s: Error reading dataset\n", fcnm);
    } 
    status += H5Sclose(memSpace);
    status += H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    return status;
}
/*

int parmt_io_openObjfnArchive( )
{

}
*/
