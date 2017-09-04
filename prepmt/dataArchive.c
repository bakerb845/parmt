#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "prepmt/prepmt_dataArchive.h"
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

/*!
 * @brief Closes the HDF5 archive.
 *
 * @param[in] h5fl    HDF5 file handle.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_dataArchive_closeArchive(const hid_t h5fl)
{
    hid_t status;
    status = H5Fclose(h5fl);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to close archive\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the archive file.  This archiver assumes the location 
 *        grid-search is purely in depth.
 *
 * @param[in] fname     Name of HDF5 archive.
 * @param[in] ndeps     Number of depths in the grid search.
 * @param[in] evla      Event latitude (degrees).
 * @param[in] evlo      Event longitude (degrees).
 * @param[in] evdps     Event depths (km).  This is an array of lenght [ndeps].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_dataArchive_createArchive(
    const char *fname, //const char *dirnm, const char *projnm,
    const int ndeps,
    const double evla,
    const double evlo,
    const double *evdps)
{
    double *evlas, *evlos;
    int ierr, nlocs;
    nlocs = ndeps;
    evlas = array_set64f(nlocs, evla, &ierr);
    evlos = array_set64f(nlocs, evlo, &ierr);
    ierr = utils_dataArchive_initialize(fname, //dirnm, projnm,
                                        nlocs, evlas, evlos, evdps);
    memory_free64f(&evlas);
    memory_free64f(&evlos);
    return ierr;
}
//============================================================================//
/*!
 * @brief Opens an archive file for read/write access.
 *
 * @param[in] fname    Name of H5 data archive to open.
 *
 * @param[out] ierr    0 indicates success.
 *
 * @result On successful exit this is the HDF5 file handle of the opened
 *         archive.
 *
 */
hid_t prepmt_dataArchive_openArchive(
    const char *fname, int *ierr) //const char *dirnm, const char *projnm, int *ierr)
{
    hid_t h5fl = 0;
/*
    char fname[PATH_MAX];
    hid_t h5fl = 0;
    *ierr = utils_dataArchive_setFileName(dirnm, projnm, fname);
    if (*ierr != 0)
    {
        fprintf(stderr, "%s: Failed to set file name\n", __func__);
        return -1;
    }
*/
    // archive doesn't exist
    if (!os_path_exists(fname))
    {
        fprintf(stderr, "%s: Archive %s doesn't exist\n", __func__, fname);
        *ierr = 1;
        return h5fl;
    }
    h5fl = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    return h5fl;
}
//============================================================================//
/*!
 * @brief Adds an observation to the HDF5 archive.
 *
 * @param[in] h5fl   HDF5 file handle.
 * @param[in] obs    SAC observation to add.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_dataArchive_addObservation(
    const hid_t h5fl, const struct sacData_struct obs)
{
    int ierr;
    ierr = utils_dataArchive_addObservation(h5fl, obs);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to add observation\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Adds a Green's function to the archive.
 *
 * @param[in] h5fl     HDF5 file handle.
 * @param[in] sac      SAC data file.  If this file does not exist in the
 *                     archive then a new observation will be added.
 * @param[in] sacGxx   Green's functions scaling the \$f m_{xx} \f$ moment 
 *                     tensor term.
 * @param[in] sacGyy   Green's functions scaling the \$f m_{yy} \f$ moment 
 *                     tensor term.
 * @param[in] sacGzz   Green's functions scaling the \$f m_{zz} \f$ moment 
 *                     tensor term.
 * @param[in] sacGxy   Green's functions scaling the \$f m_{xy} \f$ moment 
 *                     tensor term.
 * @param[in] sacGxz   Green's functions scaling the \$f m_{xz} \f$ moment 
 *                     tensor term.
 * @param[in] sacGyz   Green's functions scaling the \$f m_{yz} \f$ moment 
 *                     tensor term.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int prepmt_dataArchive_addGreensFunctions(const hid_t h5fl,
                                          const struct sacData_struct sac,
                                          const struct sacData_struct sacGxx,
                                          const struct sacData_struct sacGyy,
                                          const struct sacData_struct sacGzz,
                                          const struct sacData_struct sacGxy,
                                          const struct sacData_struct sacGxz,
                                          const struct sacData_struct sacGyz)
{
    //char varname[256];
    double evla, evlo, evdp;
    int ierr, locationID, waveformID;
    sacio_getFloatHeader(SAC_FLOAT_EVLA, sacGxx.header, &evla);
    sacio_getFloatHeader(SAC_FLOAT_EVLO, sacGxx.header, &evlo);
    sacio_getFloatHeader(SAC_FLOAT_EVDP, sacGxx.header, &evdp);
    // Get the depth index
    locationID = utils_dataArchive_getLocationID(h5fl,
                                                 evla, evlo, evdp);
    if (locationID < 0)
    {
        fprintf(stderr, "%s: Couldn't find location %f %f %f\n",
                __func__, evla, evlo, evdp);
        return -1;
    }
    waveformID = utils_dataArchive_getObservationID(h5fl, sac);
    if (waveformID < 0)
    {
        ierr = utils_dataArchive_addObservation(h5fl, sac);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Couldn't add observation\n", __func__);
            return -1;
        }
        waveformID = utils_dataArchive_getObservationID(h5fl, sac);
    }
    ierr = utils_dataArchive_addGreensFunctions(h5fl,
                                                waveformID, locationID,
                                                sacGxx, sacGyy, sacGzz,
                                                sacGxy, sacGxz, sacGyz);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Couldn't add Green's function\n", __func__);
        return -1;
    }
    return 0;
}


