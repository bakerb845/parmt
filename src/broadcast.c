#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parmt_utils.h"
#include "parmt_mtsearch.h"
#include "sacio_mpi.h"

int parmt_broadcast_parmtPolarityParms(struct parmtPolarityParms_struct *parms,
                                        const int root, const MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    if (nprocs == 1){return 0;}
    MPI_Bcast(parms->ttimesTablesDir, PATH_MAX, MPI_CHAR, root, comm);
    MPI_Bcast(parms->ttimesModel,     64,       MPI_CHAR, root, comm);
    MPI_Bcast(&parms->lcomputePolarity, 1, MPI_C_BOOL, root, comm);
    return 0;
}

int parmt_broadcast_mtSearchParms(struct parmtMtSearchParms_struct *parms,
                                  const int root, const MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    if (nprocs == 1){return 0;}
    MPI_Bcast(&parms->nb, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->ng, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->nm, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->nk, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->nt, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->ns, 1, MPI_INTEGER, root, comm);
    MPI_Bcast(&parms->betaLower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->betaUpper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->gammaLower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->gammaUpper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->kappaLower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->kappaUpper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->sigmaLower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->sigmaUpper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->thetaLower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->thetaUpper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->m0Lower, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->m0Upper, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&parms->luseLog, 1, MPI_C_BOOL, root, comm);
    return 0; 
}

int parmt_broadcast_generalParms(struct parmtGeneralParms_struct *parms,
                                 const int root, const MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    if (nprocs == 1){return 0;}
    MPI_Bcast(parms->projnm,            256,      MPI_CHAR, root, comm);
    MPI_Bcast(parms->dataFile,          PATH_MAX, MPI_CHAR, root, comm);
    MPI_Bcast(parms->resultsDir,        PATH_MAX, MPI_CHAR, root, comm);
    MPI_Bcast(parms->resultsFileSuffix, PATH_MAX, MPI_CHAR, root, comm);

    MPI_Bcast(&parms->blockSize, 1, MPI_INT, root, comm);
    MPI_Bcast(&parms->objFnType, 1, MPI_INT, root, comm);

    MPI_Bcast(&parms->defaultMaxLagTime, 1, MPI_DOUBLE, root, comm);

    MPI_Bcast(&parms->lwantLags, 1, MPI_C_BOOL, root, comm);
    return 0;
}

int parmt_broadcast_data(struct parmtData_struct *data,
                         const int root, const MPI_Comm comm)
{
    const char *fcnm = "parmt_broadcast_data\0";
    int ierr, myid, nprocs;
    size_t nwork;
    ierr = 0;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &nprocs);
    if (nprocs == 1){return 0;}
    MPI_Bcast(&data->nobs, 1, MPI_INT, root, comm);
    MPI_Bcast(&data->nlocs, 1, MPI_INT, root, comm);
    if (data->nobs == 0){return 0;}
    if (data->nlocs == 0){return 0;}
    if (myid != root)
    {
        nwork = (size_t) (data->nobs*data->nlocs);
        data->sacGxx = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        data->sacGyy = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        data->sacGzz = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        data->sacGxy = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        data->sacGxz = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        data->sacGyz = (struct sacData_struct *)
                       calloc(nwork, sizeof(struct sacData_struct));
        nwork = (size_t) (data->nobs);
        data->data = (struct sacData_struct *)
                     calloc(nwork, sizeof(struct sacData_struct));
    }
    ierr = sacio_mpi_broadcast(data->data, data->nobs, root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting data on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGxx, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gxx on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGyy, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gyy on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGzz, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gzz on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGxy, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gxy on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGxz, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gxz on process %d\n", fcnm, myid);
        return -1;
    }
    ierr = sacio_mpi_broadcast(data->sacGyz, data->nobs*data->nlocs,
                               root, comm);
    if (ierr != 0)
    {
        printf("%s: Error broadcasting Gyz on process %d\n", fcnm, myid);
        return -1;
    }
    return 0;
}
