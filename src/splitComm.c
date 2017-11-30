#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
/*!
 * @brief Splits the global communicator into a group for a parallel
 *        loop on observations, a parallel loop on locations, and
 *        a parallel loop on moment tensors.
 * 
 * @param[in] globalComm       global MPI communicator to split
 * @param[in] npInObsGroups    number of processes in an observation
 *                             group.
 * @param[in] npInLocGroups    number of processes in a location group.
 * @param[in] npInMTGroups     number of processes in a moment tensor gorup.
 *
 * @param[out] linObsComm      If true then the process is in the observation
 *                             communicator.
 * @param[out] obsComm         observation group communicator handle
 * @param[out] linLocComm      If true then the process is in the location
 *                             communicator.
 * @param[out] locComm         location group communicator handle
 * @param[out] linMTComm       If true then the process is in the moment
 *                             tensor communicator.
 * @param[out] mtComm          moment tensor group communication handle
 *
 * @result 0 indicates success
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int parmt_splitComm(const MPI_Comm globalComm,
                    const int npInObsGroups,
                    const int npInLocGroups,
                    const int npInMTGroups,
                    bool *linObsComm, MPI_Comm *obsComm,
                    bool *linLocComm, MPI_Comm *locComm,
                    bool *linMTComm,  MPI_Comm *mtComm)
{
    const char *fcnm = "parmt_splitComm\0";
    MPI_Group locGroup, locGroupInclude, mtGroup, mtGroupInclude,
              obsGroup, obsGroupInclude;
    int *ids, *locIDs, *mtIDs, *obsIDs, i, i1, i2, id, j, k,
         locCommSize, nprocs, mtCommSize, myGlobalID, obsCommSize, stride;
    // Get global communicator information and check the split can work 
    MPI_Comm_size(globalComm, &nprocs);
    MPI_Comm_rank(globalComm, &myGlobalID);
    if (nprocs < npInObsGroups*npInLocGroups*npInMTGroups)
    {
        printf("%s: Communicator is too small to split\n", fcnm);
        return -1;
    }
    if (nprocs != npInObsGroups*npInLocGroups*npInMTGroups)
    {
        printf("%s: mpi_excl not yet programmed\n", fcnm);
        printf("%s: nprocs must be evenly divisible by npLoc*npMT*npObs\n",
               fcnm);
        return -1;
    }
    // Handle base case
/*
    if (nprocs == 1)
    {
        *linObsComm = true;
        *linLocComm = true;
        *linMTComm = true;
        *obsComm = globalComm; 
        *locComm = globalComm;
        *mtComm = globalComm;
    }
*/
    // Generate the ranks definining how we split the global communicator.
    // The strategy is that the x (fast) axis of a cube defines the 
    // the MT communicator, the y (intermediate) axis of a cube defines
    // locations that each MT is in, and the z axis of a cube defines
    // the observations groups.  So for example, a 2 observations groups
    // split into 3 locations with each location having 4 moment tensor
    // search members (i.e. 2*3*4 = 24 processes in total) would have
    // numbering:
    //   observation group 1:
    //     location 1: 
    //       global PID: 0 1 2 3
    //     location 2:
    //       global PID: 4 5 6 7 
    //     location 3:
    //       global PID: 8 9 10 11
    //     the moment tensor group masters that can communicate on the
    //     location communicator are process 0, 4, and, 8.
    //   observation group 2:
    //     location 1:
    //       global PID: 12 13 14 15
    //     location 2:
    //       global PID: 16 17 18 19
    //     location 3:
    //       global PID: 20 21 22 23
    //     the moment tensor group masters that can communicate on the
    //     location communicate are processes 12, 16, and 20.
    //   the observation group masters that can communicate on the
    //   observation communicator are 0 and 12.
    *linObsComm = false;
    *linLocComm = false;
    *linMTComm = false;
    obsCommSize = 0;
    locCommSize = 0;
    mtCommSize = 0;
    obsIDs = (int *) calloc((size_t) nprocs, sizeof(int));
    locIDs = (int *) calloc((size_t) nprocs, sizeof(int));
    mtIDs  = (int *) calloc((size_t) nprocs, sizeof(int));
    for (k=0; k<npInObsGroups; k++)
    {   
        for (j=0; j<npInLocGroups; j++)
        {
            for (i=0; i<npInMTGroups; i++)
            {
                id = k*npInLocGroups*npInMTGroups + j*npInMTGroups + i;
                if (j == 0 && i == 0)
                {
                    obsIDs[obsCommSize] = id;
                    obsCommSize = obsCommSize + 1;
                    if (id == myGlobalID){*linObsComm = true;}
                }
                i1 = k*npInLocGroups*npInMTGroups;
                i2 = (k + 1)*npInLocGroups*npInMTGroups;
                if (i == 0 && myGlobalID >= i1 && myGlobalID < i2)
                {
                    locIDs[locCommSize] = id;
                    locCommSize = locCommSize + 1;
                    if (id == myGlobalID){*linLocComm = true;}
                }
                i1 = k*npInLocGroups*npInMTGroups + j*npInMTGroups;
                i2 = k*npInLocGroups*npInMTGroups + (j+1)*npInMTGroups;
                if (myGlobalID >= i1 && myGlobalID < i2)
                {
                    mtIDs[mtCommSize] = id;
                    mtCommSize = mtCommSize + 1; 
                    if (id == myGlobalID){*linMTComm = true;}
                }
            }
        }
    }
    // Create the observation group 
    MPI_Comm_group(globalComm, &obsGroup);
    MPI_Group_incl(obsGroup, obsCommSize, obsIDs, &obsGroupInclude);
    MPI_Comm_create(globalComm, obsGroupInclude, obsComm);
//if (*linObsComm){printf("%d %d\n", obsCommSize, myGlobalID);}
    // Create the location group
    MPI_Comm_group(globalComm, &locGroup);
    MPI_Group_incl(locGroup, locCommSize, locIDs, &locGroupInclude);
    MPI_Comm_create(globalComm, locGroupInclude, locComm);
//if (*linLocComm){printf("%d\n", locCommSize);}
    // Create the moment tensor group
    MPI_Comm_group(globalComm, &mtGroup);
    MPI_Group_incl(mtGroup, mtCommSize, mtIDs, &mtGroupInclude);
    MPI_Comm_create(globalComm, mtGroupInclude, mtComm);
//if (*linMTComm){printf("%d\n", mtCommSize);}
    // Free resources
    MPI_Group_free(&mtGroup);
    MPI_Group_free(&locGroup);
    MPI_Group_free(&obsGroup);
    MPI_Group_free(&mtGroupInclude);
    MPI_Group_free(&locGroupInclude);
    MPI_Group_free(&obsGroupInclude);
    free(obsIDs);
    free(locIDs);
    free(mtIDs);
    return 0;
}
