#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <mpi.h>
#include <iniparser.h>
#include "iscl/array/array.h"
#include "iscl/os/os.h"

int mkgrns_readini(const char *iniFile);

struct mkgrnsParms_struct
{
 
};

int main(int argc, char **argv)
{
    char iniFile[PATH_MAX];
    int ierr, myid, nprocs;
    const int master = 0;
    //------------------------------------------------------------------------//
    //
    // Initialize MPI and some variables
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    ierr = 0;
    if (myid == master)
    {
        memset(iniFile, 0, PATH_MAX*sizeof(char));
strcpy(iniFile, "mkgrns.ini");
        mkgrns_readini(iniFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
//============================================================================//
/*!
 * @brief Parsers the command line arguments
 */
int mkgrns_readini(const char *iniFile)
{
    const char *fcnm = "mkgrns_readini\0";
    const char *s;
    int ierr;
    dictionary *ini;
    //------------------------------------------------------------------------//
    ierr = 0;
    if (!os_path_isfile(iniFile))
    {
        printf("%s: ini %s doesn't exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);


    iniparser_freedict(ini); 
    return ierr;
}
//============================================================================//
