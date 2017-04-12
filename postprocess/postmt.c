#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include "compearth.h"
#include "parmt_utils.h"
#include "parmt_postProcess.h"

#define PROGRAM_NAME "postmt"

static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX]);
static void printUsage(void);

int main(int argc, char *argv[])
{
    struct parmtGeneralParms_struct parms;
    char iniFile[PATH_MAX];
    double *phi;
    int ierr;
    // Initialize
    memset(&parms, 0, sizeof(struct parmtGeneralParms_struct));
    // Parse the input arguments 
    ierr = parseArguments(argc, argv, iniFile);
    if (ierr != 0)
    {
        if (ierr ==-2){return 0;}
        printf("%s: Error parsing arguments\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Load the ini file - from this we should be able to deduce the archive
    ierr = parmt_utils_readGeneralParms(iniFile, &parms);
    if (ierr != 0)
    {
        printf("%s: Error reading general parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Read the archive
     
    return 0;
}

//============================================================================//
/*!
 * @brief Parses the input arguments to get the ini file anme
 */
static int parseArguments(int argc, char *argv[], char iniFile[PATH_MAX])
{
    bool linFile;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"ini_file", required_argument, 0, 'i'},
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
    return 0;
}

static void printUsage(void)
{
    printf("Usage:\n   postmt -i input_file\n\n");
    printf("Required arguments:\n");
    printf("   -i input_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("   -h displays this message\n");
    return;
}

