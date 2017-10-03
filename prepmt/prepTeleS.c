#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <limits.h>
#include <iniparser.h>
#include "prepmt/prepmt.h"
#include "iscl/iscl/iscl.h"
#include "iscl/os/os.h"
 
#define PROGRAM_NAME "xprepTeleS"
static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256]);

int main(int argc, char **argv)
{
    char iniFile[PATH_MAX], section[256];
    const char *hpulseSection = "hpulse96\0";
    const bool lisP = false;
    int ierr;
    const double dmin = PREPMT_S_MIN_TELESEISMIC_DIST;
    const double dmax = PREPMT_S_MAX_TELESEISMIC_DIST;
    iscl_init();
    ierr = parseArguments(argc, argv, iniFile, section);
    if (ierr != 0){return EXIT_FAILURE;}
    ierr = prepmt_prepData_prepTeleseismicBodyWaves(iniFile, section,
                                                    hpulseSection,lisP,
                                                    dmin, dmax);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Error preparing data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    iscl_finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Parses the command line arguments for the ini file.
 *
 * @param[in] argc      Number of arguments input to the program.
 * @param[in] argv      Command line arguments.
 *
 * @param[out] iniFile  If the result is 0 then this is the ini file.
 * @param[out] section  Section of ini file to read.  The default is
 *                      prepmt:prepTeleseismicData.
 *
 * @result 0 indicates success. \n
 *        -1 indicates the user inquired about the usage. \n
 *        -2 indicates the command line argument is invalid.
 *
 * @author Ben Baker, ISTI
 *
 */
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX], char section[256])
{
    bool linFile, lsection;
    linFile = false;
    lsection = false;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    memset(section, 0, 256*sizeof(char));
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
        c = getopt_long(argc, argv, "?hi:s:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            linFile = true;
        }
        else if (c == 's')
        {
            strcpy(section, (const char *) optarg);
            lsection = true;
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
    if (!lsection){strcpy(section, "prepmt:prepTeleseismicData\0");} 
    return 0;
}
//============================================================================//
/*!
 * @brief Utility routine to demonstrate how to call the program.
 */
static void printUsage(void)
{
    printf("Usage:\n   %s -i ini_file\n\n", PROGRAM_NAME);
    printf("Required arguments:\n");
    printf("    -i ini_file specifies the initialization file\n");
    printf("\n");
    printf("Optional arguments:\n");
    printf("    -h displays this message\n"); 
    printf("    -s section of ini file to read\n");
    return;
}
