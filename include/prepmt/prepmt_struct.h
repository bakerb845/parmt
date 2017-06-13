#ifndef PREPMT_STRUCT_H__
#define PREPMT_STRUCT_H__ 1
#include <stdbool.h>
#include <compearth.h>
#include "sacio.h"
#include "prepmt/prepmt_config.h"

//#define PREPMT_MIN_TELESEISMIC_DIST 30.0
//#define PREPMT_MAX_TELESEISMIC_DIST 95.0

enum prepmtGreens_enum
{
    G11_GRNS = 1,  /*!< Green's function that scales xx moment tensor term */
    G22_GRNS = 2,  /*!< Green's function that scales yy moment tensor term */
    G33_GRNS = 3,  /*!< Green's function that scales zz moment tensor term */
    G12_GRNS = 4,  /*!< Green's function that scales xy moment tensor term */
    G13_GRNS = 5,  /*!< Green's function that scales xz moment tensor term */
    G23_GRNS = 6   /*!< Green's function that scales yz moment tensor term */
};

struct prepmtModifyCommands_struct
{
    double targetDt;     /*!< Target sampling period (s) in decimation. */
    double cut0;         /*!< Relative cut start time. */
    double cut1;         /*!< Relative cut end time. */
    int iodva;           /*!< Output units of Hudson Green's functions:
                               0 -> displacement.
                               1 -> velocity.
                               2 -> acceleration. */
    bool ldeconvolution; /*!< If true transfer function put into deconvolution.
                              This would be appropriate for removing the
                              instrument response from data.
                              If false transfer function put into convolution.
                              This would be appropriate for simulating an
                              instrument response with Green's functions. */ 
                              
};

struct prepmtCommandsChars_struct
{
    char **cmds;
    int ncmds;
    bool lgeneric;
    char pad[3];
};

struct prepmtCommands_struct
{
    struct prepmtCommandsChars_struct *cmds;
    int nobs;
};

struct prepmtEventParms_struct
{
    double latitude;      /*!< Event latitude (degrees). */
    double longitude;     /*!< Event longitude (degrees). */
    double depth;         /*!< Event depth (km). */
    double time;          /*!< Origin epochal time (UTC seconds). */
    double m11;           /*!< Moment tensor component for xx or rr. */
    double m22;           /*!< Moment tensor component for yy or tt. */
    double m33;           /*!< Moment tensor component for zz or pp. */
    double m12;           /*!< Moment tensor component for xy or rt. */
    double m13;           /*!< Moment tensor component for xz or rp. */
    double m23;           /*!< Moment tensor component for yz or tp. */
    enum compearthCoordSystem_enum
         basis;           /*!< Moment tensor basis - NED or USE. */
    char pad[4];
};

#endif
