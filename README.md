# About

These are utilities for performing moment tensor estimation via grid-search.

# Dependencies
Below is a list of dependencies.  It is strongly recommended to use vendor BLAS/LAPACK and MPI where available.

  - CMake v2.6
  - A C11 compiler
  - [Message Passing Interface](https://www.open-mpi.org/)
  - [Intel Performance Primitives](https://software.intel.com/en-us/intel-ipp)
  - [LAPACK(E)](http://www.netlib.org/lapack/) and [(C)BLAS](http://www.netlib.org/blas/).  If using Intel then the [MKL](https://software.intel.com/en-us/mkl) should be used instead of the libraries available through a package manager.  If using IBM Power architecture then the [ESSL](https://www-03.ibm.com/systems/power/software/essl/) library should be used.
  - The high performance self describing file format [HDF5](https://support.hdfgroup.org/HDF5/).  This build will likely require [zlib](https://zlib.net/).
  - The initializing parser library [iniparser](https://github.com/ndevilla/iniparser) 
  - ISTI's scientific computing library [ISCL](https://github.com/bakerb845/libiscl).  This must have been compiled with [GeographicLib](https://geographiclib.sourceforge.io/)
  - ISTI's signals processing library [ISPL](https://github.com/bakerb845/ispl)
  - The 1D global earth [travel-times library](https://gitlab.isti.com/bbaker/libttimes) to [iaspei-tau](https://seiscode.iris.washington.edu/projects/iaspei-tau)
  - The SAC input/output library [sacio](https://gitlab.isti.com/bbaker/sacio)
  - The Computer Programs in Seismology library [libcps](https://github.com/bakerb845/libcps) 
  - The C moment tensor manipulation library [compearth](https://github.com/bakerb845/compearth)
  
# Configuration

The following is a template for configuring parmt to build against the Intel libraries.  In the source root directory one could create a script similar to

    #!/bin/sh
    /usr/bin/cmake ./ \
    -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=/path/to/C/compiler \
    -DMPI_C_COMPILER=/opt/intel/impi/2017.1.132/bin64/mpiicc \
    -DMPI_C_INCLUDE_PATH=/opt/intel/impi/2017.1.132/include64 \
    -DMPI_C_LIBRARIES=/opt/intel/impi/2017.1.132/lib64/libmpi.so \
    -DCMAKE_C_FLAGS="-g3 -O2 -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
    -DPARMT_USE_INTEL=TRUE \
    -DMKL_LIBRARY="/path/to/intel64_lin/libmkl_intel_lp64.so;/path/to/intel64_lin/libmkl_sequential.so;/path/to/intel64_lin/libmkl_core.so" \
    -DIPP_LIBRARY="/path/to/intel64_lin/libipps.so;/path/to/intel64_lin/libippvm.so;/path/to/intel64_lin/libippcore.so" \
    -DIPP_INCLUDE_DIR=/path/to/ipp/include \
    -DMKL_INCLUDE_DIR=/path/to/mkl/include \
    -DCOMPEARTH_INCLUDE_DIR=/path/tocompearth/momenttensor/c_src/include \
    -DCOMPEARTH_LIBRARY=/path/to/compearth/momenttensor/c_src/lib/libcompearth_shared.so \
    -DISCL_LIBRARY=/path/to/libiscl/lib/libiscl_shared.so \
    -DISCL_INCLUDE_DIR=/path/to/libiscl/include \
    -DISPL_LIBRARY=/path/to/ispl/lib/libispl_shared.so \
    -DISPL_INCLUDE_DIR=/path/to/ispl/include \
    -DCPS_LIBRARY=/home/path/to/libcps/lib/libcps_shared.so \
    -DCPS_INCLUDE_DIR=/path/to/libcps/include \
    -DTTIMES_LIBRARY=path/to/libttimes/lib/libttimes_shared.so \
    -DTTIMES_INCLUDE_DIR=/path/to/libttimes/include \
    -DSACIO_INCLUDE_DIR=/path/to/sacio/include \
    -DSACIO_LIBRARY=/home/path/to/lib/libsacio_shared.so \
    -DH5_C_INCLUDE_DIR=/path/to/hdf5/include \
    -DH5_LIBRARY=/path/to/hdf5/lib/libhdf5.so \
    -DINIPARSER_INCLUDE_DIR=/path/to/iniparser/src \
    -DINIPARSER_LIBRARY=/path/to/iniparser/libiniparser.a
 
