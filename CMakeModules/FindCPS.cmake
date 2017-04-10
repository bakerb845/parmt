# include(FindLibraryWithDebug)
#if (NOT CPS_INCLUDE_DIR AND CPS_LIBRARY)
  if (CPS_INCLUDE_DIR AND CPS_LIBRARY)
    set(CPS_FIND_QUIETLY TRUE)
  endif (CPS_INCLUDE_DIR AND CPS_LIBRARY)
  find_path(CPS_INCLUDE_DIR
    NAMES cps.h
    HINTS /usr/include /usr/local/include $ENV{CPSDIR}/include
  )
  find_library(CPS_LIBRARY
    NAMES cps_shared cps_static
    HINTS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 $ENV{CPSDIR}/lib
  )
  if (CPS_INCLUDE_DIR)
    message("Found CPS include")
  endif (CPS_INCLUDE_DIR)
  if (CPS_LIBRARY)
    message("Found CPS library")
  endif (CPS_LIBRARY)
  include(FindPackageHandleStandardArgs)
#endif()
find_package_handle_standard_args(CPS DEFAULT_MSG
                                  CPS_INCLUDE_DIR CPS_LIBRARY)

mark_as_advanced(CPS_INCLUDE_DIR CPS_LIBRARY)
