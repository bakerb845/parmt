# include(FindLibraryWithDebug)
#if (NOT TTIMES_INCLUDE_DIR AND TTIMES_LIBRARY)
  if (TTIMES_INCLUDE_DIR AND TTIMES_LIBRARY)
    set(TTIMES_FIND_QUIETLY TRUE)
  endif (TTIMES_INCLUDE_DIR AND TTIMES_LIBRARY)
  find_path(TTIMES_INCLUDE_DIR
    NAMES cps.h
    HINTS /usr/include /usr/local/include $ENV{TTIMESDIR}/include
  )
  find_library(TTIMES_LIBRARY
    NAMES cps_shared cps_static
    HINTS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 $ENV{TTIMESDIR}/lib
  )
  if (TTIMES_INCLUDE_DIR)
    message("Found TTIMES include")
  endif (TTIMES_INCLUDE_DIR)
  if (TTIMES_LIBRARY)
    message("Found TTIMES library")
  endif (TTIMES_LIBRARY)
  include(FindPackageHandleStandardArgs)
#endif()
find_package_handle_standard_args(TTIMES DEFAULT_MSG
                                  TTIMES_INCLUDE_DIR TTIMES_LIBRARY)

mark_as_advanced(TTIMES_INCLUDE_DIR TTIMES_LIBRARY)
