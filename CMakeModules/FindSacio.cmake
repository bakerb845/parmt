# include(FindLibraryWithDebug)
#if (NOT SACIO_INCLUDE_DIR AND SACIO_LIBRARY)
  if (SACIO_INCLUDE_DIR AND SACIO_LIBRARY)
    set(SACIO_FIND_QUIETLY TRUE)
  endif (SACIO_INCLUDE_DIR AND SACIO_LIBRARY)
  find_path(SACIO_INCLUDE_DIR
    NAMES sacio.h
    HINTS /usr/include /usr/local/include $ENV{SACIO_DIR}/include
  )
  find_library(SACIO_LIBRARY
    NAMES sacio_shared sacio_static
    HINTS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 $ENV{SACIO_DIR}/lib
  )
  if (SACIO_INCLUDE_DIR)
    message("Found sacio include")
  endif (SACIO_INCLUDE_DIR)
  if (SACIO_LIBRARY)
    message("Found sacio library")
  endif (SACIO_LIBRARY)
  include(FindPackageHandleStandardArgs)
#endif()
find_package_handle_standard_args(SACIO DEFAULT_MSG
                                  SACIO_INCLUDE_DIR SACIO_LIBRARY)

mark_as_advanced(SACIO_INCLUDE_DIR SACIO_LIBRARY)
