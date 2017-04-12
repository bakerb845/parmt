cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

PROJECT(parmt C)
enable_testing()

set(CMAKE_BUILD_TYPE Release)
# Add modules for finding requisites
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
# For testing
SET(CMAKE_TESTING_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing)
# For binaries
SET(CMAKE_BINARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/bin)
# For libraries
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

IF (PARMT_USE_INTEL)
   MESSAGE("Using Intel Libraries")
   FIND_PACKAGE(IPP REQUIRED)
   FIND_PACKAGE(MKL REQUIRED)
ELSE ()
   FIND_PACKAGE(CBLAS REQUIRED)
   FIND_PACKAGE(LAPACKE REQUIRED)
ENDIF (PARMT_USE_INTEL)
FIND_PACKAGE(ISCL REQUIRED)
FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(CPS REQUIRED)
FIND_PACKAGE(H5 REQUIRED)
FIND_PACKAGE(IniParser REQUIRED)
FIND_PACKAGE(TTimes)
FIND_PACKAGE(Sacio)
ADD_DEFINITIONS(-DPARMT_USE_ISCL)

IF (PARMT_USE_INTEL)
   INCLUDE_DIRECTORIES(
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${CMAKE_CURRENT_SOURCE_DIR}/compearth
   ${MKL_INCLUDE_DIR}
   ${ISCL_INCLUDE_DIR}
   ${IPP_INCLUDE_DIR}
   ${CPS_INCLUDE_DIR}
   ${TTIMES_INCLUDE_DIR}
   ${INIPARSER_INCLUDE_DIR}
   ${SACIO_INCLUDE_DIR}
   ${H5_C_INCLUDE_DIR}
   ${MPI_C_INCLUDE_PATH}
   )
ELSE ()

ENDIF (PARMT_USE_INTEL)

#macro for adding sources
macro (add_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else (_src ${ARGN})
            list (APPEND SRCS "${_src}")
        endif (_relPath)
    endforeach (_src ${ARGN})
    #propigate up to parent
    if (_relPath)
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif (_relPath)
endmacro (add_sources)

# Add directories to compile
add_subdirectory(include)
add_subdirectory(compearth)

add_executable(xmkgrns utils/mkgrns.c)
add_executable(parmt src/parmt.c src/synthetic.c src/splitComm.c)
add_executable(postmt postprocess/postmt.c)
add_executable(utils_test utils/unit_tests.c)
add_executable(objfn_test src/unit_tests.c)
SET (SRCS_PARMT src/broadcast.c src/discretizeMT.c src/locsearch.c
                src/magnitude.c src/mtsearch.c
                src/polarity.c src/obssearch.c
                src/invert.c src/stackedCrossCorrelation.c 
                postprocess/beachball.c
                postprocess/marginal.c
                postprocess/tt2beachball.c)
add_library(compearth_shared SHARED ${SRCS_COMPEARTH})
add_library(compearth_static STATIC ${SRCS_COMPEARTH})
add_library(parmtUtils_shared SHARED utils/ff2mtGreens.c utils/getNoiseBasis.c
            utils/readini.c io/dataArchive.c io/objfnArchive.c utils/sacGrnsToEst.c)
add_library(parmtUtils_static STATIC utils/ff2mtGreens.c utils/getNoiseBasis.c
            utils/readini.c io/dataArchive.c io/objfnArchive.c utils/sacGrnsToEst.c)
add_library(parmt_shared SHARED ${SRCS_PARMT})
add_library(parmt_static STATIC ${SRCS_PARMT})

SET(LINK_LIB_ALL ${INIPARSER_LIBRARY} ${CPS_LIBRARY} ${TTIMES_LIBRARY} ${ISCL_LIBRARY} ${MPI_C_LIBRARIES} ${H5_LIBRARY}
                   -Wl,--start-group ${MKL_LIBRARY}  -Wl,--end-group -lm -ldl -lpthread -lpng -lstdc++)

set_target_properties(xmkgrns PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(parmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(postmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(utils_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
set_target_properties(objfn_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
set_target_properties(compearth_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
set_target_properties(compearth_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
set_target_properties(parmt_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
set_target_properties(parmt_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
set_target_properties(parmtUtils_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
set_target_properties(parmtUtils_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})

TARGET_LINK_LIBRARIES(xmkgrns ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(parmt parmt_shared parmtUtils_shared compearth_shared ${SACIO_LIBRARY} ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(postmt parmt_shared parmtUtils_shared compearth_shared ${SACIO_LIBRARY} ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(utils_test parmtUtils_static ${SACIO_LIBRARY} ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(objfn_test parmt_static parmtUtils_static compearth_static
                      ${SACIO_LIBRARY} ${LINK_LIB_ALL})

ADD_TEST(NAME utils_test 
         WORKING_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY}
         COMMAND utils_test)
ADD_TEST(NAME objfn_test
         WORKING_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY}
         COMMAND objfn_test)
