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
   ${ISPL_INCLUDE_DIR}
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
ADD_SUBDIRECTORY(include)
# TODO - email carl tape and make this into its own repository
ADD_SUBDIRECTORY(compearth)

ADD_EXECUTABLE(parmt src/parmt.c src/synthetic.c)
ADD_EXECUTABLE(polarmt src/polarmt.c)
ADD_EXECUTABLE(postmt postprocess/postmt.c)
ADD_EXECUTABLE(mergemt src/mergemt.c)
ADD_EXECUTABLE(xgrnsTeleB prepmt/grnsTeleB.c)
ADD_EXECUTABLE(hspec96 prepmt/hspec96.c)
ADD_EXECUTABLE(xprepTeleP prepmt/prepTeleP.c)
##ADD_EXECUTABLE(xmkgrns utils/mkgrns.c)
#ADD_EXECUTABLE(utils_test utils/unit_tests.c)
#ADD_EXECUTABLE(objfn_test src/unit_tests.c)

SET (SRCS_PARMT src/broadcast.c src/discretizeMT.c src/locsearch.c
                src/magnitude.c src/mtsearch.c
                src/polarity.c src/obssearch.c
                src/invert.c src/stackedCrossCorrelation.c 
                src/splitComm.c
                postprocess/beachball.c
                postprocess/marginal.c
                postprocess/tt2beachball.c
                postprocess/cellSpacing.c
                postprocess/gmtHelper.c)
SET (SRCS_PREPMT prepmt/commands.c
                 prepmt/event.c
                 prepmt/hudson96.c
                 prepmt/prepData.c
                 prepmt/dataArchive.c
                 prepmt/greens.c
                 prepmt/hpulse96.c
                 prepmt/pickFile.c)
ADD_LIBRARY(compearth_shared SHARED ${SRCS_COMPEARTH})
ADD_LIBRARY(compearth_static STATIC ${SRCS_COMPEARTH})
ADD_LIBRARY(parmtUtils_shared SHARED utils/ff2mtGreens.c utils/getNoiseBasis.c
            utils/lagTime.c utils/readini.c utils/setDataOnG.c utils/weights.c
            io/dataArchive.c io/objfnArchive.c utils/sacGrnsToEst.c)
ADD_LIBRARY(parmtUtils_static STATIC utils/ff2mtGreens.c utils/getNoiseBasis.c
            utils/readini.c io/dataArchive.c io/objfnArchive.c utils/sacGrnsToEst.c)
ADD_LIBRARY(prepmt_shared SHARED ${SRCS_PREPMT})
#ADD_LIBRARY(prepmt_static STATIC ${SRC_PREPMT})
ADD_LIBRARY(parmt_shared SHARED ${SRCS_PARMT})
ADD_LIBRARY(parmt_static STATIC ${SRCS_PARMT})

SET(LINK_LIB_ALL ${INIPARSER_LIBRARY} ${CPS_LIBRARY} ${TTIMES_LIBRARY} ${SACIO_LIBRARY}
                 ${ISCL_LIBRARY} ${MPI_C_LIBRARIES} ${H5_LIBRARY}
                 -Wl,--start-group ${MKL_LIBRARY} -Wl,--end-group 
                 ${IPP_LIBRARY} -lm -ldl -lpthread -lpng -lstdc++)

SET_TARGET_PROPERTIES(xgrnsTeleB PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(hspec96 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(xprepTeleP PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
##SET_TARGET_PROPERTIES(xmkgrns PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#SET_TARGET_PROPERTIES(utils_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(objfn_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})

SET_TARGET_PROPERTIES(parmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(postmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(polarmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(mergemt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET_TARGET_PROPERTIES(compearth_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(compearth_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(parmt_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(parmt_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(parmtUtils_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(parmtUtils_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
SET_TARGET_PROPERTIES(prepmt_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
#SET_TARGET_PROPERTIES(prepmt_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})

TARGET_LINK_LIBRARIES(xgrnsTeleB
                      prepmt_shared parmt_shared parmtUtils_shared compearth_shared
                      ${LINK_LIB_ALL} ${ISPL_LIBRARY})
TARGET_LINK_LIBRARIES(hspec96
                      prepmt_shared parmt_shared parmtUtils_shared compearth_shared
                      ${LINK_LIB_ALL} ${ISPL_LIBRARY})
TARGET_LINK_LIBRARIES(xprepTeleP
                      prepmt_shared parmt_shared parmtUtils_shared compearth_shared
                      ${LINK_LIB_ALL} ${ISPL_LIBRARY})
##TARGET_LINK_LIBRARIES(xmkgrns ${LINK_LIB_ALL})
#TARGET_LINK_LIBRARIES(utils_test parmtUtils_shared ${SACIO_LIBRARY} ${LINK_LIB_ALL})
#TARGET_LINK_LIBRARIES(objfn_test parmt_shared parmtUtils_shared compearth_shared
#                      ${SACIO_LIBRARY} ${LINK_LIB_ALL}

TARGET_LINK_LIBRARIES(parmt parmt_shared parmtUtils_shared compearth_shared ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(postmt parmt_shared parmtUtils_shared compearth_shared ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(polarmt parmt_shared parmtUtils_shared compearth_shared ${LINK_LIB_ALL})
TARGET_LINK_LIBRARIES(mergemt parmt_shared parmtUtils_shared compearth_shared ${LINK_LIB_ALL})

ADD_TEST(NAME utils_test 
         WORKING_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY}
         COMMAND utils_test)
ADD_TEST(NAME objfn_test
         WORKING_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY}
         COMMAND objfn_test)
