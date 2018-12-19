cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT OFF)
SET(CTEST_TEST_TYPE Nightly)

# Begin User inputs:
set( CTEST_SITE             "avatar.scorec.rpi.edu" ) # generally the output of hostname
set( CTEST_DASHBOARD_ROOT   "/fasttmp/dibanez/cdash/trilinos" ) # writable path
set( CTEST_CMAKE_GENERATOR  "Unix Makefiles" ) # What is your compilation apps ?
set( CTEST_BUILD_CONFIGURATION  Release) # What type of build do you want ?

set( CTEST_PROJECT_NAME         "Albany" )
set( CTEST_SOURCE_NAME          repos)
set( CTEST_BUILD_NAME           "linux-gcc-${CTEST_BUILD_CONFIGURATION}")
set( CTEST_BINARY_NAME          build)

SET(PREFIX_DIR /users/ghansen)

SET (CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET (CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  FILE(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
ENDIF()
IF(NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  FILE(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
ENDIF()

configure_file(${CTEST_SCRIPT_DIRECTORY}/danCTestConfig.cmake
               ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake COPYONLY)

SET(CTEST_CMAKE_COMMAND "cmake")
SET(CTEST_COMMAND "ctest -D ${CTEST_TEST_TYPE}")
SET(CTEST_BUILD_FLAGS "-j 8")

find_program(CTEST_GIT_COMMAND NAMES git)

# Point at the public Repo
SET(Trilinos_REPOSITORY_LOCATION https://software.sandia.gov/trilinos/repositories/publicTrilinos)
SET(SCOREC_REPOSITORY_LOCATION https://github.com/SCOREC/core.git)
SET(Albany_REPOSITORY_LOCATION https://github.com/SNLComputation/Albany.git)

# Initial cache info
set( CACHE_CONTENTS "
SITE:STRING=${CTEST_SITE}
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_GENERATOR:INTERNAL=${CTEST_CMAKE_GENERATOR}
BUILD_TESTING:BOOL=OFF
PRODUCT_REPO:STRING=${Albany_REPOSITORY_LOCATION}
" )

ctest_empty_binary_directory( "${CTEST_BINARY_DIRECTORY}" )
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "${CACHE_CONTENTS}")


# Get the publicTrilinos repo

set(CTEST_CHECKOUT_COMMAND)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/publicTrilinos")
  message("Cloning Trilinos repository!")
  EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
    clone ${Trilinos_REPOSITORY_LOCATION} ${CTEST_SOURCE_DIRECTORY}/publicTrilinos
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE HAD_ERROR)
  message("out: ${_out}")
  message("err: ${_err}")
  message("res: ${HAD_ERROR}")
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot clone Trilinos repository!")
  endif()
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

# Get the SCOREC repo

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/publicTrilinos/SCOREC")
  message("Cloning SCOREC repository!")
  EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
    clone ${SCOREC_REPOSITORY_LOCATION} ${CTEST_SOURCE_DIRECTORY}/publicTrilinos/SCOREC
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE HAD_ERROR)
  message("out: ${_out}")
  message("err: ${_err}")
  message("res: ${HAD_ERROR}")
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot checkout SCOREC repository!")
  endif()
endif()

# Get Tpetra branch of Albany

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/AlbanyT")
  message("Cloning Albany repository!")
  EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
    clone -b scorec ${Albany_REPOSITORY_LOCATION} ${CTEST_SOURCE_DIRECTORY}/AlbanyT
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE HAD_ERROR)
  message("out: ${_out}")
  message("err: ${_err}")
  message("res: ${HAD_ERROR}")
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot clone Albany repository, Tpetra branch!")
  endif()
endif()

ctest_start(${CTEST_TEST_TYPE})

# Send the project structure to CDash

IF(CTEST_DO_SUBMIT)
  CTEST_SUBMIT(FILES "${CTEST_SCRIPT_DIRECTORY}/danProject.xml"
            RETURN_VALUE  HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit Albany Project.xml!")
  endif()
ENDIF()

# Update Trilinos
SET_PROPERTY (GLOBAL PROPERTY SubProject Trilinos)
SET_PROPERTY (GLOBAL PROPERTY Label Trilinos)

message("Updating Trilinos repository!")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/publicTrilinos" RETURN_VALUE count)
message("Found ${count} changed files")

IF(CTEST_DO_SUBMIT)
  CTEST_SUBMIT(PARTS Update RETURN_VALUE  HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot update Trilinos!")
  endif()
ENDIF()

# Update the SCOREC repo
SET_PROPERTY (GLOBAL PROPERTY SubProject SCOREC)
SET_PROPERTY (GLOBAL PROPERTY Label SCOREC)

message("Updating SCOREC repository!")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/publicTrilinos/SCOREC" RETURN_VALUE count)
message("Found ${count} changed files")

IF(CTEST_DO_SUBMIT)
  CTEST_SUBMIT(PARTS Update RETURN_VALUE  HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot update SCOREC!")
  endif()
ENDIF()

# Update Albany Tpetra branch
SET_PROPERTY (GLOBAL PROPERTY SubProject AlbanyTpetraBranch)
SET_PROPERTY (GLOBAL PROPERTY Label AlbanyTpetraBranch)

message("Updating Albany repository!")
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}/AlbanyT" RETURN_VALUE count)
message("Found ${count} changed files")

IF(CTEST_DO_SUBMIT)
  CTEST_SUBMIT(PARTS Update RETURN_VALUE  HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot update Albany Tpetra branch!")
  endif()
ENDIF()

# Configure the Trilinos/SCOREC build
SET_PROPERTY(GLOBAL PROPERTY SubProject Trilinos)
SET_PROPERTY(GLOBAL PROPERTY Label Trilinos)

SET(CONFIGURE_OPTIONS
  "-DTrilinos_EXTRA_REPOSITORIES:STRING=SCOREC"
  "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SOURCE_DIRECTORY}/publicTrilinos/sampleScripts/AlbanySettings.cmake"
  "-DCMAKE_BUILD_TYPE:STRING=NONE"
  "-DCMAKE_CXX_FLAGS:STRING=-O3"
  "-DCMAKE_C_FLAGS:STRING=-O3"
  "-DCMAKE_Fortran_FLAGS:STRING=-O3"
  "-DTPL_ENABLE_MPI:BOOL=ON"
  "-DMPI_BASE_DIR:PATH=${PREFIX_DIR}"
  "-DTPL_ENABLE_Matio:BOOL=OFF"
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
  "-DTrilinos_VERBOSE_CONFIGURE:BOOL=OFF"
  "-DBoost_INCLUDE_DIRS:PATH=${PREFIX_DIR}/include"
  "-DBoostAlbLib_INCLUDE_DIRS:PATH=${PREFIX_DIR}/include"
  "-DBoost_LIBRARY_DIRS:PATH=${PREFIX_DIR}/lib"
  "-DBoostAlbLib_LIBRARY_DIRS:PATH=${PREFIX_DIR}/lib"
  "-DTPL_ENABLE_Netcdf:STRING=ON"
  "-DNetcdf_INCLUDE_DIRS:PATH=${PREFIX_DIR}/parallel/include"
  "-DNetcdf_LIBRARY_DIRS:PATH=${PREFIX_DIR}/parallel/lib"
  "-DTPL_ENABLE_HDF5:STRING=ON"
  "-DHDF5_INCLUDE_DIRS:PATH=${PREFIX_DIR}/parallel/include"
  "-DHDF5_LIBRARY_DIRS:PATH=${PREFIX_DIR}/parallel/lib"
  "-DTPL_ENABLE_Zlib:STRING=ON"
  "-DZlib_INCLUDE_DIRS:PATH=${PREFIX_DIR}/include"
  "-DZlib_LIBRARY_DIRS:PATH=${PREFIX_DIR}/lib"
  "-DTPL_ENABLE_ParMETIS:STRING=ON"
  "-DParMETIS_INCLUDE_DIRS:PATH=${PREFIX_DIR}/parallel/ParMetis-4.0.3/include"
  "-DParMETIS_LIBRARY_DIRS:PATH=${PREFIX_DIR}/parallel/ParMetis-4.0.3/lib"
  "-DTrilinos_ENABLE_SCOREC:BOOL=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON"
  "-DTrilinos_ENABLE_Ifpack2:BOOL=ON"
  "-DTrilinos_ENABLE_Amesos2:BOOL=ON"
  "-DTrilinos_ENABLE_MueLu:BOOL=ON"
  "-DZoltan_ENABLE_ULONG_IDS:BOOL=ON"
  "-DTeuchos_ENABLE_COMPLEX:BOOL=OFF"
  "-DSCOREC_DISABLE_STRONG_WARNINGS:BOOL=ON"
  "-DTPL_ENABLE_Parasolid:BOOL=ON"
  "-DParasolid_INCLUDE_DIRS:PATH=/usr/local/parasolid/25.1.181"
  "-DParasolid_LIBRARY_DIRS:PATH=/usr/local/parasolid/25.1.181/shared_object"
  "-DTPL_ENABLE_SuperLU:STRING=ON"
  "-DSuperLU_INCLUDE_DIRS:PATH=${PREFIX_DIR}/SuperLU_4.3/include"
  "-DSuperLU_LIBRARY_DIRS:PATH=${PREFIX_DIR}/SuperLU_4.3/lib"
  "-DCMAKE_INSTALL_PREFIX:PATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
  "-DTrilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF"
  "-DTrilinos_ENABLE_Kokkos:BOOL=ON"
  "-DTrilinos_ENABLE_KokkosCore:BOOL=ON"
  "-DPhalanx_KOKKOS_DEVICE_TYPE:STRING="SERIAL""
  "-DPhalanx_INDEX_SIZE_TYPE:STRING="INT""
  "-DPhalanx_SHOW_DEPRECATED_WARNINGS:BOOL=OFF"
  "-DKokkos_ENABLE_Serial:BOOL=ON"
  "-DKokkos_ENABLE_OpenMP:BOOL=OFF"
  "-DKokkos_ENABLE_Pthread:BOOL=OFF"
  )

message("Configuring Trilinos repository!")
CTEST_CONFIGURE(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          SOURCE "${CTEST_SOURCE_DIRECTORY}/publicTrilinos"
          OPTIONS "${CONFIGURE_OPTIONS}"
          RETURN_VALUE HAD_ERROR
)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot configure Trilinos build!")
endif()

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Configure RETURN_VALUE  HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit Trilinos/SCOREC configure results!")
endif()
ENDIF()

# SCOREC build
SET_PROPERTY (GLOBAL PROPERTY SubProject SCOREC)
SET_PROPERTY (GLOBAL PROPERTY Label SCOREC)
SET(CTEST_BUILD_TARGET "SCOREC_libs")

MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")

CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          RETURN_VALUE  HAD_ERROR
          NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot build SCOREC!")
endif()

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Build RETURN_VALUE  HAD_ERROR)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit SCOREC build results!")
endif()
ENDIF()

# Trilinos
SET_PROPERTY (GLOBAL PROPERTY SubProject Trilinos)
SET_PROPERTY (GLOBAL PROPERTY Label Trilinos)
SET(CTEST_BUILD_TARGET all)

MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")

CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          RETURN_VALUE  HAD_ERROR
          NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
          APPEND
)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot build Trilinos!")
endif()

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Build RETURN_VALUE  HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit Trilinos build results!")
endif()
ENDIF()

execute_process(COMMAND "${CMAKE_MAKE_PROGRAM}" "install" 
  WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
  RESULT_VARIABLE makeInstallResult
  OUTPUT_VARIABLE makeInstallLog 
  ERROR_VARIABLE makeInstallLog
)

file(WRITE ${CTEST_BINARY_DIRECTORY}/makeinstall.log
  "${makeInstallLog}")

# Configure the Albany Tpetra branch build

SET_PROPERTY (GLOBAL PROPERTY SubProject AlbanyTpetraBranch)
SET_PROPERTY (GLOBAL PROPERTY Label AlbanyTpetraBranch)

# refer to doc/buildScripts/gahanse/do-cmake-albany-mpi
SET(CONFIGURE_OPTIONS
  "-DALBANY_TRILINOS_DIR:PATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
  "-DENABLE_LCM:BOOL=ON"
  "-DENABLE_SCOREC:BOOL=ON"
  "-DENABLE_SG:BOOL=OFF"
  "-DENABLE_ENSEMBLE:BOOL=OFF"
  "-DENABLE_QCAD:BOOL=OFF"
  "-DENABLE_LCM_SPECULATIVE:BOOL=OFF"
  "-DENABLE_LAME:BOOL=OFF"
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
  "-DENABLE_DEBUGGING:BOOL=OFF"
  "-DENABLE_CHECK_FPE:BOOL=OFF"
  "-DENABLE_HYDRIDE:BOOL=OFF"
  "-DENABLE_LANDICE:BOOL=OFF"
  "-DENABLE_MOR:BOOL=OFF"
  "-DENABLE_ALBANY_CI:BOOL=OFF"
  "-DENABLE_ASCR:BOOL=OFF"
  )

if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}/AlbanyT")
  FILE(MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/AlbanyT)
endif()

message("Configuring Albany!")
CTEST_CONFIGURE(
          BUILD "${CTEST_BINARY_DIRECTORY}/AlbanyT"
          SOURCE "${CTEST_SOURCE_DIRECTORY}/AlbanyT"
          OPTIONS "${CONFIGURE_OPTIONS}"
          RETURN_VALUE HAD_ERROR
          APPEND
)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot configure Albany Tpetra branch build!")
endif()

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Configure
          RETURN_VALUE  HAD_ERROR
            )

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit Albany Tpetra branch configure results!")
endif()
ENDIF()

# Build Albany Tpetra branch

message("Building Albany!")
CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}/AlbanyT"
          RETURN_VALUE  HAD_ERROR
          NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
          APPEND
)

if(HAD_ERROR)
  message(FATAL_ERROR "Cannot build Albany Tpetra branch!")
endif()

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Build RETURN_VALUE  HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit Albany Tpetra branch build results!")
endif()
ENDIF()

# Run Albany Tpetra branch tests

message("Testing Albany!")
CTEST_TEST(BUILD "${CTEST_BINARY_DIRECTORY}/AlbanyT")

IF(CTEST_DO_SUBMIT)
CTEST_SUBMIT(PARTS Test RETURN_VALUE  HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Cannot submit Albany Tpetra branch test results!")
endif()
ENDIF()

