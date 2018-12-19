
#cmake_minimum_required (VERSION 2.8)
set (CTEST_DO_SUBMIT ON)
set (CTEST_TEST_TYPE Nightly)

# What to build and test
set (DOWNLOAD_TRILINOS FALSE)
set (BUILD_TRILINOS FALSE)
set (DOWNLOAD_ALBANY TRUE) 
set (BUILD_ALBANY TRUE) 

# Begin User inputs:
set (CTEST_SITE "waterman.sandia.gov" ) # generally the output of hostname
set (CTEST_DASHBOARD_ROOT "$ENV{TEST_DIRECTORY}" ) # writable path
set (CTEST_SCRIPT_DIRECTORY "$ENV{SCRIPT_DIRECTORY}" ) # where the scripts live
set (CTEST_CMAKE_GENERATOR "Unix Makefiles" ) # What is your compilation apps ?
set (CTEST_CONFIGURATION  Release) # What type of build do you want ?

set (INITIAL_LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})

set (CTEST_PROJECT_NAME "Albany" )
set (CTEST_SOURCE_NAME repos)
set (CTEST_NAME "linux-gcc-${CTEST_BUILD_CONFIGURATION}")
set (CTEST_BINARY_NAME build)


set (CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set (CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

if (NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file (MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif ()
if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  file (MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
endif ()

configure_file (${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake
  ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake COPYONLY)

set (CTEST_NIGHTLY_START_TIME "00:00:00 UTC")
set (CTEST_CMAKE_COMMAND "cmake")
set (CTEST_COMMAND "ctest -D ${CTEST_TEST_TYPE}")
set (CTEST_FLAGS "-j32")
SET (CTEST_BUILD_FLAGS "-j32")

set (CTEST_DROP_METHOD "http")

if (CTEST_DROP_METHOD STREQUAL "http")
  set (CTEST_DROP_SITE "cdash.sandia.gov")
  set (CTEST_PROJECT_NAME "Albany")
  set (CTEST_DROP_LOCATION "/CDash-2-3-0/submit.php?project=Albany")
  set (CTEST_TRIGGER_SITE "")
  set (CTEST_DROP_SITE_CDASH TRUE)
endif ()

find_program (CTEST_GIT_COMMAND NAMES git)

set (Albany_REPOSITORY_LOCATION git@github.com:SNLComputation/Albany.git)
set (Trilinos_REPOSITORY_LOCATION git@github.com:trilinos/Trilinos.git)

#set (NVCC_WRAPPER "$ENV{jenkins_trilinos_dir}/packages/kokkos/config/nvcc_wrapper")
set (NVCC_WRAPPER ${CTEST_SCRIPT_DIRECTORY}/nvcc_wrapper_volta)
set (CUDA_MANAGED_FORCE_DEVICE_ALLOC 1)
set( CUDA_LAUNCH_BLOCKING 1)

#set(BOOST_ROOT "/home/projects/pwr8-rhel73-lsf/boost/1.60.0/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44")
#set(NETCDF_ROOT "/home/projects/pwr8-rhel73-lsf/netcdf/4.4.1/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44") 
#set(PNETCDF_ROOT "/home/projects/pwr8-rhel73-lsf/pnetcdf/1.6.1/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44") 
#set(HDF5_ROOT "/home/projects/pwr8-rhel73-lsf/hdf5/1.8.17/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44") 
#set(BLAS_ROOT "/home/projects/pwr8-rhel73-lsf/openblas/0.2.19/gcc/5.3.0")  
#set(LAPACKROOT "/home/projects/pwr8-rhel73-lsf/openblas/0.2.19/gcc/5.3.0")
#set(ZLIB_DIR "/home/projects/pwr8-rhel73-lsf/zlib/1.2.8") 

if (CLEAN_BUILD)
  # Initial cache info
  set (CACHE_CONTENTS "
  SITE:STRING=${CTEST_SITE}
  CMAKE_TYPE:STRING=Release
  CMAKE_GENERATOR:INTERNAL=${CTEST_CMAKE_GENERATOR}
  TESTING:BOOL=OFF
  PRODUCT_REPO:STRING=${Albany_REPOSITORY_LOCATION}
  " )

  ctest_empty_binary_directory( "${CTEST_BINARY_DIRECTORY}" )
  file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "${CACHE_CONTENTS}")
endif ()

if (DOWNLOAD_TRILINOS)

  set (CTEST_CHECKOUT_COMMAND)
 
  #
  # Get Trilinos
  #
  
  if (NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/Trilinos")
    execute_process (COMMAND "${CTEST_GIT_COMMAND}" 
      clone ${Trilinos_REPOSITORY_LOCATION} -b develop ${CTEST_SOURCE_DIRECTORY}/Trilinos
      OUTPUT_VARIABLE _out
      ERROR_VARIABLE _err
      RESULT_VARIABLE HAD_ERROR)
    message(STATUS "out: ${_out}")
    message(STATUS "err: ${_err}")
    message(STATUS "res: ${HAD_ERROR}")
    if (HAD_ERROR)
      message(FATAL_ERROR "Cannot clone Trilinos repository!")
    endif ()
  endif ()

endif()


if (DOWNLOAD_ALBANY)

  set (CTEST_CHECKOUT_COMMAND)
  set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
  
  #
  # Get Albany
  #

  if (NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/Albany")
    execute_process (COMMAND "${CTEST_GIT_COMMAND}" 
      clone ${Albany_REPOSITORY_LOCATION} -b master ${CTEST_SOURCE_DIRECTORY}/Albany
      OUTPUT_VARIABLE _out
      ERROR_VARIABLE _err
      RESULT_VARIABLE HAD_ERROR)
    
    message(STATUS "out: ${_out}")
    message(STATUS "err: ${_err}")
    message(STATUS "res: ${HAD_ERROR}")
    if (HAD_ERROR)
      message(FATAL_ERROR "Cannot clone Albany repository!")
    endif ()
  endif ()

  set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")


endif ()


ctest_start(${CTEST_TEST_TYPE})

#
# Send the project structure to CDash
#

#if (CTEST_DO_SUBMIT)
#  ctest_submit (FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
#    RETURN_VALUE  HAD_ERROR
#    )

#  if (HAD_ERROR)
#    message(FATAL_ERROR "Cannot submit Albany Project.xml!")
#  endif ()
#endif ()

# 
# Set the common Trilinos config options & build Trilinos
# 

if (BUILD_TRILINOS) 
  message ("ctest state: BUILD_TRILINOS")
  #
  # Configure the Trilinos/SCOREC build
  #
  set_property (GLOBAL PROPERTY SubProject IKTWatermanTrilinosCUDA)
  set_property (GLOBAL PROPERTY Label IKTWatermanTrilinosCUDA)

    #"-DCMAKE_CXX_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpicxx" 
    #"-DCMAKE_C_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpicc" 
    #"-DCMAKE_Fortran_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpifort"
    #"-DMPI_BASE_DIR:PATH=${MPI_BASE_DIR}"
  set (CONFIGURE_OPTIONS
    "-DCMAKE_INSTALL_PREFIX:PATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
    "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
    #
    "-DCMAKE_CXX_COMPILER:FILEPATH=${NVCC_WRAPPER}"
    "-DCMAKE_SKIP_RULE_DEPENDENCY=ON"
    "-DTPL_ENABLE_MPI:BOOL=ON"
    #
    "-DTPL_ENABLE_MPI:BOOL=ON"
    "-DMPI_EXEC=mpirun"
    "-DMPI_EXEC_NUMPROCS_FLAG:STRING=-n"
    #
    "-DTPL_ENABLE_BLAS:BOOL=ON"
    "-DBLAS_LIBRARY_DIRS:PATH=$ENV{BLAS_ROOT}/lib"
    "-DBLAS_LIBRARY_NAMES:STRING=blas"
    #
    "-DTPL_ENABLE_LAPACK:BOOL=ON"
    "-DLAPACK_LIBRARY_DIRS:PATH=$ENV{BLAS_ROOT}/lib"
    "-DLAPACK_LIBRARY_NAMES:STRING=lapack"
    #
    "-DTPL_ENABLE_Boost:BOOL=ON"
    "-DBoost_INCLUDE_DIRS:PATH=$ENV{BOOST_ROOT}/include"
    #
    "-DTPL_ENABLE_BoostLib:BOOL=ON"
    "-DBoostLib_INCLUDE_DIRS:PATH=$ENV{BOOST_ROOT}/include"
    "-DBoostLib_LIBRARY_DIRS:PATH=$ENV{BOOST_ROOT}/lib"
    #
    "-DTrilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF"
    "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF"
    "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
    "-DTrilinos_ENABLE_CXX11:BOOL=ON"
    "-DTrilinos_ENABLE_EXAMPLES:BOOL=OFF"
    "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=OFF"
    "-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="
    #
    "-DHAVE_INTREPID_KOKKOSCORE:BOOL=ON"
    "-DKokkos_ENABLE_Cuda_UVM:BOOL=ON"
    "-DKokkos_ENABLE_EXAMPLES:BOOL=OFF"
    "-DKokkos_ENABLE_OpenMP:BOOL=OFF"
    "-DKokkos_ENABLE_Pthread:BOOL=OFF"
    "-DKokkos_ENABLE_Serial:BOOL=ON"
    "-DKokkos_ENABLE_TESTS:BOOL=OFF"
    "-DKOKKOS_ARCH:STRING='Power9\\;Volta70'"
    "-DTPL_ENABLE_Thrust:BOOL=ON"
    "-DTPL_ENABLE_CUDA:BOOL=ON"
    "-DTPL_ENABLE_CUSPARSE:BOOL=ON"
    #
    "-DAmesos2_ENABLE_KLU2:BOOL=ON"
    "-DEpetraExt_USING_HDF5:BOOL=OFF"
    "-DIntrepid_ENABLE_TESTS:BOOL=OFF"
    "-DIntrepid2_ENABLE_TESTS:BOOL=OFF"
    "-DPhalanx_INDEX_SIZE_TYPE:STRING=UINT"
    "-DPhalanx_KOKKOS_DEVICE_TYPE:STRING=CUDA"
    "-DSacado_ENABLE_COMPLEX:BOOL=OFF"
    "-DTeuchos_ENABLE_COMPLEX:BOOL=OFF"
    "-DTpetra_ENABLE_Kokkos_Refactor:BOOL=ON"
    "-DKokkos_ENABLE_Cuda_Lambda:BOOL=ON"
    "-DTpetra_ASSUME_CUDA_AWARE_MPI:BOOL=ON"
    "-DTpetra_INST_INT_LONG_LONG:BOOL=ON"
    "-DTpetra_INST_INT_INT:BOOL=OFF"
    #
    "-DTPL_ENABLE_Matio:BOOL=OFF"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DTPL_Netcdf_INCLUDE_DIRS:PATH=$ENV{NETCDF_ROOT}/include"
    "-DTPL_Netcdf_LIBRARIES=$ENV{NETCDF_ROOT}/lib/libnetcdf.a"
    "-DTPL_Netcdf_PARALLEL:BOOL=ON"
    "-DTPL_ENABLE_Pnetcdf:STRING=ON"
    "-DTPL_Pnetcdf_INCLUDE_DIRS:PATH=$ENV{PNETCDF_ROOT}/include"
    "-DTPL_Pnetcdf_LIBRARIES=$ENV{PNETCDF_ROOT}/lib/libpnetcdf.a"
    #
    "-DTPL_ENABLE_HDF5:STRING=ON"
    "-DTPL_HDF5_INCLUDE_DIRS:PATH=$ENV{HDF5_ROOT}/include"
    "-DTPL_HDF5_LIBRARIES=$ENV{HDF5_ROOT}/lib/libhdf5_hl.a"
    "-DTrilinos_EXTRA_LINK_FLAGS:STRING='-lnetcdf -lpnetcdf -lhdf5_hl -lhdf5 -lz'"
    "-DTPL_ENABLE_X11:BOOL=OFF"
    #
    "-DTrilinos_ENABLE_Amesos2:BOOL=ON"
    "-DTrilinos_ENABLE_Amesos:BOOL=ON"
    "-DTrilinos_ENABLE_Anasazi:BOOL=ON"
    "-DTrilinos_ENABLE_AztecOO:BOOL=ON"
    "-DTrilinos_ENABLE_Belos:BOOL=ON"
    "-DTrilinos_ENABLE_EXAMPLES:BOOL=OFF"
    "-DTrilinos_ENABLE_Epetra:BOOL=ON"
    "-DTrilinos_ENABLE_EpetraExt:BOOL=ON"
    "-DTrilinos_ENABLE_Ifpack2:BOOL=ON"
    "-DTrilinos_ENABLE_Ifpack:BOOL=ON"
    "-DTrilinos_ENABLE_Intrepid:BOOL=ON"
    "-DTrilinos_ENABLE_Intrepid2:BOOL=ON"
    "-DTrilinos_ENABLE_Kokkos:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosAlgorithms:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosContainers:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosCore:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosExample:BOOL=OFF"
    "-DTrilinos_ENABLE_ML:BOOL=ON"
    "-DTrilinos_ENABLE_OpenMP:BOOL=OFF"
    "-DTrilinos_ENABLE_MiniTensor:BOOL=ON"
    "-DTrilinos_ENABLE_MueLu:BOOL=ON"
    "-DMueLu_ENABLE_Tutorial:BOOL=OFF"
    "-DTrilinos_ENABLE_NOX:BOOL=ON"
    "-DTrilinos_ENABLE_Pamgen:BOOL=ON"
    "-DTrilinos_ENABLE_Phalanx:BOOL=ON"
    "-DTrilinos_ENABLE_Piro:BOOL=ON"
    "-DTrilinos_ENABLE_Rythmos:BOOL=ON"
    "-DTrilinos_ENABLE_SEACAS:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASAprepro_lib:BOOL=OFF"
    "-DTrilinos_ENABLE_STKDoc_tests:BOOL=OFF"
    "-DTrilinos_ENABLE_STKIO:BOOL=ON"
    "-DTrilinos_ENABLE_STKMesh:BOOL=ON"
    "-DTrilinos_ENABLE_Sacado:BOOL=ON"
    "-DTrilinos_ENABLE_Shards:BOOL=ON"
    "-DTrilinos_ENABLE_Stokhos:BOOL=OFF"
    "-DTrilinos_ENABLE_Stratimikos:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
    "-DTrilinos_ENABLE_Teko:BOOL=ON"
    "-DTrilinos_ENABLE_Teuchos:BOOL=ON"
    "-DTrilinos_ENABLE_ThreadPool:BOOL=OFF"
    "-DTrilinos_ENABLE_Thyra:BOOL=ON"
    "-DTrilinos_ENABLE_Tpetra:BOOL=ON"
    "-DTrilinos_ENABLE_Zoltan2:BOOL=ON"
    "-DTrilinos_ENABLE_Zoltan:BOOL=ON"
    "-DTrilinos_ENABLE_Tempus:BOOL=ON"
    "-DTempus_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  )

  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/TriBuild")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/TriBuild)
  endif ()

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/TriBuild"
    SOURCE "$ENV{jenkins_trilinos_dir}"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE HAD_ERROR
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Configure
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Trilinos configure results!")
    endif ()
  endif ()

  if (HAD_ERROR)
    message ("Cannot configure Trilinos build!")
  endif ()

  #
  # Build the rest of Trilinos and install everything
  #

  set_property (GLOBAL PROPERTY SubProject IKTWatermanTrilinosCUDA)
  set_property (GLOBAL PROPERTY Label IKTWatermanTrilinosCUDA)
  #set (CTEST_BUILD_TARGET all)
  set (CTEST_BUILD_TARGET install)

  MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")

  CTEST_BUILD(
    BUILD "${CTEST_BINARY_DIRECTORY}/TriBuild"
    RETURN_VALUE  HAD_ERROR
    NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
    APPEND
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Build
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Trilinos build results!")
    endif ()

  endif ()

  if (HAD_ERROR)
    message ("Cannot build Trilinos!")
  endif ()

  if (BUILD_LIBS_NUM_ERRORS GREATER 0)
    message ("Encountered build errors in Trilinos build. Exiting!")
  endif ()

endif()

if (BUILD_ALBANY)

  # Configure the Albany build 
  #

  set_property (GLOBAL PROPERTY SubProject IKTWatermanAlbanyCUDA)
  set_property (GLOBAL PROPERTY Label IKTWatermanAlbanyCUDA)
  
  set (CONFIGURE_OPTIONS
    "-DALBANY_TRILINOS_DIR:FILEPATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
    "-DENABLE_DEMO_PDES:BOOL=ON"
    "-DENABLE_LANDICE:BOOL=ON"
    "-DENABLE_ALBANY_EPETRA_EXE:BOOL=OFF"
    "-DENABLE_QCAD:BOOL=ON"
    "-DENABLE_LCM:BOOL=OFF"
    "-DENABLE_AERAS:BOOL=ON"
    "-DENABLE_SG:BOOL=OFF"
    "-DENABLE_ENSEMBLE:BOOL=OFF"
    "-DENABLE_ATO:BOOL=OFF"
    "-DENABLE_MOR:BOOL=OFF"
    "-DENABLE_PERFORMANCE_TESTS:BOOL=OFF"
    "-DALBANY_LIBRARIES_ONLY=OFF"
    "-DENABLE_INSTALL:BOOL=OFF"
    "-DENABLE_KOKKOS_UNDER_DEVELOPMENT:BOOL=ON"
    "-DENABLE_TAN_FAD_TYPE:STRING=SLFad"
    "-DALBANY_TAN_SLFAD_SIZE=100"
    )
  
  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/AlbBuild")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/AlbBuild)
  endif ()

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/AlbBuild"
    SOURCE "/home/projects/albany/waterman/repos/Albany"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE HAD_ERROR
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Configure
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Albany configure results!")
    endif ()
  endif ()

  if (HAD_ERROR)
    message ("Cannot configure Albany build!")
  endif ()

  #
  # Build the rest of Albany and install everything
  #

  set_property (GLOBAL PROPERTY SubProject IKTWatermanAlbanyCUDA)
  set_property (GLOBAL PROPERTY Label IKTWatermanAlbanyCUDA)
  set (CTEST_BUILD_TARGET all)
  #set (CTEST_BUILD_TARGET install)

  MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")

  CTEST_BUILD(
    BUILD "${CTEST_BINARY_DIRECTORY}/AlbBuild"
    RETURN_VALUE  HAD_ERROR
    NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
    APPEND
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Build
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Albany build results!")
    endif ()

  endif ()

  if (HAD_ERROR)
    message ("Cannot build Albany!")
  endif ()

  if (BUILD_LIBS_NUM_ERRORS GREATER 0)
    message ("Encountered build errors in Albany build. Exiting!")
  endif ()

  #
  # Run Albany tests
  #

  set (CTEST_TEST_TIMEOUT 1500)
  CTEST_TEST (
    BUILD "${CTEST_BINARY_DIRECTORY}/AlbBuild"
    RETURN_VALUE HAD_ERROR)

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Test RETURN_VALUE S_HAD_ERROR)

    if (S_HAD_ERROR)
      message ("Cannot submit Albany test results!")
    endif ()
  endif ()


endif ()
