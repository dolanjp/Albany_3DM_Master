
# Begin User inputs:
set (CTEST_SITE "cori.nersc.gov" ) # generally the output of hostname
set (CTEST_DASHBOARD_ROOT "$ENV{TEST_DIRECTORY}" ) # writable path
set (CTEST_SCRIPT_DIRECTORY "$ENV{SCRIPT_DIRECTORY}" ) # where the scripts live
set (CTEST_CMAKE_GENERATOR "Unix Makefiles" ) # What is your compilation apps ?
set (CTEST_CONFIGURATION  Release) # What type of build do you want ?

set (INITIAL_LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})

set (CTEST_PROJECT_NAME "Albany" )
set (CTEST_SOURCE_NAME repos)
set (CTEST_NAME "cori-gcc-${CTEST_BUILD_CONFIGURATION}")
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
set (CTEST_CMAKE_COMMAND "${PREFIX_DIR}/bin/cmake")
set (CTEST_COMMAND "${PREFIX_DIR}/bin/ctest -D ${CTEST_TEST_TYPE}")
set (CTEST_FLAGS "-j16")
set (CTEST_BUILD_FLAGS "-j16")

set (CTEST_DROP_METHOD "http")

#if (CTEST_DROP_METHOD STREQUAL "http")
#  set (CTEST_DROP_SITE "my.cdash.com")
#  set (CTEST_PROJECT_NAME "Albany")
#  set (CTEST_DROP_LOCATION "/submit.php?project=Albany")
#  set (CTEST_TRIGGER_SITE "")
#  set (CTEST_DROP_SITE_CDASH TRUE)
#endif ()

find_program (CTEST_GIT_COMMAND NAMES git)
find_program (CTEST_SVN_COMMAND NAMES svn)

set (Albany_REPOSITORY_LOCATION git@github.com:SNLComputation/Albany.git)
set (Trilinos_REPOSITORY_LOCATION git@github.com:trilinos/Trilinos.git)
set (cism-piscees_REPOSITORY_LOCATION  git@github.com:E3SM-Project/cism-piscees.git)


#set (BOOST_DIR /usr/common/software/boost/1.61/hsw/gnu) 
#set (NETCDF_DIR /opt/cray/pe/netcdf-hdf5parallel/4.4.0/GNU/5.1) 

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
      clone ${Albany_REPOSITORY_LOCATION} ${CTEST_SOURCE_DIRECTORY}/Albany
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

  #
  # Get cism-piscees
  #

  if (NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/cism-piscees")
    execute_process (COMMAND "${CTEST_GIT_COMMAND}"
      clone ${cism-piscees_REPOSITORY_LOCATION} -b ali_interface ${CTEST_SOURCE_DIRECTORY}/cism-piscees
      OUTPUT_VARIABLE _out
      ERROR_VARIABLE _err
      RESULT_VARIABLE HAD_ERROR)
    message(STATUS "out: ${_out}")
    message(STATUS "err: ${_err}")
    message(STATUS "res: ${HAD_ERROR}")
    if (HAD_ERROR)
      message(FATAL_ERROR "Cannot clone cism-piscees repository!")
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
  set_property (GLOBAL PROPERTY SubProject CoriTrilinos)
  set_property (GLOBAL PROPERTY Label CoriTrilinos)


  set (CONFIGURE_OPTIONS
    "-DCMAKE_INSTALL_PREFIX:PATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
    "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
    "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
    "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF"
    "-DTrilinos_ENABLE_Fortran:BOOL=ON"
    "-DTPL_ENABLE_SuperLU:BOOL=OFF"
    "-DAmesos2_ENABLE_KLU2:BOOL=ON"
    "-DTrilinos_ASSERT_MISSING_PACKAGES=OFF"
    "-DTrilinos_ENABLE_Teuchos:BOOL=ON"
    "-DHAVE_TEUCHOS_COMM_TIMERS=ON"
    "-DTrilinos_ENABLE_Epetra:BOOL=ON"
    "-DTrilinos_ENABLE_Tpetra:BOOL=ON"
    "-DTrilinos_ENABLE_EpetraExt:BOOL=ON"
    "-DTrilinos_ENABLE_Ifpack:BOOL=ON"
    "-DTrilinos_ENABLE_Ifpack2:BOOL=ON"
    "-DTrilinos_ENABLE_AztecOO:BOOL=ON"
    "-DTrilinos_ENABLE_Amesos:BOOL=ON"
    "-DTrilinos_ENABLE_Amesos2:BOOL=ON"
    "-DTrilinos_ENABLE_Anasazi:BOOL=ON"
    "-DTrilinos_ENABLE_Belos:BOOL=ON"
    "-DTrilinos_ENABLE_Phalanx:BOOL=ON"
    "-DTrilinos_ENABLE_Kokkos:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosCore:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosContainers:BOOL=ON"
    "-DTrilinos_ENABLE_MiniTensor:BOOL=ON"
    "-DTrilinos_ENABLE_ML:BOOL=ON"
    "-DTrilinos_ENABLE_MueLu:BOOL=ON"
    "-DMueLu_ENABLE_TESTS:BOOL=OFF"
    "-DTrilinos_ENABLE_Stratimikos:BOOL=ON"
    "-DTrilinos_ENABLE_Thyra:BOOL=ON"
    "-DTrilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON"
    "-DTrilinos_ENABLE_TrilinosCouplings:BOOL=ON"
    "-DTrilinos_ENABLE_Isorropia:BOOL=ON"
    "-DTPL_FIND_SHARED_LIBS:BOOL=OFF"
    "-DBUILD_SHARED_LIBS:BOOL=OFF"
    "-DTrilinos_LINK_SEARCH_START_STATIC:BOOL=ON"
    #
    "-DTrilinos_ENABLE_Kokkos:BOOL=ON"
    "-DTrilinos_ENABLE_KokkosCore:BOOL=ON"
    "-DPhalanx_KOKKOS_DEVICE_TYPE:STRING=SERIAL"
    "-DPhalanx_INDEX_SIZE_TYPE:STRING=INT"
    "-DPhalanx_SHOW_DEPRECATED_WARNINGS:BOOL=OFF"
    "-DKokkos_ENABLE_Serial:BOOL=ON"
    "-DKokkos_ENABLE_OpenMP:BOOL=OFF"
    "-DKokkos_ENABLE_Pthread:BOOL=OFF"
    #
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    #
    "-DTPL_ENABLE_BLAS:BOOL=ON"
    "-DBLAS_LIBRARY_NAMES:STRING="
    "-DLAPACK_LIBRARY_NAMES:STRING="
    "-DTPL_ENABLE_GLM:BOOL=OFF"
    "-DTPL_ENABLE_Matio:BOOL=OFF"
    "-DTrilinos_ENABLE_Zoltan:BOOL=ON"
    "-DTrilinos_ENABLE_Zoltan2:BOOL=ON"
    "-DZoltan_ENABLE_ULONG_IDS:BOOL=ON"
    "-DZOLTAN_BUILD_ZFDRIVE:BOOL=OFF"
    "-DZoltan2_ENABLE_Experimental:BOOL=ON"
    "-DTrilinos_ENABLE_FEI:BOOL=OFF"
    #
    "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
    "-DPiro_ENABLE_TESTS:BOOL=OFF"
    "-DTrilinos_ENABLE_EXAMPLES:BOOL=OFF"
    "-DTPL_ENABLE_MPI:BOOL=ON"
    #
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=OFF"
    #
    "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    "-DTpetra_INST_INT_LONG_LONG:BOOL=ON"
    "-DTpetra_INST_INT_INT:BOOL=OFF"
    #
    "-DMPI_USE_COMPILER_WRAPPERS:BOOL=OFF"
    "-DCMAKE_CXX_COMPILER:FILEPATH=CC"
    "-DCMAKE_C_COMPILER:FILEPATH=cc"
    "-DCMAKE_Fortran_COMPILER:FILEPATH=ftn"
    "-DTrilinos_ENABLE_Fortran=ON"
    "-DCMAKE_C_FLAGS:STRING=-O3 -DREDUCE_SCATTER_BUG"
    "-DCMAKE_CXX_FLAGS:STRING=-O3 -DREDUCE_SCATTER_BUG -DBOOST_NO_HASH"
    "-DCMAKE_EXE_LINKER_FLAGS:STRING='-static -Wl,-zmuldefs'"
    "-DTrilinos_ENABLE_SHADOW_WARNINGS=OFF"
    "-DTrilinos_ENABLE_CXX11=ON"
    "-DTPL_ENABLE_Pthread:BOOL=OFF"
    "-DTPL_ENABLE_BinUtils:BOOL=OFF"
    "-DTrilinos_ENABLE_ROL:BOOL=ON"
    #
    "-DMPI_EXEC:FILEPATH=srun"
    "-DMPI_EXEC_MAX_NUMPROCS:STRING=4"
    "-DMPI_EXEC_NUMPROCS_FLAG:STRING=-n"
    #
    "-DHAVE_INTREPID_KOKKOSCORE:BOOL=ON"
    "-DTrilinos_ENABLE_Intrepid:BOOL=ON"
    "-DTrilinos_ENABLE_Intrepid2:BOOL=ON"
    "-DTrilinos_ENABLE_NOX:BOOL=ON"
    "-DTrilinos_ENABLE_Rythmos:BOOL=ON"
    "-DTrilinos_ENABLE_OptiPack:BOOL=ON"
    "-DTrilinos_ENABLE_GlobiPack:BOOL=ON"
    "-DTrilinos_ENABLE_Stokhos:BOOL=OFF"
    "-DTrilinos_ENABLE_Piro:BOOL=ON"
    "-DTrilinos_ENABLE_Teko:BOOL=ON"
    "-DTrilinos_ENABLE_Pamgen:BOOL=ON"
    "-DAnasazi_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DAztecOO_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DBelos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DIfpack_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DNOX_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DPhalanx_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=OFF"
    "-DRythmos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DStratimikos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DThyra_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DTrilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON"
    "-DBoost_LIBRARY_DIRS:FILEPATH=$ENV{BOOST_DIR}/lib"
    "-DBoost_INCLUDE_DIRS:FILEPATH=$ENV{BOOST_DIR}/include"
    "-DBoostLib_LIBRARY_DIRS:FILEPATH=$ENV{BOOST_DIR}/lib"
    "-DBoostLib_INCLUDE_DIRS:FILEPATH=$ENV{BOOST_DIR}/include"
    "-DBoostAlbLib_LIBRARY_DIRS:FILEPATH=$ENV{BOOST_DIR}/lib"
    "-DBoostAlbLib_INCLUDE_DIRS:FILEPATH=$ENV{BOOST_DIR}/include"
    "-DTPL_ENABLE_Boost:BOOL=ON"
    "-DTPL_ENABLE_BoostLib:BOOL=ON"
    "-DTPL_ENABLE_BoostAlbLib:BOOL=ON"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DNetcdf_LIBRARY_DIRS:FILEPATH=$ENV{NETCDF_DIR}/lib"
    "-DTPL_Netcdf_INCLUDE_DIRS:PATH=$ENV{NETCDF_DIR}/include"
    "-DTrilinos_ENABLE_STKIO:BOOL=ON"
    "-DTrilinos_ENABLE_STKMesh:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASExodus:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASAprepro_lib:BOOL=ON"
    "-DTrilinos_ENABLE_Shards:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASIoss:BOOL=ON"
    "-DTrilinos_ENABLE_Pamgen:BOOL=ON"
    "-DCMAKE_SKIP_INSTALL_RPATH=TRUE"
    "-DTPL_Netcdf_PARALLEL:BOOL=ON"
  )

  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/TriBuild")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/TriBuild)
  endif ()

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/TriBuild"
    SOURCE "${CTEST_SOURCE_DIRECTORY}/Trilinos"
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

  set_property (GLOBAL PROPERTY SubProject CoriTrilinos)
  set_property (GLOBAL PROPERTY Label CoriTrilinos)
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

if (BUILD_ALB_FELIX)

  # Configure the Albany build 
  # Builds FELIX only. 
  #

  set_property (GLOBAL PROPERTY SubProject CoriAlbanyFELIX)
  set_property (GLOBAL PROPERTY Label CoriAlbanyFELIX)

  set (CONFIGURE_OPTIONS
    "-DALBANY_TRILINOS_DIR:FILEPATH=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
    "-DENABLE_DEMO_PDES=OFF" 
    "-DENABLE_ASCR=OFF" 
    "-DENABLE_LCM:BOOL=OFF"
    "-DENABLE_LANDICE:BOOL=ON"
    "-DENABLE_FAD_TYPE:STRING=SFad"
    "-DALBANY_SFAD_SIZE=16"
    "-DENABLE_MPAS_INTERFACE=ON"
    "-DENABLE_CISM_INTERFACE=ON"
    "-DENABLE_GPTL:BOOL=OFF"
    "-DAlbany_BUILD_STATIC_EXE:BOOL=ON"
    "-DINSTALL_ALBANY:BOOL=ON"
    "-DCMAKE_INSTALL_PREFIX:BOOL=${CTEST_BINARY_DIRECTORY}/AlbanyFELIXInstall"
    "-DCISM_INCLUDE_DIR:FILEPATH=${CTEST_SOURCE_DIRECTORY}/cism-piscees/libdycore"
    "-DCMAKE_EXE_LINKER_FLAGS:STRING='-static -Wl,-zmuldefs'"
    "-DENABLE_STOKHOS:BOOL=OFF"
    "-DENABLE_PARAMETERS_DEPEND_ON_SOLUTION:BOOL=ON"
    )
  
  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/AlbBuild")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/AlbBuild)
  endif ()

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/AlbBuild"
    SOURCE "${CTEST_SOURCE_DIRECTORY}/Albany"
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

  set_property (GLOBAL PROPERTY SubProject CoriAlbanyFELIX)
  set_property (GLOBAL PROPERTY Label CoriAlbanyFELIX)
  #set (CTEST_BUILD_TARGET all)
  set (CTEST_BUILD_TARGET install)

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

  #CTEST_TEST(
  #  BUILD "${CTEST_BINARY_DIRECTORY}/CoriAlbanyFELIX"
    #              PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
    #              INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
    #NUMBER_FAILED  TEST_NUM_FAILED
  #  RETURN_VALUE  HAD_ERROR
  #  )

  #if (CTEST_DO_SUBMIT)
  #  ctest_submit (PARTS Test
  #    RETURN_VALUE  S_HAD_ERROR
  #    )

  #  if (S_HAD_ERROR)
  #    message(FATAL_ERROR "Cannot submit Albany test results!")
  #  endif ()
  #endif ()

  #if (HAD_ERROR)
  #	message(FATAL_ERROR "Some Albany tests failed.")
  #endif ()

endif ()

if (BUILD_CISM_PISCEES)

  # Configure the CISM-Albany build 
  #

  set_property (GLOBAL PROPERTY SubProject CoriCismAlbany)
  set_property (GLOBAL PROPERTY Label CoriCismAlbany)


  set (CONFIGURE_OPTIONS
    "-DCISM_MPI_MODE:BOOL=ON"
    "-DCISM_SERIAL_MODE:BOOL=OFF"
    "-DCISM_BUILD_CISM_DRIVER:BOOL=ON"
    "-DCISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
    #
    "-DCISM_USE_TRILINOS:BOOL=ON"
    "-DCISM_TRILINOS_DIR=${CTEST_BINARY_DIRECTORY}/TrilinosInstall"
    "-DALBANY_LANDICE_DYCORE:BOOL=ON"
    "-DALBANY_LANDICE_CTEST:BOOL=ON"
    "-DCISM_ALBANY_DIR=${CTEST_BINARY_DIRECTORY}/AlbanyFELIXInstall"
    "-DCISM_NETCDF_DIR=$ENV{NETCDF_DIR}" 
    #
    "-DCMAKE_CXX_COMPILER=CC"
    "-DCMAKE_C_COMPILER=cc"
    "-DCMAKE_Fortran_COMPILER=ftn"
    #
    "-DCMAKE_CXX_FLAGS:STRING='-O2 -static -std=c++11'" 
    "-DCMAKE_EXE_LINKER_FLAGS:STRING='-static -Wl,-zmuldefs'"
    "-DBUILD_SHARED_LIBS:BOOL=OFF"
    "-DCISM_STATIC_LINKING:BOOL=ON"
    "-DCISM_Fortran_FLAGS='-ffree-line-length-none'" 
    "-DCISM_GNU:BOOL=ON"
  )
 
  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/CoriCismAlbany")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/CoriCismAlbany)
  endif ()

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/CoriCismAlbany"
    SOURCE "${CTEST_SOURCE_DIRECTORY}/cism-piscees"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE HAD_ERROR
    APPEND
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Configure
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message(FATAL_ERROR "Cannot submit CISM-Albany configure results!")
    endif ()
  endif ()

  if (HAD_ERROR)
    message(FATAL_ERROR "Cannot configure CISM-Albany build!")
  endif ()

  #
  # Build CISM-Albany
  #

  set (CTEST_TARGET all)

  MESSAGE("\nBuilding target: '${CTEST_TARGET}' ...\n")

  CTEST_BUILD(
    BUILD "${CTEST_BINARY_DIRECTORY}/CoriCismAlbany"
    RETURN_VALUE  HAD_ERROR
    NUMBER_ERRORS  LIBS_NUM_ERRORS
    APPEND
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Build
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message(FATAL_ERROR "Cannot submit CISM-Albany build results!")
    endif ()
  endif ()

  if (HAD_ERROR)
    message(FATAL_ERROR "Cannot build CISM-Albany!")
  endif ()

  if (LIBS_NUM_ERRORS GREATER 0)
    message(FATAL_ERROR "Encountered build errors in CISM-Albany build. Exiting!")
  endif ()

  #
  # Run CISM-Albany tests
  #

  #CTEST_TEST(
  #  BUILD "${CTEST_BINARY_DIRECTORY}/CoriAlbanyFELIX"
    #              PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
    #              INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
    #NUMBER_FAILED  TEST_NUM_FAILED
  #  RETURN_VALUE  HAD_ERROR
  #  )

  #if (CTEST_DO_SUBMIT)
  #  ctest_submit (PARTS Test
  #    RETURN_VALUE  S_HAD_ERROR
  #    )

  #  if (S_HAD_ERROR)
  #    message(FATAL_ERROR "Cannot submit Albany test results!")
  #  endif ()
  #endif ()

  #if (HAD_ERROR)
  #	message(FATAL_ERROR "Some Albany tests failed.")
  #endif ()

endif ()
