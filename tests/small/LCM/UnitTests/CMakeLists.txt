# Disable these unit tests on the BGQ for now.
IF(NOT ALBANY_PARALLEL_ONLY)

# add the individual unit tests
IF(ALBANY_LCM AND LCM_TEST_EXES AND ALBANY_BGL)
  add_test(utLocalNonlinearSolver ${Albany_BINARY_DIR}/src/LCM/utLocalNonlinearSolver)
  add_test(utMiniSolvers ${Albany_BINARY_DIR}/src/LCM/utMiniSolvers)
  IF (ALBANY_ROL)
    add_test(utMiniSolversROL ${Albany_BINARY_DIR}/src/LCM/utMiniSolversROL)
  ENDIF()
  add_test(utSurfaceElement ${Albany_BINARY_DIR}/src/LCM/utSurfaceElement)
  add_test(utHeliumODEs ${Albany_BINARY_DIR}/src/LCM/utHeliumODEs)
  IF(ALBANY_LAME)
    add_test(utLameStress_elastic ${Albany_BINARY_DIR}/src/LCM/utLameStress_elastic)
  ENDIF()
ENDIF()

# create a custom target "make utest" that runs only the unit tests
add_custom_target(utest COMMAND ctest)

ENDIF()
