##*****************************************************************//
##    Albany 3.0:  Copyright 2016 Sandia Corporation               //
##    This Software is released under the BSD license detailed     //
##    in the file "license.txt" in the top-level Albany directory  //
##*****************************************************************//

SET(_3DM_DIR ${Albany_SOURCE_DIR}/src/3DM)

set(_3DM-problem-sources
  "${_3DM_DIR}/problems/_3DMThermoMechanics.cpp"
)
set(_3DM-problem-headers
  "${_3DM_DIR}/problems/_3DMThermoMechanics.hpp"
)

set(_3DM-evaluator-sources
  "${_3DM_DIR}/evaluators/_3DM_Time.cpp"
  "${_3DM_DIR}/evaluators/Energy_Dot.cpp"
  "${_3DM_DIR}/evaluators/Laser.cpp"
  "${_3DM_DIR}/evaluators/Laser_Source.cpp"
  "${_3DM_DIR}/evaluators/_3DM_Residual.cpp"
  "${_3DM_DIR}/evaluators/Phi1.cpp"
  "${_3DM_DIR}/evaluators/Phi2.cpp"
  "${_3DM_DIR}/evaluators/Psi1.cpp"
  "${_3DM_DIR}/evaluators/Psi2.cpp"
  "${_3DM_DIR}/evaluators/rho_Cp.cpp"
  "${_3DM_DIR}/evaluators/ThCond.cpp"
)
set(_3DM-evaluator-headers
  "${_3DM_DIR}/evaluators/_3DM_Time_Def.hpp"
  "${_3DM_DIR}/evaluators/_3DM_Time.hpp"
  "${_3DM_DIR}/evaluators/Energy_Dot_Def.hpp"
  "${_3DM_DIR}/evaluators/Energy_Dot.hpp"
  "${_3DM_DIR}/evaluators/Laser.hpp"
  "${_3DM_DIR}/evaluators/Laser_Source_Def.hpp"
  "${_3DM_DIR}/evaluators/Laser_Source.hpp"
  "${_3DM_DIR}/evaluators/_3DM_Residual_Def.hpp"
  "${_3DM_DIR}/evaluators/_3DM_Residual.hpp"
  "${_3DM_DIR}/evaluators/Phi1_Def.hpp"
  "${_3DM_DIR}/evaluators/Phi2_Def.hpp"
  "${_3DM_DIR}/evaluators/Phi1.hpp"
  "${_3DM_DIR}/evaluators/Phi2.hpp"
  "${_3DM_DIR}/evaluators/Psi1_Def.hpp"
  "${_3DM_DIR}/evaluators/Psi2_Def.hpp"
  "${_3DM_DIR}/evaluators/Psi1.hpp"
  "${_3DM_DIR}/evaluators/Psi2.hpp"
  "${_3DM_DIR}/evaluators/rho_Cp_Def.hpp"
  "${_3DM_DIR}/evaluators/rho_Cp.hpp"
  "${_3DM_DIR}/evaluators/ThCond_Def.hpp"
  "${_3DM_DIR}/evaluators/ThCond.hpp"
)

set(SOURCES
  ${_3DM-problem-sources}
  ${_3DM-evaluator-sources}
  )

set(HEADERS
  ${_3DM-problem-headers}
  ${_3DM-evaluator-headers}
 )

set(INCLUDE_DIRS
  ${Trilinos_INCLUDE_DIRS}
  ${Albany_SOURCE_DIR}/src
  ${Albany_SOURCE_DIR}/src/evaluators
  ${Albany_SOURCE_DIR}/src/problems 
  ${Albany_SOURCE_DIR}/src/disc/pumi
  ${Albany_SOURCE_DIR}/src/adapt
  ${_3DM_DIR}/problems
  ${_3DM_DIR}/evaluators
  ${ALBANY_PUMI_INCLUDE_DIRS} 
  ${Trilinos_TPL_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRS})
link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
add_library(_3DM ${Albany_LIBRARY_TYPE} ${SOURCES} ${HEADERS})
set_target_properties(_3DM PROPERTIES PUBLIC_HEADER "${HEADERS}")

IF ( INSTALL_ALBANY )
  install(TARGETS _3DM EXPORT albany-export
    LIBRARY DESTINATION "${LIB_INSTALL_DIR}"
    ARCHIVE DESTINATION "${LIB_INSTALL_DIR}"
    PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}")
ENDIF()
