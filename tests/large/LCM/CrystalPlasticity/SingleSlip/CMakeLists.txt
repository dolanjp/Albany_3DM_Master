add_subdirectory(SingleSlip)
add_subdirectory(SingleSlipHard)
add_subdirectory(SingleSlipSaturation)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SingleSlip.g
               ${CMAKE_CURRENT_BINARY_DIR}/SingleSlip.g COPYONLY)
