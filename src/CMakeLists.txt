# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Include directories required to build
include_directories(SYSTEM ${DOLFIN_INCLUDE_DIRS}
  ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

add_definitions(-DFENICSSOLID_GIT_COMMIT_HASH="${GIT_COMMIT_HASH}")

# Source files to build
add_library(fenics-solid-mechanics
    AnbaConstitutiveUpdate.cpp  FenicsSolidMechanics.h  QuadratureData.cpp       utils.h
    AnbaConstitutiveUpdate.h    HistoryData.cpp         QuadratureData.h        VonMises.cpp            HistoryData.h
    QuadratureFunction.cpp  UFLQuadratureFunction.cpp         VonMises.h          ConstitutiveUpdate.cpp
    PlasticityModel.cpp     QuadratureFunction.h    UFLQuadratureFunction.h     ConstitutiveUpdate.h        PlasticityModel.h
    ReturnMapping.cpp       UFLReturnMapping.cpp    DruckerPrager.cpp           PlasticityProblem.cpp   ReturnMapping.h
    UFLReturnMapping.h      DruckerPrager.h             PlasticityProblem.h     StateUpdate.h           utils.cpp)

# Link to DOLFIN library
target_link_libraries(fenics-solid-mechanics ${DOLFIN_LIBRARIES})

# BLAS/LAPACK on OSX. Should get from DOLFIN CMake config file.
#if (APPLE)
#   set(CMAKE_SHARED_LINKER_FLAGS "-framework vecLib")
#endif()

#------------------------------------------------------------------------------
# Install header files
#install(DIRECTORY ./ DESTINATION include/fenics-solid-mechanics
#  FILES_MATCHING PATTERN "*.h")

#install(TARGETS fenics-solid-mechanics EXPORT FEniCSSolidMechanicsTargets
#  RUNTIME DESTINATION ${FENICS_SOLID_MECHANICS_LIB_DIR} COMPONENT RuntimeExecutables
#  LIBRARY DESTINATION ${FENICS_SOLID_MECHANICS_LIB_DIR} COMPONENT RuntimeLibraries
#  ARCHIVE DESTINATION ${FENICS_SOLID_MECHANICS_LIB_DIR} COMPONENT Development)

#install(EXPORT FEniCSSolidMechanicsTargets
#  DESTINATION ${FENICS_SOLID_MECHANICS_SHARE_DIR}/cmake
#  COMPONENT Development)
#------------------------------------------------------------------------------
# Generate CMake config files (dolfin-config{,-version}.cmake)

#configure_file(${FENICS_SOLID_MECHANICS_CMAKE_DIR}/fenics-solid-mechanics-config.cmake.in
#  ${CMAKE_BINARY_DIR}/fenics-solid-mechanics/fenics-solid-mechanics-config.cmake @ONLY)

#install(FILES
#  ${CMAKE_BINARY_DIR}/fenics-solid-mechanics/fenics-solid-mechanics-config.cmake
#  DESTINATION ${FENICS_SOLID_MECHANICS_SHARE_DIR}/cmake
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
