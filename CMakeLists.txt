# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME demo_3D_P2_cylinder)
project(${PROJECT_NAME})

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get FEniCS Plasticity and DOLFIN configuration data
#find_package(fenics-solid-mechanics)
find_package(DOLFIN)
include(${DOLFIN_USE_FILE})

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# C++11 compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++11 HAVE_STD_CPP11)
if (HAVE_STD_CPP11)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_STD_CPP0x)
  if (HAVE_STD_CPP0x)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Include directories
include_directories(${FENICS_SOLID_MECHANICS_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_INCLUDE_DIRS}
  ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/src")
add_subdirectory(src/)

# Executable
add_executable(${PROJECT_NAME} main.cpp Plas3D.h Plas3D.cpp)

# Need to get VTK config because VTK uses advanced VTK features which
# mean it's not enough to just link to the DOLFIN target. See
# http://www.vtk.org/pipermail/vtk-developers/2013-October/014402.html
find_package(VTK HINTS ${VTK_DIR} $ENV{VTK_DIR} NO_MODULE QUIET)

# Target libraries
#target_link_libraries(${PROJECT_NAME} ${FENICS_SOLID_MECHANICS_LIBRARIES}
#  ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} fenics-solid-mechanics)

target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} fenics-solid-mechanics)


  set (COPY_FILES cylinder.xml)
foreach(loop_var ${COPY_FILES})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${loop_var} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach(loop_var)
