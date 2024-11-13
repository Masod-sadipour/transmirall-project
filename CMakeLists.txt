# CMake minimum version
cmake_minimum_required(VERSION 3.10)

# Project name
project(TransmuralLaplaceSolver)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen manually
include_directories("/usr/include/eigen3")

# Enable verbose makefile for debugging (optional)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set the installation paths for PETSc and libMesh
set(PETSC_DIR /nas/longleaf/home/masods/transmirall/external_libs/petsc_install)
set(PETSC_ARCH arch-linux2-c-opt)

# Find Boost
find_package(Boost 1.66 REQUIRED COMPONENTS system filesystem iostreams)
include_directories(${Boost_INCLUDE_DIRS})


# Set the libMesh installation path
set(LIBMESH_DIR /nas/longleaf/home/masods/transmirall/external_libs/libmesh_install)
set(CMAKE_PREFIX_PATH ${LIBMESH_DIR})
find_package(PkgConfig REQUIRED)
pkg_check_modules(LibMesh REQUIRED IMPORTED_TARGET libmesh)

# Include directories for header files
include_directories(${LibMesh_INCLUDE_DIRS})
link_directories(${LibMesh_LIBRARY_DIRS})

include_directories(${CMAKE_SOURCE_DIR}/BoundaryConditions)


# Source files
set(SOURCES
    BoundaryConditions/BCs.cpp
    src/main.cpp
    src/Poisson.cpp
)

# Add executable
add_executable(transmural_solver ${SOURCES})

# Include Eigen directories explicitly
target_include_directories(transmural_solver PRIVATE ${EIGEN3_INCLUDE_DIR})

# Link libraries
target_link_libraries(transmural_solver
    ${LibMesh_LIBRARIES}
    ${Boost_LIBRARIES}
    libmesh_getpot
)

# Set C++ compiler flags for better error reporting
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
