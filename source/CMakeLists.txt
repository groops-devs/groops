cmake_minimum_required(VERSION 3.7)

project(GROOPS LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
include_directories(${PROJECT_SOURCE_DIR})
set(SOURCES)
include(sourcesCXX.txt)
include(sourcesF77.txt)
include(sourcesF90.txt)

# Libraries
# ---------
# stdc++fs  required C++14 std::experimental::filesystem or C++17 std::filesystem
# EXPAT     required Stream-oriented XML parser library (https://libexpat.github.io/)
# BLAS      required Basic Linear Algebra Subprograms (http://www.netlib.org/blas/)
# LAPACK    required Linear Algebra PACKage (http://www.netlib.org/lapack/)
# ERFA      optional Essential Routines for Fundamental Astronomy (https://github.com/liberfa)
# Z         optional File compression (https://www.zlib.net)
# NETCDF    optional Network Common Data Form (https://www.unidata.ucar.edu/software/netcdf/)
# =========================================

find_package(EXPAT  REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
set(BASE_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXPAT_LIBRARIES} stdc++fs )

find_library(LIB_ERFA erfa)
if(LIB_ERFA)
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_ERFA})
else()
  add_definitions(-DNOLIB_ERFA)
  message(WARNING "ERFA library *NOT* found (https://github.com/liberfa). GROOPS is not able to compute Earth rotation based on IERS EOP files.")
endif()

find_library(LIB_Z z)
if(LIB_Z)
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_Z})
else()
  add_definitions(-DNOLIB_Z)
  message(WARNING "Z library *NOT* found (https://www.zlib.net). GROOPS is not able to read/write compressed *.gz files.")
endif()

find_library(LIB_NETCDF netcdf)
if(LIB_NETCDF)
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_NETCDF})
else()
  add_definitions(-DNOLIB_NETCDF)
  message(WARNING "netCDF library *NOT* found (https://www.unidata.ucar.edu/software/netcdf). GROOPS is not able to convert netCDF *.grd files.")
endif()

# =========================================

add_executable(groops ${SOURCES} ${PROJECT_SOURCE_DIR}/parallel/parallelSingle.cpp)
target_link_libraries(groops ${BASE_LIBRARIES})
set_target_properties(groops PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../bin")

# =========================================

find_package(MPI COMPONENTS CXX)
if(MPI_FOUND)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  add_executable(groopsMPI ${SOURCES} ${PROJECT_SOURCE_DIR}/parallel/parallelCluster.cpp)
  set_target_properties(groopsMPI PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../bin")
  target_link_libraries(groopsMPI ${BASE_LIBRARIES} ${MPI_CXX_LIBRARIES})
  if(MPI_COMPILE_FLAGS)
    set_target_properties(groopsMPI PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(groopsMPI PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()
else()
  message(WARNING "A GROOPS executable for MPI parallel computing will *NOT* be created.")
endif()

# =========================================
