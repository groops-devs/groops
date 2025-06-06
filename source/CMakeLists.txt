cmake_minimum_required(VERSION 3.10)

project(GROOPS LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_FLAGS "-std=legacy")
include_directories(${PROJECT_SOURCE_DIR})
include(sourcesCXX.txt)
include(sourcesF77.txt)
include(sourcesF90.txt)
add_library(groopscore OBJECT ${SOURCES})

# =========================================

# Libraries
# ---------
# stdc++fs  required C++14 std::experimental::filesystem or C++17 std::filesystem
# EXPAT     required Stream-oriented XML parser library (https://libexpat.github.io/)
# BLAS      required Basic Linear Algebra Subprograms (http://www.netlib.org/blas/)
# LAPACK    required Linear Algebra PACKage (http://www.netlib.org/lapack/)
# ERFA      optional Essential Routines for Fundamental Astronomy (https://github.com/liberfa)
# Z         optional File compression (https://www.zlib.net)
# NETCDF    optional Network Common Data Form (https://www.unidata.ucar.edu/software/netcdf/)

# External Source Files
# ---------------------
# HWM14     Horizontal Wind Model 2014 (https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/)
# NRLMSIS   NRLMSIS 2.0 thermospheric model (https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0/)
# JB2008    JB2008 thermospheric model (http://sol.spacenvironment.net/jb2008/)
# IGRF      International Geomagnetic Reference Field (https://doi.org/10.1186/s40623-015-0228-9)
# IERS      International Earth Rotation and Reference Systems Service (IERS) Conventions software collection (https://iers-conventions.obspm.fr/)

find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(EXPAT  REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

set(BASE_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXPAT_LIBRARIES} stdc++fs)

find_library(LIB_ERFA erfa)
if(LIB_ERFA AND ((NOT ${DISABLE_ERFA}) OR (NOT DEFINED DISABLE_ERFA)))
  find_path(ERFA_INCLUDE_DIR NAMES erfa.h)
  include_directories(${ERFA_INCLUDE_DIR})
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_ERFA})
else()
  add_definitions(-DGROOPS_DISABLE_ERFA)
  message(WARNING "ERFA library *NOT* found (https://github.com/liberfa). GROOPS is not able to compute Earth rotation based on IERS EOP files.")
endif()

find_library(LIB_Z z)
if(LIB_Z AND ((NOT ${DISABLE_Z}) OR (NOT DEFINED DISABLE_Z)))
  find_path(ZLIB_INCLUDE_DIR NAMES zlib.h)
  include_directories(${ZLIB_INCLUDE_DIR})
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_Z})
else()
  add_definitions(-DGROOPS_DISABLE_Z)
  message(WARNING "Z library *NOT* found (https://www.zlib.net). GROOPS is not able to read/write compressed *.gz files.")
endif()

find_library(LIB_NETCDF netcdf)
if(LIB_NETCDF AND ((NOT ${DISABLE_NETCDF}) OR (NOT DEFINED DISABLE_NETCDF)))
  find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)
  include_directories(${NETCDF_INCLUDE_DIR})
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_NETCDF})
else()
  add_definitions(-DGROOPS_DISABLE_NETCDF)
  message(WARNING "netCDF library *NOT* found (https://www.unidata.ucar.edu/software/netcdf). GROOPS is not able to convert netCDF *.grd files.")
endif()

if(${DISABLE_HWM14})
  message(WARNING "HWM14 wind model will *NOT* be compiled.")
  add_definitions(-DGROOPS_DISABLE_HWM14)
endif()

if(${DISABLE_NRLMSIS})
  message(WARNING "NRLMSIS thermospheric model will *NOT* be compiled.")
  add_definitions(-DGROOPS_DISABLE_NRLMSIS)
endif()

if(${DISABLE_JB2008})
  message(WARNING "JB2008 thermospheric model will *NOT* be compiled.")
  add_definitions(-DGROOPS_DISABLE_JB2008)
endif()

if(${DISABLE_IGRF})
  message(WARNING "IGRF magnetic field model will *NOT* be compiled.")
  add_definitions(-DGROOPS_DISABLE_IGRF)
endif()

if(${DISABLE_IERS})
  message(WARNING "IERS software routines will *NOT* be compiled.")
  add_definitions(-DGROOPS_DISABLE_IERS)
endif()

# =========================================

add_executable(groops ${PROJECT_SOURCE_DIR}/parallel/parallelSingle.cpp $<TARGET_OBJECTS:groopscore>)
target_link_libraries(groops ${BASE_LIBRARIES})

install(TARGETS groops DESTINATION bin)

# =========================================

find_package(MPI COMPONENTS CXX)
if(MPI_FOUND)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  add_executable(groopsMPI ${PROJECT_SOURCE_DIR}/parallel/parallelCluster.cpp $<TARGET_OBJECTS:groopscore>)

  target_link_libraries(groopsMPI ${BASE_LIBRARIES} ${MPI_CXX_LIBRARIES})
  if(MPI_COMPILE_FLAGS)
    set_target_properties(groopsMPI PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(groopsMPI PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()
  install(TARGETS groopsMPI DESTINATION bin)
else()
  message(WARNING "A GROOPS executable for MPI parallel computing will *NOT* be created.")
endif()

# =========================================
