cmake_minimum_required (VERSION 2.8)

project (gamil-init Fortran)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
find_package (NETCDF REQUIRED COMPONENTS Fortran)
include_directories (${NETCDF_INCLUDE_DIRS})
message (STATUS ${NETCDF_INCLUDE_DIRS})

set (sources
    "${PROJECT_SOURCE_DIR}/src/constants.F90"
    "${PROJECT_SOURCE_DIR}/src/gamil_init.F90"
    "${PROJECT_SOURCE_DIR}/src/interp.F90"
    "${PROJECT_SOURCE_DIR}/src/io_manager.F90"
    "${PROJECT_SOURCE_DIR}/src/model_gears.F90"
    "${PROJECT_SOURCE_DIR}/src/model_grids.F90"
    "${PROJECT_SOURCE_DIR}/src/model_ic.F90"
    "${PROJECT_SOURCE_DIR}/src/model_fc.F90"
    "${PROJECT_SOURCE_DIR}/src/source_data.F90"
    "${PROJECT_SOURCE_DIR}/src/utils.F90"
    "${PROJECT_SOURCE_DIR}/src/variable.F90"
)

include_directories (
    "${PROJECT_SOURCE_DIR}/src"
)

add_executable (gamil-init ${sources})
target_link_libraries (gamil-init ${NETCDF_LIBRARIES})
