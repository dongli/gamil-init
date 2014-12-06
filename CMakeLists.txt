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
    "${PROJECT_SOURCE_DIR}/src/model_bc.F90"
    "${PROJECT_SOURCE_DIR}/src/model_fc.F90"
    "${PROJECT_SOURCE_DIR}/src/source_data.F90"
    "${PROJECT_SOURCE_DIR}/src/utils.F90"
    "${PROJECT_SOURCE_DIR}/src/variable.F90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/ao.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/ao_i.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/area_ave.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/binf2c.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/cell_area.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/err_exit.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/get_curr_timestamp.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/interp_driver.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/lininterp.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/map_i.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/max_ovr.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/precision.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/regrid/regrid.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/bcgen.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/calcclim.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/calcfull.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/driver.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/output_dateinfo.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/prec.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/setup_outfile.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/solver.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/types.f90"
    "${PROJECT_SOURCE_DIR}/src/pcmdi/bcgen/wrap_nf.f90"
)

include_directories (
    "${PROJECT_SOURCE_DIR}/src"
)

add_executable (gamil-init ${sources})
target_link_libraries (gamil-init ${NETCDF_LIBRARIES})
