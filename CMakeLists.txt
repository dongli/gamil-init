cmake_minimum_required (VERSION 2.8)

project (gamil-init Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

add_subdirectory(lib/datetime)

set (sources
    src/constants.F90
    src/gamil_init.F90
    src/interp.F90
    src/io_manager.F90
    src/model_gears.F90
    src/model_grids.F90
    src/model_ic.F90
    src/model_bc.F90
    src/model_fc.F90
    src/source_data.F90
    src/utils.F90
    src/namelist_mod.F90
    src/variable.F90
    src/pcmdi/regrid/ao.f90
    src/pcmdi/regrid/ao_i.f90
    src/pcmdi/regrid/area_ave.f90
    src/pcmdi/regrid/binf2c.f90
    src/pcmdi/regrid/cell_area.f90
    src/pcmdi/regrid/err_exit.f90
    src/pcmdi/regrid/get_curr_timestamp.f90
    src/pcmdi/regrid/interp_driver.f90
    src/pcmdi/regrid/lininterp.f90
    src/pcmdi/regrid/map_i.f90
    src/pcmdi/regrid/max_ovr.f90
    src/pcmdi/regrid/precision.f90
    src/pcmdi/regrid/regrid.f90
    src/pcmdi/bcgen/bcgen.f90
    src/pcmdi/bcgen/calcclim.f90
    src/pcmdi/bcgen/calcfull.f90
    src/pcmdi/bcgen/driver.f90
    src/pcmdi/bcgen/output_dateinfo.f90
    src/pcmdi/bcgen/prec.f90
    src/pcmdi/bcgen/setup_outfile.f90
    src/pcmdi/bcgen/solver.f90
    src/pcmdi/bcgen/types.f90
    src/pcmdi/bcgen/wrap_nf.f90
)

add_executable (gamil-init ${sources})
target_link_libraries (gamil-init netcdff fortran_datetime)
