program gamil_init

    use model_grids
    use source_data
    use model_ic
    use model_bc
    use model_fc

    implicit none

    character(300) :: model_grid_file = "N/A"
    character(300) :: uvtq_data_file  = "N/A"
    character(300) :: ozone_data_file = "N/A"
    character(300) :: aero_data_file = "N/A"
    character(300) :: topo_data_file  = "N/A"
    character(300) :: sst_data_file = "N/A"
    character(300) :: ice_data_file = "N/A"

    character(300) namelist_file

    namelist /gamil_init_config/ &
        num_model_lon, num_model_lat,    &
        model_grid_type, uvtq_data_type, &
        model_grid_file, uvtq_data_file, &
        ozone_data_file, topo_data_file, &
        aero_data_file, sst_data_file,   &
        ice_data_file, &
        ! Some control parameters from modules
        model_ic_allow_extrap, &
        latmesh_B

    if (command_argument_count() /= 1) then
        call print_usage
        stop
    end if
    call get_command_argument(1,namelist_file)

    open(10, file=namelist_file, status="old")
    read(10, nml=gamil_init_config)
    close(10)

    if (model_grid_file /= "N/A") then
        call model_grids_read(model_grid_file)
    else
        call calc_model_grids(num_model_lon, num_model_lat, latmesh_B)
    end if

    if (uvtq_data_file /= "N/A") then
        call source_data_read_uvtq(uvtq_data_file)
        call model_ic_calc_topo(topo_data_file)
        call model_ic_interp
        call model_ic_write
    end if

    if (aero_data_file /= "N/A") then
        call source_data_start_aerosol(aero_data_file)
        call model_fc_interp_aerosol
        call source_data_end_aerosol
    end if

    if (sst_data_file /= "N/A") then
        call model_bc_init(sst_data_file, ice_data_file, model_ic_file_name())
        call model_bc_interp
        call model_bc_correct
    end if

contains

    subroutine print_usage

        write(*, "('------------------------------------------------------')")
        write(*, "('---------------------- GAMIL-INIT --------------------')")
        write(*, "('------------------------------------------------------')")
        write(*, "('Usage:')")
        write(*, *)
        write(*, "('  ./gamil-init <namelist>')")

    end subroutine print_usage

end program gamil_init
