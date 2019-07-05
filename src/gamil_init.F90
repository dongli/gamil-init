program gamil_init

    use namelist_mod
    use model_grids
    use source_data
    use model_ic
    use model_bc
    use model_fc

    implicit none

    character(300) namelist_file

    if (command_argument_count() /= 1) then
        call print_usage
        stop
    end if
    call get_command_argument(1, namelist_file)

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
        call model_grids_write
    end if

    if (aero_data_file /= "N/A") then
        call source_data_start_aerosol(aero_data_file)
        call model_fc_interp_aerosol
        call source_data_end_aerosol
    end if

    if (sst_data_file /= "N/A") then
        call model_bc_init(model_ic_file_name())
        call model_bc_interp
        if (correct_sst_sice) call model_bc_correct
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
