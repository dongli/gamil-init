module model_bc

    use utils
    use model_grids

    implicit none

    private

    public model_bc_init
    public model_bc_interp
    public model_bc_correct

    character(300) sst_data_file
    character(300) ice_data_file
    character(300) grid_file
    character(300) tmp_file
    character(300) clm_file

contains

    subroutine model_bc_init(sst_data_file_, ice_data_file_, grid_file_)

        character(*), intent(in) :: sst_data_file_
        character(*), intent(in) :: ice_data_file_
        character(*), intent(in) :: grid_file_

        sst_data_file = sst_data_file_
        ice_data_file = ice_data_file_
        grid_file = grid_file_
        tmp_file = "sstice.tmp.nc"
        clm_file = "sstice.clm.nc"

    end subroutine model_bc_init

    subroutine model_bc_interp

        character(30), parameter :: sub_name = "model_bc_interp"

        call notice(sub_name, "Interpolation SST and sea ice data onto model grids")

        call pcmdi_regrid(sst_data_file, ice_data_file, grid_file, tmp_file)

    end subroutine model_bc_interp

    subroutine model_bc_correct

        character(300) file_name
        integer ret

        character(30), parameter :: sub_name = "model_bc_correct"

        file_name = model_bc_file_name()

        call notice(sub_name, "Correct month mean")

        call pcmdi_bcgen(tmp_file, clm_file, file_name, &
                         1, 1870, 12, 2015, 1, 1870, 12, 2015, &
                         1, 1870, 12, 2015, 1, 1870, 12, 2015)

        call notice(sub_name, "File "//trim(file_name)//" has been generated")

        ! Remove temporal files.
        call execute_command_line("rm -rf "//trim(tmp_file), exitstat=ret)
        if (ret /= 0) then
            call report_error(sub_name, "Failed to remove "//trim(tmp_file))
        end if
        call execute_command_line("rm -rf "//trim(clm_file), exitstat=ret)
        if (ret /= 0) then
            call report_error(sub_name, "Failed to remove "//trim(clm_file))
        end if

    end subroutine model_bc_correct

    character(300) function model_bc_file_name()

        model_bc_file_name = "bc.gamil."// &
            trim(to_string(num_model_lon))//"x"// &
            trim(to_string(num_model_lat))//".nc"

    end function model_bc_file_name

end module model_bc
