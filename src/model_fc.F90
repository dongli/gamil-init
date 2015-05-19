module model_fc

    use utils
    use source_data
    use variable
    use model_grids
    use model_gears
    use io_manager

    implicit none

    type(var_list) model_aero_list

contains

    subroutine model_fc_interp_aerosol

        integer dims(3), num_sub_aero
        integer time, aero, i, j, k

        integer, parameter :: num_aero = 10
        character(30) data_aero_names(num_aero)
        character(10) model_aero_names(num_aero)
        real(8) int_coef(num_aero)

        real(8), pointer :: model_ps(:,:)
        integer, pointer :: date
        class(var), pointer :: ptr1, ptr2
        real(8), pointer :: tmp1(:,:,:)
        real(8), allocatable, target :: tmp2(:,:,:)

        character(100) output_file_name
        integer file_idx

        character(50), parameter :: sub_name = "model_fc_interp_aerosol"

        data_aero_names ( 1) = "SO4"
        model_aero_names( 1) = "MSUL_V"
        int_coef        ( 1) = 1.0d0
        data_aero_names ( 2) = "SSLT01+SSLT02+SSLT03+SSLT04"
        model_aero_names( 2) = "MSSLT_V"
        int_coef        ( 2) = 1.0d0
        data_aero_names ( 3) = "DST01"
        model_aero_names( 3) = "MDUST1_V"
        int_coef        ( 3) = 1.0d0
        data_aero_names ( 4) = "DST02"
        model_aero_names( 4) = "MDUST2_V"
        int_coef        ( 4) = 1.0d0
        data_aero_names ( 5) = "DST03"
        model_aero_names( 5) = "MDUST3_V"
        int_coef        ( 5) = 1.0d0
        data_aero_names ( 6) = "DST04"
        model_aero_names( 6) = "MDUST4_V"
        int_coef        ( 6) = 1.0d0
        data_aero_names ( 7) = "OC1"
        model_aero_names( 7) = "MOCPHO_V"
        int_coef        ( 7) = 1.0d0
        data_aero_names ( 8) = "CB1"
        model_aero_names( 8) = "MBCPHO_V"
        int_coef        ( 8) = 1.0d0
        data_aero_names ( 9) = "OC2"
        model_aero_names( 9) = "MOCPHI_V"
        int_coef        ( 9) = 1.0d0
        data_aero_names (10) = "CB2"
        model_aero_names(10) = "MBCPHI_V"
        int_coef        (10) = 1.0d0

        output_file_name = "fc.aero."// &
            trim(to_string(num_model_lon))//"x"// &
            trim(to_string(num_model_lat))//".nc"

        call io_manager_create_file(output_file_name, file_idx)

        ! Define dimensions.
        ptr1 => data_aero_list%get_var("time")
        call io_manager_def_dim(file_idx, "time", &
            data_type=ptr1%get_data_type(), units=ptr1%get_units())
        call model_grids_add_2d_dims(file_idx)
        call io_manager_add_dim(file_idx, data_aero_list%get_var("lev"))
        call io_manager_add_dim(file_idx, data_aero_list%get_var("ilev"))

        call io_manager_def_var(file_idx, data_aero_list%get_var("P0"))
        call io_manager_put_var(file_idx, data_aero_list%get_var("P0"))

        ! Add 1D auxiliary variables.
        ptr1 => data_aero_list%get_head()
        do i = 1, data_aero_list%get_num_var()
            if (ptr1%get_num_dim() /= 1 .or. &
                any(["lon ","lat ","lev ","ilev","time"] == ptr1%get_name())) then
                ptr1 => ptr1%next
                cycle
            end if
            if (any(["date   ","datesec"] == ptr1%get_name())) then
                call io_manager_def_var(file_idx, ptr1, ["time"])
            else if (any(["hyam","hybm"] == ptr1%get_name())) then
                call io_manager_def_var(file_idx, ptr1, ["lev"])
            else if (any(["hyai","hybi"] == ptr1%get_name())) then
                call io_manager_def_var(file_idx, ptr1, ["ilev"])
            end if
            call io_manager_put_var(file_idx, ptr1)
            ptr1 => ptr1%next
        end do

        ! Add aerosol data.
        dims = [num_model_lon,num_model_lat,num_aero_lev]
        do i = 1, num_aero
            num_sub_aero = count_substr(data_aero_names(i), "+")
            if (num_sub_aero == 0) then
                ptr1 => data_aero_list%get_var(data_aero_names(i))
            else
                ptr1 => data_aero_list%get_var(get_substr(data_aero_names(i), "+", 1))
            end if
            call model_aero_list%append(model_aero_names(i), &
                ptr1%get_long_name(), ptr1%get_units(), dims)
            call io_manager_def_var(file_idx, model_aero_list%get_tail(), ["lon ","lat ","lev ","time"])
        end do

        ! Add aerosol data surface pressure.
        call model_aero_list%append("PS", "surface pressure", "Pa", dims(1:2))
        call model_aero_list%get_tail_values(model_ps)
        call io_manager_def_var(file_idx, model_aero_list%get_tail(), ["lon ","lat ","time"])

        allocate(tmp2(num_aero_lon,num_aero_lat,num_aero_lev))

        do time = 1, num_aero_time
            ! ------------------------------------------------------------------
            call io_manager_put_var(file_idx, "time", aero_time(time), rec=time)
            call io_manager_put_var(file_idx, "date", aero_date(time), rec=time)
            call io_manager_put_var(file_idx, "datesec", aero_datesec(time), rec=time)
            ! ------------------------------------------------------------------
            call notice(sub_name, "Interpolate aerosol data on "// &
                trim(integer_to_string(aero_date(time))))
            ! ------------------------------------------------------------------
            call source_data_read_aerosol
            ! ------------------------------------------------------------------
            ! interpolate surface pressure onto model grids
            call model_gears_interp_h(aero_lon, aero_lat, aero_ps, &
                model_lon, model_lat, model_ps, 1)
            call io_manager_put_var(file_idx, model_aero_list%get_var("PS"), rec=time)
            ! ------------------------------------------------------------------
            ! interpolate aerosol on model grids
            do aero = 1, num_aero
                ! ==============================================================
                ! get source aerosol data
                num_sub_aero = count_substr(data_aero_names(aero), "+")+1
                if (num_sub_aero == 1) then
                    ptr1 => data_aero_list%get_var(data_aero_names(aero))
                    select type (ptr1)
                    type is (var3d_d)
                        tmp1 => ptr1%get_values()
                    end select
                else
                    tmp2 = 0.0d0
                    do i = 1, num_sub_aero
                        ptr1 => data_aero_list%get_var( &
                            get_substr(data_aero_names(aero), "+", i))
                        select type (ptr1)
                        type is (var3d_d)
                            tmp1 => ptr1%get_values()
                        end select
                        tmp2 = tmp2+tmp1
                    end do
                    tmp1 => tmp2
                end if
                ! ==============================================================
                ! Integrate source aerosol data.
                do j = 1, num_aero_lat
                    do i = 1, num_aero_lon
                        do k = 1, num_aero_lev
                            tmp2(i,j,k) = tmp1(i,j,k)*aero_dp(i,j,k)/g*int_coef(aero)
                        end do
                        do k = num_aero_lev-1, 1, -1
                            tmp2(i,j,k) = tmp2(i,j,k)+tmp2(i,j,k+1)
                        end do
                    end do
                end do
                ! ==============================================================
                ! Get model aerosol data pointer.
                ptr2 => model_aero_list%get_var(model_aero_names(aero))
                select type (ptr2)
                type is (var3d_d)
                    tmp1 => ptr2%get_values()
                end select
                ! ==============================================================
                ! Interpolate (only horizontal)
                do k = 1, num_aero_lev
                    call model_gears_interp_h(aero_lon, aero_lat, tmp2(:,:,k), &
                        model_lon, model_lat, tmp1(:,:,k), 1)
                end do
                ! ==============================================================
                call io_manager_put_var(file_idx, ptr2, rec=time)
            end do
        end do

        call io_manager_close_file(file_idx)

        call notice(sub_name, "File "//trim(output_file_name)//" is generated!")

    end subroutine model_fc_interp_aerosol

end module model_fc
