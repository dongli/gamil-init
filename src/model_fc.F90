module model_fc

    use utils
    use source_data
    use variable
    use model_grids
    use model_gears
    use io_manager

    implicit none

    type(var_list) model_aero_list
    type(var1d_d) model_aero_lev

contains

    subroutine model_fc_interp_aerosol

        integer dims(3), num_sub_aero, num_lev
        integer time, aero, i, j, k

        integer, parameter :: num_aero = 10
        character(30) data_aero_names(num_aero)
        character(10) model_aero_names(num_aero)
        real(8) int_coef(num_aero)

        real(8), pointer :: model_ps(:,:)
        integer, pointer :: date
        real(8), allocatable :: model_p(:,:,:), tmp_p(:,:,:)
        class(var), pointer :: ptr1, ptr2
        real(8), pointer :: tmp1(:,:,:), tmp2(:,:,:)
        real(8), allocatable, target :: tmp3(:,:,:), tmp4(:,:,:)

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

        ! choose the vertical levels of model aerosol data
        num_lev = num_model_lev+1
        call model_aero_lev%init("lev", "level (sigma vertical cooordinate)", "1", [num_lev])
        call model_aero_lev%set_values(model_lev_bnds)

        ! Add aerosol data.
        dims = [num_model_lon,num_model_lat,num_lev]
        do i = 1, num_aero
            num_sub_aero = count_substr(data_aero_names(i), "+")
            if (num_sub_aero == 0) then
                ptr1 => data_aero_list%get_var(data_aero_names(i))
            else
                ptr1 => data_aero_list%get_var(get_substr(data_aero_names(i), "+", 1))
            end if
            call model_aero_list%append(model_aero_names(i), &
                ptr1%get_long_name(), ptr1%get_units(), dims)
        end do

        ! Add aerosol data surface pressure.
        call model_aero_list%append("PS", "surface pressure", "Pa", dims(1:2))
        call model_aero_list%get_tail_values(model_ps)

        allocate(model_p(num_model_lon,num_model_lat,num_lev))
        allocate(tmp_p(num_model_lon,num_model_lat,num_aero_lev))
        allocate(tmp3(num_model_lon,num_model_lat,num_aero_lev))
        allocate(tmp4(num_aero_lon,num_aero_lat,num_aero_lev))

        output_file_name = "fc.aero."// &
            trim(to_string(num_model_lon))//"x"// &
            trim(to_string(num_model_lat))//"x"// &
            trim(to_string(num_model_lev))//".nc"

        call io_manager_create_file(output_file_name, file_idx)
        call io_manager_def_dim(file_idx, "time", "integer", &
            units="days since 1849-01-01")
        call model_grids_add_2d_dims(file_idx)
        call io_manager_add_dim(file_idx, model_aero_lev)

        ! Add 'date' variable.
        call model_aero_list%append("date", "date integer", "1", data_type="integer")
        call model_aero_list%get_tail_values(date)

        ptr1 => model_aero_list%get_head()
        do i = 1, model_aero_list%get_num_var()
            if (ptr1%get_name() == "PS") then
                call io_manager_def_var(file_idx, ptr1, ["lon ","lat ","time"])
            else if (ptr1%get_name() == "date") then
                call io_manager_def_var(file_idx, ptr1, ["time"])
            else
                call io_manager_def_var(file_idx, ptr1, ["lon ","lat ","lev ","time"])
            end if
            ptr1 => ptr1%next
        end do

        do time = 1, num_aero_time
            date = aero_date(time)
            ptr2 => model_aero_list%get_var("date")
            call io_manager_put_var(file_idx, ptr2, rec=time)
            ! ------------------------------------------------------------------
            call notice(sub_name, "Interpolate aerosol data on "// &
                trim(integer_to_string(aero_date(time))))
            ! ------------------------------------------------------------------
            call source_data_read_aerosol
            ! ------------------------------------------------------------------
            call io_manager_put_var(file_idx, "time", aero_time(time), rec=time)
            ! ------------------------------------------------------------------
            ! interpolate surface pressure onto model grids
            call model_gears_interp_h(aero_lon, aero_lat, aero_ps, &
                model_lon, model_lat, model_ps, 1)
            ptr2 => model_aero_list%get_var("PS")
            call io_manager_put_var(file_idx, ptr2, rec=time)
            ! ------------------------------------------------------------------
            ! calculate pressure of horizontal model grids + vertical data levels
            do k = 1, num_aero_lev
                do j = 1, num_model_lat
                    do i = 1, num_model_lon
                        tmp_p(i,j,k) = aero_hyam(k)*aero_p0+aero_hybm(k)*model_ps(i,j)
                    end do
                end do
            end do
            ! ------------------------------------------------------------------
            ! calculate pressure of model grids
            do k = 1, num_lev
                do j = 1, num_model_lat
                    do i = 1, num_model_lon
                        model_p(i,j,k) = model_pt+model_lev_bnds(k)*(model_ps(i,j)-model_pt)
                    end do
                end do
            end do
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
                    tmp4 = 0.0d0
                    do i = 1, num_sub_aero
                        ptr1 => data_aero_list%get_var( &
                            get_substr(data_aero_names(aero), "+", i))
                        select type (ptr1)
                        type is (var3d_d)
                            tmp1 => ptr1%get_values()
                        end select
                        tmp4 = tmp4+tmp1
                    end do
                    tmp1 => tmp4
                end if
                ! ==============================================================
                ! integrate source aerosol data
                do j = 1, num_aero_lat
                    do i = 1, num_aero_lon
                        do k = 1, num_aero_lev
                            tmp4(i,j,k) = tmp1(i,j,k)*aero_dp(i,j,k)/g*int_coef(aero)
                        end do
                        do k = num_aero_lev-1, 1, -1
                            tmp4(i,j,k) = tmp4(i,j,k)+tmp4(i,j,k+1)
                        end do
                    end do
                end do
                ! ==============================================================
                ! get model aerosol data pointer
                ptr2 => model_aero_list%get_var(model_aero_names(aero))
                select type (ptr2)
                type is (var3d_d)
                    tmp2 => ptr2%get_values()
                end select
                ! ==============================================================
                ! interpolate (first horizontal, then vertical)
                do k = 1, num_aero_lev
                    call model_gears_interp_h(aero_lon, aero_lat, tmp4(:,:,k), &
                        model_lon, model_lat, tmp3(:,:,k), 1)
                end do
                do j = 1, num_model_lat
                    do i = 1, num_model_lon
                        call model_gears_interp_v(tmp_p(i,j,:), tmp3(i,j,:), &
                            model_p(i,j,:), tmp2(i,j,:), 0)
                        ! constrain the boundary cells (no expolation)
                        if (model_p(i,j,1) < tmp_p(i,j,1)) then
                            tmp2(i,j,1) = tmp3(i,j,1)
                        end if
                        if (model_p(i,j,num_lev) > tmp_p(i,j,num_aero_lev)) then
                            tmp2(i,j,num_lev) = tmp3(i,j,num_aero_lev)
                        end if
                    end do
                end do
                ! ==============================================================
                call io_manager_put_var(file_idx, ptr2, rec=time)
            end do
        end do

        call io_manager_close_file(file_idx)

        call notice(sub_name, "File "//trim(output_file_name)//" is generated!")

    end subroutine model_fc_interp_aerosol

end module model_fc
