module model_ic
    
    use utils
    use model_gears
    use model_grids
    use source_data
    use constants
    use variable

    implicit none

    private

    public model_ic_calc_topo
    public model_ic_calc_ps
    public model_ic_calc_ts
    public model_ic_calc_landm
    public model_ic_interp
    public model_ic_write
    public model_ic_file_name

    type(var_list) model_vars

    real(8), pointer :: u(:,:,:), v(:,:,:), T(:,:,:), Q(:,:,:), p(:,:,:)

    real(8), pointer :: landfrac(:,:), phis(:,:), sgh(:,:), landm(:,:)

    real(8), pointer :: mslp(:,:), mslt(:,:), ps(:,:), ts(:,:)

    ! ??????????????????????????????????????????????????????????????????????????
    ! TODO: The following variables are not clear!
    real(8), pointer :: cwat(:,:,:), snowhice(:,:), tsice(:,:)
    real(8), pointer :: ts1(:,:), ts2(:,:), ts3(:,:), ts4(:,:)
    ! ??????????????????????????????????????????????????????????????????????????

contains

#include "calc_model_topo.F90"

    subroutine model_ic_calc_topo(topo_file)

        character(*), intent(in) :: topo_file

        character(30), parameter :: sub_name = "model_ic_calc_topo"

        integer ncid, dim_id, var_id, ierr

        integer num_topo_lon, num_topo_lat
        real(8), allocatable :: topo_lon(:), topo_lat(:)
        real(8), allocatable :: topo(:,:)

        call notice(sub_name, "Calculate GAMIL topography")

        ierr = nf90_open(topo_file, nf90_nowrite, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lon", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_topo_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lat", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_topo_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        allocate(topo_lon(num_topo_lon))
        allocate(topo_lat(num_topo_lat))
        allocate(topo(num_topo_lon,num_topo_lat))

        ierr = nf90_inq_varid(ncid, "lon", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, topo_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lat", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, topo_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "htopo", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, topo)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        call model_vars%append("LANDFRAC", "land fraction (0: ocean, 1: land)", "1", model_2d_dims)
        call model_vars%get_tail_values(landfrac)
        call model_vars%append("PHIS", "surface geopotential", "m2 s-2", model_2d_dims)
        call model_vars%get_tail_values(phis)
        call model_vars%append("SGH", "standard deviation of surface geopotential", "m", model_2d_dims)
        call model_vars%get_tail_values(sgh)

        ! call Bin Wang's subroutine to calculate the topography of model
        call calc_model_topo(num_topo_lon, num_topo_lat, topo_lon, topo_lat, topo, &
                             num_model_lon, num_model_lat, model_lon, model_lat, &
                             landfrac, phis, sgh)

        ! change the unit of phis to m2 s-2
        phis = phis*g

    end subroutine model_ic_calc_topo

    subroutine model_ic_calc_landm

        real(8), parameter :: coastal_check_land = 0.5
        real(8), parameter :: coastal_check_ocean = 0.1
        integer, parameter :: search_levels = 6
        real(8), parameter :: idx_dist_scale = 6.0

        integer search_levels_lat
        real(8) idx_dist_scale_lat
        real(8) tmp
        real(8) idx_dist
        integer i, j, k, m, n, ii, jj
        logical checked(num_model_lon,num_model_lat)
        logical is_coastal

        character(30), parameter :: sub_name = "model_ic_calc_landm"

        call notice(sub_name, "Calcumate GAMIL land transition mask")
        call report_warning(sub_name, &
            "The algorithm of calculation of land transition mask should be revised!")

        call model_vars%append("LANDM", "land transition mask", "1", model_2d_dims)
        call model_vars%get_tail_values(landm)

        if (.not. associated(landfrac)) then
            call report_error(sub_name, "The GAMIL topography should be calculated first")
        end if

        landm = landfrac

        do j = 1, num_model_lat
            tmp = (1.0+5.0*exp(-((90.0-abs(model_lat(j))))**2/50.0))
            search_levels_lat = tmp*search_levels
            idx_dist_scale_lat = tmp*idx_dist_scale
            do i = 1, num_model_lon
                is_coastal = .false.
                ! find the coastal grids
                if (landfrac(i,j) > coastal_check_land) then
                    landm(i,j) = min(1.0, 2.0*landfrac(i,j))
                    do n = j-1, j+1
                    do m = i-1, i+1
                        if (m > num_model_lon) then
                            ii = m-num_model_lon
                        else if (m < 1) then
                            ii = m+num_model_lon
                        else
                            ii = m
                        end if
                        jj = min(num_model_lat, max(1, n))
                        if (landfrac(ii,jj) < coastal_check_ocean) then
                            is_coastal = .true.
                        end if
                    end do
                    end do
                end if
                ! skip the non-coastal land and ocean grids
                if (.not. is_coastal) cycle
                ! here the grids should be coastal
                checked = .false.
                do k = 1, search_levels_lat
                    do n = j-k, j+k
                    do m = i-k, i+k
                        if (m > num_model_lon) then
                            ii = m-num_model_lon
                        else if (m < 1) then
                            ii = m+num_model_lon
                        else
                            ii = m
                        end if
                        jj = min(num_model_lat, max(1, n))
                        if (checked(ii,jj)) cycle
                        if (k == 1) then
                            landm(ii,jj) = max(landm(ii,jj), landm(i,j))
                        else
                            if (landm(ii,jj) < landm(i,j)) then
                                idx_dist = sqrt(dble(m-i)**2+dble(n-j)**2)
                                landm(ii,jj) = max(landm(ii,jj), &
                                    landm(i,j)*exp(-(idx_dist/idx_dist_scale_lat)**2))
                            end if
                        end if
                        checked(ii,jj) = .true.
                    end do
                    end do
                end do
            end do
        end do

    end subroutine model_ic_calc_landm

    ! --------------------------------------------------------------------------
    ! Procedure name: model_ic_calc_ps
    ! Description:
    !   Calculate the surface pressure on the model grids from the source data
    !   according to the hydrostatic approximation.
    ! Author:
    !   Li Dong, dongli@lasg.iap.ac.cn
    ! --------------------------------------------------------------------------

    subroutine model_ic_calc_ps
    
        ! reference state
        real(8) p0(num_model_lon,num_model_lat)
        real(8) t0(num_model_lon,num_model_lat)
        real(8) z0(num_model_lon,num_model_lat)
        real(8) temp1, temp2, temp3
        integer i, j

        character(30), parameter :: sub_name = "model_ic_calc_ps"

        call notice(sub_name, "Calculate GAMIL surface pressure")

        if (.not. associated(phis)) then
            call report_error(sub_name, "GAMIL topography should be calculated first")
        end if

        if (uvtq_data_type == era_interim .and. &
            .not. allocated(data_mslp)) then
            call report_error(sub_name, &
                "There is no mean sea level pressure data in source data")
        end if

        if (uvtq_data_type == model_data .and. &
            .not. allocated(data_ps) .and. .not. allocated(data_ts)) then
            call report_error(sub_name, &
                "There is no surface pressure and surface temperature data in source data")
        end if

        ! ----------------------------------------------------------------------
        temp1 = 0.006*Rd/g
        call model_vars%append("PS", "surface pressure", "Pa", model_2d_dims)
        call model_vars%get_tail_values(ps)

        ! ----------------------------------------------------------------------
        ! get the reference state
        if (uvtq_data_type == era_interim) then
            call model_gears_interp_h(data_lon, data_lat, data_mslp, &
                                      model_lon, model_lat, p0, 1)
            call model_gears_interp_h(data_lon, data_lat, data_t(:,:,num_data_lev), &
                                      model_lon, model_lat, t0, 1)
            do j = 1, num_model_lat
                do i = 1, num_model_lon
                    temp2 = p0(i,j)/data_lev(num_data_lev)
                    t0(i,j) = t0(i,j)*temp2**temp1
                end do
            end do
            z0 = 0.0d0
            ! record mean sea level pressure and temperature
            call model_vars%append("MSLP", "mean sea level pressure", "Pa", model_2d_dims)
            call model_vars%get_tail_values(mslp)
            call model_vars%append("MSLT", "mean sea level temperaure", "K", model_2d_dims)
            call model_vars%get_tail_values(mslt)
            mslp = p0
            mslt = t0
        else if (uvtq_data_type == model_data) then
            call model_gears_interp_h(data_lon, data_lat, data_ps, &
                                      model_lon, model_lat, p0, 1)
            call model_gears_interp_h(data_lon, data_lat, data_ts, &
                                      model_lon, model_lat, t0, 1)
            call model_gears_interp_h(data_lon, data_lat, data_phis, &
                                      model_lon, model_lat, z0, 1)
            z0 = z0/g
        end if

        ! ----------------------------------------------------------------------
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                temp3 = 1.0d0-0.006d0*(phis(i,j)/g-z0(i,j))/t0(i,j)
                ps(i,j) = p0(i,j)*temp3**(1.0d0/temp1)
            end do
        end do

    end subroutine model_ic_calc_ps

    ! --------------------------------------------------------------------------
    ! Procedure name: model_ic_calc_ts
    ! Description:
    !   Calculate surface temperature on the model grids.
    ! Author:
    !
    ! --------------------------------------------------------------------------

    subroutine model_ic_calc_ts

        real(8) t0(num_model_lon,num_model_lat)
        real(8) p0(num_model_lon,num_model_lat)
        real(8) temp
        integer i, j

        character(30), parameter :: sub_name = "model_ic_calc_ts"

        call notice(sub_name, "Calculate GAMIL surface temperature")

        if (.not. associated(ps)) then
            call report_warning(sub_name, "GAMIL surface pressure has not been calculated")
            call model_ic_calc_ps
        end if

        call model_vars%append("TS", "surface temperature", "K", model_2d_dims)
        call model_vars%get_tail_values(ts)

        if (uvtq_data_type == era_interim) then
            t0 = mslt
            p0 = mslp
        else if (uvtq_data_type == model_data) then
            call model_gears_interp_h(data_lon, data_lat, data_ts, &
                                      model_lon, model_lat, t0, 1)
            call model_gears_interp_h(data_lon, data_lat, data_ps, &
                                      model_lon, model_lat, p0, 1)
        end if

        temp = 0.006*Rd/g
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                ts(i,j) = t0(i,j)*(ps(i,j)/p0(i,j))**temp
            end do
        end do
    
    end subroutine model_ic_calc_ts

    subroutine model_ic_interp

        ! horizontal interpolated variables
        real(8), allocatable :: tmp(:,:,:), tmp_p(:,:,:)

        integer i, j, k

        character(30), parameter :: sub_name = "model_ic_interp"

        call notice(sub_name, "Interpolate source data onto the model grids")

        call model_ic_calc_landm
        call model_ic_calc_ps
        call model_ic_calc_ts
        ! "p" is not added into the main data
        allocate (p(num_model_lon,num_model_lat,num_model_lev))
        call model_grids_calc_p(model_hyam, model_hybm, ps, model_p0, p)

        ! intermediate grids (model horizontal grids + data vertical grids)
        allocate (tmp(num_model_lon,num_model_lat,num_data_lev))
        if (uvtq_data_type == model_data) then
            ! when interpolating from other model data, we need to calculate
            ! the pressure on the intermediate grid
            allocate (tmp_p(num_model_lon,num_model_lat,num_data_lev))
            call model_grids_calc_p(data_hyam, data_hybm, data_ps, data_p0, data_p)
            do k = 1, num_data_lev
                call model_gears_interp_h(data_lon, data_lat, data_p(:,:,k), &
                    model_lon, model_lat, tmp_p(:,:,k), 1)
            end do
        end if

        call model_vars%append("U", "zonal wind", "m s-1", model_3d_dims)
        call model_vars%get_tail_values(u)
        call model_vars%append("V", "meridional wind", "m s-1", model_3d_dims)
        call model_vars%get_tail_values(v)
        call model_vars%append("T", "temperature", "K", model_3d_dims)
        call model_vars%get_tail_values(T)
        call model_vars%append("Q", "specific humidity", "kg kg-1", model_3d_dims)
        call model_vars%get_tail_values(Q)

        do k = 1, num_data_lev
            call model_gears_interp_h(data_lon, data_lat, data_u(:,:,k), &
                                      model_lon, model_lat, tmp(:,:,k), 0)
        end do
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                if (uvtq_data_type == era_interim) then
                    call model_gears_interp_v(data_lev, tmp(i,j,:), p(i,j,:), u(i,j,:), 0)
                else if (uvtq_data_type == model_data) then
                    call model_gears_interp_v(tmp_p(i,j,:), tmp(i,j,:), p(i,j,:), u(i,j,:), 0)
                end if
            end do
        end do

        do k = 1, num_data_lev
            call model_gears_interp_h(data_lon, data_lat, data_v(:,:,k), &
                                      model_lon, model_lat, tmp(:,:,k), 0)
        end do
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                if (uvtq_data_type == era_interim) then
                    call model_gears_interp_v(data_lev, tmp(i,j,:), p(i,j,:), v(i,j,:), 0)
                else if (uvtq_data_type == model_data) then
                    call model_gears_interp_v(tmp_p(i,j,:), tmp(i,j,:), p(i,j,:), v(i,j,:), 0)
                end if
            end do
        end do

        do k = 1, num_data_lev
            call model_gears_interp_h(data_lon, data_lat, data_t(:,:,k), &
                                      model_lon, model_lat, tmp(:,:,k), 1)
        end do
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                if (uvtq_data_type == era_interim) then
                    call model_gears_interp_v(data_lev, tmp(i,j,:), p(i,j,:), T(i,j,:), 1)
                else if (uvtq_data_type == model_data) then
                    call model_gears_interp_v(tmp_p(i,j,:), tmp(i,j,:), p(i,j,:), T(i,j,:), 1)
                end if
            end do
        end do

        do k = 1, num_data_lev
            call model_gears_interp_h(data_lon, data_lat, data_q(:,:,k), &
                                      model_lon, model_lat, tmp(:,:,k), 1)
        end do
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                if (uvtq_data_type == era_interim) then
                    call model_gears_interp_v(data_lev, tmp(i,j,:), p(i,j,:), Q(i,j,:), 0)
                else if (uvtq_data_type == model_data) then
                    call model_gears_interp_v(tmp_p(i,j,:), tmp(i,j,:), p(i,j,:), Q(i,j,:), 0)
                end if
            end do
        end do

        if (uvtq_data_type == model_data) then
            call model_vars%append("CWAT", "total grid box averaged condensate amount (liquid + ice)", "kg kg-1", model_3d_dims)
            call model_vars%get_tail_values(cwat)
            call model_vars%append("SNOWHICE", "water equivalent snow depth", "m", model_2d_dims)
            call model_vars%get_tail_values(snowhice)
            call model_vars%append("TSICE", "", "K", model_2d_dims)
            call model_vars%get_tail_values(tsice)
            call model_vars%append("TS1", "1 subsoil temperature", "K", model_2d_dims)
            call model_vars%get_tail_values(ts1)
            call model_vars%append("TS2", "2 subsoil temperature", "K", model_2d_dims)
            call model_vars%get_tail_values(ts2)
            call model_vars%append("TS3", "3 subsoil temperature", "K", model_2d_dims)
            call model_vars%get_tail_values(ts3)
            call model_vars%append("TS4", "4 subsoil temperature", "K", model_2d_dims)
            call model_vars%get_tail_values(ts4)

            do k = 1, num_data_lev
                call model_gears_interp_h(data_lon, data_lat, data_cwat(:,:,k), &
                                          model_lon, model_lat, tmp(:,:,k), 1)
            end do
            do j = 1, num_model_lat
                do i = 1, num_model_lon
                    if (uvtq_data_type == era_interim) then
                        call model_gears_interp_v(data_lev, tmp(i,j,:), p(i,j,:), cwat(i,j,:), 0)
                    else if (uvtq_data_type == model_data) then
                        call model_gears_interp_v(tmp_p(i,j,:), tmp(i,j,:), p(i,j,:), cwat(i,j,:), 0)
                    end if
                end do
            end do
    
            call model_gears_interp_h(data_lon, data_lat, data_snowhice, &
                                      model_lon, model_lat, snowhice, 1)
            call model_gears_interp_h(data_lon, data_lat, data_tsice, &
                                      model_lon, model_lat, tsice, 1)
    
            call model_gears_interp_h(data_lon, data_lat, data_subts(:,:,1), &
                                      model_lon, model_lat, ts1, 1)
            call model_gears_interp_h(data_lon, data_lat, data_subts(:,:,2), &
                                      model_lon, model_lat, ts2, 1)
            call model_gears_interp_h(data_lon, data_lat, data_subts(:,:,3), &
                                      model_lon, model_lat, ts3, 1)
            call model_gears_interp_h(data_lon, data_lat, data_subts(:,:,4), &
                                      model_lon, model_lat, ts4, 1)
        end if

    end subroutine model_ic_interp

    subroutine model_ic_write

        character(300) file_name

        integer file_idx, i
        class(var), pointer :: ptr

        character(30), parameter :: sub_name = "model_ic_write"

        if (.not. associated(ts)) then
            call report_warning(sub_name, "GAMIL surface temperature has not been calculated yet")
            call model_ic_calc_ts
        end if

        file_name = model_ic_file_name()

        call io_manager_create_file(file_name, file_idx)
        call io_manager_def_dim(file_idx, "time", "integer", &
            units="days since "//trim(data_date))
        call model_grids_add_dims(file_idx)

        call io_manager_put_var(file_idx, "time", 0, rec=1)

        ptr => model_vars%get_head()
        do i = 1, model_vars%get_num_var()
            select type (ptr)
            type is (var2d_d)
                call io_manager_def_var(file_idx, ptr, ["lon ","lat ","time"])
            type is (var3d_d)
                call ptr%reshape([1,3,2])
                call io_manager_def_var(file_idx, ptr, ["lon ","lev ","lat ","time"])
            end select
            ptr => ptr%next
        end do

        ptr => model_vars%get_head()
        do i = 1, model_vars%get_num_var()
            call io_manager_put_var(file_idx, ptr, rec=1)
            ptr => ptr%next
        end do

        call io_manager_add_att(file_idx, "note", "Created by GAMIL-INIT")

        call io_manager_close_file(file_idx)

        call notice(sub_name, "File "//trim(file_name)//" has been generated")

    end subroutine model_ic_write

    character(300) function model_ic_file_name()

        model_ic_file_name = "ic.gamil."// &
            trim(to_string(num_model_lon))//"x"// &
            trim(to_string(num_model_lat))//"x"// &
            trim(to_string(num_model_lev))//"."// &
            trim(data_date)//".nc"

    end function model_ic_file_name

end module model_ic
