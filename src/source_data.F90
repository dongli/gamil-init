module source_data

    use utils
    use namelist_mod
    use variable

    implicit none

    ! --------------------------------------------------------------------------
    ! u, v, t, q data set
    integer, parameter :: era_interim = 1
    integer, parameter :: model_data  = 2
    integer num_data_lon, num_data_lat, num_data_lev
    real(8), allocatable :: data_lon(:)
    real(8), allocatable :: data_lat(:)
    real(8), allocatable :: data_lev(:)
    real(8), allocatable :: data_hyam(:)
    real(8), allocatable :: data_hybm(:)
    real(8), allocatable :: data_hyai(:)
    real(8), allocatable :: data_hybi(:)
    real(8), allocatable :: data_u(:,:,:)
    real(8), allocatable :: data_v(:,:,:)
    real(8), allocatable :: data_t(:,:,:)
    real(8), allocatable :: data_rh(:,:,:)
    real(8), allocatable :: data_q(:,:,:)
    real(8), allocatable :: data_z(:,:,:)
    real(8), allocatable :: data_p(:,:,:)
    real(8), allocatable :: data_mslp(:,:)
    real(8), allocatable :: data_ps(:,:)
    real(8), allocatable :: data_phis(:,:)
    real(8), allocatable :: data_ts(:,:)
    real(8) data_p0, data_pt

    real(8), allocatable :: data_cwat(:,:,:)
    real(8), allocatable :: data_snowhice(:,:)
    real(8), allocatable :: data_tsice(:,:)
    real(8), allocatable :: data_subts(:,:,:)
    ! data date is read from the units of time variable from file
    character(30) data_date

    ! --------------------------------------------------------------------------
    ! Ozone data set.
    integer num_ozone_lon, num_ozone_lat, num_ozone_lev, num_ozone_time
    real(8), allocatable :: data_ozone_lon(:)
    real(8), allocatable :: data_ozone_lat(:)
    real(8), allocatable :: data_ozone_lev(:)
    real(8), allocatable :: data_ozone(:,:,:,:)

    ! --------------------------------------------------------------------------
    ! Aerosol data set.
    integer aero_ncid, aero_time_step
    integer num_aero_var
    integer aero_var_ids(nf90_max_vars)
    logical aero_var_read_each_time(nf90_max_vars)
    integer num_aero_lon, num_aero_lat, num_aero_lev, num_aero_time
    real(8), pointer :: aero_p0
    real(8), pointer :: aero_lon(:)
    real(8), pointer :: aero_lat(:)
    real(8), pointer :: aero_hyam(:)
    real(8), pointer :: aero_hybm(:)
    real(8), pointer :: aero_hyai(:)
    real(8), pointer :: aero_hybi(:)
    real(8), pointer :: aero_time(:)
    integer, pointer :: aero_date(:)
    integer, pointer :: aero_datesec(:)
    real(8), pointer :: aero_ps(:,:)
    real(8), allocatable :: aero_p(:,:,:)
    real(8), allocatable :: aero_dp(:,:,:)
    type(var_list) data_aero_list

    ! --------------------------------------------------------------------------
    ! Temporary pointers.
    integer, pointer :: tmp_0d_i
    real(8), pointer :: tmp_0d_d
    integer, pointer :: tmp_1d_i(:)
    real(8), pointer :: tmp_1d_d(:)
    real(8), pointer :: tmp_2d_d(:,:)
    real(8), pointer :: tmp_3d_d(:,:,:)

contains

    subroutine source_data_read_uvtq(file_name)

        character(*), intent(in) :: file_name

        character(30), parameter :: sub_name = "source_data_read_uvtq"

        select case (uvtq_data_type)
        case (era_interim)
            call source_data_read_era_interim(file_name)
        case (model_data)
            call source_data_read_model_data(file_name)
        case default
            call report_error(sub_name, "Unknown uvtq data type!")
        end select
    
    end subroutine source_data_read_uvtq

    subroutine source_data_read_era_interim(file_name)

        character(*), intent(in) :: file_name

        integer ncid, dim_id, var_id, ierr
        integer count_2d(3), count_3d(4)
        character(30) units

        real(8) e
        integer i, j, k

        character(30), parameter :: sub_name = "source_data_read_era_interim"

        call notice(sub_name, "Read ERA-interim reanalysis data "//trim(file_name))

        ierr = nf90_open(file_name, nf90_nowrite, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lon", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lat", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lev", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ! ----------------------------------------------------------------------
        count_2d = [num_data_lon,num_data_lat,1]
        count_3d = [num_data_lon,num_data_lat,num_data_lev,1]
        allocate(data_lon(num_data_lon))
        allocate(data_lat(num_data_lat))
        allocate(data_lev(num_data_lev))
        allocate(data_u(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_v(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_t(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_rh(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_z(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_mslp(num_data_lon,num_data_lat))
        ! ----------------------------------------------------------------------
        ierr = nf90_inq_varid(ncid, "time", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_att(ncid, var_id, "units", units)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        data_date = get_substr(units, " ", 3)

        ierr = nf90_inq_varid(ncid, "lon", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lat", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lev", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_att(ncid, var_id, "units", units)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        call check_pressure_units(sub_name, data_lev, units)

        ierr = nf90_inq_varid(ncid, "u", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_u, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "v", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_v, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "T", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_t, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "rh", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_rh, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "Z", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_z, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "mslp", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_mslp, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ! ----------------------------------------------------------------------
        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ! ----------------------------------------------------------------------
        call notice(sub_name, "Calculate specific humidity in ERA-interim reanalysis data")
        ! ----------------------------------------------------------------------
        allocate(data_q(num_data_lon,num_data_lat,num_data_lev))
        do k = 1, num_data_lev
        do j = 1, num_data_lat
        do i = 1, num_data_lon
            e = data_rh(i,j,k)*0.01d0*6.11d0*10.0d0**( &
                7.5d0*(data_t(i,j,k)-273.15)/(237.3d0+(data_t(i,j,k)-273.15d0)))
            data_q(i,j,k) = 0.622d0*e/(data_lev(k)-0.378d0*e)
        end do
        end do
        end do

    end subroutine source_data_read_era_interim

    subroutine source_data_read_model_data(file_name)
    
        character(*), intent(in) :: file_name

        integer ncid, dim_id, var_id, ierr
        integer count_2d(3), count_3d(4)
        character(30) units

        integer k
        real(8), allocatable :: tmp(:,:,:)

        character(30), parameter :: sub_name = "source_data_read_model_data"

        call notice(sub_name, "Reading model data "//trim(file_name))

        ierr = nf90_open(file_name, nf90_nowrite, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ! ----------------------------------------------------------------------
        ierr = nf90_inq_dimid(ncid, "lon", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lat", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lev", dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_inquire_dimension(ncid, dim_id, len=num_data_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ! ----------------------------------------------------------------------
        count_2d = [num_data_lon,num_data_lat,1]
        count_3d = [num_data_lon,num_data_lev,num_data_lat,1]
        allocate(data_lon(num_data_lon))
        allocate(data_lat(num_data_lat))
        allocate(data_lev(num_data_lev))
        allocate(data_hyam(num_data_lev))
        allocate(data_hybm(num_data_lev))
        allocate(data_hyai(num_data_lev+1))
        allocate(data_hybi(num_data_lev+1))
        allocate(data_u(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_v(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_t(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_q(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_p(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_ps(num_data_lon,num_data_lat))
        allocate(data_ts(num_data_lon,num_data_lat))
        allocate(data_phis(num_data_lon,num_data_lat))
        allocate(data_cwat(num_data_lon,num_data_lat,num_data_lev))
        allocate(data_snowhice(num_data_lon,num_data_lat))
        allocate(data_tsice(num_data_lon,num_data_lat))
        allocate(data_subts(num_data_lon,num_data_lat,4))
        allocate(tmp(num_data_lon,num_data_lev,num_data_lat))
        ! ----------------------------------------------------------------------
        ierr = nf90_inq_varid(ncid, "time", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_att(ncid, var_id, "units", units)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        data_date = get_substr(units, " ", 3)

        ierr = nf90_inq_varid(ncid, "lon", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lat", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! TODO: Do we need "lev" any more?
        ierr = nf90_inq_varid(ncid, "lev", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! ----------------------------------------------------------------------
        ! Read in hybrid sigma-pressure coordinate parameters
        ierr = nf90_inq_varid(ncid, "P0", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_p0)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "hyam", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_hyam)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "hybm", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_hybm)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "hyai", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_hyai)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "hybi", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_hybi)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! ----------------------------------------------------------------------
        ierr = nf90_inq_varid(ncid, "U", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, tmp, [1,1,1,1], count_3d)
        do k = 1, num_data_lev
            data_u(:,:,k) = tmp(:,k,:)
        end do
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "V", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, tmp, [1,1,1,1], count_3d)
        do k = 1, num_data_lev
            data_v(:,:,k) = tmp(:,k,:)
        end do
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "T", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, tmp, [1,1,1,1], count_3d)
        do k = 1, num_data_lev
            data_t(:,:,k) = tmp(:,k,:)
        end do
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "Q", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, tmp, [1,1,1,1], count_3d)
        do k = 1, num_data_lev
            data_q(:,:,k) = tmp(:,k,:)
        end do
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "PS", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_ps, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TS", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_ts, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "PHIS", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_phis, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "CWAT", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, tmp, [1,1,1,1], count_3d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        do k = 1, num_data_lev
            data_cwat(:,:,k) = tmp(:,k,:)
        end do

        ierr = nf90_inq_varid(ncid, "SNOWHICE", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_snowhice, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TSICE", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_tsice, [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TS1", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_subts(:,:,1), [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TS2", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_subts(:,:,2), [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TS3", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_subts(:,:,3), [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "TS4", var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)
        ierr = nf90_get_var(ncid, var_id, data_subts(:,:,4), [1,1,1], count_2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        call report_warning(sub_name, "Assume the vertical coordinate of input model data is the same of output!")
        ! ----------------------------------------------------------------------
        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine source_data_read_model_data

    subroutine source_data_read_ozone(file_name)

        character(*), intent(in) :: file_name

        integer ncid, ierr
        integer time_dim_id
        integer lon_dim_id, lat_dim_id, lev_dim_id
        integer lon_var_id, lat_var_id, lev_var_id
        integer ozone_var_id

        character(30) exp_id, dim_name
        integer ndims, dimids(4), dim_len

        character(30), parameter :: sub_name = "source_data_read_ozone"

        call notice(sub_name, "Reading ozone data "//trim(file_name))

        ierr = nf90_open(file_name, nf90_nowrite, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_att(ncid, nf90_global, "experiment_id", exp_id)
        if (ierr == nf90_noerr) then
            if (exp_id == "APE") then
                call notice(sub_name, "Ozone data is for aqua-planet experiments")
            end if
        else
            call report_warning(sub_name, "Unknown experiment ID of ozone data")
        end if

        ierr = nf90_inq_varid(ncid, "ozone", ozone_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_variable(ncid, ozone_var_id, ndims=ndims, dimids=dimids)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! start here to handle time dimension 2011-06-28
        if (ndims == 2) then
            num_ozone_lon = 1

            ierr = nf90_inquire_dimension(ncid, dimids(1), name=dim_name, len=dim_len)
            call handle_netcdf_error(sub_name, __LINE__, ierr)

            if (dim_name == "lat") then
                num_ozone_lat = dim_len
            end if

            ierr = nf90_inquire_dimension(ncid, dimids(2), name=dim_name, len=dim_len)
            call handle_netcdf_error(sub_name, __LINE__, ierr)

            if (dim_name == "lev") then
                num_ozone_lev = dim_len
            end if
        else if (ndims == 3) then
            call report_error(sub_name, "3D ozone data has not been handled yet")
        end if

        allocate(data_ozone(num_ozone_lon,num_ozone_lat,num_ozone_lev,num_ozone_time))

        ierr = nf90_get_var(ncid, ozone_var_id, data_ozone)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine source_data_read_ozone

    subroutine source_data_start_aerosol(file_name)

        character(*), intent(in) :: file_name

        character(10) name
        character(100) long_name, units
        integer ierr, lon_dim_id, lat_dim_id, lev_dim_id, ilev_dim_id, time_dim_id
        integer dim_size, i, dim_id, data_type, num_dim, dim_ids(4)

        character(50), parameter :: sub_name = "source_data_start_aerosol"

        call notice(sub_name, "Start read aerosol file """//trim(file_name)//"""")

        ! Open the aerosol data file.
        ierr = nf90_open(file_name, nf90_nowrite, aero_ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! Get all the variable ID's.
        ierr = nf90_inq_varids(aero_ncid, num_aero_var, aero_var_ids)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! Get the dimension sizes.
        do i = 1, num_aero_var
            ierr = nf90_inquire_variable(aero_ncid, aero_var_ids(i), name=name)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            if (.not. any(["lon ","lat ","lev ","ilev","time"] == name)) cycle
            ierr = nf90_inq_dimid(aero_ncid, name, dim_id)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            ierr = nf90_inquire_dimension(aero_ncid, dim_id, len=dim_size)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            select case (name)
            case ("lon")
                num_aero_lon = dim_size
                ierr = nf90_inq_dimid(aero_ncid, "lon", lon_dim_id)
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            case ("lat")
                num_aero_lat = dim_size
                ierr = nf90_inq_dimid(aero_ncid, "lat", lat_dim_id)
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            case ("lev")
                num_aero_lev = dim_size
                ierr = nf90_inq_dimid(aero_ncid, "lev", lev_dim_id)
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            case ("ilev")
                ierr = nf90_inq_dimid(aero_ncid, "ilev", ilev_dim_id)
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            case ("time")
                num_aero_time = dim_size
                ierr = nf90_inq_dimid(aero_ncid, "time", time_dim_id)
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            end select
        end do

        ! Define variables.
        do i = 1, num_aero_var
            ierr = nf90_inquire_variable(aero_ncid, aero_var_ids(i), name=name, &
                xtype=data_type, ndims=num_dim, dimids=dim_ids)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            ierr = nf90_get_att(aero_ncid, aero_var_ids(i), "long_name", long_name)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            ierr = nf90_get_att(aero_ncid, aero_var_ids(i), "units", units)
            if (ierr /= nf90_noerr) then
                units = "1"
            end if
            select case (num_dim)
            case (0)
                call data_aero_list%append(name, long_name, units, data_type=data_type)
                select type (var => data_aero_list%get_var(name))
                type is (var0d_d)
                    tmp_0d_d => var%get_values()
                    ierr = nf90_get_var(aero_ncid, aero_var_ids(i), tmp_0d_d)
                    call handle_netcdf_error(sub_name, __LINE__, ierr)
                end select
                aero_var_read_each_time(i) = .false.
            case (1) ! 1D variable
                if (dim_ids(1) == lon_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lon], data_type=data_type)
                else if (dim_ids(1) == lat_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lat], data_type=data_type)
                else if (dim_ids(1) == lev_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lev], data_type=data_type)
                else if (dim_ids(1) == ilev_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lev+1], data_type=data_type)
                else if (dim_ids(1) == time_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_time], data_type=data_type)
                end if
                select type (var => data_aero_list%get_var(name))
                type is (var1d_i)
                    tmp_1d_i => var%get_values()
                    ierr = nf90_get_var(aero_ncid, aero_var_ids(i), tmp_1d_i)
                    call handle_netcdf_error(sub_name, __LINE__, ierr)
                type is (var1d_d)
                    tmp_1d_d => var%get_values()
                    ierr = nf90_get_var(aero_ncid, aero_var_ids(i), tmp_1d_d)
                    call handle_netcdf_error(sub_name, __LINE__, ierr)
                end select
                aero_var_read_each_time(i) = .false.
            case (3) ! 2D variable
                if (dim_ids(1) == lon_dim_id .and. dim_ids(2) == lat_dim_id .and. dim_ids(3) == time_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lon,num_aero_lat], data_type=data_type)
                end if
                aero_var_read_each_time(i) = .true.
            case (4) ! 3D variable
                if (dim_ids(1) == lon_dim_id .and. dim_ids(2) == lat_dim_id .and. dim_ids(3) == lev_dim_id .and. dim_ids(4) == time_dim_id) then
                    call data_aero_list%append(name, long_name, units, [num_aero_lon,num_aero_lat,num_aero_lev], data_type=data_type)
                end if
                aero_var_read_each_time(i) = .true.
            end select
        end do

        call data_aero_list%get_values("P0", aero_p0)
        call data_aero_list%get_values("time", aero_time)
        call data_aero_list%get_values("date", aero_date)
        call data_aero_list%get_values("datesec", aero_datesec)
        call data_aero_list%get_values("lon", aero_lon)
        call data_aero_list%get_values("lat", aero_lat)
        call data_aero_list%get_values("hyam", aero_hyam)
        call data_aero_list%get_values("hybm", aero_hybm)
        call data_aero_list%get_values("hyai", aero_hyai)
        call data_aero_list%get_values("hybi", aero_hybi)
        call data_aero_list%get_values("PS", aero_ps)

        ! ----------------------------------------------------------------------
        ! Create pressure and pressure thickness variables.
        allocate(aero_p(num_aero_lon,num_aero_lat,num_aero_lev))
        allocate(aero_dp(num_aero_lon,num_aero_lat,num_aero_lev))

        aero_time_step = 0

    end subroutine source_data_start_aerosol

    subroutine source_data_read_aerosol

        character(10) name
        integer ierr, start2d(3), count2d(3),  start3d(4), count3d(4)
        integer data_type, num_dim, dim_ids(4)
        integer i, j, k

        character(50), parameter :: sub_name = "source_data_read_aerosol"

        aero_time_step = aero_time_step+1

        start2d = [1,1,aero_time_step]
        count2d = [num_aero_lon,num_aero_lat,1]
        start3d = [1,1,1,aero_time_step]
        count3d = [num_aero_lon,num_aero_lat,num_aero_lev,1]

        do i = 1, num_aero_var
            if (.not. aero_var_read_each_time(i)) cycle
            ierr = nf90_inquire_variable(aero_ncid, aero_var_ids(i), name=name, &
                xtype=data_type, ndims=num_dim, dimids=dim_ids)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            select case (num_dim)
            case (3) ! 2D variable
                select type (var => data_aero_list%get_var(name))
                type is (var2d_d)
                    tmp_2d_d => var%get_values()
                    ierr = nf90_get_var(aero_ncid, aero_var_ids(i), tmp_2d_d, start2d, count2d)
                    call handle_netcdf_error(sub_name, __LINE__, ierr)
                end select
            case (4) ! 3D variable
                select type (var => data_aero_list%get_var(name))
                type is (var3d_d)
                    tmp_3d_d => var%get_values()
                    ierr = nf90_get_var(aero_ncid, aero_var_ids(i), tmp_3d_d, start3d, count3d)
                    call handle_netcdf_error(sub_name, __LINE__, ierr)
                end select
            end select
        end do

        ! Calculate each level's pressure.
        do k = 1, num_aero_lev
            do j = 1, num_aero_lat
                do i = 1, num_aero_lon
                    aero_p(i,j,k) = aero_hyam(k)*aero_p0+aero_hybm(k)*aero_ps(i,j)
                end do
            end do
        end do

        ! Calculate each level's pressure thickness.
        do k = 1, num_aero_lev
            do j = 1, num_aero_lat
                do i = 1, num_aero_lon
                    aero_dp(i,j,k) = aero_p0*(aero_hyai(k+1)-aero_hyai(k))+aero_ps(i,j)*(aero_hybi(k+1)-aero_hybi(k))
                end do
            end do
        end do

    end subroutine source_data_read_aerosol

    subroutine source_data_end_aerosol

        integer ierr

        character(50), parameter :: sub_name = "source_data_end_aerosol"

        ierr = nf90_close(aero_ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        call data_aero_list%final

        deallocate(aero_p)
        deallocate(aero_dp)

    end subroutine source_data_end_aerosol

end module source_data
