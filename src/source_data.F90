module source_data

    use utils
    use variable

    implicit none

    ! --------------------------------------------------------------------------
    ! u, v, t, q data set
    integer, parameter :: era_interim = 1
    integer, parameter :: model_data  = 2
    integer :: uvtq_data_type = -1
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
    ! ozone data set
    integer num_ozone_lon, num_ozone_lat, num_ozone_lev, num_ozone_time
    real(8), allocatable :: data_ozone_lon(:)
    real(8), allocatable :: data_ozone_lat(:)
    real(8), allocatable :: data_ozone_lev(:)
    real(8), allocatable :: data_ozone(:,:,:,:)

    ! --------------------------------------------------------------------------
    ! aerosol data set
    integer aero_ncid, aero_time_step, aero_ps_var_id
    integer, allocatable :: aero_var_id(:)
    integer num_aero_lon, num_aero_lat, num_aero_lev, num_aero_time
    real(8) aero_p0
    real(8), allocatable :: aero_lon(:)
    real(8), allocatable :: aero_lat(:)
    real(8), allocatable :: aero_hyam(:)
    real(8), allocatable :: aero_hybm(:)
    real(8), allocatable :: aero_hyai(:)
    real(8), allocatable :: aero_hybi(:)
    integer, allocatable :: aero_time(:)
    integer, allocatable :: aero_date(:)
    real(8), allocatable :: aero_ps(:,:)
    real(8), allocatable :: aero_p(:,:,:)
    real(8), allocatable :: aero_dp(:,:,:)
    type(var_list) data_aero_list

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

        integer ierr
        integer time_dim_id, lon_dim_id, lat_dim_id, lev_dim_id
        integer time_var_id, date_var_id, lon_var_id, lat_var_id
        integer p0_var_id, hyam_var_id, hybm_var_id, hyai_var_id, hybi_var_id

        integer, parameter :: num_aero_var = 13
        character(6) :: aero_var_names(num_aero_var) = ["SO4   ", &
            "SSLT01","SSLT02","SSLT03","SSLT04", &
            "DST01 ","DST02 ","DST03 ","DST04 ", &
            "OC1   ","CB1   ","OC2   ","CB2   "]
        character(30) long_name, units
        integer dims2d(2), dims3d(3), i

        character(50), parameter :: sub_name = "source_data_start_aerosol"

        call notice(sub_name, "Start read aerosol file """//trim(file_name)//"""")

        ierr = nf90_open(file_name, nf90_nowrite, aero_ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(aero_ncid, "time", time_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(aero_ncid, time_dim_id, len=num_aero_time)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(aero_ncid, "lon", lon_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(aero_ncid, lon_dim_id, len=num_aero_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(aero_ncid, "lat", lat_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(aero_ncid, lat_dim_id, len=num_aero_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(aero_ncid, "lev", lev_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(aero_ncid, lev_dim_id, len=num_aero_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        allocate(aero_time(num_aero_time))
        allocate(aero_date(num_aero_time))
        allocate(aero_lon(num_aero_lon))
        allocate(aero_lat(num_aero_lat))
        allocate(aero_hyam(num_aero_lev))
        allocate(aero_hybm(num_aero_lev))
        allocate(aero_hyai(num_aero_lev+1))
        allocate(aero_hybi(num_aero_lev+1))

        ierr = nf90_inq_varid(aero_ncid, "time", time_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, time_var_id, aero_time)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "date", date_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, date_var_id, aero_date)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "lon", lon_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, lon_var_id, aero_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "lat", lat_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, lat_var_id, aero_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "P0", p0_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, p0_var_id, aero_p0)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "hyam", hyam_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, hyam_var_id, aero_hyam)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "hybm", hybm_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, hybm_var_id, aero_hybm)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "hyai", hyai_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, hyai_var_id, aero_hyai)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(aero_ncid, "hybi", hybi_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(aero_ncid, hybi_var_id, aero_hybi)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! ----------------------------------------------------------------------
        ! inquire surface pressure id
        allocate(aero_ps(num_aero_lon,num_aero_lat))
        allocate(aero_p(num_aero_lon,num_aero_lat,num_aero_lev))
        allocate(aero_dp(num_aero_lon,num_aero_lat,num_aero_lev))

        ierr = nf90_inq_varid(aero_ncid, "PS", aero_ps_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ! ----------------------------------------------------------------------
        ! inquire aerosol varible id and create aerosol variable list
        allocate(aero_var_id(num_aero_var))
        dims2d = [num_aero_lon,num_aero_lat]
        dims3d = [num_aero_lon,num_aero_lat,num_aero_lev]
        do i = 1, num_aero_var
            ierr = nf90_inq_varid(aero_ncid, aero_var_names(i), aero_var_id(i))
            call handle_netcdf_error(sub_name, __LINE__, ierr)

            ierr = nf90_get_att(aero_ncid, aero_var_id(i), "long_name", long_name)
            call handle_netcdf_error(sub_name, __LINE__, ierr)

            ierr = nf90_get_att(aero_ncid, aero_var_id(i), "units", units)
            call handle_netcdf_error(sub_name, __LINE__, ierr)

            call data_aero_list%append(aero_var_names(i), long_name, units, dims3d)
        end do

        aero_time_step = 0

    end subroutine source_data_start_aerosol

    subroutine source_data_read_aerosol

        integer ierr, start2d(3), count2d(3),  start3d(4), count3d(4)

        integer i, j, k
        real(8), pointer :: tmp(:,:,:)
        class(var), pointer :: ptr

        integer ncid
        integer lon_dim_id, lat_dim_id, lev_dim_id
        integer lon_var_id, lat_var_id, p_var_id

        character(50), parameter :: sub_name = "source_data_read_aerosol"

        aero_time_step = aero_time_step+1

        start2d = [1,1,aero_time_step]
        count2d = [num_aero_lon,num_aero_lat,1]
        start3d = [1,1,1,aero_time_step]
        count3d = [num_aero_lon,num_aero_lat,num_aero_lev,1]

        ierr = nf90_get_var(aero_ncid, aero_ps_var_id, aero_ps, start2d, count2d)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ptr => data_aero_list%get_head()
        do i = 1, data_aero_list%get_num_var()
            select type (ptr)
            type is (var3d)
                tmp => ptr%get_values()
            end select
            ierr = nf90_get_var(aero_ncid, aero_var_id(i), tmp, start3d, count3d)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
            ptr => ptr%next
        end do

        ! calculate each level's pressure
        do k = 1, num_aero_lev
            do j = 1, num_aero_lat
                do i = 1, num_aero_lon
                    aero_p(i,j,k) = aero_hyam(k)*aero_p0+aero_hybm(k)*aero_ps(i,j)
                end do
            end do
        end do

        ! calculate each level's pressure thickness
        do k = 1, num_aero_lev
            do j = 1, num_aero_lat
                do i = 1, num_aero_lon
                    aero_dp(i,j,k) = aero_ps(i,j)*(aero_hybi(k+1)-aero_hybi(k))
                end do
            end do
        end do

        ierr = nf90_create("check_p.nc", nf90_clobber, ncid)
        ierr = nf90_def_dim(ncid, "lon", num_aero_lon, lon_dim_id)
        ierr = nf90_def_dim(ncid, "lat", num_aero_lat, lat_dim_id)
        ierr = nf90_def_dim(ncid, "lev", num_aero_lev, lev_dim_id)
        ierr = nf90_def_var(ncid, "lon", nf90_double, lon_dim_id, lon_var_id)
        ierr = nf90_def_var(ncid, "lat", nf90_double, lat_dim_id, lat_var_id)
        ierr = nf90_def_var(ncid, "p", nf90_double, [lon_dim_id,lat_dim_id,lev_dim_id], p_var_id)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, lon_var_id, aero_lon)
        ierr = nf90_put_var(ncid, lat_var_id, aero_lat)
        ierr = nf90_put_var(ncid, p_var_id, aero_p)
        ierr = nf90_close(ncid)

    end subroutine source_data_read_aerosol

    subroutine source_data_end_aerosol

        integer ierr

        character(50), parameter :: sub_name = "source_data_end_aerosol"

        ierr = nf90_close(aero_ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine source_data_end_aerosol

end module source_data
