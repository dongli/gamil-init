module model_grids

    use utils
    use namelist_mod
    use variable
    use io_manager

    implicit none

    integer num_model_lev

    integer model_2d_dims(2), model_3d_dims(3)

    type(var_list) model_dims
    type(var_list) model_dims_aux

    ! Horizontal lat-lon mesh center grid coordinates
    real(8), pointer :: model_lon(:)
    real(8), pointer :: model_lat(:)
    real(8), pointer :: model_lon_bnds(:)
    real(8), pointer :: model_lat_bnds(:)

    integer, parameter :: classic_sigma_pressure = 1
    integer, parameter :: hybrid_sigma_pressure = 2
    integer :: model_vertical_coordinate_type = 0

    ! Classic sigma-pressure coordinate parameters
    real(8), pointer :: model_lev(:)
    real(8), pointer :: model_lev_bnds(:)
    real(8), pointer :: model_pt

    ! Hybrid sigma-pressure coordinate parameters
    real(8), pointer :: model_hyai(:)
    real(8), pointer :: model_hybi(:)
    real(8), pointer :: model_hyam(:)
    real(8), pointer :: model_hybm(:)
    real(8), pointer :: model_p0

    interface model_grids_calc_p
        module procedure model_grids_calc_p_sigma
        module procedure model_grids_calc_p_hybrid
    end interface model_grids_calc_p

contains

#include "calc_model_grids.F90"

    subroutine model_grids_read(file_name)

        character(*), intent(in) :: file_name

        integer ncid
        integer lon_dim_id, lat_dim_id, lev_dim_id
        integer lon_var_id, lat_var_id, lev_var_id, lev_bnds_var_id
        integer pmtop_var_id
        integer ierr

        character(*), parameter :: sub_name = "model_grids_read"

        call notice(sub_name, "Read GAMIL grid dimension information")

        ierr = nf90_open(file_name, nf90_nowrite, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lon", lon_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lat", lat_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_dimid(ncid, "lev", lev_dim_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(ncid, lon_dim_id, len=num_model_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(ncid, lat_dim_id, len=num_model_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inquire_dimension(ncid, lev_dim_id, len=num_model_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        call model_dims%append("lon", "longitude", "degrees_east", [num_model_lon])
        call model_dims%get_tail_values(model_lon)
        call model_dims%append("lat", "latitude", "degrees_north", [num_model_lat])
        call model_dims%get_tail_values(model_lat)
        call model_dims%append("lev", "level (sigma vertical coordinate)", "1", [num_model_lev])
        call model_dims%get_tail_values(model_lev)
        call model_dims%append("lev_bnds", "level bounds", "1", [num_model_lev+1])
        call model_dims%get_tail_values(model_lev_bnds)

        ierr = nf90_inq_varid(ncid, "lon", lon_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lat", lat_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lev", lev_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "lev_bnds", lev_bnds_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(ncid, "pmtop", pmtop_var_id)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(ncid, lon_var_id, model_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(ncid, lat_var_id, model_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(ncid, lev_var_id, model_lev)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(ncid, lev_bnds_var_id, model_lev_bnds)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_get_var(ncid, pmtop_var_id, model_pt)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        model_2d_dims = [num_model_lon,num_model_lat]
        model_3d_dims = [num_model_lon,num_model_lat,num_model_lev]

    end subroutine model_grids_read

    subroutine model_grids_write

        integer ncid, ierr, i, j, k
        integer grid_rank_dimid
        integer grid_size_dimid
        integer grid_corners_dimid
        integer grid_dims_varid
        integer grid_imask_varid
        integer grid_center_lon_varid
        integer grid_center_lat_varid
        integer grid_corner_lon_varid
        integer grid_corner_lat_varid

        integer, allocatable :: grid_imask(:)
        real(8), allocatable :: grid_center_lon(:)
        real(8), allocatable :: grid_center_lat(:)
        real(8), allocatable :: grid_corner_lon(:,:)
        real(8), allocatable :: grid_corner_lat(:,:)

        character(30) file_name
        character(*), parameter :: sub_name = "model_grids_write"

        file_name = "grid.gamil." // trim(to_string(num_model_lon)) // "x" // trim(to_string(num_model_lat)) // ".nc"

        ! Write out a SCRIP format data for generating mapping files.
        ierr = nf90_create(file_name, nf90_clobber, ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_dim(ncid, "grid_rank", 2, grid_rank_dimid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_dim(ncid, "grid_size", num_model_lon * num_model_lat, grid_size_dimid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_dim(ncid, "grid_corners", 4, grid_corners_dimid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_dims", nf90_int, [grid_rank_dimid], grid_dims_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_imask", nf90_int, [grid_size_dimid], grid_imask_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(ncid, grid_imask_varid, "units", "unitless")
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_center_lon", nf90_double, [grid_size_dimid], grid_center_lon_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(ncid, grid_center_lon_varid, "units", "degrees")
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_center_lat", nf90_double, [grid_size_dimid], grid_center_lat_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(ncid, grid_center_lat_varid, "units", "degrees")
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_corner_lon", nf90_double, [grid_corners_dimid,grid_size_dimid], grid_corner_lon_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(ncid, grid_corner_lon_varid, "units", "degrees")
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(ncid, "grid_corner_lat", nf90_double, [grid_corners_dimid,grid_size_dimid], grid_corner_lat_varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(ncid, grid_corner_lat_varid, "units", "degrees")
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_var(ncid, grid_dims_varid, [num_model_lon, num_model_lat])
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        allocate(grid_imask(num_model_lon * num_model_lat))
        allocate(grid_center_lon(num_model_lon * num_model_lat))
        allocate(grid_center_lat(num_model_lon * num_model_lat))
        allocate(grid_corner_lon(4,num_model_lon * num_model_lat))
        allocate(grid_corner_lat(4,num_model_lon * num_model_lat))

        grid_imask(:) = 1
        k = 1
        do j = 1, num_model_lat
            do i = 1, num_model_lon
                grid_center_lon(k) = model_lon(i)
                grid_center_lat(k) = model_lat(j)
                grid_corner_lon(:,k) = [model_lon_bnds(i),model_lon_bnds(i+1),model_lon_bnds(i+1),model_lon_bnds(i)]
                grid_corner_lat(:,k) = [model_lat_bnds(j),model_lat_bnds(j),model_lat_bnds(j+1),model_lat_bnds(j+1)]
                k = k + 1
            end do
        end do

        ierr = nf90_put_var(ncid, grid_imask_varid, grid_imask)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_var(ncid, grid_center_lon_varid, grid_center_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_var(ncid, grid_center_lat_varid, grid_center_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_var(ncid, grid_corner_lon_varid, grid_corner_lon)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_var(ncid, grid_corner_lat_varid, grid_corner_lat)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_close(ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        deallocate(grid_imask)
        deallocate(grid_center_lon)
        deallocate(grid_center_lat)
        deallocate(grid_corner_lon)
        deallocate(grid_corner_lat)

        call notice(sub_name, "File "//trim(file_name)//" has been generated")

    end subroutine model_grids_write

    subroutine model_grids_add_dims(file_idx)

        integer, intent(in) :: file_idx

        integer i
        class(var), pointer :: ptr

        ptr => model_dims%get_head()
        do i = 1, model_dims%get_num_var()
            call io_manager_add_dim(file_idx, ptr)
            ptr => ptr%next
        end do
        ptr => model_dims_aux%get_head()
        do i = 1, model_dims_aux%get_num_var()
            if (ptr%get_name() == "hyam" .or. ptr%get_name() == "hybm") then
                call io_manager_add_var(file_idx, ptr, ["lev"])
            else if (ptr%get_name() == "hyai" .or. ptr%get_name() == "hybi") then
                call io_manager_add_var(file_idx, ptr, ["ilev"])
            else
                call io_manager_add_var(file_idx, ptr)
            end if
            ptr => ptr%next
        end do

    end subroutine model_grids_add_dims

    subroutine model_grids_add_2d_dims(file_idx)

        integer, intent(in) :: file_idx

        integer i
        class(var), pointer :: ptr

        ptr => model_dims%get_head()
        do i = 1, model_dims%get_num_var()
            if (ptr%get_name() == "lon" .or. ptr%get_name() == "lat") then
                call io_manager_add_dim(file_idx, ptr)
            end if
            ptr => ptr%next
        end do

    end subroutine model_grids_add_2d_dims

    subroutine model_grids_calc_p_sigma(sig, ps, pt, p)

        real(8), intent(in) :: sig(:), ps(:,:), pt
        real(8), intent(out) :: p(:,:,:)

        integer dims(2), i, j, k
        integer num_lon, num_lat, num_lev

        dims = shape(ps)
        num_lev = size(sig)
        num_lon = dims(1)
        num_lat = dims(2)

        do k = 1, num_lev
            do j = 1, num_lat
                do i = 1, num_lon
                    p(i,j,k) = pt+sig(k)*(ps(i,j)-pt)
                end do
            end do
        end do

    end subroutine model_grids_calc_p_sigma

    subroutine model_grids_calc_p_hybrid(hya, hyb, ps, p0, p)

        real(8), intent(in) :: hya(:), hyb(:), ps(:,:), p0
        real(8), intent(out) :: p(:,:,:)

        integer dims(3), i, j, k
        integer num_lon, num_lat, num_lev

        dims = shape(p)
        num_lon = dims(1)
        num_lat = dims(2)
        num_lev = dims(3)

        do k = 1, num_lev
            do j = 1, num_lat
                do i = 1, num_lon
                    p(i,j,k) = hya(k)*p0+hyb(k)*ps(i,j)
                end do
            end do
        end do

    end subroutine model_grids_calc_p_hybrid

end module model_grids
