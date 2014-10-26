module model_grids

    use utils
    use variable
    use io_manager

    implicit none

    integer, parameter :: equal_interval_grid = 1
    integer, parameter :: even_area_grid = 2
    integer :: model_grid_type = even_area_grid

    integer num_model_lon, num_model_lat
    integer num_model_lev

    integer model_2d_dims(2), model_3d_dims(3)

    type(var_list) model_dims
    type(var_list) model_dims_aux

    real(8), pointer :: model_lon(:)
    real(8), pointer :: model_lat(:)
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

        character(30), parameter :: sub_name = "model_grids_read"

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

    subroutine model_grids_write(file_name)

        character(*), intent(in) :: file_name

        integer file_idx

        call io_manager_create_file(file_name, file_idx)
        call model_grids_add_dims(file_idx)
        call io_manager_close_file(file_idx)

    end subroutine model_grids_write

    subroutine model_grids_add_dims(file_idx)

        integer, intent(in) :: file_idx

        integer i
        class(var), pointer :: ptr

        ! Set lev and ilev from the equivalent hybrid coefficients.
        model_lev = model_hybm
        model_lev_bnds = model_hybi

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
        ptr => model_dims_aux%get_head()
        do i = 1, model_dims_aux%get_num_var()
            call io_manager_add_var(file_idx, ptr)
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
                    p(i,j,k) = hya(k)*ps(i,j)+hyb(k)*p0
                end do
            end do
        end do

    end subroutine model_grids_calc_p_hybrid

end module model_grids
