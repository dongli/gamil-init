module io_manager

    use utils
    use variable

    implicit none

    integer, parameter :: num_max_file = 100
    integer, parameter :: num_max_dim = 10
    integer, parameter :: num_max_var = 1000

    integer :: num_file = 0

    type dim_ptr
        type(var1d_d), pointer :: ptr
    end type dim_ptr

    type file_info
        character(300) name
        integer ncid
        integer :: num_dim = 0
        integer dimids(num_max_dim)
        integer :: num_var = 0
        integer var_ids(num_max_var)
        type(dim_ptr) dims(num_max_dim)
    end type file_info

    type(file_info), target :: file_infos(num_max_file)

    interface io_manager_put_var
        module procedure io_manager_put_var_good
        module procedure io_manager_put_var_int_0d
        module procedure io_manager_put_var_double_0d
        module procedure io_manager_put_var_float_0d
    end interface io_manager_put_var

    interface io_manager_add_att
        module procedure io_manager_add_var_att
        module procedure io_manager_add_file_att
    end interface io_manager_add_att

contains

    subroutine io_manager_create_file(file_name, file_idx)

        character(*), intent(in) :: file_name
        integer, intent(out) :: file_idx

        type(file_info), pointer :: file
        integer ierr

        character(50), parameter :: sub_name = "io_manager_create_file"

        if (num_file > num_max_file) then
            call report_error(sub_name, "Exceed maximum file number!")
        end if

        num_file = num_file+1
        file => file_infos(num_file)
        file%name = file_name

        ierr = nf90_create(file_name, nf90_clobber, file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        file_idx = num_file

        call notice(sub_name, "File "//trim(file_name)//" is created.")

    end subroutine io_manager_create_file

    subroutine io_manager_add_dim(file_idx, dim)

        integer, intent(in) :: file_idx
        class(var), intent(in), target :: dim

        type(file_info), pointer :: file
        integer, pointer :: dimid
        integer varid, ierr

        character(50), parameter :: sub_name = "io_manager_add_dim"

        call get_file_info(sub_name, file_idx, file)

        file%num_dim = file%num_dim+1
        if (file%num_dim > num_max_dim) then
            call report_error(sub_name, "Exceed maximum dimension number!")
        end if
        dimid => file%dimids(file%num_dim)
        select type (dim)
        type is (var1d_d)
            file%dims(file%num_dim)%ptr => dim
        class default
            call report_error(sub_name, "Dimension variable must be 1D!")
        end select

        ierr = nf90_redef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_dim(file%ncid, dim%get_name(), dim%get_dims(), dimid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_def_var(file%ncid, dim%get_name(), nf90_double, [dimid], varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, varid, "long_name", dim%get_long_name())
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, varid, "units", dim%get_units())
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        select type (dim)
        type is (var1d_d)
            ierr = nf90_put_var(file%ncid, varid, dim%get_values())
            call handle_netcdf_error(sub_name, __LINE__, ierr)
        end select

    end subroutine io_manager_add_dim

    subroutine io_manager_def_var(file_idx, vari, dim_names)

        integer, intent(in) :: file_idx
        class(var), intent(in) :: vari
        character(*), intent(in), optional :: dim_names(:)

        type(file_info), pointer :: file
        integer dimids(4), i
        integer num_dim, ierr
        integer, pointer :: varid
        integer data_type

        character(50), parameter :: sub_name = "io_manager_def_var"

        call get_file_info(sub_name, file_idx, file)

        file%num_var = file%num_var+1
        if (file%num_var > num_max_var) then
            call report_error(sub_name, "Exceed maximum variable number!")
        end if
        varid => file%var_ids(file%num_var)

        if (present(dim_names)) then
            num_dim = size(dim_names)
            dimids = -1
            do i = 1, num_dim
                ierr = nf90_inq_dimid(file%ncid, dim_names(i), dimids(i))
                call handle_netcdf_error(sub_name, __LINE__, ierr)
            end do
        end if

        ierr = nf90_redef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        select case (vari%get_data_type())
        case ("double")
            data_type = nf90_double
        case ("integer")
            data_type = nf90_int
        end select

        if (present(dim_names)) then
            ierr = nf90_def_var(file%ncid, vari%get_name(), data_type, &
                dimids(1:num_dim), varid)
        else
            ierr = nf90_def_var(file%ncid, vari%get_name(), data_type, varid)
        end if
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, varid, "long_name", vari%get_long_name())
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, varid, "units", vari%get_units())
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_def_var

    subroutine io_manager_put_var_good(file_idx, vari, rec)

        integer, intent(in) :: file_idx
        class(var), intent(in) :: vari
        integer, intent(in), optional :: rec

        type(file_info), pointer :: file
        integer ierr, varid, num_dim, dimids(4)
        integer start(4), count(4), i

        character(50), parameter :: sub_name = "io_manager_put_var_good"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_inq_varid(file%ncid, vari%get_name(), varid)

        dimids = -1
        ierr = nf90_inquire_variable(file%ncid, varid, ndims=num_dim, dimids=dimids)

        start = 1
        do i = 1, num_dim
            ierr = nf90_inquire_dimension(file%ncid, dimids(i), len=count(i))
        end do
        if (present(rec)) then
            if (num_dim /= 0) then
                start(num_dim) = rec
                count(num_dim) = 1
            end if
        end if

        select type (vari)
        type is (var0d_i)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start)
        type is (var0d_d)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start)
        type is (var1d_d)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start, count)
        type is (var1d_i)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start, count)
        type is (var2d_d)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start, count)
        type is (var3d_d)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start, count)
        type is (var4d_d)
            ierr = nf90_put_var(file%ncid, varid, vari%get_values(), start, count)
        end select

    end subroutine io_manager_put_var_good

    subroutine io_manager_add_var(file_idx, vari, dim_names)

        integer, intent(in) :: file_idx
        class(var), intent(in) :: vari
        character(*), intent(in), optional :: dim_names(:)

        if (present(dim_names)) then
            call io_manager_def_var(file_idx, vari, dim_names)
        else
            call io_manager_def_var(file_idx, vari)
        end if

        call io_manager_put_var(file_idx, vari)

    end subroutine io_manager_add_var

    subroutine io_manager_close_file(file_idx)

        integer, intent(in) :: file_idx

        type(file_info), pointer :: file
        integer ierr

        character(50), parameter :: sub_name = "io_manager_close_file"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_close(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_close_file

    subroutine io_manager_def_dim(file_idx, dim_name, var_type, &
        dim_size, long_name, units)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: dim_name
        character(*), intent(in) :: var_type
        integer, intent(in), optional :: dim_size
        character(*), intent(in), optional :: long_name, units

        type(file_info), pointer :: file
        integer, pointer :: dimid
        integer varid, ierr

        character(50), parameter :: sub_name = "io_manager_def_dim"

        call get_file_info(sub_name, file_idx, file)

        file%num_dim = file%num_dim+1
        if (file%num_dim > num_max_dim) then
            call report_error(sub_name, "Exceed maximum dimension number!")
        end if
        dimid => file%dimids(file%num_dim)

        ierr = nf90_redef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        if (present(dim_size)) then
            ierr = nf90_def_dim(file%ncid, dim_name, dim_size, dimid)
        else
            ierr = nf90_def_dim(file%ncid, dim_name, nf90_unlimited, dimid)
        end if
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        select case (var_type)
        case ("integer")
            ierr = nf90_def_var(file%ncid, dim_name, nf90_int, dimid, varid)
        case ("real(8)")
            ierr = nf90_def_var(file%ncid, dim_name, nf90_double, dimid, varid)
        case ("real(4)")
            ierr = nf90_def_var(file%ncid, dim_name, nf90_float, dimid, varid)
        end select
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        if (present(long_name)) then
            ierr = nf90_put_att(file%ncid, varid, "long_name", long_name)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
        end if

        if (present(units)) then
            ierr = nf90_put_att(file%ncid, varid, "units", units)
            call handle_netcdf_error(sub_name, __LINE__, ierr)
        end if

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_def_dim

    subroutine io_manager_put_var_int_0d(file_idx, var_name, value, rec)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: var_name
        integer, intent(in) :: value
        integer, intent(in), optional :: rec

        type(file_info), pointer :: file
        integer varid, ierr

        character(50), parameter :: sub_name = "io_manager_put_var_int_0d"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_inq_varid(file%ncid, var_name, varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        if (present(rec)) then
            ierr = nf90_put_var(file%ncid, varid, value, [rec])
        else
            ierr = nf90_put_var(file%ncid, varid, value)
        end if
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_put_var_int_0d

    subroutine io_manager_put_var_double_0d(file_idx, var_name, value, rec)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: var_name
        real(8), intent(in) :: value
        integer, intent(in), optional :: rec

        type(file_info), pointer :: file
        integer varid, ierr

        character(50), parameter :: sub_name = "io_manager_put_var_double_0d"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_inq_varid(file%ncid, var_name, varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        if (present(rec)) then
            ierr = nf90_put_var(file%ncid, varid, value, [rec])
        else
            ierr = nf90_put_var(file%ncid, varid, value)
        end if
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_put_var_double_0d

    subroutine io_manager_put_var_float_0d(file_idx, var_name, value, rec)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: var_name
        real(4), intent(in) :: value
        integer, intent(in), optional :: rec

        type(file_info), pointer :: file
        integer varid, ierr

        character(50), parameter :: sub_name = "io_manager_put_var_float_0d"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_inq_varid(file%ncid, var_name, varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        if (present(rec)) then
            ierr = nf90_put_var(file%ncid, varid, value, [rec])
        else
            ierr = nf90_put_var(file%ncid, varid, value)
        end if
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_put_var_float_0d

    subroutine io_manager_add_var_att(file_idx, var_name, att_name, value)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: var_name, att_name, value

        type(file_info), pointer :: file
        integer ierr, varid

        character(50), parameter :: sub_name = "io_manager_add_var_att"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_redef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_inq_varid(file%ncid, var_name, varid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, varid, att_name, value)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_add_var_att

    subroutine io_manager_add_file_att(file_idx, att_name, value)

        integer, intent(in) :: file_idx
        character(*), intent(in) :: att_name, value

        type(file_info), pointer :: file
        integer ierr

        character(50), parameter :: sub_name = "io_manager_add_file_att"

        call get_file_info(sub_name, file_idx, file)

        ierr = nf90_redef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_put_att(file%ncid, nf90_global, att_name, value)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

        ierr = nf90_enddef(file%ncid)
        call handle_netcdf_error(sub_name, __LINE__, ierr)

    end subroutine io_manager_add_file_att

    subroutine get_file_info(sub_name, file_idx, file)

        character(*), intent(in) :: sub_name
        integer, intent(in) :: file_idx
        type(file_info), pointer :: file

        if (file_idx > num_max_file) then
            call report_error(sub_name, "Exceed maximum file number!")
        end if
        file => file_infos(file_idx)

    end subroutine get_file_info

end module io_manager
