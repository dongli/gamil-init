module utils

    use constants
    use netcdf

    interface to_string
        module procedure integer_to_string
        module procedure real_to_string
    end interface to_string

contains

    subroutine notice(sub_name, message)

        character(*), intent(in) :: sub_name
        character(*), intent(in) :: message

        write(*, "('[Notice]: ', A, ': ', A)") trim(sub_name), trim(message)

    end subroutine notice

    subroutine report_error(sub_name, message)
    
        character(*), intent(in) :: sub_name
        character(*), intent(in) :: message

        write(*, "('[Error]: ', A, ': ', A)") trim(sub_name), trim(message)
        stop

    end subroutine report_error

    subroutine report_warning(sub_name, message)
    
        character(*), intent(in) :: sub_name
        character(*), intent(in) :: message

        write(*, "('[Warning]: ', A, ': ', A)") trim(sub_name), trim(message)

    end subroutine report_warning

    subroutine handle_netcdf_error(sub_name, line, ierr)

        character(*), intent(in) :: sub_name
        integer, intent(in) :: line
        integer, intent(in) :: ierr

        if (ierr /= nf90_noerr) then
            write(*, "('[Error]: ', A, ':', I3, ': ', A)") &
                trim(sub_name), line, trim(nf90_strerror(ierr))//"."
            stop
        end if

    end subroutine handle_netcdf_error

    subroutine check_pressure_units(sub_name, p, units)

        character(*), intent(in) :: sub_name
        real(8), intent(inout) :: p(:)
        character(*), intent(inout) :: units

        if (units == "hPa" .or. units == "millibars") then
            call notice(sub_name, "Change the pressure units from "// &
                trim(units)//" to Pa")
            units = "Pa"
            p = p*100.0
        end if

    end subroutine check_pressure_units

    subroutine check_file_exist(sub_name, file_name)

        character(*), intent(in) :: sub_name
        character(*), intent(in) :: file_name

        logical alive

        inquire(file=file_name, exist=alive)
        if (.not. alive) then
            call report_error(sub_name, &
                "File "//trim(file_name)//" does not exist!")
        end if

    end subroutine check_file_exist

    ! --------------------------------------------------------------------------
    ! string procedures
    integer function count_substr(string, substr) result (res)

        character(*), intent(in) :: string, substr

        integer i, j, k, n

        i = 1
        j = 0
        n = len_trim(string)

        res = 0
        do
            k = index(string(i:n), substr)
            j = j+k
            if (k /= 0) then
                i = j+1
                res = res+1
            else
                return
            end if
        end do

    end function count_substr

    character(100) function get_substr(string, delim, idx) result (res)

        character(*), intent(in) :: string, delim
        integer, intent(in) :: idx

        integer i1, i2, j, c, c0, n

        c0 = count_substr(string, delim)
        if (idx > c0+1) then
            res = "[Error]: get_substr: Index exceeds!"
            return
        end if

        i1 = 1; i2 = 1; j = 0; c = 0
        n = len_trim(string)

        do
            j = j+index(string(i2:n), delim)
            i2 = j+1
            c = c+1
            if (c == idx) then
                res = string(i1:i2-2)
                return
            else if (c == c0) then
                res = string(i2:n)
                return
            else
                i1 = i2
            end if
        end do

    end function get_substr

    character(100) function integer_to_string(n, format) result (res)

        integer, intent(in) :: n
        character(*), intent(in), optional :: format

        if (present(format)) then
            write(res, "("//format//")") n
        else
            write(res, "(I100)") n
        end if

        res = adjustl(res)

    end function integer_to_string

    character(100) function real_to_string(n, format) result (res)

        real(8), intent(in) :: n
        character(*), intent(in) :: format

        write(res, format) n

        res = adjustl(res)

    end function real_to_string

end module utils
