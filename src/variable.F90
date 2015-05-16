! ------------------------------------------------------------------------------
! Description:
!
!   This module defines "var" type, which is used to bind variable and its
!   metadata together.
!
! Author:
!
!   - Li Dong <dongli@lasg.iap.ac.cn>
! ------------------------------------------------------------------------------

module variable

    implicit none

    private

    public var, var_list
    public var0d_d, var1d_d, var2d_d, var3d_d, var4d_d
    public var0d_i, var1d_i

    ! --------------------------------------------------------------------------
    integer, parameter :: name_str_len = 8
    integer, parameter :: long_name_str_len = 255
    integer, parameter :: units_str_len = 30
    integer, parameter :: data_type_str_len = 10

    ! --------------------------------------------------------------------------
    type var
        character(name_str_len), private :: name
        character(long_name_str_len), private :: long_name
        character(units_str_len), private :: units
        class(var), pointer :: prev, next
    contains
        procedure :: init => var_init
        procedure :: get_name => var_get_name
        procedure :: get_long_name => var_get_long_name
        procedure :: get_units => var_get_units
        procedure :: get_dims => var_get_dims
        procedure :: get_data_type => var_get_data_type
        procedure :: reshape => var_reshape
        ! ======================================================================
        procedure :: var_set_values_1d
        generic :: set_values => var_set_values_1d
    end type var

    type, extends(var) :: var0d_d
        real(8), pointer, private :: values
    contains
        procedure :: get_values => var0d_d_get_values
        procedure :: get_data_type => var0d_d_get_data_type
    end type var0d_d

    type, extends(var) :: var0d_i
        integer, pointer, private :: values
    contains
        procedure :: get_values => var0d_i_get_values
        procedure :: get_data_type => var0d_i_get_data_type
    end type var0d_i

    type, extends(var) :: var1d_d
        integer, private :: dims(1)
        real(8), pointer, private :: values(:)
    contains
        procedure :: get_values => var1d_d_get_values
        procedure :: get_data_type => var1d_d_get_data_type
    end type var1d_d

    type, extends(var) :: var1d_i
        integer, private :: dims(1)
        integer, pointer, private :: values(:)
    contains
        procedure :: get_values => var1d_i_get_values
        procedure :: get_data_type => var1d_i_get_data_type
    end type var1d_i

    type, extends(var) :: var2d_d
        integer, private :: dims(2)
        real(8), pointer, private :: values(:,:)
    contains
        procedure :: get_values => var2d_d_get_values
        procedure :: get_data_type => var2d_d_get_data_type
    end type var2d_d

    type, extends(var) :: var3d_d
        integer, private :: dims(3)
        real(8), pointer, private :: values(:,:,:)
    contains
        procedure :: get_values => var3d_d_get_values
        procedure :: get_data_type => var3d_d_get_data_type
    end type var3d_d

    type, extends(var) :: var4d_d
        integer, private :: dims(4)
        real(8), pointer, private :: values(:,:,:,:)
    contains
        procedure :: get_values => var4d_d_get_values
        procedure :: get_data_type => var4d_d_get_data_type
    end type var4d_d

    ! --------------------------------------------------------------------------
    type var_list
        class(var), pointer, private :: head => null()
        class(var), pointer, private :: tail => null()
        integer, private :: num_var
    contains
        procedure :: append => var_list_append
        procedure :: get_num_var => var_list_get_num_var
        procedure :: get_var => var_list_get_var
        procedure :: get_head => var_list_get_head
        procedure :: get_tail => var_list_get_tail
        ! ======================================================================
        procedure :: var_list_get_tail_values_0d_d
        procedure :: var_list_get_tail_values_0d_i
        procedure :: var_list_get_tail_values_1d_d
        procedure :: var_list_get_tail_values_1d_i
        procedure :: var_list_get_tail_values_2d_d
        procedure :: var_list_get_tail_values_3d_d
        procedure :: var_list_get_tail_values_4d_d
        generic :: get_tail_values => var_list_get_tail_values_0d_d, &
                                      var_list_get_tail_values_0d_i, &
                                      var_list_get_tail_values_1d_d, &
                                      var_list_get_tail_values_1d_i, &
                                      var_list_get_tail_values_2d_d, &
                                      var_list_get_tail_values_3d_d, &
                                      var_list_get_tail_values_4d_d
    end type var_list

contains

    subroutine var_init(this, name, long_name, units, dims)

        class(var), intent(out) :: this
        character(*), intent(in) :: name, long_name, units
        integer, intent(in), optional :: dims(:)

        this%name = name
        this%long_name = long_name
        this%units = units

        if (present(dims)) then
            if (size(dims) > 4) then
                write(*, "('[Error]: var_init: Dimension number exceeds 4!')")
                stop
            end if
        end if

        select type (this)
        type is (var0d_d)
            allocate (this%values)
        type is (var0d_i)
            allocate (this%values)
        type is (var1d_d)
            if (size(dims) /= 1) then
                write(*, "('[Error]: var_init: Dimension number /= 1!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1)))
        type is (var1d_i)
            if (size(dims) /= 1) then
                write(*, "('[Error]: var_init: Dimension number /= 1!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1)))
        type is (var2d_d)
            if (size(dims) /= 2) then
                write(*, "('[Error]: var_init: Dimension number /= 2!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1),dims(2)))
        type is (var3d_d)
            if (size(dims) /= 3) then
                write(*, "('[Error]: var_init: Dimension number /= 3!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1),dims(2),dims(3)))
        type is (var4d_d)
            if (size(dims) /= 4) then
                write(*, "('[Error]: var_init: Dimension number /= 4!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1),dims(2),dims(3),dims(4)))
        end select

    end subroutine var_init

    character(name_str_len) function var_get_name(this) result (name)

        class(var), intent(in) :: this

        name = this%name

    end function var_get_name

    character(long_name_str_len) function var_get_long_name(this) result (long_name)

        class(var), intent(in) :: this

        long_name = this%long_name

    end function var_get_long_name

    character(units_str_len) function var_get_units(this) result (units)

        class(var), intent(in) :: this

        units = this%units

    end function var_get_units

    integer function var_get_dims(this, dim_idx) result (res)

        class(var), intent(in) :: this
        integer, intent(in), optional :: dim_idx

        select type (this)
        type is (var0d_d)
            res = 1
        type is (var0d_i)
            res = 1
        type is (var1d_d)
            res = this%dims(1)
        type is (var1d_i)
            res = this%dims(1)
        type is (var2d_d)
            res = this%dims(dim_idx)
        type is (var3d_d)
            res = this%dims(dim_idx)
        type is (var4d_d)
            res = this%dims(dim_idx)
        end select

    end function var_get_dims

    function var_get_data_type(this) result (res)

        class(var), intent(in) :: this
        character(data_type_str_len) res

        select type (this)
        type is (var0d_d)
            res = this%get_data_type()
        type is (var0d_i)
            res = this%get_data_type()
        type is (var1d_d)
            res = this%get_data_type()
        type is (var1d_i)
            res = this%get_data_type()
        type is (var2d_d)
            res = this%get_data_type()
        type is (var3d_d)
            res = this%get_data_type()
        type is (var4d_d)
            res = this%get_data_type()
        end select

    end function var_get_data_type

    subroutine var_reshape(this, order)

        class(var), intent(inout) :: this
        integer, intent(in) :: order(:)

        integer dims2d(2), dims3d(3), dims4d(4)
        real(8), allocatable :: tmp2d(:,:), tmp3d(:,:,:), tmp4d(:,:,:,:)

        select type (this)
        type is (var3d_d)
            dims3d(1) = this%get_dims(order(1))
            dims3d(2) = this%get_dims(order(2))
            dims3d(3) = this%get_dims(order(3))
            allocate (tmp3d(dims3d(1),dims3d(2),dims3d(3)))
            tmp3d = reshape(this%values, dims3d, order=order)
            deallocate (this%values)
            allocate (this%values(dims3d(1),dims3d(2),dims3d(3)))
            this%values = tmp3d
        end select

    end subroutine var_reshape

    function var0d_d_get_values(this) result (res)

        class(var0d_d), intent(in) :: this
        real(8), pointer :: res

        res => this%values

    end function var0d_d_get_values

    function var0d_d_get_data_type(this) result (res)

        class(var0d_d), intent(in) :: this
        character(data_type_str_len) res

        res = "double"

    end function var0d_d_get_data_type

    function var0d_i_get_values(this) result (res)

        class(var0d_i), intent(in) :: this
        integer, pointer :: res

        res => this%values

    end function var0d_i_get_values

    function var0d_i_get_data_type(this) result (res)

        class(var0d_i), intent(in) :: this
        character(data_type_str_len) res

        res = "integer"

    end function var0d_i_get_data_type

    function var1d_d_get_values(this) result (res)

        class(var1d_d), intent(in) :: this
        real(8), pointer :: res(:)

        res => this%values

    end function var1d_d_get_values

    function var1d_d_get_data_type(this) result (res)

        class(var1d_d), intent(in) :: this
        character(data_type_str_len) res

        res = "double"

    end function var1d_d_get_data_type

    subroutine var_set_values_1d(this, values)

        class(var), intent(inout) :: this
        real(8), intent(in) :: values(:)

        select type (this)
        type is (var1d_d)
            this%values = values
        class default
            write(*, "('[Error]: var_set_values_1d: Variable is not 1d!')")
            stop
        end select

    end subroutine var_set_values_1d

    function var1d_i_get_values(this) result (res)

        class(var1d_i), intent(in) :: this
        integer, pointer :: res(:)

        res => this%values

    end function var1d_i_get_values

    function var1d_i_get_data_type(this) result (res)

        class(var1d_i), intent(in) :: this
        character(data_type_str_len) res

        res = "integer"

    end function var1d_i_get_data_type

    function var2d_d_get_values(this) result (res)

        class(var2d_d), intent(in) :: this
        real(8), pointer :: res(:,:)

        res => this%values

    end function var2d_d_get_values

    function var2d_d_get_data_type(this) result (res)

        class(var2d_d), intent(in) :: this
        character(data_type_str_len) res

        res = "double"

    end function var2d_d_get_data_type

    function var3d_d_get_values(this) result (res)

        class(var3d_d), intent(in) :: this
        real(8), pointer :: res(:,:,:)

        res => this%values

    end function var3d_d_get_values

    function var3d_d_get_data_type(this) result (res)

        class(var3d_d), intent(in) :: this
        character(data_type_str_len) res

        res = "double"

    end function var3d_d_get_data_type

    function var4d_d_get_values(this) result (res)

        class(var4d_d), intent(in) :: this
        real(8), pointer :: res(:,:,:,:)

        res => this%values

    end function var4d_d_get_values

    function var4d_d_get_data_type(this) result (res)

        class(var4d_d), intent(in) :: this
        character(data_type_str_len) res

        res = "double"

    end function var4d_d_get_data_type

    subroutine var_list_append(this, name, long_name, units, dims, pointer, data_type)

        class(var_list), intent(inout) :: this
        character(*), intent(in) :: name, long_name, units
        integer, intent(in), optional :: dims(:)
        class(var), intent(out), optional, pointer :: pointer
        character(*), intent(in), optional :: data_type

        class(var), pointer :: ptr
        character(data_type_str_len) data_type_

        type(var0d_d) var0d_d_type
        type(var0d_i) var0d_i_type
        type(var1d_d) var1d_d_type
        type(var1d_i) var1d_i_type
        type(var2d_d) var2d_d_type
        type(var3d_d) var3d_d_type
        type(var4d_d) var4d_d_type

        if (present(data_type)) then
            data_type_ = data_type
        else
            data_type_ = "double"
        end if

        if (.not. associated(this%head)) then
            if (.not. present(dims)) then
                select case (data_type_)
                case ("double")
                    allocate (this%head, source=var0d_d_type)
                case ("integer")
                    allocate (this%head, source=var0d_i_type)
                end select
            else
                select case (size(dims))
                case (1)
                    select case (data_type_)
                    case ("double")
                        allocate (this%head, source=var1d_d_type)
                    case ("integer")
                        allocate (this%head, source=var1d_i_type)
                    end select
                case (2)
                    allocate (this%head, source=var2d_d_type)
                case (3)
                    allocate (this%head, source=var3d_d_type)
                case (4)
                    allocate (this%head, source=var4d_d_type)
                end select
            end if
            this%tail => this%head
            nullify(this%head%prev)
            nullify(this%head%next)
            this%num_var = 1
        else
            if (.not. present(dims)) then
                select case (data_type_)
                case ("double")
                    allocate (this%tail%next, source=var0d_d_type)
                case ("integer")
                    allocate (this%tail%next, source=var0d_i_type)
                end select
            else
                select case (size(dims))
                case (1)
                    select case (data_type_)
                    case ("double")
                        allocate (this%tail%next, source=var1d_d_type)
                    case ("integer")
                        allocate (this%tail%next, source=var1d_i_type)
                    end select
                case (2)
                    allocate (this%tail%next, source=var2d_d_type)
                case (3)
                    allocate (this%tail%next, source=var3d_d_type)
                case (4)
                    allocate (this%tail%next, source=var4d_d_type)
                end select
            end if
            ptr => this%tail
            this%tail => this%tail%next
            this%tail%prev => ptr
            ptr%next => this%tail
            this%num_var = this%num_var+1
        end if

        if (present(dims)) then
            call this%tail%init(name, long_name, units, dims)
        else
            call this%tail%init(name, long_name, units)
        end if

        if (present(pointer)) then
            pointer => this%tail
        end if

    end subroutine var_list_append

    pure integer function var_list_get_num_var(this) result (num_var)

        class(var_list), intent(in) :: this

        num_var = this%num_var

    end function var_list_get_num_var

    function var_list_get_var(this, name) result (res)

        class(var_list), intent(in) :: this
        character(*), intent(in) :: name
        class(var), pointer :: res

        integer i
        class(var), pointer :: ptr

        ptr => this%get_head()
        do i = 1, this%get_num_var()
            if (ptr%get_name() == name) then
                res => ptr
                return
            end if
            ptr => ptr%next
        end do
        write(*, "('[Error]: var_list_get_var: Unknown variable ', A)") name
        write(*, "('[Debug]: var_list_get_var: Available variables:')")
        ptr => this%get_head()
        do i = 1, this%get_num_var()
            print *, ptr%get_name()
            ptr => ptr%next
        end do
        stop

    end function var_list_get_var

    function var_list_get_head(this) result (head)

        class(var_list), intent(in) :: this
        class(var), pointer :: head

        head => this%head

    end function var_list_get_head

    function var_list_get_tail(this) result (tail)

        class(var_list), intent(in) :: this
        class(var), pointer :: tail

        tail => this%tail

    end function var_list_get_tail

    subroutine var_list_get_tail_values_0d_d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values

        select type (tail => this%tail)
        type is (var0d_d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_0d_d

    subroutine var_list_get_tail_values_0d_i(this, values)

        class(var_list), intent(in) :: this
        integer, pointer :: values

        select type (tail => this%tail)
        type is (var0d_i)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_0d_i

    subroutine var_list_get_tail_values_1d_d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:)

        select type (tail => this%tail)
        type is (var1d_d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_1d_d

    subroutine var_list_get_tail_values_1d_i(this, values)

        class(var_list), intent(in) :: this
        integer, pointer :: values(:)

        select type (tail => this%tail)
        type is (var1d_i)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_1d_i

    subroutine var_list_get_tail_values_2d_d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:)

        select type (tail => this%tail)
        type is (var2d_d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_2d_d

    subroutine var_list_get_tail_values_3d_d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:,:)

        select type (tail => this%tail)
        type is (var3d_d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_3d_d

    subroutine var_list_get_tail_values_4d_d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:,:,:)

        select type (tail => this%tail)
        type is (var4d_d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_4d_d

end module variable
