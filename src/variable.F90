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

    public var, var0d, var1d, var2d, var3d, var4d, var_list

    ! --------------------------------------------------------------------------
    integer, parameter :: name_str_len = 8
    integer, parameter :: long_name_str_len = 255
    integer, parameter :: units_str_len = 30

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
        procedure :: reshape => var_reshape
        ! ======================================================================
        procedure :: var_set_values_1d
        generic :: set_values => var_set_values_1d
    end type var

    type, extends(var) :: var0d
        real(8), pointer, private :: values
    contains
        procedure :: get_values => var0d_get_values
    end type var0d

    type, extends(var) :: var1d
        integer, private :: dims(1)
        real(8), pointer, private :: values(:)
    contains
        procedure :: get_values => var1d_get_values
    end type var1d

    type, extends(var) :: var2d
        integer, private :: dims(2)
        real(8), pointer, private :: values(:,:)
    contains
        procedure :: get_values => var2d_get_values
    end type var2d

    type, extends(var) :: var3d
        integer, private :: dims(3)
        real(8), pointer, private :: values(:,:,:)
    contains
        procedure :: get_values => var3d_get_values
    end type var3d

    type, extends(var) :: var4d
        integer, private :: dims(4)
        real(8), pointer, private :: values(:,:,:,:)
    contains
        procedure :: get_values => var4d_get_values
    end type var4d

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
        procedure :: var_list_get_tail_values_0d
        procedure :: var_list_get_tail_values_1d
        procedure :: var_list_get_tail_values_2d
        procedure :: var_list_get_tail_values_3d
        procedure :: var_list_get_tail_values_4d
        generic :: get_tail_values => var_list_get_tail_values_0d, &
                                      var_list_get_tail_values_1d, &
                                      var_list_get_tail_values_2d, &
                                      var_list_get_tail_values_3d, &
                                      var_list_get_tail_values_4d
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
        type is (var0d)
            allocate (this%values)
        type is (var1d)
            if (size(dims) /= 1) then
                write(*, "('[Error]: var_init: Dimension number /= 1!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1)))
        type is (var2d)
            if (size(dims) /= 2) then
                write(*, "('[Error]: var_init: Dimension number /= 2!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1),dims(2)))
        type is (var3d)
            if (size(dims) /= 3) then
                write(*, "('[Error]: var_init: Dimension number /= 3!')")
                stop
            end if
            this%dims = dims
            allocate (this%values(dims(1),dims(2),dims(3)))
        type is (var4d)
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
        type is (var0d)
            res = 1
        type is (var1d)
            res = this%dims(1)
        type is (var2d)
            res = this%dims(dim_idx)
        type is (var3d)
            res = this%dims(dim_idx)
        type is (var4d)
            res = this%dims(dim_idx)
        end select

    end function var_get_dims

    subroutine var_reshape(this, order)

        class(var), intent(inout) :: this
        integer, intent(in) :: order(:)

        integer dims2d(2), dims3d(3), dims4d(4)
        real(8), allocatable :: tmp2d(:,:), tmp3d(:,:,:), tmp4d(:,:,:,:)

        select type (this)
        type is (var3d)
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

    function var0d_get_values(this) result (res)

        class(var0d), intent(in) :: this
        real(8), pointer :: res

        select type (this)
        type is (var0d)
            res => this%values
        end select

    end function var0d_get_values

    function var1d_get_values(this) result (res)

        class(var1d), intent(in) :: this
        real(8), pointer :: res(:)

        res => this%values

    end function var1d_get_values

    subroutine var_set_values_1d(this, values)

        class(var), intent(inout) :: this
        real(8), intent(in) :: values(:)

        select type (this)
        type is (var1d)
            this%values = values
        class default
            write(*, "('[Error]: var_set_values_1d: Variable is not 1d!')")
            stop
        end select

    end subroutine var_set_values_1d

    function var2d_get_values(this) result (res)

        class(var2d), intent(in) :: this
        real(8), pointer :: res(:,:)

        res => this%values

    end function var2d_get_values

    function var3d_get_values(this) result (res)

        class(var3d), intent(in) :: this
        real(8), pointer :: res(:,:,:)

        res => this%values

    end function var3d_get_values

    function var4d_get_values(this) result (res)

        class(var4d), intent(in) :: this
        real(8), pointer :: res(:,:,:,:)

        res => this%values

    end function var4d_get_values

    subroutine var_list_append(this, name, long_name, units, dims, pointer)

        class(var_list), intent(inout) :: this
        character(*), intent(in) :: name, long_name, units
        integer, intent(in), optional :: dims(:)
        class(var), intent(out), optional, pointer :: pointer

        class(var), pointer :: ptr

        type(var0d) var0d_type
        type(var1d) var1d_type
        type(var2d) var2d_type
        type(var3d) var3d_type
        type(var4d) var4d_type

        if (.not. associated(this%head)) then
            if (.not. present(dims)) then
                allocate (this%head, source=var0d_type)
            else
                select case (size(dims))
                case (1)
                    allocate (this%head, source=var1d_type)
                case (2)
                    allocate (this%head, source=var2d_type)
                case (3)
                    allocate (this%head, source=var3d_type)
                case (4)
                    allocate (this%head, source=var4d_type)
                end select
            end if
            this%tail => this%head
            nullify(this%head%prev)
            nullify(this%head%next)
            this%num_var = 1
        else
            if (.not. present(dims)) then
                allocate (this%tail%next, source=var0d_type)
            else
                select case (size(dims))
                case (1)
                    allocate (this%tail%next, source=var1d_type)
                case (2)
                    allocate (this%tail%next, source=var2d_type)
                case (3)
                    allocate (this%tail%next, source=var3d_type)
                case (4)
                    allocate (this%tail%next, source=var4d_type)
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

    subroutine var_list_get_tail_values_0d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values

        select type (tail => this%tail)
        type is (var0d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_0d

    subroutine var_list_get_tail_values_1d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:)

        select type (tail => this%tail)
        type is (var1d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_1d

    subroutine var_list_get_tail_values_2d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:)

        select type (tail => this%tail)
        type is (var2d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_2d

    subroutine var_list_get_tail_values_3d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:,:)

        select type (tail => this%tail)
        type is (var3d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_3d

    subroutine var_list_get_tail_values_4d(this, values)

        class(var_list), intent(in) :: this
        real(8), pointer :: values(:,:,:,:)

        select type (tail => this%tail)
        type is (var4d)
            values => tail%values
        end select

    end subroutine var_list_get_tail_values_4d

end module variable
