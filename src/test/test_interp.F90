program test_interp

    use interp

    real(8), parameter :: PI = 4.0d0*atan(1.0d0)
    real(8), parameter :: Rad2Deg = 180.0d0/PI

    real(8) dx
    real(8) x1(10), y1(10)
    real(8) x2(10), y2(10)

    integer i

    dx = 1.0d0/10
    do i = 1, 10
        x1(i) = 0.5d0*dx+(i-1)*dx
        y1(i) = sin(2.0d0*PI*x1(i))
    end do

    do i = 1, 10
        x2(i) = (i-1)*dx
    end do

    call interp_linear(x1, y1, x2, y2)

    open (10, file="test_interp.dat")
    do i = 1, 10
        write(10, *) x1(i), y1(i), x2(i), y2(i)
    end do
    close(10)

end program test_interp
