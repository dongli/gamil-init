subroutine calc_model_topo(num_lon_topo, num_lat_topo, lon_topo, lat_topo, topo, &
                           num_lon_model, num_lat_model, lon_model, lat_model, &
                           landfrac, phis, sgh)

    implicit none

    integer, intent(in) :: num_lon_topo, num_lat_topo
    real(8), intent(in) :: lon_topo(num_lon_topo), lat_topo(num_lat_topo)
    real(8), intent(in) :: topo(num_lon_topo,num_lat_topo)
    integer, intent(in) :: num_lon_model, num_lat_model
    real(8), intent(in) :: lon_model(num_lon_model), lat_model(num_lat_model)
    real(8), intent(inout) :: landfrac(num_lon_model,num_lat_model)
    real(8), intent(inout) :: phis(num_lon_model,num_lat_model)
    real(8), intent(inout) :: sgh(num_lon_model,num_lat_model)

    integer ix1(num_lon_model)
    integer ix2(num_lon_model)
    integer jy1(num_lat_model)
    integer jy2(num_lat_model)

    real(8) lon1, lon2
    real(8) lat1, lat2, dlat
    real(8) phis0, landfrac0, sgh0, num_topo_grid, num_topo_grid0

    integer i, j, k, l, i1, i2, j1, j2

    ! To obtain the indices of east-west boundary for boxes of model grids

    do i = 1, num_lon_model
        if (i == 1) then
            lon1 = (lon_model(i) + lon_model(num_lon_model) - 360.0d0) * 0.5d0
        else
            lon1 = (lon_model(i) + lon_model(i-1)) * 0.5d0
        end if

        if (lon1 < lon_topo(1)) then
            lon1 = lon1 + 360.0d0
        end if

        if (i == num_lon_model) then
            lon2 = (lon_model(i) + lon_model(1) + 360.0d0) * 0.5d0
        else
            lon2 = (lon_model(i) + lon_model(i+1)) * 0.5d0
        end if

        if (lon2 > lon_topo(num_lon_topo)) then
            lon2 = lon2 - 360.0d0
        end if

        i1 = 0
        i2 = 0
        do k = 1, num_lon_topo
            if (lon_topo(k) >= lon1 .and. i1 == 0) then
                i1 = k
            end if
            if (lon_topo(k) >= lon2 .and. i2 == 0) then
                i2 = k - 1
            end if
            if (i1 /= 0 .and. i2 /= 0) then
                exit
            end if
        end do

        ix1(i) = i1
        ix2(i) = i2
    end do

    ! To obtain the indices of south-north boundary for boxes of model grids

    do j = 2, num_lat_model - 1
        lat1 = (lat_model(j) + lat_model(j-1))*0.5d0
        lat2 = (lat_model(j) + lat_model(j+1))*0.5d0

        dlat = (lat2 - lat1) * 0.5d0

        lat1 = lat_model(j) - dlat
        lat2 = lat_model(j) + dlat

        j1 = 0
        j2 = 0
        do k = 1, num_lat_topo
            if (lat_topo(k) >= lat1 .and. j1 == 0) then
                j1 = k
            end if

            if (lat_topo(k) >= lat2 .and. j2 == 0) then
                j2 = k - 1
            end if

            if (j1 /= 0 .and. j2 /= 0) then
                exit
            end if
        end do

        jy1(j) = j1
        jy2(j) = j2
    end do

    jy1(1) = 1
    jy2(1) = jy1(2) - 1

    jy1(num_lat_model) = jy2(num_lat_model-1) + 1
    jy2(num_lat_model) = num_lat_topo
 
    ! To calculate the averaged terrain of the model grid box as the model terrain,
    ! the corresponding standard deviation, and the land-sea mask

    ! For the region excluding the poles
    do j = 2, num_lat_model - 1
        j1 = jy1(j)
        j2 = jy2(j)

       if (ix1(1) > ix2(1)) then
            ! Split into two parts
            i1 = 1
            i2 = ix2(1)
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(1,j)      = phis0
            landfrac(1,j)  = landfrac0
            sgh(1,j)       = sgh0
            num_topo_grid0 = num_topo_grid

            i1 = ix1(1)
            i2 = num_lon_topo
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(1,j)      = phis(1,j)     + phis0
            landfrac(1,j)  = landfrac(1,j) + landfrac0
            sgh(1,j)       = sgh(1,j)      + sgh0
            num_topo_grid0 = num_topo_grid + num_topo_grid0

            phis(1,j)     = phis(1,j)      / num_topo_grid0
            landfrac(1,j) = landfrac(1,j)  / num_topo_grid0
            sgh(1,j)      = dsqrt(sgh(1,j) / num_topo_grid0)
        else
            i1 = ix1(1)
            i2 = ix2(1)
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(1,j)     = phis0      / num_topo_grid
            landfrac(1,j) = landfrac0  / num_topo_grid
            sgh(1,j)      = dsqrt(sgh0 / num_topo_grid)
        end if

        do i = 2, num_lon_model - 1
            i1 = ix1(i)
            i2 = ix2(i)
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(i,j)     = phis0      / num_topo_grid
            landfrac(i,j) = landfrac0  / num_topo_grid
            sgh(i,j)      = dsqrt(sgh0 / num_topo_grid)
        end do

        if (ix1(num_lon_model) > ix2(num_lon_model)) then
            i1 = 1
            i2 = ix2(num_lon_model)
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(num_lon_model,j)     = phis0
            landfrac(num_lon_model,j) = landfrac0
            sgh(num_lon_model,j)      = sgh0
            num_topo_grid0            = num_topo_grid

            i1 = ix1(1)
            i2 = num_lon_topo
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(num_lon_model,j)     = phis(num_lon_model,j)     + phis0
            landfrac(num_lon_model,j) = landfrac(num_lon_model,j) + landfrac0
            sgh(num_lon_model,j)      = sgh(num_lon_model,j)      + sgh0
            num_topo_grid0            = num_topo_grid + num_topo_grid0

            phis(num_lon_model,j)     = phis(num_lon_model,j)      / num_topo_grid0
            landfrac(num_lon_model,j) = landfrac(num_lon_model,j)  / num_topo_grid0
            sgh(num_lon_model,j)      = dsqrt(sgh(num_lon_model,j) / num_topo_grid0)
        else
            i1 = ix1(num_lon_model)
            i2 = ix2(num_lon_model)
            call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
            phis(num_lon_model,j)     = phis0      / num_topo_grid
            landfrac(num_lon_model,j) = landfrac0  / num_topo_grid
            sgh(num_lon_model,j)      = dsqrt(sgh0 / num_topo_grid)
        end if
    end do

    ! For the poles

    do j=1, num_lat_model, num_lat_model - 1
        j1 = jy1(j)
        j2 = jy2(j)
    
        i1 = 1
        i2 = num_lon_topo
        call termask(phis0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
        phis0     = phis0     / num_topo_grid
        landfrac0 = landfrac0 / num_topo_grid
        sgh0      = sgh0      / num_topo_grid

        do i = 1, num_lon_model
            phis(i,j)     = phis0
            landfrac(i,j) = landfrac0
            sgh(i,j)      = dsqrt(sgh0)
        end do
    end do

    return

end subroutine calc_model_topo

subroutine termask(phis, sgh, landfrac, num_topo_grid, i1, i2, j1, j2, topo, num_lon_topo, num_lat_topo)
    integer, intent(in) :: i1, i2, j1, j2, num_lon_topo, num_lat_topo
    real(8), intent(in) :: topo(num_lon_topo,num_lat_topo)

    real(8), intent(out) :: phis, sgh, landfrac, num_topo_grid
    real(8) osmm, oshs, stdgs

    integer k, l

    phis          = 0.0d0
    landfrac      = 0.0d0
    num_topo_grid = 0.0d0
    do k = j1, j2
        do l = i1, i2
            if (topo(l,k) > 0.0d0) then
                phis     = phis     + topo(l,k)
                landfrac = landfrac + 1.0d0
            end if
            num_topo_grid = num_topo_grid + 1.0d0
        end do
    end do
    osmm = phis / num_topo_grid

    sgh = 0.0d0
    do k = j1, j2
        do l = i1, i2
            oshs = topo(l,k)
            if (oshs <= 0.0d0) oshs = 0.0
            stdgs = oshs - osmm
            sgh = sgh + stdgs * stdgs
        end do
    end do

    return

end subroutine termask
