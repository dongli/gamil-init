! ------------------------------------------------------------------------------
! Description:
!
!   This module encloses several routines that is useful for processing data
!   relating to model.
!
! Authors:
!
!   Li Dong - dongli@lasg.iap.ac.cn
! ------------------------------------------------------------------------------

module model_gears

    use interp

    implicit none

    private

    public model_gears_interp_h
    public model_gears_interp_v

contains

    subroutine model_gears_interp_h(lon1, lat1, data1, lon2, lat2, data2, flag)

        real(8), intent(in) :: lon1(:), lat1(:), data1(:,:), lon2(:), lat2(:)
        real(8), intent(out) :: data2(:,:)
        integer, intent(in) :: flag

        integer size(1), num_lon1, num_lat1, num_lon2, num_lat2

        size = shape(lon1); num_lon1 = size(1)
        size = shape(lat1); num_lat1 = size(1)
        size = shape(lon2); num_lon2 = size(1)
        size = shape(lat2); num_lat2 = size(1)

        call interp_bilinear(lon1, lat1, data1, lon2, lat2, data2, 360.0d0)

        ! hanele pole grids
        if (flag == 0) then
            data2(:,1) = 0.0
            data2(:,num_lat2) = 0.0
        else if (flag == 1) then
            data2(:,1) = sum(data1(:,1))/num_lon1
            data2(:,num_lat2) = sum(data1(:,num_lat1))/num_lon1
        end if

    end subroutine model_gears_interp_h

    subroutine model_gears_interp_v(lev1, data1, lev2, data2, flag)
 
        real(8), intent(in) :: lev1(:), data1(:), lev2(:)
        real(8), intent(out) :: data2(:)
        integer, intent(in) :: flag

        if (flag == 0) then
            call interp_linear(lev1, data1, lev2, data2)
        else
            call interp_log_linear(lev1, data1, lev2, data2)
        end if

    end subroutine model_gears_interp_v

end module model_gears
