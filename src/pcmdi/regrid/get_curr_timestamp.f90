subroutine get_curr_timestamp(time)
! return timestamp formatted as "YYYY-MM-DD HH:MM:SS"   
   character(len=19), intent(out) :: time
   integer :: t(8)
   call date_and_time(values=t)
   write(time,'(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') t(1),'-',t(2),'-',t(3),' ',&
                                                         t(5),':',t(6),':',t(7)
end subroutine get_curr_timestamp
