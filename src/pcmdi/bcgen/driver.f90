subroutine pcmdi_bcgen(infil, outfilclim, outfilamip, mon1, iyr1, monn, iyrn, &
                       mon1rd, iyr1rd, monnrd, iyrnrd, mon1clm, iyr1clm, &
                       monnclm, iyrnclm, mon1out, iyr1out, monnout, iyrnout)
!------------------------------------------------------------------------------------
!
! Purpose: Wrapper program for routines (provided by Karl Taylor of PCMDI) which modify 
!          mid-month values of SST and ice concentration to preserve monthly means upon 
!          linear time interpolation.
!
! Method: Read in a namelist containing the relevant variables.  Check for validity,
!         then call driving subroutine.
!
!------------------------------------------------------------------------------------
   implicit none

   include 'netcdf.inc'

   character(*), intent(in) :: infil      ! input filename.
   character(*), intent(in) :: outfilclim ! output climatology file
   character(*), intent(in) :: outfilamip ! output AMIP-style file
   integer, intent(in) :: mon1            ! start month of period of interest plus buffer
   integer, intent(in) :: iyr1            ! start year of period of interest plus buffer
   integer, intent(in) :: monn            ! end month of period of interest plus buffer
   integer, intent(in) :: iyrn            ! end year of period of interest plus buffer
   integer, intent(in) :: mon1rd          ! start month to read from input data
   integer, intent(in) :: iyr1rd          ! start year to read from input data
   integer, intent(in) :: monnrd          ! end month to read from input data
   integer, intent(in) :: iyrnrd          ! end year to read from input data
   integer, intent(in) :: mon1clm         ! start month to use for climatology
   integer, intent(in) :: iyr1clm         ! start year to use for climatology
   integer, intent(in) :: monnclm         ! end month to use for climatology
   integer, intent(in) :: iyrnclm         ! end year to use for climatology
   integer, intent(in) :: mon1out         ! start month written to output file (default mon1rd)
   integer, intent(in) :: iyr1out         ! start year written to output file (default iyr1rd)
   integer, intent(in) :: monnout         ! end month written to output file (default monnrd)
   integer, intent(in) :: iyrnout         ! end year written to output file (default iyrnrd)

   integer :: nmon                ! total number of months of period of interest plus buffer
   integer :: nlat = -1           ! number of latitudes (e.g. 64 for typical T42)
   integer :: nlon = -1           ! number of longitudes (e.g. 128 for typical T42)
   logical :: oldttcalc = .false. ! true => bfb agreement with original code

   character(len=256) :: string           ! temporary character variable

   character(len=19)  :: cur_timestamp
   character(len=1024) :: prev_history = ' ' ! history attribute from input file
   character(len=1024) :: history = ' '      ! history attribute for output files
!
! netcdf info
!
   integer :: ncidin = -1         ! input file handle
   integer :: londimid = -1       ! longitude dimension id
   integer :: latdimid = -1       ! latitude dimension id
!
! Check that all required input items were specified in the namelist
!
   call verify_input (mon1, 'mon1', 'start month of period of interest plus buffer')
   call verify_input (iyr1, 'iyr1', 'start year of period of interest plus buffer')
   call verify_input (monn, 'monn', 'end month of period of interest plus buffer')
   call verify_input (iyrn, 'iyrn', 'end year of period of interest plus buffer')
   call verify_input (mon1rd, 'mon1rd', 'start month of input data')
   call verify_input (iyr1rd, 'iyr1rd', 'start year of input data')
   call verify_input (monnrd, 'monnrd', 'end month of input data')
   call verify_input (iyrnrd, 'iyrnrd', 'end year of input data')
   call verify_input (mon1clm, 'mon1clm', 'start month of output climatology')
   call verify_input (iyr1clm, 'iyr1clm', 'start year of output climatology')
   call verify_input (monnclm, 'monnclm', 'end month of output climatology')
   call verify_input (iyrnclm, 'iyrnclm', 'end year of output climatology')
!
! Check that all specified input values are valid
!
   call verify_monthindx (mon1, 'mon1')
   call verify_monthindx (monn, 'monn')
   call verify_monthindx (mon1rd, 'mon1rd')
   call verify_monthindx (monnrd, 'monnrd')
   call verify_monthindx (mon1clm, 'mon1clm')
   call verify_monthindx (monnclm, 'monnclm')
!
! Check that input dates are valid with respect to each other
!
   if ((mon1clm + 12*iyr1clm) < (mon1rd + 12*iyr1rd)) then
      write(6,*)'mon1clm, iyr1clm=', mon1clm, iyr1clm
      write(6,*)'mon1rd,  iyr1rd= ', mon1rd, iyr1rd
      call err_exit ('Dates for data to be read must bracket climatology')
   end if

   if (mon1rd + 12*iyr1rd > monnrd + 12*(iyrnrd-1)) then
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('must read at least 1 year of data total')
   end if

   if ((monnclm + 12*iyrnclm) > (monnrd + 12*iyrnrd)) then
      write(6,*)'monnclm, iyrnclm=', monnclm, iyrnclm
      write(6,*)'monnrd,  iyrnrd= ', monnrd, iyrnrd
      call err_exit ('End date for climatology exceeds end date of data to be read')
   end if

   if (iyr1rd > iyrnrd .or. iyr1rd == iyrnrd .and. mon1rd > mon1rd) then
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('Start date for input data exceeds end date')
   end if

   if (iyr1 > iyrn .or. iyr1 == iyrn .and. mon1 > monn) then
      write(6,*)'mon1, iyr1=', mon1, iyr1
      write(6,*)'monn, iyrn=', monn, iyrn
      call err_exit ('Start date for output data exceeds end date')
   end if

   if (.not. (monn == mon1-1 .or. mon1 == 1 .and. monn == 12)) then
      write(6,*)'mon1, monn=', mon1, monn
      call err_exit ('Period to be treated must be an integral number of years')
   end if

   if ((mon1 + 12*iyr1) > (mon1rd + 12*iyr1rd)) then
      write(6,*)'mon1,   iyr1=  ', mon1, iyr1
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      call err_exit ('Start date of data to be read must not preceed period of interest')
   end if

   if ((monn + 12*iyrn) < (monnrd + 12*iyrnrd)) then
      write(6,*)'monn,   iyrn=  ', monn, iyrn
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('End date of data to be read must not follow period of interest')
   end if

   if (mod((monn+12-mon1+1),12) /= 0) then
      write(6,*)'mon1, monn=',mon1, monn
      string = 'error in time specifications: only integral number of yrs allowed'
      call err_exit (string)
   end if
!
! Ensure that climatology period is as expected (1982-2001).  To enable other
! averaging periods, just comment out the following bit of code
!
!   if (mon1clm /= 1  .or. iyr1clm /= 1982 .or. &
!       monnclm /= 12 .or. iyrnclm /= 2001) then
!      write(6,*)'Climatological averaging period is not as expected (1982-2001)'
!      write(6,*)'If you REALLY want to change the averaging period, delete the'
!      write(6,*)'appropriate err_exit call from driver.f90'
!      call err_exit ('pcmdisst')
!   end if
!
! Calculate derived variables
!
   nmon = (iyrn - iyr1)*12 + (monn - mon1 + 1)

   if (monn - mon1 + 12*(iyrn - iyr1)+1 /= nmon) then
      string = 'error in time specifications: parameter nmon must be consistent with '// &
               'first and last months specified. Check nmon, mon1, iyr1, monn, iyrn'
      call err_exit (string)
   end if
!
! Check validity of output dates
!
   call verify_monthindx (mon1out, 'mon1out')
   call verify_monthindx (monnout, 'monnout')

   if (iyr1out > iyrnout .or. iyr1out == iyrnout .and. mon1out > mon1out) then
      call err_exit ('Start date for output exceeds end date')
   end if
!
! Calculate derived variables
!
   nmon = (iyrn - iyr1)*12 + (monn - mon1 + 1)
!
! Open input file and obtain grid size (nlon and nlat)
!
   call wrap_nf_open (infil, NF_NOWRITE, ncidin)
   call wrap_nf_inq_dimid (ncidin, 'lon', londimid)
   call wrap_nf_inq_dimid (ncidin, 'lat', latdimid)
   call wrap_nf_inq_dimlen (ncidin, londimid, nlon)
   call wrap_nf_inq_dimlen (ncidin, latdimid, nlat)
!
! Add to or define history attribute.
!
   if (nf_get_att_text (ncidin, NF_GLOBAL, 'history', prev_history) /= NF_NOERR) then
      write(6,*)'nf_get_att_text() failed for history attribute'
   end if
   
   call get_curr_timestamp(cur_timestamp)
   if (len_trim(prev_history) == 0) then
      history = cur_timestamp
   else
      history = trim(prev_history) // char(10) // cur_timestamp
   end if

   write(6,*)'Grid size is nlon,nlat=', nlon, nlat
!
! Call Karl's main program or the derivative code.  Either should give the same answers
!
   call bcgen (mon1, iyr1, monn, iyrn, mon1rd,                &
               iyr1rd, monnrd, iyrnrd, mon1clm, iyr1clm,      &
               monnclm, iyrnclm, nlat, nlon, mon1out,         &
               iyr1out, monnout, iyrnout, ncidin, outfilclim, &
               outfilamip, nmon, oldttcalc, history)
   write(6,*)'Done writing output files ', trim(outfilclim), ' and ', trim(outfilamip)
end subroutine pcmdi_bcgen

subroutine verify_monthindx (indx, varname)
!------------------------------------------------------------------------------------
!
! Purpose: Check validity of input month index.
!
!------------------------------------------------------------------------------------
   implicit none
!
! Input arguments
!
   integer, intent(in) :: indx
   character(len=*), intent(in) :: varname

   if (indx < 1 .or. indx > 12) then
      write(6,*) varname, '=', indx, ' is an invalid month index'
      stop 999
   end if

   return
end subroutine verify_monthindx

subroutine verify_input (ivar, ivarname, string)
!------------------------------------------------------------------------------------
!
! Purpose: Check that a required input variable was actually set.  If not, print an
!          error msg and exit.
!
! Method: All namelist input integer variables must have a non-negative value.  Since
!         they are initialized to -1, a check on < 0 effectively checks whether they 
!         were correctly set in the namelist.
!------------------------------------------------------------------------------------
   implicit none
!
! Input arguments
!
   integer, intent(in) :: ivar
   character(len=*), intent(in) :: ivarname
   character(len=*), intent(in) :: string

   if (ivar < 0) then
      write(6,*) ivarname, ' must be set in the namelist to a non-negative value'
      write(6,*) 'verbose description of this variable: ', trim(string)
      stop 999
   end if

   return
end subroutine verify_input
