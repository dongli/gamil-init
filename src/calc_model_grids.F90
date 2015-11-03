subroutine calc_model_grids(NX, NY)

    use utils
    use constants

    implicit none

    integer, intent(in) :: NX, NY

    character(50), parameter :: sub_name = "calc_model_grids"

    real(8) DX       ! ZONAL STEPSIZE
    real(8) DY       ! MERIDIONAL STEPSIZE
    real(8) YTHU(NY) ! LATITUDE AT THE NORNAL MERIDIONAL GRID Yj
    real(8) YTHV(NY) ! LATITUDE AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
    real(8) WTGU(NY) ! AREA-WEIGHTING AT THE NORNAL MERIDIONAL GRID Yj
    real(8) WTGV(NY) ! AREA-WEIGHTING AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
!
!     4) WORKING VARIABLES: U,S,DA,A2,U2,A,B,J,M1
!
    real(8) AS ! AREA SIZE WITH RESPECT TO LATITUDE
    real(8) B  ! A CONTROL PARAMETER RELATED TO THE MERIDIONAL RESOLUTION
!
! The resolution formula:
!
!   DY = (180-10B)/(NY-1) (in degree) PI(1-B/18)/(NY-1) (in arc)
!
!   when B=2.0, NY=41, DY=4 degree; when B=2.0, NY=81, DY=2 degree
!
    real(8) A  ! A  = B/(0.5*PI) A CONTROL PARAMETER RELATED TO THE AREA-SIZE COMPUTING
    real(8) A2 ! A2 = A*2
    real(8) DA ! DA = 1/A
    real(8) S  ! S  = AS(PI/2)   THE TOTAL AREA-SIZE OF THE WHOLE REGION
    real(8) S1 ! S1 = AS(-PI/3)	 THE AREA-SIZE AT THE POINT THETA=-PI/3
    real(8) S2 ! S2 = AS(PI/3)	 THE AREA-SIZE AT THE POINT THETA=PI/3
    real(8) U1 ! U1 = 1+2B/3     NEEDED WHEN CALCULATE YTHU, YTHV
    real(8) U2 ! U2 = (1-B/3)^2	 NEEDED WHEN CALCULATE YTHU, YTHV
    integer I, J
    integer M1 ! M1 = NY-1
!
! Description of the even-area method developed by Bin Wang in 2000:
!
!   Under the weighting w(theta)=1.0-a(|theta|-PI/3) when |theta|>=PI/3
!                       w(theta)=1.0                 when |theta|< PI/3
!
!   The area size AS is:
!     AS(theta)=[(1+a*theta+a*PI/3)^2-(1-b/3)^2]/(2a)  when -PI/2<=theta<=-PI/3
!     AS(theta)=AS(-PI/3)+(theta+PI/3)                 when -PI/3< theta<= PI/3
!     AS(theta)=AS(PI/3)+[1-(1-a*theta+a*PI/3)^2]/(2a) when  PI/3< theta<= PI/2
!     where a= b/(0.5*Pi), 0<b<1, b=1.0-0.25*sqrt(25/(m1-20))
!
!   Suppose the total area-size of the whole region S is partitioned into NY-1
!   equal small area: AS(theta(j+1))-AS(theta(j))=DY=constant, then theta(j)
!   can be calculated according to the formula of AS. Obviously,
!   theta(j+1)-theta(j) will not be a constant when they are not in the interval
!   [-PI/3, PI/3]. Especially, closer to poles theta(j) is, bigger
!   theta(j+1)-theta(j) will become. In this way, the physical stepsizes in the
!   polar regions increase and the computational stability becomes better.
!   Note that: the physical mesh is not even, but the computing mesh is, which
!   makes the meridional discretization easy.
!
    ! ??????????????????????????????????????????????????????????????????????????
    ! Temporal codes:
    integer ncid, dim_id, var_id, ierr, k
    real(8), parameter :: cam_ref_ps = 1.0d5 ! Reference surface pressure (Pa)
    ! ??????????????????????????????????????????????????????????????????????????

    ! --------------------------------------------------------------------------
    ! Set grid longitudes
    num_model_lon = NX
    call model_dims%append("lon", "longitude", "degrees_east", [num_model_lon])
    call model_dims%get_tail_values(model_lon)
    call model_dims%append("lon_bnds", "longitude bounds", "degrees_east", [num_model_lon])
    call model_dims%get_tail_values(model_lon_bnds)
    DX = PI*2.0d0/NX
    do i = 1, NX
        model_lon(i) = DX*(i-1)
        model_lon_bnds(i) = model_lon(i)+DX*0.5d0
    end do
    model_lon = model_lon*Rad2Deg
    model_lon_bnds = model_lon_bnds*Rad2Deg

    ! --------------------------------------------------------------------------
    ! Set grid latitude.
    num_model_lat = NY
    call model_dims%append("lat", "latitude", "degrees_north", [num_model_lat])
    call model_dims%get_tail_values(model_lat)
    call model_dims%append("lat_bnds", "latitude", "degrees_north", [num_model_lat])
    call model_dims%get_tail_values(model_lat_bnds)

    if (model_grid_type == equal_interval_grid) then
        call notice(sub_name, "GAMIL grid is equal-interval grid")
        DY = PI/(NY-1)
        do j = 1, NY
            model_lat(j) = -PI*0.5+DY*(j-1)
            model_lat_bnds(j) = model_lat(j)+DY*0.5d0
        end do
        model_lat = model_lat*Rad2Deg
        model_lat_bnds = model_lat_bnds*Rad2Deg
    else if (model_grid_type == even_area_grid) then
        call notice(sub_name, "GAMIL grid is even-area grid")
        M1= NY-1
        B = 2
        A  = B*2.0D0/PI
        A2 = A*2.0D0
        DA = 1.0D0/A
        S  = PI*(1.0D0-B/18.0D0)
        S1 = PI*(1.0D0-B/6.0D0)/6.0D0
        S2 = PI*(5.0D0-B/6.0D0)/6.0D0
        U1 = 1.0D0+2.0D0*B/3.0D0
        U2 = (1.0D0-B/3.0D0)*(1.0D0-B/3.0D0)
        DY = S/DFLOAT(M1)
        DO J = 0, M1
            AS = DY*DFLOAT(J)
            IF (AS.LE.S1) THEN
                YTHU(J+1) = (DSQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
                IF (YTHU(J+1).LT.0.0) YTHU(J+1)=0.0D0
                WTGU(J+1) = 1.0D0-A*(DABS(YTHU(J+1)-PI*0.5D0)-PI/3.0D0)
            ELSE IF (AS.LE.S2) THEN
                YTHU(J+1) = AS-S1-PI/3.0D0+PI*0.5D0
                WTGU(J+1)=1.0D0
            ELSE
                YTHU(J+1)=(U1-DSQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
                WTGU(J+1)=1.0D0-A*(DABS(YTHU(J+1)-PI*0.5D0)-PI/3.0D0)
            END IF
            AS = DY*(DFLOAT(J)+0.5D0)
            IF (AS.LE.S1) THEN
                YTHV(J+1) = (DSQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
                WTGV(J+1) = 1.0D0-A*(DABS(YTHV(J+1)-PI*0.5D0)-PI/3.0D0)
            ELSE IF (AS.LE.S2) THEN
                YTHV(J+1) = AS-S1-PI/3.0+PI*0.5
                WTGV(J+1) = 1.0D0
            ELSE IF (J.LT.M1) THEN
                YTHV(J+1)=(U1-DSQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
                WTGV(J+1)=1.0D0-A*(DABS(YTHV(J+1)-PI*0.5D0)-PI/3.0D0)
            END IF
        END DO

        YTHV(NY) = PI
        WTGV(NY) = 1.0D0-A*(PI*0.5-PI/3.0)

        do j = 1, NY
            model_lat(j) = YTHU(J)*Rad2Deg-90.0d0
            model_lat_bnds(j) = YTHV(J)*Rad2Deg-90.0d0
        end do
    end if

    model_lat(1)  = max(-90.0d0, model_lat(1))
    model_lat(NY) = min(90.0d0, model_lat(NY))
    if (abs(model_lat(NY)-90.0d0) < 1.0e-6) then
        model_lat(NY) = 90.0d0
    end if

    ! --------------------------------------------------------------------------
    ! Set vertical grids.
    ! ??????????????????????????????????????????????????????????????????????????
    ! TODO: For the time being, we can just read the sigma coordinates from file
    !       but we should use more elegant way to set them.
    call report_warning(sub_name, "Read vertical sigma coordinates from file!")

    call check_file_exist(sub_name, "./hybpara.26.nc")

    ierr = nf90_open("./hybpara.26.nc", nf90_nowrite, ncid)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_inq_dimid(ncid, "lev", dim_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_inquire_dimension(ncid, dim_id, len=num_model_lev)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    call model_dims_aux%append("P0", "reference pressure", "Pa")
    call model_dims_aux%get_tail_values(model_p0)
    call model_dims_aux%append("hyam", "hybrid A coefficient at layer midpoints", "1", [num_model_lev])
    call model_dims_aux%get_tail_values(model_hyam)
    call model_dims_aux%append("hybm", "hybrid B coefficient at layer midpoints", "1", [num_model_lev])
    call model_dims_aux%get_tail_values(model_hybm)
    call model_dims_aux%append("hyai", "hybrid A coefficient at layer interfaces", "1", [num_model_lev+1])
    call model_dims_aux%get_tail_values(model_hyai)
    call model_dims_aux%append("hybi", "hybrid B coefficient at layer interfaces", "1", [num_model_lev+1])
    call model_dims_aux%get_tail_values(model_hybi)

    ierr = nf90_inq_varid(ncid, "P0", var_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_get_var(ncid, var_id, model_p0)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_inq_varid(ncid, "hyam", var_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_get_var(ncid, var_id, model_hyam)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_inq_varid(ncid, "hybm", var_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_get_var(ncid, var_id, model_hybm)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_inq_varid(ncid, "hyai", var_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_get_var(ncid, var_id, model_hyai)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_inq_varid(ncid, "hybi", var_id)
    call handle_netcdf_error(sub_name, __LINE__, ierr)
    ierr = nf90_get_var(ncid, var_id, model_hybi)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ierr = nf90_close(ncid)
    call handle_netcdf_error(sub_name, __LINE__, ierr)

    ! NOTICE: Although the vertical coordinate data in the file show itself as
    !         hybrid sigma-pressure coordinate, it is classic sigma-pressure
    !         coordinate in GAMIL currently!
    model_vertical_coordinate_type = classic_sigma_pressure
    call model_dims%append("lev", "", "", [num_model_lev])
    call model_dims%get_tail_values(model_lev)
    call model_dims%append("ilev", "", "", [num_model_lev+1])
    call model_dims%get_tail_values(model_lev_bnds)
    call model_dims_aux%append("pmtop", "model top pressure", "Pa")
    call model_dims_aux%get_tail_values(model_pt)
    ! Set the classic sigma-pressure coordinates from the read hybrid ones.
    model_pt = model_hyai(1)*model_p0
    model_lev_bnds(1) = 0.0d0
    do k = 1, num_model_lev
        model_lev_bnds(k+1) = (model_hyai(k+1)*cam_ref_ps+ &
                               model_hybi(k+1)*model_p0-model_pt)/ &
                              (model_p0-model_pt)
        model_lev(k) = 0.5d0*(model_lev_bnds(k)+model_lev_bnds(k+1))
    end do
    ! ??????????????????????????????????????????????????????????????????????????

    model_2d_dims = [num_model_lon,num_model_lat]
    model_3d_dims = [num_model_lon,num_model_lat,num_model_lev]

end subroutine calc_model_grids
