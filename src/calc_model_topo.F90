subroutine calc_model_topo(NLONH, NLATH, ALONH, ALATH, OSH, &
                           NLONM, NLATM, ALONM, ALATM, INDLS, &
                           OSM, STDOSM)

    implicit none

    integer, intent(in) :: NLONH, NLATH
    real(8), intent(in) :: ALONH(NLONH), ALATH(NLATH)
    real(8), intent(in) :: OSH(NLONH,NLATH)
    integer, intent(in) :: NLONM, NLATM
    real(8), intent(in) :: ALONM(NLONM), ALATM(NLATM)
    real(8), intent(inout) :: INDLS(NLONM,NLATM)
    real(8), intent(inout) :: OSM(NLONM,NLATM)
    real(8), intent(inout) :: STDOSM(NLONM,NLATM)

    integer ixa(NLONM),ixb(NLONM),jya(NLATM),jyb(NLATM)

    real*8 xm,xa,xb,xh,ym,ya,yb,yh,dy
    real*8 osms,oshs,stdgs,stdms,anum,anuma,inds

    integer i, j, k, l,ia,ib,ja,jb

!   To obtain the indices of east-west boundary for boxes of model grids

    do i=1,NLONM
       xm=ALONM(i)

       if (i.eq.1) then
          xa=(xm+ALONM(NLONM)-360.0d0)*0.5d0
       else
          xa=(xm+ALONM(i-1))*0.5d0
       endif

       if (xa.lt.ALONH(1)) then
          xa=xa+360.0d0
       endif

       if (i.eq.NLONM) then
          xb=(xm+ALONM(1)+360.0d0)*0.5d0
       else
          xb=(xm+ALONM(i+1))*0.5d0
       endif

       if (xb.gt.ALONH(NLONH)) then
          xb=xb-360.0d0
       endif

       ia=0
       ib=0
       do k=1,NLONH
	  xh=ALONH(k)
	  if ((xh.ge.xa).and.(ia.eq.0)) then
	     ia=k
	  endif

	  if ((xh.ge.xb).and.(ib.eq.0)) then
	     ib=k-1
         endif
       enddo 

       ixa(i)=ia
       ixb(i)=ib
    enddo

!   do i=1,NLONM
!      print *,i,ixa(i),ixb(i)
!      print *,alonh(ixa(i)),alonm(i),alonh(ixb(i))
!   enddo
!   stop

!   To obtain the indices of south-north boundary for boxes of model grids

    do j=2,NLATM-1
       ym=ALATM(j)

       ya=(ym+ALATM(j-1))*0.5d0
       yb=(ym+ALATM(j+1))*0.5d0

       dy=(yb-ya)*0.5d0

       ya=ym-dy
       yb=ym+dy

       ja=0
       jb=0
       do k=1,NLATH
	  yh=ALATH(k)
	  if ((yh.ge.ya).and.(ja.eq.0)) then
	     ja=k
	  endif

	  if ((yh.ge.yb).and.(jb.eq.0)) then
	     jb=k-1
	  endif
       enddo 

       jya(j)=ja
       jyb(j)=jb
    enddo

    jya(1)=1
    jyb(1)=jya(2)-1

    jya(NLATM)=jyb(NLATM-1)+1
    jyb(NLATM)=NLATH
 
!   do j=1,NLATM
!      print *,j,jya(j),jyb(j)
!      print *,alath(jya(j)),alatm(j),alath(jyb(j))
!   enddo
!   stop

!   To calculate the averaged terrain of the model grid box as the model terrain,
!   the corresponding standard deviation, and the land-sea mask

!   For the region excluding the poles

    do j=2,NLATM-1
       ja=jya(j)
       jb=jyb(j)

       if (ixa(1).gt.ixb(1)) then
          ia=1
          ib=ixb(1)
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(1,j)=osms
          INDLS(1,j)=inds
          STDOSM(1,j)=stdms
          anuma=anum

          ia=ixa(1)
          ib=NLONH
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(1,j)=osms+OSM(1,j)
          INDLS(1,j)=inds+INDLS(1,j)
          STDOSM(1,j)=stdms+STDOSM(1,j)
          anuma=anum+anuma

          OSM(1,j)=OSM(1,j)/anuma
          INDLS(1,j)=INDLS(1,j)/anuma
          STDOSM(1,j)=dsqrt(STDOSM(1,j)/anuma)
       else
	  ia=ixa(1)
          ib=ixb(1)
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(1,j)=osms/anum
          INDLS(1,j)=inds/anum
          STDOSM(1,j)=dsqrt(stdms/anum)
       endif

       do i=2,NLONM-1
	  ia=ixa(i)
          ib=ixb(i)
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(i,j)=osms/anum
          INDLS(i,j)=inds/anum
          STDOSM(i,j)=dsqrt(stdms/anum)
       enddo

       if (ixa(NLONM).gt.ixb(NLONM)) then
          ia=1
          ib=ixb(NLONM)
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(NLONM,j)=osms
          INDLS(NLONM,j)=inds
          STDOSM(NLONM,j)=stdms
          anuma=anum

          ia=ixa(1)
          ib=NLONH
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(NLONM,j)=osms+OSM(NLONM,j)
          INDLS(NLONM,j)=inds+INDLS(NLONM,j)
          STDOSM(NLONM,j)=stdms+STDOSM(NLONM,j)
          anuma=anum+anuma

          OSM(NLONM,j)=OSM(NLONM,j)/anuma
          INDLS(NLONM,j)=INDLS(NLONM,j)/anuma
          STDOSM(NLONM,j)=dsqrt(STDOSM(NLONM,j)/anuma)
       else
	  ia=ixa(NLONM)
          ib=ixb(NLONM)
          call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
          OSM(NLONM,j)=osms/anum
          INDLS(NLONM,j)=inds/anum
          STDOSM(NLONM,j)=dsqrt(stdms/anum)
       endif
    enddo 

!   For the poles

    do j=1,NLATM,NLATM-1
       ja=jya(j)
       jb=jyb(j)
    
       ia=1
       ib=NLONH
       call termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
       osms=osms/anum
       inds=inds/anum
       stdms=stdms/anum

       do i=1,NLONM
          OSM(i,j)=osms
          INDLS(i,j)=inds
          STDOSM(i,j)=dsqrt(stdms)
       enddo
    enddo 

!   do j=1,NLATM
!      do i=1,NLONM
!         print *,j,i,real(OSM(i,j)),real(INDLS(i,j)),real(STDOSM(i,j))
!      enddo
!   enddo
!   stop

    return

end subroutine calc_model_topo

subroutine termask(osms,stdms,inds,anum,ia,ib,ja,jb,OSH,NLONH,NLATH)
    integer, intent(in) :: NLONH, NLATH
    real(8), intent(in) :: OSH(NLONH,NLATH)

    real*8 osms,osmm,oshs,stdgs,stdms,anum,inds

    integer k, l,ia,ib,ja,jb

    osms=0.0d0
    inds=0.0d0
    anum=0.0d0
    do k=ja,jb
       do l=ia,ib
          if (OSH(l,k).gt.0.0d0) then
             osms=osms+OSH(l,k)
             inds=inds+1.0d0
          endif
          anum=anum+1.0d0
       enddo
    enddo
    osmm=osms/anum

    stdms=0.0d0
    do k=ja,jb
       do l=ia,ib
          oshs=OSH(l,k)
          if (oshs.le.0.0d0) oshs=0.0
          stdgs=oshs-osmm
          stdms=stdms+stdgs*stdgs
       enddo
    enddo

    return

end subroutine termask
