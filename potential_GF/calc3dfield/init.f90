subroutine coord_init
use common_val
implicit none
integer :: idb, jdb, kdb
!
epsx = 127d0*real(NNX)/real(NNX)
epsy = 127d0*real(NNY)/real(NNY)
epsz = 252d0*real(NNZ)/real(NNZ)
!
alphax = (para_x0*(real(nx)*0.5d0)**nbeki &
     &  - real(NNX)*0.5d0*para_ix0**nbeki) &
     &  /(para_ix0*(real(nx)*0.5d0) &
     & *((real(nx)*0.5d0)**(nbeki-1) &
     &  -(para_ix0**(nbeki-1))))
betax =  (real(NNX)*0.5d0*para_ix0 &
     &  - real(nx)*0.5d0*para_x0) &
     &  /(para_ix0*(real(nx)*0.5d0) &
     & *((real(nx)*0.5d0)**(nbeki-1) &
     &  -(para_ix0**(nbeki-1))))
alphay = (para_y0*(real(ny)*0.5d0)**nbeki &
     &  - real(NNY)*0.5d0*para_jy0**nbeki) &
     &  /(para_jy0*(real(ny)*0.5d0) &
     & *((real(ny)*0.5d0)**(nbeki-1) &
     &  -(para_jy0**(nbeki-1))))
betay =  (real(NNY)*0.5d0*para_jy0 &
     &  - real(ny)*0.5d0*para_y0) &
     &  /(para_jy0*(real(ny)*0.5d0) &
     & *((real(ny)*0.5d0)**(nbeki-1) &
     &  -(para_jy0**(nbeki-1))))
alphaz = (para_z0*real(nz)**nbeki &
     &  - real(NNZ)*para_kz0**nbeki) &
     &  /(para_kz0*real(nz) &
     & *(real(nz)**(nbeki-1) &
     &  -(para_kz0**(nbeki-1))))
betaz =  (real(NNZ)*para_kz0 &
     &  - real(nz)*para_z0) &
     &  /(para_kz0*real(nz) &
     & *(real(nz)**(nbeki-1) &
     &  -(para_kz0**(nbeki-1))))
!
delx = real(NNX)/real(nx)
!delx = xl/(nx-nx0)
do idb = nx0, nx
!       xc(idb) = ((idb/rnx) - 0.5d0)*NNX
   xc(idb) = alphax*(idb-nx/2) + betax*(idb-nx/2)**nbeki
!   xc(idb) = (idb-nx/2)*delx
enddo
!
!dely = yl/(ny-ny0)
dely = real(NNY)/real(ny)
do jdb = ny0, ny
!       yc(jdb) = ((jdb/rny) - 0.5d0)*NNY
   yc(jdb) = alphay*(jdb-ny/2) + betay*(jdb-ny/2)**nbeki
!   yc(jdb) = (jdb-ny/2)*dely
enddo
!
!delz = zl/(nz-nz0)
delz = real(NNZ)/real(nz)
do kdb = nz0, nz
!       zc(kdb) = (kdb/rnz)*NNZ
   zc(kdb) = alphaz*kdb + betaz*kdb**nbeki
!   zc(kdb) = kdb*delz
enddo
!
end subroutine coord_init
!!!
!!!
!!!
subroutine boundary_condition
!---
!--- prepare bottom condition of spheromak type
!---
use common_val
implicit none
!integer :: i,j,k
real(8) :: rr=0d0, th=0d0
real(8) :: r0=1d0
real(8) :: yy0=at
real(8) :: ycd(ny0:ny)=0d0
real(8) :: xcd(nx0:nx)=0d0
real(8) :: rrd=0d0
real(8),allocatable :: Br(:,:), Bth(:,:), Bph(:,:)
integer,parameter :: test_ps = 0
integer,parameter :: mask_num = 0

allocate(Br(nx0:nx,ny0:ny))
allocate(Bth(nx0:nx,ny0:ny))
allocate(Bph(nx0:nx,ny0:ny))
Br=0d0; Bth=0d0; Bph=0d0

if(test_ps.eq.0) then
!----- Bz dist sech field ------
do i = nx0, nx
do j = ny0, ny/2
   Bz_init(i,j) = tanh(yc(j)/dtany) &
        & / cosh(sqrt((xc(i)/decx)**2 &
        & + ((yc(j)+lpilsec)/decy)**2))
end do
do j = ny/2, ny
   Bz_init(i,j) = tanh(yc(j)/dtany) &
        & / cosh(sqrt((xc(i)/decx)**2 &
        & + ((yc(j)-lpilsec)/decy)**2))
end do
!--- mask ambient field to increase kappa value
! if(mask_num.eq.1) then
!    if(yc(j).ge.0.0) then
!       Bz_init(i,j) = Bz_init(i,j) &
!            &       * (1.0 - tanh((yc(j) - pmask)/dmask))
!    else if(yc(j).le.0.0) then
!       Bz_init(i,j) = Bz_init(i,j) &
!            &       * (1.0 + tanh((yc(j) + pmask)/dmask))
!    end if
! end if
! enddo
enddo
else if(test_ps.eq.1) then
!----- Bz dist by point source -------------
do i = nx0, nx
!   xcd(i) = xc(i) * 0.25d0
!   xcd(i) = xc(i) * 0.125d0
   xcd(i) = xc(i)
do j = ny0, ny
   ycd(j) = yc(j)
!   ycd(j) = yc(j) * 0.25d0
!-----
   ! if(yc(j).gt.lpil) then
   !    ycd(j) = yc(j) * 0.5d0 + 0.5d0 * lpil
   ! else if(yc(j).lt.-lpil) then
   !    ycd(j) = yc(j) * 0.5d0 - 0.5d0 * lpil
   ! else
   !    ycd(j) = yc(j)
   ! end if
   Bz_init(i,j) = qps*dep/(xcd(i)**2+(ycd(j)-lpil)**2+dep**2)**(1.5) &
        &       - qps*dep/(xcd(i)**2+(ycd(j)+lpil)**2+dep**2)**(1.5)
!--- mask ambient field to increase kappa value
   ! if(ycd(j).ge.0.0) then
   !    Bz_init(i,j) = Bz_init(i,j) &
   !         &       * (1.0 - tanh((ycd(j) - 1.5)/0.3))
   ! else if(ycd(j).le.0.0) then
   !    Bz_init(i,j) = Bz_init(i,j) &
   !         &       * (1.0 + tanh((ycd(j) + 1.5)/0.3))
   ! end if
end do
end do
else
!----- Bz of spheromack field---------------
!+++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++
if(flatten_switch.eq.0) then
do i = nx0, nx
do j = ny0, ny
!	rr = (xc(i)**2 + yc(j)**2)**(0.5)
	rr = (xc(i)**2 + yc(j)**2 + zcut**2)**(0.5)
!--- because Bphi_init_yz diverge for rtmp*(-1)
	if(rr==0) then
		Bz_init(i,j) = 0d0
!--- vanish around magnetic field
	else if(rr .gt. a_suf) then
		Bz_init(i,j) = 0d0
	else
		th = asin(yc(j)/rr)
!                th = asin(xc(i)/rr)
!--------------- normal spheromak field
                ! Br(i,j)  = -2d0 * Bs0 * cos(th) &
                !      &* (sin(alp*rr) - alp*rr * cos(alp*rr)) / ((alp*rr)**3)
                ! Bth(i,j) = Bs0 * sin(th) &
                !      &* (sin(alp*rr) - alp*rr * cos(alp*rr) &
                !      & + 2d0*(alp*rr)*cos(alp*rr) &
                !      & + ((alp*rr)**2 - 2d0)*sin(alp*rr)) &
                !      & / ((alp*rr)**3)
		Bph(i,j) = Bs0 * sin(th) &
                     &* (sin(alp*rr) - alp*rr * cos(alp*rr)) / ((alp*rr)**2)
!		Bz_init(i,j) = Bs0 * sin(th) *(sin(alp*rr) - alp*rr * cos(alp*rr)) / ((alp*rr)**2)
!--------------- decay for exp(-rr*rr)
!		Bz_init(i,j) = exp(-rr*rr/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-xc*xc)
!		Bz_init(i,j) = exp(-xc(i)**2/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-yc*yc)
!		Bz_init(i,j) = exp(-yc(i)**2/(yy0*yy0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
	endif
enddo
enddo
!Bz_init = (Br**2 + Bth**2 + Bph**2)**0.5d0
Bz_init = Bph
!+++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++
!--- extending to x and y direction
!--- 'xcd' and 'ycd' are coordinats 
!---   for normal Bz dist.
!+++++++++++++++++++++++++++++++++++++++++++
else if(flatten_switch.eq.1) then
do i = nx0, nx
   xcd(i) = xc(i) * 0.5d0
!   xcd(i) = xc(i) * 0.25d0
enddo
do j = ny0, ny
   ycd(j) = yc(j) * 1d0
enddo
!
do i = nx0, nx
do j = ny0, ny
	rr = (xc(i)**2 + yc(j)**2)**(0.5)
!	rrd = (xcd(i)**2 + yc(j)**2)**(0.5)
	rrd = (xcd(i)**2 + ycd(j)**2)**(0.5)
!--- because Bphi_init_yz diverge for rtmp*(-1)
	if(rrd==0) then
		Bz_init(i,j) = 0d0
!--- vanish around magnetic field
	else if(rrd .gt. a_suf) then
		Bz_init(i,j) = 0d0
	else
!		th = asin(yc(j)/rrd)
		th = asin(ycd(j)/rrd)
!--------------- decay for exp(-rr*rr)
!		Bz_init(i,j) = exp(-rr*rr/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-xc*xc)
!		Bz_init(i,j) = exp(-xc(i)**2/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-yc*yc)
!		Bz_init(i,j) = exp(-yc(i)**2/(yy0*yy0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- flattend spheromak field
!		Bz_init(i,j) = Bs0 * sin(th) *(sin(alp*rrd) - alp*rrd * cos(alp*rrd)) / ((alp*rrd)**2)
		Bz_init(i,j) =(Bs0 * sin(th) *(sin(alp*rrd) - alp*rrd * cos(alp*rrd)) / ((alp*rrd)**2))&
       & * exp(-yc(j)**2/(ydecay**2))
	endif
enddo
enddo
!
!+++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++
else if(flatten_switch.eq.2) then
do j = ny0, ny
if(yc(j).ge.0) then
!   ycd(j) = yc(j)*(at*(2.0-tanh(yc(j))) - at + 1)
!   ycd(j) = yc(j)*(at*(2.0-tanh(yc(j)))/3d0 -(at/3d0-1))
   ycd(j) = yc(j)*(at*(2.0-tanh(2.0*yc(j)/3.0))/3.0 -(at/3.0-1))
else
!   ycd(j) = yc(j)*(at*(2.0+tanh(yc(j))) - at + 1)
!   ycd(j) = yc(j)*(at*(2.0+tanh(yc(j)))/3d0 -(at/3d0-1))
   ycd(j) = yc(j)*(at*(2.0+tanh(2.0*yc(j)/3.0))/3.0 -(at/3.0-1))
endif
enddo
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
do i = nx0, nx
do j = ny0, ny
	rr = (xc(i)**2 + yc(j)**2)**(0.5)
	rrd = (xc(i)**2 + ycd(j)**2)**(0.5)
!--- because Bphi_init_yz diverge for rtmp*(-1)
	if(rrd==0) then
		Bz_init(i,j) = 0d0
!--- vanish around magnetic field
	else if(rrd .gt. a_suf) then
		Bz_init(i,j) = 0d0
	else
		th = asin(ycd(j)/rrd)
!--------------- decay for exp(-rr*rr)
!		Bz_init(i,j) = exp(-rr*rr/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-xc*xc)
!		Bz_init(i,j) = exp(-xc(i)**2/(r0*r0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- decay for exp(-yc*yc)
!		Bz_init(i,j) = exp(-yc(i)**2/(yy0*yy0)) * Bs0 * sin(th) *(sin(rr) - rr * cos(rr)) / (rr**2)
!--------------- flattend spheromak field no.2
		Bz_init(i,j) = Bs0 * sin(th) *(sin(alp*rrd) - alp*rrd * cos(alp*rrd)) / ((alp*rrd)**2)
!		Bz_init(i,j) = Bz_init(i,j)*exp(-yc(j)**2/(yy0*yy0))
	endif
enddo
enddo
!
else
print *, "flatten_switch is abnormal value!!!"
print *, "flatten_switch = ", flatten_switch
stop
endif
!
end if
!--- normalize Bz
Bz_init(:,:) = 1d0*Bz_init(:,:)/maxval(abs(Bz_init(:,:)))
!
end
