program test_foot_point
implicit none
real(8),allocatable :: inifx(:,:), inify(:,:), inifz(:,:)
real(8),allocatable :: endfxp(:,:), endfyp(:,:), endfzp(:,:)
real(8),allocatable :: endfxm(:,:), endfym(:,:), endfzm(:,:)
real(8),allocatable :: leng_ini(:,:), leng(:,:)
real(8),allocatable :: leng_ini_p(:,:), leng_ini_m(:,:)
real(8),allocatable :: leng_p(:,:), leng_m(:,:)
real(8),allocatable :: diff_L(:,:,:)
integer :: ts, i, j, tsini, tsfin
integer,allocatable :: ckpm(:,:)
integer :: nxmin, nxmax, nx, nymin, nymax, ny, nzmin, nzmax
!----- parameter ---------
real(8),parameter ::r_em=0.4
real(8),parameter :: xcen = 0.0, ycen = 0.0
real(8),parameter :: r_ar = 6.0
real(8),parameter :: z_aj = 0.01
!real(8),parameter :: z_aj = 0.0235
!-------------------------
integer,allocatable :: nopp(:,:), nopm(:,:)
character :: cts*3

!--- get info ---
open(100,file="../../../../code_output/info",status="old")
read(100,*) tsini; read(100,*) tsfin
read(100,*) nxmin; read(100,*) nxmax
read(100,*) nymin; read(100,*) nymax
read(100,*) nzmin; read(100,*) nzmax
close(100)
nx = nxmax - nxmin; ny = nymax - nymin

allocate(inifx(0:nx,0:ny)); inifx=0
allocate(inify(0:nx,0:ny)); inify=0
allocate(inifz(0:nx,0:ny)); inifz=0
allocate(endfxp(0:nx,0:ny)); endfxp=0
allocate(endfyp(0:nx,0:ny)); endfyp=0
allocate(endfzp(0:nx,0:ny)); endfzp=0
allocate(endfxm(0:nx,0:ny)); endfxm=0
allocate(endfym(0:nx,0:ny)); endfym=0
allocate(endfzm(0:nx,0:ny)); endfzm=0
allocate(leng_ini(0:nx,0:ny)); leng_ini=0
allocate(leng(0:nx,0:ny)); leng=0
allocate(leng_ini_p(0:nx,0:ny)); leng_ini_p=0
allocate(leng_p(0:nx,0:ny)); leng_p=0
allocate(leng_ini_m(0:nx,0:ny)); leng_ini_m=0
allocate(leng_m(0:nx,0:ny)); leng_m=0
allocate(diff_L(0:nx,0:ny,tsini:tsfin)); diff_L=0
allocate(nopp(0:nx,0:ny)); nopp=0
allocate(nopm(0:nx,0:ny)); nopm=0
allocate(ckpm(0:nx,0:ny)); ckpm=0

!--- open seed file
!
open(21,file="seedall.txt",status="replace")

!--- initial foot point
!
ts=tsini
write(cts,"(I3.3)") ts
write(*,*) cts

open(70,file="../../Bfootp1"//cts, status="old")
do i = 0, nx
do j = 0, ny
read(70,*) inifx(i,j), inify(i,j), inifz(i,j)
read(70,*) endfxp(i,j), endfyp(i,j), endfzp(i,j)
read(70,*) nopp(i,j)
enddo 
enddo
close(70)

open(71,file="../../Bfootm1"//cts, status="old")
do i = 0, nx
do j = 0, ny
read(71,*) inifx(i,j), inify(i,j), inifz(i,j)
read(71,*) endfxm(i,j), endfym(i,j), endfzm(i,j)
read(71,*) nopm(i,j)
enddo 
enddo
close(71)

do i = 0, nx
do j = 0, ny
leng_ini_p(i,j) = sqrt((endfxp(i,j) - inifx(i,j))**2 &
     &          +      (endfyp(i,j) - inify(i,j))**2 &
     &          +      (endfzp(i,j) - inifz(i,j))**2)
leng_ini_m(i,j) = sqrt((endfxm(i,j) - inifx(i,j))**2 &
     &          +      (endfym(i,j) - inify(i,j))**2 &
     &          +      (endfzm(i,j) - inifz(i,j))**2)
if(leng_ini_p(i,j) .gt. leng_ini_m(i,j)) then
   leng_ini(i,j) = leng_ini_p(i,j)
else
   leng_ini(i,j) = leng_ini_m(i,j)
end if
!
diff_L(i,j,tsini) = leng_ini(i,j) - leng_ini(i,j)
enddo
enddo

open(20,file="diff_L"//cts,status="replace")
do i = 0, nx
do j = 0, ny
write(20,*) diff_L(i,j,tsini)
enddo
enddo
close(20)

!write(21,'(f10.5,2x,f10.5,2x,f10.5,i3)') xcen, ycen, inifz(i,j)+z_aj, tsini-9

!!---------------------------------
!!=================================
!------ begin do loop 
!------ calculate footpoint and distance from init position
!
do ts = tsini+1, tsfin
write(cts,'(I3.3)') ts
write(*,*) "cts =", cts

open(70,file="../../Bfootp1"//cts, status="old")
do i = 0, nx
do j = 0, ny
read(70,*) inifx(i,j), inify(i,j), inifz(i,j)
read(70,*) endfxp(i,j), endfyp(i,j), endfzp(i,j)
read(70,*) nopp(i,j)
enddo 
enddo
close(70)

open(71,file="../../Bfootm1"//cts, status="old")
do i = 0, nx
do j = 0, ny
read(71,*) inifx(i,j), inify(i,j), inifz(i,j)
read(71,*) endfxm(i,j), endfym(i,j), endfzm(i,j)
read(71,*) nopm(i,j)
enddo 
enddo
close(71)

do i = 0, nx
do j = 0, ny
leng_p(i,j) = sqrt((endfxp(i,j) - inifx(i,j))**2 &
     &      +      (endfyp(i,j) - inify(i,j))**2 &
     &      +      (endfzp(i,j) - inifz(i,j))**2)
leng_m(i,j) = sqrt((endfxm(i,j) - inifx(i,j))**2 &
     &      +      (endfym(i,j) - inify(i,j))**2 &
     &      +      (endfzm(i,j) - inifz(i,j))**2)
if(leng_p(i,j) .ge. leng_m(i,j)) then
   leng(i,j) = leng_p(i,j)
   ckpm(i,j) = +1
else
   leng(i,j) = leng_m(i,j)
   ckpm(i,j) = -1
end if
!
! "ckpm" means the direction of magnetic field
!
if(((ckpm(i,j).eq.+1).and.(endfzp(i,j).le.0)).or.&
  &((ckpm(i,j).eq.-1).and.(endfzm(i,j).le.0))) then
   diff_L(i,j,ts) = leng(i,j) - leng_ini(i,j)
else
   diff_L(i,j,ts) = 0d0
endif
enddo
enddo

open(20,file="diff_L"//cts,status="replace")
do i = 0, nx
do j = 0, ny
write(20,*) diff_L(i,j,ts)
enddo
enddo
close(20)

do i = 0, nx
do j = 0, ny
!-----
! change of footpoint position is larger than
! the radius of emerging flux?
!-----
if(diff_L(i,j,ts).gt.r_em) then
!-----
! footpoint does not connect to the outside of AR?
!-----
if(sqrt((inifx(i,j)-xcen)**2+(inify(i,j)-ycen)**2).lt.r_ar) then
!-----
! footpoint does not connect to inner side of emerging flux?
!-----
if(sqrt((inifx(i,j)-xcen)**2+(inify(i,j)-ycen)**2).ge.r_em) then
!-----
! initial footpoint is in region 1 or 3?
!-----
if(((inifx(i,j).le.xcen).and.(inify(i,j).le.ycen))&
     &.or.&
  &((inifx(i,j).ge.xcen).and.(inify(i,j).ge.ycen))) then
!-----
! initial footpoint is inner region of 2 times 
! of the radius of emerging flux for y direction?
!-----
!if((inify(i,j)-ycen).le.2d0*r_em) then
!-----
! the direction of Bz is plus?
!-----
if(ckpm(i,j).eq.+1) then
   !-----
   ! end of footpoint connects to outside of emerging flux?
   !-----
   if(sqrt((endfxp(i,j)-xcen)**2+(endfyp(i,j)-ycen)**2).ge.r_em) then
   !-----
   ! end of footpoint is in region 1 or 3?
   !-----
   if(((endfxp(i,j).le.xcen).and.(endfyp(i,j).le.ycen))&
        &.or.&
     &((endfxp(i,j).ge.xcen).and.(endfyp(i,j).ge.ycen))) then
   !-----
   ! end footpoint is inner region of 2 times 
   ! of the radius of emerging flux for y direction?
   !-----
!   if((endfyp(i,j)-ycen).le.2d0*r_em) then
      write(21,'(f10.5,2x,f10.5,2x,f10.5,3x,i4)') inifx(i,j), inify(i,j), inifz(i,j)+z_aj, ts-tsini
!   end if
   endif
   endif
!-----
! the direction of Bz is minus?
!-----
else if(ckpm(i,j).eq.-1) then
   !-----
   ! end of footpoint connects to outside of emerging flux?
   !-----
   if(sqrt((endfxm(i,j)-xcen)**2+(endfym(i,j)-ycen)**2).ge.r_em) then
   !-----
   ! end of footpoint is in region 1 or 3?
   !-----
   if(((endfxm(i,j).le.xcen).and.(endfym(i,j).le.ycen))&
        &.or.&
     &((endfxm(i,j).ge.xcen).and.(endfym(i,j).ge.ycen))) then
   !-----
   ! end footpoint is inner region of 2 times 
   ! of the radius of emerging flux for y direction?
   !-----
!   if((endfym(i,j)-ycen).le.2d0*r_em) then
      write(21,'(f10.5,2x,f10.5,2x,f10.5,3x,i4)') inifx(i,j), inify(i,j), inifz(i,j)+z_aj, ts-tsini
!   end if
   endif
   endif
end if
!
!end if
end if
!
end if
end if
end if
enddo
enddo

end do

close(21)

end program test_foot_point
