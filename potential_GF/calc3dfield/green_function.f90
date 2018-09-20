subroutine green_function_method
use common_val
implicit none
!integer :: i,j,k
integer,parameter :: ompnum=28 !--- number for open mp parallel
integer,parameter :: klower=27
!--- switching parameter for only limited calculation --- 1 -> yes
integer,parameter :: calc_limited = 0
!--- switching parameter under the xy-symmetry condition --- 1 -> yes
integer,parameter :: symmetry_xy = 1
!--- params for extended coordinates
integer,parameter :: nx0d = nx0-3, nx1d = nx+3
integer,parameter :: ny0d = ny0-3, ny1d = ny+3
integer,parameter :: nz0d = nz0, nz1d = nz
real(8) :: xcd(nx0d:nx1d)=0d0, ycd(ny0d:ny1d)=0d0, zcd(nz0d:nz1d)=0d0
integer :: outn = 4
integer :: t1, t2, t_rate, t_max
! for integration by alternately grid system
real(8),allocatable :: xcme(:),ycme(:), Bzme(:,:)
integer,parameter :: nx0me=nx0d, nx1me=nx1d+1
integer,parameter :: ny0me=ny0d, ny1me=ny1d+1
! for integration by alternately grid system, more fine case
real(8),allocatable :: xcme2(:),ycme2(:), Bzme2(:,:)
integer,parameter :: nx0me2=2*nx0d, nx1me2=2*nx1d+1
integer,parameter :: ny0me2=2*ny0d, ny1me2=2*ny1d+1

!---------------------------------------------------------
!-- Make extended coordinates xcd, ycd, zcd
!-- for better calculation around xy-boundary
do i = nx0d, nx1d
   xcd(i) = alphax*(i-nx/2) + betax*(i-nx/2)**nbeki
!   xcd(i) = (i-nx/2)*delx      !!!--- for debug uniform grid
end do
do j = ny0d, ny1d
   ycd(j) = alphay*(j-ny/2) + betay*(j-ny/2)**nbeki
!   ycd(j) = (j-ny/2)*dely      !!!--- for debug uniform grid
end do
do k = nz0d, nz1d
   zcd(k) = alphaz*k + betaz*k**nbeki
!   zcd(k) = k*delz         !!!--- for debug uniform grid
end do
!--- output for check of grid system -----
open(101,file=dirname//"/xcd.dat",status="replace")
do i = nx0d, nx1d
   write(101,*) xcd(i)
enddo
close(101)
open(102,file=dirname//"/ycd.dat",status="replace")
do j = ny0d, ny1d
   write(102,*) ycd(j)
enddo
close(102)
open(103,file=dirname//"/zcd.dat",status="replace")
do k = nz0d, nz1d
   write(103,*) zcd(k)
enddo
close(103)
!-------------------------------------------------------
allocate(xcme(nx0me:nx1me)); allocate(ycme(ny0me:ny1me))
allocate(Bzme(nx0me:nx1me,ny0me:ny1me))
!
xcme(nx0me) = xcd(nx0d) - 0.5d0*(xcd(nx0d+1)-xcd(nx0d))
xcme(nx1me) = xcd(nx1d) + 0.5d0*(xcd(nx1d)-xcd(nx1d-1))
do i = nx0me+1, nx1me-1
!   xcme(i) = 0.5d0*(xc(i) + xc(i+1))
   xcme(i) = 0.5d0*(xcd(i-1) + xcd(i))
enddo
ycme(ny0me) = ycd(ny0d) - 0.5d0*(ycd(ny0d+1)-ycd(ny0d))
ycme(ny1me) = ycd(ny1d) + 0.5d0*(ycd(ny1d)-ycd(ny1d-1))
do j = ny0me+1, ny1me-1
!   ycme(j) = 0.5d0*(yc(j) + yc(j+1))
   ycme(j) = 0.5d0*(ycd(j-1) + ycd(j))
enddo
Bzme=0d0
do i = nx0me, nx1me
! do j = ny0, ny/2
!    Bzme(i,j) = tanh(ycme(j)/dtany) &
!         & /cosh(sqrt((xcme(i)/decx)**2 &
!         & + ((ycme(j)+lpilsec)/decy)**2 ))
! end do
! do j = 1+ny/2, ny+1
!       Bzme(i,j) = tanh(ycme(j)/dtany) &
!         & /cosh(sqrt((xcme(i)/decx)**2 &
!         & + ((ycme(j)-lpilsec)/decy)**2 ))
! end do
do j = ny0me, ny1me
!   Bzme(i,j) = qps*dep*((xcme(i)**2+(ycme(j)-lpil)**2+dep**2)**(-1.5)) &
!        &    - qps*dep*((xcme(i)**2+(ycme(j)+lpil)**2+dep**2)**(-1.5))
   if(ycme(j).le.0d0) then
      Bzme(i,j) = tanh(ycme(j)/dtany) &
           & /cosh(sqrt((xcme(i)/decx)**2 &
           & + ((ycme(j)+lpilsec)/decy)**2 ))
   else
      Bzme(i,j) = tanh(ycme(j)/dtany) &
           & /cosh(sqrt((xcme(i)/decx)**2 &
           & + ((ycme(j)-lpilsec)/decy)**2 ))
   end if
end do
end do
Bzme = 1d0*Bzme/maxval(abs(Bzme))
!--- check grid system -----
open(71,file=dirname//"/xcme.dat",status="replace")
do i = nx0me, nx1me
   write(71,*) xcme(i)
enddo
close(71)
!------ more fine grid system for integration of green's function ----
allocate(xcme2(nx0me2:nx1me2)); allocate(ycme2(ny0me2:ny1me2))
allocate(Bzme2(nx0me2:nx1me2,ny0me2:ny1me2))
do i = nx0d, nx1d
   xcme2(2*i  ) = 0.5d0*(xcd(i)+xcme(i))
   xcme2(2*i+1) = 0.5d0*(xcd(i)+xcme(i+1))
end do
do j = ny0d, ny1d
   ycme2(2*j  ) = 0.5d0*(ycd(j)+ycme(j))
   ycme2(2*j+1) = 0.5d0*(ycd(j)+ycme(j+1))
end do
Bzme2 = 0d0
do i = nx0me2, nx1me2
do j = ny0me2, ny1me2
!   Bzme2(i,j) = qps*dep*((xcme2(i)**2+(ycme2(j)-lpil)**2+dep**2)**(-1.5)) &
!        &     - qps*dep*((xcme2(i)**2+(ycme2(j)+lpil)**2+dep**2)**(-1.5))
   if(ycme2(j).le.0d0) then
      Bzme2(i,j) = tanh(ycme2(j)/dtany) &
           & /cosh(sqrt((xcme2(i)/decx)**2 &
           & + ((ycme2(j)+lpilsec)/decy)**2 ))
   else
      Bzme2(i,j) = tanh(ycme2(j)/dtany) &
           & /cosh(sqrt((xcme2(i)/decx)**2 &
           & + ((ycme2(j)-lpilsec)/decy)**2 ))
   end if
end do
end do
Bzme2 = 1d0*Bzme2/maxval(abs(Bzme2))
!--- check grid system -----
open(81,file=dirname//"/xcme2.dat",status="replace")
do i = nx0me2, nx1me2
   write(81,*) xcme2(i)
enddo
close(81)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print *, "=========="
print *, "calculate potential field"
!
call system_clock(t1)
!------------------------------------------------------
call integration_greenfunc_omp_lower
call integration_greenfunc_omp_upper
!   call integration_greenfunc_omp_lower(ompnum)
!   call integration_greenfunc_omp_upper(ompnum)
!------------------------------------------------------
!----- copy Bxyz(nx0:nx/2,ny0:ny/2,:) -->> other region
!--- 3 to 2
do i = nx0, nx/2
do j = ny/2, ny
   Bx(i,j,:) = -Bx(i,ny-j,:)
   By(i,j,:) =  By(i,ny-j,:)
   Bz(i,j,:) = -Bz(i,ny-j,:)
enddo
enddo
!--- 3 to 4
do i = nx/2, nx
do j = ny0, ny/2
   Bx(i,j,:) = -Bx(nx-i,j,:)
   By(i,j,:) =  By(nx-i,j,:)
   Bz(i,j,:) =  Bz(nx-i,j,:)
enddo
enddo
!--- 3 to 1
do i = nx/2, nx
do j = ny/2, ny
   Bx(i,j,:) =  Bx(nx-i,ny-j,:)
   By(i,j,:) =  By(nx-i,ny-j,:)
   Bz(i,j,:) = -Bz(nx-i,ny-j,:)
enddo
enddo
!
!--- Bx
! Bx(nx0:nx/2,ny:ny/2 ,:) = -Bx(nx0:nx/2,ny0+2:ny/2,:) !!! 3 to 2
! Bx(nx:nx/2 ,ny0:ny/2,:) = -Bx(nx0+2:nx/2,ny0:ny/2,:) !!! 3 to 4
! Bx(nx:nx/2 ,ny:ny/2 ,:) =  Bx(nx0+2:nx/2,ny0+2:ny/2,:) !!! 3 to 1
! !--- By
! By(nx0:nx/2,ny:ny/2 ,:) =  By(nx0:nx/2,ny0+2:ny/2,:) !!! 3 to 2
! By(nx:nx/2 ,ny0:ny/2,:) =  By(nx0+2:nx/2,ny0:ny/2,:) !!! 3 to 4
! By(nx:nx/2 ,ny:ny/2 ,:) =  By(nx0+2:nx/2,ny0+2:ny/2,:) !!! 3 to 1
! !--- Bz
! Bz(nx0:nx/2,ny:ny/2 ,:) = -Bz(nx0:nx/2,ny0+2:ny/2,:) !!! 3 to 2
! Bz(nx:nx/2 ,ny0:ny/2,:) =  Bz(nx0+2:nx/2,ny0:ny/2,:) !!! 3 to 4
! Bz(nx:nx/2 ,ny:ny/2 ,:) = -Bz(nx0+2:nx/2,ny0+2:ny/2,:) !!! 3 to 1
!------------------------------------------------------
call system_clock(t2,t_rate,t_max)
write(*,*) "t1, t2, t_rate:", t1, t2, t_rate
write(*,*) "calc_time:", dble(t2 - t1)/dble(t_rate)
!
Bz(:,:,nz0) = Bz_init(:,:)
!
print *, "...finish"
print *, "-----------"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!
!subroutine integration_greenfunc_omp_lower(ompnum,klower)
subroutine integration_greenfunc_omp_lower
implicit none
real(8) :: bxdct(ny0me2:ny1me2)=0d0
real(8) :: bydct(ny0me2:ny1me2)=0d0
real(8) :: bzdct(ny0me2:ny1me2)=0d0
real(8) :: Totx = 0d0, Toty = 0d0, Totz = 0d0
real(8) :: rz2 = 0d0
integer :: i0=0, j0=0
integer :: nxbeg=0, nxfin=0
!$omp parallel num_threads(ompnum)
!$omp do private(rz2,Totx,Toty,Totz,k,j,i,j0,i0,bxdct,bydct,bzdct,nxbeg,nxfin)
do k = nz0, klower
   rz2 = (zc(k) - zc(nz0))**2
!do j = ny0, ny
do j = ny0, ny/2
if(calc_limited.eq.1) then
!---- minimum test calculation for check of decay index ----
if((j.ge.(ny/2)-2).and.(j.le.(ny/2)+2)) then
   nxbeg = nx0; nxfin = nx/2
else
   nxbeg = (nx/2)-2; nxfin = (nx/2)+2
end if
!-----------------------------------------------------------
else
!   nxbeg = nx0; nxfin = nx
   nxbeg = nx0; nxfin = nx/2
end if
do i = nxbeg, nxfin
!
!-----
! integration by trapezoidal rule
!-----
   Totx = 0d0; Toty = 0d0; Totz = 0d0
   bxdct = 0d0; bydct = 0d0; bzdct = 0d0
   do j0 = ny0me2, ny1me2
      do i0 = nx0me2, nx1me2-1
         bxdct(j0) = bxdct(j0) + &
              & (Bzme2(i0,j0)*(xc(i)-xcme2(i0)) &
              &  * ((xc(i)-xcme2(i0))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5) &
              & +Bzme2(i0+1,j0)*(xc(i)-xcme2(i0+1)) &
              &  * ((xc(i)-xcme2(i0+1))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5)) &
              & *(xcme2(i0+1)-xcme2(i0))
         bydct(j0) = bydct(j0) + &
              & (Bzme2(i0,j0) &
              &  * ((xc(i)-xcme2(i0))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5) &
              & +Bzme2(i0+1,j0) &
              &  * ((xc(i)-xcme2(i0+1))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5)) &
              & *(xcme2(i0+1)-xcme2(i0))*(yc(j)-ycme2(j0))
         bzdct(j0) = bzdct(j0) + &
              & (Bzme2(i0,j0) &
              &  * ((xc(i)-xcme2(i0))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5) &
              & +Bzme2(i0+1,j0) &
              &  * ((xc(i)-xcme2(i0+1))**2+(yc(j)-ycme2(j0))**2+rz2)**(-1.5)) &
              & *(xcme2(i0+1)-xcme2(i0))*(zc(k)-zc(nz0))
      enddo
   enddo
   do j0 = ny0me2, ny1me2-1
      Totx = Totx + &
           & (bxdct(j0+1) + bxdct(j0))*(ycme2(j0+1)-ycme2(j0))
      Toty = Toty + &
           & (bydct(j0+1) + bydct(j0))*(ycme2(j0+1)-ycme2(j0))
      Totz = Totz + &
           & (bzdct(j0+1) + bzdct(j0))*(ycme2(j0+1)-ycme2(j0))
   enddo
!!! due to of trapezodal rule and coef(greenfunc), for fast calclation
   Bx(i,j,k) = 0.25d0*Totx/(2d0*pi)
   By(i,j,k) = 0.25d0*Toty/(2d0*pi)
   Bz(i,j,k) = 0.25d0*Totz/(2d0*pi)
!
end do
end do
   if(mod(k,outn).eq.0) then
      print *, "--------"
      print *, "now k =", k
   endif
enddo
!$omp end do
!$omp end parallel
end subroutine integration_greenfunc_omp_lower
!
!
subroutine integration_greenfunc_omp_upper
implicit none
real(8) :: bxdct(ny0me:ny1me)=0d0
real(8) :: bydct(ny0me:ny1me)=0d0
real(8) :: bzdct(ny0me:ny1me)=0d0
real(8) :: Totx = 0d0, Toty = 0d0, Totz = 0d0
real(8) :: rz2 = 0d0
integer :: i0=0, j0=0
integer :: nxbeg=0, nxfin=0
!subroutine integration_greenfunc_omp_upper(ompnum,klower)
!$omp parallel num_threads(ompnum)
!$omp do private(rz2,Totx,Toty,Totz,k,j,i,j0,i0,bxdct,bydct,bzdct,nxbeg,nxfin)
do k = klower+1, nz
   rz2 = (zc(k) - zc(nz0))**2
!do j = ny0, ny
do j = ny0, ny/2
if(calc_limited.eq.1) then
!---- minimum test calculation for check of decay index ----
if((j.ge.(ny/2)-2).and.(j.le.(ny/2)+2)) then
   nxbeg = nx0; nxfin = nx/2
else
   nxbeg = (nx/2)-2; nxfin = (nx/2)+2
end if
!-----------------------------------------------------------
else
!   nxbeg = nx0; nxfin = nx
   nxbeg = nx0; nxfin = nx/2
end if
do i = nxbeg, nxfin
!
!-----
! integration by trapezoidal rule
!-----
   Totx = 0d0; Toty = 0d0; Totz = 0d0
   bxdct = 0d0; bydct = 0d0; bzdct = 0d0
   do j0 = ny0me, ny1me
      do i0 = nx0me, nx1me-1
         bxdct(j0) = bxdct(j0) + &
              & (Bzme(i0,j0)*(xc(i)-xcme(i0)) &
              &  * ((xc(i)-xcme(i0))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5) &
              & +Bzme(i0+1,j0)*(xc(i)-xcme(i0+1)) &
              &  * ((xc(i)-xcme(i0+1))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5)) &
              & *(xcme(i0+1)-xcme(i0))
         bydct(j0) = bydct(j0) + &
              & (Bzme(i0,j0) &
              &  * ((xc(i)-xcme(i0))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5) &
              & +Bzme(i0+1,j0) &
              &  * ((xc(i)-xcme(i0+1))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5)) &
              & *(xcme(i0+1)-xcme(i0))*(yc(j)-ycme(j0))
         bzdct(j0) = bzdct(j0) + &
              & (Bzme(i0,j0) &
              &  * ((xc(i)-xcme(i0))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5) &
              & +Bzme(i0+1,j0) &
              &  * ((xc(i)-xcme(i0+1))**2+(yc(j)-ycme(j0))**2+rz2)**(-1.5)) &
              & *(xcme(i0+1)-xcme(i0))*(zc(k)-zc(nz0))
      enddo
   enddo
   do j0 = ny0me, ny1me-1
      Totx = Totx + &
           & (bxdct(j0+1) + bxdct(j0))*(ycme(j0+1)-ycme(j0))
      Toty = Toty + &
           & (bydct(j0+1) + bydct(j0))*(ycme(j0+1)-ycme(j0))
      Totz = Totz + &
           & (bzdct(j0+1) + bzdct(j0))*(ycme(j0+1)-ycme(j0))
   enddo
!!! due to of trapezodal rule and coef(greenfunc), for fast calclation
   Bx(i,j,k) = 0.25d0*Totx/(2d0*pi)
   By(i,j,k) = 0.25d0*Toty/(2d0*pi)
   Bz(i,j,k) = 0.25d0*Totz/(2d0*pi)
!
end do
end do
   if(mod(k,outn).eq.0) then
      print *, "--------"
      print *, "now k =", k
   endif
enddo
!$omp end do
!$omp end parallel
end subroutine integration_greenfunc_omp_upper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine green_function_method
