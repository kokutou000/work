subroutine calc_velocity
use common_val
implicit none
!integer :: i,j,k
integer :: im, ip, jm, jp
real(8) :: dx, dy
real(8) :: dxp, dxm, dyp, dym
real(8) :: maxvv=0d0
real(8) :: maxbz=0d0

Vx_init=0d0
Vy_init=0d0
Vz_init=0d0
ck1Bz = 0d0
ck2Bz = 0d0
dx = xc(1) - xc(0)
dy = yc(1) - yc(0)
!
!
do i = nx0+1, nx-1
dxp = xc(i+1) - xc(i)
dxm = xc(i) - xc(i-1)
do j = ny0+1, ny-1
dyp = yc(j+1) - yc(j)
dym = yc(j) - yc(j-1)
im = i - 1
ip = i + 1
jm = j - 1
jp = j + 1
!if(i.eq.0) then
!	im = nx-1
!	ip = 1
!else if(i.eq.nx) then
!	im = nx-1
!	ip = 1
!endif
!if(j.eq.0) then
!	jm = ny-1
!	jp = 1
!else if(j.eq.ny) then
!	jm = ny-1
!	jp = 1
!endif
!
Vx_init(i,j) = tanh(Bz_init(i,j))*&
     &       ( Bz_init(i,jp)*dym*dym &
     &       + Bz_init(i,j)*(dyp*dyp-dym*dym)&
     &       - Bz_init(i,jm)*dyp*dyp)&
     &       /(dyp*dym*(dyp+dym))
Vy_init(i,j) =-tanh(Bz_init(i,j))*&
     &       ( Bz_init(ip,j)*dxm*dxm &
     &       + Bz_init(i,j)*(dxp*dxp-dxm*dxm)&
     &       - Bz_init(im,j)*dxp*dxp)&
     &       /(dxp*dxm*(dxp+dxm))
!Vx_init(i,j) = tanh(Bz_init(i,j))*(Bz_init(i,jp)-Bz_init(i,jm))/dy
!Vy_init(i,j) =-tanh(Bz_init(i,j))*(Bz_init(ip,j)-Bz_init(im,j))/dx
!
!if(Bz_init(i,j).ne.0) then
!Vx_init(i,j) = (Bz_init(i,j)/abs(Bz_init(i,j)))*(Bz_init(i,jp)-Bz_init(i,jm))/dy
!Vy_init(i,j) =-(Bz_init(i,j)/abs(Bz_init(i,j)))*(Bz_init(ip,j)-Bz_init(im,j))/dx
!else
!Vx_init(i,j) = (Bz_init(i,jp)-Bz_init(i,jm))/dy
!Vy_init(i,j) =-(Bz_init(ip,j)-Bz_init(im,j))/dx
!end if

enddo
enddo
!---
!----- mask velocity field not to twist strongly around center of Bz
!---
maxbz = maxval(abs(Bz_init))
!Vx_init = Vx_init * 0.5 * (1.0 - tanh(8.0*(abs(Bz_init/maxbz)-0.5)))
!Vy_init = Vy_init * 0.5 * (1.0 - tanh(8.0*(abs(Bz_init/maxbz)-0.5)))
!Vx_init = Vx_init * 0.5 * (1.0 - tanh(15.0*(abs(Bz_init/maxbz)-0.7)))
!Vy_init = Vy_init * 0.5 * (1.0 - tanh(15.0*(abs(Bz_init/maxbz)-0.7)))
!Vx_init = Vx_init * 0.5 * (1.0 - tanh(20.0*(abs(Bz_init/maxbz)-0.8)))
!Vy_init = Vy_init * 0.5 * (1.0 - tanh(20.0*(abs(Bz_init/maxbz)-0.8)))
! Vx_init = Vx_init * 0.5 * (1.0 - tanh((abs(Bz_init/maxbz)-0.95)/0.025))
! Vy_init = Vy_init * 0.5 * (1.0 - tanh((abs(Bz_init/maxbz)-0.95)/0.025))
maxvv = maxval(sqrt(Vx_init**2+Vy_init**2))
Vx_init = Vx_init/maxvv
Vy_init = Vy_init/maxvv
Vx_init = Vx_init*0.01d0
Vy_init = Vy_init*0.01d0
!Vx_init = Vx_init*0.1d0
!Vy_init = Vy_init*0.1d0

!----------
!----------
!----------
! do i = nx0, nx
! do j = ny0, ny
! im = i - 1
! ip = i + 1
! jm = j - 1
! jp = j + 1
! if(i.eq.nx0) then
! 	im = nx-1
! 	ip = 0
! else if(i.eq.nx) then
! 	im = nx-1
! 	ip = 0
! endif
! if(j.eq.ny0) then
! 	jm = ny-1
! 	jp = 0
! else if(j.eq.ny) then
! 	jm = ny-1
! 	jp = 0
! endif

! ck1Bz(i,j) = Vx_init(i,j)*(Bz_init(ip,j)-Bz_init(im,j))/dx &
!      &   + Vy_init(i,j)*(Bz_init(i,jp)-Bz_init(i,jm))/dy
! ck2Bz(i,j) = Bz_init(i,j)*((Vx_init(ip,j)-Vx_init(im,j))/dx &
!      &   + (Vy_init(i,jp)-Vy_init(i,jm))/dy)
! !delBzinit(i,j) = ck1Bz(i,j) + ck2Bz(i,j)

! enddo
! enddo

end
