subroutine calc_current_divb
use common_val
implicit none
real(8),allocatable :: cx(:,:,:), cy(:,:,:), cz(:,:,:)
real(8),allocatable :: divb(:,:,:)
! for debug divb calculation
real(8),allocatable :: dbxdx(:,:,:), dbydy(:,:,:), dbzdz(:,:,:)
real(8) :: dxp, dxm, dyp, dym, dzp, dzm
!
allocate(cx(nx0:nx,ny0:ny,nz0:nz)); cx=0d0
allocate(cy(nx0:nx,ny0:ny,nz0:nz)); cy=0d0
allocate(cz(nx0:nx,ny0:ny,nz0:nz)); cz=0d0
allocate(divb(nx0:nx,ny0:ny,nz0:nz)); divb=0d0
allocate(dbxdx(nx0:nx,ny0:ny,nz0:nz)); dbxdx=0d0
allocate(dbydy(nx0:nx,ny0:ny,nz0:nz)); dbydy=0d0
allocate(dbzdz(nx0:nx,ny0:ny,nz0:nz)); dbzdz=0d0
dxp=0d0; dxm=0d0; dyp=0d0; dym=0d0; dzp=0d0; dzm=0d0
!
print *, "calculate electric current for check of B"
cx = 0d0
cy = 0d0
cz = 0d0
divb = 0d0
dbxdx = 0d0; dbydy = 0d0; dbzdz = 0d0
do i = nx0+1, nx-1
   dxp = xc(i+1) - xc(i)
   dxm = xc(i) - xc(i-1)
do j = ny0+1, ny-1
   dyp = yc(j+1) - yc(j)
   dym = yc(j) - yc(j-1)
do k = nz0, nz-1
   if(k.ne.0) then
      dzp = zc(k+1) - zc(k)
      dzm = zc(k) - zc(k-1)
   endif
if(k.eq.nz0) then
cx(i,j,k) = (dym*dym*bz(i,j+1,k) + &
     & (dyp*dyp-dym*dym)*bz(i,j,k) - &
     & dyp*dyp*bz(i,j-1,k)) / &
     & (dyp*dym*(dyp+dym)) - &
!     & (by(i,j,k+1)-by(i,j,k))/(zc(k+1)-zc(k))
     & ((by(i,j,k+2)*(zc(k+1)-zc(k))**2) &
     & -(by(i,j,k+1)*(zc(k+2)-zc(k))**2) &
     & +(by(i,j,k)*((zc(k+2)-zc(k))**2-(zc(k+1)-zc(k))**2))) &
     & /((zc(k+1)-zc(k))*(zc(k+2)-zc(k))*(zc(k+1)-zc(k+2)))
!cy(i,j,k) = (bx(i,j,k+1)-bx(i,j,k))/(zc(k+1)-zc(k)) - &
cy(i,j,k) = ((bx(i,j,k+2)*(zc(k+1)-zc(k))**2) &
     & -(bx(i,j,k+1)*(zc(k+2)-zc(k))**2) &
     & +(bx(i,j,k)*((zc(k+2)-zc(k))**2-(zc(k+1)-zc(k))**2))) &
     & /((zc(k+1)-zc(k))*(zc(k+2)-zc(k))*(zc(k+1)-zc(k+2))) - &
     & (dxm*dxm*bz(i+1,j,k) + &
     & (dxp*dxp-dxm*dxm)*bz(i,j,k) - &
     & dxp*dxp*bz(i-1,j,k)) / &
     & (dxp*dxm*(dxp+dxm))
dbzdz(i,j,k) = (bz(i,j,k+1)-bz(i,j,k))/(zc(k+1)-zc(k))
else
cx(i,j,k) = (dym*dym*bz(i,j+1,k) + &
     & (dyp*dyp-dym*dym)*bz(i,j,k) - &
     & dyp*dyp*bz(i,j-1,k)) / &
     & (dyp*dym*(dyp+dym)) - &
     & (dzm*dzm*by(i,j,k+1) + &
     & (dzp*dzp-dzm*dzm)*by(i,j,k) - &
     & dzp*dzp*by(i,j,k-1)) / &
     & (dzp*dzm*(dzp+dzm))
cy(i,j,k) = (dzm*dzm*bx(i,j,k+1) + &
     & (dzp*dzp-dzm*dzm)*bx(i,j,k) - &
     & dzp*dzp*bx(i,j,k-1)) / &
     & (dzp*dzm*(dzp+dzm)) - &
     & (dxm*dxm*bz(i+1,j,k) + &
     & (dxp*dxp-dxm*dxm)*bz(i,j,k) - &
     & dxp*dxp*bz(i-1,j,k)) / &
     & (dxp*dxm*(dxp+dxm))
! dbzdz(i,j,k) = (dzm*dzm*bz(i,j,k+1) + &
!      & (dzp*dzp-dzm*dzm)*bz(i,j,k) - &
!      & dzp*dzp*bz(i,j,k-1)) / &
!      & (dzp*dzm*(dzp+dzm))
dbzdz(i,j,k) = (bz(i,j,k+1)*dzm)/(dzp*(dzp+dzm)) &
     & + (bz(i,j,k)*dzp)/(dzm*(dzp+dzm)) &
     & -((bz(i,j,k)*dzm)/(dzp*(dzm+dzp)) &
     & + (bz(i,j,k-1)*dzp)/(dzm*(dzp+dzm)))
endif
cz(i,j,k) = (dxm*dxm*by(i+1,j,k) + &
     & (dxp*dxp-dxm*dxm)*by(i,j,k) - &
     & dxp*dxp*by(i-1,j,k)) / &
     & (dxp*dxm*(dxp+dxm)) - &
     & (dym*dym*bx(i,j+1,k) + &
     & (dyp*dyp-dym*dym)*bx(i,j,k) - &
     & dyp*dyp*bx(i,j-1,k)) / &
     & (dyp*dym*(dyp+dym))
dbxdx(i,j,k) = (dxm*dxm*bx(i+1,j,k) + &
     & (dxp*dxp-dxm*dxm)*bx(i,j,k) - &
     & dxp*dxp*bx(i-1,j,k)) / &
     & (dxp*dxm*(dxp+dxm))
dbydy(i,j,k) = (dym*dym*by(i,j+1,k) + &
     & (dyp*dyp-dym*dym)*by(i,j,k) - &
     & dyp*dyp*by(i,j-1,k)) / &
     & (dyp*dym*(dyp+dym))
enddo
enddo
enddo
divb = dbxdx + dbydy + dbzdz
!
print *, "======================="
print *, "OUTPUT ELECTRIC CURRENT"
print *, "======================="
open(56,file=dirname//"/cx_debug",status="replace",form="unformatted")
write(56) cx
close(56)
open(57,file=dirname//"/cy_debug",status="replace",form="unformatted")
write(57) cy
close(57)
open(58,file=dirname//"/cz_debug",status="replace",form="unformatted")
write(58) cz
close(58)
open(61,file=dirname//"/divb_debug",status="replace",form="unformatted")
write(61) divb
close(61)
open(62,file=dirname//"/dbxdx_debug",status="replace",form="unformatted")
write(62) dbxdx
close(62)
open(63,file=dirname//"/dbydy_debug",status="replace",form="unformatted")
write(63) dbydy
close(63)
open(64,file=dirname//"/dbzdz_debug",status="replace",form="unformatted")
write(64) dbzdz
close(64)
!
end subroutine calc_current_divb
