subroutine sint(xi,yi,zi,sx,sy,sz,dcell,iout,rotbsx,rotbsy,rotbsz)
implicit none
      include "common.f90"
      real(8),intent(inout) :: xi
      real(8),intent(in) :: yi,zi
!      real(8),intent(in) :: xi,yi,zi
      real(8),intent(out) :: sx,sy,sz,dcell
      real(8),intent(out) :: rotbsx,rotbsy,rotbsz
      integer,intent(inout) :: iout
!
      integer :: ii, ix0, ix1, iy0, iy1, iz0, iz1
      real(8) :: x, y, z, dx, dy, dz
      real(8) :: w000, w001, w010, w100, w011, w101, w110, w111
      real(8) :: wx0, wy0, wz0, wx1, wy1, wz1
      real(8) :: celldx, celldy, celldz
      integer :: find_cell_index
!...............................
      iout = 0
      celldx = xc(1)-xc(0)
      celldy = yc(1)-yc(0)
      celldz = zc(1)-zc(0)
!.............................. check the location inside the domain
!---------- x dir is periodic boundary
!      if(xi.lt.xc(0)) then
!!         iout = -100
!!         return
!          xi = xi + (xc(nxl)-xc(0))
!---------- This x-boundary is periodic, then we don't have to stop calculation at xc(nxl-1)
!---------- but we should connect xc(nxl) <-> xc(0)
!      else if(xi.gt.xc(nxl)) then
!!      else if(xi.gt.xc(nxl-1)) then
!!         iout = 100
!!         return
!         xi = xi + (xc(0)- xc(nxl))
!      end if
!
!---------- in the case of x dir is fixed boundary
      if(xi.lt.xc(0)) then
         iout = -100
         return
!      else if(xi.gt.xc(nxl)) then
      else if(xi.gt.xc(nxl-1)) then
         iout = 100
         return
      end if
!-------------------------------------------------
      if(yi.lt.yc(0)) then
         iout = -101
         return
      else if(yi.gt.yc(nyl-1)) then
!      else if(yi.gt.yc(nyl)) then
         iout = 101
         return
      end if
!
      if(zi.lt.zc(0)) then
         iout = -102
         return
      else if(zi.gt.zc(nzl-1)) then
!      else if(zi.gt.zc(nzl)) then
         iout = 102
         return
      end if
!
      x = xi
      y = yi
      z = zi
!
! >>>> Find the cell in which the point is
!
! [X]
         ix0 = find_cell_index(xc,nxl,x)
         if(ix0.lt.0.or.ix0.gt.nxl) then
            write(*,*) 'out of range (x)',x,xc(0),xc(nxl),ix0,nxl 
         end if
         ix1 = ix0+1
         celldx = xc(ix1)-xc(ix0)
         wx1= (x-xc(ix0))/celldx
         wx0= (xc(ix1)-x)/celldx
! [Y]
         iy0 = find_cell_index(yc,nyl,y)
         if(iy0.lt.0.or.iy0.gt.nyl) then
            write(*,*) 'out of range (y)',y
         end if
         iy1 = iy0+1
         celldy = yc(iy1)-yc(iy0)
         wy1= (y-yc(iy0))/celldy
         wy0= (yc(iy1)-y)/celldy
! [Z]
         iz0 = find_cell_index(zc,nzl,z)
         if(iz0.lt.0.or.iz0.gt.nzl) then
            write(*,*) 'out of range (z)',z
         end if
         iz1 = iz0+1
         celldz = zc(iz1)-zc(iz0)
         wz1= (z-zc(iz0))/celldz
         wz0= (zc(iz1)-z)/celldz
!
         dcell = min(celldx,celldy,celldz)
!
! >>> interpolation (linear)
!
         w000 = wx0*wy0*wz0
         w001 = wx0*wy0*wz1
         w010 = wx0*wy1*wz0
         w011 = wx0*wy1*wz1
         w100 = wx1*wy0*wz0
         w101 = wx1*wy0*wz1
         w110 = wx1*wy1*wz0
         w111 = wx1*wy1*wz1
!
         sx = w000*bx(ix0,iy0,iz0) &
             +w001*bx(ix0,iy0,iz1) &
             +w010*bx(ix0,iy1,iz0) &
             +w011*bx(ix0,iy1,iz1) &
             +w100*bx(ix1,iy0,iz0) &
             +w101*bx(ix1,iy0,iz1) &
             +w110*bx(ix1,iy1,iz0) &
             +w111*bx(ix1,iy1,iz1)
!
         sy = w000*by(ix0,iy0,iz0) &
             +w001*by(ix0,iy0,iz1) &
             +w010*by(ix0,iy1,iz0) &
             +w011*by(ix0,iy1,iz1) &
             +w100*by(ix1,iy0,iz0) &
             +w101*by(ix1,iy0,iz1) &
             +w110*by(ix1,iy1,iz0) &
             +w111*by(ix1,iy1,iz1)
!
         sz = w000*bz(ix0,iy0,iz0) &
             +w001*bz(ix0,iy0,iz1) &
             +w010*bz(ix0,iy1,iz0) &
             +w011*bz(ix0,iy1,iz1) &
             +w100*bz(ix1,iy0,iz0) &
             +w101*bz(ix1,iy0,iz1) &
             +w110*bz(ix1,iy1,iz0) &
             +w111*bz(ix1,iy1,iz1)
!
         rotbsx = w000*rotbx(ix0,iy0,iz0) &
             +w001*rotbx(ix0,iy0,iz1) &
             +w010*rotbx(ix0,iy1,iz0) &
             +w011*rotbx(ix0,iy1,iz1) &
             +w100*rotbx(ix1,iy0,iz0) &
             +w101*rotbx(ix1,iy0,iz1) &
             +w110*rotbx(ix1,iy1,iz0) &
             +w111*rotbx(ix1,iy1,iz1)
!
         rotbsy = w000*rotby(ix0,iy0,iz0) &
             +w001*rotby(ix0,iy0,iz1) &
             +w010*rotby(ix0,iy1,iz0) &
             +w011*rotby(ix0,iy1,iz1) &
             +w100*rotby(ix1,iy0,iz0) &
             +w101*rotby(ix1,iy0,iz1) &
             +w110*rotby(ix1,iy1,iz0) &
             +w111*rotby(ix1,iy1,iz1)
!
         rotbsz = w000*rotbz(ix0,iy0,iz0) &
             +w001*rotbz(ix0,iy0,iz1) &
             +w010*rotbz(ix0,iy1,iz0) &
             +w011*rotbz(ix0,iy1,iz1) &
             +w100*rotbz(ix1,iy0,iz0) &
             +w101*rotbz(ix1,iy0,iz1) &
             +w110*rotbz(ix1,iy1,iz0) &
             +w111*rotbz(ix1,iy1,iz1)
!
         if(ix0.eq.nxl/2.and.sx.le.0) then
            xx0 = x
            yx0 = y
            zx0 = z
         else if(ix1.eq.nxl/2.and.sx.ge.0) then
            xx0 = x
            yx0 = y
            zx0 = z
         endif
!
      return
end subroutine sint
!!$integer function find_cell_index(a,n,x)
!!$real(8),intent(in) :: a(n)
!!$integer,intent(in) :: n
!!$real(8),intent(in) :: x
!!$integer :: iminl(1)
!!$iminl = minloc((a(:)-x)**2)
!!$if(a(iminl(1)).ge.x) then
!!$   find_cell_index = iminl(1)-1
!!$else
!!$   find_cell_index = iminl(1)
!!$end if
!!$end function find_cell_index
integer function find_cell_index(a,n,x)
real(8),intent(in) :: a(n+1)
integer,intent(in) :: n
real(8),intent(in) :: x
integer :: iminl(1)
iminl = minloc((a(:)-x)**2)
!if(a(iminl(1)).ge.x.and.a(1)/=x) then
!if(a(iminl(1)).ge.x) then
if(a(iminl(1)).gt.x) then
   find_cell_index = iminl(1)-1
else
   find_cell_index = iminl(1)
end if
!!!!!!! debug
! if(a(find_cell_index).gt.x.or.a(find_cell_index+1).lt.x) then
! write(*,*) 'find_cell_index ',find_cell_index,x,a(find_cell_index),a(find_cell_index+1)
! end if
!!!!!!!
find_cell_index=find_cell_index-1
end function find_cell_index


