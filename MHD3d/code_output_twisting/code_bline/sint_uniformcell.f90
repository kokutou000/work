subroutine sint(xi,yi,zi,sx,sy,sz,dcell,iout)
implicit none
      include "common.f90"
      real(8),intent(in) :: xi,yi,zi
      real(8),intent(out) :: sx,sy,sz,dcell
      integer,intent(inout) :: iout
!
      integer :: ii, ix0, ix1, iy0, iy1, iz0, iz1
      real(8) :: x, y, z, dx, dy, dz
      real(8) :: w000, w001, w010, w100, w011, w101, w110, w111
      real(8) :: wx0, wy0, wz0, wx1, wy1, wz1
      real(8) :: celldx, celldy, celldz
!...............................
      iout = 0
      celldx = xc(1)-xc(0)
      celldy = yc(1)-yc(0)
      celldz = zc(1)-zc(0)
!.............................. check the location inside the domain
      if(xi.lt.xc(0)) then
         iout = -100
         return
      else if(xi.gt.xc(nxl)) then
         iout = 100
         return
      end if
!
      if(yi.lt.yc(0)) then
         iout = -101
         return
      else if(yi.gt.yc(nyl)) then
         iout = 101
         return
      end if
!
      if(zi.lt.zc(0)) then
         iout = -102
         return
      else if(zi.gt.zc(nzl)) then
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
         ix0 = (x-xc(0))/celldx
         if(ix0.lt.0.or.ix0.gt.nxl) then
            write(*,*) 'out of range (x)',x,xc(0),xc(nxl),ix0,nxl 
         end if
         ix1 = ix0+1
         wx1= (x-xc(ix0))/celldx
         wx0= (xc(ix1)-x)/celldx
! [Y]
         iy0 = (y-yc(0))/celldy
         if(iy0.lt.0.or.iy0.gt.nyl) then
            write(*,*) 'out of range (y)',y
         end if
         iy1 = iy0+1
         wy1= (y-yc(iy0))/celldy
         wy0= (yc(iy1)-y)/celldy
! [Z]
         iz0 = (z-zc(0))/celldz
         if(iz0.lt.0.or.iz0.gt.nzl) then
            write(*,*) 'out of range (z)',z
         end if
         iz1 = iz0+1
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
      return
end subroutine sint



