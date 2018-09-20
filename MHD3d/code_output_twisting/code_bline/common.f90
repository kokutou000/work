include "paramete.f90"
!
      common /comcq/ &
       c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
      real(8) :: c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
!
      common /cmbfld/ bx(0:nxl,0:nyl,0:nzl), &
                      by(0:nxl,0:nyl,0:nzl), &
                      bz(0:nxl,0:nyl,0:nzl)
      real(4) bx, by, bz
!
      common /rotb/ rotbx(0:nxl,0:nyl,0:nzl), &
		     rotby(0:nxl,0:nyl,0:nzl), &
		     rotbz(0:nxl,0:nyl,0:nzl)
      real(4) rotbx, rotby, rotbz
!========== rotbx,rotby,rotbz -> rot B
      common /yzplane/ xx0, yx0, zx0
      real(8) xx0, yx0, zx0
!
      common /alpha_xy/ alpxy(0:nxl,0:nyl), &
                        alpxy1(0:nxl,0:nyl), &
                        alpxy2(0:nxl,0:nyl)
      real(4) alpxy, alpxy1, alpxy2
!
      common /twist_xy/ tw_xy(0:nxl,0:nyl)
      real(4) tw_xy
!
      common /combl/ width,rbl,gbl,bbl,wbl
      real*4 width,rbl,gbl,bbl,wbl
!
      common /cmxyz/ xl,yl,zl,                      &
                     xc(0:nxl),yc(0:nyl),zc(0:nzl), &
                     xc0(0:nx),yc0(0:ny),zc0(0:nz)
      real(8) xl,yl,zl,xc,yc,zc,xc0,yc0,zc0
!
       common /cometch/ in_dir, out_dir
      character*12,dimension(10) :: in_dir
      character*2 :: out_dir
!
      namelist /nlist00/ xl, yl, zl




