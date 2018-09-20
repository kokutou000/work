 subroutine preprg
 implicit none
     include "common.f90"
      integer,parameter :: mbl=500000
!      integer,parameter :: mbl=100000
      common /blxy/ blx(mbl),bly(mbl),blz(mbl)
      real(8) :: blx, bly, blz
      real(8) :: blx0,bly0,blz0 
      real(8) :: qbx, qby, qbz, dbx, dby, dbz
      real(8) :: sbx, sby, sbz, dcell, bb, dbl
      real(8) :: jsbx, jsby, jsbz
      real(8),dimension(mbl) :: x,y,z
      integer :: iout, idir, ipcell, isol, nop, l0, l1
      integer :: icell=1,isci=2
      integer :: iperix=0, iperiy=0
!          >>> iperio = 1, periodic line plotted
      integer :: i, j, k
      integer :: iptw, jptw
!
      real(8) :: tw=0d0
      real(8) :: tw2 = 0d0
!      real(8) :: tw3 = 0d0
      integer :: rightleft
!	>>> rightleft is checker whether idir directs right or left
!		right(positive) -->> 1, left(negative) -->> -1
      integer :: bzidir
!       >>> bzidir is check that we should output footpoints bz > 0 or bz < 0

!
      c1 = 0.5d0
      c2 = 1.0d0 - sqrt(0.5d0)
      c3 = 1.0d0 + sqrt(0.5d0)
      c4 = 1.0d0 / 6.0d0
      cq1 = -2.0d0
      cq2 =  1.0d0 - 3.0d0 * c2
      cs2 =  2.0d0 * c2
      cq3 =  1.0d0 - 3.0d0 * c3
      cs3 =  2.0d0 * c3
!
      return
!
! ================
      entry calbl(blx0,bly0,blz0,isol,idir,iptw,jptw)
!                          isol = 1   --> eular
!                          otherwise  --> RKG
! ================
!
         if(abs(idir).ne.1) then
            write(6,*) '>>WARNING calbl: idir', idir
            return
         end if
         ipcell = idir*icell
!=== idir is "1" or "-1" -> ipcell = "1" or "-1"
!
         i = 1
         xx0 = 100
!
         blx(i) = blx0
         bly(i) = bly0
         blz(i) = blz0
!
!=== initialize rightleft=0, which means on y-z plane
!	rightleft=0
 1       continue
!
      if(isol.eq.1) then
! ++++++++++++++ Eular method
!
         call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,dcell,iout,jsbx,jsby,jsbz)
!
        if(i==1) then
                if(sbz>=0) bzidir=1
                if(sbz<0) bzidir=-1
        endif
!
         i = i+1
         if(iout.ne.0.or.i.gt.mbl) then
!            write(6,*) ' Warning: iout ',iout, i, blx(i),bly(i),blz(i)
            go to 10
         end if
!
         bb = sqrt(sbx**2+sby**2+sbz**2)
!	bb = sbx**2+sby**2+sbz**2
!         sbx = sbx/bb
!         sby = sby/bb
!         sbz = sbz/bb
!	if(i==2) rightleft = nint(sbx/abs(sbx)*(idir/abs(idir)))
         dbl = dcell/ipcell
!=== dcell is min(celldx,celldy,celldz)
!=== so, dbl is minimum celld?(xyz) over 1 or -1
!
!write(*,*) "debug1", i, mbl
         if(i.gt.mbl) go to 10
!write(*,*) "debug2"
!
         qbx = dbl*sbx/bb
         qby = dbl*sby/bb
         qbz = dbl*sbz/bb
!
         blx(i) = blx(i-1) + qbx
         bly(i) = bly(i-1) + qby
         blz(i) = blz(i-1) + qbz
!!
	 tw2 = tw2 + (jsbx*qbx + jsby*qby + jsbz*qbz)/bb
         tw = tw + sqrt(qbx**2+qby**2+qbz**2)*(jsbx*sbx + jsby*sby + jsbz*sbz)/(bb*bb)
!         tw3 = tw3 + sqrt(jsbx*qbx + jsby*qby + jsbz*qbz)/bb
      else
! +++++++++++++++ RKG
! ++++++++++++++ step 1
!
         call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,dcell,iout,jsbx,jsby,jsbz)

! ---- normalization ----
         bb = sqrt(sbx**2+sby**2+sbz**2)
         sbx = sbx/bb
         sby = sby/bb
         sbz = sbz/bb
! -----------------------

         dbl = dcell/ipcell
!         if(blz(i).lt.0.5) dbl=idir*xl/200
!
         i = i+1
!
         if(iout.ne.0.or.i.gt.mbl) go to 10
!
         qbx = dbl*sbx
         qby = dbl*sby 
         qbz = dbl*sbz
!
         blx(i) = blx(i-1) + c1*qbx
         bly(i) = bly(i-1) + c1*qby
         blz(i) = blz(i-1) + c1*qbz
!
! ++++++++++++++ step 2
!
         call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,dcell,iout,jsbx,jsby,jsbz)

! ---- normalization ----
         bb = sqrt(sbx**2+sby**2+sbz**2)
         sbx = sbx/bb
         sby = sby/bb
         sbz = sbz/bb
! -----------------------

         dbl = dcell/ipcell
!
         if(iout.ne.0) go to 10
!
         dbx = dbl*sbx
         dby = dbl*sby
         dbz = dbl*sbz
!
         blx(i) = blx(i) + c2 * (dbx -qbx)
         bly(i) = bly(i) + c2 * (dby -qby)
         blz(i) = blz(i) + c2 * (dbz -qbz)
!
         qbx = cq2 * qbx + cs2 * dbx
         qby = cq2 * qby + cs2 * dby
         qbz = cq2 * qbz + cs2 * dbz
!
!
! ++++++++++++++ step 3
!
         call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,dcell,iout,jsbx,jsby,jsbz)

! ---- normalization ----
         bb = sqrt(sbx**2+sby**2+sbz**2)
         sbx = sbx/bb
         sby = sby/bb
         sbz = sbz/bb
! -----------------------

         dbl = dcell/ipcell
!
         if(iout.ne.0) go to 10
!
         dbx = dbl*sbx
         dby = dbl*sby
         dbz = dbl*sbz
!
         blx(i) = blx(i) + c3 * (dbx -qbx)
         bly(i) = bly(i) + c3 * (dby -qby)
         blz(i) = blz(i) + c3 * (dbz -qbz)         
!
         qbx = cq3 * qbx + cs3 * dbx
         qby = cq3 * qby + cs3 * dby
         qbz = cq3 * qbz + cs3 * dbz
!
!
! ++++++++++++++ step 4
!
         call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,dcell,iout,jsbx,jsby,jsbz)

! ---- normalization ----
         bb = sqrt(sbx**2+sby**2+sbz**2)
         sbx = sbx/bb
         sby = sby/bb
         sbz = sbz/bb
! -----------------------

         dbl = dcell/ipcell
!
         if(iout.ne.0) go to 10
!
         dbx = dbl*sbx
         dby = dbl*sby
         dbz = dbl*sbz
!
         blx(i) = blx(i) + c4 * (dbx - 2.*qbx)
         bly(i) = bly(i) + c4 * (dby - 2.*qby)
         blz(i) = blz(i) + c4 * (dbz - 2.*qbz)
!
! ++++++++++++++++++++++ RKG ends
      end if
!
        go to 1
 10     nop = i - 1
!	if(idir==1.and.mod(nop,200)==0) print *, "output blplot, nop = ", nop, tw
	tw = tw / (4d0 * pi)
	if(idir==-1) tw2 = tw2 * -1d0
        tw_xy(iptw,jptw) = tw + tw_xy(iptw,jptw)
!        call blplot(nop,blx(1),bly(1),blz(1),idir,tw,blx(nop),bly(nop),blz(nop))
!--- debug check
if(iptw.eq.NX/2 .and. jptw.eq.9*NY/16) then
write(*,*) "iptw,jptw,nop,idir", iptw,jptw,nop, idir
write(*,*) "tw,tw2,twxy", tw,tw2,tw_xy(iptw,jptw)
write(*,*) "bz(z=0)", bz(iptw,jptw,0)
write(*,*) "iout", iout
endif
!----- output foot point data
if(idir.eq.1) then
   write(70,*) blx(1), bly(1), blz(1)
   write(70,*) blx(nop), bly(nop), blz(nop)
   write(70,*) nop
else if(idir.eq.-1) then
   write(71,*) blx(1), bly(1), blz(1)
   write(71,*) blx(nop), bly(nop), blz(nop)
   write(71,*) nop
end if
!----- output connection from footpoints
!----- to yz plane
if(bzidir.eq.idir) then
   if(xx0.ne.100) then
      write(60,*) blx(1), bly(1), blz(1)
      write(60,*) xx0, yx0, zx0
   endif
endif
! ------ Initializing twist number = 0d0
	tw = 0d0
        tw2 = 0d0
      return
end subroutine preprg
