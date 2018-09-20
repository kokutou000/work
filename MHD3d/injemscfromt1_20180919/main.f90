! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
program main
! -----------------------------------------------------------
!                      PROGRAM MAIN
! -----------------------------------------------------------
use common
use job
use mpi
use mpiut
use pset
use iset
use rkg
use mhd
use out
implicit none
integer :: istep, iset_err
! ==================== set 0 to 3d arrays =================
call pset__zero_set

! ==================== read NAMELIST ======================
call pset__namelist

! ==================== MPI init  ==========================
call mpiut__init
cmyrank = '.'//chari4(myrank)

! ==================== set parameters =====================
call pset__init ! set the parameters

! ========== debug ==========
!!$     write(*,*) '::main: nproc = ',nproc
!!$     write(*,*) '::main: myrank, rank_left, rank_right, rank_up, rank_bottom, index_x, index_y '
!!$     write(*,*)  myrank,rank_left,rank_right,rank_up, rank_down, index_x, index_y
! ===========================

call iset__initial(iset_err) ! set the initial condition 
!! call iset__model  !!!! uniform magnetic field for debug
if(iset_err /= 0) go to 900

call mhd__sub

call rkg__set_coefficient 
call pset__dt
call pset__resistivity

if(start_from_initial) call out__disk

! ===================== MAIN LOOP ==========================

do while(job__is_fine())

     nloop = nloop + 1
     atime = atime + dtstep


     do istep = 1, 4

       call mhd__sub
!call ckNaNsub
       call mhd__delta
!call ckNaNdelta
       !! if(atime.le.atime_pmotion) then
       !!    call mhd__forcing
       !! end if

       if(istep.eq.1) call rkg__prog1
       if(istep.eq.2) call rkg__prog2
       if(istep.eq.3) call rkg__prog3
       if(istep.eq.4) call rkg__prog4

     end do
!
! check NaN for debug
!
!call check_NaN
!------------------
if(mod(nloop,nloop_output) == 0) call out__disk
!--debug
!--
!call debug_out
!---
!
     if(mod(nloop,5) == 0) then
             call pset__dt
             call pset__resistivity
             if((index_x.eq.0).and.(index_y.eq.0)) then
                write(*,*) "Now loop count:", nloop, ", TIME: ", atime 
             end if
     end if
!#############
!    write(*,*) '#DEBUG ',nproc
!    write(*,*) 'min max of ro(0)',minval(ro(0,:,:)),maxval(ro(0,:,:))
!    write(*,*) 'min max of pr(0)',minval(pr(0,:,:)),maxval(pr(0,:,:))
!    write(*,*) 'min max of vz(0)',minval(vz(0,:,:)),maxval(vz(0,:,:))
!#############

end do

! =================== out restart file =======================
call out__restart
! =================== close files ============================
900 call job__finalize

stop
end program main
!
!
!
!
!--- debug for check of NaN ---
subroutine check_NaN
use common
use mpiut
implicit none
integer :: ick, jck, kck
integer :: hanteiNaN = 0

  do ick = -1, NX
  do jck = -1, NY
  do kck = 0, NZ
     if(vx(kck,jck,ick).ne.vx(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN vx !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, vx(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(vy(kck,jck,ick).ne.vy(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN vy !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, vy(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(vz(kck,jck,ick).ne.vz(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN vz !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, vz(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(bx(kck,jck,ick).ne.bx(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN bx !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, bx(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(by(kck,jck,ick).ne.by(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN by !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, by(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(bz(kck,jck,ick).ne.bz(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN bz !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, bz(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(ro(kck,jck,ick).ne.ro(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN ro !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, ro(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
  end do
  end do
  end do
!
  if(hanteiNaN .eq. 1) then
     open(99,file="NaNcheck"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & bx, by, bz, vx, vy, vz, ro
     close(99)
     open(99,file="NaNcksub"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & cx, cy, cz, ex, ey, ez
     close(99)
  end if
!
  call mpiut__barrier
!
  if(hanteiNaN .eq. 1) then
     stop
  endif

end subroutine check_NaN
!
!
!
subroutine ckNaNsub
use common
use mpiut
implicit none
integer :: ick, jck, kck
integer :: hanteiNaN = 0

  do ick = -1, NX
  do jck = -1, NY
  do kck = 0, NZ
     if(cx(kck,jck,ick).ne.cx(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN cx !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, cx(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(cy(kck,jck,ick).ne.cy(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN cy !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, cy(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(cz(kck,jck,ick).ne.cz(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN cz !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, cz(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(ex(kck,jck,ick).ne.ex(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN ex !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, ex(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(ey(kck,jck,ick).ne.ey(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN ey !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, ey(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(ez(kck,jck,ick).ne.ez(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN ez !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, ez(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
  end do
  end do
  end do
!
  if(hanteiNaN .eq. 1) then
     open(99,file="NaNcheck"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & bx, by, bz, vx, vy, vz, ro
     close(99)
     open(99,file="NaNcksub"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & cx, cy, cz, ex, ey, ez
     close(99)
  end if
!
  call mpiut__barrier
!
  if(hanteiNaN .eq. 1) then
     stop
  endif

end subroutine ckNaNsub
!
!
!
subroutine ckNaNdelta
use common
use mpiut
implicit none
integer :: ick, jck, kck
integer :: hanteiNaN = 0

  do ick = -1, NX
  do jck = -1, NY
  do kck = 0, NZ
     if(dvx(kck,jck,ick).ne.dvx(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dvx !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dvx(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dvy(kck,jck,ick).ne.dvy(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dvy !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dvy(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dvz(kck,jck,ick).ne.dvz(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dvz !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dvz(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dbx(kck,jck,ick).ne.dbx(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dbx !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dbx(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dby(kck,jck,ick).ne.dby(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dby !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dby(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dbz(kck,jck,ick).ne.dbz(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dbz !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dbz(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
     if(dro(kck,jck,ick).ne.dro(kck,jck,ick)) then
        write(6, *) "===== !!! find NaN dro !!!  ====="
        write(6,'(i3,1x,i3,1x,i3,1x,f10.5,1x,i2,1x,i2,1x,i6)') &
             &kck, jck, ick, dro(kck,jck,ick), index_x, index_y, nloop
        hanteiNaN=1
     endif
  end do
  end do
  end do
!
  if(hanteiNaN .eq. 1) then
     open(99,file="NaNcheck"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & bx, by, bz, vx, vy, vz, ro
     close(99)
     open(99,file="NaNcksub"//cmyrank,form="unformatted",status="replace")
     write(99) iwrite, nloop, atime, dtstep, &
          & cx, cy, cz, ex, ey, ez
     close(99)
  end if
!
  call mpiut__barrier
!
  if(hanteiNaN .eq. 1) then
     stop
  endif

end subroutine ckNaNdelta
!
!
subroutine debug_out
use common
use mpiut

if(nloop.eq.14160) then
   open(98,file="DEBUG160"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14161) then
   open(98,file="DEBUG161"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14162) then
   open(98,file="DEBUG162"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14163) then
   open(98,file="DEBUG163"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14164) then
   open(98,file="DEBUG164"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14165) then
   open(98,file="DEBUG165"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14166) then
   open(98,file="DEBUG166"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14167) then
   open(98,file="DEBUG167"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14168) then
   open(98,file="DEBUG168"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
!
if(nloop.eq.14169) then
   open(98,file="DEBUG169"//cmyrank,form="unformatted")
   write(98) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
   close(98)
   call mpiut__barrier
endif
end subroutine debug_out
!
