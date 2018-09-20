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
integer :: ick, jck, kck
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

       call mhd__delta

       !! if(atime.le.atime_pmotion) then
       !!    call mhd__forcing
       !! end if

       if(istep.eq.1) call rkg__prog1
       if(istep.eq.2) call rkg__prog2
       if(istep.eq.3) call rkg__prog3
       if(istep.eq.4) call rkg__prog4

     end do
!
!--- debug for check of NaN ---
!      if(hanteiNaN.eq.0) then
!      do ick = 0, NX
!      do jck = 0, NY
!      do kck = 0, NZ-1
!	if(vx(kck,jck,ick).ne.vx(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "vx", vx(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!	if(vy(kck,jck,ick).ne.vy(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "vy", vy(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!	if(vz(kck,jck,ick).ne.vz(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "vz", vz(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!	if(bx(kck,jck,ick).ne.bx(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "bx", bx(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!	if(by(kck,jck,ick).ne.by(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "by", by(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!	if(bz(kck,jck,ick).ne.bz(kck,jck,ick)) then
!		write(6,'(i3,1x,i3,1x,i3,1x,A3,1x,f10.5,1x,i2,1x,i2,1x,i4)') &
!&kck, jck, ick, "bz", bz(kck,jck,ick), index_x, index_y, nloop
!		hanteiNaN=1
!	endif
!---
!      end do
!      end do
!      end do
!end if
!------------------
if(mod(nloop,nloop_output) == 0) call out__disk
     if(mod(nloop,5) == 0) then
             call pset__dt
             call pset__resistivity
     end if
     if(mod(nloop,100) == 0) then
        if(index_x == 0 .and.index_y==0) write(*,*) "Now loop count:", nloop, ", TIME: ", atime 
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

















