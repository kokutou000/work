! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
program main
! -----------------------------------------------------------
!                      PROGRAM MAIN
! -----------------------------------------------------------
! --- test for git 
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
!------------------
if(mod(nloop,nloop_output) == 0) call out__disk
     if(mod(nloop,5) == 0) then
             call pset__dt
             call pset__resistivity
     end if
!--- comment out to calculate fast
     ! if(mod(nloop,100) == 0) then
     !    if(index_x == 0 .and.index_y==0) write(*,*) "Now loop count:", nloop, ", TIME: ", atime 
     ! end if
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

















