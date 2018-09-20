! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
module iset
! -----------------------------------------------------------
!                    MODULE iset
! -----------------------------------------------------------
use common
implicit none

! === 2D vector field in 2D space for initial equilibrium ===
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_bx
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_by
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_bz
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_vx
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_vy
real(DP), dimension(-1:NY2D,0:NZ2D) :: vec2_vz

! === eigenfunction for the perturbation
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_bx
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_by
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_bz
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_vx
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_vy
complex(DP), dimension(-1:NY2D,0:NZ2D) :: cvec2_vz
complex(DP) :: expkx ! exp(i pi2 m x/lx)

! === for perturbation phase
real(DP) :: rand ! random number
integer :: irand ! source integer for random number

integer :: mmode, mode, i,j,k,iw
real(8) :: bz0max

contains
! -----------------------------------------------------------
subroutine iset__initial(iset_err)
integer, INTENT(OUT) :: iset_err
!!!! iset_err = 0 : fine, 100 : grid check fail


! ------------------------------ reset main array fields

    vx(:,:,:) = 0.0
    vy(:,:,:) = 0.0
    vz(:,:,:) = 0.0
    bx(:,:,:) = 0.0
    by(:,:,:) = 0.0
    bz(:,:,:) = 0.0
    ro(:,:,:) = 0.0
    pr(:,:,:) = 0.0
    dvx(:,:,:) = 0.0 
    dvy(:,:,:) = 0.0 
    dvz(:,:,:) = 0.0
    dbx(:,:,:) = 0.0 
    dby(:,:,:) = 0.0 
    dbz(:,:,:) = 0.0
    dro(:,:,:) = 0.0 
    dpr(:,:,:) = 0.0
    qvx(:,:,:) = 0.0 
    qvy(:,:,:) = 0.0 
    qvz(:,:,:) = 0.0
    qbx(:,:,:) = 0.0 
    qby(:,:,:) = 0.0 
    qbz(:,:,:) = 0.0
    qro(:,:,:) = 0.0 
    qpr(:,:,:) = 0.0
    ex (:,:,:) = 0.0 
    ey (:,:,:) = 0.0 
    ez (:,:,:) = 0.0
    cx (:,:,:) = 0.0 
    cy (:,:,:) = 0.0 
    cz (:,:,:) = 0.0
    rvx(:,:,:) = 0.0 
    rvy(:,:,:) = 0.0
    rvz(:,:,:) = 0.0
    qx (:,:,:) = 0.0 
    qy (:,:,:) = 0.0
    qz (:,:,:) = 0.0
    tmp(:,:,:) = 0.0
    phi(:,:,:) = 0.0
    dphi(:,:,:) = 0.0
    qphi(:,:,:) = 0.0

  if(start_from_initial) then
   !!-----------------------------------------------------
   !! The initial B is given by FILE_NLFF and 
   !! B&V on the bottom boundary are given by FILE_2D_XY
   !    call iset__equilibrium_nlff_vtvc
   !!-----------------------------------------------------
   !! The initial B&V are given by FILE_2D_FIELD that is 
   !! the (Y,Z) data.
   !!#### check grid number ####
   !   if(NY*nproc_y.ne.NY2D) then
   !     write(*,*) '##ERR:iset__initial: NY*nproc_y .ne. NY2D',NY,nproc_y,NY2D
   !    iset_err = 100
   !    return
   !  end if
   !!###########################
   !   call iset__equilibrium_file2d
   !!-----------------------------------------------------
   !! The initial B is given by the LFFF
   !call iset__equilibrium_lfff
   !!-----------------------------------------------------
   !! add instability eigenfunction from FILE_EIGENMODE 
   !! as perturbation
   !    if(add_perturbation) then
   !     call iset__perturbation  
   !    end if
   !!-----------------------------------------------------
   !! The initial condition is given by potential field based on spheromak sunspot  
   call iset__init_file3d
   !
   !!-----------------------------------------------------

      nloop_end = nloop_incmax
  
    else
 
      call iset__restart
      nloop_end = nloop + nloop_incmax
!------------------------------------------------------
!----- get initvx, vy arrays
 open(FILE_3D_INIT,file=trim(cfile_3d_init), &
!!      convert='big_endian', &  ! for Absoft_f90
      form='unformatted')
    read(FILE_3D_INIT) initbx3d, initby3d, initbz3d, &
                     & initvx2d, initvy2d, initvz2d
 close(FILE_3D_INIT)
!----- reset density too
! ro = sqrt(bx**2+by**2+bz**2)
!   do i = 0, NX
!   do j = 0, NY
!   do k = 0, NZ
!   if(ro(k,j,i).le.ro_min) ro(k,j,i) = ro_min
!   enddo
!   enddo
!   enddo
!----- set visc3d again
 visc3d(:,:,:) = visc
 do i = -1, NX
 do j = -1, NY
 do k =  0, NZ
    visc3d(k,j,i)=visc*(dx(i)*dy(j)*dz(k)/(0.02023*0.01723**2))**(1.0/3.0)
 end do
 end do
 end do
!----------------------------------------------------
  end if

  iset_err = 0

  return

end subroutine iset__initial
! -----------------------------------------------------------
subroutine iset__equilibrium_file2d
use common
integer :: j2d
!!
 open(FILE_2D_FIELD,file=trim(cfile_2d_field), &
!!      convert='big_endian', &  ! for Absoft_f90
      form='unformatted')
    read(FILE_2D_FIELD) vec2_bx, vec2_by, vec2_bz, &
                        vec2_vx, vec2_vy, vec2_vz
 close(FILE_2D_FIELD)

! === copy to 3d vector
  do i = -1, NX
  do j = -1, NY
     j2d = index_y*NY + j
  do k =  0, NZ
     bx(k,j,i) = vec2_bx(j2d,k)
     by(k,j,i) = vec2_by(j2d,k)
     bz(k,j,i) = vec2_bz(j2d,k)

     vx(k,j,i) = vec2_vx(j2d,k)
     vy(k,j,i) = vec2_vy(j2d,k)
     vz(k,j,i) = vec2_vz(j2d,k)
  end do
  end do
  end do

! === reset the velocity in the whole volume ===
  if(reset_velocity) then
     vx(:,:,:) = 0.0d0
     vy(:,:,:) = 0.0d0
     vz(:,:,:) = 0.0d0
  end if

! === reset density and pressure ===
  ro(:,:,:) = ro_init
  pr(:,:,:) = pr_init

  return

end subroutine iset__equilibrium_file2d
! -----------------------------------------------------------
! -----------------------------------------------------------
subroutine iset__init_file3d
use common
integer :: j3d, i3d
!!
 open(FILE_3D_INIT,file=trim(cfile_3d_init), &
!!      convert='big_endian', &  ! for Absoft_f90
      form='unformatted')
    read(FILE_3D_INIT) initbx3d, initby3d, initbz3d, &
                     & initvx2d, initvy2d, initvz2d
 close(FILE_3D_INIT)
! === decrease initial velocity
! initvx2d = initvx2d/maxval(abs(initvx2d))
! initvy2d = initvy2d/maxval(abs(initvy2d))
! initvx2d = initvx2d*0.1d0
! initvy2d = initvy2d*0.1d0
! initbx3d = initbx3d*0.1d0
! initby3d = initby3d*0.1d0
! initbz3d = initbz3d*0.1d0
! === reset the velocity in the whole volume ===
  vx(:,:,:) = 0.0d0
  vy(:,:,:) = 0.0d0
  vz(:,:,:) = 0.0d0

 ! === copy to 3d vector
  do i = -1, NX
     i3d = index_x*NX + i
  do j = -1, NY
     j3d = index_y*NY + j
  do k =  0, NZ
    bx(k,j,i) = initbx3d(i3d,j3d,k)
    by(k,j,i) = initby3d(i3d,j3d,k)
    bz(k,j,i) = initbz3d(i3d,j3d,k)
!    if(k.eq.0) then
!       vx(k,j,i) = initvx2d(i3d,j3d)
!       vy(k,j,i) = initvy2d(i3d,j3d)
!    end if
  end do
  end do
  end do
!--- mask velocity on z=0 for some value
!   bz0max = maxval(initbz3d(:,:,0))
!   do i = -1, NX
!      i3d = index_x*NX + i
!      do j = -1, NY
!         j3d = index_y*NY + j
!         k = 0
! !
!         vx(k,j,i) = vx(k,j,i) * 0.5*(1d0-tanh(10*(abs(bz(0,j,i)/bz0max) - 0.6)))
!         vy(k,j,i) = vy(k,j,i) * 0.5*(1d0-tanh(10*(abs(bz(0,j,i)/bz0max) - 0.6)))
! !
! !        if(abs(bz(0,j,i))>bz0max*0.8) then
! !           vx(k,j,i) = 0d0
! !           vy(k,j,i) = 0d0
! !        end if
! !
!      end do
!   end do
!---
! === reset the velocity in the whole volume ===
!  vx(:,:,:) = 0.0d0
!  vy(:,:,:) = 0.0d0
!  vz(:,:,:) = 0.0d0

! === reset density and pressure ===
  ro(:,:,:) = ro_init
  pr(:,:,:) = pr_init
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! === for constant alfven wave ---20161223
 ro(:,:,:) = sqrt(bx(:,:,:)**2+by(:,:,:)**2+bz(:,:,:)**2)
 do k = 0, NZ
 do j = 0, NY
 do i = 0, NX
    ! if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. &
    !      & (abs(zc(k)).le.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i)*(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) + &
    !         & sqrt(ro(k,j,i))*(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)
    ! else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
    !      & (abs(zc(k)).le.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i)*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) + &
    !         & sqrt(ro(k,j,i))*(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)
    ! else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. & 
    !      & (abs(zc(k)).ge.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i)*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) + &
    !         & sqrt(ro(k,j,i))*(zc(k)-0.9*zl)/(0.1*zl)
    ! else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
    !      & (abs(zc(k)).le.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i)* &
    !         & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) + &
    !         & sqrt(ro(k,j,i))* (1.0-&
    !         & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)))
    ! else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. & 
    !      & (abs(zc(k)).ge.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i)* &
    !         & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) + &
    !         & sqrt(ro(k,j,i))* (1.0-&
    !         & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
    ! else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
    !      & (abs(zc(k)).ge.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i) &
    !         &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) &
    !         &+sqrt(ro(k,j,i)) * (1.0 &
    !         &-(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
    ! else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
    !      & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
    !      & (abs(zc(k)).ge.0.9d0*zl)) then
    !    ro(k,j,i) = ro(k,j,i) &
    !         &*(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) &
    !         &+sqrt(ro(k,j,i)) * (1.0 &
    !         &-(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
    !         &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
    !         &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
    ! end if
 if(ro(k,j,i).le.ro_min) ro(k,j,i) = ro_min
 enddo
 enddo
 enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!----- set large viscosity near boundary to imitate open condition
 visc3d(:,:,:) = visc
 do i = -1, NX
 do j = -1, NY
 do k =  0, NZ
    visc3d(k,j,i)=visc*(dx(i)*dy(j)*dz(k)/(0.02023*0.01723**2))**(1.0/3.0)
 end do
 end do
 end do
!---------------------------------------

  return

end subroutine iset__init_file3d
! -----------------------------------------------------------
subroutine iset__equilibrium_lfff
use common
real(DP) :: bx_lfff, by_lfff, bz_lfff
real(DP) :: phase
!!
! === LFFF ===
!
  do k =  0, NZ
  do j = -1, NY
     phase   = lfff_k*(yc(j)-lfff_y0)

     bx_lfff = lfff_b0*lfff_alpha/lfff_k     &
              *cos(phase)*exp(-lfff_kk*zc(k))
     by_lfff =-lfff_b0*lfff_kk/lfff_k        &
              *cos(phase)*exp(-lfff_kk*zc(k))
     bz_lfff = lfff_b0                       &
              *sin(phase)*exp(-lfff_kk*zc(k))

  do i = -1, NX
     bx(k,j,i) = bx_lfff
     by(k,j,i) = by_lfff
     bz(k,j,i) = bz_lfff

     vx(k,j,i) = 0.0
     vy(k,j,i) = 0.0
     vz(k,j,i) = 0.0
  end do

  end do
  end do

! === set density and pressure ===

  ro(:,:,:) = ro_init
  pr(:,:,:) = pr_init

  return

end subroutine iset__equilibrium_lfff
! -----------------------------------------------------------
subroutine iset__equilibrium_nlff_vtvc
! set the NLFF to B[x,y,z], and 
!     the converging and twisting flow to V[x,y,z]
!:::::::::::::::::::::::::::::::::::::::::::20071102(KusanoK)
use common
real(DP), dimension(0:NX_NLFF,0:NY_NLFF,0:NZ) :: bx_nlff, by_nlff, bz_nlff
real(DP), dimension(0:NX_NLFF,0:NY_NLFF) :: bx2d, by2d, bz2d
real(DP), dimension(0:NX_NLFF,0:NY_NLFF) :: vx2d, vy2d, vz2d
real(DP), dimension(0:NX_NLFF,0:NY_NLFF) :: vx2d_c, vy2d_c, vz2d_c
real(DP), dimension(0:NX_NLFF,0:NY_NLFF) :: vx2d_t, vy2d_t, vz2d_t
real(DP) :: vel_0, vel_max
integer :: i_nlff, j_nlff

 open(FILE_NLFF,file=trim(cfile_nlff), &
!!       convert='big_endian', &  ! for Absoft_f90
       form='unformatted')
     read(FILE_NLFF) bx_nlff, by_nlff, bz_nlff
 close(FILE_NLFF)

 open(FILE_2D_XY,file=trim(cfile_nlff_2d), &
      form='unformatted')
    read(FILE_2D_XY) bx2d, by2d, bz2d
    read(FILE_2D_XY) vx2d_t, vy2d_t, vz2d_t
    read(FILE_2D_XY) vx2d_c, vy2d_c, vz2d_c
 close(FILE_2D_XY)

! make the forcing field

  vx2d(:,:) = av_c*vx2d_c(:,:)+av_t*vx2d_t(:,:)
  vy2d(:,:) = av_c*vy2d_c(:,:)+av_t*vy2d_t(:,:)
  vz2d(:,:) = av_c*vz2d_c(:,:)+av_t*vz2d_t(:,:)

! normalization of velocity

  vel_max = sqrt(maxval(vx2d(:,:)**2+vy2d(:,:)**2+vz2d(:,:)**2))
  vel_0 = AMPP/vel_max

  vx2d(:,:) = vx2d(:,:)*vel_0
  vy2d(:,:) = vy2d(:,:)*vel_0
  vz2d(:,:) = vz2d(:,:)*vel_0

! === reset the velocity in the whole volume ===

     vx(:,:,:) = 0.0d0
     vy(:,:,:) = 0.0d0
     vz(:,:,:) = 0.0d0

! === copy to 3d vector
  do i = -1, NX
     i_nlff = index_x*NX + i
     if(i_nlff.eq.-1) i_nlff=NX_NLFF-1
     if(i_nlff.gt.NX_NLFF) write(FILE_SYSOUT,*) '**WARNING NX_NLFF<i_nlff',NX_NLFF,i_nlff

  do j = -1, NY
     j_nlff = index_y*NY + j
     if(j_nlff.eq.-1) j_nlff=NY_NLFF-1
     if(j_nlff.gt.NY_NLFF) write(FILE_SYSOUT,*) '**WARNING NY_NLFF<j_nlff',NY_NLFF,j_nlff

   do k =  0, NZ
     bx(k,j,i) = bx_nlff(i_nlff,j_nlff,k)
     by(k,j,i) = by_nlff(i_nlff,j_nlff,k)
     bz(k,j,i) = bz_nlff(i_nlff,j_nlff,k)
   end do

! photospheric velocity
!  vx(0,j,i) = vx2d(i_nlff,j_nlff)
!  vy(0,j,i) = vy2d(i_nlff,j_nlff)
!  vz(0,j,i) = vz2d(i_nlff,j_nlff)

! forcing field  (introduced in NL3DpwE_mpi_02g_nlff: 20071122)
   force_x(j,i) = vx2d(i_nlff,j_nlff)
   force_y(j,i) = vy2d(i_nlff,j_nlff)
   force_z(j,i) = vz2d(i_nlff,j_nlff)

  end do
  end do

! === reset density and pressure ===
  ro(:,:,:) = ro_init
  pr(:,:,:) = pr_init

  return

end subroutine iset__equilibrium_nlff_vtvc
! -----------------------------------------------------------
subroutine iset__reset_pmotion

vx(0,:,:) = 0.0
vy(0,:,:) = 0.0
vz(0,:,:) = 0.0

return
end subroutine iset__reset_pmotion
! -----------------------------------------------------------
subroutine iset__model
use common

  if(start_from_initial) then

     bx(:,:,:) = 1.0d0
     by(:,:,:) = 0.0d0
     bz(:,:,:) = 0.0d0
     vx(:,:,:) = 0.0d0
     vy(:,:,:) = 0.0d0
     vz(:,:,:) = 0.0d0
     ro(:,:,:) = 1.0d0
     pr(:,:,:) = 0.0d0

      nloop_end = nloop_incmax
  
    else
 
      call iset__restart
      nloop_end = nloop + nloop_incmax

  end if

  return

end subroutine iset__model
! -----------------------------------------------------------
subroutine iset__perturbation
use common
integer :: j2d
! ---------------------------------
! read eigenfunction
! ---------------------------------
open(FILE_EIGENMODE,file=trim(cfile_eigenmode), &
!!     convert='big_endian', ! for Absoft_f90
     form='unformatted')
  read(FILE_EIGENMODE,err=800) mode, cvec2_vx, cvec2_vy, cvec2_vz, &
                                     cvec2_bx, cvec2_by, cvec2_bz
close(FILE_EIGENMODE)
! ---------------------------------
  irand = 12345
 
  do mmode = mmode1, mmode2

    call urand1(1,rand,irand)

      do i = -1, NX
      do j = -1, NY
         j2d = index_y*NY + j
      do k =  0, NZ
        expkx = ampp*exp(IUNIT*PI2*(mmode*xc(i)/xl+rand))

         vx(k,j,i) = vx(k,j,i) + cvec2_vx(j2d,k)*expkx
         vy(k,j,i) = vy(k,j,i) + cvec2_vy(j2d,k)*expkx
         vz(k,j,i) = vz(k,j,i) + cvec2_vz(j2d,k)*expkx
         bx(k,j,i) = bx(k,j,i) + cvec2_bx(j2d,k)*expkx
         by(k,j,i) = by(k,j,i) + cvec2_by(j2d,k)*expkx
         bz(k,j,i) = bz(k,j,i) + cvec2_bz(j2d,k)*expkx
      end do
      end do
      end do

   end do
        
  return
800 write(FILE_SYSOUT,*) ' ### ERR: File Not Found in iset__perturbation '
    stop
  
end subroutine iset__perturbation

! -----------------------------------------------------------
subroutine iset__restart
use common

 open(FILE_RESTART,file=trim(cfile_restart)//cmyrank, &
!!      convert='big_endian', ! for Absoft_f90
      form='unformatted')
  read(FILE_RESTART,err=801) iw, nloop, atime, dtstep, &
                     bx, by, bz, vx, vy, vz, ro, pr, phi
 close(FILE_RESTART)

  iwrite = iw

  return

801 write(FILE_SYSOUT,*) ' ### ERR: File Not Found in iset__restart '
    stop

end subroutine iset__restart

end module iset












