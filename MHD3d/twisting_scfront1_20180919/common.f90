! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
module common
! -----------------------------------------------------------
!                    MODULE COMMON
! -----------------------------------------------------------
use constants
implicit none

! === define of basic strucutre ===
type vector_3d_
   real(DP), dimension(0:NZ,-1:NY,-1:NX) :: x
   real(DP), dimension(0:NZ,-1:NY,-1:NX) :: y
   real(DP), dimension(0:NZ,-1:NY,-1:NX) :: z
end type vector_3d_

! === 3D vector fields ===
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  vx,  vy,  vz
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  bx,  by,  bz
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  ro,  pr
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dvx, dvy, dvz
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dbx, dby, dbz
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dro, dpr
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: qvx, qvy, qvz
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: qbx, qby, qbz
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: qro, qpr
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  ex,  ey,  ez
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  cx,  cy,  cz
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: rvx, rvy, rvz ! ro*v[x,y,z]
real(DP), dimension(0:NZ,-1:NY,-1:NX) ::  qx,  qy,  qz ! heat flux
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: tmp ! temperature
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: phi ! variable for divb = 0
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dphi, qphi

! === 3d vector field in 3d space for initial and lateral condition ===
real(DP),dimension(-1:NX3D,-1:NY3D,0:NZ3D) :: initbx3d
real(DP),dimension(-1:NX3D,-1:NY3D,0:NZ3D) :: initby3d
real(DP),dimension(-1:NX3D,-1:NY3D,0:NZ3D) :: initbz3d

! === initial and lateral velocity on bottom ===
real(DP),dimension(-1:NX3D,-1:NY3D) :: initvx2d
real(DP),dimension(-1:NX3D,-1:NY3D) :: initvy2d
real(DP),dimension(-1:NX3D,-1:NY3D) :: initvz2d

! === forcing field  ===
real(DP), dimension(-1:NY,-1:NX) :: force_x, force_y, force_z

! === scaler field in 3D Space ===
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: res ! resistivity
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: visc3d ! viscosicy 3d

! === geometrical parameters ===
real(DP) :: xl       ! size of the x coordinate
real(DP) :: yl       ! size of the y coordinate
real(DP) :: zl       ! size of the z coordinate
real(DP) :: epsx     ! grid packing parameter for x
real(DP) :: epsy     ! grid packing parameter for y
real(DP) :: epsz     ! grid packing parameter for z

real(DP), dimension(-1:NX) :: xc, dx, ddx, dx2, ddx2
real(DP), dimension(-1:NY) :: yc, dy, ddy, dy2, ddy2
real(DP), dimension( 0:NZ) :: zc, dz, ddz, dz2, ddz2
real(DP), dimension(-1:NX,-1:1) :: d1x, d2x
real(DP), dimension(-1:NY,-1:1) :: d1y, d2y
real(DP), dimension(-1:NZ,-1:1) :: d1z, d2z
real(DP), dimension(0:NZ,-1:NY,-1:NX) :: d2xyz0
! --------------------(memo 2005.03.01)--------------------------
! laplacian f = f(k+1,,)*d2z(k,1) + f(k-1,,)*d2z(j,-1) &
!             + f(,j+1,)*d2y(j,1) + f(,j-1,)*d2y(j,-1) &
!             + f(,,i+1)*d2x(i,1) + f(,,i-1)*d2z(i,-1) &
!             + f(k,j,i)*d2xyz0(k,j,i)
! where d2xyz0 = d2z(0)+d2y(0)+d2x(0)
! ---------------------------------------------------------------
real(DP), dimension(-1:1) :: d1y_0, d2y_0
real(DP), dimension(-1:1) :: d1y_ny,d2y_ny
real(DP), dimension(-1:1) :: d1x_0
real(DP), dimension(-1:1) :: d1x_nx

! === diffusion parameters ===
real(DP) :: eta, eta1, cc0, visc, diff_ro, akappa, iohm_heat

! === physical paramter ===
real(DP) :: gamma, gamma1

! === paramter for time and timestep ===
real(DP) :: atime, dtstep, cfr
real(DP) :: atime_end=160d0   ! finish Time of calc
real(DP) :: time_ramp_init=10d0   ! Time of initial incresing state
real(DP) :: time_ramp_end=10d0   ! Time of end decresing state
real(DP) :: time_end_tw=160d0   ! Start time of end decresing state
integer :: run_number   ! run number
character(3) :: crun_number ! character for run number
integer :: nloop        ! current loop number
integer :: nloop_output ! output interval step
integer :: nloop_incmax ! maximum calculation loop in the current task
integer :: nloop_end    ! last loop number
!integer :: iwrite  !  index of the data on disk
integer :: iwrite=0  !  index of the data on disk
character(4) :: cwrite

! === parameter for divb=0
real(DP) :: ch2, cp2 ! coefficient for divb cleaning

! === parameter for check of NaN
integer :: hanteiNaN=0

! === parameter for initial condition ===
integer :: mmode1, mmode2 ! Fourier mode number (min & max) for the perturbation 
logical :: start_from_initial 
logical :: reset_velocity 
logical :: add_perturbation
real(DP) :: ampp ! amplitude perturbation
real(DP) :: av_c, av_t ! the normalized fraction of Vc and Vt
real(DP) :: atime_pmotion 
! the photospheric motion is imposed only for atime<atime_pmotion, however,
! if atime_pmotion.le.0.0 then the photospheric motion is not terminated for ever.
! === parameters for flux emerging
real(DP) :: eflux_x0, eflux_y0, eflux_z0, & ! Coordinate of torus center
            eflux_ph, & ! Horizontal angle of emerging flux
            eflux_vz, & ! Emerging speed
            eflux_t0, eflux_t1, & ! Time the emerging starts and stops.
            eflux_b0, eflux_r0, eflux_rm, eflux_a0

! === parameters for lfff ===
real(DP) :: lfff_b0, lfff_alpha, lfff_kly_pi, lfff_y0, &
            lfff_k, lfff_kk

! === parameters for density and pressure
real(DP) :: ro_init, pr_init
real(DP) :: ro_min, pr_min

! === parameter for mpi ===
integer :: nproc, nproc_x, nproc_y, MPI_COMM_CART
integer :: myrank, root=0
integer :: rank_right, rank_left ! neighbor region for x-coordinate
integer :: rank_up, rank_down    ! neighbor region for y-coordinate 
integer :: index_x, index_y ! coordinate index for each region
character(5) :: cmyrank

! === file name ===
character(100) :: cfile_sysout
character(100) :: cfile_run_number
character(100) :: cfile_coordinate_x
character(100) :: cfile_coordinate_y
character(100) :: cfile_coordinate_z
character(100) :: cfile_coordinate_z_all
character(100) :: cfile_eigenmode
character(100) :: cfile_output_list
character(100) :: cfile_time_list
character(100) :: cfile_nlff
character(100) :: cfile_nlff_2d
character(100) :: cfile_2d_field
character(100) :: cfile_3d_field
character(100) :: cfile_3d_init
character(100) :: cfile_restart
character(100) :: cfile_slice_xz

namelist /nlist00/ xl, yl, zl
namelist /nlist01/ eta, eta1, cc0, iohm_heat
namelist /nlist02/ visc, diff_ro, akappa, gamma
namelist /nlist03/ cfr
namelist /nlist04/ nloop_incmax, &
                   nloop_output
namelist /nlist05/ epsx, epsy, epsz
namelist /nlist06/ reset_velocity,     &
                   add_perturbation,   &
                   ro_init, pr_init,   &
                   ro_min,  pr_min,    &
                   mmode1, mmode2, ampp
namelist /nlist06a/ av_c, av_t, atime_pmotion
namelist /nlist07/ cfile_sysout,       &
                   cfile_run_number,   &
                   cfile_coordinate_x, &
                   cfile_coordinate_y, &
                   cfile_coordinate_z, &
                   cfile_coordinate_z_all, &
                   cfile_eigenmode,    &
                   cfile_output_list,  &
                   cfile_time_list,    &
                   cfile_nlff,         &
                   cfile_nlff_2d,      &
                   cfile_2d_field,     &
                   cfile_3d_field,     &
                   cfile_3d_init,      &
                   cfile_restart,      &
                   cfile_slice_xz
namelist /nlist08/ nproc_x, nproc_y
namelist /nlist_ef/                   &
        eflux_x0, eflux_y0, eflux_z0, & ! Coordinate of torus center
        eflux_ph, & ! angle of emeging flux vector 
                    ! (0 correponds x-direction)
        eflux_vz, eflux_t0, eflux_t1, & 
                    ! Emerging speed, times of start and end
        eflux_b0, eflux_r0, eflux_rm, eflux_a0
namelist /nlist_lfff/  &
        lfff_b0, lfff_alpha, lfff_kly_pi, lfff_y0
namelist /nlist_divbeq0/  &
        ch2, cp2
!
contains
!-----------------------------------------------------------------------
      function chari2(in)
      character(3) chari2
      integer, INTENT(IN) :: in
! to make character of length 2 corresponding to 
! the integer 'in' with pedding.
!
      if(in.gt.99) then
         write(*,*) ' ++ warning: chari3, in > 999:',in
      end if
      write(chari2,'(I2.2)') in
      return
      end function chari2
!-----------------------------------------------------------------------
      function chari3(in)
      character(3) chari3
      integer, INTENT(IN) :: in
! to make character of length 3 corresponding to 
! the integer 'in' with pedding.
!
      if(in.gt.999) then
         write(*,*) ' ++ warning: chari3, in > 999:',in
      end if
      write(chari3,'(I3.3)') in
      return
      end function chari3
!-----------------------------------------------------------------------
      function chari4(in)
      character(4) chari4
      integer, INTENT(IN) :: in
! to make character of length 4 corresponding to 
! the integer 'in' with pedding.
!
      if(in.gt.9999) then
         write(*,*) ' ++ warning: chari3, in > 999:',in
      end if
      write(chari4,'(I4.4)') in
      return
      end function chari4

end module common






