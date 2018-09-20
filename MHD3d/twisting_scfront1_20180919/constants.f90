! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
module constants
! -----------------------------------------------------------
!                      MODULE CONSTANTS
! -----------------------------------------------------------
implicit none

! === grid number ===
!MEMO! (256=32*8, 256=64*4)
  integer, parameter :: NX=64
  integer, parameter :: NY=64
  integer, parameter :: NZ=256

! === grid number for the 3D initial state ===
  integer, parameter :: NX3D=NX*4
  integer, parameter :: NY3D=NY*4
  integer, parameter :: NZ3D=NZ

! === grid number for the 2D initial state ===
  integer, parameter :: NY2D=512
  integer, parameter :: NZ2D=NZ

! === grid number for NLFF ===
  integer, parameter :: NX_NLFF=512
  integer, parameter :: NY_NLFF=256

! === sub-grid number ===
  integer, parameter :: nxm1=NX-1
  integer, parameter :: nxp1=NX+1
  integer, parameter :: nym1=NY-1
  integer, parameter :: nyp1=NY+1
  integer, parameter :: nzm1=NZ-1
  integer, parameter :: nzp1=NZ+1

! === precision type ===
integer, parameter :: DP = kind(1.0d0)
integer, parameter :: SP = kind(1.0)
integer, parameter :: DPC = kind((1.0d0,1.0d0))

! === mathematical constants ===
real(DP), parameter :: PI = 3.1415926535897932384626433_DP
real(DP), parameter :: PI2 = 2*PI
complex(DPC), parameter :: IUNIT = (0.0_DP, 1.0_DP)

! === I/O file number
integer, parameter :: FILE_SYSOUT       = 06 ! sysout file
integer, parameter :: FILE_RUN_NUMBER   = 09 ! number list calculated so far
integer, parameter :: FILE_NAMELIST     = 10 ! namelist
integer, parameter :: FILE_COORDINATE_X = 11 ! coordinate x with '(6e15.5)'
integer, parameter :: FILE_COORDINATE_Y = 12 ! coordinate y with '(6e15.5)'
integer, parameter :: FILE_COORDINATE_Z = 13 ! coordinate z with '(6e15.5)'
integer, parameter :: FILE_COORDINATE_Z_ALL = 14 ! coordinate z with '(6e15.5)'
                                             ! given by NLFF model.
integer, parameter :: FILE_EIGENMODE    = 21 ! eigenmode perturbation
integer, parameter :: FILE_OUTPUT_LIST  = 50 ! output list
integer, parameter :: FILE_TIME_LIST    = 60 ! sequential list (time,loop,etc.)
integer, parameter :: FILE_NLFF         = 62 ! NLFF (3d whole domain data)
integer, parameter :: FILE_2D_FIELD     = 65 ! initial 2d equilibrium
integer, parameter :: FILE_2D_XY        = 66 ! 2d (x,y) B&V field on the bottom
integer, parameter :: FILE_3D_INIT      = 67 ! initial 3d initial condition
integer, parameter :: FILE_3D_FIELD     = 70 ! 3d results 
integer, parameter :: FILE_RESTART      = 80 ! 3d last result for restart
integer, parameter :: FILE_SLICE_XZ     = 90 ! 2d slice data

end module constants
