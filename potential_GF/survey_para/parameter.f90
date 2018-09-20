module parameter_all
!--- BOX Length for each dimension
!integer,parameter :: NNX=20, NNY=20, NNZ=20
!integer,parameter :: NNX=4, NNY=2, NNZ=4
!integer,parameter :: NNX=8, NNY=8, NNZ=8
!integer,parameter :: NNX=10, NNY=10, NNZ=10
!integer,parameter :: NNX=128, NNY=128, NNZ=128
integer,parameter :: NNX=50, NNY=50, NNZ=50
!integer,parameter :: NNX=100, NNY=100, NNZ=2   !!! debug for openmp
!---
!integer,parameter :: nx = 16, ny = 16, nz = 256  !!! for debug1
!integer,parameter :: nx = 128, ny = 128, nz = 64  !!! for debug2 about 2min
!integer,parameter :: nx = 256, ny = 256, nz = 128 !!! about 30min
integer,parameter :: nx = 512, ny = 512, nz = 256
!integer,parameter :: nx = 512, ny = 512, nz = 8  !!! debug for openmp
!integer,parameter :: nx = 576, ny = 576, nz = 320
!---
integer,parameter :: nx0 = -1, ny0 = -1, nz0 = 0
integer,parameter :: nx0m1 = -2, ny0m1 = -2, nz0m1 = -1
real(8),parameter :: rnx = nx*1d0, rny = ny*1d0, rnz = nz*1d0
real(8),parameter :: Bs0 = 1d0
real(8),parameter :: at = 4.493409458
!real(8),parameter :: a_suf = 0.5d0
real(8),parameter :: a_suf = 2.0d0
real(8),parameter :: alp = at/a_suf
real(8),parameter :: zcut = a_suf*0.6d0
real(8),parameter :: ydecay = a_suf*0.5
real(8),parameter :: pi=3.14159265358979329d0,pi2=6.2831853071795864d0
integer,parameter :: flatten_switch=1
!----- para for point source
real(8),parameter :: dep = 0.50d0
real(8),parameter :: lpil = 0.25d0
real(8),parameter :: qps = 1.0d0
!----- para for sech field
real(8),parameter :: decx = 6.0d0
real(8),parameter :: decy = 0.3d0
real(8),parameter :: lpilsec = 0.3d0
real(8),parameter :: dtany = 0.05d0
real(8),parameter :: xmask = 20.0
!----- para mask field
real(8),parameter :: pmask = 1.5d0
real(8),parameter :: dmask = 0.3d0
!--- for 100-100-100 512x512x256
real(8),parameter :: nbeki=11
real(8),parameter :: para_x0=8.0, para_ix0=224.0
real(8),parameter :: para_y0=8.0, para_jy0=224.0
real(8),parameter :: para_z0=8.0, para_kz0=205.0
!real(8),parameter :: para_z0=8.0, para_kz0=210.0
!real(8),parameter :: para_z0=1.0, para_kz0=6.0 !!! debug for openmp
!---- for 150-150-150 576x576x288
!real(8),parameter :: nbeki=11
!real(8),parameter :: para_x0=4.0, para_ix0=48.0
!real(8),parameter :: para_y0=4.0, para_jy0=48.0
!real(8),parameter :: para_z0=6.0, para_kz0=48.0
!--- for debug 100-100-100 128x128x64
!real(8),parameter :: nbeki=11
!real(8),parameter :: para_x0=6.0, para_ix0=48.0
!real(8),parameter :: para_y0=6.0, para_jy0=48.0
!real(8),parameter :: para_z0=8.0, para_kz0=48.0
end
