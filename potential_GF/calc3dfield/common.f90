module common_val
use parameter_all
implicit none
!
!
!real(8),dimension(nx0:nx,ny0:ny,nz0:nz) :: Bx, By, Bz
real(8),allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
!
!
!real(8),dimension(nx0m1:nx+1,ny0m1:ny+1,nz0m1:nz+1) :: Psi0
!
!
real(8),dimension(nx0:nx,ny0:ny) :: Bz_init
real(8),dimension(nx0:nx,ny0:ny) :: Vx_init, Vy_init, Vz_init
real(8),dimension(nx0:nx,ny0:ny) :: ck1Bz, ck2Bz
!
real(8),dimension(nx0:nx) :: xc
real(8),dimension(ny0:ny) :: yc
real(8),dimension(nz0:nz) :: zc
!
character(8) :: date0="0"
character(10) :: time0="0"
character(25) :: dirname="0"
!
integer :: i,j,k
!
real(8) :: delx=0d0, dely=0d0, delz=0d0
real(8) :: epsx=0d0, epsy=0d0, epsz=0d0
!
real(8) :: alphax, alphay, alphaz
real(8) :: betax, betay, betaz
!
end module common_val
