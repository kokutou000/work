program calc_vel2d
implicit none
integer,parameter :: NNX=50, NNY=50, NNZ=50
integer,parameter :: nx=512, ny=512, nz=256
integer :: i,j,k
real(8),allocatable :: bx3d(:,:,:), by3d(:,:,:), bz3d(:,:,:)
real(8),allocatable :: vx2d(:,:), vy2d(:,:), vz2d(:,:)
real(8),allocatable :: vxca(:,:), vyca(:,:), vzca(:,:)
real(8) :: xc(-1:nx), yc(-1:ny), zc(0:nz)
real(8) :: delx, dely, delz, epsz
real(8) :: dbdy = 0d0, dbdx = 0d0
real(8),parameter :: pi=3.1415926535
real(8) :: pi2, xi0, yj0, zk0
real(8) :: dxp, dxm, dyp, dym
integer :: ip, im, jp, jm
real(8) :: maxvv=0d0, maxbz=0d0
!character,parameter :: dirname*100="uniform_for_xdim"
!character,parameter :: dirname*100="hiroiv"
character,parameter :: dirname*100="testtmp"
!character,parameter :: dirname*100="nomask"
!character,parameter :: dirname*100="notanh"
allocate(bx3d(-1:nx,-1:ny,0:nz))
allocate(by3d(-1:nx,-1:ny,0:nz))
allocate(bz3d(-1:nx,-1:ny,0:nz))
allocate(vx2d(-1:nx,-1:ny))
allocate(vy2d(-1:nx,-1:ny))
allocate(vz2d(-1:nx,-1:ny))
allocate(vxca(-1:nx,-1:ny))
allocate(vyca(-1:nx,-1:ny))
allocate(vzca(-1:nx,-1:ny))

write(*,*) "READ FILE B3D_init"
open(10,file="../B3D_init",form="unformatted",status="old")
read(10) bx3d, by3d, bz3d, vx2d, vy2d, vz2d
close(10)
!
open(21,file="../coord.xgc",status="old")
do i = -1, nx
   read(21,'(f10.5)') xc(i)
end do
close(21)
open(22,file="../coord.ygc",status="old")
do j = -1, ny
   read(22,'(f10.5)') yc(j)
end do
close(22)
open(23,file="../coord.zgc",status="old")
do k = 0, nz
   read(23,'(f10.5)') zc(k)
end do
close(23)

! pi2 = pi*2.0d0
! epsz = 252d0
! delx = real(NNX)/real(nx)
! dely = real(NNY)/real(ny)
! delz = real(NNZ)/real(nz)
! do i = -1, nx
!    xi0 = delx * i
!    xc(i) = (xi0 - NNX*0.5d0)*0.03d0 + &
!         & ((4.0*(xi0 - NNX*0.5d0))**11.0)*62.08/(0.5d0*nx)**11.0
! end do
! do j = -1, ny
!    yj0 = dely * j
!    yc(j) = (yj0 - NNY*0.5d0)*0.01d0 + &
!         & ((4.0*(yj0 - NNY*0.5d0))**11.0)*63.36/(0.5*ny)**11.0
! end do
! do k = 0, nz
!    zk0 = delz * k
!    zc(k) = zk0 + epsz*sin(pi*zk0/NNZ + pi)/pi2
! end do

!write(*,*) vx2d(256,200:312)
vxca = 0.0d0
vyca = 0.0d0
vzca = vz2d
print *, bz3d(256,256:300,0)

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

do i = 0, nx-1
dxp = xc(i+1) - xc(i)
dxm = xc(i) - xc(i-1)
im = i-1; ip = i+1
do j = 0, ny-1
dyp = yc(j+1) - yc(j)
dym = yc(j) - yc(j-1)
jm = j-1; jp = j+1
!
! dbdy = ( bz3d(i,jp,0)*dym*dym &
!      &    + bz3d(i,j,0)*(dyp*dyp-dym*dym) &
!      &    - bz3d(i,jm,0)*dyp*dyp) &
!      &    /(dyp*dym*(dyp+dym))
! dbdx = ( bz3d(ip,j,0)*dxm*dxm &
!      &    + bz3d(i,j,0)*(dxp*dxp-dxm*dxm) &
!      &    - bz3d(im,j,0)*dxp*dxp) &
!      &    /(dxp*dxm*(dxp+dxm))
! vxca(i,j) = tanh(bz3d(i,j,0))*&
!      &     (tanh((bz3d(i,j,0)+0.8d0)/0.1d0)+1d0)*&
!      &    (-tanh((bz3d(i,j,0)-0.8d0)/0.1d0)+1d0)*&
!      &     dbdy/sqrt(dbdx**2+dbdy**2)
! vyca(i,j) =-tanh(bz3d(i,j,0))*&
!      &     (tanh((bz3d(i,j,0)+0.8d0)/0.1d0)+1d0)*&
!      &    (-tanh((bz3d(i,j,0)-0.8d0)/0.1d0)+1d0)*&
!      &     dbdx/sqrt(dbdx**2+dbdy**2)
if(bz3d(i,j,0).ne.0.0) then
!vxca(i,j) = bz3d(i,j,0)*&
!vxca(i,j) = &
vxca(i,j) = tanh(bz3d(i,j,0)/0.50d0)*&
     &    ( bz3d(i,jp,0)*dym*dym &
     &    + bz3d(i,j,0)*(dyp*dyp-dym*dym) &
     &    - bz3d(i,jm,0)*dyp*dyp) &
     &    /(dyp*dym*(dyp+dym))
!vyca(i,j) =-bz3d(i,j,0)*&
!vyca(i,j) =-&
vyca(i,j) =-tanh(bz3d(i,j,0)/0.50d0)*&
     &    ( bz3d(ip,j,0)*dxm*dxm &
     &    + bz3d(i,j,0)*(dxp*dxp-dxm*dxm) &
     &    - bz3d(im,j,0)*dxp*dxp) &
     &    /(dxp*dxm*(dxp+dxm))
else
vxca(i,j) = 0.0
vyca(i,j) = 0.0
end if
enddo
enddo
write(*,*) vxca(256,200:312)
!
!maxbz = maxval(abs(bz3d(:,:,0)))
!----- original vfield
!vxca = vxca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.8)/0.05))
!vyca = vyca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.8)/0.05))
!----- nearPIL0
!vxca = vxca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.4)/0.2))
!vyca = vyca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.4)/0.2))
!----- nearPIL1
!vxca = vxca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.3)/0.6))
!vyca = vyca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.3)/0.6))
!----- nearPIL1.1
!vxca = vxca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.3)/0.6))
!vyca = vyca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.3)/0.6))
!vxca = vxca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.95)/0.025))
!vyca = vyca * 0.5 * (1.0 - tanh((abs(bz3d(:,:,0)/maxbz)-0.95)/0.025))
!----- nomask
!vxca = vxca
!vyca = vyca
maxvv = maxval(sqrt(vxca**2+vyca**2))
vxca = vxca/maxvv
vyca = vyca/maxvv
vxca = vxca*0.01d0
vyca = vyca*0.01d0
!
!--- reset velocity
! vxca = 0d0
! vyca = 0d0
! !--- calc uniform v for x direction
! do j = 0, ny-1
! vxca(:,j) = 0.01*tanh(yc(j)/0.1)&
!      &*(1.0-tanh((yc(j)-1.0)/0.2))*0.5&
!      &*(1.0+tanh((yc(j)+1.0)/0.2))*0.5
! end do
!

open(20,file="output/"//trim(dirname)//"/vxca",status="replace",form="unformatted")
write(20) vxca
close(20)
open(20,file="output/"//trim(dirname)//"/vyca",status="replace",form="unformatted")
write(20) vyca
close(20)
open(20,file="output/"//trim(dirname)//"/vzca",status="replace",form="unformatted")
write(20) vzca
close(20)
open(20,file="output/"//trim(dirname)//"/vx2d",status="replace",form="unformatted")
write(20) vx2d
close(20)
open(20,file="output/"//trim(dirname)//"/vy2d",status="replace",form="unformatted")
write(20) vy2d
close(20)

write(*,*) vxca(256,200:312)

write(*,*) "OUTPUT FILE B3D_init"
open(30,file="output/"//trim(dirname)//"/B3D_init",form="unformatted")
write(30) bx3d, by3d, bz3d, vxca, vyca, vz2d
close(30)

end program calc_vel2d
