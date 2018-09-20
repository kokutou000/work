program calc_decay_index
implicit none
integer,parameter :: nnx = 512, nny = 512, nnz = 256
integer,parameter :: nx0 = -1, ny0 = -1, nz0 = 0
integer,parameter :: ycut = 256
integer,parameter :: xcut = 256
integer :: i,j,k
real(8),allocatable :: By(:,:,:), Bx(:,:,:), Bz(:,:,:)
real(8),allocatable :: BB(:,:,:)
real(8) :: xc(nx0:nnx), yc(ny0:nny), zc(nz0:nnz)
real(8) :: nc(nx0:nnx,ny0:nny,nz0:nnz)=0d0
allocate(Bx(nx0:nnx,ny0:nny,nz0:nnz))
allocate(By(nx0:nnx,ny0:nny,nz0:nnz))
allocate(Bz(nx0:nnx,ny0:nny,nz0:nnz))
allocate(BB(nx0:nnx,ny0:nny,nz0:nnz))

print *, "test"
open(11,file="../By8bin_3d_init.dat",status="old",form="unformatted")
read(11) By
close(11)
open(12,file="../Bx8bin_3d_init.dat",status="old",form="unformatted")
read(12) Bx
close(12)
open(13,file="../Bz8bin_3d_init.dat",status="old",form="unformatted")
read(13) Bz
close(13)
open(21,file="../coord.xgc",status="old")
open(22,file="../coord.ygc",status="old")
open(23,file="../coord.zgc",status="old")
do i = nx0, nnx
read(21,'(f10.5)') xc(i)
enddo
do j = ny0, nny
read(22,'(f10.5)') yc(j)
enddo
do k = nz0, nnz
read(23,'(f10.5)') zc(k)
enddo
close(21)
close(22)
close(23)
!
BB = sqrt(Bx**2 + By**2 + Bz**2)
print *, "test2"
! nc = 0d0
! do i = nx0, nnx
! do k = nz0+1, nnz-1
!    nc(i,k) = - zc(k)*(by(i,ycut,k+1)-by(i,ycut,k-1)) &
!         &  / (2d0*0.5d0*(zc(k+1)-zc(k-1))*by(i,ycut,k))
! end do
! end do
nc = 0d0
do i = nx0, nnx
do j = ny0, nny
do k = nz0+1, nnz-1
   nc(i,j,k) = - zc(k)*(BB(i,j,k+1) - BB(i,j,k-1)) &
        & / (2d0*0.5d0*(zc(k+1)-zc(k-1))*BB(i,j,k))
enddo
enddo
enddo
!
print *, "test3"
!
open(110,file="bin_nc",status="replace",form="unformatted")
write(110) nc
close(110)
!
open(111,file="n_yz.dat",status="replace")
do j = ny0, nny
do k = nz0, nnz
write(111,'(f10.5)') nc(xcut,j,k)
enddo
enddo
close(111)
!
open(112,file="n_zx.dat",status="replace")
do i = nx0, nnx
do k = nz0, nnz
write(112,'(f10.5)') nc(i,ycut,k)
enddo
enddo
close(112)
!
open(113,file="n_at_center.dat",status="replace")
do k = nz0, nnz
write(113,'(f10.5,2x,f10.5)') zc(k), nc(xcut,ycut,k)
enddo
close(113)
end program calc_decay_index
