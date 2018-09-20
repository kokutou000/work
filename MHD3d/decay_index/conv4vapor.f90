program convert_n_yz_for_vapor
implicit none
integer,parameter :: nx1=512, ny1=512, nz1=256
integer,parameter :: nx0=-1, ny0=-1, nz0=0
real(8) :: n_bin(nx0:nx1,ny0:ny1,nz0:nz1)
integer,parameter :: nnx0=32, nnx1=480
integer,parameter :: nny0=32, nny1=480
integer,parameter :: nnz0=0, nnz1=205
!real(4) :: n_yz_vapor4(nny0:nny1,nnz0:nnz1)
real(4) :: n_3d_vapor4(nnx0:nnx1,nny0:nny1,nnz0:nnz1)
integer :: i,j,k

open(10,file="bin_nc",status="old",form="unformatted")
read(10) n_bin
close(10)

do i = nnx0, nnx1
do j = nny0, nny1
do k = nnz0, nnz1
   n_3d_vapor4(i,j,k) = n_bin(i,j,k)
end do
end do
end do

open(20,file="vapor_n_3d",status="replace",form="unformatted")
write(20) n_3d_vapor4
close(20)

end program convert_n_yz_for_vapor
