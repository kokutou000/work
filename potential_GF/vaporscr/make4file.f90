program make_bin4_files
implicit none
integer,parameter :: nx = 512, ny = 512, nz = 256
integer,parameter :: nx0 = -1, ny0 = -1, nz0 = 0
real(8),allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
allocate(Bx(nx0:nx,ny0:ny,nz0:nz)); Bx=0d0
allocate(By(nx0:nx,ny0:ny,nz0:nz)); By=0d0
allocate(Bz(nx0:nx,ny0:ny,nz0:nz)); Bz=0d0

!
open(21,file="../Bx8bin_3d_init.dat",status="old",form="unformatted")
read(21) Bx
close(21)
open(22,file="../By8bin_3d_init.dat",status="old",form="unformatted")
read(22) By
close(22)
open(23,file="../Bz8bin_3d_init.dat",status="old",form="unformatted")
read(23) Bz
close(23)
!
open(31,file="Bx4bin3d",status="replace",form="unformatted")
write(31) real(Bx)
close(31)
open(32,file="By4bin3d",status="replace",form="unformatted")
write(32) real(By)
close(32)
open(33,file="Bz4bin3d",status="replace",form="unformatted")
write(33) real(Bz)
close(33)
!
end program make_bin4_files
