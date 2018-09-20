program convert_little_to_big
implicit none
integer,parameter :: nx1=512, ny1=512, nz1=256
integer,parameter :: nx0=-1, ny0=-1, nz0=0
real(8),allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:)
real(8),allocatable :: vx(:,:), vy(:,:), vz(:,:)
allocate(bx(nx0:nx1,ny0:ny1,nz0:nz1)); bx=0d0
allocate(by(nx0:nx1,ny0:ny1,nz0:nz1)); by=0d0
allocate(bz(nx0:nx1,ny0:ny1,nz0:nz1)); bz=0d0
allocate(vx(nx0:nx1,ny0:ny1)); vx=0d0
allocate(vy(nx0:nx1,ny0:ny1)); vy=0d0
allocate(vz(nx0:nx1,ny0:ny1)); vz=0d0

write(*,*) "READ FILE B3D_init: LITTLE_ENDIAN"
open(10,file="../vchange/B3D_init", status="old", form="unformatted")
read(10) bx, by, bz, vx, vy, vz
close(10)

write(*,*) "WRITE FILE B3D_init: BIG_ENDIAN"
open(11,file="B3D_init_BIG", status="replace", form="unformatted",convert="BIG_ENDIAN")
write(11) bx, by, bz, vx, vy, vz
close(11)

end program convert_little_to_big
