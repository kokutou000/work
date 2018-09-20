program calculate_dbzdx_and_dbzdy
implicit none
integer,parameter :: nx = 512, ny = 512, nz = 256
integer,parameter :: nx0 = -1, ny0 = -1, nz0 = 0
integer :: i,j,k, num
real(8) :: xc(nx0:nx), yc(ny0:ny), zc(nz0:nz)
real(8),allocatable :: dbzdx(:,:), dbzdy(:,:)
real(8),allocatable :: dbzdytanhbz(:,:)
real(8) :: dxp, dxm, dyp, dym
integer :: nbz ! number of array of Bz 0.0 to 1.0
integer :: maxnbz, minnbz
real(8),allocatable :: aarr(:), barr(:)
real(8),allocatable :: dbzdycal(:,:)
real(8),allocatable :: newf(:,:)
integer :: tmp0
real(8),allocatable :: vxtmp(:,:)
real(8),allocatable :: vytmp(:,:)
real(8) :: vvmax=0d0
real(8),allocatable :: bx3d(:,:,:), by3d(:,:,:), bz3d(:,:,:)
real(8),allocatable :: vx2d(:,:), vy2d(:,:), vz2d(:,:)
real(8),allocatable :: bxd0(:,:,:), byd0(:,:,:), bzd0(:,:,:)

allocate(bx3d(nx0:nx,ny0:ny,nz0:nz)); bx3d=0d0
allocate(by3d(nx0:nx,ny0:ny,nz0:nz)); by3d=0d0
allocate(bz3d(nx0:nx,ny0:ny,nz0:nz)); bz3d=0d0
allocate(bxd0(nx0:nx,ny0:ny,nz0:nz)); bxd0=0d0
allocate(byd0(nx0:nx,ny0:ny,nz0:nz)); byd0=0d0
allocate(bzd0(nx0:nx,ny0:ny,nz0:nz)); bzd0=0d0
allocate(vx2d(nx0:nx,ny0:ny)); vx2d=0d0
allocate(vy2d(nx0:nx,ny0:ny)); vy2d=0d0
allocate(vz2d(nx0:nx,ny0:ny)); vz2d=0d0

allocate(dbzdx(nx0:nx,ny0:ny)); dbzdx=0d0
allocate(dbzdy(nx0:nx,ny0:ny)); dbzdy=0d0
allocate(dbzdytanhbz(nx0:nx,ny0:ny)); dbzdytanhbz=0d0
allocate(dbzdycal(nx0:nx,ny0:ny)); dbzdycal=0d0
allocate(newf(nx0:nx,ny0:ny)); newf=0d0
allocate(vxtmp(nx0:nx,ny0:ny)); vxtmp=0d0
allocate(vytmp(nx0:nx,ny0:ny)); vytmp=0d0

write(*,*) "READ FILE B3D_init"
open(10,file="../B3D_init",form="unformatted",status="old")
read(10) bx3d, by3d, bz3d, vx2d, vy2d, vz2d
close(10)
!
open(11,file="../Bx8bin_3d_init.dat",form="unformatted",status="old")
open(12,file="../By8bin_3d_init.dat",form="unformatted",status="old")
open(13,file="../Bz8bin_3d_init.dat",form="unformatted",status="old")
read(11) bxd0; read(12) byd0; read(13) bzd0
close(11); close(12); close(13)
!
open(21,file="../coord.xgc",status="old")
open(22,file="../coord.ygc",status="old")
open(23,file="../coord.zgc",status="old")
do i = nx0, nx
   read(21,'(f10.5)') xc(i)
enddo
do j = ny0 ,ny
   read(22,'(f10.5)') yc(j)
enddo
do k = nz0 ,nz
   read(23,'(f10.5)') zc(k)
enddo
close(21);close(22);close(23)
!
do i = nx0+1, nx-1
   dxp = xc(i+1)-xc(i)
   dxm = xc(i)-xc(i-1)
do j = ny0+1, ny-1
   dyp = yc(j+1)-yc(j)
   dym = yc(j)-yc(j-1)
   dbzdx(i,j) = (bzd0(i+1,j,0)*dxm*dxm &
             & + bzd0(i  ,j,0)*(dxp*dxp-dxm*dxm) & 
             & - bzd0(i-1,j,0)*dxp*dxp) &
             & /(dxp*dxm*(dxp+dxm))
   dbzdy(i,j) = (bzd0(i,j+1,0)*dym*dym &
             & + bzd0(i,j  ,0)*(dyp*dyp-dym*dym) & 
             & - bzd0(i,j-1,0)*dyp*dyp) &
             & /(dyp*dym*(dyp+dym))
   dbzdytanhbz(i,j) = dbzdy(i,j)*tanh(bzd0(i,j,0))
enddo
enddo
!
open(41,file="dbzdx.dat",status="replace")
open(42,file="dbzdy.dat",status="replace")
open(43,file="dbzdytanhbz.dat",status="replace")
do i = nx0, nx
do j = ny0, ny
write(41,'(f10.6,2x,f10.6,2x,f12.8)') xc(i), yc(j), dbzdx(i,j)
write(42,'(f10.6,2x,f10.6,2x,f12.8)') xc(i), yc(j), dbzdy(i,j)
write(43,'(f10.6,2x,f10.6,2x,f12.8)') xc(i), yc(j), dbzdytanhbz(i,j)
end do
end do
close(41); close(42); close(43)
!
open(51,file="bz_vs_dbzdy.dat",status="replace")
do j = ny/2, ny
write(51,'(f12.8,2x,f12.8,2x,f12.8)') yc(j), bzd0(256,j,0), dbzdy(256,j)
end do
close(51)
!
maxnbz = maxloc(bzd0(nx/2,:,0),1); minnbz = minloc(abs(bzd0(nx/2,:,0)),1)
print *, maxnbz, minnbz, maxval(bzd0(nx/2,:,0)), minval(abs(bzd0(nx/2,:,0)))
nbz = maxnbz - minnbz + 1
allocate(aarr(1:nbz)); allocate(barr(1:nbz))
aarr = 0d0; barr = 0d0
!
do num = 1, nbz-1
   aarr(num) = (dbzdy(256,num+(ny/2)) - dbzdy(256,num+(ny/2)-1))&
        &/(bzd0(256,num+(ny/2),0)-bzd0(256,num+(ny/2)-1,0))
   barr(num) = (dbzdy(256,num+(ny/2)-1)*bzd0(256,num+(ny/2),0)&
        &-dbzdy(256,num+(ny/2))*bzd0(256,num+(ny/2)-1,0))&
        &/(bzd0(256,num+(ny/2),0)-bzd0(256,num+(ny/2)-1,0))
enddo
open(52,file="interpolation1d.dat",status="replace")
do num = 1, nbz-1
write(52,'(f12.6,2x,f12.6,2x,f12.6)') bzd0(256,num+(ny/2)-1,0), aarr(num), barr(num)
end do
close(52)
!
open(200,file="debug.dat",status="replace")
do i = nx0, nx
do j = ny0, ny
   tmp0 = 0
   do num = 1, nbz-1
      if((abs(bzd0(i,j,0)) .ge. abs(bzd0(256,num+(ny/2)-1,0)))&
           &.and. (abs(bzd0(i,j,0)) .lt. abs(bzd0(256,num+(ny/2),0))) ) then
         dbzdycal(i,j) = aarr(num)*abs(bzd0(i,j,0)) + barr(num)
!         print *, "debug", bzd0(256,num+(ny/2)-1,0), bzd0(256,num+(ny/2),0),dbzdycal(i,j)
         tmp0 = 1
      else if(abs(bzd0(i,j,0)).eq.abs(bzd0(256,256+nbz-1,0))) then
         dbzdycal(i,j) = 1d0
         tmp0 = 1
!       else
!          dbzdycal(i,j) = 0d0
!  !        print *, "value error and stop program!"
! !         print *, "@(i,j)=",i,j
!          write(200,*) i,j,bzd0(i,j,0)
! !         stop
      end if
      if(tmp0.eq.1) exit
      if(num.eq.nbz-1) write(200,*) i,j,bzd0(i,j,0)
   end do
enddo
enddo
close(200)
!
open(101,file="dbzdycal.dat",status="replace")
do i = nx0, nx
do j = ny0, ny
write(101,*) bzd0(i,j,0), dbzdycal(i,j)
enddo
enddo
close(101)
!
do i = nx0, nx
do j = ny0, ny
!newf(i,j) = 2d0*dbzdycal(i,j) / (1-tanh((abs(bzd0(i,j,0))-0.95)/0.05))
!newf(i,j) = dbzdycal(i,j) + 0.5d0*(1d0+tanh((abs(bzd0(i,j,0))-0.98)/0.02))
!newf(i,j) = dbzdycal(i,j) / (0.5d0*(1d0-tanh((abs(bzd0(i,j,0))-0.95)/0.01))&
!& + dbzdycal(i,j) * 0.5d0 * (1d0+tanh((abs(bzd0(i,j,0))-0.95)/0.01)))
newf(i,j) = dbzdycal(i,j) / (0.5d0*(1d0-tanh((abs(bzd0(i,j,0))-0.95)/0.01))&
& + dbzdycal(i,j) * 0.5d0 * (1d0+tanh((abs(bzd0(i,j,0))-0.95)/0.01)))
enddo
enddo
open(102,file="newfunction.dat",status="replace")
do i = nx0, nx
do j = ny0, ny
write(102,*) bzd0(i,j,0), newf(i,j)
enddo
enddo
close(102)
!
do i = nx0, nx
do j = ny0, ny
   vxtmp(i,j) = tanh(bzd0(i,j,0))*dbzdy(i,j) / newf(i,j)
   vytmp(i,j) =-tanh(bzd0(i,j,0))*dbzdx(i,j) / newf(i,j)
enddo
enddo
vvmax = maxval(sqrt(vxtmp**2 + vytmp**2))
vxtmp = 0.01 * vxtmp / vvmax
vytmp = 0.01 * vytmp / vvmax
open(111,file="vxtmp.dat",status="replace")
open(112,file="vytmp.dat",status="replace")
do i = nx0, nx
do j = ny0, ny
write(111,*) vxtmp(i,j)
write(112,*) vytmp(i,j)
enddo
enddo
close(111); close(112)
!
vx2d = vxtmp; vy2d = vytmp
write(*,*) "OUTPUT FILE B3D_init"
open(200,file="B3D_init",form="unformatted")
!write(200) bx3d, by3d, bz3d, vx2d, vy2d, vz2d
write(200) bxd0, byd0, bzd0, vx2d, vy2d, vz2d
close(200)
end program calculate_dbzdx_and_dbzdy
