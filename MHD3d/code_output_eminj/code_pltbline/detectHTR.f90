program calc_flux_from_high_tw_region
implicit none
integer :: ts=0, i, j, k, tsini, tsfin
integer,parameter :: nx = 512, ny = 512, nz = 256
integer :: nnx,nny,nnz,nxmin,nxmax,nymin,nymax,nzmin,nzmax
real(8),parameter :: tw_crit = -0.5
real(8),parameter :: r_em = 0.4
real(4),allocatable :: tw2d(:,:)
real(4),allocatable :: bz3d(:,:,:)
real(8),allocatable :: flux(:), sumFRy(:), sumFRz(:), aveFRy(:), aveFRz(:)
real(8),allocatable :: rawxp(:), rawyp(:), rawzp(:)
integer,allocatable :: count(:)
character*3 :: cts="000"
character*100 :: dirbinfile="../../../../3dbin/"
!--- get info ---
open(100,file="../../../../code_output/info",status="old")
read(100,*) tsini; read(100,*) tsfin
read(100,*) nxmin; read(100,*) nxmax
read(100,*) nymin; read(100,*) nymax
read(100,*) nzmin; read(100,*) nzmax
close(100)
nnx = nxmax - nxmin; nny = nymax - nymin; nnz = nzmax - nzmin

allocate(tw2d(0:nnx,0:nny)); tw2d=0d0
allocate(bz3d(0:nnx,0:nny,0:nnz)); bz3d=0d0
allocate(flux(tsini:tsfin)); flux=0d0
allocate(rawxp(0:nnx)); allocate(rawyp(0:nny)); allocate(rawzp(0:nnz))
allocate(sumFRy(tsini:tsfin)); sumFRy=0d0
allocate(sumFRz(tsini:tsfin)); sumFRz=0d0
allocate(aveFRy(tsini:tsfin)); aveFRy=0d0
allocate(aveFRz(tsini:tsfin)); aveFRz=0d0
allocate(count(tsini:tsfin)); count=0
open(11,file=trim(dirbinfile)//"coord.xgc",status="old")
open(12,file=trim(dirbinfile)//"coord.ygc",status="old")
open(13,file=trim(dirbinfile)//"coord.zgc",status="old")
read(11,'(e25.16)') rawxp
read(12,'(e25.16)') rawyp
read(13,'(e25.16)') rawzp
close(11)
close(12)
close(13)

do ts = tsini, tsfin
   write(cts,'(I3.3)') ts
   write(*,*) "now time step: ", cts
!
   open(31,file="../../Tw2dbin."//cts,form="unformatted",status="old")
   read(31) tw2d
   close(31)
   open(32,file="../../../../3dbin/Bzbin_3d_R."//cts,form="unformatted",status="old")
   read(32) Bz3d
   close(32)
!---
   do j = 0, nny
   do i = 0, nnx
      if(sqrt(rawxp(i)**2+rawyp(j)**2).ge.r_em) then
      if(tw2d(i,j).le.tw_crit) then
         flux(ts) = flux(ts) &
              &+ abs(bz3d(i,j,0))*0.5d0&       !!! for 2-time foot points
              &*(rawyp(j+1)-rawyp(j-1))*0.5d0& !!! for half length
              &*(rawxp(i+1)-rawxp(i-1))*0.5d0  !!! for half length
!         count(ts) = count(ts) + 1
!         sumFRy(ts) = sumFRy(ts) + rawyp(j)
!         sumFRz(ts) = sumFRz(ts) + rawzp(k)
      end if
      end if
   end do
   end do
!   aveFRy(ts) = sumFRy(ts)/real(count(ts))
!   aveFRz(ts) = sumFRz(ts)/real(count(ts))
!
end do
!
open(201,file="hightwflux.txt")
do ts = tsini, tsfin
   write(201,*) flux(ts)
end do
close(201)
!
!open(210,file="posiFR.txt")
!do ts = tsini, tsfin
!   write(210,'(i5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)') count(ts), aveFRy(ts), aveFRz(ts), sumFRy(ts), sumFRz(ts)
!end do
!close(210)
!
end program calc_flux_from_high_tw_region
