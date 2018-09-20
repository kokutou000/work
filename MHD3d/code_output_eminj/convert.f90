program file_for_vapor3d
implicit none
integer :: i,j,k
integer,parameter :: nxall = 512
integer,parameter :: nyall = 512
integer,parameter :: nzall = 256
!integer,parameter :: nxmin = 32, nxmax = 480 ! -8<x<8
!integer,parameter :: nymin = 32, nymax = 480 ! -8<y<8
!integer,parameter :: nzmin = 0,  nzmax = 205 !   0<z<8
integer :: tsini, tsfin
integer :: nxmin, nxmax, nymin, nymax, nzmin, nzmax
integer :: ts = 0
character*3 :: cts = "000"
real(4),allocatable :: bxall(:,:,:), byall(:,:,:), bzall(:,:,:), &
     & vxall(:,:,:), vyall(:,:,:), vzall(:,:,:), roall(:,:,:)
real(4),dimension(0:nxall,0:nyall) :: bz2dall
real(8),dimension(0:nxall) :: xcall
real(8),dimension(0:nyall) :: ycall
real(8),dimension(0:nzall) :: zcall
character*100 :: dirbinfile="../binfile/"
!
!--- 1d
real(4),allocatable :: xc(:), yc(:), zc(:)
real(4),allocatable :: engb(:), engv(:)
real(4),allocatable :: engbx(:), engby(:), engbz(:)
real(4),allocatable :: tmpeb2(:,:,:), tmpev2(:,:,:)
real(4),allocatable :: tmpeb1(:,:), tmpev1(:,:)
real(4),allocatable :: tmpeb(:), tmpev(:)
!
real(4),allocatable :: testv_norho(:)
!--- 2d
real(4),allocatable :: bz2d(:,:)
!--- 3d
real(4),allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:), &
     & vx(:,:,:), vy(:,:,:), vz(:,:,:), ro(:,:,:)
!--- get info ---
open(99,file="../code_output/info",status="old")
read(99,*) tsini;read(99,*) tsfin
read(99,*) nxmin;read(99,*) nxmax
read(99,*) nymin;read(99,*) nymax
read(99,*) nzmin;read(99,*) nzmax
close(99)
print *, tsini, tsfin, nxmin, nxmax, nymin, nymax, nzmin, nzmax
!----------------
!
allocate(bxall(0:nxall,0:nyall,0:nzall))
allocate(byall(0:nxall,0:nyall,0:nzall))
allocate(bzall(0:nxall,0:nyall,0:nzall))
allocate(vxall(0:nxall,0:nyall,0:nzall))
allocate(vyall(0:nxall,0:nyall,0:nzall))
allocate(vzall(0:nxall,0:nyall,0:nzall))
allocate(roall(0:nxall,0:nyall,0:nzall))
!
allocate(xc(nxmin:nxmax));allocate(yc(nymin:nymax));allocate(zc(nzmin:nzmax))
allocate(engb(tsini:tsfin)); allocate(engv(tsini:tsfin))
allocate(engbx(tsini:tsfin)); allocate(engby(tsini:tsfin))
allocate(engbz(tsini:tsfin))
allocate(tmpeb2(nymin:nymax,nzmin:nzmax,tsini:tsfin))
allocate(tmpev2(nymin:nymax,nzmin:nzmax,tsini:tsfin))
allocate(tmpeb1(nzmin:nzmax,tsini:tsfin))
allocate(tmpev1(nzmin:nzmax,tsini:tsfin))
allocate(tmpeb(tsini:tsfin))
allocate(tmpev(tsini:tsfin))
allocate(testv_norho(tsini:tsfin))
allocate(bz2d(nxmin:nxmax,nymin:nymax))
allocate(bx(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(by(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(bz(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(vx(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(vy(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(vz(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
allocate(ro(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
engb = 0d0; engv = 0d0
engbx = 0d0; engby = 0d0; engbz = 0d0
tmpeb2 = 0d0; tmpev2 = 0d0
tmpeb1 = 0d0; tmpev1 = 0d0
tmpeb = 0d0; tmpev = 0d0
testv_norho = 0d0

dirbinfile=trim(dirbinfile)

!--- coordinate
open(11,file=trim(dirbinfile)//"coord.xgc",status="old")
read(11,'(e25.16)') xcall
close(11)
open(12,file=trim(dirbinfile)//"coord.ygc",status="old")
read(12,'(e25.16)') ycall
close(12)
open(13,file=trim(dirbinfile)//"coord.zgc",status="old")
read(13,'(e25.16)') zcall
close(13)
!
xc = xcall(nxmin:nxmax)
yc = ycall(nymin:nymax)
zc = zcall(nzmin:nzmax)
!
open(11,file="coord.xgc",status="replace")
write(11,'(e25.16)') xc
close(11)
open(12,file="coord.ygc",status="replace")
write(12,'(e25.16)') yc
close(12)
open(13,file="coord.zgc",status="replace")
write(13,'(e25.16)') zc
close(13)
!---------------
do ts = tsini, tsfin
   write(cts,'(I3.3)') ts
   write(*,*) "now time step: ", cts
!--- 2dvals
   open(21,file=trim(dirbinfile)//"Bz2d."//cts,form="unformatted",status="old")
   read(21) bz2dall
   close(21)
!
   bz2d(nxmin:nxmax,nymin:nymax) = bz2dall(nxmin:nxmax,nymin:nymax)
!
   open(21,file="Bz2d_R."//cts,form="unformatted",status="replace")
   write(21) bz2d
   close(21)
!---------------
!--- 3dvals
   open(31,file=trim(dirbinfile)//"Bxbin_3d."//cts,form="unformatted",status="old")
   open(32,file=trim(dirbinfile)//"Bybin_3d."//cts,form="unformatted",status="old")
   open(33,file=trim(dirbinfile)//"Bzbin_3d."//cts,form="unformatted",status="old")
   open(34,file=trim(dirbinfile)//"Vxbin_3d."//cts,form="unformatted",status="old")
   open(35,file=trim(dirbinfile)//"Vybin_3d."//cts,form="unformatted",status="old")
   open(36,file=trim(dirbinfile)//"Vzbin_3d."//cts,form="unformatted",status="old")
   open(37,file=trim(dirbinfile)//"Robin_3d."//cts,form="unformatted",status="old")
   read(31) bxall
   read(32) byall
   read(33) bzall
   read(34) vxall
   read(35) vyall
   read(36) vzall
   read(37) roall
   close(31)
   close(32)
   close(33)
   close(34)
   close(35)
   close(36)
   close(37)
!
   bx(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & bxall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   by(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & byall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   bz(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & bzall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   vx(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & vxall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   vy(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & vyall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   vz(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & vzall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   ro(nxmin:nxmax,nymin:nymax,nzmin:nzmax) = &
        & roall(nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!
   open(31,file="Bxbin_3d_R."//cts,form="unformatted", status="replace")
   open(32,file="Bybin_3d_R."//cts,form="unformatted", status="replace")
   open(33,file="Bzbin_3d_R."//cts,form="unformatted", status="replace")
   open(34,file="Vxbin_3d_R."//cts,form="unformatted", status="replace")
   open(35,file="Vybin_3d_R."//cts,form="unformatted", status="replace")
   open(36,file="Vzbin_3d_R."//cts,form="unformatted", status="replace")
   open(37,file="Robin_3d_R."//cts,form="unformatted", status="replace")
   write(31) bx
   write(32) by
   write(33) bz
   write(34) vx
   write(35) vy
   write(36) vz
   write(37) ro
   close(31)
   close(32)
   close(33)
   close(34)
   close(35)
   close(36)
   close(37)
!
!--- integration for z>0
   do k = nzmin+1, nzmax
   do j = nymin, nymax
   do i = nxmin, nxmax
      engb(ts) = engb(ts) + 0.5d0*(bx(i,j,k)*bx(i,j,k)+by(i,j,k)*by(i,j,k)+bz(i,j,k)*bz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k-1))*0.5d0
      engv(ts) = engv(ts) + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k))*ro(i,j,k) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k-1))*0.5d0
      testv_norho(ts) = testv_norho(ts) + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k-1))*0.5d0
!      engbx(ts) = 0.5d0*(bx(i,j,k)**2)
!      engby(ts) = 0.5d0*(by(i,j,k)**2)
!      engbz(ts) = 0.5d0*(bz(i,j,k)**2)
   end do
   end do
   end do
!--- integration for z=0
   k = nzmin
   do j = nymin, nymax
   do i = nxmin, nxmax
      engb(ts) = engb(ts) + 0.5d0*(bx(i,j,k)*bx(i,j,k)+by(i,j,k)*by(i,j,k)+bz(i,j,k)*bz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k))*0.5d0
      engv(ts) = engv(ts) + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k))*ro(i,j,k) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k))*0.5d0
      testv_norho(ts) = testv_norho(ts) + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k))*0.5d0
   end do
   end do
!--- calc by 2-bungi hou
!--- integraion for z>0
   do k = nzmin+1, nzmax
   do j = nymin, nymax
   do i = nxmin, nxmax
      tmpeb2(j,k,ts) = tmpeb2(j,k,ts) + 0.5d0*(bx(i,j,k)*bx(i,j,k)+by(i,j,k)*by(i,j,k)+bz(i,j,k)*bz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k-1))*0.5d0
      tmpev2(j,k,ts) = tmpev2(j,k,ts) &
           & + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k))*ro(i,j,k) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k-1))*0.5d0
   enddo   
   enddo
   enddo
!--- integraion for z=0
   k = nzmin
   do j = nymin, nymax
   do i = nxmin, nxmax
      tmpeb2(j,k,ts) = tmpeb2(j,k,ts) + 0.5d0*(bx(i,j,k)*bx(i,j,k)+by(i,j,k)*by(i,j,k)+bz(i,j,k)*bz(i,j,k)) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k))*0.5d0
      tmpev2(j,k,ts) = tmpev2(j,k,ts) &
           & + 0.5d0*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k))*ro(i,j,k) &
           &*(xcall(i+1)-xcall(i-1))*0.5d0 &
           &*(ycall(j+1)-ycall(j-1))*0.5d0 &
           &*(zcall(k+1)-zcall(k))*0.5d0
   enddo   
   enddo
!--- integraion for y dir
   do k = nzmin, nzmax
   do j = nymin, nymax
      tmpeb1(k,ts) = tmpeb1(k,ts) + tmpeb2(j,k,ts)
      tmpev1(k,ts) = tmpev1(k,ts) + tmpev2(j,k,ts)
   enddo
   enddo
!--- integration for z dir
   do k = nzmin, nzmax
      tmpeb(ts) = tmpeb(ts) + tmpeb1(k,ts)
      tmpev(ts) = tmpev(ts) + tmpev1(k,ts)
   enddo
!---
!--- end do loop for time step
end do
!-----
open(41,file="engb_R",status="replace")
open(42,file="engv_R",status="replace")
do ts = tsini, tsfin
!   write(41,*) ts, engb(ts)
!   write(42,*) ts, engv(ts)
   write(41,*) ts, tmpeb(ts)
   write(42,*) ts, tmpev(ts)
end do
close(41)
close(42)
print *, engb
print *, engv
print *, tmpeb
print *, tmpev
!--- output test calc v*v *dxdydz 
!--- so kinetic energy without no rho term
open(51,file="checkv_R",status="replace")
do ts = tsini, tsfin
write(51,*) ts, testv_norho(ts)
end do
close(51)
end program file_for_vapor3d
