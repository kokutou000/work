program main
  implicit none
  integer,parameter :: nx=64
  integer,parameter :: ny=64
  integer,parameter :: nz=256
  integer,parameter :: nproc_x=8, nproc_y=8
  integer,parameter :: nx3d=nx*nproc_x
  integer,parameter :: ny3d=ny*nproc_y
  integer :: tsini, tsbeg1, tsfin1, ts_re
  integer :: i, j, k, l, m
  integer :: nx0, nxf, ny0, nyf
  integer :: nproc, ts
  integer :: iwrite, nloop
  integer :: nxmin,nxmax,nymin,nymax,nzmin,nzmax
  real(4) :: atime, dtstep
  real(8) :: a8time, d8tstep
  real(4),dimension(0:nz,-1:ny,-1:nx) :: bx,by,bz,vx,vy,vz,ro,phi
  real(4),allocatable :: bx3d(:,:,:), by3d(:,:,:), bz3d(:,:,:), &
       & vx3d(:,:,:), vy3d(:,:,:), vz3d(:,:,:), ro3d(:,:,:), phi3d(:,:,:)
!  real(4),dimension(0:nx3d,0:ny3d,0:nz) :: bx3d, by3d, bz3d, &
!					& vx3d, vy3d, vz3d, &
!                                        & ro3d
  real(8),dimension(0:nz,-1:ny,-1:nx) :: bx8,by8,bz8,vx8,vy8,vz8,ro8,pr8,phi8
  character*100 :: dirname1="../../"
  character*100 :: dirname2="../../../data/newbz_20180725/sh0.1_cont/"
!  character*100 :: dirname_ini="../../../../../../make_shearedfield/data/relax_20171102/sh17_newbztf250/"
  character*100 :: dirname_ini=""
!  character*100 :: dirname_ini="../../../../../../makeshear_modified20180612/data/exec_dibvclean/sh1.1_0825relax/"
  

  character*100 :: fileb3d="B3D+"
  character*100 :: filere="RESTART."
  character(4) :: chr_process_number
  character(3) :: chr_timestep
  real(8),dimension(0:nx3d) :: xc
  real(8),dimension(0:ny3d) :: yc
  real(8),dimension(0:nz) :: zc
  real(8),dimension(-1:nx) :: xc0
  real(8),dimension(-1:ny) :: yc0
 
  integer :: forread(1:5), ixx

!--- get info ---
open(99,file="../code_output/info",status="old")
read(99,*) tsini; read(99,*) ts_re
read(99,*) nxmin; read(99,*) nxmax
read(99,*) nymin; read(99,*) nymax
read(99,*) nzmin; read(99,*) nzmax
read(99,'(a)') dirname_ini
close(99)
tsbeg1 = tsini + 1; tsfin1 = ts_re - 1
dirname_ini = trim(dirname_ini)
!+++++++++++
allocate(bx3d(0:nx3d,0:ny3d,0:nz)); allocate(by3d(0:nx3d,0:ny3d,0:nz))
allocate(bz3d(0:nx3d,0:ny3d,0:nz)); allocate(vx3d(0:nx3d,0:ny3d,0:nz))
allocate(vy3d(0:nx3d,0:ny3d,0:nz)); allocate(vz3d(0:nx3d,0:ny3d,0:nz))
allocate(ro3d(0:nx3d,0:ny3d,0:nz)); allocate(phi3d(0:nx3d,0:ny3d,0:nz))
!
ts = tsini
   call make_merge_file_8x8(ts,dirname_ini,filere,0)

do ts = tsbeg1, tsfin1
   call make_merge_file_8x8(ts,dirname1,fileb3d,1)
enddo
!
! do ts = tsbeg2, tsfin2
!    call make_merge_file_8x8(ts,dirname2,fileb3d)
! end do

ts = ts_re
   call make_merge_file_8x8(ts,dirname1,filere,1)
!-----------
!-----------
  !--- y dim
  ! open(40,file=dirname//"COORDINATE_y.0000", status="old")
  ! read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  ! read(40,'(e25.16)') yc0
  ! close(40)
  ! do i = 0, nx
  !    yc(i) = yc0(i)
  ! enddo
  ! open(40,file=dirname//"COORDINATE_y.0007", status="old")
  ! read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  ! read(40,'(e25.16)') yc0
  ! close(40)
  ! do i = 0, ny
  !    yc(i+ny) = yc0(i)
  ! enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0000", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0001", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0002", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*2) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0003", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*3) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0004", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*4) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0005", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*5) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0006", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*6) = yc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_y.0007", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') yc0
  close(40)
  do i = 0, ny
     yc(i+ny*7) = yc0(i)
  enddo

  !--- x dim
  open(40,file=trim(dirname1)//"COORDINATE_x.0000", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0008", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0016", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+2*nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0024", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+3*nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0032", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+4*nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0040", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+5*nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0048", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+6*nx) = xc0(i)
  enddo
  open(40,file=trim(dirname1)//"COORDINATE_x.0056", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') xc0
  close(40)
  do i = 0, nx
     xc(i+7*nx) = xc0(i)
  enddo

  !--- z dim
  open(40,file=trim(dirname1)//"COORDINATE_z.0000", status="old",convert="BIG_ENDIAN")
  read(40,*) forread(1), forread(2), forread(3), forread(4), forread(5)
  read(40,'(e25.16)') zc
  close(40)
  !--- output coordinate
  open(50,file="coord.xgc",status="replace")
  write(50,'(e25.16)') xc(0:nx3d)
  close(50)
  open(51,file="coord.ygc",status="replace")
  write(51,'(e25.16)') yc(0:ny3d)
  close(51)
  open(52,file="coord.zgc",status="replace")
  write(52,'(e25.16)') zc
  close(52)
  !
!
!  
!-----------
!  end do
!-----------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
subroutine make_merge_file_8x8(sub_ts,sub_dirname,sub_fname,check_endian)
integer,intent(in) :: sub_ts
character*100,intent(in) :: sub_dirname
character*100,intent(in) :: sub_fname
integer :: check_endian ! 0 --> LITTLE, 1 --> BIG

  print *, "file number:", sub_ts
!+++++++++++
  write(chr_timestep,'(i3.3)') sub_ts
  nproc = nproc_x * nproc_y
!-----------
  do i = 0, nproc-1
!-----------
  write(chr_process_number,'(i4.4)') i
  if(trim(sub_fname).eq."B3D+") then
  if(check_endian.eq.0) then
     open(10,file=trim(sub_dirname)//trim(sub_fname)//chr_timestep//"."//chr_process_number, form="unformatted", &
          & status="old")
  else if(check_endian.eq.1) then
     open(10,file=trim(sub_dirname)//trim(sub_fname)//chr_timestep//"."//chr_process_number, form="unformatted", &
          & status="old",convert="BIG_ENDIAN")
  else
     write(*,*) "ERROR CHECK_ENDIAN VALUE!!"
     stop
  end if
     read(10) iwrite, nloop, atime, dtstep, &
          & bx, by, bz, vx, vy, vz, ro, phi
     close(10)
  else if(trim(sub_fname).eq."RESTART.") then
  if(check_endian.eq.0) then
     open(10,file=trim(sub_dirname)//trim(sub_fname)//chr_process_number, form="unformatted", status="old")
  else if(check_endian.eq.1) then
     open(10,file=trim(sub_dirname)//trim(sub_fname)//chr_process_number, form="unformatted", status="old",convert="BIG_ENDIAN")
  else
     write(*,*) "ERROR CHECK_ENDIAN VALUE!!"
     stop
  end if
     read(10) iwrite, nloop, a8time, d8tstep, &
          & bx8, by8, bz8, vx8, vy8, vz8, ro8, pr8, phi8
     close(10)
     bx = bx8; by = by8; bz = bz8
     vx = vx8; vy = vy8; vz = vz8
     ro = ro8; phi = phi8
  end if
  
  nx0 = (i/8) * nx
  nxf = ((i/8) + 1) * nx - 1
  ny0 = mod(i,8) * ny
  nyf = (mod(i,8)+1) * ny - 1
  write(*,*) "debug", i, nx0, nxf, ny0, nyf

  do l = nx0, nxf
  do m = ny0, nyf
  bx3d(l,m,:) = bx(:,m-ny0,l-nx0)
  by3d(l,m,:) = by(:,m-ny0,l-nx0)
  bz3d(l,m,:) = bz(:,m-ny0,l-nx0)
  vx3d(l,m,:) = vx(:,m-ny0,l-nx0)
  vy3d(l,m,:) = vy(:,m-ny0,l-nx0)
  vz3d(l,m,:) = vz(:,m-ny0,l-nx0)
  ro3d(l,m,:) = ro(:,m-ny0,l-nx0)
  phi3d(l,m,:) = phi(:,m-ny0,l-nx0)
  enddo
!  write(*,*) l, nx0, nxf, m, ny0, nyf
  enddo
!-----------
  end do
!--debug
write(*,*) "max bz3d =", maxval(bz3d)
write(*,*) "min bz3d =", minval(bz3d)
!-----------
!  open(20,file=dirname//"Bxbin_3d", status="replace", form="unformatted")
!  open(21,file=dirname//"Bybin_3d", status="replace", form="unformatted")
!  open(22,file=dirname//"Bzbin_3d", status="replace", form="unformatted")
!  open(23,file=dirname//"Vxbin_3d", status="replace", form="unformatted")
!  open(24,file=dirname//"Vybin_3d", status="replace", form="unformatted")
!  open(25,file=dirname//"Vzbin_3d", status="replace", form="unformatted")
  open(30,file="Bxbin_3d."//chr_timestep, status="replace", form="unformatted")
  open(31,file="Bybin_3d."//chr_timestep, status="replace", form="unformatted")
  open(32,file="Bzbin_3d."//chr_timestep, status="replace", form="unformatted")
  open(33,file="Vxbin_3d."//chr_timestep, status="replace", form="unformatted")
  open(34,file="Vybin_3d."//chr_timestep, status="replace", form="unformatted")
  open(35,file="Vzbin_3d."//chr_timestep, status="replace", form="unformatted")
  open(36,file="Robin_3d."//chr_timestep, status="replace", form="unformatted")
  open(37,file="Phibin_3d."//chr_timestep, status="replace", form="unformatted")
  write(30) bx3d; write(31) by3d; write(32) bz3d
  write(33) vx3d; write(34) vy3d; write(35) vz3d
  write(36) ro3d; write(37) phi3d
  close(30); close(31); close(32)
  close(33); close(34); close(35)
  close(36); close(37)
!++++++++++++
  open(41,file="Bz2d."//chr_timestep,status="replace",form="unformatted")
  write(41) bz3d(:,:,0)
  close(41)
!++++++++++++  
end subroutine make_merge_file_8x8
!
end program main
