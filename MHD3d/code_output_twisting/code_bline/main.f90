program bline
      include "common.f90"
      integer :: i, j, k, numbl, isol=1, idir, line, dim_lines, line_each_dim
      integer :: rand_seed=12345, iout
      integer,parameter :: num_rand=3
	integer :: jnum

! === file name ===
character(100) :: cfile_sysout
character(100) :: cfile_run_number
character(100) :: cfile_coordinate_x
character(100) :: cfile_coordinate_y
character(100) :: cfile_coordinate_z
character(100) :: cfile_coordinate_z_all
character(100) :: cfile_eigenmode
character(100) :: cfile_output_list
character(100) :: cfile_time_list
character(100) :: cfile_nlff
character(100) :: cfile_nlff_2d
character(100) :: cfile_2d_field
character(100) :: cfile_3d_field
character(100) :: cfile_3d_init
character(100) :: cfile_initial
character(100) :: cfile_restart
character(100) :: cfile_slice_xz

namelist /nlist07/ cfile_sysout,       &
                   cfile_run_number,   &
                   cfile_coordinate_x, &
                   cfile_coordinate_y, &
                   cfile_coordinate_z, &
                   cfile_coordinate_z_all, &
                   cfile_eigenmode,    &
                   cfile_output_list,  &
                   cfile_time_list,    &
!                   cfile_nlff,         &
!                   cfile_nlff_2d,      &
                   cfile_2d_field,     &
                   cfile_3d_field,     &
                   cfile_3d_init,      &
!                   cfile_initial,      &
                   cfile_restart,      &
                   cfile_slice_xz

!
! isol = 1: Euler,  2: RKG
!
! filed line starting points are determined as
!  a point for dim_lines=0
!  equi-interval poins on a line for dim_lines = 1
!  grid points on a plane for dim_lines = 2
!  grid points in a cube for dim_lines = 3
      character clp*3,filen*11
      real(4),dimension(0:NX,0:NY,0:NZ) :: a3d
      real(8) :: r_start(3), x0, y0, z0
      real(8) :: dxl, dyl, dzl, del_x, del_y, del_z
      real(8) :: bx0, by0, bz0
      real(8) :: r_0(3),r_1(3), r_a(3),r_b(3),r_c(3)
      real(8) :: r_0a(3),r_0b(3),r_0c(3)
      real(8) :: alpha, beta, gamma
      real(8) :: del_alpha, del_beta, del_gamma
      real(8) :: randn(num_rand), dcell, bz_abs_max, edge
      integer :: nx_0,nx_1,ny_0,ny_1
!
! ++++ read Input ++++
!
      open(7,file='NAMELIST',status='OLD',err=99)
      read(7,nlist00)
      read(7,nlist07)
      close(7)
      write(*,*) '###DEBUG 1'
!
! +++ open Bline file +++
!
!      open(20,file='Blines',status='UNKNOWN')
!
!
! ++++ coordinate data input ++++
!
!      open(8,file='../binfile/coord.xgc',status='OLD')
      open(8,file='../3dbin/coord.xgc',status='OLD')
        do i = 0, NX
          read(8,*) xc0(i)
        end do
      close(8)
!
!      open(8,file='../binfile/coord.ygc',status='OLD')
      open(8,file='../3dbin/coord.ygc',status='OLD')
        do j = 0, NY
         read(8,*) yc0(j)
        end do
      close(8)
!
!      open(8,file='../binfile/coord.zgc',status='OLD')
      open(8,file='../3dbin/coord.zgc',status='OLD')
        do k = 0, NZ
         read(8,*) zc0(k)
        end do
!
      do i = 0, nxl
         xc(i) = xc0(nx0+i*leapx)
      end do
!
      do j = 0, nyl
         yc(j) = yc0(ny0+j*leapy)
      end do
!
      do k = 0, nzl
         zc(k) = zc0(nz0+k*leapz)
      end do
!
! >>>>>>>>> read field data
!
!      open(50,file='readN',status='OLD')
!      read(50,*) j
!      close(50)
!	do jnum = 1, 10      !----- for debug or test
	do jnum = 49, 82
      write(clp,'(I3.3)') jnum
        write(*,*) "-------------" 
	write(*,*) "clp =", clp
!
! .................read 3d field
!     
!      open(10,file=trim(cfile_3d_field)//'_R.'//clp//'.BX',form='UNFORMATTED',status='OLD')
!      open(10,file="../binfile/Bxbin_3d."//clp,form='UNFORMATTED',status='OLD')
      open(10,file="../3dbin/Bxbin_3d_R."//clp,form='UNFORMATTED',status='OLD')
      read(10) a3d
      close(10)
      call reduced_array(a3d,bx)

!      open(10,file=trim(cfile_3d_field)//'_R.'//clp//'.BY',form='UNFORMATTED',status='OLD')
!      open(10,file="../binfile/Bybin_3d."//clp,form='UNFORMATTED',status='OLD')
      open(10,file="../3dbin/Bybin_3d_R."//clp,form='UNFORMATTED',status='OLD')
      read(10) a3d
      close(10)
      call reduced_array(a3d,by)

!      open(10,file=trim(cfile_3d_field)//'_R.'//clp//'.BZ',form='UNFORMATTED',status='OLD')
!      open(10,file="../binfile/Bzbin_3d."//clp,form='UNFORMATTED',status='OLD')
      open(10,file="../3dbin/Bzbin_3d_R."//clp,form='UNFORMATTED',status='OLD')
      read(10) a3d
      close(10)
 write(*,*) '#DEBUG BZ min max ',minval(a3d),maxval(a3d)
      call reduced_array(a3d,bz)

!-------------------- prepare rot B
	call calc_rotb_grid
	write(*,*) "output file ROT_B"
	open(28,file="file_output/Jx3dbin."//clp,form="unformatted",status="UNKNOWN")
	write(28) rotbx
	close(28)
	open(29,file="file_output/Jy3dbin."//clp,form="unformatted",status="UNKNOWN")
	write(29) rotby
	close(29)
	open(30,file="file_output/Jz3dbin."//clp,form="unformatted",status="UNKNOWN")
	write(30) rotbz
        close(30)
 	open(36,file="file_output/JJ3dbin."//clp,form="unformatted",status="UNKNOWN")
	write(36) sqrt(rotbx**2+rotby**2+rotbz**2)
        close(36)
write(*,*) "bx(60,60,0) =", bx(60,60,0)
write(*,*) "by(60,60,0) =", by(60,60,0)
write(*,*) "bz(60,60,0) =", bz(60,60,0)
write(*,*) "max bz =", maxval(bz)
write(*,*) "min bz =", minval(bz)
write(*,*) "rotbx(60,60,0) =", rotbx(60,60,0)
write(*,*) "rotby(60,60,0) =", rotby(60,60,0)
write(*,*) "rotbz(60,60,0) =", rotbz(60,60,0)
!write(*,*) "min rotbz on z=0 =", minval(rotbz(:,:,0))
!write(*,*) "max rotbz on z=0 =", maxval(rotbz(:,:,0))
!write(*,*) "min rotbz location on z=0", maxloc(rotbx(:,:,0),2)
!
	open(31,file="file_output/alp2dbin."//clp,form="unformatted",status="UNKNOWN")
	write(31) alpxy
	close(31)
!
! 	open(32,file="file_output/alpz1_2dbin."//clp,form="unformatted",status="UNKNOWN")
! 	write(32) alpxy1
! 	close(32)
! !
! 	open(33,file="file_output/alpz2_2dbin."//clp,form="unformatted",status="UNKNOWN")
! 	write(33) alpxy2
! 	close(33)
write(*,*) "alpha_xy(60,60,0) =", alpxy(60,60)
write(*,*) "min aplha_xy =", minval(alpxy)
write(*,*) "max aplha_xy =", maxval(alpxy)
!write(*,*) "min rotbz location on z=0", minloc(alpxy,2)
!
write(*,*) "alpha_xy(60,60,1) =", alpxy1(60,60)
write(*,*) "min aplha_xy on z=1 =", minval(alpxy1)
write(*,*) "max aplha_xy on z=1 =", maxval(alpxy1)
!write(*,*) "min rotbz location on z=1", minloc(alpxy1,2)
!-------------------- open twist file on yz-plane
!	open(40,file="file_output/twist_yzplane."//clp,status='UNKNOWN')
!-------------------- open init foot point and final foot point file for posi(r) and nega(l)
!	open(50,file="file_output/Bfootyz_right"//clp,status="UNKNOWN")
!	open(55,file="file_output/Bfootyz_left"//clp,status="UNKNOWN")
!--------- open init and final foot point file for idir
open(70,file="file_output/Bfootp1"//clp,status="UNKNOWN")
open(71,file="file_output/Bfootm1"//clp,status="UNKNOWN")
!--------- open file init and yz-plane point information
open(60,file="file_output/Byzp"//clp,status="UNKNOWN")
! ...................Initial Point

      open(10,file='Param_bline',status='OLD')
 500  read(10,*,end=999,err=999) dim_lines,numbl,rbl,gbl,bbl,wbl
      write(*,*) '## Param_bline :',dim_lines,numbl,rbl,gbl,bbl,wbl
      if(dim_lines.lt.0) go to 999
!
! >>>>>>>> filed line tracing
! 
       call preprg
!
!## POINT ##
  if(dim_lines.eq.0) then

      read(10,*,end=999) r_start(1:3)

        idir = 1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
        idir = -1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)

!## LINE ##         
   else if(dim_lines.eq.1) then

      read(10,*,end=999) r_0(1:3), r_0a(1:3)

      do line = 1, numbl

         r_start(:) = r_0(:) + (line-1)*r_0a(:)

        idir = 1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
        idir = -1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)

      end do

!## PLANE ##         
   else if(dim_lines.eq.2) then

      read(10,*,end=999) r_0(1:3), r_1(1:3)
      write(*,*) r_0(1:3), r_1(1:3)
!
      r_0(1) = 0.5d0; r_0(2) = 1d0; r_0(3) = 1d0
      r_1(1) = 0.5d0; r_1(2) = 0d0; r_0(3) = 0d0

!      line_each_dim = max(2, int(real(numbl)**0.5))

!	line_each_dim = 1025
!      del_x =(r_1(1)-r_0(1))/(line_each_dim-1)
!	line_each_dim = 257
!      del_y =(r_1(2)-r_0(2))/(line_each_dim-1)
!	line_each_dim = 513
!      del_z =(r_1(3)-r_0(3))/(line_each_dim-1)
!=====
!      do i = 1, line_each_dim+512
!	print *, "NOW calbl step i =", i
!      do j = 1, line_each_dim-256
!
!         r_start(1) = r_0(1) + del_x*(i-1) 
!         r_start(2) = r_0(2) + del_y*(j-1) 
!         r_start(3) = r_0(3) + del_z*(i-1)
!
!         write(*,*) '*r_start ',r_start
!
!        idir = 1
!        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
!        idir = -1
!        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
!
!      end do
!      end do
!=====
	do i = 0, NX
!	  if(mod(i,64)==0) print *, "===== NOW calbl x step i =", i
	do j = 0, NY
	  r_start(1) = xc(i)
	  r_start(2) = yc(j)
	  r_start(3) = zc(o)
!	  write(*,*) "r_start ", r_start
	  idir = 1
	  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,i,j)
	  idir = -1
	  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,i,j)
	enddo
!if(mod(i,NX/8).eq.0) write(*,*) "NOW i =", i, xc(i)
	enddo
!----
	print *, "output file twist on 'z=0 plane'", clp
	open(32,file="file_output/Tw2dbin."//clp,form="unformatted",status="UNKNOWN")
	write(32) tw_xy
	close(32)
write(*,*) "min twist_xy =", minval(tw_xy)
write(*,*) "max twist_xy =", maxval(tw_xy)
!--- initializing twist_xy
tw_xy=0d0
!----
!## CUBE ##         
   else if(dim_lines.eq.3) then

      read(10,*,end=999) r_0(1:3), r_0a(1:3), r_0b(1:3), r_0c(1:3)

      line_each_dim = max(2, int(real(numbl)**(1.0d0/3.0d0)))

!      r_0a(:) = r_a(:)-r_0(:)
!      r_0b(:) = r_b(:)-r_0(:)
!      r_0c(:) = r_c(:)-r_0(:)

!      del_alpha = 1.0d0/(line_each_dim-1)
!      del_beta  = 1.0d0/(line_each_dim-1)
!      del_gamma = 1.0d0/(line_each_dim-1)
       del_alpha = 1.0d0
       del_beta  = 1.0d0
       del_gamma = 1.0d0

      do i = 1, line_each_dim
         alpha = del_alpha*(i-1)

      do j = 1, line_each_dim
         beta  = del_beta*(j-1)

      do k = 1, line_each_dim
         gamma  = del_gamma*(k-1)

         r_start(:) = r_0(:) + alpha*r_0a(:) + beta*r_0b(:) + gamma*r_0c(:)

        idir = 1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
        idir = -1
        call calbl(r_start(1),r_start(2),r_start(3),isol,idir)

      end do
      end do
      end do

!## PLANE at random ##         
   else if(dim_lines.eq.10) then
   read(10,*,end=999,err=999) nx_0,nx_1,ny_0,ny_1
   write(*,*) nx_0,nx_1,ny_0,ny_1

      bz_abs_max = maxval(abs(bz(nx_0:nx_1,ny_0:ny_1,0)))
      line = 0
      do while (line<numbl)
         
         call urand(num_rand,randn,rand_seed)

         edge=0.00
         x0 = randn(1)*(xc(nx_1)-xc(nx_0))*(1.0d0-edge*2) &
                      +(xc(nx_0)+edge)
         y0 = randn(2)*(yc(ny_1)-yc(ny_0))*(1.0d0-edge*2) &
                      +(yc(ny_0)+edge)
         z0 = 2.0d-3
         call sint(x0,y0,z0,bx0,by0,bz0,dcell,iout)
!write(*,*) line,x0,y0,abs(bz0)/bz_abs_max,randn(3) !DEBUG

         bz0_norm = abs(bz0)/bz_abs_max
         if(randn(3).le.bz0_norm) then
           r_start(1) = x0
           r_start(2) = y0
           r_start(3) = z0
!           write(*,*) i,j,r_start(:),alpha,beta !DEBUG

          if(bz0.gt.0.0) then
             idir = 1
             call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
           else
             idir = -1
             call calbl(r_start(1),r_start(2),r_start(3),isol,idir)
          end if

           line = line + 1

         end if

      end do

   end if

      go to 500

 999  continue
!
      close(10)
!      close(20)
	close(40)
	close(50)
	close(55)
        close(70)
        close(71)
	enddo
      stop
 99   write(6,*) ' Input ERR '
!      close(20)
	close(40)
	close(50)
	close(55)
        close(70)
        close(71)
      stop
contains
!
subroutine reduced_array(a,b)
real(4),dimension(0:NX,0:NY,0:NZ) :: a
real(4),dimension(0:nxl,0:nyl,0:nzl) :: b
integer :: i, j, k
!
do k = 0, nzl
do j = 0, nyl
do i = 0, nxl
   b(i,j,k) = a(nx0+i*leapx,ny0+j*leapy,nz0+k*leapz)
end do
end do
end do
!
end subroutine reduced_array
subroutine urand(n,x,ir)
!************************************************************************
!* UNIFORM RANDOM NUMBER GENERATOR (MIXED CONGRUENTIAL METHOD)          *
!*     PORTABLE BUT SLOW.  THE PERIOD IS ONLY 1664501.                  *
!* PARAMETERS                                                           *
!*   (1) N      (I) THE NUMBER OF RANDOM NUMBERS TO BE GENERATED        *
!*                  (INPUT)                                             *
!*   (2) X      (D) UNIFORM RANDOM NUMBERS (OUTPUT)                     *
!*   (3) IR     (I) THE INITIAL SEED  (INPUT)                           *
!*                  THE SEED FOR THE NEXT CALL (OUTPUT)                 *
!* COPYRIGHT: Y. OYANAGI, JUNE 30, 1989  V.1                            *
!************************************************************************
!*
       DOUBLE PRECISION X(N), INVM
       PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
       PARAMETER (INVM = 1.0D0 / M)
!*PAREMETER CHECK
      IF( N .LE. 0) THEN
       WRITE(6,*) '(SUBR.URAND1) PARAMETER ERROR. N = ', N
       WRITE(6,*) 'RETURN WITH NO FURTHER CALCULATION.'
       RETURN
      END IF
      IF( IR .LT. 0 .OR. IR .GE. M) THEN
       WRITE(6,*) '(SUBR.URAND1) WARNING. IR = ', IR
      END IF
!*MAIN LOOP
      DO 10 I = 1, N
       IR = MOD( LAMBDA * IR + MU, M)
       X(I) = IR * INVM
   10 CONTINUE
      RETURN
 end subroutine urand

end program bline
