! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
!  Files:
!      RUN_NUMBER: if not exist, the calculation starts from the
!                  initial condition.
!                  if exits, read the last run_number from this file.
!  2005/06/09 new subroutines pset__zero_set, cp_double_to_d3d
!  2005/06/09 ddx2,ddy2,ddz2 multiplyed by two
module pset
! -----------------------------------------------------------
!                    MODULE pset
! -----------------------------------------------------------
use common
use mpiut
use job

implicit none

integer :: index_x_check, index_y_check

public :: pset__zero_set, &
          pset__init, &
          pset__resistivity, &
          pset__dt, &
          pset__namelist

private :: pset__x_grid, &
           pset__y_grid, &
           pset__z_grid, &
           cp_double_to_d3d

contains

! -----------------------------------------------------------
subroutine pset__zero_set

call cp_double_to_d3d(0.0d0,vx)
call cp_double_to_d3d(0.0d0,vy)
call cp_double_to_d3d(0.0d0,vz)
call cp_double_to_d3d(0.0d0,bx)
call cp_double_to_d3d(0.0d0,by)
call cp_double_to_d3d(0.0d0,bz)
call cp_double_to_d3d(0.0d0,ro)
call cp_double_to_d3d(0.0d0,pr)

call cp_double_to_d3d(0.0d0,dvx)
call cp_double_to_d3d(0.0d0,dvy)
call cp_double_to_d3d(0.0d0,dvz)
call cp_double_to_d3d(0.0d0,dbx)
call cp_double_to_d3d(0.0d0,dby)
call cp_double_to_d3d(0.0d0,dbz)
call cp_double_to_d3d(0.0d0,dro)
call cp_double_to_d3d(0.0d0,dpr)

call cp_double_to_d3d(0.0d0,qvx)
call cp_double_to_d3d(0.0d0,qvy)
call cp_double_to_d3d(0.0d0,qvz)
call cp_double_to_d3d(0.0d0,qbx)
call cp_double_to_d3d(0.0d0,qby)
call cp_double_to_d3d(0.0d0,qbz)
call cp_double_to_d3d(0.0d0,qro)
call cp_double_to_d3d(0.0d0,qpr)

call cp_double_to_d3d(0.0d0,ex)
call cp_double_to_d3d(0.0d0,ey)
call cp_double_to_d3d(0.0d0,ez)
call cp_double_to_d3d(0.0d0,cx)
call cp_double_to_d3d(0.0d0,cy)
call cp_double_to_d3d(0.0d0,cz)

call cp_double_to_d3d(0.0d0,rvx)
call cp_double_to_d3d(0.0d0,rvy)
call cp_double_to_d3d(0.0d0,rvz)
call cp_double_to_d3d(0.0d0,qx)
call cp_double_to_d3d(0.0d0,qy)
call cp_double_to_d3d(0.0d0,qz)

call cp_double_to_d3d(0.0d0,tmp)
call cp_double_to_d3d(0.0d0,res)
call cp_double_to_d3d(0.0d0,visc3d)
call cp_double_to_d3d(0.0d0,phi)

end subroutine pset__zero_set
! -----------------------------------------------------------
subroutine pset__init
integer :: isend, ierr
integer :: i,j,k

! ------------------------------<< sysout file >>
  open(FILE_SYSOUT,file=trim(cfile_sysout)//cmyrank, &
       form='formatted')
! ------------------------------<< run number >>
  
  if(myrank == root) then

    open(FILE_RUN_NUMBER,file=trim(cfile_run_number),form='formatted')
       run_number = 0
          write(*,*) '::>> pset__init: run_number(before read)',run_number
       do
          write(*,*) '::>> pset__init: run_number(before read 2)',run_number
          read(FILE_RUN_NUMBER,*,end=101) run_number
          write(*,*) '::>> pset__init: run_number(read)',run_number
       end do
101    run_number = run_number + 1
!
       write(FILE_RUN_NUMBER,*) run_number   ! wrtie the current run_number 
    close(FILE_RUN_NUMBER)

  end if

! ------------------------------<<send run_number to other nodes>>

    isend = run_number
    call mpi_bcast(isend,1,MPI_INTEGER, &
                   root,MPI_COMM_CART,ierr)
    run_number = isend

!========= DEBUG
    write(*,*) '::>> pset__init: run_number',run_number
!=========

    crun_number = '_'//chari2(run_number)

! ================ start from initial ==================
    if(run_number == 1) then
         start_from_initial = .true.
     else
         start_from_initial = .false.
    end if
! ======================================================

! ------------------------------<< open coordinate file >>

 open(FILE_COORDINATE_X,file=trim(cfile_coordinate_x)//cmyrank, &
      form='formatted')
 open(FILE_COORDINATE_Y,file=trim(cfile_coordinate_y)//cmyrank, &
      form='formatted')
 open(FILE_COORDINATE_Z,file=trim(cfile_coordinate_z)//cmyrank, &
      form='formatted')

! ------------------------------<< output file >>

if(myrank == root) then
 open(FILE_TIME_LIST,   file=trim(cfile_time_list)//crun_number, &
      form='formatted')
end if

 open(FILE_OUTPUT_LIST, file=trim(cfile_output_list)//crun_number//cmyrank, &
      form='formatted')

 write(FILE_OUTPUT_LIST,*) 
 write(FILE_OUTPUT_LIST,*) ' ++ 3D zero-beta nonlinear MHD MODEL ++'
 write(FILE_OUTPUT_LIST,*) ' NX, NY, NZ = ',nx,ny,nz
 write(FILE_OUTPUT_LIST,*) 

 write(FILE_OUTPUT_LIST,nml=nlist00)
 write(FILE_OUTPUT_LIST,nml=nlist01)
 write(FILE_OUTPUT_LIST,nml=nlist02)
 write(FILE_OUTPUT_LIST,nml=nlist03)
 write(FILE_OUTPUT_LIST,nml=nlist04)
 write(FILE_OUTPUT_LIST,nml=nlist05)
 write(FILE_OUTPUT_LIST,nml=nlist06)
!! write(FILE_OUTPUT_LIST,nml=nlist06a) 
 write(FILE_OUTPUT_LIST,nml=nlist07)
 write(FILE_OUTPUT_LIST,nml=nlist08)
 write(FILE_OUTPUT_LIST,nml=nlist_ef)
 write(FILE_OUTPUT_LIST,nml=nlist_lfff)

! -------------------------------<< set parameters for lfff >>
  lfff_k  = lfff_kly_pi/yl*PI
  lfff_kk = sqrt(lfff_k**2-lfff_alpha**2)

! -------------------------------<< set grid coordinate >>

  call pset__x_grid
  call pset__y_grid
  call pset__z_grid
!
! .... coefficient at the center for the laplacian operator

  do i = -1, NX
  do j = -1, NY
  do k =  0, NZ
     d2xyz0(k,j,i) = d2x(i,0)+d2y(j,0)+d2z(k,0)
  end do
  end do
  end do

! -------------------------------<< set physical parameter >>

  gamma1 = gamma - 1.0d0

end subroutine pset__init

! -----------------------------------------------------------
subroutine pset__namelist
!
! ++++ input parameter ++++
!
open(FILE_NAMELIST,file='NAMELIST')

      read(FILE_NAMELIST,nml=nlist00)
      read(FILE_NAMELIST,nml=nlist01)
      read(FILE_NAMELIST,nml=nlist02)
      read(FILE_NAMELIST,nml=nlist03)
      read(FILE_NAMELIST,nml=nlist04)
      read(FILE_NAMELIST,nml=nlist05)
      read(FILE_NAMELIST,nml=nlist06)
!!      read(FILE_NAMELIST,nml=nlist06a)  !! for the forcing to converge
      read(FILE_NAMELIST,nml=nlist07)
      read(FILE_NAMELIST,nml=nlist08)
      read(FILE_NAMELIST,nml=nlist_ef)
      read(FILE_NAMELIST,nml=nlist_lfff)
      read(FILE_NAMELIST,nml=nlist_divbeq0)

close(FILE_NAMELIST)

! ==============================WARNING============================
if(nloop_output .gt. nloop_incmax) then
   write(FILE_SYSOUT) ':: WARNING : nloop_output > nloop_incmax ', &
   nloop_output, nloop_incmax
end if

end subroutine pset__namelist

! -----------------------------------------------------------
subroutine pset__x_grid

real(DP) :: delx, xi0, dxmin, dxmax, dxchg, bunbo
real(DP) :: xleft ! uniform x coordinate at the left end
real(DP) :: xright ! uniform x coordinate at the right end
real(DP) :: xlength ! uniform x length included in this rank
integer :: i
real(DP) :: para_x0, para_ix0
real(DP) :: nbeki
real(DP) :: alphax, betax
para_x0 = 8.0 
para_ix0 = 224.0
nbeki = 11
!
alphax = (para_x0*(real(nx3d)*0.5d0)**nbeki &
     &  - real(xl)*0.5d0*para_ix0**nbeki) &
     &  /(para_ix0*(real(nx3d)*0.5d0) &
     & *((real(nx3d)*0.5d0)**(nbeki-1) &
     &  -(para_ix0**(nbeki-1))))
betax =  (real(xl)*0.5d0*para_ix0 &
     &  - real(nx3d)*0.5d0*para_x0) &
     &  /(para_ix0*(real(nx3d)*0.5d0) &
     & *((real(nx3d)*0.5d0)**(nbeki-1) &
     &  -(para_ix0**(nbeki-1))))
!

   xlength = xl / nproc_x
    
       delx = xlength / nx
       xleft = xlength * index_x
       xright= xleft + xlength
!#### debug ####
! write(*,*) '##pset__x_grid index_x, xleft,xlength', index_x, xleft,xlength
!###############

      do i = -1, nx
         xi0 = xleft + delx*i
!         xc(i) = xi0 + epsx*xl/(PI2)*sin(PI2*xi0/xl)
!         xc(i) = xi0 + epsx/(PI2)*sin(PI2*xi0/xl)
!         xc(i) = xi0 * 0.03d0 &
!              & + ((4.0*xi0)**11.0)*62.08/(0.5d0*nx)**11.0 &
!         xc(i) = (xi0 - xl*0.5d0)*0.03d0 &
!              & + ((4.0*(xi0 - xl*0.5d0))**11.0)*62.08/(0.5d0*nx3d)**11.0 &
!              & + xl*0.5d0
         xc(i) = alphax*(i+nx*index_x-nx3d/2) + betax*(i+nx*index_x-nx3d/2)**nbeki
!         xc(i) = xi0
      end do
!
      do i =-1, nxm1
         dx(i) = xc(i+1) - xc(i)
        ddx(i) = 1.0d0/dx(i)
        if(i.eq.-1) then
           dxmin = dx(i)
           dxmax = dx(i)
         else
           dxmin = min(dx(i), dxmin)
           dxmax = max(dx(i), dxmax)
        end if
        if(i.eq.0) then
           dxchg = abs(dx(i)-dx(i-1))/dx(i)
         else if(i.gt.0) then
           dxchg = max(dxchg, abs(dx(i)-dx(i-1))/dx(i))
        end if
      end do
!
      dx(nx) = dx(nxm1)
      ddx(nx)= ddx(nxm1)
!
      do i = 0, nx
         dx2(i) = dx(i-1) + dx(i)
        ddx2(i) = 2.0d0/dx2(i)
!
         bunbo = 1.0d0/(dx(i-1)*dx(i)*dx2(i))
         d1x(i,-1) = -dx(i)**2*bunbo
         d1x(i, 0) = (dx(i)**2-dx(i-1)**2)*bunbo
         d1x(i, 1) =  dx(i-1)**2*bunbo
!
         d2x(i,-1) = 2*dx (i  )*bunbo
         d2x(i, 0) =-2*dx2(i  )*bunbo
         d2x(i, 1) = 2*dx (i-1)*bunbo

      end do
!
   if(index_x == 0) then
      i = 0
!      i = -1
!
         dx2(i) = dx(i)+dx(i+1)
!        ddx2(i) = 2.0d0/dx2(i)

         bunbo = 1./(dx(i)*dx(i+1)*dx2(i+1))
         d1x_0(-1) = -dx (i+1)*(dx(i+1)+2*dx(i))*bunbo
         d1x_0( 0) =  dx2(i+1)**2*bunbo
         d1x_0( 1) = -dx (i)**2*bunbo
!
!         d2x_0(-1) = 2*dx (i+1)*bunbo
!         d2x_0( 0) =-2*dx2(i+1)*bunbo
!         d2x_0( 1) = 2*dx (i  )*bunbo
   end if
!
   if(index_x == nproc_x-1) then
      i = nx
!
         dx2(i) = dx(i-1)+dx(i)
!        ddx2(i) = 2.0d0/dx2(i)

         bunbo = 1./(dx(i-1)*dx(i-2)*dx2(i-1))
         d1x_nx(-1) =  dx (i-1)**2*bunbo
         d1x_nx( 0) = -dx2(i-1)**2*bunbo
         d1x_nx( 1) =  dx (i-2)*(dx(i-2)+2*dx(i-1))*bunbo
!
!         d2x_nx(-1) = 2*dx (i-1)*bunbo
!         d2x_nx( 0) =-2*dx2(i-1)*bunbo
!         d2x_nx( 1) = 2*dx (i-2)*bunbo         
    end if

! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_x_grid : ',dxmin, dxmax, dxchg
! -----------------------------------------------------------

! ------------- write coordinate file X ---------------------
!  if(run_number == 1) then
      write(FILE_COORDINATE_X,*) index_x, index_y, nx, ny, nz
      write(FILE_COORDINATE_X,'(e25.16)') xc
!  else
!     read(FILE_COORDINATE_X,*) index_x_check, index_y_check
!     call job__check_equal(index_x, index_x_check, '**X**')
!     call job__check_equal(index_y, index_y_check, '**Y**')
!  end if 
! -----------------------------------------------------------
end subroutine pset__x_grid

! -----------------------------------------------------------
subroutine pset__y_grid

implicit none
real(DP) :: dely, yj0, dymin, dymax, dychg, bunbo
real(DP) :: ybottom ! uniform y coordinate at the bottom
real(DP) :: ytop ! uniform y coordinate at the top
real(DP) :: ylength ! uniform y length included in this rank
integer :: j
real(DP) :: para_y0, para_jy0
real(DP) :: nbeki
real(DP) :: alphay, betay
para_y0 = 8.0
para_jy0 = 224.0
nbeki = 11
!
alphay = (para_y0*(real(ny3d)*0.5d0)**nbeki &
     &  - real(yl)*0.5d0*para_jy0**nbeki) &
     &  /(para_jy0*(real(ny3d)*0.5d0) &
     & *((real(ny3d)*0.5d0)**(nbeki-1) &
     &  -(para_jy0**(nbeki-1))))
betay =  (real(yl)*0.5d0*para_jy0 &
     &  - real(ny3d)*0.5d0*para_y0) &
     &  /(para_jy0*(real(ny3d)*0.5d0) &
     & *((real(ny3d)*0.5d0)**(nbeki-1) &
     &  -(para_jy0**(nbeki-1))))
!
   ylength = yl / nproc_y

      dely = ylength / ny
      ybottom = ylength * index_y
      ytop    = ybottom + ylength
      do j = -1, ny
         yj0 = ybottom + dely*j
!         yc(j) =yj0 + epsy*yl/(PI2)*sin(PI2*yj0/yl)
!         yc(j) =yj0 + epsy/(PI2)*sin(PI2*yj0/yl)
!         yc(j) = (yj0 - yl*0.5d0) * 0.01d0 &
!              & + ((4.0*(yj0 - yl*0.5d0))**11.0)*63.36/(0.5*ny3d)**11.0 &
!              & + yl * 0.5d0
         yc(j) = alphay*(j+ny*index_y-ny3d/2) + betay*(j+ny*index_y-ny3d/2)**nbeki
!         yc(j) = yj0
      end do
!
      do j =-1, nym1
         dy(j) = yc(j+1) - yc(j)
        ddy(j) = 1.0d0/dy(j)
        if(j.eq.-1) then
           dymin = dy(j)
           dymax = dy(j)
         else
           dymin = min(dy(j), dymin)
           dymax = max(dy(j), dymax)
        end if
        if(j.eq.0) then
           dychg = abs(dy(j)-dy(j-1))/dy(j)
         else if(j.gt.0) then
           dychg = max(dychg, abs(dy(j)-dy(j-1))/dy(j))
        end if
       end do
!
      dy(ny) = dy(nym1)
      ddy(ny)= ddy(nym1)
!
      do j = 0, nym1
         dy2(j) = dy(j-1) + dy(j)
        ddy2(j) = 2.0d0/dy2(j)
!
         bunbo = 1.0d0/(dy(j-1)*dy(j)*dy2(j))
         d1y(j,-1) = -dy(j)**2*bunbo
         d1y(j, 0) = (dy(j)**2-dy(j-1)**2)*bunbo
         d1y(j, 1) =  dy(j-1)**2*bunbo
!
         d2y(j,-1) = 2*dy (j  )*bunbo
         d2y(j, 0) =-2*dy2(j  )*bunbo
         d2y(j, 1) = 2*dy (j-1)*bunbo
      end do
!
   if(index_y == 0) then
      j = 0
!      j = -1
!
         dy2(j) = dy(j)+dy(j+1)
        ddy2(j) = 2.0d0/dy2(j)

         bunbo = 1./(dy(j)*dy(j+1)*dy2(j+1))
         d1y_0(-1) = -dy (j+1)*(dy(j+1)+2*dy(j))*bunbo
         d1y_0( 0) =  dy2(j+1)**2*bunbo
         d1y_0( 1) = -dy (j)**2*bunbo
!
         d2y_0(-1) = 2*dy (j+1)*bunbo
         d2y_0( 0) =-2*dy2(j+1)*bunbo
         d2y_0( 1) = 2*dy (j  )*bunbo
   end if
!
   if(index_y == nproc_y-1) then
      j = ny
!
         dy2(j) = dy(j-1)+dy(j)
        ddy2(j) = 2.0d0/dy2(j)

         bunbo = 1./(dy(j-1)*dy(j-2)*dy2(j-1))
         d1y_ny(-1) =  dy (j-1)**2*bunbo
         d1y_ny( 0) = -dy2(j-1)**2*bunbo
         d1y_ny( 1) =  dy (j-2)*(dy(j-2)+2*dy(j-1))*bunbo
!
         d2y_ny(-1) = 2*dy (j-1)*bunbo
         d2y_ny( 0) =-2*dy2(j-1)*bunbo
         d2y_ny( 1) = 2*dy (j-2)*bunbo         
    end if
! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_y_grid : ',dymin, dymax, dychg
! -----------------------------------------------------------

! ------------- write coordinate file Y ---------------------
!  if(run_number == 1) then
      write(FILE_COORDINATE_Y,*) index_x, index_y, nx, ny, nz
      write(FILE_COORDINATE_Y,'(e25.16)') yc
!  else
!     read(FILE_COORDINATE_Y,*) index_x_check, index_y_check
!     call job__check_equal(index_x, index_x_check, '**X**')
!     call job__check_equal(index_y, index_y_check, '**y**')
!  end if 
! -----------------------------------------------------------
end subroutine pset__y_grid

! -----------------------------------------------------------
subroutine pset__z_grid

implicit none
real(DP) :: fraction_fine_grid, &
            z1, delz, zca, zcb, zcc, zcd, &
            zk, dzmin, dzmax, dzchg, &
            bunbo
real(DP) :: zk0
integer :: k
integer :: SWITCH_COORDZ_FILE = 0  ! 0: zc(:) unform grid
                                   ! 1: zc(:) determined by the file(coord.zgc)
real(DP) :: para_z0, para_kz0
real(DP) :: nbeki
real(DP) :: alphaz, betaz
para_z0 = 8.0
para_kz0 = 205.0
nbeki = 11
!
alphaz = (para_z0*real(nz3d)**nbeki &
     &  - real(zl)*para_kz0**nbeki) &
     &  /(para_kz0*real(nz3d) &
     & *(real(nz3d)**(nbeki-1) &
     &  -(para_kz0**(nbeki-1))))
betaz =  (real(zl)*para_kz0 &
     &  - real(nz3d)*para_z0) &
     &  /(para_kz0*real(nz3d) &
     & *(real(nz3d)**(nbeki-1) &
     &  -(para_kz0**(nbeki-1))))
!

!------set coordinate array---------
!! read coord.zgc 
if(SWITCH_COORDZ_FILE.eq.1) then
    open(FILE_COORDINATE_Z_ALL, &
         file=trim(cfile_coordinate_z_all), &
         form='formatted')
    read(FILE_COORDINATE_Z_ALL,'(e25.16)') zc
    close(FILE_COORDINATE_Z_ALL)
else
      delz = zl/NZ
      do k = 0, NZ
!         zc(k) = delz*k
         zk0 = delz*k
!         zc(k) = zk0 + epsz/(PI2)*sin(PI*zk0/zl+PI)
         zc(k) = alphaz*k + betaz*k**nbeki
      end do
end if
!-----------------------------------
!
      do k = 0, nzm1
         dz(k) = zc(k+1) - zc(k)
        ddz(k) = 1.0d0/dz(k)
        if(k.eq.0) then
           dzmin = dz(k)
           dzmax = dz(k)
         else
           dzmin = min(dz(k), dzmin)
           dzmax = max(dz(k), dzmax)
        end if
        if(k.eq.1) then
           dzchg = abs(dz(k)-dz(k-1))/dz(k)
         else if(k.gt.1) then
           dzchg = max(dzchg, abs(dz(k)-dz(k-1))/dz(k))
        end if
      end do

      dz(nz) = dz(nzm1)
      ddz(nz) = ddz(nzm1)
!
      do k = 1, nzm1
         dz2(k) = dz(k-1) + dz(k)
        ddz2(k) = 2.0d0/dz2(k)
!
         bunbo = 1./(dz(k)*dz(k-1)*dz2(k))
         d1z(k,-1) = -dz(k)**2*bunbo
         d1z(k, 0) = (dz(k)**2-dz(k-1)**2)*bunbo
         d1z(k, 1) =  dz(k-1)**2*bunbo
!
         d2z(k,-1) = 2*dz (k  )*bunbo
         d2z(k, 0) =-2*dz2(k  )*bunbo
         d2z(k, 1) = 2*dz (k-1)*bunbo
      end do
!
!
      k = 0
!
         dz2(k) = dz(k)
        ddz2(k) = 2.0d0/dz2(k)

         bunbo = 1./(dz(k)*dz(k+1)*dz2(k+1))
         d1z(k,-1) = -dz (k+1)*(dz(k+1)+2*dz(k))*bunbo
         d1z(k, 0) =  dz2(k+1)**2*bunbo
         d1z(k, 1) = -dz (k)**2*bunbo
!
         d2z(k,-1) = 2*dz (k+1)*bunbo
         d2z(k, 0) =-2*dz2(k+1)*bunbo
         d2z(k, 1) = 2*dz (k  )*bunbo
!
      k = nz
!
         dz2(k) = dz(k-1)+dz(k)
        ddz2(k) = 2.0d0/dz2(k)

         bunbo = 1./(dz(k-1)*dz(k-2)*dz2(k-1))
         d1z(k,-1) =  dz (k-1)**2*bunbo
         d1z(k, 0) = -dz2(k-1)**2*bunbo
         d1z(k, 1) =  dz (k-2)*(dz(k-2)+2*dz(k-1))*bunbo
!
         d2z(k,-1) = 2*dz (k-1)*bunbo
         d2z(k, 0) =-2*dz2(k-1)*bunbo
         d2z(k, 1) = 2*dz (k-2)*bunbo         

! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_z_grid : ',dzmin, dzmax, dzchg
! -----------------------------------------------------------

! ------------- write coordinate file Z ---------------------
      write(FILE_COORDINATE_Z,*) index_x, index_y, nx, ny, nz
      write(FILE_COORDINATE_Z,'(e25.16)') zc
! -----------------------------------------------------------
end subroutine pset__z_grid

! -----------------------------------------------------------
subroutine pset__resistivity

implicit none
real(DP) :: cc
integer :: i,j,k

      do i =-1, NX
      do j =-1, NY
      do k = 0, NZ
         cc = sqrt(cx(k,j,i)**2 &
                  +cy(k,j,i)**2 &
                  +cz(k,j,i)**2)    ! current density
         if(cc.gt.cc0) then
            res(k,j,i) = eta + eta1*((cc-cc0)/cc0)**2 ! anomalous resistivity
            write(6,'(A,2x,3i,2x,3i,2x,3i,2x,12f)') "over cc0", k,j,i,cc
         else
            res(k,j,i) = eta ! uniform resistivity
         end if
      end do
      end do
      end do

! ----------------------------top & bottom
       res(0,:,:) = 0.0d0
!       res(0,:,:) = eta
       res(nz,:,:)= 0.0d0
! ----------------------------

! ----------------------------lateral
       if(index_y.eq.0)        res(:,0,:) = 0.0
       if(index_y.eq.nproc_y-1)  res(:,ny,:)= 0.0
! ----------------------------
! ----------------------------lateral also for x
       if(index_x.eq.0)        res(:,:,0) = 0.0
       if(index_x.eq.nproc_x-1)  res(:,:,nx)= 0.0
! ----------------------------
end subroutine pset__resistivity
! -----------------------------------------------------------
subroutine pset__dt
!  set the time step and apply the limiter to density and 
!  and pressure.
! 2005/05/30 acoustic wave is taken into account.
!
implicit none
real(DP) :: fast_wave, diffusion, ttx, tty, ttz 
real(DP) :: dt, dt_diff
integer :: i,j,k
real(DP) :: ttnux, ttnuy, ttnuz
!
  dt = 1.0d0/NY

      do i = 0, NX
      do j = 0, NY
      do k = 0, NZ-1
!
! === density & pressure limiter ============================
         if(ro(k,j,i).le.ro_min) ro(k,j,i) = ro_min
         if(pr(k,j,i).le.pr_min) pr(k,j,i) = pr_min
! ===========================================================
!
         fast_wave = sqrt((gamma*pr(k,j,i)+ &
                           bx(k,j,i)**2+ &
                           by(k,j,i)**2+ &
                           bz(k,j,i)**2) &
                          /ro(k,j,i))
!         fast_wave = 1.0d0
!
!         diffusion = max(visc3d(k,j,i), res(k,j,i), diff_ro, gamma1*akappa/ro(k,j,i))
         diffusion = max(visc, res(k,j,i), diff_ro)
         dt_diff   = min(dx(i),dy(j),dz(k))**2/diffusion/2
!         dt_diff = 1.0d0
!
           ttx = dx(i) / ( abs(vx(k,j,i)) + fast_wave )
           tty = dy(j) / ( abs(vy(k,j,i)) + fast_wave ) 
           ttz = dz(k) / ( abs(vz(k,j,i)) + fast_wave ) 
!
           ttnux = dx(i)*dx(i) / (2.0*visc3d(k,j,i)) 
           ttnuy = dy(j)*dy(j) / (2.0*visc3d(k,j,i)) 
           ttnuz = dz(k)*dz(k) / (2.0*visc3d(k,j,i)) 
!
!           dt = min(dt,ttx,tty,ttz,dt_diff)
           dt = min(dt,ttx,tty,ttz,dt_diff,ttnux,ttnuy,ttnuz)
!---
!if(abs(vx(k,j,i)).gt.1.0d0) write(*,*) k, j, i, "vx", vx(k,j,i), index_x, index_y
!if(abs(vy(k,j,i)).gt.1.0d0) write(*,*) k, j, i, "vy", vy(k,j,i), index_x, index_y
!if(abs(vz(k,j,i)).gt.1.0d0) write(*,*) k, j, i, "vz", vz(k,j,i), index_x, index_y
!---
!if(hanteiNaN.eq.0) then
!	if(vx(k,j,i).ne.vx(k,j,i)) then
!		write(*,*) k, j, i, "vx", vx(k,j,i), index_x, index_y
!		hanteiNaN=1
!	endif
!if(vy(k,j,i).ne.vy(k,j,i)) write(*,*) k, j, i, "vy", vy(k,j,i), index_x, index_y
!	if(vx(k,j,i).ne.vx(k,j,i)) then
!		write(*,*) k, j, i, "vx", vx(k,j,i), index_x, index_y
!		hanteiNaN=1
!	endif
!if(vz(k,j,i).ne.vz(k,j,i)) write(*,*) k, j, i, "vz", vz(k,j,i), index_x, index_y
!	if(vx(k,j,i).ne.vx(k,j,i)) then
!		write(*,*) k, j, i, "vx", vx(k,j,i), index_x, index_y
!		hanteiNaN=1
!	endif
!end if
!---
      end do
      end do
      end do

      dtstep = dt*cfr

! ##### debug #####
!      write(*,*) '## pset__dt: dtstep, myrank ', dtstep, myrank,index_x,index_y
! #################

      call mpiut__min(dtstep)

      return

end subroutine pset__dt
! -----------------------------------------------------------
subroutine cp_double_to_d3d(d,d3d)
real(DP),intent(in) :: d
real(DP),intent(out) :: d3d(0:NZ,-1:NY,-1:NX)

d3d(:,:,:) = d

return
end subroutine cp_double_to_d3d
! -----------------------------------------------------------
end module pset






