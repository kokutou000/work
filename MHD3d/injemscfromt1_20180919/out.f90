! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
module out
! -----------------------------------------------------------
!                    MODULE out
! -----------------------------------------------------------
use common
use mpiut
use mpi
implicit none
real(DP), dimension(2) :: eng, sendbuf, recvbuf

integer,private :: j

contains

! -----------------------------------------------------------
subroutine out__disk
integer :: ierr
!
real(DP), dimension(-1:NX,0:NZ) :: vx_xz
real(DP), dimension(-1:NX,0:NZ) :: vy_xz
real(DP), dimension(-1:NX,0:NZ) :: vz_xz
real(DP), dimension(-1:NX,0:NZ) :: bx_xz
real(DP), dimension(-1:NX,0:NZ) :: by_xz
real(DP), dimension(-1:NX,0:NZ) :: bz_xz
real(DP), dimension(-1:NX,0:NZ) :: cx_xz
real(DP), dimension(-1:NX,0:NZ) :: cy_xz
real(DP), dimension(-1:NX,0:NZ) :: cz_xz
real(DP), dimension(-1:NX,0:NZ) :: ro_xz
real(DP), dimension(-1:NX,0:NZ) :: pr_xz

! #### DEBUG ####
logical :: file_3d_created = .true.
! ###############

! ----------------------
      iwrite = iwrite + 1   ! sequential data number
      cwrite = '+'//chari3(iwrite)
! ----------------------
! write 3d data 
    if(file_3d_created) then
      open(FILE_3D_FIELD,file=trim(cfile_3d_field)//cwrite//cmyrank, &
       form='unformatted')

!#### 20110822 #####
!      write(FILE_3D_FIELD) iwrite,nloop,atime,dtstep, &
!                           bx,by,bz,vx,vy,vz,ro,pr
!#### 20110822 #####
      write(FILE_3D_FIELD) iwrite,nloop,real(atime),real(dtstep), &
                           real(bx),real(by),real(bz), &
                           real(vx),real(vy),real(vz), &
                           real(ro)
!###################
      close(FILE_3D_FIELD)
    end if
! ----------------------
! total energy
      call out__energy(eng(1),eng(2))  ! calculate total energy
!
      sendbuf(:) = eng(:)
      call MPI_REDUCE(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_CART,ierr)

    if(myrank == root) then

      write(FILE_TIME_LIST,'(i3,i8,5e13.6)') iwrite, & 
        nloop, atime, dtstep, recvbuf(1), recvbuf(2), &
        recvbuf(1)+recvbuf(2)

    end if
! ----------------------
! slice data
!  j=0
!  if(mod(nproc_y,2)==0 .and. index_y+1.eq.nproc_y/2) j=ny
!  if(mod(nproc_y,2)==1 .and. index_y  .eq.nproc_y/2) j=ny/2
!  if(j.ne.0) then
!      call out__slice_j(bx,bx_xz,j)
!      call out__slice_j(by,by_xz,j)
!      call out__slice_j(bz,bz_xz,j)
!      call out__slice_j(vx,vx_xz,j)
!      call out__slice_j(vy,vy_xz,j)
!      call out__slice_j(vz,vz_xz,j)
!      call out__slice_j(ro,ro_xz,j)
!      call out__slice_j(pr,pr_xz,j)
!      call out__slice_j(cx,cx_xz,j)
!      call out__slice_j(cy,cy_xz,j)
!      call out__slice_j(cz,cz_xz,j)
!!
!      open(FILE_SLICE_XZ,file=trim(cfile_slice_xz)//cwrite//cmyrank, &
!      form='unformatted')
!
!      write(FILE_SLICE_XZ) iwrite,bx_xz,by_xz,bz_xz, &
!                                  vx_xz,vy_xz,vz_xz, &
!                                  ro_xz,pr_xz, &
!                                  cx_xz,cy_xz,cz_xz
!      close(FILE_SLICE_XZ)
!   end if
!
 return

end subroutine out__disk
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine out__energy(engb,engv)
real(DP), INTENT(OUT) :: engb, engv
integer :: i,j,k
!
      engb = 0.0d0
      engv = 0.0d0
! ------------------ boundaries (top and bottom)
      do i = 0, NX-1
      do j = 0, NY-1
      do k = 0, NZ, NZ
         engb = engb + (bx(k,j,i)**2 &
                       +by(k,j,i)**2 &
                       +bz(k,j,i)**2)*dx(i)*dy(j)*dz(k)
         engv = engv + (vx(k,j,i)**2 &
                       +vy(k,j,i)**2 &
                       +vz(k,j,i)**2)*dx(i)*dy(j)*dz(k)
      end do
      end do
      end do
!
         engb = 0.5d0 * engb
         engv = 0.5d0 * engv
!
! ----------------- integration in the  domain
      do i = 0, NX-1
      do j = 0, NY-1
      do k = 1, NZ-1
         engb = engb + (bx(k,j,i)**2 &
                       +by(k,j,i)**2 &
                       +bz(k,j,i)**2)*dx(i)*dy(j)*dz(k)
         engv = engv + (vx(k,j,i)**2 &
                       +vy(k,j,i)**2 &
                       +vz(k,j,i)**2)*dx(i)*dy(j)*dz(k)
      end do
      end do
      end do
!
! ------------------ Eb = (1/2) int B^2 dv
! ------------------ Ev = (1/2) int V^2 dv
      engb = 0.5d0 * engb 
      engv = 0.5d0 * engv 
!

      return
end subroutine out__energy
! -------------------------------------------------------
subroutine out__slice_j(d3d,d2d,j_slice)
real(DP), INTENT(IN), dimension (0:NZ,-1:NY,-1:NX) :: d3d
real(DP), INTENT(OUT), dimension(-1:NX,0:NZ) :: d2d
integer,  INTENT(IN) :: j_slice
integer :: i, k
!
      do i = -1, nx
      do k =  0, nz
        d2d(i,k) = d3d(k,j_slice,i)
      end do
      end do

      return
end subroutine out__slice_j
! -------------------------------------------------------
subroutine out__restart

 open(FILE_RESTART,file=trim(cfile_restart)//cmyrank, &
!!      convert='big_endian', & ! for Absoft_f90
      form='unformatted')

  write(FILE_RESTART) iwrite, nloop, atime, dtstep, &
                      bx, by, bz, vx, vy, vz, ro, pr
 close(FILE_RESTART)

end subroutine out__restart

end module out
