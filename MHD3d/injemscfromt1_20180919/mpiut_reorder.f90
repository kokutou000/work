! ===========================================================
!      3D zero-beta nonlinear MHD Model 
! ===========================================================
module mpiut
! -----------------------------------------------------------
!                    MODULE mpiut
! -----------------------------------------------------------
use common
use mpi

implicit none

integer :: ierr
integer :: dims(2), coords(2)
logical :: periods(2), reorder

contains

! -----------------------------------------------------------
subroutine mpiut__init

  call mpi_init(ierr)

!-------------- Original communicator --------
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

! ###### debug #######
  write(*,*) ' #### communicator(WORLD), nproc, myrank ',MPI_COMM_WORLD,nproc,myrank
! ####################
!-------------- Cartesion communicator -------

  dims(1) = nproc_x
  dims(2) = nproc_y

  periods(1) = .true.
  periods(2) = .false.

  reorder = .true.

  call mpi_cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,MPI_COMM_CART,ierr)

  call mpi_comm_size(MPI_COMM_CART,nproc,ierr)
  call mpi_comm_rank(MPI_COMM_CART,myrank,ierr)
!
  if(nproc .ne. nproc_x*nproc_y) then
     write(*,*) '::ERR:mpiut__init, nproc not equal to nproc_x*nproc_y', &
                nproc, nproc_x, nproc_y
     call mpiut__finalize
     stop
  end if

! ++++++++++++++ My Coordinate ++++++++++++++
  call MPI_CART_COORDS(MPI_COMM_CART, myrank, 2, coords, ierr)
  index_x = coords(1)
  index_y = coords(2)

! ++++++++++++++ Coordinate of neighbors ++++++++++++++
  call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, rank_left, rank_right, ierr)
  call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, rank_down, rank_up, ierr)

! ###### debug #######
  write(*,*) ' #### communicator(CART), nproc, myrank ',MPI_COMM_CART,nproc,myrank
  write(*,*) ' #### index(x,y), neighbors ',index_x, index_y, rank_left, rank_right, rank_down, rank_up
  write(*,*)
! ####################

end subroutine mpiut__init
! -----------------------------------------------------------
subroutine mpiut__finalize

  call mpi_finalize(ierr)

end subroutine mpiut__finalize
! -----------------------------------------------------------
subroutine mpiut__barrier

  call MPI_BARRIER(MPI_COMM_CART,ierr)

end subroutine mpiut__barrier
! -----------------------------------------------------------
subroutine mpiut__min(d)

  real(DP) :: sendbuf, recvbuf, d

  sendbuf = d

  call MPI_REDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                  root, MPI_COMM_CART, ierr)

  call MPI_BCAST(recvbuf, 1, MPI_DOUBLE_PRECISION, &
                 root, MPI_COMM_CART, ierr)

  d = recvbuf

end subroutine mpiut__min
! -----------------------------------------------------------
subroutine mpiut__exchange_sub_x

  real(DP), dimension(0:NZ,-1:NY,6) :: sendbuf, recvbuf
  integer :: sendcount = (NY+2)*(NZ+1)*6
  integer :: recvcount = (NY+2)*(NZ+1)*6
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------ from left to right

  sendbuf(:,:,1) = cx(:,:,NX-1)
  sendbuf(:,:,2) = cy(:,:,NX-1)
  sendbuf(:,:,3) = cz(:,:,NX-1)
  sendbuf(:,:,4) = qx(:,:,NX-1)
  sendbuf(:,:,5) = qy(:,:,NX-1)
  sendbuf(:,:,6) = qz(:,:,NX-1)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_right, 0, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_left, 0, &
                    MPI_COMM_CART, status, ierr)

  cx(:,:,-1) = recvbuf(:,:,1)
  cy(:,:,-1) = recvbuf(:,:,2)
  cz(:,:,-1) = recvbuf(:,:,3)
  qx(:,:,-1) = recvbuf(:,:,4)
  qy(:,:,-1) = recvbuf(:,:,5)
  qz(:,:,-1) = recvbuf(:,:,6)

! ------------ from right to left

  sendbuf(:,:,1) = cx(:,:,0)
  sendbuf(:,:,2) = cy(:,:,0)
  sendbuf(:,:,3) = cz(:,:,0)
  sendbuf(:,:,4) = qx(:,:,0)
  sendbuf(:,:,5) = qy(:,:,0)
  sendbuf(:,:,6) = qz(:,:,0)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_left, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_right, 1, &
                    MPI_COMM_CART, status, ierr)

  cx(:,:,NX) = recvbuf(:,:,1)
  cy(:,:,NX) = recvbuf(:,:,2)
  cz(:,:,NX) = recvbuf(:,:,3)
  qx(:,:,NX) = recvbuf(:,:,4)
  qy(:,:,NX) = recvbuf(:,:,5)
  qz(:,:,NX) = recvbuf(:,:,6)

return
end subroutine mpiut__exchange_sub_x
! -----------------------------------------------------------
subroutine mpiut__exchange_sub_y

  real(DP), dimension(0:NZ,-1:NX,6) :: sendbuf, recvbuf
  integer :: sendcount = (NX+2)*(NZ+1)*6
  integer :: recvcount = (NX+2)*(NZ+1)*6
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------ from down to up

    sendbuf(:,:,1) = cx(:,NY-1,:)
    sendbuf(:,:,2) = cy(:,NY-1,:)
    sendbuf(:,:,3) = cz(:,NY-1,:)
    sendbuf(:,:,4) = qx(:,NY-1,:)
    sendbuf(:,:,5) = qy(:,NY-1,:)
    sendbuf(:,:,6) = qz(:,NY-1,:)

    call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 0, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 0, &
                    MPI_COMM_CART, status, ierr)

! #### not periodic boundary condition ####
  if(index_y.ne.0) then
    cx(:,-1,:) = recvbuf(:,:,1)
    cy(:,-1,:) = recvbuf(:,:,2)
    cz(:,-1,:) = recvbuf(:,:,3)
    qx(:,-1,:) = recvbuf(:,:,4)
    qy(:,-1,:) = recvbuf(:,:,5)
    qz(:,-1,:) = recvbuf(:,:,6)
  end if

! ------------ from up to down

  sendbuf(:,:,1) = cx(:,0,:)
  sendbuf(:,:,2) = cy(:,0,:)
  sendbuf(:,:,3) = cz(:,0,:)
  sendbuf(:,:,4) = qx(:,0,:)
  sendbuf(:,:,5) = qy(:,0,:)
  sendbuf(:,:,6) = qz(:,0,:)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 1, &
                    MPI_COMM_CART, status, ierr)

! #### not periodic boundary condition ####
  if(index_y.ne.nproc_y-1) then
   cx(:,NY,:) = recvbuf(:,:,1)
   cy(:,NY,:) = recvbuf(:,:,2)
   cz(:,NY,:) = recvbuf(:,:,3)
   qx(:,NY,:) = recvbuf(:,:,4)
   qy(:,NY,:) = recvbuf(:,:,5)
   qz(:,NY,:) = recvbuf(:,:,6)
  end if

return
end subroutine mpiut__exchange_sub_y
! -----------------------------------------------------------
subroutine mpiut__exchange_d_x

  real(DP), dimension(0:NZ,-1:NY,8) :: sendbuf, recvbuf
  integer :: sendcount = (NY+2)*(NZ+1)*8
  integer :: recvcount = (NY+2)*(NZ+1)*8
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------ from left to right

  sendbuf(:,:,1) = dvx(:,:,NX-1)
  sendbuf(:,:,2) = dvy(:,:,NX-1)
  sendbuf(:,:,3) = dvz(:,:,NX-1)
  sendbuf(:,:,4) = dbx(:,:,NX-1)
  sendbuf(:,:,5) = dby(:,:,NX-1)
  sendbuf(:,:,6) = dbz(:,:,NX-1)
  sendbuf(:,:,7) = dro(:,:,NX-1)
  sendbuf(:,:,8) = dpr(:,:,NX-1)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_right, 2, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_left, 2, &
                    MPI_COMM_CART, status, ierr)

  dvx(:,:,-1) = recvbuf(:,:,1)
  dvy(:,:,-1) = recvbuf(:,:,2)
  dvz(:,:,-1) = recvbuf(:,:,3)
  dbx(:,:,-1) = recvbuf(:,:,4)
  dby(:,:,-1) = recvbuf(:,:,5)
  dbz(:,:,-1) = recvbuf(:,:,6)
  dro(:,:,-1) = recvbuf(:,:,7)
  dpr(:,:,-1) = recvbuf(:,:,8)

! ------------ from right to left

  sendbuf(:,:,1) = dvx(:,:,0)
  sendbuf(:,:,2) = dvy(:,:,0)
  sendbuf(:,:,3) = dvz(:,:,0)
  sendbuf(:,:,4) = dbx(:,:,0)
  sendbuf(:,:,5) = dby(:,:,0)
  sendbuf(:,:,6) = dbz(:,:,0)
  sendbuf(:,:,7) = dro(:,:,0)
  sendbuf(:,:,8) = dpr(:,:,0)


  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_left, 3, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_right, 3, &
                    MPI_COMM_CART, status, ierr)

  dvx(:,:,NX) = recvbuf(:,:,1)
  dvy(:,:,NX) = recvbuf(:,:,2)
  dvz(:,:,NX) = recvbuf(:,:,3)
  dbx(:,:,NX) = recvbuf(:,:,4)
  dby(:,:,NX) = recvbuf(:,:,5)
  dbz(:,:,NX) = recvbuf(:,:,6)
  dro(:,:,NX) = recvbuf(:,:,7)
  dpr(:,:,NX) = recvbuf(:,:,8)

return
end subroutine mpiut__exchange_d_x
subroutine mpiut__exchange_d_y

  real(DP), dimension(0:NZ,-1:NX,8) :: sendbuf, recvbuf
  integer :: sendcount = (NX+2)*(NZ+1)*8
  integer :: recvcount = (NX+2)*(NZ+1)*8
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------ from down to up

  sendbuf(:,:,1) = dvx(:,NY-1,:)
  sendbuf(:,:,2) = dvy(:,NY-1,:)
  sendbuf(:,:,3) = dvz(:,NY-1,:)
  sendbuf(:,:,4) = dbx(:,NY-1,:)
  sendbuf(:,:,5) = dby(:,NY-1,:)
  sendbuf(:,:,6) = dbz(:,NY-1,:)
  sendbuf(:,:,7) = dro(:,NY-1,:)
  sendbuf(:,:,8) = dpr(:,NY-1,:)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 2, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 2, &
                    MPI_COMM_CART, status, ierr)

! ######################################################
! #### fixed boundary condition on the lateral wall ####
! ######################################################

  if(index_y.ne.0) then
    dvx(:,-1,:) = recvbuf(:,:,1)
    dvy(:,-1,:) = recvbuf(:,:,2)
    dvz(:,-1,:) = recvbuf(:,:,3)
    dbx(:,-1,:) = recvbuf(:,:,4)
    dby(:,-1,:) = recvbuf(:,:,5)
    dbz(:,-1,:) = recvbuf(:,:,6)
    dro(:,-1,:) = recvbuf(:,:,7)
    dpr(:,-1,:) = recvbuf(:,:,8)
  end if

! ------------ from right to left

  sendbuf(:,:,1) = dvx(:,0,:)
  sendbuf(:,:,2) = dvy(:,0,:)
  sendbuf(:,:,3) = dvz(:,0,:)
  sendbuf(:,:,4) = dbx(:,0,:)
  sendbuf(:,:,5) = dby(:,0,:)
  sendbuf(:,:,6) = dbz(:,0,:)
  sendbuf(:,:,7) = dro(:,0,:)
  sendbuf(:,:,8) = dpr(:,0,:)

  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 3, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 3, &
                    MPI_COMM_CART, status, ierr)

! ######################################################
! #### fixed boundary condition on the lateral wall ####
! ######################################################

  if(index_y.ne.nproc_y-1) then
    dvx(:,NY,:) = recvbuf(:,:,1)
    dvy(:,NY,:) = recvbuf(:,:,2)
    dvz(:,NY,:) = recvbuf(:,:,3)
    dbx(:,NY,:) = recvbuf(:,:,4)
    dby(:,NY,:) = recvbuf(:,:,5)
    dbz(:,NY,:) = recvbuf(:,:,6)
    dro(:,NY,:) = recvbuf(:,:,7)
    dpr(:,NY,:) = recvbuf(:,:,8)
  end if

return
end subroutine mpiut__exchange_d_y
end module mpiut




