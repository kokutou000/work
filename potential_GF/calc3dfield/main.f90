program main_prog
use common_val
use vals_out
implicit none
!integer :: i,j,k
character :: date0time0*14
!
allocate(Bx(nx0:nx,ny0:ny,nz0:nz))
allocate(By(nx0:nx,ny0:ny,nz0:nz))
allocate(Bz(nx0:nx,ny0:ny,nz0:nz))
!

call date_and_time(DATE=date0,TIME=time0)
date0time0 = date0//time0
print *, "date0time0: ", date0time0
dirname="output_dat/"//date0time0
print *, "dirname =", dirname
call system('mkdir '//dirname)
!
call out_param
!-----
call coord_init
!
call out_coord
!-----

call boundary_condition
call calc_velocity
call green_function_method
call calc_current_divb

write(*,*) "debug1"
!-----
call out_values

end
