module vals_out
use common_val
implicit none

contains
!===========================================
subroutine out_param
!---- check vals
print *, "nx,ny,nz:", nx,ny,nz
print *, "at:", at
print *, "a_surface:", a_suf
print *, "alp:", alp
print *, "Now date:", date0
print *, "Now time:", time0
print *, "para_x0", para_x0
print *, "para_ix0", para_ix0
print *, "para_y0", para_y0
print *, "para_jy0", para_jy0
print *, "para_z0", para_z0
print *, "para_kz0", para_kz0
!--- output parameters
open(15,file=dirname//"/param.dat",status="replace")
write(15,*) "date0:", date0
write(15,*) "time0:", time0
write(15,*) "NNX, NNY, NNZ:", NNX, NNY, NNZ
write(15,*) "nx,ny,nz:", nx,ny,nz
write(15,*) "nx0,ny0,nz0", nx0, ny0, nz0
write(15,*) "rnx,rny,rnz:", rnx,rny,rnz
write(15,*) "Bs0:", Bs0
write(15,*) "at:", at
write(15,*) "a_surface:", a_suf
write(15,*) "alp:", alp
write(15,*) "flatten switch:", flatten_switch
write(15,*) "ydecay:", ydecay
write(15,*) "para_x0", para_x0
write(15,*) "para_ix0", para_ix0
write(15,*) "para_y0", para_y0
write(15,*) "para_jy0", para_jy0
write(15,*) "para_z0", para_z0
write(15,*) "para_kz0", para_kz0
write(15,*) "depth of source", dep
write(15,*) "distance from PIL", lpil
write(15,*) "strength of point source", qps
write(15,*) "widthx sech", decx
write(15,*) "widthy sech", decy
write(15,*) "l from PIL for sech field", lpilsec
write(15,*) "width tanh-ydir around center", dtany
write(15,*) "mask begin position", pmask
write(15,*) "mask width of tanh", dmask
close(15)
end subroutine out_param
!===========================================
subroutine out_coord
integer :: idb2,jdb2,kdb2

print *, dirname
open(10,file=dirname//"/coord.xgc",status="replace")
do idb2 = nx0, nx
write(10,'(f10.5)') xc(idb2)
enddo
close(10)
!
open(11,file=dirname//"/coord.ygc",status="replace")
do jdb2 = ny0, ny
write(11,'(f10.5)') yc(jdb2)
enddo
close(11)
!
open(12,file=dirname//"/coord.zgc",status="replace")
do kdb2 = nz0, nz
write(12,'(f10.5)') zc(kdb2)
enddo
close(12)
!
end subroutine out_coord
!=======================================
subroutine out_values
!integer :: i,j,

write(*,*) "tekitou"
!print *, dirname
open(22,file=dirname//"/Bz_init.dat",status="replace")
	do i = nx0, nx
	do j = ny0, ny
	write(22,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Bz_init(i,j)
	enddo
	enddo
close(22)
!-----
write(*,*) "debug2"
open(23,file=dirname//"/Vx_init.dat",status="replace")
open(24,file=dirname//"/Vy_init.dat",status="replace")
	do i = nx0, nx
	do j = ny0, ny
	write(23,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Vx_init(i,j)
	write(24,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Vy_init(i,j)
	enddo
	enddo
close(23)
close(24)
!-----
!Vxbin_init = Vx_init
!Vybin_init = Vy_init
!Vzbin_init = Vz_init
!open(25,file=dirname//"/Vxbin_init.dat",status="replace",form="unformatted")
!open(26,file=dirname//"/Vybin_init.dat",status="replace",form="unformatted")
!open(27,file=dirname//"/Vybin_init.dat",status="replace",form="unformatted")
!	write(25) Vxbin_init
!	write(26) Vybin_init
!	write(27) Vzbin_init
!close(25)
!close(26)
!close(27)
!-----
!Bxbin = Bx
!Bybin = By
!Bzbin = Bz
!Bx4bin = Bx
!By4bin = By
!Bz4bin = Bz
!open(33,file=dirname//"/Bxbin_3d_init.dat",status="replace",form="unformatted")
!open(34,file=dirname//"/Bybin_3d_init.dat",status="replace",form="unformatted")
!open(35,file=dirname//"/Bzbin_3d_init.dat",status="replace",form="unformatted")
open(36,file=dirname//"/Bx8bin_3d_init.dat",status="replace",form="unformatted")
open(37,file=dirname//"/By8bin_3d_init.dat",status="replace",form="unformatted")
open(38,file=dirname//"/Bz8bin_3d_init.dat",status="replace",form="unformatted")
!	write(33) Bxbin
!	write(34) Bybin
!	write(35) Bzbin
	write(36) Bx
	write(37) By
	write(38) Bz
!close(33)
!close(34)
!close(35)
close(36)
close(37)
close(38)
!
!open(40,file=dirname//"/delBzinit.dat",status="replace")
open(41,file=dirname//"/ck1bz.dat",status="replace")
open(42,file=dirname//"/ck2bz.dat",status="replace")
	do i = nx0, nx
	do j = ny0, ny
!	write(40,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), delBzinit(i,j)
	write(41,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), ck1Bz(i,j)
	write(42,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), ck2Bz(i,j)
	enddo
	enddo
!close(40)
close(41)
close(42)
!
! open(51,file=dirname//"/Psi0_z=0.dat",status="replace")
! !open(52,file=dirname//"/ckbx_z=0.dat",status="replace")
! !open(53,file=dirname//"/ckby_z=0.dat",status="replace")
! !open(54,file=dirname//"/ckbz_z=0.dat",status="replace")
!         do i = nx0, nx
!         do j = ny0, ny
! 	write(51,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Psi0(i,j,0)
! !        write(52,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Bx(i,j,0)
! !        write(53,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), By(i,j,0)
! !        write(54,'(f10.5,3x,f10.5,7x,f12.8)') xc(i), yc(j), Bz(i,j,0)
!         enddo
!         enddo
! close(51)
! !close(52)
! !close(53)
! !close(54)
!=============
open(60,file=dirname//"/B3D_init",status="replace",form="unformatted")
write(60) Bx, By, Bz, Vx_init, Vy_init, Vz_init
close(60)
!=============
!-----
!Vx4bin_init = Vx_init
!Vy4bin_init = Vy_init
!Vz4bin_init = Vz_init
open(70,file=dirname//"/Vx4bin_init.dat",status="replace",form="unformatted")
open(71,file=dirname//"/Vy4bin_init.dat",status="replace",form="unformatted")
open(72,file=dirname//"/Vz4bin_init.dat",status="replace",form="unformatted")
	write(70) real(Vx_init)
	write(71) real(Vy_init)
	write(72) real(Vz_init)
close(70)
close(71)
close(72)
!
end subroutine out_values
!================================
end module vals_out
