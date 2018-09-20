subroutine out_tw
implicit none
	include "common.f90"
!	integer :: itw, jtw!, ktw
	print *, "output file ROT_B on 'z=0 plane'"
	open(30,file="ROT_B",form="unformatted",status="UNKNOWN")
	write(30) twz(:,:,0)
	close(30)
!	do itw = 0, nxl
!	  do jtw = 0, nyl
!            write(30,*) twx(itw,jtw,1), twy(itw,jtw,1), twz(itw,jtw,1)
!            write(30,*) twx(itw,jtw,nzl/2), twy(itw,jtw,nzl/2), twz(itw,jtw,nzl/2)
!	    write(30,*) bx(itw,jtw,1), by(itw,jtw,1), bz(itw,jtw,1)
!	    do ktw = 0, nzl
!		write(30,*) twx(itw,jtw,ktw), twy(itw,jtw,ktw), twz(itw,jtw,ktw)
!	    enddo
!	  enddo
!	enddo
!	close(30)
end subroutine
