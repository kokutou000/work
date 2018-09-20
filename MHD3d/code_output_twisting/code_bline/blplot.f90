      subroutine blplot(np,x,y,z,idir,tw,xfin,yfin,zfin)
      include "common.f90"
      real*8 x(np),y(np),z(np)
      real(8) :: xfin,yfin,zfin 
      real(8) :: tw
	integer :: tmp
!
!      print *, "np=",np,"rbl=",rbl,"gbl",gbl
!      print *, "bbl",bbl,"wbl",wbl,"idir",idir
!       write(20,*) np,rbl,gbl,bbl,wbl,idir
!       do i = 1, np
!        write(20,*) x(i),y(i),z(i)
!       end do
!
!       write(40,*) x(1), y(1), tw, idir
	if(rl==1) then
		write(50,*) x(1), y(1), z(1)
		write(50,*) xfin, yfin, zfin
	else if(rl==-1) then
		write(55,*) x(1), y(1), z(1)
		write(55,*) xfin, yfin, zfin
	else
!		print *, "right left become -->>", rl
		tmp = 0
	endif
       write(40,*) y(1), z(1), tw, idir, np

      return
      end subroutine blplot
