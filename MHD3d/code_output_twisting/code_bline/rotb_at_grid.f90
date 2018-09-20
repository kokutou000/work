subroutine calc_rotb_grid
implicit none
	include "common.f90"
	integer :: itw, jtw, ktw
	real(8) :: dzp1, dzp2, dzm1, dzm2
	real(8) :: bb2
!----- bb2 -> sqrt(B*B) = B
	real(8) :: dxp, dxm, dyp, dym, dzp, dzm
!----- dxm -> delta x minus
!----- dxp -> delta x plus

!	do itw = 0, nxl
!	  do jtw = 0, nyl
!	    do ktw = 0, nzl
!----- rotbx, rotby, rotbz -> rot B / B
!----- so, = alpha_x, alpha_y, alpha_z
	rotbx=0d0; rotby=0d0; rotbz=0d0
	print *, "nxl, nyl, nzl =", nxl, nyl, nzl
	do itw = 1, nxl-1
	  do jtw = 1, nyl-1
	    do ktw = 1, nzl-1
		dxp = xc(itw+1)-xc(itw)
		dxm = xc(itw)-xc(itw-1)
		dyp = yc(jtw+1)-yc(jtw)
		dym = yc(jtw)-yc(jtw-1)
		dzp = zc(ktw+1)-zc(ktw)
		dzm = zc(ktw)-zc(ktw-1)
!		if(itw==nxl/2.and.ktw==nzl/2.and.jtw==nyl/2) print *, dxp, dxm, dyp, dym, dzp, dzm
		rotbx(itw,jtw,ktw) = &
		& (bz(itw,jtw+1,ktw)*dym*dym&
		&+ bz(itw,jtw  ,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&-(bz(itw,jtw  ,ktw)*dym*dym&
		&+ bz(itw,jtw-1,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&-(by(itw,jtw,ktw+1)*dzm*dzm&
		&+ by(itw,jtw,ktw  )*dzp*dzp&
		&)/(dzp*dzm*(dzp+dzm))&
		&+(by(itw,jtw,ktw  )*dzm*dzm&
		&+ by(itw,jtw,ktw-1)*dzp*dzp&
		&)/(dzp*dzm*(dzp+dzm))
!		
!		& (dym*(bz(itw,jtw+1,ktw)-bz(itw,jtw  ,ktw))/(dyp*(dyp+dym))&
!		&+ dyp*(bz(itw,jtw  ,ktw)-bz(itw,jtw-1,ktw))/(dym*(dyp+dym)))&
!		&-(dzm*(by(itw,jtw,ktw+1)-by(itw,jtw,ktw  ))/(dzp*(dzp+dzm))&
!		&+ dzp*(by(itw,jtw,ktw  )-by(itw,jtw,ktw-1))/(dzm*(dzp+dzm)))
!		
!		& (dym*dym*(bz(itw,jtw+1,ktw)-bz(itw,jtw  ,ktw))&
!		&+ dyp*dyp*(bz(itw,jtw  ,ktw)-bz(itw,jtw-1,ktw))&
!		&)/(dyp*dym*(dyp+dym))&
!		&-(dzm*dzm*(by(itw,jtw,ktw+1)-by(itw,jtw,ktw  ))&
!		&+ dzp*dzp*(by(itw,jtw,ktw  )-by(itw,jtw,ktw-1))&
!		&)/(dzp*dzm*(dzp+dzm))
		!
		rotby(itw,jtw,ktw) = &
		& (bx(itw,jtw,ktw+1)*dzm*dzm&
		&+ bx(itw,jtw,ktw  )*dzp*dzp&
		&)/(dzp*dzm*(dzp+dzm))&
		&-(bx(itw,jtw,ktw  )*dzm*dzm&
		&+ bx(itw,jtw,ktw-1)*dzp*dzp&
		&)/(dzp*dzm*(dzp+dzm))&
		&-(bz(itw+1,jtw,ktw)*dxm*dxm&
		&+ bz(itw  ,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&+(bz(itw  ,jtw,ktw)*dxm*dxm&
		&+ bz(itw-1,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))
!
!		& (dzm*(bx(itw,jtw,ktw+1)-bx(itw,jtw,ktw  ))/(dzp*(dzp+dzm))&
!		&+ dzp*(bx(itw,jtw,ktw  )-bx(itw,jtw,ktw-1))/(dzm*(dzp+dzm)))&
!		&-(dxm*(bz(itw+1,jtw,ktw)-bz(itw  ,jtw,ktw))/(dxp*(dxp+dxm))&
!		&+ dxp*(bz(itw  ,jtw,ktw)-bz(itw-1,jtw,ktw))/(dxm*(dxp+dxm)))
!
!		& (dzm*dzm*(bx(itw,jtw,ktw+1)-bx(itw,jtw,ktw  ))&
!		&+ dzp*dzp*(bx(itw,jtw,ktw  )-bx(itw,jtw,ktw-1))&
!		&)/(dzp*dzm*(dzp+dzm))&
!		&-(dxm*dxm*(bz(itw+1,jtw,ktw)-bz(itw  ,jtw,ktw))&
!		&+ dxp*dxp*(bz(itw  ,jtw,ktw)-bz(itw-1,jtw,ktw))&
!		&)/(dxp*dxm*(dxp+dxm))
		!
		rotbz(itw,jtw,ktw) = &
		& (by(itw+1,jtw,ktw)*dxm*dxm&
		&+ by(itw  ,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&-(by(itw  ,jtw,ktw)*dxm*dxm&
		&+ by(itw-1,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&-(bx(itw,jtw+1,ktw)*dym*dym&
		&+ bx(itw,jtw  ,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&+(bx(itw,jtw  ,ktw)*dym*dym&
		&+ bx(itw,jtw-1,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))
!
!		& (dxm*(by(itw+1,jtw,ktw)-by(itw  ,jtw,ktw))/(dxp*(dxp+dxm))&
!		&+ dxp*(by(itw  ,jtw,ktw)-by(itw-1,jtw,ktw))/(dxm*(dxp+dxm)))&
!		&-(dym*(bx(itw,jtw+1,ktw)-bx(itw,jtw  ,ktw))/(dyp*(dyp+dym))&
!		&+ dyp*(bx(itw,jtw  ,ktw)-bx(itw,jtw-1,ktw))/(dym*(dyp+dym)))
!
!		& (dxm*dxm*(by(itw+1,jtw,ktw)-by(itw  ,jtw,ktw))&
!		&+ dxp*dxp*(by(itw  ,jtw,ktw)-by(itw-1,jtw,ktw))&
!		&)/(dxp*dxm*(dxp+dxm))&
!		&-(dym*dym*(bx(itw,jtw+1,ktw)-bx(itw,jtw  ,ktw))&
!		&+ dyp*dyp*(bx(itw,jtw  ,ktw)-bx(itw,jtw-1,ktw))&
!		&)/(dyp*dym*(dyp+dym))
!		bb2 = sqrt(bx(itw,jtw,ktw)**2+by(itw,jtw,ktw)**2+bz(itw,jtw,ktw)**2)
!		if(itw==nxl/2.and.ktw==nzl/2.and.jtw==nyl/2) print *, rotbx(itw,jtw,ktw), rotby(itw,jtw,ktw), rotbz(itw,jtw,ktw)
!		if(itw==nxl/2.and.ktw==nzl/2.and.jtw==nyl/2) print *, bb2
!		rotbx(itw,jtw,ktw) = rotbx(itw,jtw,ktw)/bb2
!		rotby(itw,jtw,ktw) = rotby(itw,jtw,ktw)/bb2
!		rotbz(itw,jtw,ktw) = rotbz(itw,jtw,ktw)/bb2
!		if(itw==nxl/2.and.ktw==nzl/2.and.jtw==nyl/2) print *, rotbx(itw,jtw,ktw), rotby(itw,jtw,ktw), rotbz(itw,jtw,ktw)
	    end do
	  end do
	end do
!----- calculate rot(B)_z on bottom and ceiling -----
	do itw = 1, nxl-1
	do jtw = 1, nyl-1
	dxp = xc(itw+1)-xc(itw)
	dxm = xc(itw)-xc(itw-1)
	dyp = yc(jtw+1)-yc(jtw)
	dym = yc(jtw)-yc(jtw-1)
	!
	ktw = 0
	dzp1 = zc(ktw+1) - zc(ktw)
	dzp2 = zc(ktw+2) - zc(ktw)
		!
		rotbx(itw,jtw,ktw) = &
		& (bz(itw,jtw+1,ktw)*dym*dym&
		&+ bz(itw,jtw  ,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&-(bz(itw,jtw  ,ktw)*dym*dym&
		&+ bz(itw,jtw-1,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&-(by(itw,jtw,ktw+1)*dzp2*dzp2&
		&+ by(itw,jtw,ktw  )*dzp1*dzp1&
		&)/(dzp1*dzp2*(dzp2-dzp1))&
		&+(by(itw,jtw,ktw  )*dzp2*dzp2&
		&+ by(itw,jtw,ktw+2)*dzp1*dzp1&
		&)/(dzp1*dzp2*(dzp2-dzp1))
!
!		& (dym*(bz(itw,jtw+1,ktw)-bz(itw,jtw  ,ktw))/(dyp*(dyp+dym))&
!		&+ dyp*(bz(itw,jtw  ,ktw)-bz(itw,jtw-1,ktw))/(dym*(dyp+dym)))&
!		&-(dzp2*(by(itw,jtw,ktw+1)-by(itw,jtw,ktw  ))/(dzp1*(dzp2-dzp1))&
!		&+ dzp1*(by(itw,jtw,ktw  )-by(itw,jtw,ktw+2))/(dzp2*(dzp2-dzp1)))
!
!		& (dym*dym*(bz(itw,jtw+1,ktw)-bz(itw,jtw  ,ktw))&
!		&+ dyp*dyp*(bz(itw,jtw  ,ktw)-bz(itw,jtw-1,ktw))&
!		&)/(dyp*dym*(dyp+dym))&
!		&-(dzp2*dzp2*(by(itw,jtw,ktw+1)-by(itw,jtw,ktw  ))&
!		&+ dzp1*dzp1*(by(itw,jtw,ktw  )-by(itw,jtw,ktw+2))&
!		&)/(dzp1*dzp2*(dzp2-dzp1))
		!
		rotby(itw,jtw,ktw) = &
		& (bx(itw,jtw,ktw+1)*dzp2*dzp2&
		&+ bx(itw,jtw,ktw  )*dzp1*dzp1&
		&)/(dzp1*dzp2*(dzp2-dzp1))&
		&-(bx(itw,jtw,ktw  )*dzp2*dzp2&
		&+ bx(itw,jtw,ktw+2)*dzp1*dzp1&
		&)/(dzp1*dzp2*(dzp2-dzp1))&
		&-(bz(itw+1,jtw,ktw)*dxm*dxm&
		&+ bz(itw  ,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&+(bz(itw  ,jtw,ktw)*dxm*dxm&
		&+ bz(itw-1,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))
!
!		& (dzp2*(bx(itw,jtw,ktw+1)-bx(itw,jtw,ktw  ))/(dzp1*(dzp2-dzp1))&
!		&+ dzp1*(bx(itw,jtw,ktw  )-bx(itw,jtw,ktw+2))/(dzp2*(dzp2-dzp1)))&
!		&-(dxm*(bz(itw+1,jtw,ktw)-bz(itw  ,jtw,ktw))/(dxp*(dxp+dxm))&
!		&+ dxp*(bz(itw  ,jtw,ktw)-bz(itw-1,jtw,ktw))/(dxm*(dxp+dxm)))
!
!		& (dzp2*dzp2*(bx(itw,jtw,ktw+1)-bx(itw,jtw,ktw  ))&
!		&+ dzp1*dzp1*(bx(itw,jtw,ktw  )-bx(itw,jtw,ktw+2))&
!		&)/(dzp1*dzp2*(dzp2-dzp1))&
!		&-(dxm*dxm*(bz(itw+1,jtw,ktw)-bz(itw  ,jtw,ktw))&
!		&+ dxp*dxp*(bz(itw  ,jtw,ktw)-bz(itw-1,jtw,ktw))&
!		&)/(dxp*dxm*(dxp+dxm))
		!
		rotbz(itw,jtw,ktw) = &
		& (by(itw+1,jtw,ktw)*dxm*dxm&
		&+ by(itw  ,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&-(by(itw  ,jtw,ktw)*dxm*dxm&
		&+ by(itw-1,jtw,ktw)*dxp*dxp&
		&)/(dxp*dxm*(dxp+dxm))&
		&-(bx(itw,jtw+1,ktw)*dym*dym&
		&+ bx(itw,jtw  ,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))&
		&+(bx(itw,jtw  ,ktw)*dym*dym&
		&+ bx(itw,jtw-1,ktw)*dyp*dyp&
		&)/(dyp*dym*(dyp+dym))
!
!		& (dxm*(by(itw+1,jtw,ktw)-by(itw  ,jtw,ktw))/(dxp*(dxp+dxm))&
!		&+ dxp*(by(itw  ,jtw,ktw)-by(itw-1,jtw,ktw))/(dxm*(dxp+dxm)))&
!		&-(dym*(bx(itw,jtw+1,ktw)-bx(itw,jtw  ,ktw))/(dyp*(dyp+dym))&
!		&+ dyp*(bx(itw,jtw  ,ktw)-bx(itw,jtw-1,ktw))/(dym*(dyp+dym)))
!
!		& (dxm*dxm*(by(itw+1,jtw,ktw)-by(itw  ,jtw,ktw))&
!		&+ dxp*dxp*(by(itw  ,jtw,ktw)-by(itw-1,jtw,ktw))&
!		&)/(dxp*dxm*(dxp+dxm))&
!		&-(dym*dym*(bx(itw,jtw+1,ktw)-bx(itw,jtw  ,ktw))&
!		&+ dyp*dyp*(bx(itw,jtw  ,ktw)-bx(itw,jtw-1,ktw))&
!		&)/(dyp*dym*(dyp+dym))
	!
	ktw = nzl
	dzm1 = zc(ktw) - zc(ktw-1)
	dzm2 = zc(ktw) - zc(ktw-2)
		!
		rotbx(itw,jtw,ktw) = &
		& (dym*dym*(bz(itw,jtw+1,ktw)-bz(itw,jtw  ,ktw))&
		&+ dyp*dyp*(bz(itw,jtw  ,ktw)-bz(itw,jtw-1,ktw))&
		&)/(dyp*dym*(dyp+dym))&
		&-(dzm2*dzm2*(by(itw,jtw,ktw-1)-by(itw,jtw,ktw  ))&
		&+ dzm1*dzm1*(by(itw,jtw,ktw  )-by(itw,jtw,ktw-2))&
		&)/(dzm1*dzm2*(dzm1-dzm2))
		!
		rotby(itw,jtw,ktw) = &
		& (dzm2*dzm2*(bx(itw,jtw,ktw-1)-bx(itw,jtw,ktw  ))&
		&+ dzm1*dzm1*(bx(itw,jtw,ktw  )-bx(itw,jtw,ktw-2))&
		&)/(dzm1*dzm2*(dzm1-dzm2))&
		&-(dxm*dxm*(bz(itw+1,jtw,ktw)-bz(itw  ,jtw,ktw))&
		&+ dxp*dxp*(bz(itw  ,jtw,ktw)-bz(itw-1,jtw,ktw))&
		&)/(dxp*dxm*(dxp+dxm))
		!
		rotbz(itw,jtw,ktw) = &
		& (dxm*dxm*(by(itw+1,jtw,ktw)-by(itw  ,jtw,ktw))&
		&+ dxp*dxp*(by(itw  ,jtw,ktw)-by(itw-1,jtw,ktw))&
		&)/(dxp*dxm*(dxp+dxm))&
		&-(dym*dym*(bx(itw,jtw+1,ktw)-bx(itw,jtw  ,ktw))&
		&+ dyp*dyp*(bx(itw,jtw  ,ktw)-bx(itw,jtw-1,ktw))&
		&)/(dyp*dym*(dyp+dym))
	end do
	end do
!
!----- calculate force free alpha on z = 0
!
	do itw = 0, nxl
	do jtw = 0, nyl
		alpxy(itw,jtw) = &
		&(rotbx(itw,jtw,0)*bx(itw,jtw,0)&
		&+rotby(itw,jtw,0)*by(itw,jtw,0)&
		&+rotbz(itw,jtw,0)*bz(itw,jtw,0))&
		&/(bx(itw,jtw,0)**2+by(itw,jtw,0)**2+bz(itw,jtw,0)**2)
		!
		alpxy1(itw,jtw) = &
		&(rotbx(itw,jtw,1)*bx(itw,jtw,1)&
		&+rotby(itw,jtw,1)*by(itw,jtw,1)&
		&+rotbz(itw,jtw,1)*bz(itw,jtw,1))&
		&/(bx(itw,jtw,1)**2+by(itw,jtw,1)**2+bz(itw,jtw,1)**2)
		!
		alpxy2(itw,jtw) = &
		&(rotbx(itw,jtw,2)*bx(itw,jtw,2)&
		&+rotby(itw,jtw,2)*by(itw,jtw,2)&
		&+rotbz(itw,jtw,2)*bz(itw,jtw,2))&
		&/(bx(itw,jtw,2)**2+by(itw,jtw,2)**2+bz(itw,jtw,2)**2)
	end do
	end do
	
return

end subroutine
