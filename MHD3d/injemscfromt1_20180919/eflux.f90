! ===========================================================
!     Subroutine for the setting the boundary condition in flux emerging 
!     process NL3DpwE_mpi_02h_nlff by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
!                    SUBROUTINE eflux__mag
!                   elliptical cross-section
! -----------------------------------------------------------
! &NLIST_EF
! eflux_x0 = 0.25, ! X coordinate of emerging flux center
! eflux_y0 = 0.25, ! Y coordinate of emerging flux center
! eflux_z0 =-0.1,  ! Z coordinate of emerging flux center on the initial 
! eflux_ph =45.0,  ! The azimuthal angle of horizontal field of emerging flux
!                  ! f.g. eflux_ph=0 means B_t directed to the X-axis.
! eflux_vz = 0.1,  ! The emerging speed.
! eflux_t0 = 0.0,  ! The time the emerging starts.
! eflux_t1 = 0.1,  ! The time the emerging stops.
! eflux_b0 = 1.0,  ! The intensity of the emerging magnetic field.
! eflux_r0 = 0.1,  ! The radius of emerging flux.
! eflux_a0 = 0.1,  ! The width of emerging flux.
! eflux_rm = 0.01, ! The innermost radius of emerging flux
!         so that 
!         if c_r < eflux_r0 and |c_z| < eflux_a0, then
!                    |B| = eflux_b0*(1.0-(c_z/eflux_a0)**2)* &
!                          eflux_r0/max(c_r,eflux_rm)
!         else       |B| = 0.0,
!         where (c_r,c_z) is the cylilndrical coordinate, in which
!         the origin is located at the emerging flux center.
! &END
! -----------------------------------------------------------
!
subroutine eflux__mag(bx2d,by2d,bz2d,vz2d,time_eflux)
use common
!use test_eflux_subroutines

implicit none

real(DP), dimension(-1:NY,-1:NX), intent(OUT) :: bx2d, by2d, bz2d, vz2d
real(DP) :: time_eflux, eflux_z
real(DP) :: xx, yy, zz, c_r, c_t, c_z, b_r, b_t, b_z, b_xx, b_yy, b_zz
real(DP) :: phi_c, theta, b0, ellipse
integer :: i, j

phi_c = eflux_ph/180.0*PI+PI/2

bx2d(:,:) = 0.0
by2d(:,:) = 0.0
bz2d(:,:) = 0.0
vz2d(:,:) = 0.0

eflux_z = eflux_z0 + eflux_vz*min(eflux_t1-eflux_t0, &
                                  max(0.0,time_eflux-eflux_t0))
!DEBUG
! write(*,*) '## eflux_z=', eflux_z,' atime=',atime,' phi_c=',phi_c
!#####

do i = -1, NX
do j = -1, NY

  ! convert to the rectangle coordinate, in which the origin is located at
  ! the torus center.
  xx = cos(phi_c)*(xc(i)-eflux_x0)+sin(phi_c)*(yc(j)-eflux_y0)
  yy =-sin(phi_c)*(xc(i)-eflux_x0)+cos(phi_c)*(yc(j)-eflux_y0)
  zz = 0.0d0-eflux_z

  ! convert to the cylindrival coordinate
  c_r = sqrt(yy**2+zz**2)
 if(c_r.ne.0.0) then
  c_t = atan2(zz,yy)
 else
  c_t = 0.0
 end if
  c_z = xx

  ellipse = (c_r/eflux_r0)**2 + (c_z/eflux_a0)**2 

  ! emerging velocity

   if(time_eflux.ge.eflux_t0.AND.time_eflux.le.eflux_t1) then
    if(ellipse.le.1.0)  then 

       vz2d(j,i) = eflux_vz

    end if
   end if

  ! magnetic field of emerging flux on the cylindrical coordinate
   if(ellipse.le.1.0.AND.c_r.ne.0.0)  then
      b0 = 1.0d0
    else
      b0 = 0.0
   end if

      b_r = 0.0
      b_t = b0*eflux_b0
      b_z = 0.0

  ! magnetic field on the rectangle coordinate originated at the torus center
  ! theta = atan2(zz,yy)
  b_xx = b_z
  b_yy = b_r*cos(c_t)-b_t*sin(c_t)
  b_zz = b_r*sin(c_t)+b_t*cos(c_t)

  ! magnetic field on the global rectangle coordinate
  bx2d(j,i) = b_xx*cos(phi_c)-b_yy*sin(phi_c)
  by2d(j,i) = b_xx*sin(phi_c)+b_yy*cos(phi_c)
  bz2d(j,i) = b_zz

end do
end do

end subroutine eflux__mag
