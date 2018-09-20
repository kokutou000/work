! ===========================================================
!     Subroutine for the setting the boundary condition in flux emerging 
!     process NL3DpwE_mpi_02h_nlff by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
!                    SUBROUTINE eflux__mag
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
!!! use test_eflux_subroutines

implicit none

real(DP), dimension(-1:NY,-1:NX), intent(OUT) :: bx2d, by2d, bz2d, vz2d
real(DP) :: time_eflux, eflux_z
real(DP) :: xx, yy, zz, c_r, c_t, c_z, b_r, b_t, b_z, b_xx, b_yy, b_zz
real(DP) :: phi_c, theta, b0
integer :: i, j

phi_c = eflux_ph/180.0*PI+PI/2
!DEBUG
!write(*,*) '## eflux_z=', eflux_z,' atime=',atime,' phi_c=',phi_c
eflux_z = eflux_z0 + eflux_vz*time_eflux

do i = -1, NX
do j = -1, NY

  ! convert to the rectangle coordinate, in which the origin is located at
  ! the torus center.
  xx = cos(phi_c)*(xc(i)-eflux_x0)+sin(phi_c)*(yc(j)-eflux_y0)
  yy =-sin(phi_c)*(xc(i)-eflux_x0)+cos(phi_c)*(yc(j)-eflux_y0)
  zz = 0.0d0-eflux_z

  ! convert to the cylindrival coordinate
  c_r = sqrt(yy**2+zz**2)
  c_t = atan2(zz,yy)
  c_z = xx

  ! emerging velocity
   if(c_r.le.eflux_r0.AND.abs(c_z).le.eflux_a0)  then

      vz2d(j,i) = eflux_vz

   else if(c_r.gt.eflux_r0.AND.c_r.le.eflux_r0*2 &
      .AND.abs(c_z).le.eflux_a0+eflux_a0/eflux_r0*(c_r-eflux_r0))  then

      vz2d(j,i) = eflux_vz*(2*eflux_r0-c_r)/eflux_r0

   else if(abs(c_z).gt.eflux_a0.AND.abs(c_z).le.eflux_a0*2 &
      .AND.c_r.le.eflux_r0+eflux_r0/eflux_a0*(abs(c_z)-eflux_a0)) then

      vz2d(j,i) = eflux_vz*(2*eflux_a0-abs(c_z))/eflux_a0

   else
 
      vz2d(j,i) = 0.0

   end if


  ! magnetic field of emerging flux on the cylindrical coordinate
   if(c_r.le.eflux_r0.AND.abs(c_z).le.eflux_a0)  then
      b0 = eflux_b0*(1.0-(c_z/eflux_a0)**2)
    else
      b0 = 0.0
   end if

      b_r = 0.0
      b_t = b0*eflux_r0/max(c_r,eflux_rm)
      b_z = 0.0

  ! magnetic field on the rectangle coordinate originated at the torus center
  theta = atan2(zz,yy)
  b_xx = b_z
  b_yy = b_r*cos(theta)-b_t*sin(theta)
  b_zz = b_r*sin(theta)+b_t*cos(theta)

  ! magnetic field on the global rectangle coordinate
  bx2d(j,i) = b_xx*cos(phi_c)-b_yy*sin(phi_c)
  by2d(j,i) = b_xx*sin(phi_c)+b_yy*cos(phi_c)
  bz2d(j,i) = b_zz

end do
end do

end subroutine eflux__mag
