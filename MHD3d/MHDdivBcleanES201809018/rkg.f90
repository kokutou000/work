! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
module rkg
! -----------------------------------------------------------
!                    MODULE rkg
! -----------------------------------------------------------
use common
implicit none
! === coefficient parameters for runge-kutta-gill scheme
real(DP) :: c1, c2, c3, c4
real(DP) :: cq1, cq2, cq3, cs2, cs3
integer :: i,j,k

public
contains
!
! -----------------------------------------------------------
subroutine rkg__set_coefficient
!
      c1 = 0.5d0
      c2 = 1.0d0 - sqrt(0.5d0)
      c3 = 1.0d0 + sqrt(0.5d0)
      c4 = 1.0d0 / 6.0d0
      cq1 = -2.0d0
      cq2 =  1.0d0 - 3.0d0 * c2
      cs2 =  2.0d0 * c2
      cq3 =  1.0d0 - 3.0d0 * c3
      cs3 =  2.0d0 * c3

!
return
end subroutine rkg__set_coefficient
! -----------------------------------------------------------
subroutine rkg__prog1
!
       vx(:,:,:) = vx(:,:,:) + c1 * dvx(:,:,:)
       vy(:,:,:) = vy(:,:,:) + c1 * dvy(:,:,:)
       vz(:,:,:) = vz(:,:,:) + c1 * dvz(:,:,:)
       bx(:,:,:) = bx(:,:,:) + c1 * dbx(:,:,:)
       by(:,:,:) = by(:,:,:) + c1 * dby(:,:,:)
       bz(:,:,:) = bz(:,:,:) + c1 * dbz(:,:,:)
       ro(:,:,:) = ro(:,:,:) + c1 * dro(:,:,:)
       pr(:,:,:) = pr(:,:,:) + c1 * dpr(:,:,:)
       phi(:,:,:) = phi(:,:,:) + c1 * dphi(:,:,:)
!
       qvx(:,:,:) = dvx(:,:,:)
       qvy(:,:,:) = dvy(:,:,:)
       qvz(:,:,:) = dvz(:,:,:)
       qbx(:,:,:) = dbx(:,:,:)
       qby(:,:,:) = dby(:,:,:)
       qbz(:,:,:) = dbz(:,:,:)
       qro(:,:,:) = dro(:,:,:)
       qpr(:,:,:) = dpr(:,:,:)
       qphi(:,:,:) = dphi(:,:,:)
!
      return
end subroutine rkg__prog1
! -----------------------------------------------------------
subroutine rkg__prog2
!
       vx(:,:,:) = vx(:,:,:) + c2 * ( dvx(:,:,:) - qvx(:,:,:) )
       vy(:,:,:) = vy(:,:,:) + c2 * ( dvy(:,:,:) - qvy(:,:,:) )
       vz(:,:,:) = vz(:,:,:) + c2 * ( dvz(:,:,:) - qvz(:,:,:) )
       bx(:,:,:) = bx(:,:,:) + c2 * ( dbx(:,:,:) - qbx(:,:,:) )
       by(:,:,:) = by(:,:,:) + c2 * ( dby(:,:,:) - qby(:,:,:) )
       bz(:,:,:) = bz(:,:,:) + c2 * ( dbz(:,:,:) - qbz(:,:,:) )
       ro(:,:,:) = ro(:,:,:) + c2 * ( dro(:,:,:) - qro(:,:,:) )
       pr(:,:,:) = pr(:,:,:) + c2 * ( dpr(:,:,:) - qpr(:,:,:) )
       phi(:,:,:) = phi(:,:,:) + c2 * ( dphi(:,:,:) - qphi(:,:,:) )
!
       qvx(:,:,:) = cq2 * qvx(:,:,:) + cs2 * dvx(:,:,:)
       qvy(:,:,:) = cq2 * qvy(:,:,:) + cs2 * dvy(:,:,:)
       qvz(:,:,:) = cq2 * qvz(:,:,:) + cs2 * dvz(:,:,:)
       qbx(:,:,:) = cq2 * qbx(:,:,:) + cs2 * dbx(:,:,:)
       qby(:,:,:) = cq2 * qby(:,:,:) + cs2 * dby(:,:,:)
       qbz(:,:,:) = cq2 * qbz(:,:,:) + cs2 * dbz(:,:,:)
       qro(:,:,:) = cq2 * qro(:,:,:) + cs2 * dro(:,:,:)
       qpr(:,:,:) = cq2 * qpr(:,:,:) + cs2 * dpr(:,:,:)
       qphi(:,:,:) = cq2 * qphi(:,:,:) + cs2 * dphi(:,:,:)
!
      return
end subroutine rkg__prog2
! -----------------------------------------------------------
subroutine rkg__prog3
!
       vx(:,:,:) = vx(:,:,:) + c3 * ( dvx(:,:,:) - qvx(:,:,:) )
       vy(:,:,:) = vy(:,:,:) + c3 * ( dvy(:,:,:) - qvy(:,:,:) )
       vz(:,:,:) = vz(:,:,:) + c3 * ( dvz(:,:,:) - qvz(:,:,:) )
       bx(:,:,:) = bx(:,:,:) + c3 * ( dbx(:,:,:) - qbx(:,:,:) )
       by(:,:,:) = by(:,:,:) + c3 * ( dby(:,:,:) - qby(:,:,:) )
       bz(:,:,:) = bz(:,:,:) + c3 * ( dbz(:,:,:) - qbz(:,:,:) )
       ro(:,:,:) = ro(:,:,:) + c3 * ( dro(:,:,:) - qro(:,:,:) )
       pr(:,:,:) = pr(:,:,:) + c3 * ( dpr(:,:,:) - qpr(:,:,:) )
       phi(:,:,:) = phi(:,:,:) + c3 * ( dphi(:,:,:) - qphi(:,:,:) )
!
       qvx(:,:,:) = cq3 * qvx(:,:,:) + cs3 * dvx(:,:,:)
       qvy(:,:,:) = cq3 * qvy(:,:,:) + cs3 * dvy(:,:,:)
       qvz(:,:,:) = cq3 * qvz(:,:,:) + cs3 * dvz(:,:,:)
       qbx(:,:,:) = cq3 * qbx(:,:,:) + cs3 * dbx(:,:,:)
       qby(:,:,:) = cq3 * qby(:,:,:) + cs3 * dby(:,:,:)
       qbz(:,:,:) = cq3 * qbz(:,:,:) + cs3 * dbz(:,:,:)
       qro(:,:,:) = cq3 * qro(:,:,:) + cs3 * dro(:,:,:)
       qpr(:,:,:) = cq3 * qpr(:,:,:) + cs3 * dpr(:,:,:)
       qphi(:,:,:) = cq3 * qphi(:,:,:) + cs3 * dphi(:,:,:)
!
      return
end subroutine rkg__prog3
! -----------------------------------------------------------
subroutine rkg__prog4
!
       vx(:,:,:) = vx(:,:,:) + c4 * ( dvx(:,:,:) - 2. * qvx(:,:,:) )
       vy(:,:,:) = vy(:,:,:) + c4 * ( dvy(:,:,:) - 2. * qvy(:,:,:) )
       vz(:,:,:) = vz(:,:,:) + c4 * ( dvz(:,:,:) - 2. * qvz(:,:,:) )
       bx(:,:,:) = bx(:,:,:) + c4 * ( dbx(:,:,:) - 2. * qbx(:,:,:) )
       by(:,:,:) = by(:,:,:) + c4 * ( dby(:,:,:) - 2. * qby(:,:,:) )
       bz(:,:,:) = bz(:,:,:) + c4 * ( dbz(:,:,:) - 2. * qbz(:,:,:) )
       ro(:,:,:) = ro(:,:,:) + c4 * ( dro(:,:,:) - 2. * qro(:,:,:) )
       pr(:,:,:) = pr(:,:,:) + c4 * ( dpr(:,:,:) - 2. * qpr(:,:,:) )
       phi(:,:,:) = phi(:,:,:) + c4 * ( dphi(:,:,:) - 2. * qphi(:,:,:) )
!
      return
end subroutine rkg__prog4

! -----------------------------------------------------------
subroutine rkg__write(message)
use common
character(5),intent(IN) :: message

write(*,*) '======= '//trim(message)//' ========'
write(*,*) 'rkg::',bx(NX/2,NY/2,0:2)
return
end subroutine rkg__write

end module rkg




