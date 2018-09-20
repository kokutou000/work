! ===========================================================
!      3D finite-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
!
! 2005/06/03 debug index in the coefficients of density diffusion
! 2005/06/06 itne initial set for heat flux (qx,qy,qz) added. 
! 2005/06/09 qx,qy,qz are defined at half-integer mesh
! 2007/02/14 eflux included.
! ===========================================================
module mhd
! -----------------------------------------------------------
!                    MODULE mhd
! -----------------------------------------------------------
use common
use mpiut
implicit none
integer :: i,j,k,kp,km,jp,jm,ip,im
integer :: jpp,jmm,kpp,kmm,ipp,imm
integer :: ix0, jy0
real(DP) :: roi
real(DP), dimension(-1:NY,-1:NX) :: eflux_bx2d, eflux_by2d, eflux_bz2d, eflux_vz2d

contains
! -----------------------------------------------------------
subroutine mhd__delta
use common

! ####################################################
!~~~~~~~~ debug 2016/12/26
!if(index_x.eq.0.and.index_y.eq.1) then
!k=1
!j=2
!i=0
!ip=i+1
!im=i-1
!jp=j+1
!jm=j-1
!kp=k+1
!km=k-1
!write(*,*) "+++++++++++"
!write(*,*) vx(k,j,i), vx(k ,j ,ip), d1x(i,+1)
!write(*,*) vx(k ,j ,i ), d1x(i, 0)
!write(*,*) vx(k ,j ,im), d1x(i,-1)
!write(*,*) vy(k,j,i), vx(k ,jp,i ), d1y(j,+1)
!write(*,*) vx(k ,j ,i ), d1y(j, 0)
!write(*,*) vx(k ,jm,i ), d1y(j,-1)
!write(*,*) vz(k,j,i), vx(kp,j ,i ), d1z(k,+1)
!write(*,*) vx(k ,j ,i ), d1z(k, 0)
!write(*,*) vx(km,j ,i ), d1z(k,-1)
!write(*,*) 1.0d0/ro(k,j,i)
!write(*,*) cy(k ,j ,i), bz(k, j, i)
!write(*,*) cz(k ,j ,i), by(k, j, i)
!write(*,*) pr(k ,j ,ip), d1x(i,+1)
!write(*,*) pr(k ,j ,i ), d1x(i, 0)
!write(*,*) pr(k ,j ,im), d1x(i,-1)
!write(*,*) vx(k ,j ,ip), d2x(i, 1)
!write(*,*) vx(k ,j ,im), d2x(i,-1)
!write(*,*) vx(k ,jp,i ), d2y(j, 1)
!write(*,*) vx(k ,jm,i ), d2y(j,-1)
!write(*,*) vx(kp,j ,i ), d2z(k, 1)
!write(*,*) vx(km,j, i ), d2z(k,-1)
!write(*,*) vx(k ,j ,i ), d2xyz0(k,j,i)
!write(*,*) visc
!end if
!~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~ debug 2017/02/02 ~~~~~~~
!if(nloop .eq. 80) then
!write(*,*) "11111"
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) "++++++++++++++"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!end if
!
if(index_x.eq.0) then
   ix0 = 1
else
   ix0 = 0
endif
if(index_y.eq.0) then
   jy0 = 1
else
   jy0 = 0
endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i = ix0, NX-1
         ip=i+1
         im=i-1
      do j = jy0, NY-1
         jp=j+1
         jm=j-1
      do k = 1, NZ-1
         kp=k+1
         km=k-1
         roi = 1.0d0/ro(k,j,i)
! ..................................[VX]
         dvx(k,j,i) = dtstep *                          &
          (-(vx(k,j,i)*(vx(k ,j ,ip)*d1x(i,+1)          &
                       +vx(k ,j ,i )*d1x(i, 0)          &
                       +vx(k ,j ,im)*d1x(i,-1))         &
            +vy(k,j,i)*(vx(k ,jp,i )*d1y(j,+1)          &
                       +vx(k ,j ,i )*d1y(j, 0)          &
                       +vx(k ,jm,i )*d1y(j,-1))         &
            +vz(k,j,i)*(vx(kp,j ,i )*d1z(k,+1)          &
                       +vx(k ,j ,i )*d1z(k, 0)          &
                       +vx(km,j ,i )*d1z(k,-1)))        &
           +roi*                                        &
            ((cy(k ,j ,i)*bz(k, j, i)                   &
             -cz(k ,j ,i)*by(k, j, i))                  &
            -(pr(k ,j ,ip)*d1x(i,+1)                    &
             +pr(k ,j ,i )*d1x(i, 0)                    &
             +pr(k ,j ,im)*d1x(i,-1)))                  &
           +( vx(k ,j ,ip)*d2x(i, 1)                    &
             +vx(k ,j ,im)*d2x(i,-1)                    &
             +vx(k ,jp,i )*d2y(j, 1)                    &
             +vx(k ,jm,i )*d2y(j,-1)                    &
             +vx(kp,j ,i )*d2z(k, 1)                    &
             +vx(km,j, i )*d2z(k,-1)                    &
             +vx(k ,j ,i )*d2xyz0(k,j,i)                &
            )*visc3d(k,j,i))                                     
! ..................................[VY]
         dvy(k,j,i) = dtstep *                          &
          (-(vx(k,j,i)*(vy(k ,j ,ip)*d1x(i,+1)          &
                       +vy(k ,j ,i )*d1x(i, 0)          &
                       +vy(k ,j ,im)*d1x(i,-1))         &
            +vy(k,j,i)*(vy(k ,jp,i )*d1y(j,+1)          &
                       +vy(k ,j ,i )*d1y(j, 0)          &
                       +vy(k ,jm,i )*d1y(j,-1))         &
            +vz(k,j,i)*(vy(kp,j ,i )*d1z(k,+1)          &
                       +vy(k ,j ,i )*d1z(k, 0)          &
                       +vy(km,j ,i )*d1z(k,-1)))        &
           +roi*                                        &
            ((cz(k ,j ,i)*bx(k, j, i)                   &
             -cx(k ,j ,i)*bz(k, j, i))                  &
            -(pr(k ,jp,i )*d1y(j,+1)                    &
             +pr(k ,j ,i )*d1y(j, 0)                    &
             +pr(k ,jm,i )*d1y(j,-1)))                  &
           +( vy(k ,j ,ip)*d2x(i, 1)                    &
             +vy(k ,j ,im)*d2x(i,-1)                    &
             +vy(k ,jp,i )*d2y(j, 1)                    &
             +vy(k ,jm,i )*d2y(j,-1)                    &
             +vy(kp,j ,i )*d2z(k, 1)                    &
             +vy(km,j ,i )*d2z(k,-1)                    &
             +vy(k ,j ,i )*d2xyz0(k,j,i)                &
            )*visc3d(k,j,i))
! ..................................[VZ]
         dvz(k,j,i) = dtstep *                          &
          (-(vx(k,j,i)*(vz(k ,j ,ip)*d1x(i,+1)          &
                       +vz(k ,j ,i )*d1x(i, 0)          &
                       +vz(k ,j ,im)*d1x(i,-1))         &
            +vy(k,j,i)*(vz(k ,jp,i )*d1y(j,+1)          &
                       +vz(k ,j ,i )*d1y(j, 0)          &
                       +vz(k ,jm,i )*d1y(j,-1))         &
            +vz(k,j,i)*(vz(kp,j ,i )*d1z(k,+1)          &
                       +vz(k ,j ,i )*d1z(k, 0)          &
                       +vz(km,j ,i )*d1z(k,-1)))        &
           +roi*                                        &
            ((cx(k ,j ,i )*by(k, j, i)                  &
             -cy(k ,j ,i )*bx(k, j, i))                 &
            -(pr(kp,j ,i )*d1z(k,+1)                    &
             +pr(k ,j ,i )*d1z(k, 0)                    &
             +pr(km,j ,i )*d1z(k,-1)))                  &
           +( vz(k ,j ,ip)*d2x(i, 1)                    &
             +vz(k ,j ,im)*d2x(i,-1)                    &
             +vz(k ,jp,i )*d2y(j, 1)                    &
             +vz(k ,jm,i )*d2y(j,-1)                    &
             +vz(kp,j ,i )*d2z(k, 1)                    &
             +vz(km,j ,i )*d2z(k,-1)                    &
             +vz(k ,j ,i )*d2xyz0(k,j,i)                &
           )*visc3d(k,j,i))
! ..................................[BX]
         dbx(k,j,i) =-dtstep *           &
           ((ez(k ,jp,i )*d1y(j,+1)      &
            +ez(k ,j ,i )*d1y(j, 0)      &
            +ez(k ,jm,i )*d1y(j,-1))     &
           -(ey(kp,j ,i )*d1z(k,+1)      &
            +ey(k ,j ,i )*d1z(k, 0)      &
            +ey(km,j ,i )*d1z(k,-1))     &
           +(phi(k,j,ip)*d1x(i,+1)       &
            +phi(k,j,i )*d1x(i, 0)       &
            +phi(k,j,im)*d1x(i,-1))      &
           )
! ..................................[BY]
         dby(k,j,i) =-dtstep *           &
           ((ex(kp,j ,i )*d1z(k,+1)      &
            +ex(k ,j ,i )*d1z(k, 0)      &
            +ex(km,j ,i )*d1z(k,-1))     &
           -(ez(k ,j ,ip)*d1x(i,+1)      &
            +ez(k ,j ,i )*d1x(i, 0)      &
            +ez(k ,j ,im)*d1x(i,-1))     &
           +(phi(k,jp,i)*d1y(j,+1)       &
            +phi(k,j ,i)*d1y(j, 0)       &
            +phi(k,jm,i)*d1y(j,-1))      &
           )
! ..................................[BZ]
         dbz(k,j,i) =-dtstep *           &
           ((ey(k ,j ,ip)*d1x(i,+1)      &
            +ey(k ,j ,i )*d1x(i, 0)      &
            +ey(k ,j ,im)*d1x(i,-1))     &
           -(ex(k ,jp,i )*d1y(j,+1)      &
            +ex(k ,j ,i )*d1y(j, 0)      &
            +ex(k ,jm,i )*d1y(j,-1))     &
           +(phi(kp,j,i)*d1z(k,+1)       &
            +phi(k ,j,i)*d1z(k, 0)       &
            +phi(km,j,i)*d1z(k,-1))      &
           )
! ..................................[PHI]
         dphi(k,j,i) = -dtstep * ch2 *  &
             (bx(k,j,ip)*d1x(i,+1)      &
             +bx(k,j,i )*d1x(i, 0)      &
             +bx(k,j,im)*d1x(i,-1)      &
             +by(k,jp,i)*d1y(j,+1)      &
             +by(k,j ,i)*d1y(j, 0)      &
             +by(k,jm,i)*d1y(j,-1)      &
             +bz(kp,j,i)*d1z(k,+1)      &
             +bz(k ,j,i)*d1z(k, 0)      &
             +bz(km,j,i)*d1z(k,-1)      &
             +phi(k,j,i)/cp2       &
             )
! ..................................[RO]
         dro(k,j,i) = 0.0       !2008.02.16 (02j)

!!       dro(k,j,i) = dtstep *           & 
!!       (-((rvx(k ,j ,ip)*d1x(i,+1)     & 
!!          +rvx(k ,j ,i )*d1x(i, 0)     &
!!          +rvx(k ,j ,im)*d1x(i,-1))    &
!!         +(rvy(k ,jp,i )*d1y(j,+1)     &
!!          +rvy(k ,j ,i )*d1y(j, 0)     &
!!          +rvy(k ,jm,i )*d1y(j,-1))    &
!!         +(rvz(kp,j ,i )*d1z(k,+1)     &
!!          +rvz(k ,j ,i )*d1z(k, 0)     &
!!          +rvz(km,j ,i )*d1z(k,-1)))   &
!!        +diff_ro*                      &
!!          (ro (k ,j ,ip)*d2x(i,+1)     &
!!          +ro (k ,j ,im)*d2x(i,-1)     &
!!          +ro (k ,jp,i )*d2y(j,+1)     &
!!          +ro (k ,jm,i )*d2y(j,-1)     &
!!          +ro (kp,j ,i )*d2z(k,+1)     &
!!          +ro (km,j ,i )*d2z(k,-1)     &
!!          +ro (k ,j ,i )*d2xyz0(k,j,i)))
! ..................................[PR]
         dpr(k,j,i) = 0.0        !2008.02.16 (02j)

!!       dpr(k,j,i) = dtstep *                     &
!!        (-(vx(k,j,i)*(pr(k ,j ,ip)*d1x(i,+1)     &
!!                     +pr(k ,j ,i )*d1x(i, 0)     &
!!                     +pr(k ,j ,im)*d1x(i,-1))    &
!!          +vy(k,j,i)*(pr(k ,jp,i )*d1y(j,+1)     &
!!                     +pr(k ,j ,i )*d1y(j, 0)     &
!!                     +pr(k ,jm,i )*d1y(j,-1))    &
!!          +vz(k,j,i)*(pr(kp,j ,i )*d1z(k,+1)     &
!!                     +pr(k ,j ,i )*d1z(k, 0)     &
!!                     +pr(km,j ,i )*d1z(k,-1)))   &
!!         -gamma*pr(k,j,i)*                       &
!!                     (vx(k ,j ,ip)*d1x(i,+1)     & 
!!                     +vx(k ,j ,i )*d1x(i, 0)     &
!!                     +vx(k ,j ,im)*d1x(i,-1)     &
!!                     +vy(k ,jp,i )*d1y(j,+1)     &
!!                     +vy(k ,j ,i )*d1y(j, 0)     &
!!                     +vy(k ,jm,i )*d1y(j,-1)     &
!!                     +vz(kp,j ,i )*d1z(k,+1)     &
!!                     +vz(k ,j ,i )*d1z(k, 0)     &
!!                     +vz(km,j ,i )*d1z(k,-1))    &
!!         -gamma1*   ((qx(k ,j ,i )               &
!!                     -qx(k ,j ,im))*ddx2(i)      &
!!                    +(qy(k ,j ,i )               &
!!                     -qy(k ,jm,i ))*ddy2(j)      &
!!                    +(qz(k ,j ,i )               &
!!                     -qz(km,j ,i ))*ddz2(k))     &
!!         +gamma1*res(k,j,i)*iohm_heat*           &
!!                     (cx(k,j,i)**2               &
!!                     +cy(k,j,i)**2               &
!!                     +cz(k,j,i)**2)              &
!!         )
      end do
      end do
      end do
!~~~~~~~~ debug 2016/12/26 ~~~~~~~
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) "???????????"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~ debug 2017/02/02 ~~~~~~~
!if(nloop .eq. 80) then
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "2222222"
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) "++++++++++++++"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!end if
!
! >>>>>>>> TOP & BOTTOM boundary
!
!      do i = 0, nxm1
      do i = ix0, nxm1
         ip=i+1
         im=i-1
!      do j = 0, nym1
      do j = jy0, nym1
         jp=j+1
         jm=j-1
!
         k =0
         kp =k+1
         kpp=k+2
!
         dvx(k,j,i) = 0.0d0
         dvy(k,j,i) = 0.0d0
         dvz(k,j,i) = 0.0d0
!
         dbx(k,j,i) =-dtstep *        &
           ((ez(k ,jp,i )*d1y(j,+1)   &
            +ez(k ,j ,i )*d1y(j, 0)   &
            +ez(k ,jm,i )*d1y(j,-1))  &
           -(ey(kpp,j ,i )*d1z(k,+1)  &
            +ey(kp,j ,i  )*d1z(k, 0)  &
            +ey(k ,j ,i  )*d1z(k,-1)) &
           +(phi(k,j,ip)*d1x(i,+1)       &
            +phi(k,j,i )*d1x(i, 0)       &
            +phi(k,j,im)*d1x(i,-1))      &
           )
!
         dby(k,j,i) =-dtstep *        &
           ((ex(kpp,j ,i )*d1z(k,+1)  &
            +ex(kp,j ,i  )*d1z(k, 0)  &
            +ex(k ,j ,i  )*d1z(k,-1)) &
           -(ez(k ,j ,ip)*d1x(i,+1)   &
            +ez(k ,j ,i )*d1x(i, 0)   &
            +ez(k ,j ,im)*d1x(i,-1))  &
           +(phi(k,jp,i)*d1y(j,+1)       &
            +phi(k,j ,i)*d1y(j, 0)       &
            +phi(k,jm,i)*d1y(j,-1))      &
           )
!
        dbz(k,j,i) =-dtstep *        &
          ((ey(k ,j ,ip)*d1x(i,+1)   &
           +ey(k ,j ,i )*d1x(i, 0)   &
           +ey(k ,j ,im)*d1x(i,-1))  &
          -(ex(k ,jp,i )*d1y(j,+1)   &
           +ex(k ,j ,i )*d1y(j, 0)   &
           +ex(k ,jm,i )*d1y(j,-1))  &
           +(phi(kpp,j,i)*d1z(k,+1)       &
            +phi(kp ,j,i)*d1z(k, 0)       &
            +phi(k  ,j,i)*d1z(k,-1))      &
          )
!
        dphi(k,j,i) = 0d0
        ! dphi(k,j,i) =-((ey(k,j,ip)*d1x(i,+1) &
        !      &        +ey(k,j,i )*d1x(i, 0) &
        !      &        +ey(k,j,im)*d1x(i,-1))&
        !      &       -(ey(k,jp,i)*d1y(j,+1) &
        !      &        +ey(k,j ,i)*d1y(j, 0) &
        !      &        +ey(k,jm,i)*d1y(j,-1))&
        !      &        +phi(kpp,j,i)*d1z(k,+1)&
        !      &        +phi(kp ,j,i)*d1z(k, 0))&
        !      &       /(d1z(k,-1))&
        !      &        -phi(k,j,i)
!
! ####### density fixed #########
         dro(k,j,i) = 0.0
! ###############################
!!$         dro(k,j,i) = dtstep *           & 
!!$         (-( rvz(kpp,j ,i )*d1z(k,+1)    &
!!$            +rvz(kp ,j ,i )*d1z(k, 0)    &
!!$            +rvz(k  ,j ,i )*d1z(k,-1))   &
!!$          +(diff_ro*                     &
!!$            (ro (k ,j ,ip)*d2x(i,+1)     &
!!$            +ro (k ,j ,i )*d2x(i, 0)     &
!!$            +ro (k ,j ,im)*d2x(i,-1)     &
!!$            +ro (k ,jp,i )*d2y(j,+1)     &
!!$            +ro (k ,j ,i )*d2y(j, 0)     &
!!$            +ro (k ,jm,i )*d2y(j,-1)     &
!!$            +ro (kpp,j ,i )*d2z(k,+1)    &
!!$            +ro (kp ,j ,i )*d2z(k, 0)    &
!!$            +ro (k  ,j ,i )*d2z(k,-1))))
!
! ####### pressure fixed #########
         dpr(k,j,i) = 0.0
! ################################
!!$         dpr(k,j,i) = dtstep *                     &
!!$          (-gamma*pr(k,j,i)*                       &
!!$                       (vz(kpp,j ,i )*d1z(k,+1)    &
!!$                       +vz(kp ,j ,i )*d1z(k, 0)    &
!!$                       +vz(k  ,j ,i )*d1z(k,-1))   &
!!$           -gamma1*   ((qx(k ,j ,i )               &
!!$                       -qx(k ,j ,im))*ddx2(i)      &
!!$                      +(qy(k ,j ,i )               &
!!$                       -qy(k ,jm,i ))*ddy2(j)      &
!!$                      +(qz(k ,j ,i )               &
!!$                                    )*ddz2(k)))
! ################################
!
         k =nz
         km =k-1
         kmm=k-2
!
         dvx(k,j,i) = 0.0d0
         dvy(k,j,i) = 0.0d0
         dvz(k,j,i) = 0.0d0
!
         dbx(k,j,i) =-dtstep *        &
           ((ez(k ,jp,i )*d1y(j,+1)   &
            +ez(k ,j ,i )*d1y(j, 0)   &
            +ez(k ,jm,i )*d1y(j,-1))  &
           -(ey(k ,j ,i )*d1z(k,+1)   &
            +ey(km,j ,i )*d1z(k, 0)   &
            +ey(kmm,j ,i)*d1z(k,-1))  &
           +(phi(k,j,ip)*d1x(i,+1)       &
            +phi(k,j,i )*d1x(i, 0)       &
            +phi(k,j,im)*d1x(i,-1))      &
           )
!
         dby(k,j,i) =-dtstep *        &
           ((ex(k ,j ,i )*d1z(k,+1)   &
            +ex(km,j ,i )*d1z(k, 0)   &
            +ex(kmm,j,i )*d1z(k,-1))  &
           -(ez(k ,j ,ip)*d1x(i,+1)   &
            +ez(k ,j ,i )*d1x(i, 0)   &
            +ez(k ,j ,im)*d1x(i,-1))  &
           +(phi(k,jp,i)*d1y(j,+1)       &
            +phi(k,j ,i)*d1y(j, 0)       &
            +phi(k,jm,i)*d1y(j,-1))      &
           )
!
         dbz(k,j,i) =-dtstep *        &
           ((ey(k ,j ,ip)*d1x(i,+1)   &
            +ey(k ,j ,i )*d1x(i, 0)   &
            +ey(k ,j ,im)*d1x(i,-1))  &
           -(ex(k ,jp,i )*d1y(j,+1)   &
            +ex(k ,j ,i )*d1y(j, 0)   &
            +ex(k ,jm,i )*d1y(j,-1))  &
           +(phi(k  ,j,i)*d1z(k,+1)       &
            +phi(km ,j,i)*d1z(k, 0)       &
            +phi(kmm,j,i)*d1z(k,-1))      &
           )
!
         dphi(k,j,i) = 0d0
        ! dphi(k,j,i) =-((ey(k,j,ip)*d1x(i,+1) &
        !      &        +ey(k,j,i )*d1x(i, 0) &
        !      &        +ey(k,j,im)*d1x(i,-1))&
        !      &       -(ey(k,jp,i)*d1y(j,+1) &
        !      &        +ey(k,j ,i)*d1y(j, 0) &
        !      &        +ey(k,jm,i)*d1y(j,-1))&
        !      &        +phi(km ,j,i)*d1z(k, 0) &
        !      &        +phi(kmm,j,i)*d1z(k,-1))&
        !      &       /(d1z(k,+1))&
        !      &        -phi(k,j,i)
!
! ####### density fixed #########
         dro(k,j,i) = 0.0
! ###############################
!!$         dro(k,j,i) = dtstep *           & 
!!$          (-(rvz(k  ,j ,i )*d1z(k,+1)    &
!!$            +rvz(km ,j ,i )*d1z(k, 0)    &
!!$            +rvz(kmm,j ,i )*d1z(k,-1))   &
!!$           +diff_ro*                     &
!!$            (ro (k ,j ,ip)*d2x(i,+1)     &
!!$            +ro (k ,j ,i )*d2x(i, 0)     &
!!$            +ro (k ,j ,im)*d2x(i,-1)     &
!!$            +ro (k ,jp,i )*d2y(j,+1)     &
!!$            +ro (k ,j ,i )*d2y(j, 0)     &
!!$            +ro (k ,jm,i )*d2y(j,-1)     &
!!$            +ro (k  ,j ,i )*d2z(k,+1)    &
!!$            +ro (km ,j ,i )*d2z(k, 0)    &
!!$            +ro (kmm,j ,i )*d2z(k,-1)))
!
! ####### pressure fixed #########
         dpr(k,j,i) = 0.0
! ################################
!!$         dpr(k,j,i) = dtstep *                     &
!!$          (-gamma*pr(k,j,i)*                       &
!!$                       (vz(k  ,j ,i )*d1z(k,+1)    &
!!$                       +vz(km ,j ,i )*d1z(k, 0)    &
!!$                       +vz(kmm,j ,i )*d1z(k,-1))   &
!!$           -gamma1*   ((qx(k ,j ,i )               &
!!$                       -qx(k ,j ,im))*ddx2(i)      &
!!$                      +(qy(k ,j ,i )               &
!!$                       -qy(k ,jm,i ))*ddy2(j)      &
!!$                      +(                           &
!!$                       -qz(km,j ,i ))*ddz2(k)))
! ################################
!
      end do
      end do
!
!~~~~~~~~ debug 2016/12/26 ~~~~~~~
!if(index_x == 0 .and.index_y==0) then
!write(*,*) '**1', index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) '++2', index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) '--3', index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) '==4', index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!end if
!if(nloop .eq. 80) then
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "333333333"
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) "++++++++++++++"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ########################################
! ##### data exchage for Y direction #####
! ########################################

      if(nproc_y /= 0) then
       call mpiut__exchange_d_y
      end if

! ########################################
!~~~~~~~~ debug 2016/12/26 ~~~~~~~
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) ":::::::::::::::::"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!if(nloop .eq. 80) then
!if(index_x == 0 .and.index_y==0) then
!write(*,*) "444444444"
!write(*,*) "**1", index_y, dvx(1,nx-1,-1), dvx(1,nx-1,0), dvx(1,nx-1,1)
!write(*,*) "++2", index_y, dvx(1,nx,-1), dvx(1,nx,0), dvx(1,nx,1)
!write(*,*) "++++++++++++++"
!end if
!if(index_x == 0 .and.index_y==1) then
!write(*,*) "--3", index_y, dvx(1,-1,-1), dvx(1,-1,0), dvx(1,-1,1)
!write(*,*) "==4", index_y, dvx(1,0,-1), dvx(1,0,0), dvx(1,0,1)
!write(*,*) "~~~~~~~~~~~~~~~~~"
!end if
!end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ###### y=0 boundry #####
!

   if(index_y == 0) then

         j = 0
!         j = -1
         jp=j+1
         jpp=j+2

!      do i = 1, nxm1
      do i = ix0, nxm1
         ip=i+1
         im=i-1
      do k = 1, nzm1
         kp=k+1
         km=k-1

!
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
!
         dbx(k,j,i) =-dtstep *        &
           ((ez(k ,jpp,i)*d1y_0(+1)   &
            +ez(k ,jp ,i)*d1y_0( 0)   &
            +ez(k ,j  ,i)*d1y_0(-1))  &
           -(ey(kp,j ,i )*d1z(k,+1)   &
            +ey(k ,j ,i )*d1z(k, 0)   &
            +ey(km,j ,i )*d1z(k,-1))  &
           +(phi(k,j,ip)*d1x(i,+1)       &
            +phi(k,j,i )*d1x(i, 0)       &
            +phi(k,j,im)*d1x(i,-1))      &
           )
!
!         dby(k,j,i) = 0.0
         dby(k,j,i) =-dtstep *           &
           ((ex(kp,j ,i )*d1z(k,+1)      &
            +ex(k ,j ,i )*d1z(k, 0)      &
            +ex(km,j ,i )*d1z(k,-1))     &
           -(ez(k ,j ,ip)*d1x(i,+1)      &
            +ez(k ,j ,i )*d1x(i, 0)      &
            +ez(k ,j ,im)*d1x(i,-1))     &
           +(phi(k,jpp,i)*d1y_0(+1)       &
            +phi(k,jp ,i)*d1y_0( 0)       &
            +phi(k,j  ,i)*d1y_0(-1))      &
           )
!
         dbz(k,j,i) =-dtstep *        &
           ((ey(k,j ,ip)*d1x(i,+1)    &
            +ey(k,j ,i )*d1x(i, 0)    &
            +ey(k,j ,im)*d1x(i,-1))   &
           -(ex(k,jpp,i)*d1y_0(+1)    &
            +ex(k,jp ,i)*d1y_0( 0)    &
            +ex(k,j  ,i)*d1y_0(-1))   &
           +(phi(kp,j,i)*d1z(k,+1)       &
            +phi(k ,j,i)*d1z(k, 0)       &
            +phi(km,j,i)*d1z(k,-1))      &
           )

         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0

         dphi(k,j,i) = 0d0
         ! dphi(k,j,i) =-((ex(kp,j ,i )*d1z(k,+1)&
         !     & +ex(k ,j ,i )*d1z(k, 0)      &
         !     & +ex(km,j ,i )*d1z(k,-1))     &
         !     & -(ez(k ,j ,ip)*d1x(i,+1)      &
         !     & +ez(k ,j ,i )*d1x(i, 0)      &
         !     & +ez(k ,j ,im)*d1x(i,-1))     &
         !     & +(phi(k,jpp,i)*d1y_0(+1)      &
         !     & +phi(k,jp ,i)*d1y_0( 0))     &
         !     & )/(d1y_0(-1)) &
         !     & - phi(k,j,i)

      end do
      end do

   end if
!
! ###### y=Ly boundry #####
!
   if(index_y == nproc_y-1) then

         j  = ny
!         j  =ny-1
         jm =j-1
         jmm=j-2
!
!      do i = 1, nxm1
      do i = ix0, nxm1
         ip=i+1
         im=i-1
      do k = 1, nzm1
         kp=k+1
         km=k-1

         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
!
         dbx(k,j,i) =-dtstep *        &
           ((ez(k,j  ,i )*d1y_ny(+1)   &
            +ez(k,jm ,i )*d1y_ny( 0)   &
            +ez(k,jmm,i )*d1y_ny(-1))  &
           -(ey(kp,j ,i )*d1z(k,+1)   &
            +ey(k ,j ,i )*d1z(k, 0)   &
            +ey(km,j ,i )*d1z(k,-1))  &
           +(phi(k,j,ip)*d1x(i,+1)       &
            +phi(k,j,i )*d1x(i, 0)       &
            +phi(k,j,im)*d1x(i,-1))      &
           )        
!
!         dby(k,j,i) = 0.0
         dby(k,j,i) =-dtstep *           &
           ((ex(kp,j ,i )*d1z(k,+1)      &
            +ex(k ,j ,i )*d1z(k, 0)      &
            +ex(km,j ,i )*d1z(k,-1))     &
           -(ez(k ,j ,ip)*d1x(i,+1)      &
            +ez(k ,j ,i )*d1x(i, 0)      &
            +ez(k ,j ,im)*d1x(i,-1))     &
           +(phi(k,j  ,i)*d1y_ny(+1)       &
            +phi(k,jm ,i)*d1y_ny( 0)       &
            +phi(k,jmm,i)*d1y_ny(-1))      &
           )
!
         dbz(k,j,i) =-dtstep *        &
           ((ey(k ,j ,ip)*d1x(i,+1)   &
            +ey(k ,j ,i )*d1x(i, 0)   &
            +ey(k ,j ,im)*d1x(i,-1))  &
           -(ex(k ,j ,i )*d1y_ny(+1)   &
            +ex(k ,jm,i )*d1y_ny( 0)   &
            +ex(k ,jmm,i)*d1y_ny(-1))  &
           +(phi(kp,j,i)*d1z(k,+1)       &
            +phi(k ,j,i)*d1z(k, 0)       &
            +phi(km,j,i)*d1z(k,-1))      &
           )

         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0

         dphi(k,j,i) = 0d0
         ! dphi(k,j,i) =-((ex(kp,j ,i )*d1z(k,+1)      &
         !      & +ex(k ,j ,i )*d1z(k, 0)      &
         !      & +ex(km,j ,i )*d1z(k,-1))     &
         !      & -(ez(k ,j ,ip)*d1x(i,+1)      &
         !      & +ez(k ,j ,i )*d1x(i, 0)      &
         !      & +ez(k ,j ,im)*d1x(i,-1))     &
         !      & +phi(k,jm ,i)*d1y_ny( 0)       &
         !      & +phi(k,jmm,i)*d1y_ny(-1)      &
         !      & )/(d1y_ny(+1)) &
         !      & - phi(k,j,i)

       end do
       end do

   end if
!
!######### boundary condition for x #########
!
!      if(nproc_x == 1) then
!         call mhd__period(dvx)
!         call mhd__period(dvy)
!         call mhd__period(dvz)
!         call mhd__period(dbx)
!         call mhd__period(dby)
!         call mhd__period(dbz)
!         call mhd__period(dro)
!         call mhd__period(dpr)
!       else
!         call mpiut__exchange_d_x
!      end if
!!
!
! ########################################
! ##### data exchage for X direction #####
! ########################################

      if(nproc_x /= 0) then
       call mpiut__exchange_d_x
      end if

! ########################################
!
! ###### x=0 boundry #####
!

   if(index_x == 0) then

!         j = 0
!         jp=j+1
!         jpp=j+2
      i = 0
!      i = -1
      ip=i+1
      ipp=i+2
      
!      do i = 0, nxm1
!         ip=i+1
!         im=i-1
!      do j = 1, nym1
      do j = jy0, nym1
         jp=j+1
         jm=j-1
      do k = 1, nzm1
         kp=k+1
         km=k-1

!
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
!
!         dbx(k,j,i) = 0.0         
         dbx(k,j,i) =-dtstep *           &
           ((ez(k ,jp,i )*d1y(j,+1)      &
            +ez(k ,j ,i )*d1y(j, 0)      &
            +ez(k ,jm,i )*d1y(j,-1))     &
           -(ey(kp,j ,i )*d1z(k,+1)      &
            +ey(k ,j ,i )*d1z(k, 0)      &
            +ey(km,j ,i )*d1z(k,-1))     &
           +(phi(k,j,ipp)*d1x_0(+1)      &
            +phi(k,j,ip )*d1x_0( 0)      &
            +phi(k,j,i  )*d1x_0(-1))     &
            )
!!
!            dby(k,j,i) = 0.0
         dby(k,j,i) =-dtstep *        &
           ((ex(kp,j ,i )*d1z(k,+1)   &
            +ex(k ,j ,i )*d1z(k, 0)   &
            +ex(km,j ,i )*d1z(k,-1))  &
           -(ez(k ,j,ipp)*d1x_0(+1)   &
            +ez(k ,j ,ip)*d1x_0( 0)   &
            +ez(k ,j ,i )*d1x_0(-1))  &
           +(phi(k,jp,i)*d1y(j,+1)       &
            +phi(k,j ,i)*d1y(j, 0)       &
            +phi(k,jm,i)*d1y(j,-1))      &
           )
!!
         dbz(k,j,i) =-dtstep *        &
!           ((ey(k,j ,ip)*d1x(i,+1)    &
!            +ey(k,j ,i )*d1x(i, 0)    &
!            +ey(k,j ,im)*d1x(i,-1))   &
!           -(ex(k,jpp,i)*d1y_0(+1)    &
!            +ex(k,jp ,i)*d1y_0( 0)    &
!            +ex(k,j  ,i)*d1y_0(-1))   &
           ((ey(k,j,ipp)*d1x_0(+1)    &
            +ey(k,j ,ip)*d1x_0( 0)    &
            +ey(k,j  ,i)*d1x_0(-1))   &
           -(ex(k,jp ,i)*d1y(j,+1)    &
            +ex(k,j  ,i)*d1y(j, 0)    &
            +ex(k,jm ,i)*d1y(j,-1))   &         
           +(phi(kp,j,i)*d1z(k,+1)       &
            +phi(k ,j,i)*d1z(k, 0)       &
            +phi(km,j,i)*d1z(k,-1))      &
           )

         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0

         dphi(k,j,i) = 0d0
         ! dphi(k,j,i) =-((ez(k ,jp,i )*d1y(j,+1)      &
         !      & +ez(k ,j ,i )*d1y(j, 0)      &
         !      & +ez(k ,jm,i )*d1y(j,-1))     &
         !      & -(ey(kp,j ,i )*d1z(k,+1)      &
         !      & +ey(k ,j ,i )*d1z(k, 0)      &
         !      & +ey(km,j ,i )*d1z(k,-1))     &
         !      & +(phi(k,j,ipp)*d1x_0(+1)      &
         !      & +phi(k,j,ip )*d1x_0( 0))      &
         !      & )/(d1x_0(-1))&
         !      & - phi(k,j,i)

      end do
      end do

   end if
!
! ###### x=Lx boundry #####
!
   if(index_x == nproc_x-1) then

!         j  = ny
!         jm =j-1
!         jmm=j-2
         i  = nx
!         i  = nx-1
         im =i-1
         imm=i-2
!!
      do j = jy0, nym1
         jp=j+1
         jm=j-1
      do k = 1, nzm1
         kp=k+1
         km=k-1

         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
!!
!         dbx(k,j,i) = 0.0
         dbx(k,j,i) =-dtstep *           &
           ((ez(k ,jp,i )*d1y(j,+1)      &
            +ez(k ,j ,i )*d1y(j, 0)      &
            +ez(k ,jm,i )*d1y(j,-1))     &
           -(ey(kp,j ,i )*d1z(k,+1)      &
            +ey(k ,j ,i )*d1z(k, 0)      &
            +ey(km,j ,i )*d1z(k,-1))     &
           +(phi(k,j,i  )*d1x_nx(+1)      &
            +phi(k,j,im )*d1x_nx( 0)      &
            +phi(k,j,imm)*d1x_nx(-1))     &
            )
!!
!         dby(k,j,i) = 0.0
         dby(k,j,i) =-dtstep *        &
           ((ex(kp,j ,i )*d1z(k,+1)   &
            +ex(k ,j ,i )*d1z(k, 0)   &
            +ex(km,j ,i )*d1z(k,-1))  &
           -(ez(k ,j ,i )*d1x_nx(+1)   &
            +ez(k ,j ,im)*d1x_nx( 0)   &
            +ez(k ,j,imm)*d1x_nx(-1))  &
           +(phi(k,jp,i)*d1y(j,+1)       &
            +phi(k,j ,i)*d1y(j, 0)       &
            +phi(k,jm,i)*d1y(j,-1))      &
           )
!!
!         dbz(k,j,i) =-dtstep *        &
!           ((ey(k ,j ,ip)*d1x(i,+1)   &
!            +ey(k ,j ,i )*d1x(i, 0)   &
!            +ey(k ,j ,im)*d1x(i,-1))  &
!           -(ex(k ,j ,i )*d1y_ny(+1)   &
!            +ex(k ,jm,i )*d1y_ny( 0)   &
!            +ex(k ,jmm,i)*d1y_ny(-1))  &
!           )
         dbz(k,j,i) =-dtstep *        &
           ((ey(k ,j  ,i)*d1x_nx(+1)   &
            +ey(k ,j ,im)*d1x_nx( 0)   &
            +ey(k ,j,imm)*d1x_nx(-1))  &
           -(ex(k ,jp,i )*d1y(j,+1)   &
            +ex(k ,j ,i )*d1y(j, 0)   &
            +ex(k ,jm,i )*d1y(j,-1))  &
           +(phi(kp,j,i)*d1z(k,+1)       &
            +phi(k ,j,i)*d1z(k, 0)       &
            +phi(km,j,i)*d1z(k,-1))      &
           )
!!
         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0

         dphi(k,j,i) = 0d0
         ! dphi(k,j,i) = -((ez(k ,jp,i )*d1y(j,+1)      &
         !      & +ez(k ,j ,i )*d1y(j, 0)      &
         !      & +ez(k ,jm,i )*d1y(j,-1))     &
         !      & -(ey(kp,j ,i )*d1z(k,+1)      &
         !      & +ey(k ,j ,i )*d1z(k, 0)      &
         !      & +ey(km,j ,i )*d1z(k,-1))     &
         !      & +phi(k,j,im )*d1x_nx( 0)      &
         !      & +phi(k,j,imm)*d1x_nx(-1)     &
         !      & )/(d1x_nx(+1))&
         !      & - phi(k,j,i)

       end do
       end do

   end if
!
!
   if((index_x.eq.0).and.(index_y.eq.0)) then
      i = 0
      j = 0
      do k = 0, nz
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
         dbx(k,j,i) = 0.0
         dby(k,j,i) = 0.0
         dbz(k,j,i) = 0.0
         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0
         dphi(k,j,i) = 0.0
      enddo
   end if
   if((index_x.eq.0).and.(index_y.eq.nproc_y-1)) then
      i = 0; j = ny
      do k = 0, nz
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
         dbx(k,j,i) = 0.0
         dby(k,j,i) = 0.0
         dbz(k,j,i) = 0.0
         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0
         dphi(k,j,i) = 0.0
      enddo
   end if
   if((index_x.eq.nproc_x-1).and.(index_y.eq.0)) then
      i = nx; j = 0
      do k = 0, nz
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
         dbx(k,j,i) = 0.0
         dby(k,j,i) = 0.0
         dbz(k,j,i) = 0.0
         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0
         dphi(k,j,i) = 0.0
      enddo
   end if
   if((index_x.eq.nproc_x-1).and.(index_y.eq.nproc_y-1)) then
      i = nx; j = ny
      do k = 0, nz
         dvx(k,j,i) = 0.0
         dvy(k,j,i) = 0.0
         dvz(k,j,i) = 0.0
         dbx(k,j,i) = 0.0
         dby(k,j,i) = 0.0
         dbz(k,j,i) = 0.0
         dro(k,j,i) = 0.0
         dpr(k,j,i) = 0.0
         dphi(k,j,i) = 0.0
      enddo
   end if
!
!
     return
end subroutine mhd__delta
! -----------------------------------------------------------
! -----------------------------------------------------------
! -----------------------------------------------------------
! -----------------------------------------------------------
! -----------------------------------------------------------
subroutine mhd__sub
! 2005/06/02 limiter introduced
! 2005/06/09 qx,qy,qz defined at half integer mech
!
!   |-------|-------|-------|-------| <= major variables defined
!   i-2     i-1     i       i+1     i+2
!       +       +       +       +     <= q[x,y,z] defined
!       i-2     i-1     i       i+1   index for q[x,y,z]
!   qx is defined at the point (i+1/2,j,k)
!   qy is defined at the point (i,j+1/2,k)
!   qz is defined at the point (i,j,k+1/2)
!   f.g. qx(1,1,1) defind at x=(xc(1)+xc(2))/2,y=yc(1),z=zc(1)
!   f.g. qy(1,1,1) defind at x=xc(1),y=(yc(1)+yc(2))/2,z=zc(1)
!   f.g. qz(1,1,1) defind at x=xc(1),y=yc(1),z=(zc(1)+zc(1))/2
!
use common
integer :: j3d, i3d
integer :: ix0, jy0
!
! ++++ emerging flux ++++
!
!     call eflux__mag(eflux_bx2d,eflux_by2d,eflux_bz2d,eflux_vz2d,atime)
!     vx(0,:,:) = 0.0
!     vy(0,:,:) = 0.0
!!     vz(0,:,:) = eflux_vz2d(:,:)
!     vz(0,:,:) = 0.0
! === copy to 3d vector for lateral condition on z=0 for velocity
! === magnetic field is not fixed condition!! 
! === but almost fixed because dB/dt~0
  do i = -1, NX
     i3d = index_x*NX + i
  do j = -1, NY
     j3d = index_y*NY + j
!!     if((i3d.eq.-1).or.(j3d.eq.-1)) then
!!        bx(0,j,i) = 0d0
!!        by(0,j,i) = 0d0
!!        bz(0,j,i) = 0d0
!!        vx(0,j,i) = 0d0
!!        vy(0,j,i) = 0d0
!!     else
!!        bx(0,j,i) = initbx3d(i3d,j3d,0)
!!        by(0,j,i) = initby3d(i3d,j3d,0)
!!        bz(0,j,i) = initbz3d(i3d,j3d,0)
     ! if(atime.le.time_ramp_init) then
     !    vx(0,j,i) = initvx2d(i3d,j3d)*(atime/time_ramp_init)
     !    vy(0,j,i) = initvy2d(i3d,j3d)*(atime/time_ramp_init)
     ! else if( (atime.gt.time_ramp_init) &
     !   &.and. (atime.le.time_end_tw) ) then
     !    vx(0,j,i) = initvx2d(i3d,j3d)
     !    vy(0,j,i) = initvy2d(i3d,j3d)
     ! else if( (atime.gt.time_end_tw) &
     !   &.and. (atime.le.(time_end_tw + time_ramp_end)) ) then
     !    vx(0,j,i) = initvx2d(i3d,j3d) &
     !         & * (1d0+(-atime+time_end_tw)/time_ramp_init)
     !    vy(0,j,i) = initvy2d(i3d,j3d) &
     !         & * (1d0+(-atime+time_end_tw)/time_ramp_init)
     ! else
        vx(0,j,i) = 0d0
        vy(0,j,i) = 0d0
     ! end if
!!     end if
  end do
  end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!--- rho = B*B for constant alfven speed --- 20161222
  ro(:,:,:) = sqrt(bx(:,:,:)**2+by(:,:,:)**2+bz(:,:,:)**2)
!--- to avoid NaN, upper and side density make dense --- 20170829
 ! do k = 0, NZ
 ! do j = 0, NY
 ! do i = 0, NX
 !    if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. &
 !         & (abs(zc(k)).le.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i)*(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) + &
 !            & sqrt(ro(k,j,i))*(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)
 !    else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
 !         & (abs(zc(k)).le.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i)*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) + &
 !            & sqrt(ro(k,j,i))*(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)
 !    else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. & 
 !         & (abs(zc(k)).ge.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i)*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) + &
 !            & sqrt(ro(k,j,i))*(zc(k)-0.9*zl)/(0.1*zl)
 !    else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
 !         & (abs(zc(k)).le.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i)* &
 !            & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) + &
 !            & sqrt(ro(k,j,i))* (1.0-&
 !            & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)))
 !    else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).le.0.4d0*yl).and. & 
 !         & (abs(zc(k)).ge.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i)* &
 !            & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) + &
 !            & sqrt(ro(k,j,i))* (1.0-&
 !            & (1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
 !    else if((abs(xc(i)-0.5*xl).le.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
 !         & (abs(zc(k)).ge.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i) &
 !            &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) &
 !            &+sqrt(ro(k,j,i)) * (1.0 &
 !            &-(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
 !    else if((abs(xc(i)-0.5*xl).ge.0.4d0*xl).and. &
 !         & (abs(yc(j)-0.5*yl).ge.0.4d0*yl).and. & 
 !         & (abs(zc(k)).ge.0.9d0*zl)) then
 !       ro(k,j,i) = ro(k,j,i) &
 !            &*(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)) &
 !            &+sqrt(ro(k,j,i)) * (1.0 &
 !            &-(1.0-(abs(xc(i)-0.5*xl)-0.4*xl)/(0.1*xl)) &
 !            &*(1.0-(abs(yc(j)-0.5*yl)-0.4*yl)/(0.1*yl)) &
 !            &*(1.0-(zc(k)-0.9*zl)/(0.1*zl)))
 !    end if
 ! end do
 ! end do
 ! end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!--- upper boundary condition for velocity -- open boundary
!
  ! k = nz
  ! vx(k,:,:) = (vx(k-1,:,:)*(zc(k)-zc(k-2))**2  &
  !           & -vx(k-2,:,:)*(zc(k)-zc(k-1))**2) &
  !           & /( (zc(k)-zc(k-2))**2 - (zc(k)-zc(k-1))**2 )
  ! vy(k,:,:) = (vy(k-1,:,:)*(zc(k)-zc(k-2))**2  &
  !           & -vy(k-2,:,:)*(zc(k)-zc(k-1))**2) &
  !           & /( (zc(k)-zc(k-2))**2 - (zc(k)-zc(k-1))**2 )
  ! vz(k,:,:) = (vz(k-1,:,:)*(zc(k)-zc(k-2))**2  &
  !           & -vz(k-2,:,:)*(zc(k)-zc(k-1))**2) &
  !           & /( (zc(k)-zc(k-2))**2 - (zc(k)-zc(k-1))**2 )
!
!--- lateral boundary condition for velocity -- open boundary
!
     ! if(index_x == 0) then
     !  vx(:,:,0) = (vx(:,:,1)*(xc(0)-xc(2))**2  &
     !            & -vx(:,:,2)*(xc(0)-xc(1))**2) &
     !            & /( (xc(0)-xc(2))**2 - (xc(0)-xc(1))**2 )
     !  vy(:,:,0) = (vy(:,:,1)*(xc(0)-xc(2))**2  &
     !            & -vy(:,:,2)*(xc(0)-xc(1))**2) &
     !            & /( (xc(0)-xc(2))**2 - (xc(0)-xc(1))**2 )
     !  vz(:,:,0) = (vz(:,:,1)*(xc(0)-xc(2))**2  &
     !            & -vz(:,:,2)*(xc(0)-xc(1))**2) &
     !            & /( (xc(0)-xc(2))**2 - (xc(0)-xc(1))**2 )
     !  vx(:,:,-1) = vx(:,:,0)
     !  vy(:,:,-1) = vy(:,:,0)
     !  vz(:,:,-1) = vz(:,:,0)
     ! end if

     ! if(index_x == nproc_x-1) then
     !   vx(:,:,nx) = (vx(:,:,nx-1)*(xc(nx)-xc(nx-2))**2  &
     !              & -vx(:,:,nx-2)*(xc(nx)-xc(nx-1))**2) &
     !              & /( (xc(nx)-xc(nx-2))**2 - (xc(nx)-xc(nx-1))**2 )
     !   vy(:,:,nx) = (vy(:,:,nx-1)*(xc(nx)-xc(nx-2))**2  &
     !              & -vy(:,:,nx-2)*(xc(nx)-xc(nx-1))**2) &
     !              & /( (xc(nx)-xc(nx-2))**2 - (xc(nx)-xc(nx-1))**2 )
     !   vz(:,:,nx) = (vz(:,:,nx-1)*(xc(nx)-xc(nx-2))**2  &
     !              & -vz(:,:,nx-2)*(xc(nx)-xc(nx-1))**2) &
     !              & /( (xc(nx)-xc(nx-2))**2 - (xc(nx)-xc(nx-1))**2 )
     ! end if
!
     ! if(index_y == 0) then
     !  vx(:,0,:) = (vx(:,1,:)*(yc(0)-yc(2))**2  &
     !            & -vx(:,2,:)*(yc(0)-yc(1))**2) &
     !            & /( (yc(0)-yc(2))**2 - (yc(0)-yc(1))**2 )
     !  vy(:,0,:) = (vy(:,1,:)*(yc(0)-yc(2))**2  &
     !            & -vy(:,2,:)*(yc(0)-yc(1))**2) &
     !            & /( (yc(0)-yc(2))**2 - (yc(0)-yc(1))**2 )
     !  vz(:,0,:) = (vz(:,1,:)*(yc(0)-yc(2))**2  &
     !            & -vz(:,2,:)*(yc(0)-yc(1))**2) &
     !            & /( (yc(0)-yc(2))**2 - (yc(0)-yc(1))**2 )
     !  vx(:,-1,:) = vx(:,0,:)
     !  vy(:,-1,:) = vy(:,0,:)
     !  vz(:,-1,:) = vz(:,0,:)
     ! end if

     ! if(index_y== nproc_y-1) then
     !   vx(:,nx,:) = (vx(:,nx-1,:)*(yc(nx)-yc(nx-2))**2  &
     !              & -vx(:,nx-2,:)*(yc(nx)-yc(nx-1))**2) &
     !              & /( (yc(nx)-yc(nx-2))**2 - (yc(nx)-yc(nx-1))**2 )
     !   vy(:,nx,:) = (vy(:,nx-1,:)*(yc(nx)-yc(nx-2))**2  &
     !              & -vy(:,nx-2,:)*(yc(nx)-yc(nx-1))**2) &
     !              & /( (yc(nx)-yc(nx-2))**2 - (yc(nx)-yc(nx-1))**2 )
     !   vz(:,nx,:) = (vz(:,nx-1,:)*(yc(nx)-yc(nx-2))**2  &
     !              & -vz(:,nx-2,:)*(yc(nx)-yc(nx-1))**2) &
     !              & /( (yc(nx)-yc(nx-2))**2 - (yc(nx)-yc(nx-1))**2 )
     ! end if

!--- now emerging flux is assumed zero
     eflux_bx2d=0d0
     eflux_by2d=0d0
     eflux_bz2d=0d0
!
! ++++ limiter for density and pressure ++++
!
      ro(:,:,:) = max(ro(:,:,:),ro_min)
      pr(:,:,:) = max(pr(:,:,:),pr_min)
!
! ++++ momentum (ro*v) ++++

      rvx(:,:,:) = ro(:,:,:)*vx(:,:,:)
      rvy(:,:,:) = ro(:,:,:)*vy(:,:,:)
      rvz(:,:,:) = ro(:,:,:)*vz(:,:,:)

! ++++ temperature ++++

      tmp(:,:,:) = pr(:,:,:)/ro(:,:,:)
!
! ++++ heat flux ++++

      do i = 0, nxm1
         ip = i+1
         im = i-1
      do j = 0, nym1
         jp = j+1
         jm = j-1
      do k = 0, nzm1
         kp = k+1
         km = k-1
!
         qx(k,j,i) =-akappa*(tmp(k,j,ip)-tmp(k,j,i))*ddx(i)
         qy(k,j,i) =-akappa*(tmp(k,jp,i)-tmp(k,j,i))*ddy(j)
         qz(k,j,i) =-akappa*(tmp(kp,j,i)-tmp(k,j,i))*ddz(k)
!
      end do
      end do
      end do
!
! ++++ current ++++
!
if(index_x.eq.0) then
   ix0 = 1
else
   ix0 = 0
endif
if(index_y.eq.0) then
   jy0 = 1
else
   jy0 = 0
endif
!
      do i = ix0, nxm1
         ip = i+1
         im = i-1
      do j = jy0, nym1
         jp = j+1
         jm = j-1
      do k = 1, nzm1
         kp = k+1
         km = k-1
!
         cx(k,j,i) = (bz(k ,jp,i )*d1y(j,+1)  &
                     +bz(k, j ,i )*d1y(j, 0)  &
                     +bz(k ,jm,i )*d1y(j,-1)) &
                    -(by(kp,j ,i )*d1z(k,+1)  &
                     +by(k ,j ,i )*d1z(k, 0)  &
                     +by(km,j ,i )*d1z(k,-1))  

         cy(k,j,i) = (bx(kp,j ,i )*d1z(k,+1)  &
                     +bx(k ,j ,i )*d1z(k, 0)  &
                     +bx(km,j ,i )*d1z(k,-1)) &
                    -(bz(k ,j ,ip)*d1x(i,+1)  &
                     +bz(k ,j ,i )*d1x(i, 0)  &
                     +bz(k ,j ,im)*d1x(i,-1))  

         cz(k,j,i) = (by(k ,j ,ip)*d1x(i,+1)  &
                     +by(k ,j ,i )*d1x(i, 0)  &
                     +by(k ,j ,im)*d1x(i,-1)) &
                    -(bx(k ,jp,i )*d1y(j,+1)  &
                     +bx(k ,j ,i )*d1y(j, 0)  &
                     +bx(k ,jm,i )*d1y(j,-1))  
!
      end do
      end do
      end do
!
! +++ electric current at the top & bottom boundary +++

      do i = ix0, nxm1
         ip = i+1
         im = i-1
      do j = jy0, nym1
         jp = j+1
         jm = j-1
!
         k  = 0
         kp = k+1
         kpp= k+2
!
         cx(k,j,i) = (bz(k ,jp,i )*d1y(j,+1)  &
                     +bz(k , j,i )*d1y(j, 0)  &
                     +bz(k ,jm,i )*d1y(j,-1)) &
                    -(by(kpp,j,i )*d1z(k,+1)  &
                     +by(kp, j,i )*d1z(k, 0)  &
                     +by(k , j,i )*d1z(k,-1))  

         cy(k,j,i) = (bx(kpp,j,i )*d1z(k,+1)  &
                     +bx(kp,j ,i )*d1z(k, 0)  &
                     +bx(k ,j ,i )*d1z(k,-1)) &
                    -(bz(k ,j ,ip)*d1x(i,+1)  &
                     +bz(k ,j ,i )*d1x(i, 0)  &
                     +bz(k ,j ,im)*d1x(i,-1))   

         cz(k,j,i) = (by(k ,j ,ip)*d1x(i,+1)  &
                     +by(k ,j ,i )*d1x(i, 0)  &
                     +by(k ,j ,im)*d1x(i,-1)) &
                    -(bx(k ,jp,i )*d1y(j,+1)  &
                     +bx(k ,j ,i )*d1y(j, 0)  &
                     +bx(k ,jm,i )*d1y(j,-1))   
!
         k = nz
         km = k-1
         kmm= k-2
!
         cx(k,j,i) = (bz(k ,jp,i )*d1y(j,+1)  &
                     +bz(k, j, i )*d1y(j, 0)  &
                     +bz(k ,jm,i )*d1y(j,-1)) &
                    -(by(k ,j ,i )*d1z(k,+1)  &
                     +by(km,j ,i )*d1z(k, 0)  &
                     +by(kmm,j,i )*d1z(k,-1))  

         cy(k,j,i) = (bx(k ,j ,i )*d1z(k,+1)  &
                     +bx(km,j ,i )*d1z(k, 0)  &
                     +bx(kmm,j,i )*d1z(k,-1)) &
                    -(bz(k ,j ,ip)*d1x(i,+1)  &
                     +bz(k ,j ,i )*d1x(i, 0)  &
                     +bz(k ,j ,im)*d1x(i,-1))  

         cz(k,j,i) = (by(k ,j ,ip)*d1x(i,+1)  &
                     +by(k ,j ,i )*d1x(i, 0)  &
                     +by(k ,j ,im)*d1x(i,-1)) &
                    -(bx(k ,jp,i )*d1y(j,+1)  &
                     +bx(k, j ,i )*d1y(j, 0)  &
                     +bx(k ,jm,i )*d1y(j,-1))
!
      end do
      end do
!
! ########################################
! ##### data exchage for Y direction #####
! ########################################

     if(nproc_y /= 1) then
        call mpiut__exchange_sub_y
     end if

! ########################################
! ##### data exchage for X direction #####
! ########################################

!     if(nproc_x == 1) then
!        call mhd__period(cx)
!        call mhd__period(cy)
!        call mhd__period(cz)
!        call mhd__period(qx)
!        call mhd__period(qy)
!        call mhd__period(qz)
!      else
!        call mpiut__exchange_sub_x
!     end if
!~~~~~~~~~~~~~~~~~ test 20161225
     if(nproc_x /= 1) then
        call mpiut__exchange_sub_x
     end if
!~~~~~~~~~~~~~~~~~
! -------------------------------------------------------------
!
! ++++ electric field ++++
!
       ex(:,:,:) = res(:,:,:)*cx(:,:,:) &
                  -(vy(:,:,:)*bz(:,:,:)-vz(:,:,:)*by(:,:,:))
       ey(:,:,:) = res(:,:,:)*cy(:,:,:) &
                  -(vz(:,:,:)*bx(:,:,:)-vx(:,:,:)*bz(:,:,:))
       ez(:,:,:) = res(:,:,:)*cz(:,:,:) &
                  -(vx(:,:,:)*by(:,:,:)-vy(:,:,:)*bx(:,:,:))
!
! ++++ emerging flux ++++
!
!       ex(0,:,:) = -(vy(0,:,:)*eflux_bz2d(:,:) &
!                    -vz(0,:,:)*eflux_by2d(:,:))
!       ey(0,:,:) = -(vz(0,:,:)*eflux_bx2d(:,:) &
!                    -vx(0,:,:)*eflux_bz2d(:,:))
!       ez(0,:,:) = -(vx(0,:,:)*eflux_by2d(:,:) &
!                    -vy(0,:,:)*eflux_bx2d(:,:))
!
! ++++ lateral boundary condition ++++

     if(index_y == 0) then
      ex(:,0,:) = 0.0
      ey(:,0,:) = 0.0
      ez(:,0,:) = 0.0
       ex(:,-1,:) = 0.0
       ey(:,-1,:) = 0.0
       ez(:,-1,:) = 0.0
     end if

     if(index_y == nproc_y-1) then
       ex(:,NY,:) = 0.0
       ey(:,NY,:) = 0.0
       ez(:,NY,:) = 0.0
     end if
!
!--- for debug
!     if(index_y == nproc_y-1) then
!       ex(:,NY-1,:) = 0.0
!       ey(:,NY-1,:) = 0.0
!       ez(:,NY-1,:) = 0.0
!     end if
! ++++ lateral boundary condition for x ++++

     if(index_x == 0) then
      ex(:,:,0) = 0.0
      ey(:,:,0) = 0.0
      ez(:,:,0) = 0.0
       ex(:,:,-1) = 0.0
       ey(:,:,-1) = 0.0
       ez(:,:,-1) = 0.0
     end if

     if(index_x == nproc_x-1) then
       ex(:,:,NX) = 0.0
       ey(:,:,NX) = 0.0
       ez(:,:,NX) = 0.0
     end if
!
!--- for debug
!     if(index_x == nproc_x-1) then
!       ex(:,:,NX-1) = 0.0
!       ey(:,:,NX-1) = 0.0
!       ez(:,:,NX-1) = 0.0
!     end if
!-----------------------------
!--- bottom and top boundary condition for phi
     do i = ix0, nxm1
        ip = i+1; im = i-1
     do j = jy0, nym1
        jp = j+1; jm = j-1
!
        k = 0; kp = k+1; kpp = k+2
        phi(k,j,i) =-((ey(k,j,ip)*d1x(i,+1) &
             &        +ey(k,j,i )*d1x(i, 0) &
             &        +ey(k,j,im)*d1x(i,-1))&
             &       -(ex(k,jp,i)*d1y(j,+1) &
             &        +ex(k,j ,i)*d1y(j, 0) &
             &        +ex(k,jm,i)*d1y(j,-1))&
             &        +phi(kpp,j,i)*d1z(k,+1)&
             &        +phi(kp ,j,i)*d1z(k, 0))&
             &       /(d1z(k,-1))
        k = nz; km = k-1; kmm = k-2
        phi(k,j,i) =-((ey(k,j,ip)*d1x(i,+1) &
             &        +ey(k,j,i )*d1x(i, 0) &
             &        +ey(k,j,im)*d1x(i,-1))&
             &       -(ex(k,jp,i)*d1y(j,+1) &
             &        +ex(k,j ,i)*d1y(j, 0) &
             &        +ex(k,jm,i)*d1y(j,-1))&
             &        +phi(km ,j,i)*d1z(k, 0) &
             &        +phi(kmm,j,i)*d1z(k,-1))&
             &       /(d1z(k,+1))
     end do
     end do
!--- side boundary (y) condition for phi
     if(index_y == 0) then
     j = 0; jp = j+1; jpp = j+2
     do i = ix0, nxm1
        ip = i+1; im = i-1
     do k = 1, nzm1
        kp = k+1; km = k-1
        phi(k,j,i) =-          &
           ((ex(kp,j ,i )*d1z(k,+1)      &
            +ex(k ,j ,i )*d1z(k, 0)      &
            +ex(km,j ,i )*d1z(k,-1))     &
           -(ez(k ,j ,ip)*d1x(i,+1)      &
            +ez(k ,j ,i )*d1x(i, 0)      &
            +ez(k ,j ,im)*d1x(i,-1))     &
            +phi(k,jpp,i)*d1y_0(+1)       &
            +phi(k,jp ,i)*d1y_0( 0)       &
           )/(d1y_0(-1))
     enddo
     enddo
     end if
!
     if(index_y == nproc_y - 1) then
     j = ny; jm = j-1; jmm = j-2
     do i = ix0, nxm1
        ip = i+1; im = i-1
     do k = 1, nzm1
        kp = k+1; km = k-1
        phi(k,j,i) =-           &
           ((ex(kp,j ,i )*d1z(k,+1)      &
            +ex(k ,j ,i )*d1z(k, 0)      &
            +ex(km,j ,i )*d1z(k,-1))     &
           -(ez(k ,j ,ip)*d1x(i,+1)      &
            +ez(k ,j ,i )*d1x(i, 0)      &
            +ez(k ,j ,im)*d1x(i,-1))     &
            +phi(k,jm ,i)*d1y_ny( 0)       &
            +phi(k,jmm,i)*d1y_ny(-1)      &
           )/(d1y_ny(+1))
     enddo
     enddo
     end if
!--- side boundary (x) condition for phi
     if(index_x == 0) then
     i = 0; ip = i+1; ipp = i+2
     do j = jy0, nym1
        jp = j+1; jm = j-1
     do k = 1, nzm1
        kp = k+1; km = k-1
        phi(k,j,i) =-           &
           ((ez(k ,jp,i )*d1y(j,+1)      &
            +ez(k ,j ,i )*d1y(j, 0)      &
            +ez(k ,jm,i )*d1y(j,-1))     &
           -(ey(kp,j ,i )*d1z(k,+1)      &
            +ey(k ,j ,i )*d1z(k, 0)      &
            +ey(km,j ,i )*d1z(k,-1))     &
            +phi(k,j,ipp)*d1x_0(+1)      &
            +phi(k,j,ip )*d1x_0( 0)      &
            )/(d1x_0(-1))
     enddo
     enddo
     endif
!
     if(index_x == nproc_x - 1) then
     i = nx; im = i-1; imm = i-2
     do j = jy0, nym1
        jp = j+1; jm = j-1
     do k = 1, nzm1
        kp = k+1; km = k-1
        phi(k,j,i) =-           &
           ((ez(k ,jp,i )*d1y(j,+1)      &
            +ez(k ,j ,i )*d1y(j, 0)      &
            +ez(k ,jm,i )*d1y(j,-1))     &
           -(ey(kp,j ,i )*d1z(k,+1)      &
            +ey(k ,j ,i )*d1z(k, 0)      &
            +ey(km,j ,i )*d1z(k,-1))     &
            +phi(k,j,im )*d1x_nx( 0)      &
            +phi(k,j,imm)*d1x_nx(-1)     &
            )/(d1x_nx(+1))
     enddo
     enddo
     end if
!-------------------
!-------------------
     if(nproc_y /= 1) then
        call mpiut__exchange_phi_y
     end if
     if(nproc_x /= 1) then
        call mpiut__exchange_phi_x
     end if
!-------------------
      return
end subroutine mhd__sub
! -----------------------------------------------------------
subroutine mhd__period(a)
use constants
real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(INOUT) :: a
!
         a(:,:,NX) = a(:,:,0)
         a(:,:,-1) = a(:,:,NX-1)
!
      return
end subroutine mhd__period
! -----------------------------------------------------------
subroutine mhd__forcing
! add the external forcing
! introduced at NL3DpwE_mpi_02g_nlff on 20071122

   dvx(1,:,:) = dvx(1,:,:)+dtstep*force_x(:,:)
   dvy(1,:,:) = dvy(1,:,:)+dtstep*force_y(:,:)
   dvz(1,:,:) = dvz(1,:,:)+dtstep*force_z(:,:)

end subroutine mhd__forcing

end module mhd






