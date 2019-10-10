#include "cppdefs.h"
      MODULE wdrift_floats_mod
#if defined FLOAT_OIL && (defined BULK_FLUXES || defined BULK_FLUXES2D)
!
!
!================================================== Dmitry Dukhovskoy  ===
!
!=======================================================================
!                                                                      !
!  These routines compute wind drift linearly added to oil U,V fields  !
!   locations: wnudg                                                   !
!                                                                      !
!  Reference: Samuels et al., 1982 wind drift model with wind-dependt  !
!  deflection angle   + our BOEM report, 2017                          !
!                                                                      !
! Compute wind drift with wind-dependent deflection angle              !
! Wind-drag (Cw) coefficient is fixed at the moment                    ! 
! In the future, needs to be a function of the upper-most              ! 
! layer thickness (inversly proportional), thinner ocean layer         ! 
! has more wind effect already in the ocean current and                ! 
! wind correction has to be reduced - discussed in our                 ! 
! BOEM report (WAMOST 2017, Remote Sensing Assessment of Surface Oil   !
! Transport and Fate during Spills in the Gulf of Mexico               !
!                                                                      ! 
! Oil surface drift is considered as superposition of                  ! 
! advective velocities and turbulent diffusion                         ! 
! Advective velocities are liniear combination of ocean                ! 
! surface currents (ROMS) + wind drift (+ stokes drift - not           ! 
! implemented here but easy to add using the wind)                     ! 
! Vectors:                                                             ! 
! Ua = Uroms + Uwind                                                   ! 
!                                                                      ! 
! Uwind(Udrft, Vdrft) = Cwind * THETA * |W10|                          ! 
!    where Cwind is air drag coefficient                               ! 
!          THETA is a unit vector directed at deflection angle (theta) !
!                from the wind                                         ! 
!          THETA = exp(i*theta)                                        ! 
!          W10 is wind speed at 10 m                                   ! 
!                                                                      ! 
!  Wind-dependent deflection angle (degrees):                          !
!   theta = -25*exp(-1e-8*|W10|/(nu*g))                                !
!   nu - kinematic viscosity of sea water                              !
!                                                                      ! 
! Kinematic viscosity of sea water is a constant here,                 ! 
! possible improvement: compare with the model value                   ! 
!                                                                      ! 
!=======================================================================
!
      implicit none
!
!
      CONTAINS
!***********************************************************************
      SUBROUTINE wdrift_floats (ng, LBi, UBi, LBj, UBj,                 &
     &                          Istr, Iend, Jstr, Jend,                 &
     &                          nfindx, l,                              &
# ifdef OIL_DEBUG
     &                          Uwind, Vwind, track, ldebug)  
# else
     &                          Uwind, Vwind, track)  
# endif
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_grid
      USE mod_scalars
      USE mod_ncparam
      USE nrutil, ONLY : gasdev
!      USE mod_forces

!      USE interp_floats_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
# endif

!
!  Imported variable declarations.
!
!      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: ng, l
      integer, intent(in) :: LBi, UBi, LBj, UBj
!      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Istr, Iend, Jstr, Jend 
      integer, intent(in) :: nfindx
# ifdef OIL_DEBUG
      logical, intent(in) :: ldebug  ! float # to report in log
!      logical, intent(in) :: Predictor
# endif
!
# ifdef ASSUMED_SHAPE
!      logical, intent(in) :: bounded(:)
!      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
!      logical, intent(in) :: bounded(Nfloats(ng))
!      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, icc
      integer :: Ir, Jr, i1, i2, j1, j2

!      real(r8), parameter :: Cwnd=0.035_r8 ! constant for now - make f(layer thickness)
      real(r8), parameter :: nuk=1.05e-6_r8  ! kinem visc of sea water, constant
!      real(r8), parameter :: dfsn=2.0_r8  ! m2/s eddy coeff approx diffusiv random walk

!      real(r8), dimension(Istr:Iend,Jstr:Jend) :: Udrft
!      real(r8), dimension(Istr:Iend,Jstr:Jend) :: Vdrft
!      real(r8) :: wnudg(Lstr:Lend)  ! 

      real(r8) :: Cwnd, dfsn
      real(r8) :: ixold, iyold, ixnew, iynew
      real(r8) :: zoil, w10, theta
      real(r8) :: p1, p2, q1, q2
      real(r8) :: udrflt, vdrflt, udrindx, vdrindx
      real(r8) :: sgm, nrnd, ixdfs, iydfs
      real(r8) :: zb11, zb12, zb21, zb22, zflt0, zflt1
      real(r8) :: rm11, rm12, rm21, rm22
      real(r8) :: uw11, uw12, uw21, uw22
      real(r8) :: vw11, vw12, vw21, vw22
      real(r8) :: udrft11, udrft12, udrft21, udrft22
      real(r8) :: vdrft11, vdrft12, vdrft21, vdrft22

!# include "set_bounds.h"

      Cwnd=CWndOil ! later implement variable coeff 
      dfsn=DfsOil
!
!-----------------------------------------------------------------------
!  Compute wind drift and adjust oil float position
!-----------------------------------------------------------------------
!

! Interpolate wind into float position
      Ir=FLOOR(track(ixgrd,nfindx,l))
      Jr=FLOOR(track(iygrd,nfindx,l))
!
      i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
      i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
      j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
      j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
! wind U component
      uw11=Uwind(i1,j1)
      uw12=Uwind(i1,j2)
      uw21=Uwind(i2,j1)
      uw22=Uwind(i2,j2)
! Wind V component
      vw11=Vwind(i1,j1)
      vw12=Vwind(i1,j2)
      vw21=Vwind(i2,j1)
      vw22=Vwind(i2,j2)
!
! Calculate wind drift at i,j locations around the float
!
      w10=sqrt(uw11*uw11+vw11*vw11)
      theta=-25.0_r8*exp(-1.0e-8_r8*w10**3/(9.8_r8*nuk))
      theta=theta*pi/180.0_r8
      udrft11=Cwnd*(cos(theta)*uw11-sin(theta)*vw11) ! rotate by theta angle and scale by Cwnd     
      vdrft11=Cwnd*(sin(theta)*uw11+cos(theta)*vw11)
      
      w10=sqrt(uw12*uw12+vw12*vw12)
      theta=-25.0_r8*exp(-1.0e-8_r8*w10**3/(9.8_r8*nuk))
      theta=theta*pi/180.0_r8
      udrft12=Cwnd*(cos(theta)*uw12-sin(theta)*vw12) ! rotate by theta angle and scale by Cwnd     
      vdrft12=Cwnd*(sin(theta)*uw12+cos(theta)*vw12)
      
      w10=sqrt(uw21*uw21+vw21*vw21)
      theta=-25.0_r8*exp(-1.0e-8_r8*w10**3/(9.8_r8*nuk))
      theta=theta*pi/180.0_r8
      udrft21=Cwnd*(cos(theta)*uw21-sin(theta)*vw21) ! rotate by theta angle and scale by Cwnd     
      vdrft21=Cwnd*(sin(theta)*uw21+cos(theta)*vw21)
      
      w10=sqrt(uw22*uw22+vw22*vw22)
      theta=-25.0_r8*exp(-1.0e-8_r8*w10**3/(9.8_r8*nuk))
      theta=theta*pi/180.0_r8
      udrft22=Cwnd*(cos(theta)*uw22-sin(theta)*vw22) ! rotate by theta angle and scale by Cwnd     
      vdrft22=Cwnd*(sin(theta)*uw22+cos(theta)*vw22)
      

      p2=REAL(i2-i1,r8)*(track(ixgrd,nfindx,l)-REAL(i1,r8))
      q2=REAL(j2-j1,r8)*(track(iygrd,nfindx,l)-REAL(j1,r8))
      p1=1.0_r8-p2
      q1=1.0_r8-q2

      udrflt=p1*q1*udrft11+p2*q1*udrft21+p1*q2*udrft12+p2*q2*udrft22
      vdrflt=p1*q1*vdrft11+p2*q1*vdrft21+p1*q2*vdrft12+p2*q2*vdrft22
!
!      udrflt=p1*q1*Udrft(i1,j1)+p2*q1*Udrft(i2,j1)+                     &
!     &       p1*q2*Udrft(i1,j2)+p2*q2*Udrft(i2,j2)
!      vdrflt=p1*q1*Vdrft(i1,j1)+p2*q1*Vdrft(i2,j1)+                     &
!     &       p1*q2*Vdrft(i1,j2)+p2*q2*Vdrft(i2,j2)
!
!
! Interpolate depth into float position
! Land mask: make -1 = ocean, +1 = land
!
      rm11=-2.0_r8*GRID(ng)%rmask(i1,j1)+1.0_r8 
      rm12=-2.0_r8*GRID(ng)%rmask(i1,j2)+1.0_r8 
      rm21=-2.0_r8*GRID(ng)%rmask(i2,j1)+1.0_r8 
      rm22=-2.0_r8*GRID(ng)%rmask(i2,j2)+1.0_r8 
!
      zb11=rm11*GRID(ng)%h(i1,j1)
      zb12=rm12*GRID(ng)%h(i1,j2)
      zb21=rm21*GRID(ng)%h(i2,j1)
      zb22=rm22*GRID(ng)%h(i2,j2)
      
      zflt0=p1*q1*zb11+p2*q1*zb21+p1*q2*zb12+p2*q2*zb22

!
! Wind drift in index/sec
!
      udrindx=udrflt*GRID(ng)%pm(i1,j1)
      vdrindx=vdrflt*GRID(ng)%pn(i1,j1)
!
! Add random walk - diffusion
! xt=x(t-1)+u*dt+/-sqrt(2D*dt)*Xr
! yt=y(t-1)+v*dt+/-sqrt(2D*dt)*Yr
! Xr,Yr - random numbr Gauss. (0,1)
! D - diffusion, dt - time scale
! Define random walk as u=U*zt*cos(2*pi*zt),
! v=U*zt*sin(2*pi*zt), zt - normally distributed on [-1,1]
! or simply: u=U*randn; v=V*randn;
! See paper: Lonin, 1999 for example
! 
      sgm=sqrt(2.0_r8*dfsn*dt(ng))
! Keep random number < 3 stdev
      nrnd=0.0_r8
      DO icc=1,1000 
        CALL gasdev(nrnd)
        IF (abs(nrnd).le.3.0_r8) EXIT 
      END DO
      ixdfs=sgm*nrnd*GRID(ng)%pm(i1,j1)

      nrnd=0.0_r8
      DO icc=1,1000 
        CALL gasdev(nrnd)
        IF (abs(nrnd).le.3.0_r8) EXIT 
      END DO
      iydfs=sgm*nrnd*GRID(ng)%pn(i1,j1)

!
! Adjust drift location due to wind drift
! and random walk
!
      ixold=track(ixgrd,nfindx,l)
      iyold=track(iygrd,nfindx,l)
      ixnew=track(ixgrd,nfindx,l)+dt(ng)*udrindx+ixdfs
      iynew=track(iygrd,nfindx,l)+dt(ng)*vdrindx+iydfs
!
! Check if new location is onland 
! if onland - do not move the float to the 
! new position
!
      Ir=FLOOR(ixnew)
      Jr=FLOOR(iynew)
!
      i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
      i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
      j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
      j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
      p2=REAL(i2-i1,r8)*(ixnew-REAL(i1,r8))
      q2=REAL(j2-j1,r8)*(iynew-REAL(j1,r8))
      p1=1.0_r8-p2
      q1=1.0_r8-q2
!
! Interpolate depth into new float position
! Land mask: make -1 = ocean, +1 = land
!
      rm11=-2.0_r8*GRID(ng)%rmask(i1,j1)+1.0_r8 
      rm12=-2.0_r8*GRID(ng)%rmask(i1,j2)+1.0_r8 
      rm21=-2.0_r8*GRID(ng)%rmask(i2,j1)+1.0_r8 
      rm22=-2.0_r8*GRID(ng)%rmask(i2,j2)+1.0_r8 

      zb11=rm11*GRID(ng)%h(i1,j1)
      zb12=rm12*GRID(ng)%h(i1,j2)
      zb21=rm21*GRID(ng)%h(i2,j1)
      zb22=rm22*GRID(ng)%h(i2,j2)

      zflt1=p1*q1*zb11+p2*q1*zb21+p1*q2*zb12+p2*q2*zb22

      IF (zflt1.ge.-1.0e-3_r8) THEN
        ixnew=ixold
        iynew=iyold
      ENDIF

      track(ixgrd,nfindx,l)=ixnew
      track(iygrd,nfindx,l)=iynew

# ifdef OIL_DEBUG
      IF (ldebug) THEN
        print*,' ===============  WIND DRIFT   ========= '
        print*,':: i1=',i1,' j1=',j1,' Ir=',Ir,' Jr=',Jr
        print*,':: Wind-x=',Uwind(i1,j1),' Wind-y=',Vwind(i1,j1)
        print*,':: Udrft=',udrft11,' Vdrft=',vdrft11
        print*,':: Interp to float Udrft=',udrflt,' Vdrft=',vdrflt
        print*,':: Float drift index/sec Uix=',udrindx,' Viy=',vdrindx
        print*,':: RandWalk Difs: dX=',ixdfs,' dY=',iydfs
        print*,':: Old float location i=',ixold,' j=',iyold
        print*,':: Adjusted float location i=',track(ixgrd,nfindx,l),   &
     &           'j=',track(iygrd,nfindx,l)
        print*,':: Depth old location=',zflt0,' new=',zflt1
        print*,'  =============  END WIND DRIFT ========= '
      ENDIF
# endif

      RETURN
      END SUBROUTINE wdrift_floats

#endif
      END MODULE wdrift_floats_mod
