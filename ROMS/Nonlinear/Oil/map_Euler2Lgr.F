#include "cppdefs.h"

      MODULE euler2lagr_mod
#if defined FLOAT_OIL && defined OIL_EULR
!==============================================  Dmitry S Dukhovskoy ===
!                                                                      !
!=======================================================================
!                                                                      !
!  This routine maps oil variables from Eulerian -> Lagrangian coord   !
!  Tis is needed after oil concentration field has been changed by     !
!  biological/ sedimentation processes                                 !
!                                                                      !
!  It is assumed that biodegr. process changes oil component structure !
!  because aromatics and saturates are consumed at a faster rate       !
!  asphaltines and raisins consumed less, thus oil partilce becomes    !
!  heavier and smaller. Need to find these changes                     !
!                                                                      !
!  I'm using the following approach:                                   !
! For every grid cell with oil                                         !
! update oil weight fraction for Nocmp components -> this changes oil  !
!                                                    particle density  !
! During bio degradation, oil particle typically becomes               !
!    heavier as lighter  components got consumed faster                !
!                                                                      !
! update oil size (Doil) : this is done such that the overall Mass of  !
! Oil in grid cell                                                     !
! with updated oil density  = Coil(new)*Vgrid_cell                     !
!                                                                      !
! Looking for a size reduction coefficient alfa                        !
! voil_prtcl(new)=voil_prtcl(old)*alfa                                 !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: euler2lagr

      CONTAINS

!
!***********************************************************************
      SUBROUTINE euler2lagr(ng, tile, Lstr, Lend, ntt)
!***********************************************************************
! ntt - time level
!
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_grid
!      USE mod_stepping
      USE mod_oil_Eulvar
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstr, Lend, ntt
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL euler2lagr_tile(ng, tile, Lstr, Lend,                        &
     &                     LBi, UBi, LBj, UBj, ntt,                     &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     OIL3D(ng) % Coil,                            &
     &                     OIL3D(ng) % Coil0,                           &
     &                     OIL3D(ng) % NFLT3D,                          &
     &                     DRIFTER(ng) % track,                         &
     &                     DRIFTER(ng) % bounded,                       &
     &                     DRIFTER(ng) % wfroil0,                       &
     &                     DRIFTER(ng) % szoil0,                        &
     &                     GRID(ng) % omn,                              &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)

      RETURN
      END SUBROUTINE euler2lagr

!
!***********************************************************************
      SUBROUTINE euler2lagr_tile(ng, tile, Lstr, Lend,                  &
     &                           LBi, UBi, LBj, UBj, ntt,               &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Coil1, Coil0, NFLT3D, track,           &
     &                           bounded, wfroil0, szoil0, omn,         &
     &                           Hz, z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_stepping
      USE oil_floats_mod, ONLY : oil_density
!      USE map_Lgr2Euler,   ONLY: lagr2euler_tile

      implicit none

!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile, Lstr, Lend
      integer, intent(in) :: LBi, UBi, LBj, UBj, ntt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: track(:,0:,:)
      real(r8), intent(inout) :: wfroil0(:,:)
      real(r8), intent(inout) :: szoil0(:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)

      logical, intent(in) :: bounded(:)

      integer, intent(in) :: NFLT3D(LBi:,LBj:,:)
      real(r8), intent(in) :: Coil1(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Coil0(LBi:,LBj:,:,:)
#else
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
      real(r8), intent(inout) :: wfroil0(Nfloats(ng),Nocmp)
      real(r8), intent(inout) :: szoil0(Nfloats(ng))
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))

      logical, intent(in) ::  bounded(Nfloats(ng))

      integer, intent(in) ::  NFLT3D(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Coil1(LBi:UBi,LBj:UBj,N(ng),Nocmp)
      real(r8), intent(in) :: Coil0(LBi:,LBj:,:,:)
#endif

!
!  Local variables
!
      integer :: i, j, k, ic, l, nfgrd, il, nff1
      integer :: Lindx(Lstr:Lend,3)         ! grid cell indices where floats falls in 
      integer :: lfltgrd(Lstr-Lend+1)       ! # of floats in 1 grid cell
!
      real(r8) :: Zstr, Zend, Xstr, Xend, Ystr, Yend     
      real(r8) :: xflt, yflt, zflt, rhoo0, rhoo1
      real(r8) :: frsat, frarom, rr, wsum
      real(r8) :: dmm, alfa, vo0, vo1, dmm1

      real(r8) :: GrdVol(IminS:ImaxS,JminS:JmaxS,N(ng))
!      real(r8) :: MoilGrd(IminS:ImaxS,JminS:JmaxS,N(ng),Nocmp)
      real(r8) :: rCoil(Nocmp)                ! fraction of oil reduced by biodegr / sedimentation
      real(r8) :: wfr0(Nocmp), wfr1(Nocmp)    ! old/new weight fractions
      real(r8) :: Moil0, Moil1                ! oil mass in 1 grid cell
      real(r8) :: moil_prt0                   ! oil mass in 1 oil particle
      real(r8) :: Nprt_flt          ! # of oil prt in grd cell: Moil=sum(flts)[rhoo*voil*Nprtcl]
      real(r8) :: dsum
      real(r8) :: wfr_min(Nocmp), wfr_max(Nocmp) ! for checking
!      real(r8) :: tsrfo
      real(r8), dimension(Lstr:Lend) :: Voil0, Voil1, Doil0, Doil1,     &
     &                                  Rhoil0, Rhoil1
      
      logical, dimension(Lstr:Lend) :: my_thread
      logical :: lchck

#ifdef OIL_MAP_DEBUG
      integer :: ifchck, cc
      integer :: nrcoil(Nocmp), brcoil(Nocmp)
      real(r8) :: mnrcoil(Nocmp), mxrcoil(Nocmp)
#endif

#include "set_bounds.h"
!
      DO k=1,N(ng)
        GrdVol(IminS:ImaxS,JminS:JmaxS,k)=                              &
     &    omn(IminS:ImaxS,JminS:JmaxS)*Hz(IminS:ImaxS,JminS:JmaxS,k)
      ENDDO

! Floats only that belong to the thread:
      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8

!
!     Mapping Eulerian -> Lagrangian
!
      
!
!  Collocate floats with Eulerian grid cells
!
      Lindx=0
      DO l=Lstr, Lend
        my_thread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,ntt,l)).and.                           &
     &      (track(ixgrd,ntt,l).lt.Xend).and.                           &
     &      (Ystr.le.track(iygrd,ntt,l)).and.                           &
     &      (track(iygrd,ntt,l).lt.Yend)) THEN
          my_thread(l)=.TRUE.
        ELSE
          CYCLE
        ENDIF

        IF (my_thread(l).and.bounded(l)) THEN
          xflt=track(ixgrd,ntt,l)
          yflt=track(iygrd,ntt,l)
          zflt=track(izgrd,ntt,l)
!
          i=nint(xflt)
          j=nint(yflt)
          k=nint(zflt)
          IF (k.lt.1) k=1
          IF (k.gt.N(ng)) k=N(ng)
          IF (i.lt.BOUNDS(ng)%Istr(MyRank)) i=BOUNDS(ng)%Istr(MyRank)
          IF (i.gt.BOUNDS(ng)%Iend(MyRank)) i=BOUNDS(ng)%Iend(MyRank)
          IF (j.lt.BOUNDS(ng)%Jstr(MyRank)) j=BOUNDS(ng)%Jstr(MyRank)
          IF (j.gt.BOUNDS(ng)%Jend(MyRank)) j=BOUNDS(ng)%Jend(MyRank)
          Lindx(l,1)=i
          Lindx(l,2)=j
          Lindx(l,3)=k
        ENDIF
      ENDDO
!
!  For every grid cell with oil:
!
#ifdef OIL_MAP_DEBUG
      ifchck=4
      lchck=.TRUE.
      nrcoil=0
      brcoil=0
      mnrcoil=1.0_r8
      mxrcoil=0.0_r8
#endif
      wfr_min=1.e6_r8
      wfr_max=-1.0_r8
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            nfgrd=NFLT3D(i,j,k)
            IF (nfgrd.eq.0) CYCLE
! Check
            dsum = 0.0_r8
            DO ic=1,Nocmp
              dsum = dsum+Coil1(i,j,k,ic)
            ENDDO
            IF (dsum.lt.1.e-30_r8) THEN
              print*,'map_Eulr2Lgr: ??? ERR: 0 oil grid conc Coil1??', &
     &        'i,j,k=',i,j,k,' Coil1=',Coil1(i,j,k,:),  &
     &        ' Coil0=',Coil0(i,j,k,:),' Nflts in grid=',nfgrd
            ENDIF
!
! Find component reduction during bio/sedimentation processes
!
            DO ic=1,Nocmp
              rCoil(ic)=Coil1(i,j,k,ic)/Coil0(i,j,k,ic)
#ifdef OIL_MAP_DEBUG
              IF (abs(rCoil(ic)-1.0_r8).gt.1.0e-15_r8) THEN
                nrcoil(ic)=nrcoil(ic)+1
                IF (rCoil(ic).lt.mnrcoil(ic)) mnrcoil(ic)=rCoil(ic)
                IF (rCoil(ic).gt.mxrcoil(ic)) mxrcoil(ic)=rCoil(ic)
              ENDIF
              IF ((rCoil(ic)-1.0_r8).gt.1.0e-15) THEN
                brcoil(ic)=brcoil(ic)+1
              ENDIF
#endif
            ENDDO
!
! Old/New mass of oil in the grid cell
! 
            Moil0=0.0_r8
            Moil1=0.0_r8
            DO ic=1,Nocmp
              Moil0=Moil0+Coil0(i,j,k,ic)*GrdVol(i,j,k)
              Moil1=Moil1+Coil1(i,j,k,ic)*GrdVol(i,j,k)
            ENDDO
!
! Find all floats in the grid cell & check with NFLT
! passed from the Lagr2eulr mapping
!
            nff1=0
            lfltgrd=0
            DO l=Lstr,Lend
              IF (my_thread(l).and.bounded(l)) THEN
                IF (Lindx(l,1).eq.i .and.                               &
     &              Lindx(l,2).eq.j .and.                               &
     &              Lindx(l,3).eq.k) THEN
                  nff1=nff1+1
                  lfltgrd(nff1)=l
                ENDIF
              ENDIF
            ENDDO

            IF (nff1.ne.nfgrd) THEN
              print*,'*** ERRR: map_eulr2lagr: Inconsistent # floats' 
              print*,'*** ERR: i=',i,'j=',j,'k=',k,                     &
     &               'eulr2lagr  N=',nff1,'lagr2eulr N=',nfgrd
            ENDIF
!
!  Create new/old arrays 
!
            Doil0=0.0_r8
            Doil1=0.0_r8
            Rhoil0=0.0_r8
            Rhoil1=0.0_r8
            Nprt_flt=0.0_r8            
            DO il=1,nfgrd          
              l=lfltgrd(il)
              Doil0(l)=track(isizo,ntt,l)
              Rhoil0(l)=track(iroil,ntt,l)
            ENDDO
!
! For each float in the grid cell
! Update info only until oil float hits the surface
!           
            DO il=1,nfgrd
              l=lfltgrd(il)             
!              tsrfo=track(isrfo,ntt,l)
!
! (1) Update oil weight fraction after different oil components
! have been consumed in biology or deposited in sediment modules
! Use the constrain that sum(w(ic))=1 to figure out all weights
!
!     Weight fraction of oil components SAR+A:
              wfr0(1)=track(isats,ntt,l)
              wfr0(2)=track(iarom,ntt,l)
              wfr0(3)=1.0_r8-(wfr0(1)+wfr0(2))

              wsum=0.0_r8
              rr=0.0_r8
              DO ic=1,Nocmp
                rr=rCoil(ic)
                wfr1(ic)=rr*wfr0(ic)
                wsum=wsum+wfr1(ic)
              ENDDO

              IF (wsum.le.1.e-10) THEN
                print*,'map_Euler2Lagr: **ERR: ZERO wsum=',wsum, &
     &          ' rCoil=',rCoil,  &
     &          ' wfr0=',wfr0,' l=',l
              ENDIF
 
!
              wfr1=wfr1/wsum

              wsum=0.0_r8
              DO ic=1,Nocmp
                wsum=wsum+wfr1(ic)
              ENDDO
          
              IF (abs(wsum-1.0_r8)>1.e-10_r8) THEN
                print*,' *** ERR: map_eulr2lagr total weight is not 1'
              ENDIF

! Checking weight fractions of oil components
              DO ic=1,Nocmp
                dmm1=wfr_min(ic)
                dmm1=min(dmm1,wfr1(ic))
                wfr_min(ic)=dmm1

                dmm1=wfr_max(ic)
                dmm1=max(dmm1,wfr1(ic))
                wfr_max(ic)=dmm1
              ENDDO
# ifdef OIL_MAP_DEBUG
              IF (lchck .and. il.eq.ifchck) THEN
                print*,'eulr2lgr: i,j,k=',i,j,k,'float=',l,             &
     &                 'N flts in grid=',nfgrd,                         &
     &                 ' reductn rCoil=',rCoil,                         &
     &                 'old oil w.fract:',wfr0,                         &
     &                 ' new oil w.fract:',wfr1
              ENDIF
# endif
! Update oil component weight fraction:
!              IF (tsrfo<1.e-3_r8) THEN
              track(isats,ntt,l)=wfr1(1)
              track(iarom,ntt,l)=wfr1(2)
!              ENDIF
!
! Update initial oil component mass fraction 
! for evaporation module
              DO ic=1,Nocmp
                wfroil0(l,ic)=wfr1(ic)
              ENDDO
!
! (2) New Oil particle density using updated weight fraction:
! assuming all particles within the grid cell have
! equal reduction of oil components (given by new weight fractions)
! RhoOilComp - densities of the oil components, input info
! Update oil particle density for each float in the grid cell
!              dmm=sum(wfr1/RhoOilComp)
!              Rhoil1(l)=1.0_r8/dmm           ! mean oil particle density, in the float
               CALL oil_density(Rhoil1(l),wfr1,RhoOilComp)  ! mean oil particle density, in the float

              IF (Rhoil1(l).lt.500.0_r8) THEN
                print*,'eulr2lgr: ** ERR Rho<500, RhoOilComp=',         &
     &                 RhoOilComp
                print*,'eulr2lgr: ** ERR old rho=',Rhoil0(l),           &
     &                 ' newrho=',Rhoil1(l),'wfr1=',wfr1
              ENDIF
# ifdef OIL_MAP_DEBUG
              IF (lchck .and. il.eq.ifchck) THEN
!                print*,'eulr2lgr: Rho Oil Comp=',RhoOilComp
                print*,'eulr2lgr: old rho=',Rhoil0(l),                  &
     &                 ' new rho=',Rhoil1(l),' Mass fractn wfr=',wfr1
              ENDIF
# endif
            ENDDO                            ! floats in the grid cell

! (3) Find new (reduced) oil size
! Doil(new) = alfa*Doil(old)              
! Using the overall mass of oil (in the grid cell) M(new)
! sum(rhooil(new)*voil(new)*Nprtcl)=M(ned)
! Nprtcl - # of oil particles in 1 float
! this does not change (at least during biodegradation over 1 step)
! voil=4/3*pi*(Doil/2)**3
! Then alfa= M(new)/sum(rhooil(new)*voil(old)*Nprtcl)
! Note that oil mass within the grid cell can be found as
! vol. integral of Coil1(i,j,k) dv = Moil1  and
! sum(floats)[rho_oil(float)*vol(float)*Nprtcl] = Moil1
!
! Find Nprtcl in the float
            dmm=0.0_r8
            DO il=1,nfgrd
              l=lfltgrd(il)
              Voil0(l)=4.0_r8/3.0_r8*pi*(Doil0(l)/2.0_r8)**3  ! oil particle volume
              moil_prt0=Rhoil0(l)*Voil0(l)
              dmm=dmm+moil_prt0
            ENDDO
            Nprt_flt=Moil0/dmm

!# ifdef OIL_MAP_DEBUG
!            IF (lchck) THEN
!              DO ic=1,3
!                print*,'eulr2lgr: Coil0=',Coil0(i,j,k,ic),              & 
!     &                 'Coil1=',Coil1(i,j,k,ic)
!              ENDDO
!              print*,'Volume Grid cell=',GrdVol(i,j,k)
!              print*,'eulr2lgr: In grid cell Mass oil=',Moil0,          &
!     &                'Mass oil floats',dmm,' Nmb prtlcs=',Nprt_flt
!            ENDIF
!# endif

! 
! Calculate sum(rho(new)*voil(old)*Nprtcl - summation over all floats in
! the grid cell
            dmm=0.0_r8
            DO il=1,nfgrd
              l=lfltgrd(il)             
              dmm=dmm+Rhoil1(l)*Voil0(l)*Nprt_flt
            ENDDO
            alfa = Moil1/dmm

            IF (alfa>1.0_r8+1.e-10_r8) THEN
              print*,'eulr2lgr: *** ERR alfa>1',alfa,'nfgrd=',nfgrd,    &
     &               'Moil1=',Moil1,'Moil0=',Moil0,' Moil1=',Moil1,     &
     &               'NprtFlt=',Nprt_flt,'GrdVol=',GrdVol(i,j,k)
              DO il=1,nfgrd
                l=lfltgrd(il)
                print*,'eulr2lgr: ERR alfa>1 l=',l,'Rhoil1=',Rhoil1(l), &
     &           'Voil0=',Voil0(l), 'Doil0',Doil0(l),                   &
     &                 'Rhoil0=',Rhoil0(l)
!                DO ic=1,Nocmp
!                  Moil0=Moil0+Coil0(i,j,k,ic)*GrdVol(i,j,k)
!                  Moil1=Moil1+Coil1(i,j,k,ic)*GrdVol(i,j,k)
!                ENDDO
                DO ic=1,3
                  print*,'eulr2lgr: ERR Coil0(ijk,:)=',Coil0(i,j,k,ic)
                  print*,'eulr2lgr: ERR Coil1(ijk,:)=',Coil1(i,j,k,ic)
                ENDDO
                print*,'   ------  '
              ENDDO
              print*,'=========================='
            ENDIF
            
            DO il=1,nfgrd
              l=lfltgrd(il)
              vo0=Voil0(l)

              IF (abs(alfa-1.0_r8).gt.1.0e-20_r8) THEN
                vo1=alfa*vo0
                Doil1(l)=2.0_r8*(3.0_r8/4.0_r8*vo1/pi)**(1.0_r8/3.0_r8)
              ELSE
                Doil1(l)=Doil0(l)
              ENDIF
                 
# ifdef OIL_MAP_DEBUG
              IF (lchck .and. il.eq.ifchck) THEN
                print*,'eulr2lgr: l=',l,',alfa=',alfa,                  &
     &                 'Doil0=',Doil0(l),' Doil1=',Doil1(l)     
                lchck=.FALSE.
              ENDIF
# endif
            ENDDO
!
! Update
!
            DO il=1,nfgrd
              l=lfltgrd(il)
!              tsrfo=track(isrfo,ntt,l)
!              IF (tsrfo.lt.1.e-3) THEN
                track(iroil,ntt,l)=Rhoil1(l)
                track(isizo,ntt,l)=Doil1(l)
! 
! Update initial oil droplet size for evaporation module
                szoil0(l)=Doil1(l)
!              ENDIF
            ENDDO

          ENDDO
        ENDDO
      ENDDO
!
# ifdef OIL_MAP_DEBUG
      IF (maxval(nrcoil).gt.0) THEN
        print*,'eulr2lgr: Rank',MyRank,' # Cnew/Cold < 1 =',nrcoil
        print*,'eulr2lgr: Rank',MyRank,' min Cnew/Cold   =',mnrcoil
        print*,'eulr2lgr: Rank',MyRank,' max Cnew/Cold   =',mxrcoil
      ENDIF
! This shouldn't happen:
      IF (maxval(brcoil).gt.0) THEN
        print*,'eulr2lgr: **ERR MyRank',MyRank, &
     &         ' # Cnew/Cold > 1 =',brcoil
      ENDIF
# endif

      DO ic=1,Nocmp
# ifdef OIL_DEBUG
        IF (wfr_max(ic).gt.0.0_r8)                                      &
     &    print*,'map_Euler2Lgr: compound=',ic,'New min/max(wfr)=',     &
     &          wfr_min(ic),wfr_max(ic),' MyRank=',MyRank
# endif
! zero wfr(1) and wfr(2) possible if bio /sed are coupled
! possibly evaportation could lead to very low saturates/aromatics
        IF (wfr_min(ic).lt.1.0e-20) THEN
          print*,'map_Euler2Lgr: CAUTION wfr1(ic) ZERO: ic=',ic,        &
     &            ' wfr=',wfr_min(ic)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE euler2lagr_tile
#endif
      END MODULE euler2lagr_mod

