#include "cppdefs.h"

      MODULE lagr2euler_mod
#if defined FLOAT_OIL && defined OIL_EULR
!==============================================  Dmitry S Dukhovskoy ===
!                                                                      !
!=======================================================================
!                                                                      !
!  This routine maps oil variables from Lagrangain -> Eulerian coord   !
!                                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: lagr2euler

      CONTAINS

!
!***********************************************************************
      SUBROUTINE lagr2euler(ng, tile, Lstr, Lend)
!***********************************************************************
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
      integer, intent(in) :: ng, tile, Lstr, Lend
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL lagr2euler_tile(ng, tile, Lstr, Lend,                        &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     OIL3D(ng) % Coil,                            &
     &                     OIL3D(ng) % Doil,                            &
     &                     OIL3D(ng) % NFLT3D,                          &
     &                     DRIFTER(ng) % track,                         &
     &                     DRIFTER(ng) % bounded,                       &
     &                     GRID(ng) % omn,                              &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              & 
     &                     GRID(ng) % z_w)

      OIL3D(ng) % Coil0 = OIL3D(ng) % Coil

      RETURN
      END SUBROUTINE lagr2euler

!
!***********************************************************************
      SUBROUTINE lagr2euler_tile(ng, tile, Lstr, Lend,                  &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Coil, Doil, NFLT3D, track,             &
     &                           bounded, omn, Hz, z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_stepping

      implicit none

!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile, Lstr, Lend
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: track(:,0:,:)    
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)

      logical, intent(in) :: bounded(:)

      integer, intent(inout) :: NFLT3D(LBi:,LBj:,:)
      real(r8), intent(inout) :: Coil(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Doil(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng)) 
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))

      logical, intent(in) :: bounded(Nfloats(ng))
      
      integer, intent(inout) :: NFLT3D(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Coil(LBi:UBi,LBj:UBj,N(ng),Nocmp)
      real(r8), intent(inout) :: Doil(LBi:UBi,LBj:UBj,N(ng))
#endif
!
!  Local variables
!
      integer :: i, j, k, ic, l, nfgrd
      integer :: mxNflt, sumFlt  ! for debugging
!
      real(r8) :: Zstr, Zend, Xstr, Xend, Ystr, Yend
      real(r8) :: xflt, yflt, zflt, rhoo
      real(r8) :: frsat, frarom

      real(r8) :: GrdVol(IminS:ImaxS,JminS:JmaxS,N(ng))
      real(r8) :: MoilGrd(IminS:ImaxS,JminS:JmaxS,N(ng),Nocmp)
      real(r8) :: wfroil(Nocmp)
      real(r8) :: mxDoil, mxDoil_lgr ! debugging max oil size in Euler. & Lagr grid

      logical, dimension(Lstr:Lend) :: my_thread
      logical :: lchck

#include "set_bounds.h"
!
!      lchck = .TRUE.

      NFLT3D=0
      MoilGrd=0.0_r8
      GrdVol=0.0_r8
      mxDoil=0.0_r8
      mxDoil_lgr=0.0_r8
      Coil=0.0_r8
      Doil=0.0_r8

!# ifdef OIL_MAP_DEBUG      
!      IF (MyRank.eq.MyMaster) THEN
!        print*,'My Rank =', MyRank, 'tile=',tile
!        print*,'map_Lagr2Euler.F: LBi=',LBi,' UBi=',UBi,                &
!     &         ' LBj=',LBj,' UBj=',UBj
!        print*,'  IminS=',IminS,' ImaxS=',ImaxS,' JminS=',JminS,        &
!     &         '  JmaxS=',JmaxS,'N=',N(ng)
!        print*,' Jstr=',Jstr,' Jend=',Jend,' Istr=',Istr,' Iend=',Iend
!        print*,'max NFLT', maxval(NFLT3D)
!        print*,'Size I NFLT3D=',size(NFLT3D,1),' should be',UBi-LBi+1
!        print*,'Size J NFLT3D=',size(NFLT3D,2),' should be',UBj-LBj+1
!        print*,'Size K NFLT3D=',size(NFLT3D,3),' should be',N(ng)
!# endif
!
! Count # of floats in every grid cell
! Then calculate oil mass in every gird cell
! and finally oil conc = M oil / grid cell
! Do this for all components 
      DO k=1,N(ng)
        GrdVol(IminS:ImaxS,JminS:JmaxS,k)=                              &
     &    omn(IminS:ImaxS,JminS:JmaxS)*Hz(IminS:ImaxS,JminS:JmaxS,k)
      ENDDO

! Floats only that belong to the thread:
      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8

      mxNflt=0
      DO l=Lstr,Lend
        my_thread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nfp1(ng),l)).and.                      &
     &      (track(ixgrd,nfp1(ng),l).lt.Xend).and.                      &
     &      (Ystr.le.track(iygrd,nfp1(ng),l)).and.                      &
     &      (track(iygrd,nfp1(ng),l).lt.Yend)) THEN
          my_thread(l)=.TRUE.
        ELSE
          CYCLE
        ENDIF

        IF (my_thread(l).and.bounded(l)) THEN
          xflt=track(ixgrd,nfp1(ng),l)
          yflt=track(iygrd,nfp1(ng),l)
          zflt=track(izgrd,nfp1(ng),l)
          mxNflt=mxNflt+1
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

          NFLT3D(i,j,k)=NFLT3D(i,j,k)+1
          Doil(i,j,k)=Doil(i,j,k)+track(isizo,nfp1(ng),l)
          rhoo=track(iroil,nfp1(ng),l)
          mxDoil_lgr=max(mxDoil_lgr,track(isizo,nfp1(ng),l))

!     
!     This part of the code will need to be modified 
!     if varying Nocmp = # of oil components (fractions) is needed
!     Currently it is hard-coded for Nocmp=3 
!     SAR+A components in agreement
!     with oil_plume.h and Modules/mod_floats.F
!     
!     Weight fraction of oil components SAR+A:
          wfroil(1)=track(isats,nfp1(ng),l)
          wfroil(2)=track(iarom,nfp1(ng),l)
          wfroil(3)=1.0_r8-(wfroil(1)+wfroil(2))

!     Oil mass in grd cell by components
          DO ic=1,Nocmp
             MoilGrd(i,j,k,ic)=MoilGrd(i,j,k,ic)+                       &
     &                         VoilFlt*rhoo*wfroil(ic)
          ENDDO

!#if defined OIL_MAP_DEBUG
!          IF (i.eq.141 .and. j.eq.79 .and. k.eq.4) THEN
!            print*,'### i,j,k: ',i,j,k, 'NFLT=',NFLT3D(i,j,k),          &
!     &           ' Doil=',Doil(i,j,k),'Lagr Doil=',track(isizo,nfp1(ng),l)
!            print*,'### Moil1=',MoilGrd(i,j,k,1),' rhoo=',rhoo,         &
!     &          'wfroil=',wfroil(1)
!            print*,'### Moil2=',MoilGrd(i,j,k,2),' rhoo=',rhoo,         &
!     &          'wfroil=',wfroil(2)
!            print*,'### Moil3=',MoilGrd(i,j,k,3),' rhoo=',rhoo,         &
!     &          'wfroil=',wfroil(3)
!          ENDIF
!#endif

        ENDIF
      ENDDO

!# ifdef OIL_MAP_DEBUG      
!      IF (MyRank.eq.MyMaster) THEN
!      IF (lchck .and. maxval(NFLT3D)>0) THEN
!        print*,'My Rank =', MyRank, 'maxNFLT=',maxval(NFLT3D)
!      ENDIF
!# endif
!

!      mxNflt=0
      sumFlt=0
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            nfgrd=NFLT3D(i,j,k)
!            mxNflt=max(mxNflt,nfgrd)
            sumFlt=sumFlt+nfgrd
            IF (nfgrd.eq.0) THEN
              IF (Doil(i,j,k).ge.1.0e-20_r8) THEN
                print*,'map_Lagr2Eulr: ** ERR: nflt=',nfgrd,            &
     &                 ' Doil=',Doil(i,j,k),' i,j,k=',i,j,k
              ENDIF
              CYCLE
            ENDIF
            Doil(i,j,k)=Doil(i,j,k)/REAL(nfgrd) ! mean oil size, gr cell
            mxDoil=max(mxDoil,Doil(i,j,k))
!            
            DO ic=1,Nocmp
              Coil(i,j,k,ic)=MoilGrd(i,j,k,ic)/GrdVol(i,j,k)
            ENDDO
!
            IF (NFLT3D(i,j,k).ge.1000 .or. mxDoil.gt.mxDoil_lgr) THEN ! mxDoil (cell mean) should be <= mxDoil_lgr 
              print*,'*** ERR  My rank',MyRank,                         &
     &               'nfgrd=',nfgrd
              print*,'i=',i,'j=',j,'k=',k,'Nflts=',nfgrd
              print*,'Max Doil Eulr=',mxDoil,                           &
     &               ' Max Doil Lagr=',mxDoil_lgr
              print*,'   *********  '
            ENDIF

          ENDDO
        ENDDO
      ENDDO

#ifdef OIL_MAP_DEBUG
      IF (sumFlt>0 .and. mxDoil.gt.1.0e-1) THEN  
        print*,'lagr2euler_tile end, CHECK: Eulr MaxDoil=',mxDoil,      &
     &         'Lagr MaxDoil=',mxDoil_lgr, 'Eulr # Floats=',sumFlt,     &
     &         'Lagr # Floats=',mxNflt
      ENDIF
#endif
      RETURN
      END SUBROUTINE lagr2euler_tile

#endif
      END MODULE lagr2euler_mod 
