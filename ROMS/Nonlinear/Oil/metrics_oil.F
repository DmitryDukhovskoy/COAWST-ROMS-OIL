#include "cppdefs.h"
      MODULE metrics_oil_mod
#if defined FLOAT_OIL && defined OIL_EULR
!
!=======================================================================
!                                                                      !
!========================================== Dmitry S Dukhovskoy      ===
!                                                                      !
!     Computes diagnsotics of oil 3D fields and reports                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: metrics_oil

      CONTAINS
!
!***********************************************************************
      SUBROUTINE metrics_oil (ng, tile, model, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_floats
      USE mod_stepping
      USE mod_oil_Eulvar
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: Lstr, Lend
!
!  Local variable declaration
!
#include "tile.h"
!
!      print*,'metrics_oil, tile=',tile

      CALL metrics_oil_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Lstr, Lend, model,                         &
     &                       OIL3D(ng) % Coil,                          &
     &                       OIL3D(ng) % Doil,                          &
     &                       OIL3D(ng) % NFLT3D,                        &
     &                       nf(ng),                                    &
     &                       DRIFTER(ng) % track,                       &
     &                       DRIFTER(ng) % bounded)

      RETURN
      END SUBROUTINE metrics_oil
!
!***********************************************************************
      SUBROUTINE metrics_oil_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Lstr, Lend, model,                   & 
     &                             Coil, Doil, NFLT3D, nf,              &
     &                             track, bounded)                   
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_iounits
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif

      implicit none

!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lstr, Lend, model, nf
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: track(:,0:,:)

      logical, intent(in) :: bounded(:)

      integer, intent(inout) :: NFLT3D(LBi:,LBj:,:)
      real(r8), intent(inout) :: Coil(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Doil(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng))

      logical, intent(in) :: bounded(Nfloats(ng))

      integer, intent(inout) :: NFLT3D(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Coil(LBi:UBi,LBj:UBj,N(ng),Nocmp)
      real(r8), intent(inout) :: Doil(LBi:UBi,LBj:UBj,N(ng))
#endif
! 
! Local Variables
!
      integer :: nfltLagr, nfltEulr, nff
      integer :: i, j, k, ic, l, Nflt2
#ifdef DISTRIBUTE
      character (len=3), dimension(6) :: op_handle
#endif
      real(r8), dimension(6) :: buffer
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8) :: minC, maxC, avC, minD, maxD, avD
      real(r8) :: dmm, cmm
      logical, dimension(Lstr:Lend) :: my_thread

#include "set_bounds.h"

      buffer    = 0.0
      buffer(3) = 1.0e6
      buffer(5) = 1.0e6 

!  Calculate total number of active floats in Lagr. space
! There are 2 ways to count floats
! Keep them both as an example
! (1) sum across the tiles - for DISTRIBUTED memory configuration
      Xstr=REAL(Istr,r8)-0.5_r8
      Xend=REAL(Iend,r8)+0.5_r8
      Ystr=REAL(Jstr,r8)-0.5_r8
      Yend=REAL(Jend,r8)+0.5_r8
      DO l=Lstr,Lend
        my_thread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          my_thread(l)=.TRUE.
        ENDIF


        IF (my_thread(l) .and. bounded(l)) THEN
          buffer(1)=buffer(1)+1.0
        ENDIF
      ENDDO
!
! (2) since bounded() array is transmitted across the tiles
! simply sum over 1 tile
      Nflt2=0
#ifdef DISTRIBUTE
      IF (Master) THEN
#endif
        DO l=Lstr,Lend
          IF (bounded(l)) Nflt2=Nflt2+1
        ENDDO
#ifdef DISTRIBUTE
      ENDIF
#endif

!  Calculate total number of active floats in Eulerian space (3D array)
!  min, max, mean Coil and Doil
      minD=100.0_r8
      minC=1.e6_r8
      maxD=0.0_r8
      maxC=0.0_r8
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            nff=NFLT3D(i,j,k)
# if defined OIL_DEBUG
            IF (nff.gt.0 .and. Doil(i,j,k).le.1.0e-20_r8) THEN
                print*,'metrics_oil: *** ERR minD=',dmm,'i,j,k=',i,j,k  
            ENDIF
# endif
            IF (nff.eq.0) CYCLE

            buffer(2)=buffer(2)+REAL(nff)
            dmm=Doil(i,j,k)
            IF (nff>0) THEN
              minD=min(minD,dmm)
              maxD=max(maxD,dmm)
              IF (dmm.le.1.0e-20_r8 .or. minD.le.1.0e-20_r8) THEN
                print*,'metrics_oil: *** ERR dmm=',dmm,'minD=',minD,    &
     &                 'i,j,k=',i,j,k  
              ENDIF
            ENDIF
            buffer(3)=minD
            buffer(4)=maxD

            cmm=Coil(i,j,k,1) ! saturates
            IF (nff>0) THEN
              minC=min(minC,cmm)
              maxC=max(maxC,cmm)
            ENDIF
            buffer(5)=minC
            buffer(6)=maxC
          ENDDO
        ENDDO
      ENDDO


#ifdef DISTRIBUTE
      op_handle(1)='SUM'
      op_handle(2)='SUM'
      op_handle(3)='MIN'
      op_handle(4)='MAX'
      op_handle(5)='MIN'
      op_handle(6)='MAX'

      CALL mp_reduce (ng, model, 6, buffer, op_handle) 

      nfltLagr=buffer(1)
      nfltEulr=buffer(2)
      minD=buffer(3)
      maxD=buffer(4)
      minC=buffer(5)
      maxC=buffer(6)
#else
      nfltLagr=buffer(1)
      nfltEulr=buffer(2)
#endif      

#ifdef OIL_DEBUG
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
        write(stdout,10) nfltLagr, Nflt2, nfltEulr, minD, maxD,         &
     &                   minC, maxC
# ifdef DISTRIBUTE
      ENDIF
# endif
#endif

 10   FORMAT(/,'    Lagrangian -> Eulerian Mapping:             ',/,    &
     &         '================================',/,                    &
     &         'Released floats: (1) Lagr=',i6,' (2) Lagr=',i6,         &
     &         '  N floats Mapped to Eulerian=',i6,/,                   &
     &         'min/max Coil(kg/m3)=',2e10.3,                           &
     &         ' min/max Doil(m)=',2e10.3,/)


      RETURN
      END SUBROUTINE metrics_oil_tile
#endif
      END MODULE metrics_oil_mod

