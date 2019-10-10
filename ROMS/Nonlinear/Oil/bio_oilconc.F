#include "cppdefs.h"

      MODULE bio_oil_mod
#if defined FLOAT_OIL && defined OIL_EULR && defined OIL_BIO
!==============================================  Dmitry S Dukhovskoy ===
!                                                                      !
!=======================================================================
!                                                                      !
!  This routine mimics biodegradataion process by simply reducing      !
!  the concentration of individual oil components (saturates,          !
!  aromatics, asphaltines+raisins)                                     !
!  i.e. acts as a sink term                                            !
!  This should be replaced by a bio-model that will simulate oil       !
!  consumption by microbes. This will result in reducing (at different !
!  rate) the concentration of light - heavy components of oil          !
!  leading to changes in oil particles' densities and size             !
!  which is implemented in map_Eulr2Lagr.F                             !
!                                                                      !
!  Bio-degradation is turned on only until oil particles hits          !
!  surface, thus limited to layers below 5 surface layers              !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: bio_degradation

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bio_degradation(ng, tile)
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
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      
      CALL bio_degradation_tile(ng, LBi, UBi, LBj, UBj,                 &
     &                           OIL3D(ng) % Coil)

      RETURN
      END SUBROUTINE bio_degradation

!***********************************************************************
      SUBROUTINE bio_degradation_tile(ng, LBi, UBi, LBj, UBj,           &
     &                                Coil1)
!***********************************************************************
      USE mod_param
      USE mod_floats
!      USE mod_scalars
!      USE mod_grid
!      USE mod_parallel

      implicit none

!
!  Imported variable declarations
!
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Coil1(LBi:,LBj:,:,:)
#else
      real(r8), intent(inout) :: Coil1(LBi:UBi,LBj:UBj,N(ng),Nocmp)
#endif
       
!
!  Local variables
!
      integer :: ic
      real(r8) :: rr(Nocmp)    ! reduction rate of oil conc by oil components
!      real(r8) :: zoil
!#include "set_bounds.h"
!
! In subsurface layers only
! Oil componenet reduction is 
! Ci(t)=(1-rr)Ci(t-1)
!
      rr(1)=1.e-4_r8
      rr(2)=1.e-6_r8
      rr(3)=1.e-12_r8

      DO ic=1,Nocmp
        Coil1(:,:,1:N(ng)-5,ic)=Coil1(:,:,1:N(ng)-5,ic)*(1.0_r8-rr(ic))
      ENDDO

      RETURN
      END SUBROUTINE bio_degradation_tile
# endif
      END MODULE bio_oil_mod

