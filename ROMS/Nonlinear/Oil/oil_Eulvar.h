! Subroutines for allocating and initializing
! 3D arrays for keeping oil variables
! in Eularian coordinates for
! providing input fields to the Biological and sediments modules
!
!=======================================================================
!                                                                      !
!=======================================================================
!                                                                      !
! Allocate and initialize arrays                                       !
!                                                                      !
!  COil - Oil Content/concentr in a grid cell, kg/m3 3D array          !
!         Coil(:,:,:,1:Nocmp) => Oil content/conc by oil components    !
!         Coil = Mass oil/Vgrid by oil components                      !
!                                                                      !
!  COil0 - Oil concentration field before bio/sedimentation            !
!          Needed for mapping back to Lagrangian frame                 !
!                                                                      !
!  Doil - Average Oil droplet size in 3D grid cells                    !
!                                                                      !
!  NFLT3D - Number of oil Lagr. floats in 3D grid cells                !
!                                                                      !
!                                                                      !
!=======================================================================

      USE mod_kinds
      USE mod_floats
!
      implicit none
!
      TYPE T_OIL3D
!
! Oil Conc in grid cells, 3D array by components
! Mean droplet size
!
      real(r8), pointer :: Coil(:,:,:,:)
      real(r8), pointer :: Coil0(:,:,:,:)
      real(r8), pointer :: Doil(:,:,:)
      integer, pointer :: NFLT3D(:,:,:)
!
      END TYPE T_OIL3D

      TYPE (T_OIL3D), allocatable :: OIL3D(:)

      CONTAINS
!
      SUBROUTINE allocate_oil_Eulvar (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates 3D arrays for Oil variables                  !
!  mapped onto  Eularian coordinates from Lagrangian                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_parallel  !! 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OIL3D(Ngrids) )
!
      allocate ( OIL3D(ng) % Coil(LBi:UBi,LBj:UBj,N(ng),Nocmp) )
      allocate ( OIL3D(ng) % Coil0(LBi:UBi,LBj:UBj,N(ng),Nocmp) )
      allocate ( OIL3D(ng) % Doil(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OIL3D(ng) % NFLT3D(LBi:UBi,LBj:UBj,N(ng)) ) 
!
! Debug
      IF (MyRank.eq.MyMaster) THEN
        print*,'Allocating: My Rank =', MyRank
        print*,'  LBi=',LBi,' UBi=',UBi,' LBj=',LBj,' UBj=',UBj
        print*,'Size I NFLT3D=',size( OIL3d(ng) % NFLT3D,1 )
        print*,'Size J NFLT3D=',size( OIL3d(ng) % NFLT3D,2 )
        print*,'Size K NFLT3D=',size( OIL3d(ng) % NFLT3D,3 )
      ENDIF
! end debug
!
      RETURN
      END SUBROUTINE allocate_oil_Eulvar

      SUBROUTINE initialize_oil_Eulvar (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the  shared     !
!  arrays  across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k

      real(r8), parameter :: IniVal = 0.0_r8
#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize oil structure variables, only nonlinear model state
!-----------------------------------------------------------------------
!
!
!      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO itrc=1,Nocmp
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                OIL3D(ng) % Coil(i,j,k,itrc) = IniVal
                OIL3D(ng) % Coil0(i,j,k,itrc) = IniVal
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              OIL3D(ng) % Doil(i,j,k) = IniVal
              OIL3D(ng) % NFLT3D(i,j,k) = 0
            ENDDO
          ENDDO
        ENDDO
!     ENDIF


      RETURN
      END SUBROUTINE initialize_oil_Eulvar



