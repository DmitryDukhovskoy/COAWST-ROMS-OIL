!
!svn $Id: sedflocs_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sediment opa Model Kernel Variables:                               !
!                                                                      !
#if defined SEDIMENT && defined SED_OPA
!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,itauc) => mean critical erosion stress  !
#endif
!                                                                      !
!=======================================================================

      USE mod_kinds

      implicit none


      TYPE T_SEDOPA

#if defined SEDIMENT && defined SED_OPA
      real(r8), pointer :: f_diam(:)
      real(r8), pointer :: f_vol(:)
      real(r8), pointer :: f_rho(:)
      real(r8), pointer :: f_mass(:)

      real(r8), pointer :: oil_diam(:)
      real(r8), pointer :: oil_vol(:)
      real(r8), pointer :: oil_rho(:)
      real(r8), pointer :: oil_mass(:)

      real(r8), pointer :: opa_diam(:)
      real(r8), pointer :: opa_vol(:)
      real(r8), pointer :: opa_rho(:)
      real(r8), pointer :: opa_mass(:)

      real(r8),pointer :: opa_g1(:,:,:)
      real(r8),pointer :: oil_l1(:,:)
      real(r8),pointer :: sed_l1(:,:)
      real(r8),pointer :: opa_coll_prob1(:,:)
      real(r8),pointer :: opa_oil1(:,:,:)

      real(r8),pointer :: opa_g2(:,:,:)
      real(r8),pointer :: opa_l2(:,:)
      real(r8),pointer :: oil_l2(:,:)
      real(r8),pointer :: opa_coll_prob2(:,:)
      real(r8),pointer :: opa_oil2(:,:,:)

      real(r8),pointer :: opa_g3(:,:,:)
      real(r8),pointer :: opa_l3(:,:)
      real(r8),pointer :: sed_l3(:,:)
      real(r8),pointer :: opa_coll_prob3(:,:)

      real(r8),pointer :: sed_g4(:,:,:)
      real(r8),pointer :: sed_l4(:,:)
      real(r8),pointer :: opa_coll_prob4(:,:)
#endif
      END TYPE T_SEDOPA

      TYPE (T_SEDOPA), allocatable :: SEDOPA(:)

      CONTAINS

      SUBROUTINE allocate_sedopa (ng, LBi, UBi, LBj, UBj)

!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj

      integer, parameter :: NSED=3
      integer, parameter :: NOIL=1
      integer, parameter :: NOPA=4
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SEDOPA(Ngrids) )
!
!  Nonlinear model state.

#if defined SEDIMENT && defined SED_OPA
      allocate (SEDOPA(ng) % f_diam(NSED))
      allocate (SEDOPA(ng) % f_vol(NSED))
      allocate (SEDOPA(ng) % f_rho(NSED))
      allocate (SEDOPA(ng) % f_mass(0:NSED+1))

      allocate (SEDOPA(ng) % oil_diam(NOIL))
      allocate (SEDOPA(ng) % oil_vol(NOIL))
      allocate (SEDOPA(ng) % oil_rho(NOIL))
      allocate (SEDOPA(ng) % oil_mass(NOIL))

      allocate (SEDOPA(ng) % opa_diam(NOPA))
      allocate (SEDOPA(ng) % opa_vol(NOPA))
      allocate (SEDOPA(ng) % opa_rho(NOPA))
      allocate (SEDOPA(ng) % opa_mass(0:NOPA+1))

      allocate (SEDOPA(ng) %  opa_g1(NOIL,NSED,NOPA))
      allocate (SEDOPA(ng) %  oil_l1(NOIL,NSED))
      allocate (SEDOPA(ng) %  sed_l1(NOIL,NSED))
      allocate (SEDOPA(ng) %  opa_coll_prob1(NOIL,NSED))
      allocate (SEDOPA(ng) %  opa_oil1(NOIL,NSED,NOPA))

      allocate (SEDOPA(ng) %  opa_g2(NOIL,NOPA,NOPA))
      allocate (SEDOPA(ng) %  opa_l2(NOIL,NOPA))
      allocate (SEDOPA(ng) %  oil_l2(NOIL,NOPA))
      allocate (SEDOPA(ng) %  opa_coll_prob2(NOIL,NOPA))
      allocate (SEDOPA(ng) %  opa_oil2(NOIL,NOPA,NOPA))

      allocate (SEDOPA(ng) %  opa_g3(NSED,NOPA,NOPA))
      allocate (SEDOPA(ng) %  opa_l3(NSED,NOPA))
      allocate (SEDOPA(ng) %  sed_l3(NSED,NOPA))
      allocate (SEDOPA(ng) %  opa_coll_prob3(NSED,NOPA))

      allocate (SEDOPA(ng) %  sed_g4(NSED,NSED,NSED))
      allocate (SEDOPA(ng) %  sed_l4(NSED,NSED))
      allocate (SEDOPA(ng) %  opa_coll_prob4(NSED,NSED))

#endif

      RETURN
      END SUBROUTINE allocate_sedopa

      
      SUBROUTINE initialize_sedopa (ng, tile, model)

!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer, parameter :: NSED=3
      integer, parameter :: NOIL=1
      integer, parameter :: NOPA=4
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
!  Initialize sediment structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        CALL initialize_sedopa_param (ng, tile,                       &
     &                      SEDOPA(ng) %  f_mass,                       & 
     &                      SEDOPA(ng) %  f_diam,                       &
     &                      SEDOPA(ng) %  oil_mass,                     &
     &                      SEDOPA(ng) %  oil_diam,                     &
     &                      SEDOPA(ng) %  opa_mass,                     &
     &                      SEDOPA(ng) %  opa_diam,                     &
     &                      SEDOPA(ng) %  opa_g1,                       &
     &                      SEDOPA(ng) %  oil_l1,                       &
     &                      SEDOPA(ng) %  sed_l1,                       &
     &                      SEDOPA(ng) %  opa_coll_prob1,               &
     &                      SEDOPA(ng) %  opa_oil1,                     &
     &                      SEDOPA(ng) %  opa_g2,                       &
     &                      SEDOPA(ng) %  opa_l2,                       &
     &                      SEDOPA(ng) %  oil_l2,                       &
     &                      SEDOPA(ng) %  opa_coll_prob2,               &
     &                      SEDOPA(ng) %  opa_oil2,                     &
     &                      SEDOPA(ng) %  opa_g3,                       &
     &                      SEDOPA(ng) %  opa_l3,                       &
     &                      SEDOPA(ng) %  sed_l3,                       &
     &                      SEDOPA(ng) %  opa_coll_prob3,               &
     &                      SEDOPA(ng) %  sed_g4,                       &
     &                      SEDOPA(ng) %  sed_l4,                       &
     &                      SEDOPA(ng) %  opa_coll_prob4)
      END IF
!
      RETURN
      END SUBROUTINE initialize_sedopa
!
!***********************************************************************
      SUBROUTINE initialize_sedopa_param (ng, tile,                     &
     &           f_mass,f_diam,oil_mass,oil_diam,opa_mass,opa_diam,     &
     &           opa_g1,oil_l1,sed_l1,opa_coll_prob1,opa_oil1,          &
     &           opa_g2,opa_l2,oil_l2,opa_coll_prob2,opa_oil2,          &
     &           opa_g3,opa_l3,sed_l3,opa_coll_prob3,                   &
     &           sed_g4,sed_l4,opa_coll_prob4)
!***********************************************************************

      USE mod_param
      USE mod_scalars
      USE mod_sediment

      implicit none

!  Imported variable declarations.
!

      integer, intent(in) :: ng, tile
      integer, parameter  :: NSED=3
      integer, parameter  :: NOIL=1
      integer, parameter  :: NOPA=4

      real(r8), intent(inout)  :: f_mass(0:NSED+1)
      real(r8), intent(inout)  :: f_diam(NSED)
      real(r8) :: f_vol(NSED)
      real(r8) :: f_rho(NSED)
      real(r8), intent(inout)  :: oil_mass(NOIL)
      real(r8), intent(inout)  :: oil_diam(NOIL)
      real(r8) :: oil_vol(NOIL)
      real(r8) :: oil_rho(NOIL)
      real(r8), intent(inout)  :: opa_mass(0:NOPA+1)
      real(r8), intent(inout)  :: opa_diam(NOPA)
      real(r8) :: opa_vol(NOPA)
      real(r8) :: opa_rho(NOPA)
      
      real(r8), intent(inout)  :: opa_g1(NOIL,NSED,NOPA)
      real(r8), intent(inout)  :: oil_l1(NOIL,NSED)
      real(r8), intent(inout)  :: sed_l1(NOIL,NSED)
      real(r8), intent(inout)  :: opa_coll_prob1(NOIL,NSED)
      real(r8), intent(inout)  :: opa_oil1(NOIL,NSED,NOPA)

      real(r8), intent(inout)  :: opa_g2(NOIL,NOPA,NOPA)
      real(r8), intent(inout)  :: opa_l2(NOIL,NOPA)
      real(r8), intent(inout)  :: oil_l2(NOIL,NOPA)
      real(r8), intent(inout)  :: opa_coll_prob2(NOIL,NOPA)
      real(r8), intent(inout)  :: opa_oil2(NOIL,NOPA,NOPA)

      real(r8), intent(inout)  :: opa_g3(NSED,NOPA,NOPA)
      real(r8), intent(inout)  :: opa_l3(NSED,NOPA)
      real(r8), intent(inout)  :: sed_l3(NSED,NOPA)
      real(r8), intent(inout)  :: opa_coll_prob3(NSED,NOPA)

      real(r8), intent(inout)  :: sed_g4(NSED,NSED,NSED)
      real(r8), intent(inout)  :: sed_l4(NSED,NSED)
      real(r8), intent(inout)  :: opa_coll_prob4(NSED,NSED)

      real(r8), parameter :: IniVal = 0.0_r8 
!------------------------------------------------------------------------
! Initialize OPA structure variables.
!------------------------------------------------------------------------
!
!
             SEDOPA(ng) %  f_mass(:) = IniVal        
             SEDOPA(ng) %  f_diam(:) = IniVal        
             SEDOPA(ng) %  oil_mass(:) = IniVal
             SEDOPA(ng) %  oil_diam(:) = IniVal
             SEDOPA(ng) %  opa_mass(:) = IniVal
             SEDOPA(ng) %  opa_diam(:) = IniVal
             SEDOPA(ng) %  opa_g1(:,:,:) = IniVal
             SEDOPA(ng) %  oil_l1(:,:)   = IniVal
             SEDOPA(ng) %  sed_l1(:,:)   = IniVal
             SEDOPA(ng) %  opa_coll_prob1(:,:) = IniVal
             SEDOPA(ng) %  opa_oil1(:,:,:)  = IniVal
             SEDOPA(ng) %  opa_g2(:,:,:) = IniVal
             SEDOPA(ng) %  opa_l2(:,:) = IniVal
             SEDOPA(ng) %  oil_l2(:,:) = IniVal
             SEDOPA(ng) %  opa_coll_prob2(:,:) = IniVal
             SEDOPA(ng) %  opa_oil2(:,:,:) = IniVal
             SEDOPA(ng) %  opa_g3(:,:,:)   = IniVal
             SEDOPA(ng) %  opa_l3(:,:) = IniVal
             SEDOPA(ng) %  sed_l3(:,:) = IniVal
             SEDOPA(ng) %  opa_coll_prob3(:,:) = IniVal
             SEDOPA(ng) %  sed_g4(:,:,:) = IniVal
             SEDOPA(ng) %  sed_l4(:,:) = IniVal
             SEDOPA(ng) %  opa_coll_prob4(:,:) = IniVal

     RETURN
     END SUBROUTINE initialize_sedopa_param
