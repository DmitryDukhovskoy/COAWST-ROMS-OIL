#include "cppdefs.h"
      MODULE oil_floats_mod
#if defined NONLINEAR && defined FLOATS && defined FLOAT_OIL
!
! D. Dukhovskoy FSU COAPS 
!=======================================================================
!                                                                      !
!  These routine computes plume dynamics for oil Lagrangian floats.    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: oil_floats, oil_param_initialize, oil_density

      CONTAINS

# ifdef FLOAT_OIL
#  include <oil_plume.h>
# endif

#endif

      END MODULE oil_floats_mod
