/*
 * **        oil_Eulvar_def.h
 * ************************************************* Dmitry Dukhovskoy  ***
 * **                                                                    **
 * ************************************************************************
 * **                                                                    **
 * **  Defines oil model input parameters in output NetCDF files.        **
 * **  Saved with 3D Eulerian oil fields 
 * **  It is included in routine "def_info.F".                           **
 * **                                                                    **
 * ************************************************************************
 * */

!
!  Define parameters of Eulerina oil model fields
!  
! Not finished: May include information about
! initial Doil, densities by components, other parameters 
! 

      Vinfo( 1)='Doil0'
      Vinfo( 2)='Mean size of oil particles, Gamma distr'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/doildim/), Aval, Vinfo, ncname,               &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      Vinfo( 1)='Coil'
      Vinfo( 2)='Oil concentration by oil components'
      Vinfo( 3)='kg m-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/coildim/), Aval, Vinfo, ncname,               &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN



