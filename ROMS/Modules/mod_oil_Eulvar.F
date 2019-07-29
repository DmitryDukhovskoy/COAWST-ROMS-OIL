#include "cppdefs.h"
      MODULE mod_oil_Eulvar
/*
*************************************************** D Dukhovskoy      **
** This module declares  parameters of oil fields on Eulerian grid    **
** The fields are derived from Lagrangian oil floats                  **
**                                                                    **
** Note that all the *.h files are located in ROMS/Nonlinear/Oil      **
** and included within  <...>  to allow the user to customize any of  **
** them in the project directory  while keeping the distributed code  **
** intact (check the build script for details).                       **
**                                                                    **
**   *  Internal model parameters declarations:                       **
**                                                                    **
**      oil_Eulvar.h                                                  **
**                                                                    **
**   *  Code to define input model parameters in all output           **
**      NetCDF files:                                                 **
**                                                                    **
**      oil_Eulvar_def.h                                              **
**                                                                    **
************************************************************************
*/
#if defined FLOATS && defined FLOAT_OIL
# include <oil_Eulvar.h>
#endif

      END MODULE mod_oil_Eulvar
