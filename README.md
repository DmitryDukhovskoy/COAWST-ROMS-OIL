Coupled ROMS-OIL-SEDOPA-BIOCHEMICAL model 
Can be integrated in online (regular) or offline modes
For offline mode, a climatology file from an online integration should be prepared
where all needed fields have to be saved: u,v,, ... and input fileds (tracers) 
for bio (NO3, NH4)  or sediment

Description of the CSOMIO system is here:
Dukhovskoy et al., 2021. Development of the CSOMIO Coupled Ocean-Oil-Sediment- Biology Model.
Frontiers in Marine Science
https://www.frontiersin.org/articles/10.3389/fmars.2021.629299/full

Description of offline simulation is in this paper:
Kristen M. Thyng et al. 2020; Performance of offline passive tracer advection 
in the Regional Ocean Modeling System (ROMS; v3.6, revision 904);
https://gmd.copernicus.org/articles/14/391/2021/
 
# ROMS-OIL
Oil plume model coupled with ROMS

Oil plume model utilizes Lagrangian float option in ROMS with added vertical velocity determined from the oil droplet sizes, oil component wegith fraction, density. At the surface, oil droplets are subject to weathering (evaporation) and wind drift (optional). Oil droplet size is assigned by generating random number from Gamma distribution. Two alogirhtms calculating vertical velocity of oil droplets are available and defined by users (two-equation or integrated). 

The oil module resides in the directory ROMS/Nonlinear/Oil. However, there are multiple changes throughout the ROMS code to incroporate the oil module, oil array passing, and I/O. Several new script files have been added in the ROMS directories. Thus it is recommended to copy whole ROMS directory into COAWST. 

To compile Oil module, one needs to add directories with the oil module and includes/headers in makefile script:
ifdef USE_ROMS
 includes +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Oil \
		ROMS/Nonlinear/Sediment \
...


ifdef USE_ROMS
 modules  +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Oil \
...

Note the floats include file has been modified to add input information for the oil model:
e.g., more <expt>_floats.in

NFLOATS == 4000
Lfloats == T
Fprint == T
FRREC == 0
Nocmp == 3
OilComp == saturates aromatics res_asph
OilDens == 800.0 850.0 1030.0
DoilMn == 350.0e-6
DoilMin == 0.19e-6
FlxVoilDay == 3000.0
CWndOil == 0.01
DfsOil == 20.0

POS = G, C, T,   N,     Ft0,    Fx0,      Fy0,     Fz0,    Fdt,   Fdx,   Fdy,   Fdz

      1  1  1 2000  0.d0   -88.367d0 28.737d0  -1400.d0  0.03d0 0.d0   0.d0  0.d0
      1  1  1 2000  0.d0   -88.369d0 28.739d0  -1400.d0  0.03d0 0.d0   0.d0  0.d0


#ROMS-OIL-SEDIMENT-BIO ONLINE
Change definitions in the <expt>.h
To turn the oil model on with Lagrangian -> Euler mapping of the oil fields
two-equaltion algorithm for the oil plume is implemented

#define FLOATS\
#define FLOAT_OIL\
#define FLOAT_VWALK\
#undef WOIL_INTEGRATED\
#undef OIL_DEBUG\
#define OIL_EULR\
  
wind drift is automaticly "on" when bulk fluxes are defined unless CWndOil =0 in the <expt>_floats.in

#define  BULK_FLUXES\
  
To turn off the oil model and ROMS as a stand-alone model:

#undef FLOATS\
#undef FLOAT_OIL\
#undef FLOAT_VWALK\
#undef WOIL_INTEGRATED\
#undef OIL_DEBUG\
#undef OIL_EULR\

Do not turn OIL_DEBUG on, this will produce large log file. OIL_EULR can be turned off if mapping is not needed. 


# ROMS-OIL-SEDIMENT-BIO OFFLINE
See Thyng et al. paper on how to prepare offline simulation
https://gmd.copernicus.org/articles/14/391/2021/

Change definitions in the <expt>.h file to allow offline with/without Sediment and Bio
see /Projects/oil_sedbiofR/oil_sedbiofR.h

The code has been run for several simulations but it has not been thoroughly tested
Report any bugs to ddukhovskoy@fsu.edu. 


