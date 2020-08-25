/*
** svn $Id: ias.h 830 2017-01-24 21:21:11Z arango $
*******************************************************************************
** Copyright (c) 2002-2018 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for NGOM 1/25 DEG resolution.
** CSOMIO coupled OIL-BIO GENOME-SEDOPA
** To run a fully coupled system:
** define : OIL_BIO OIL_SEDIMENT FLOATS FLOAT_OIL
**          OIL_EULR should be automatically defined 
**          with any of SEDIMENT or BIO_CSOMIO turned on
**          
** Application flag:   oil_sedbiof
** Input script:       ocean_oil_sedbiof.in
*/

#undef WRF_MODEL
#undef SWAN_MODEL
#define ROMS_MODEL

/* Oil */
#define OUT_DOUBLE

/*
*******-----------------------------------------------------------------------------
*******  Offline settings
*******-----------------------------------------------------------------------------
******/
#define OFFLINE
#ifdef OFFLINE
# define OCLIMATOLOGY
# undef ANA_AKTCLIMA  /* need this for offline */


#define MIXCLIMATOLOGY
#define AKXCLIMATOLOGY
#define AKSCLIMATOLOGY
/* for offline, turn off forcings (bulk forcing undefined elsewhere in file)
 * All forcing is coming in through climatology from online case. */
/*#define ANA_INITIAL*/

/* #define ANA_VMIX */

/*
******-----------------------------------------------------------------------------
******  Adding offline passive tracers
*****-----------------------------------------------------------------------------
****/
#undef OFFLINE_TPASSIVE  
# ifdef OFFLINE_TPASSIVE
#  define T_PASSIVE
#  define ANA_BPFLUX
#  define ANA_SPFLUX
# endif
#endif


/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state settings.
**-----------------------------------------------------------------------------
*/
#define  AVERAGES               /*Write out time-averaged data*/
#undef  AVERAGES_FLUXES
#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV               /*mixing of momentum along constant s surfaces*/
#define SPLINES_VDIFF
#define SPLINES_VVISC
#undef TS_HSIMT
#undef TS_U3HADVECTION
#undef TS_C4VADVECTION
#define TS_MPDATA
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define MASKING
#undef  SRELAXATION            /*salinity relaxation as freshwater flux*/

#undef LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define  GLS_MIXING
#ifdef GLS_MIXING
# undef  LMD_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
#endif

#undef BIOLOGY
#undef BIO_CSOMIO
#define OIL_BIO  /*coupled csomio oil model */
#if defined OIL_BIO
#  define BIOLOGY
#  define BIO_CSOMIO
#endif

#ifdef BIO_CSOMIO
# undef ANA_BIOLOGY
# undef CARBON
# undef DENITRIFICATION
# define OXYGEN
# undef BIO_SEDIMENT
# undef DIAGNOSTICS_BIO
# define ANA_BPFLUX
# define ANA_SPFLUX
#endif


#undef SEDIMENT
#define OIL_SEDIMENT
#ifdef OIL_SEDIMENT
# define SEDIMENT
#endif
#ifdef  SEDIMENT
# define SED_TOY
# define SUSPLOAD
# undef  BEDLOAD_SOULSBY
# undef  BEDLOAD_MPM
# undef  SED_DENS
# undef COHESIVE_BED
# define NONCOHESIVE_BED2
# undef  SED_MORPH
# undef RIVER_SEDIMENT
# ifdef SOLVE3D
#  define ANA_SEDIMENT
#  define ANA_BPFLUX
#  define ANA_BSFLUX
#  define ANA_BTFLUX
#  define ANA_SPFLUX
#  define ANA_SRFLUX
#  define ANA_SSFLUX
#  define ANA_STFLUX
# endif
#endif


#undef  BULK_FLUXES
#ifdef OFFLINE
# undef BULK_FLUXES
#endif
#ifdef BULK_FLUXES
# undef  QCORRECTION
# undef  LONGWAVE
# define LONGWAVE_OUT
# define SPECIFIC_HUMIDITY
# define SOLAR_SOURCE
# undef CLOUDS
#else
# undef  QCORRECTION
# undef  SOLAR_SOURCE
# undef  DIURNAL_SRFLUX
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX
#undef ANA_PERTURB

#undef FORWARD_WRITE
#define OUT_DOUBLE
#undef FORWARD_READ
#undef FORWARD_MIXING

/*
**----------------------------------------------
**  Adding FLoats and oil and Euler-Lagr mapping
**----------------------------------------------
**/
#define FLOATS
#if defined FLOATS && defined OFFLINE
# define OFFLINE_FLOATS
# undef OFFLINE_FLOATS_LATLON
# define FLOAT_VWALK
# ifdef FLOAT_VWALK
#  define AKTCLIMATOLOGY
# endif
#endif

#define FLOAT_OIL
#if !defined FLOATS
# undef FLOAT_OIL
#endif
#if defined FLOAT_OIL
# undef WOIL_INTEGRATED
# define OIL_DEBUG
# define OIL_MAP_DEBUG
# undef OIL_BIO_DEBUG
# define OIL_EULR 
# define OIL_WDRIFT      /* Surface wind drift - need wind fields */
#endif  /* FLOAT_OIL */
/* If oil is on - need T and S passed from climatology file */
#if defined OFFLINE && defined FLOAT_OIL
# define ATCLIMATOLOGY
#endif

#if defined OIL_BIO || defined OIL_SEDIMENT
# define OIL_EULR  /* EULER->/<- Lagrangian mapping */ 
#endif
#if !defined OIL_BIO
# undef OIL_BIO_DEBUG
#endif
#if defined OIL_SEDIMENT
# define SED_OPA
#endif
#if !defined FLOAT_OIL
# undef OIL_EULR
# undef OIL_BIO
# undef OIL_SEDIMENT
#endif

/*
**-----------------------------------------------------------------------------
**  Variational Data Assimilation.
**-----------------------------------------------------------------------------
*/

#ifdef NORMALIZATION
# undef  MULTIPLE_TLM
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  W4DVAR
# undef  R_SYMMETRY
# define CORRELATION
# undef  CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# undef  BALANCE_OPERATOR
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#if defined IS4DVAR || defined IS4DVAR_OLD
# undef  MULTIPLE_TLM
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  W4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# undef  CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# undef  BALANCE_OPERATOR
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef W4DVAR
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  IS4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# define CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# define RPM_RELAXATION
# undef  TLM_CHECK
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef W4DPSAS
# undef  AVERAGES
# undef  AVOID_ADJOINT
# undef  IS4DVAR
# undef  R_SYMMETRY
# undef  CORRELATION
# define CONVOLVE
# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  TLM_CHECK
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

#ifdef SANITY_CHECK
# define FULL_GRID
# define FORWARD_READ
# define FORWARD_WRITE
# define FORWARD_MIXING
# define OUT_DOUBLE
# define ANA_PERTURB
# define ANA_INITIAL
#endif

