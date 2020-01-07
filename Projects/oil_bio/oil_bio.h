/*
** svn $Id: ias.h 830 2017-01-24 21:21:11Z arango $
*******************************************************************************
** Copyright (c) 2002-2018 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for NGOM 1/25 DEG resolution.
**
** Application flag:   oil_bio
** Input script:       ocean_oil_bio.in
*/

#undef WRF_MODEL
#undef SWAN_MODEL
#define ROMS_MODEL

#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */
#undef  VERIFICATION            /* NL Observation Verification Driver */
#undef  NORMALIZATION           /* Background error Covariance Normalization */
#undef  AD_SENSITIVITY          /* Adjoint Sensitivity Driver */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state settings.
**-----------------------------------------------------------------------------
*/
#undef  AVERAGES               /*Write out time-averaged data*/
#undef  AVERAGES_FLUXES
#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV               /*mixing of momentum along constant s surfaces*/
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define UV_U3HADVECTION
#undef TS_U3HADVECTION
#define TS_MPDATA
#undef TS_C4VADVECTION
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

#define BIOLOGY
#define BIO_CSOMIO

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
#ifdef  SEDIMENT
# define SUSPLOAD
# undef  BEDLOAD_SOULSBY
# undef  BEDLOAD_MPM
# undef  SED_DENS
# undef COHESIVE_BED
# undef  SED_MORPH
# define RIVER_SEDIMENT
# ifdef SOLVE3D
#  undef ANA_SEDIMENT
#  define ANA_BPFLUX
#  define ANA_BSFLUX
#  define ANA_BTFLUX
#  define ANA_SPFLUX
#  define ANA_SRFLUX
#  define ANA_SSFLUX
#  define ANA_STFLUX
# endif
#endif

#define  BULK_FLUXES
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
#endif

#undef ANA_BSFLUX
#undef ANA_BTFLUX
#undef ANA_PERTURB

#define FORWARD_WRITE
#undef OUT_DOUBLE
#undef FORWARD_READ
#undef FORWARD_MIXING

#define FLOATS
#define FLOAT_VWALK
#define FLOAT_OIL
#undef WOIL_INTEGRATED
#define OIL_DEBUG
#define OIL_MAP_DEBUG
#undef OIL_BIO
#undef OIL_SEDIMENT
#if defined OIL_BIO || defined OIL_SEDIMENT
#  define OIL_EULR  /* EULER->/<- Lagrangian mapping */ 
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
