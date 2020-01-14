/*
** svn $Id: fennel_var.h 830 2017-01-24 21:21:11Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2018 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Fennel et al. (2006) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iLDeN)')
                idTvar(iLDeN)=varid
              CASE ('idTvar(iSDeN)')
                idTvar(iSDeN)=varid
              CASE ('idTvar(iChlo)')
                idTvar(iChlo)=varid
# ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
              CASE ('idTvar(iLDeC)')
                idTvar(iLDeC)=varid
              CASE ('idTvar(iSDeC)')
                idTvar(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
# endif
/* XCHEN MMMM */
# ifdef ADDZX
              CASE ('idTvar(iZx)')
                idTvar(iZx)=varid
# endif
/* XCHEN WWWW */

/*XCHEN OILBUGS MMMM*/
/*# ifdef OILBUGS */
/* DDMITRY */
# ifdef OIL_BIO
                CASE ('idTvar(iSbug1)')
                  idTvar(iSbug1)=varid
                CASE ('idTvar(iSbug2)')
                  idTvar(iSbug2)=varid
                CASE ('idTvar(iSbug3)')
                  idTvar(iSbug3)=varid
                CASE ('idTvar(iSbug4)')
                  idTvar(iSbug4)=varid
                CASE ('idTvar(iLbug1)')
                  idTvar(iLbug1)=varid
                CASE ('idTvar(iLbug2)')
                  idTvar(iLbug2)=varid
                CASE ('idTvar(iLbug3)')
                  idTvar(iLbug3)=varid
                CASE ('idTvar(iLbug4)')
                  idTvar(iLbug4)=varid
                CASE ('idTvar(iSat)')
                  idTvar(iSat)=varid
                CASE ('idTvar(iArmt)')
                  idTvar(iArmt)=varid
                CASE ('idTvar(iResin)')
                  idTvar(iResin)=varid
# endif
/*XCHEN OILBUGS WWWW*/

/*
**  Adjoint sensitivity state biological tracers.
*/

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iNH4_)')
                idTads(iNH4_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iLDeN)')
                idTads(iLDeN)=varid
              CASE ('idTads(iSDeN)')
                idTads(iSDeN)=varid
              CASE ('idTads(iChlo)')
                idTads(iChlo)=varid
# ifdef CARBON
              CASE ('idTads(iTIC_)')
                idTads(iTIC_)=varid
              CASE ('idTads(iTAlk)')
                idTads(iTAlk)=varid
              CASE ('idTads(iLDeC)')
                idTads(iLDeC)=varid
              CASE ('idTads(iSDeC)')
                idTads(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
# endif
/* XCHEN MMMM */
# ifdef ADDZX
              CASE ('idTads(iZx)')
                idTads(iZx)=varid
# endif
/* XCHEN WWWW */
#endif

/*XCHEN OILBUGS MMMM*/
/*# ifdef OILBUGS */ /* DDMITRY */
# ifdef OIL_BIO 
                CASE ('idTads(iSbug1)')
                  idTads(iSbug1)=varid
                CASE ('idTads(iSbug2)')
                  idTads(iSbug2)=varid
                CASE ('idTads(iSbug3)')
                  idTads(iSbug3)=varid
                CASE ('idTads(iSbug4)')
                  idTads(iSbug4)=varid
                CASE ('idTads(iLbug1)')
                  idTads(iLbug1)=varid
                CASE ('idTads(iLbug2)')
                  idTads(iLbug2)=varid
                CASE ('idTads(iLbug3)')
                  idTads(iLbug3)=varid
                CASE ('idTads(iLbug4)')
                  idTads(iLbug4)=varid
                CASE ('idTads(iSat)')
                  idTads(iSat)=varid
                CASE ('idTads(iArmt)')
                  idTads(iArmt)=varid
                CASE ('idTads(iResin)')
                  idTads(iResin)=varid
# endif
/*XCHEN OILBUGS WWWW*/

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid

              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid

              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid

              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid

              CASE ('idTbry(iwest,iSDeN)')
                idTbry(iwest,iSDeN)=varid
              CASE ('idTbry(ieast,iSDeN)')
                idTbry(ieast,iSDeN)=varid
              CASE ('idTbry(isouth,iSDeN)')
                idTbry(isouth,iSDeN)=varid
              CASE ('idTbry(inorth,iSDeN)')
                idTbry(inorth,iSDeN)=varid

              CASE ('idTbry(iwest,iLDeN)')
                idTbry(iwest,iLDeN)=varid
              CASE ('idTbry(ieast,iLDeN)')
                idTbry(ieast,iLDeN)=varid
              CASE ('idTbry(isouth,iLDeN)')
                idTbry(isouth,iLDeN)=varid
              CASE ('idTbry(inorth,iLDeN)')
                idTbry(inorth,iLDeN)=varid

              CASE ('idTbry(iwest,iChlo)')
                idTbry(iwest,iChlo)=varid
              CASE ('idTbry(ieast,iChlo)')
                idTbry(ieast,iChlo)=varid
              CASE ('idTbry(isouth,iChlo)')
                idTbry(isouth,iChlo)=varid
              CASE ('idTbry(inorth,iChlo)')
                idTbry(inorth,iChlo)=varid

#ifdef CARBON
              CASE ('idTbry(iwest,iSDeC)')
                idTbry(iwest,iSDeC)=varid
              CASE ('idTbry(ieast,iSDeC)')
                idTbry(ieast,iSDeC)=varid
              CASE ('idTbry(isouth,iSDeC)')
                idTbry(isouth,iSDeC)=varid
              CASE ('idTbry(inorth,iSDeC)')
                idTbry(inorth,iSDeC)=varid

              CASE ('idTbry(iwest,iLDeC)')
                idTbry(iwest,iLDeC)=varid
              CASE ('idTbry(ieast,iLDeC)')
                idTbry(ieast,iLDeC)=varid
              CASE ('idTbry(isouth,iLDeC)')
                idTbry(isouth,iLDeC)=varid
              CASE ('idTbry(inorth,iLDeC)')
                idTbry(inorth,iLDeC)=varid

              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid

              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
#endif
#ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
#endif
/* XCHEN MMMM */
#ifdef ADDZX
              CASE ('idTbry(iwest,iZx)')
                idTbry(iwest,iZx)=varid
              CASE ('idTbry(ieast,iZx)')
                idTbry(ieast,iZx)=varid
              CASE ('idTbry(isouth,iZx)')
                idTbry(isouth,iZx)=varid
              CASE ('idTbry(inorth,iZx)')
                idTbry(inorth,iZx)=varid
#endif
/* XCHEN WWWW */

/*XCHEN OILBUGS MMMM*/ 
/*DDMITRY */
# ifdef OIL_BIO 
            CASE ('idTbry(iwest,iSbug1)')
              idTbry(iwest,iSbug1)=varid
            CASE ('idTbry(ieast,iSbug1)')
              idTbry(ieast,iSbug1)=varid
            CASE ('idTbry(isouth,iSbug1)')
              idTbry(isouth,iSbug1)=varid
            CASE ('idTbry(inorth,iSbug1)')
              idTbry(inorth,iSbug1)=varid

            CASE ('idTbry(iwest,iSbug2)')
              idTbry(iwest,iSbug2)=varid
            CASE ('idTbry(ieast,iSbug2)')
              idTbry(ieast,iSbug2)=varid
            CASE ('idTbry(isouth,iSbug2)')
              idTbry(isouth,iSbug2)=varid
            CASE ('idTbry(inorth,iSbug2)')
              idTbry(inorth,iSbug2)=varid

            CASE ('idTbry(iwest,iSbug3)')
              idTbry(iwest,iSbug3)=varid
            CASE ('idTbry(ieast,iSbug3)')
              idTbry(ieast,iSbug3)=varid
            CASE ('idTbry(isouth,iSbug3)')
              idTbry(isouth,iSbug3)=varid
            CASE ('idTbry(inorth,iSbug3)')
              idTbry(inorth,iSbug3)=varid

            CASE ('idTbry(iwest,iSbug4)')
              idTbry(iwest,iSbug4)=varid
            CASE ('idTbry(ieast,iSbug4)')
              idTbry(ieast,iSbug4)=varid
            CASE ('idTbry(isouth,iSbug4)')
              idTbry(isouth,iSbug4)=varid
            CASE ('idTbry(inorth,iSbug4)')
              idTbry(inorth,iSbug4)=varid

            CASE ('idTbry(iwest,iLbug1)')
              idTbry(iwest,iLbug1)=varid
            CASE ('idTbry(ieast,iLbug1)')
              idTbry(ieast,iLbug1)=varid
            CASE ('idTbry(isouth,iLbug1)')
              idTbry(isouth,iLbug1)=varid
            CASE ('idTbry(inorth,iLbug1)')
              idTbry(inorth,iLbug1)=varid

            CASE ('idTbry(iwest,iLbug2)')
              idTbry(iwest,iLbug2)=varid
            CASE ('idTbry(ieast,iLbug2)')
              idTbry(ieast,iLbug2)=varid
            CASE ('idTbry(isouth,iLbug2)')
              idTbry(isouth,iLbug2)=varid
            CASE ('idTbry(inorth,iLbug2)')
              idTbry(inorth,iLbug2)=varid

            CASE ('idTbry(iwest,iLbug3)')
              idTbry(iwest,iLbug3)=varid
            CASE ('idTbry(ieast,iLbug3)')
              idTbry(ieast,iLbug3)=varid
            CASE ('idTbry(isouth,iLbug3)')
              idTbry(isouth,iLbug3)=varid
            CASE ('idTbry(inorth,iLbug3)')
              idTbry(inorth,iLbug3)=varid

            CASE ('idTbry(iwest,iLbug4)')
              idTbry(iwest,iLbug4)=varid
            CASE ('idTbry(ieast,iLbug4)')
              idTbry(ieast,iLbug4)=varid
            CASE ('idTbry(isouth,iLbug4)')
              idTbry(isouth,iLbug4)=varid
            CASE ('idTbry(inorth,iLbug4)')
              idTbry(inorth,iLbug4)=varid

            CASE ('idTbry(iwest,iSat)')
              idTbry(iwest,iSat)=varid
            CASE ('idTbry(ieast,iSat)')
              idTbry(ieast,iSat)=varid
            CASE ('idTbry(isouth,iSat)')
              idTbry(isouth,iSat)=varid
            CASE ('idTbry(inorth,iSat)')
              idTbry(inorth,iSat)=varid

            CASE ('idTbry(iwest,iArmt)')
              idTbry(iwest,iArmt)=varid
            CASE ('idTbry(ieast,iArmt)')
              idTbry(ieast,iArmt)=varid
            CASE ('idTbry(isouth,iArmt)')
              idTbry(isouth,iArmt)=varid
            CASE ('idTbry(inorth,iArmt)')
              idTbry(inorth,iArmt)=varid

            CASE ('idTbry(iwest,iResin)')
              idTbry(iwest,iResin)=varid
            CASE ('idTbry(ieast,iResin)')
              idTbry(ieast,iResin)=varid
            CASE ('idTbry(isouth,iResin)')
              idTbry(isouth,iResin)=varid
            CASE ('idTbry(inorth,iResin)')
              idTbry(inorth,iResin)=varid
# endif
/*XCHEN OILBUGS WWWW*/

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iLDeN)')
                idRtrc(iLDeN)=varid
              CASE ('idRtrc(iSDeN)')
                idRtrc(iSDeN)=varid
              CASE ('idRtrc(iChlo)')
                idRtrc(iChlo)=varid
#ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
              CASE ('idRtrc(iLDeC)')
                idRtrc(iLDeC)=varid
              CASE ('idRtrc(iSDeC)')
                idRtrc(iSDeC)=varid
#endif
#ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
#endif
/* XCHEN MMMM */
#ifdef ADDZX
              CASE ('idRtrc(iZx)')
                idRtrc(iZx)=varid
#endif
/* XCHEN WWWW */

/*XCHEN OILBUGS MMMM*/
/*# ifdef OILBUGS*/
/* DDMITRY */
# ifdef OIL_BIO

              CASE ('idRtrc(iSbug1)')
                idRtrc(iSbug1)=varid
              CASE ('idRtrc(iSbug2)')
                idRtrc(iSbug2)=varid
              CASE ('idRtrc(iSbug3)')
                idRtrc(iSbug3)=varid
              CASE ('idRtrc(iSbug4)')
                idRtrc(iSbug4)=varid
              CASE ('idRtrc(iLbug1)')
                idRtrc(iLbug1)=varid
              CASE ('idRtrc(iLbug2)')
                idRtrc(iLbug2)=varid
              CASE ('idRtrc(iLbug3)')
                idRtrc(iLbug3)=varid
              CASE ('idRtrc(iLbug4)')
                idRtrc(iLbug4)=varid
              CASE ('idRtrc(iSat)')
                idRtrc(iSat)=varid
              CASE ('idRtrc(iArmt)')
                idRtrc(iArmt)=varid
              CASE ('idRtrc(iResin)')
                idRtrc(iResin)=varid
# endif
/*XCHEN OILBUGS WWWW*/

#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/
# ifdef DENITRIFICATION
              CASE ('iDbio2(iDNIT)')
                iDbio2(iDNIT)=varid
# endif
# ifdef CARBON
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(ipCO2)')
                iDbio2(ipCO2)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
# endif
              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
#endif
