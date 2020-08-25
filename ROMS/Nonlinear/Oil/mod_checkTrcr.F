#include "cppdefs.h"
      MODULE check_tracer_mod
#if defined OIL_DEBUG
!
!=============================================
! Check tracer arrays 
!=============================================
!
      implicit none

      PRIVATE
      PUBLIC :: check_tracer

      CONTAINS
!
!***********************************************************************
      SUBROUTINE check_tracer (ng,stl)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_parallel

      integer, intent(in) :: ng
      character(LEN=*), intent(in) :: stl

      integer :: lfl, ii, jj, kk, inn
      real(r8) :: dmm, tmin, tmax, smin, smax
      real(r8) :: tmn1,tmx1,tmn2,tmx2,tmn3,tmx3

!        DO ii=Istr(MyRank),Iend(MyRank)
!        DO jj=Jstr(MyRank),Jend(MyRank)
        DO inn=1,NT(ng)
          tmx1=-1.e30
          tmn1=1.e30
          tmx2=-1.e30
          tmn2=1.e30
          tmx3=-1.e30
          tmn3=1.e30
          DO kk=1,N(ng)
            dmm=minval(OCEAN(ng)%t(:,:,kk,1,inn))
            tmn1=minval((/tmn1,dmm/))
            dmm=maxval(OCEAN(ng)%t(:,:,kk,1,inn))
            tmx1=maxval((/tmx1,dmm/))
            dmm=minval(OCEAN(ng)%t(:,:,kk,2,inn))
            tmn2=minval((/tmn1,dmm/))
            dmm=maxval(OCEAN(ng)%t(:,:,kk,2,inn))
            tmx2=maxval((/tmx1,dmm/))
            dmm=minval(OCEAN(ng)%t(:,:,kk,3,inn))
            tmn3=minval((/tmn1,dmm/))
            dmm=maxval(OCEAN(ng)%t(:,:,kk,3,inn))
            tmx3=maxval((/tmx1,dmm/))
          ENDDO
          print*,stl,': check_tracer, MyRank =',MyRank,' Trc=',inn,     &
     &           'time1: min/max=',tmn1,tmx1,                           &
     &           'time2: min/max=',tmn2,tmx2,                           &
     &           'time3: min/max=',tmn3,tmx3
        ENDDO
      END SUBROUTINE check_tracer

#endif
      END MODULE check_tracer_mod
