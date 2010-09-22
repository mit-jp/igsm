
#include "ctrparam.h"

!	This subroutine checks model settings
!
      subroutine check_settings
      common/libset/nlev_l,nlat_l,nlon_l

c       Checking of consistency of model settings
c
       LMO=N_LEVOCEAN
       nlev=N_LEV
       nlat=N_LAT
       nlon0=N_LON0
       print *,'From check_settings',nlev,nlat,nlon0
       if(nlat.eq.24.and.nlon0.ne.36)then
         print *,'From model'
         print *,' N_LAT inconsistent with N_LON0'
         print *,N_LAT,N_LON0
         stop
       endif
       if(nlat.eq.46.and.nlon0.ne.72)then
         print *,'From model'
         print *,' N_LAT inconsistent with N_LON0'
         print *,N_LAT,N_LON0
         stop
       endif

       call check_lib
       if(nlev.ne.nlev_l)then
         print *,' N_LEV inconsistent with N_LEV_LIB'
         print *,'N_LEV_MODEL=',nlev,' N_LEV_LIB=',nlev_l
         stop
       endif
       if(nlat.ne.nlat_l)then
         print *,' N_LAT inconsistent with N_LAT_LIB'
         print *,'N_LAT_MODEL=',nlat,' N_LAT_LIB=',nlat_l
         stop
       endif
       print *,'nlev=',nlev
C
C   RUNS WITH CHEMISTRY
C
#if ( defined CPL_CHEM )
!
        print *,' Climate model with interactive chemistry'
!
#  if ( !defined PREDICTED_GASES )
c
       print *, 'PREDICTED_GASES undefined'
       stop
c
#  endif
#  if ( !defined PREDICTED_AEROSOL )
c
       print *, 'PREDICTED_AEROSOL undefined'
       stop
c
#  endif
#  if ( !defined PREDICTED_BC )
c
       print *, 'PREDICTED_BC undefined'
       stop
c
#  endif
#  if ( defined O3_FORCED) || (defined SVI_ALBEDO)
       print *,' O3_FORCED and SVI_ALBEDO should be undefined'
       print *,O3_FORCED,SVI_ALBEDO 
       stop
#  endif
c
#  if ( defined CPL_OCEANCO2 )
      print *,' Run with ocean chemistry'
#  if ( defined ML_2D )
       if(LMO.ne.11) then
         print *,'Wrong LMO for OCM model',LMO
         stop
       endif
#  endif
#  endif
#  if (defined CPL_META )
       print *,'Run with META model'
#  endif
C
C   RUNS WITHOUT CHEMISTRY
C
#else 
       print *,'Climate model without chemistry'

#  if ( defined PREDICTED_BC )
       print *,'Run with chahges in black carbon predicted by'
       print *,' climate-chemistry model READGHG should be 1 or 2'
      
#  endif
c
c
#  if (defined CPL_META )
c
       print *, 'CPL_META should be also undefined'
       print *,  CPL_META 
       stop
#  endif
#  if ( defined PREDICTED_GASES )
       print *,'Run with chahges in GHGs predicted by'
       print *,' climate-chemistry model READGHG should be 1 or 2'
#    if ( defined O3_FORCED) || (defined SVI_ALBEDO)
       print *,' O3_FORCED and SVI_ALBEDO should be undefined'
       print *,O3_FORCED,SVI_ALBEDO 
       stop
#    endif
#  endif
#  if ( defined PREDICTED_AEROSOL )
       print *,'Run with chahges in aerosol predicted by'
       print *,' climate-chemistry model READGHG should be  1 or 2'
#    if ( defined O3_FORCED) || (defined SVI_ALBEDO)
       print *,' O3_FORCED and SVI_ALBEDO should be undefined'
       print *,O3_FORCED,SVI_ALBEDO 
       stop
#    endif
#  endif
#  if ( defined O3_FORCED) 
       print *,' Run with prescribed changes in ozone'
c      if(NLEN.ne.11) then
c       print *,'NLEN should be 11'
c       print *,'NLEV=',N_LEV
c       stop
c      endif
       print *,'Make sure that OBSFOR = true'
#  endif
#  if (defined SVI_ALBEDO)
       print *,' Run with prescribed changes in aerosol (surf. albedo)'
c      if(NLEN.ne.11) then
c       print *,'NLEN should be 11'
c       print *,'NLEV=',N_LEV
c       stop
c      endif
       print *,'Make sure that OBSFOR = true'
#  endif
c
#endif
#if ( defined CPL_TEM )
        print *,' Run with TEM'
#endif
#if ( defined CPL_OCEANCO2 )
      print *,' Run with ocean chemistry'
#if ( defined ML_2D )
       if(LMO.ne.11) then
         print *,'Wrong LMO for OCM model',LMO
         stop
       endif
#endif
#endif
#if ( defined CPL_NEM )
        print *,' Run with NEM'
#endif
      return
      end
