C $Header: /home/igsm-repository/ATM/check_lib.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	This subroutine saves library settings
!
      subroutine check_lib
      common /libset/nlev,nlat,nlon0

c
       nlev=N_LEV
       nlat=N_LAT
       nlon0=N_LON0
       print *,' FROM LIB',nlev,nlat,nlon0
       if(nalt.eq.24.and.nlon1.ne.36)then
         print *,'From Clim Lib'
         print *,' N_LAT inconsistent with N_LON0'
         print *,nlat,nlon0
         stop
       endif
       if(nalt.eq.46.and.nlon1.ne.72)then
         print *,'From Clim Lib'
         print *,' N_LAT inconsistent with N_LON0'
         print *,nlat,nlon1
         stop
       endif
c      print *,'nlev=',nlev
C
      return
      end
