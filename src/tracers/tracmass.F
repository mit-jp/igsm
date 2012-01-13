
#include "ctrparam.h"

!	============================================================
!
!	TRACAIRMASS.F:  Subroutine for calculating air mass
!			in MIT Global Chemistry Model
!
!	------------------------------------------------------------
!
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!				of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------
!	122095	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	rev.
!	091811	A     S   	rev.
!
!	==========================================================

!	=========================
	Subroutine tracairmass(p)
!	=========================

#include "chem_para"
#include "chem_com"
!#include "BD2G04.COM"

          dimension p(nlon,nlat)

!	----------------------------------------------------------


c-------------------------------
c  Calculate air mass:
c
!      print *,'From TRACAIEMASS'
!      print *,'P(1,23)=',p(1,23)
      
      i= 1
      do k=1,nlev
      do j=1,nlat
        airmass(i,j,k) = airmass0(i,j,k)*p(i,j)
      end do
      end do


      return
       end

!	================================
	Subroutine tracmass0
!	================================

c==================================================================c
c                                                                  c
c   CHEMMASS.F:  Subroutine for calculating total amount (mass)   c
c                   of tracers in MIT Global Climate Model       c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    August 8, 1995                             c
c   Last Revised on:    Septrmber 14, 2011                             c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"
#include "tracer_com.h"


!	--------------------------------------------------------


c-------------------------------
c  Calculate total amount of tracer:
c
!     print *,'From TRACERSMASS0 ntracers=',ntracers
       i = 1
      do n=1,ntracers
!        print *,airmass(i,38,7),tracers(i,38,7,n)
       tracmass(n)= 0.0
       do k=1,nlev
       do j=1,nlat
!        print *,k,j,airmass(i,j,k)
!    &          , tracers(i,j,k,n)
        tracmass(n)= tracmass(n)+
     &          + airmass(i,j,k)
     &          * tracers(i,j,k,n)
       end do
       end do
      end do



      return
       end

!	=======================================
	Subroutine tracmassadj
!	=======================================

c==================================================================c
c                                                                  c
c   TRACMASSADJ.F:  Subroutine for conpensating mass loss during     c
c		    simulation based on unified mapping 	   c
c                   in MIT Global Climate Model 	           c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    September 15, 1995                         c
c   Last Revised on:    September 14, 2011                         c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"
#include "tracer_com.h"


!	-------------------------------------------------


c---------------------------
c  Readjust tracer's mass:
c
      do n=1,ntracers
       tracmass1(n)= 0.0
       i = 1
       do k=1,nlev
       do j=1,nlat
        tracmass1(n)= tracmass1(n)+
     &          + airmass(i,j,k)
     &          *  max(0.0,tracers(i,j,k,n))
       end do
       end do


       if(tracmass1(n).gt.0.0) then
        xratio = tracmass(n)/tracmass1(n)
       else
         xratio = 1.0
       endif

       do k=1,nlev
       do j=1,nlat

        tracers(i,j,k,n)= tracers(i,j,k,n) * xratio

       end do
       end do

      end do  ! ntracers

      return
       end
	
!	=============================================
	Subroutine tracmassloss
!	=============================================

c==================================================================c
c                                                                  c
c   TRACMASSLOSS.F:  Subroutine                                    c
c		     take away lossed mass directly according      c
c			to the tropospheric life time              c
c                   in MIT Global Chemistry Model 	           c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    September 15, 1995                         c
c   Last Revised on:    September 14, 2011                         c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"
#include "tracer_com.h"


!	----------------------------------------------------


c---------------------------
c

      do n=1,ntracers
!      if( trtype(n).eq.1.or.trtype(n).eq.2 ) then
       if( treftime(n).gt.0.0) then

        i=1
        do k=1,nlev
        do j=1,nlat

         tracers(i,j,k,n)= tracers(i,j,k,n) * 
     &     (1. - 1./(treftime(n)*24.) )  ! treftime is in days

        end do
        end do

       endif
      end do


      return
       end

