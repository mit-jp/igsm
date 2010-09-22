C $Header: /home/igsm-repository/CHEM/chemtrop0.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	============================================================
!
!	CHEMTROP0.F:	Interface for subroutine CHEMTROP.F
!                             of MIT Global Chemistry Model
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
!	052300	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	092801	Chien Wang	add bc and oc
!	093001	Chien Wang	add S(VI) RH dependency
!	051904	Chien Wang	rev.
!
!	==========================================================

!
      subroutine chemtrop0(ifss, pT, qv, dtr, nloop)
!     =============================================

#include "chem_para"
#include "chem_com"

#include "BD2G04.COM"

	common U,V,T,P,Q
	
      dimension pT  (nlon,nlat,nlev)
      dimension Temp(nlon,nlat,nlev)
      dimension qv  (nlon,nlat,nlev)
      dimension den (nlon,nlat,nlev)
      dimension rh  (nlon,nlat,nlev)

      dimension tmp_co  (nlon,nlat,nlev)
      dimension tmp_ch4 (nlon,nlat,nlev)
      dimension tmp_o3  (nlon,nlat,nlev)
      dimension tmp_svi (nlon,nlat,nlev)
      dimension tmp_no  (nlon,nlat,nlev)
      dimension tmp_no2 (nlon,nlat,nlev)
      dimension tmp_nv  (nlon,nlat,nlev)
      dimension tmp_ch2o(nlon,nlat,nlev)

!	--------------------------------------------

#if ( defined CPL_CHEM )

 	ktrop = n_tropopause

c === 032697
c ===	add diagnostic procedure:
c
	do k=1,ktrop
	do j=1,nlat
	  tmp_co  (1,j,k) = co   (1,j,k)
	  tmp_ch4 (1,j,k) = ch4  (1,j,k)
	  tmp_o3  (1,j,k) = o3   (1,j,k)
	  tmp_svi (1,j,k) = h2so4(1,j,k)
	  tmp_no  (1,j,k) = xno  (1,j,k)
	  tmp_no2 (1,j,k) = xno2 (1,j,k)
	  tmp_nv  (1,j,k) = hno3 (1,j,k)
	  tmp_ch2o(1,j,k) = ch2o (1,j,k)
	enddo
	enddo

c---------
c Note the T from
c  main.f is a fraction of potential temprerature
c
      do k = 1, nlev
      do j = 1,n2dh 
	airptmp1    = (sig(k)*p(1,j) + 10.0)
        Temp(1,j,k) = T(1,j,k)*airptmp1**0.286
        den(1,j,k)  = airptmp1/(2.87*Temp(1,j,k))
	rh (1,j,k)  = 3.80/airpress(k)
     &      *exp(17.67*(Temp(1,j,k) - 273.15)
     &		/(Temp(1,j,k) - 29.65))
	rh (1,j,k)  = qv(1,j,k)/rh(1,j,k)*100.0
      end do
      end do

c     do 2 ntime =1,nloop

        call chemtrop(dtr, 0, ktrop, Temp, qv, den)

c2      continue

c === 032697
c ===	add diagnostic procedure:
c
	do k=1,ktrop
	do j=1,nlat
	  photo_co   (1,j,k) = photo_co  (1,j,k)
     &	     + (co   (1,j,k) -   tmp_co  (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_ch4  (1,j,k) = photo_ch4 (1,j,k)
     &	     + (ch4  (1,j,k) -   tmp_ch4 (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_o3   (1,j,k) = photo_o3  (1,j,k)
     &	     + (o3   (1,j,k) -   tmp_o3  (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_svi  (1,j,k) = photo_svi (1,j,k)
     &	     + (h2so4(1,j,k) -   tmp_svi (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_no   (1,j,k) = photo_no  (1,j,k)
     &	     + (xno  (1,j,k) -   tmp_no  (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_no2  (1,j,k) = photo_no2 (1,j,k)
     &	     + (xno2 (1,j,k) -   tmp_no2 (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_nv   (1,j,k) = photo_nv  (1,j,k)
     &	     + (hno3 (1,j,k) -   tmp_nv  (1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	  photo_ch2o (1,j,k) = photo_ch2o(1,j,k)
     &	     + (ch2o (1,j,k) -   tmp_ch2o(1,j,k))
     &	     *airmass(1,j,k)*1.e-18			!TGspecies
	enddo
	enddo


	i = 1

	do j=1,nlat
	  sviod(i,j,nlev ) = 0.0
	  sviod(i,j,nlev1) = 0.0
	  bcod (i,j,nlev ) = 0.0
	  bcod (i,j,nlev1) = 0.0
	  ocod (i,j,nlev ) = 0.0
	  ocod (i,j,nlev1) = 0.0
	end do

	do k=nlev1,1,-1
	do j=1,nlat
! =====
! Calculate optical depth of S(VI) aerosols:
!   ref. Charlson et al., 1992
!
	! Qex*f(rh) = 5.0*1.7 for rh = 80%
	!qex_svi = 8.5e-6/dxyp(j)
	!
	! === add frh based on calculated rh
	!
	if ( rh(i,j,k) .le. 60.0 ) then
	  frh = 1.0
	else if ( rh(i,j,k) .ge. 80.0 ) then
	  frh = 2.8
	else
	  frh = rh(i,j,k)
	  frh = -9.2906106183 
     &		+ frh*(0.52570211505
     &		+ frh*(-0.0089285760691+5.0877212432e-05*frh))
	end if
	qex_svi = 5.0e-6*frh/dxyp(j)
     &	

! === bc
	! Qex*f(rh) = 9.0*1.0 (550 micron)
	qex_bc = 8.0e-6/dxyp(j)	! normal
	!qex_bc = 11.0e-6/dxyp(j)	! high

! === oc
	! Qex*f(rh) = 6.8*1.0 (550 micron), assume rh factor
	qex_oc = 6.8e-6/dxyp(j)

	  sviod(i,j,k) = sviod(i,j,k+1) 
     &                 + airmass(i,j,k)*h2so4(i,j,k)*qex_svi
!         if( j.eq.33.and.k.eq.1) then 
!          print *,airmass(i,j,k),rh(i,j,k),qex_svi
!          print *,h2so4(i,j,k),sviod(i,j,k+1),sviod(i,j,k)
!         endif
	  bcod(i,j,k) = bcod(i,j,k+1) 
     &                 + airmass(i,j,k)*bcarbon(i,j,k)*qex_bc
	  ocod(i,j,k) = ocod(i,j,k+1) 
     &                 + airmass(i,j,k)*ocarbon(i,j,k)*qex_oc
	end do
	end do
#endif

       return
        end

