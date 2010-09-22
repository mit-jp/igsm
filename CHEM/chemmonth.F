C $Header: /home/igsm-repository/CHEM/chemmonth.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	============================================================
!
!	CHEMMONTH1.F:   Subroutine for accumulating tracer mixing
!				ratios prepared for monthly averaging
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
!	052000	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	091901	Chien Wang	fix "write(180,*),"
!	092001	Chine Wang	add bc and oc
!	051804	Chien Wang	rev.
!
!	==========================================================

!	=====================
 	subroutine chemmonth1
!	=====================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

!	------------------------------------------------------

#if ( defined CPL_CHEM )

      do i=1,n3d
        cfc11m(i,1,1) = cfc11m(i,1,1) + cfc11(i,1,1)
        cfc12m(i,1,1) = cfc12m(i,1,1) + cfc12(i,1,1)
        xn2om (i,1,1) = xn2om (i,1,1) + xn2o (i,1,1)
        o3m   (i,1,1) = o3m   (i,1,1) + o3   (i,1,1)
        com   (i,1,1) = com   (i,1,1) + co   (i,1,1)
        zco2m (i,1,1) = zco2m (i,1,1) + zco2 (i,1,1)
        hoxm  (i,1,1) = hoxm  (i,1,1) 
c     &                + atomh (i,1,1)
     &                + ho    (i,1,1)
c     &                + ho2   (i,1,1)
        znox          = xno   (i,1,1) + xno2 (i,1,1)
        xnoxm (i,1,1) = xnoxm (i,1,1) + znox
        xnoym (i,1,1) = xnoym (i,1,1) + znox
     &                + xn2o5 (i,1,1)
     &                + hno3  (i,1,1)
        ch4m  (i,1,1) = ch4m  (i,1,1) + ch4(i,1,1)
        so2m  (i,1,1) = so2m  (i,1,1) + so2(i,1,1)
        h2so4m(i,1,1) = h2so4m(i,1,1) + h2so4(i,1,1)
        sviodm(i,1,1) = sviodm(i,1,1) + sviod(i,1,1)

#ifdef INC_3GASES
        ! === if hfc, pfc, and sf6 are included:
	! === 032698
        hfc134am(i,1,1) = hfc134am(i,1,1) + hfc134a(i,1,1)
        pfcm    (i,1,1) = pfcm    (i,1,1) + pfc(i,1,1)
        sf6m    (i,1,1) = sf6m    (i,1,1) + sf6(i,1,1)
	! ===
#endif
	bcm     (i,1,1) = bcm     (i,1,1) + bcarbon(i,1,1)
	ocm     (i,1,1) = ocm     (i,1,1) + ocarbon(i,1,1)
        bcodm   (i,1,1) = bcodm   (i,1,1) + bcod   (i,1,1)
        ocodm   (i,1,1) = ocodm   (i,1,1) + ocod   (i,1,1)
      end do

      monthstep=monthstep+1

#endif

      return
       end

!	=====================
	subroutine chemmonth2
!	=====================

!	=============================================================
!
!	CHEMMONTH2.F:   Subroutine for calculating monthly averaged
!				mixing ratios of tracers
!	-------------------------------------------------------------
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!			    of Global Change
!	Last Revised:	June 29, 1999
!
!	=============================================================

#include "chem_para"
#include "chem_const1"
#include "chem_com"
#include "BD2G04.COM"

#if ( defined CPL_META )
#include "chem_meta"
#endif

!	---------------------------------------------------------

#if ( defined CPL_CHEM )

	haha=1./float(monthstep)

c calculate monthly averaged values:

      do i=1,n3d
        cfc11m(i,1,1) = cfc11m(i,1,1)*haha
        cfc12m(i,1,1) = cfc12m(i,1,1)*haha
        xn2om (i,1,1) = xn2om (i,1,1)*haha
        o3m   (i,1,1) = o3m   (i,1,1)*haha
        com   (i,1,1) = com   (i,1,1)*haha
        zco2m (i,1,1) = zco2m (i,1,1)*haha
        hoxm  (i,1,1) = hoxm  (i,1,1)*haha
        xnoxm (i,1,1) = xnoxm (i,1,1)*haha
        xnoym (i,1,1) = xnoym (i,1,1)*haha
        ch4m  (i,1,1) = ch4m  (i,1,1)*haha
        so2m  (i,1,1) = so2m  (i,1,1)*haha
        h2so4m(i,1,1) = h2so4m(i,1,1)*haha
        sviodm(i,1,1) = sviodm(i,1,1)*haha

#ifdef INC_3GASES
	! === if hfc, pfc, and sf6 are included:
	! === 032698:
        hfc134am(i,1,1) = hfc134am(i,1,1)*haha
        pfcm    (i,1,1) = pfcm    (i,1,1)*haha
        sf6m    (i,1,1) = sf6m    (i,1,1)*haha
#endif
        bcm 	(i,1,1) = bcm     (i,1,1)*haha
        ocm 	(i,1,1) = ocm     (i,1,1)*haha
        bcodm   (i,1,1) = bcodm   (i,1,1)*haha
        ocodm   (i,1,1) = ocodm   (i,1,1)*haha
      end do

 	write(169)cfc11m
 	write(169)cfc12m
	write(169)xn2om
	write(169)o3m
	write(169)com
	write(169)zco2m
	write(169)hoxm
	write(169)xnoxm
	write(169)xnoym
	write(169)ch4m
 	write(169)so2m
 	write(169)h2so4m
	write(169)sviodm
	write(169)bcm
	write(169)ocm
	write(169)bcodm
	write(169)ocodm

#ifdef INC_3GASES
	! === if hfc, pfc, and sf6 are included:
	! === 032698:
  	write(179)hfc134am
	write(179)pfcm
	write(179)sf6m
#endif

c === 032697
c ===	add diagnostic output:
c
	write(177)photo_co
	write(177)photo_ch4
	write(177)photo_o3
	write(177)photo_svi
	write(177)photo_no
	write(177)photo_no2
	write(177)photo_nv
	write(177)photo_ch2o

#if ( defined CPL_META )
!
! === 020999
! === monthly avaraged meta model results
!
	haha = 1./float(nstep_meta)
	do j=1,nlat
	do ntype=1,3
	do i=1,meta_nvar
	  results_mon(i,ntype,j) = results_mon (i,ntype,j)*haha
	end do
	end do
	end do
	nstep_meta = 0

	write(181)results_mon
	
	do j=1,nlat
	do ntype=1,3
	do i=1,meta_nvar
	  results_mon(i,ntype,j) = 0.0
	end do
	end do
	end do
#endif
	
c calculate tropospheric mass-averaged mixing ratios:

!
! --- NOTE: Currently for N_LEV == 9 & N_LEV == 11
!	troposphere is defined from 1 to 7, therefore
!	no cpp control is applied here
!	
      tropmass = 0.0
      do j=1,nlat
      do k=1,n_tropopause
        tropmass = tropmass + airmass(1,j,k)
      enddo
      enddo
	
      globalmass = tropmass
      do j=1,nlat
      do k=n_tropopause+1,nlev
        globalmass = globalmass + airmass(1,j,k)
      enddo
      enddo

      tropmass     = 28.97296245/tropmass
      globalmass_m = 1./globalmass
      globalmass   = 28.97296245/globalmass

      cfc11global = 0.0
      cfc12global = 0.0
      xn2oglobal  = 0.0
      o3global    = 0.0
      zco2trop    = 0.0
      coglobal    = 0.0
      ch4global   = 0.0
      bcglobal    = 0.0
      ocglobal    = 0.0

#ifdef INC_3GASES
! === 032698:
      hfc134aglobal = 0.0
      pfcglobal     = 0.0
      sf6global     = 0.0
#endif

      do k=1,n_tropopause
      do j=1,nlat
        cfc11global = cfc11global
     &              + airmass(1,j,k)
     &               *cfc11m (1,j,k)
        cfc12global = cfc12global
     &              + airmass(1,j,k)
     &               *cfc12m (1,j,k)
        xn2oglobal  = xn2oglobal
     &              + airmass(1,j,k)
     &               *xn2om  (1,j,k)
        o3global    = o3global
     &              + airmass(1,j,k)
     &               *o3m    (1,j,k)
        zco2trop    = zco2trop
     &              + airmass(1,j,k)
     &               *zco2m  (1,j,k)
        coglobal    = coglobal
     &              + airmass(1,j,k)
     &               *com    (1,j,k)
        ch4global   = ch4global
     &              + airmass(1,j,k)
     &               *ch4m   (1,j,k)
        bcglobal    = bcglobal
     &              + airmass(1,j,k)
     &               *bcm    (1,j,k)
        ocglobal    = ocglobal
     &              + airmass(1,j,k)
     &               *ocm    (1,j,k)

#ifdef INC_3GASES
	! === if hfc, pfc, and sf6 are included:
	! === 032698:
        hfc134aglobal = hfc134aglobal
     &                + airmass(1,j,k)
     &                 *hfc134am (1,j,k)
        pfcglobal     = pfcglobal
     &                + airmass(1,j,k)
     &                 *pfcm (1,j,k)
        sf6global     = sf6global
     &                + airmass(1,j,k)
     &                 *sf6m (1,j,k)
#endif
      enddo
      enddo

	zco2global = zco2trop
      do k=n_tropopause+1,nlev
      do j=1,nlat
        zco2global  = zco2global
     &              + airmass(1,j,k)
     &               *zco2m  (1,j,k)
        bcglobal    = bcglobal
     &              + airmass(1,j,k)
     &               *bcm    (1,j,k)
        ocglobal    = ocglobal
     &              + airmass(1,j,k)
     &               *ocm    (1,j,k)
      enddo
      enddo

        cfc11global = cfc11global*tropmass/137.3675*1.e3 !pptv
        cfc12global = cfc12global*tropmass/120.9054*1.e3 !pptv
        xn2oglobal  = xn2oglobal *tropmass/44.0000       !ppbv
        o3global    = o3global   *tropmass/47.9982       !ppbv
        zco2trop    = zco2trop   *tropmass/44.0098*1.e-3 !ppmv
        zco2global  = zco2global *globalmass/44.0098*1.e-3 !ppmv
        coglobal    = coglobal   *tropmass/28.0104       !ppbv
        ch4global   = ch4global  *tropmass/16.0426*1.e-3 !ppmv
	bcglobal    = bcglobal   *globalmass_m*1.e3	 !pptm
	ocglobal    = ocglobal   *globalmass_m*1.e3	 !pptm

c      write(176,101)cfc11global,cfc12global,xn2oglobal,
c     &             o3global,zco2trop,zco2global,coglobal,ch4global
c101   format(8f11.3)

c 020196:
!monthnumber = (myyear - 1)*12 + mymonth
	monthnumber = (iyearchem - 1)*12 + mymonth

      write(176,101)monthnumber,cfc11global,cfc12global,xn2oglobal,
     &             o3global,zco2trop,zco2global,coglobal,ch4global,
     &		   bcglobal,ocglobal
101   format(i6,10f10.3)

#ifdef INC_3GASES
	! ===== if hfc, pfc, and sf6 are included:
	! === 032698:
        hfc134aglobal = hfc134aglobal
     &		       *tropmass/awHFC134a*1.e3	!pptv
        pfcglobal     = pfcglobal
     &		       *tropmass/awPFmethane*1.e3	!pptv
        sf6global     = sf6global
     &		       *tropmass/awSF6*1.e3		!pptv

	write(180,102)monthnumber,
     &		   hfc134aglobal,pfcglobal,sf6global
102	format(i6,3f10.3)
#endif

      do i=1,n3d
        cfc11m(i,1,1) = 0.0
        cfc12m(i,1,1) = 0.0
        xn2om (i,1,1) = 0.0
        o3m   (i,1,1) = 0.0
        com   (i,1,1) = 0.0
        zco2m (i,1,1) = 0.0
        hoxm  (i,1,1) = 0.0
        xnoxm (i,1,1) = 0.0
        xnoym (i,1,1) = 0.0
        ch4m  (i,1,1) = 0.0
        so2m  (i,1,1) = 0.0
        h2so4m(i,1,1) = 0.0
        sviodm(i,1,1) = 0.0

#ifdef INC_3GASES
	! === if hfc, pfc, and sf6 are included:
	! === 032698:
        hfc134am(i,1,1) = 0.0
        pfcm    (i,1,1) = 0.0
        sf6m    (i,1,1) = 0.0
#endif

        bcm 	(i,1,1) = 0.0
        ocm 	(i,1,1) = 0.0
        bcodm 	(i,1,1) = 0.0
        ocodm 	(i,1,1) = 0.0

! === 032697:
        photo_co  (i,1,1) = 0.0
        photo_ch4 (i,1,1) = 0.0
        photo_o3  (i,1,1) = 0.0
        photo_svi (i,1,1) = 0.0
        photo_no  (i,1,1) = 0.0
        photo_no2 (i,1,1) = 0.0
        photo_nv  (i,1,1) = 0.0
        photo_ch2o(i,1,1) = 0.0
      enddo

      monthstep=0

#endif

      return
       end

     
