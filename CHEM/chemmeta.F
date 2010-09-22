C $Header: /home/igsm-repository/CHEM/chemmeta.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	============================================================
!
!	CHEMMETA.F:   Interface between CliChem and META Models.
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
!	062999	Chien Wang	using time evolution data
!					for urban partition
!	120799	Chien Wang	use real temperature in meta
!	052200	Chien Wang	use urban daily emission data
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	back to the v58
!
!	==========================================================

!	===================
	subroutine chemmeta
!	===================

	parameter (ktop = 2)
	
#include "chem_para"
#include "chem_const1"
	
#include "chem_com"
#include "chem_meta"	
#include "BD2G04.COM"

	common U,V,T,P,Q
	
	real Temp	
	real tmp(meta_nvar), conv0, conv, x11(ktop), xxx
	real r_so2,     r_co,      r_nox,    r_voc
	real ymass_no,  ymass_no2, ymass_o3, ymass_hno3, 
     &       ymass_n2o5,ymass_co,  ymass_ch2o
     	real yflux_no,  yflux_no2, yflux_o3, yflux_hno3, 
     &       yflux_n2o5,yflux_co,  yflux_ch2o
	real urban_area,rural_area,total_area

	integer ind_lat(nlat), id_lat, jid
	
! --------------------------------------------
!  rlati:	|lat| in degree				(0:65)
!  rtemp:	deviation from normal temperature 	(-10:10)
!  rclou:	tenth of mean cloud coverage		(0:10)
!  rmixi:	maximum mixing layer top in meter	(200:2000)
!  rtime:	residence time 200km/U in second	(21600:259200)
!  rso2:	emission of SO2 in kg/km^2/day 
!  rco:		emission of SO2 in kg/km^2/day
!  rnox:	emission of SO2 in kg/km^2/day
!  rvoc:	deviation from total daily VOC emission	(-0.7:0.7)
!  raqino
!  raqinoo
!  raqivoc
!  raqiozo
!  raqisox
!  results(meta_nvar)
!
! --- Order of results and others:
! ---
! ---   1:	NO   mass
! ---   2:	NO   flux
! ---   3:	NO2  mass
! ---   4:	NO2  flux
! ---   5:	O3   mass
! ---   6:	O3   flux
! ---   7:	HONO flux
! ---   8:	HNO3 flux
! ---   9:	N2O5 flux
! ---  10:	CO   mass
! ---  11:	CO   flux
! ---  12:	HCHO mass
! ---  13:	HCHO flux
! ---  14:	ALD2 mass
! ---  15:	MEK  mass
! ---  16:	PAN  flux
! ---  17:	ALKA mass
! ---  18:	ETHE mass
! ---  19:	ALKE mass
! ---  20:	TOLU mass
! ---  21:	AROM mass
! ---  22:	SO2  mass
! ---  23:	SO2  flux
! ---  24:	SO3  mass
! ---  25:	SO3  flux
! ---  26:	O3   max
! ---  27:	O3   peak hour
! ---  Daily averaged mass
! ---  28:	NO
! ---  29:	NO2
! ---  30:	O3
! ---  31:	HONO
! ---  32:	HNO3
! ---  33:	HNO4
! ---  34:	N2O5
! ---  35:	NO3
! ---  36:	HO2
! ---  37:	CO
! ---  38:	HCHO
! ---  39:	ALD2
! ---  40:	MEK
! ---  41:	MGLY
! ---  42:	PAN
! ---  43:	RO2
! ---  44:	MCO3
! ---  45:	ALKN
! ---  46:	ALKA
! ---  47:	ETHE
! ---  48:	ALKE
! ---  49:	TOLU
! ---  50:	AROM
! ---  51:	DIAL
! ---  52:	CRES
! ---  53:	NPHE
! ---  54:	H2O2
! ---  55:	MEOH
! ---  56:	ETOH
! ---  57:	SO2
! ---  58:	SO3
! ----------------------------------------------------
!
#if ( defined CPL_CHEM ) && ( defined CPL_META )

	! --- All metamodel calls and calculations are applied 
	! ---	to latitudes ranged from -39 to 63 degree
	! ---	or j=8,20 in GACM grid index, though arrays 
	! ---	are still indexed from 1 to nlat.
	id_lat = 0
	do j=1,nlat
	  ind_lat(j) = 0
	end do

#if ( N_LAT == 24 )
	do j=8,20
#endif
#if ( N_LAT == 46 )
	do j=14,39
#endif
	  if ( n_total_urban(j,myyear) .ne. 0 ) then
		id_lat          = id_lat + 1
		ind_lat(id_lat) = j
	  endif
	end do

	do 10 jid = 1,id_lat
	  j = ind_lat(jid)

	  do ntype = 1,3
	  do iii   = 1,meta_nvar
		results_meta(iii,ntype,j) = 0.0
	  end do
	  end do

	  cldcvr = 0.0
	  do k=2,6
	    xxx = chem_cldss(1,j,k) + chem_cldmc(1,j,k)
     &          - chem_cldss(1,j,k) * chem_cldmc(1,j,k)
!	    if(xxx.gt.1.0)   xxx    = 1.0
	    if(xxx.gt.0.9)   xxx    = 0.9
	    if(xxx.lt.0.0)   xxx    = 0.0
	    if(xxx.gt.cldcvr)cldcvr = xxx
	  end do

	  rtemp(j)  = 0.0

	  ! --- cloud coverage in tenth	  
	  rclou(j)  = cldcvr*10.0
	  
	  ! --- m, top of layer 2
	  rmixi(j)  = 1500.0			
	  
	  ! --- 21600 < time < 259200, fixed size of 200 km
	  rtime(j)  = 2.e5/max(0.771605, abs(pvv(1,j,1)))
	  if(rtime(j).lt. 21600.0) rtime(j) =  21600.0
	  if(rtime(j).gt.259200.0) rtime(j) = 259200.0

	  urban_area = float(n_total_urban(j,myyear))*4.e10
	  total_area = 1./dxyp(j)
	  urban_area = urban_area*total_area
	  rural_area = 1.0 - urban_area
	  if (rural_area .le. 0.0) then
		rural_area = 0.0
		urban_area = 1.0
	  end if
	  
	  tmass = 0.0                !Total air mass in kg
	  do k=1,ktop
	  	tmass = tmass + airmass(1,j,k)
	  enddo
          tmass = 1./tmass

	  ymass_no   = 0.0
	  ymass_no2  = 0.0
	  ymass_o3   = 0.0
	  ymass_hno3 = 0.0
	  ymass_n2o5 = 0.0
	  ymass_co   = 0.0
	  ymass_ch2o = 0.0
	  yflux_no   = 0.0
	  yflux_no2  = 0.0
	  yflux_o3   = 0.0
	  yflux_hno3 = 0.0
	  yflux_n2o5 = 0.0
	  yflux_co   = 0.0
	  yflux_ch2o = 0.0
	  
	  ! --- convert daily emission for a 95x95 km^2 "core"
	  ! 		from 10^-9 kg to kg/km^2
	  ! 		conv0 = 1.e-9/(95x95)/n
	  ! 
	  conv0      = 1.10803e-13
     &           /(float(n_total_urban(j,myyear)))	

	do 20 ntype=1,3	! 3 different types

	 if (n_urban(ntype,j,myyear).ne.0) then

	  ! --- decide the emission strength of 
	  ! ---  different types of cities
	  if (n_total_urban(j,myyear).lt.4)then
	 	conv = conv0*urban_beta_1(ntype)
	  else if (n_total_urban(j,myyear).lt.20)then
	 	conv = conv0*urban_beta_2(ntype)
	  else
	 	conv = conv0*urban_beta_3(ntype)
	  end if
	  	
!	  r_so2 = edailyso2(1,j,myyear)*alpha_so2(j)
	  r_so2 = edailyusox(1,j,myyear)
     &		* conv
     
!     	  r_co  = edailyco (1,j,myyear)*alpha_co (j,myyear)
	  r_co  = edailyuco (1,j,myyear)
     &		* conv
     
!     	  r_nox = edailynox(1,j,myyear)*alpha_nox(j,myyear)
	  r_nox = edailyunox(1,j,myyear)
     &		* conv
     
!     	  r_voc = 0.0
	  r_voc = edailyunmv(1,j,myyear)
     &		* conv

	  Temp = T(1,j,1)*airpress(1)**0.286

	  ! --- jday =julian day [1,365]
	  xxx = float(jday)
	  		       
	  call metamodel( Temp, j_date, 
     &			  rlati(j),rtemp(j),rclou(j),
     &			  rmixi(j),rtime(j),
     & 			  r_so2,   r_co,   r_nox,  r_voc,     
     &			  raqino (ntype,j),raqinoo(ntype,j),
     &    		  raqivoc(ntype,j),raqiozo(ntype,j), 
     &			  raqisox(ntype,j),tmp)

!
!	if(myyear.eq.1.and.mymonth.eq.8)then
!	  print *,"ntype = ",ntype, "j = ",j
!	  print *,"AQINO = ",raqino(ntype,j),"AQINO2 = ",raqinoo(ntype,j)
!	  print *,"AQIVOC = ",raqivoc(ntype,j),"AQIOZO = ",raqiozo(ntype,j)
!	  print *,"AQISOX = ",raqisox(ntype,j)
!	endif
!
	
	  ! --- convert mass and flux from  kg/km^2 to kg
	  do iii=1,25
	    results_meta(iii,ntype,j)= 
     &			tmp(iii)*n_urban(ntype,j,myyear)*4.e4
     &		       *( 1.0 + sin(3.1415926			 
     &		        *(xxx - xc_meta(iii))/w_meta(iii))
     &				             /c_meta(iii) )
     &			             
     	  end do
	  results_meta(26,ntype,j) = tmp(26)
     &		       *( 1.0 + sin(3.1415926
     &		        *(xxx - xc_meta(26))/w_meta(26))
     &				            /c_meta(26) )	  
	  results_meta(27,ntype,j) = tmp(27)
	  
	  do iii=28,meta_nvar
	    results_meta(iii,ntype,j)= 
     &			tmp(iii)*n_urban(ntype,j,myyear)*4.e4
     &		  *( 1.0 + sin(3.1415926			 
     &		   *(xxx - xc_meta(iii))/w_meta(iii))/c_meta(iii) )
     	  end do


	  ! --- all results should be positive definite
	  do iii=1,meta_nvar
	    if(results_meta(iii,ntype,j).le.0.0)
     &         results_meta(iii,ntype,j) = 0.0
	  end do
	
	  ! --- convert flux to kg (mass is already in kg)
	  !  Note: CO and NO fluxes have been included 
	  !	in chememission.F, they should not be recounted here
	  !  Note also: ymass mixing with the grid model is no longer
	  !	needed so that both meta and grid calculates its own
	  !	concentration, only connection is the yflux
	  !	This has been tested via eppa02 May 2004. 
	  !	Chien Wang 062304
	  !
!	  ymass_no   = ymass_no  
!     &		     + results_meta(1,ntype,j)
!	  ymass_no2  = ymass_no2
!     &		     + results_meta(3,ntype,j)
     	  yflux_no2  = yflux_no2
     &		     + results_meta(4,ntype,j)
!	  ymass_o3   = ymass_o3
!     &		     + results_meta(5,ntype,j)
	  yflux_o3   = yflux_o3
     &		     + results_meta(6,ntype,j)
	  yflux_hno3 = yflux_hno3
     &		     + results_meta(8,ntype,j)
	  yflux_n2o5 = yflux_n2o5
     &		     + results_meta(9,ntype,j)
!     	  ymass_co   = ymass_co
!     &		     + results_meta(10,ntype,j)
!	  ymass_ch2o = ymass_ch2o
!     &		     + results_meta(12,ntype,j)
	  yflux_ch2o = yflux_ch2o
     &		     + results_meta(13,ntype,j)
	 end if
	 
20	continue

	  ! ---
	  ! --- incorperating meta results into mixing ratios:
	  ! ---
	  ! --- NO
	  do k=1,ktop
		x11(k) = xno (1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_no,yflux_no,x11)
	  do k=1,ktop
		xno(1,j,k) = x11(k)
	  enddo
		  
	  ! --- NO2
	  do k=1,ktop
		x11(k) = xno2(1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_no2,yflux_no2,x11)
	  do k=1,ktop
		xno2(1,j,k) = x11(k)
	  enddo
		  
	  ! --- O3
	  do k=1,ktop
		x11(k) = o3 (1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_o3,yflux_o3,x11)
	  do k=1,ktop
		o3 (1,j,k) = x11(k)
	  enddo

	  ! --- HNO3
	  do k=1,ktop
		x11(k) = hno3 (1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_hno3,yflux_hno3,x11)
	  do k=1,ktop
		hno3(1,j,k) = x11(k)
	  enddo

	  ! --- N2O5
	  do k=1,ktop
		x11(k) = xn2o5(1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_n2o5,yflux_n2o5,x11)
	  do k=1,ktop
		xn2o5(1,j,k) = x11(k)
	  enddo

	  ! --- CO
	  do k=1,ktop
		x11(k) = co(1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_co,yflux_co,x11)
	  do k=1,ktop
		co(1,j,k) = x11(k)
	  enddo
	  
	  ! --- HCHO
	  do k=1,ktop
		x11(k) = ch2o(1,j,k)
	  enddo
	  call chemmeta_mass(1,j,ktop,tmass,
     &			     urban_area,rural_area,
     &		             ymass_ch2o,yflux_ch2o,x11)
	  do k=1,ktop
		ch2o(1,j,k) = x11(k)
	  enddo
		    	
	  ! ===
	  ! === convert mass into mole fraction in ppb
	  ! ===
	  conv0 = tmass/urban_area*28.97296245*1.e9
	  do ntype=1,3
	    if(n_urban(ntype,j,myyear).gt.0)then
	   	conv = conv0*float(n_total_urban(j,myyear))
     &		            /float(n_urban(ntype,j,myyear))
		! ---   1:	NO   mass
	  	results_meta(1,ntype,j)  = 
     &					 results_meta(1,ntype,j)
     &				         /awNO*conv
		! ---   3:	NO2  mass
		results_meta(3,ntype,j)  = 
     &					 results_meta(3,ntype,j)
     &				         /awNO2*conv
		! ---   5:	O3   mass
		results_meta(5,ntype,j)  = 
     &					 results_meta(5,ntype,j)
     &				    	 /awO3*conv
     		! ---  10:	CO mass
		results_meta(10,ntype,j) = 
     &					 results_meta(10,ntype,j)
     &				         /awCO*conv
		! ---  12:	HCHO mass
		results_meta(12,ntype,j) = 
     &					 results_meta(12,ntype,j)
     &				         /awCH2O*conv
		! ---  14:	ALD2 mass
		results_meta(14,ntype,j) = 
     &					 results_meta(14,ntype,j)
     &				         /awALD2*conv
		! ---  15:	MEK  mass
		results_meta(15,ntype,j) = 
     &					 results_meta(15,ntype,j)
     &				         /awMEK*conv
		! ---  17:	ALKA mass
		results_meta(17,ntype,j) = 
     &					 results_meta(17,ntype,j)
     &				         /awALKA*conv
		! ---  18:	ETHE mass
		results_meta(18,ntype,j) = 
     &					 results_meta(18,ntype,j)
     &				         /awETHE*conv
		! ---  19:	ALKE mass
		results_meta(19,ntype,j) = 
     &					 results_meta(19,ntype,j)
     &				         /awALKE*conv
		! ---  20:	TOLU mass
		results_meta(20,ntype,j) = 
     &					 results_meta(20,ntype,j)
     &				    	 /awTOLU*conv
		! ---  21:	AROM mass
		results_meta(21,ntype,j) = 
     &				         results_meta(21,ntype,j)
     &				         /awAROM*conv
		! ---  22:	SO2  mass
		results_meta(22,ntype,j) = 
     &					 results_meta(22,ntype,j)
     &				         /awSO2*conv
		! ---  24:	SO3  mass
		results_meta(24,ntype,j) = 
     &					 results_meta(24,ntype,j)
     &				         /awSO3*conv
     		! ---  28:	NO daily-mean mass
		results_meta(28,ntype,j) = 
     &					 results_meta(28,ntype,j)
     &				         /awNO*conv
		! ---  29:	NO2 daily-mean mass
		results_meta(29,ntype,j) = 
     &					 results_meta(29,ntype,j)
     &				         /awNO2*conv
		! ---  30:	O3 daily-mean mass
		results_meta(30,ntype,j) = 
     &					 results_meta(30,ntype,j)
     &				    	 /awO3*conv
		! ---  31:	HONO daily-mean mass
		results_meta(31,ntype,j) = 
     &				         results_meta(31,ntype,j)
     &				         /awHONO*conv
		! ---  32:	HNO3 daily-mean mass
		results_meta(32,ntype,j) = 
     &					 results_meta(32,ntype,j)
     &				         /awHNO3*conv
		! ---  33:	HNO4 daily-mean mass
		results_meta(33,ntype,j) = 
     &					 results_meta(33,ntype,j)
     &				         /awHNO4*conv
     		! ---  34:	N2O5 daily-mean mass
		results_meta(34,ntype,j) = 
     &					 results_meta(34,ntype,j)
     &				         /awN2O5*conv
		! ---  35:	NO3 daily-mean mass
		results_meta(35,ntype,j) = 
     &					 results_meta(35,ntype,j)
     &				         /awNO3*conv
		! ---  36:	HO2 daily-mean mass
		results_meta(36,ntype,j) = 
     &					 results_meta(36,ntype,j)
     &				    	 /awHO2*conv
		! ---  37:	CO daily-mean mass
		results_meta(37,ntype,j) = 
     &				         results_meta(37,ntype,j)
     &				         /awCO*conv
		! ---  38:	HCHO daily-mean mass
		results_meta(38,ntype,j) = 
     &					 results_meta(38,ntype,j)
     &				         /awCH2O*conv
		! ---  39:	ALD2 daily-mean  mass
		results_meta(39,ntype,j) = 
     &					 results_meta(39,ntype,j)
     &				         /awALD2*conv
     		! ---  40:	MEK daily-mean mass
		results_meta(40,ntype,j) = 
     &					 results_meta(40,ntype,j)
     &				         /awMEK*conv
		! ---  42:	PAN daily-mean mass
		results_meta(42,ntype,j) = 
     &					 results_meta(42,ntype,j)
     &				         /awPAN*conv
		! ---  46:	ALKA daily-mean mass
		results_meta(46,ntype,j) = 
     &					 results_meta(46,ntype,j)
     &				    	 /awALKA*conv
		! ---  47:	ETHE daily-mean mass
		results_meta(47,ntype,j) = 
     &				         results_meta(47,ntype,j)
     &				         /awETHE*conv
		! ---  48:	ALKE daily-mean mass
		results_meta(48,ntype,j) = 
     &					 results_meta(48,ntype,j)
     &				         /awALKE*conv
		! ---  49:	TOLU daily-mean mass
		results_meta(49,ntype,j) = 
     &					 results_meta(49,ntype,j)
     &				         /awTOLU*conv
     		! ---  50:	AROM daily-mean mass
		results_meta(50,ntype,j) = 
     &					 results_meta(50,ntype,j)
     &				         /awAROM*conv
		! ---  54:	H2O2 daily-mean mass
		results_meta(54,ntype,j) = 
     &					 results_meta(54,ntype,j)
     &				         /awH2O2*conv
		! ---  57:	SO2 daily-mean mass
		results_meta(57,ntype,j) = 
     &					 results_meta(57,ntype,j)
     &				    	 /awSO2*conv
		! ---  58:	SO3 daily-mean mass
		results_meta(58,ntype,j) = 
     &				         results_meta(58,ntype,j)
     &				         /awSO3*conv
	    else
	  	results_meta(1, ntype,j) = 0.0
		results_meta(3, ntype,j) = 0.0 
		results_meta(5, ntype,j) = 0.0
		results_meta(10,ntype,j) = 0.0
		results_meta(12,ntype,j) = 0.0
		results_meta(14,ntype,j) = 0.0
		results_meta(15,ntype,j) = 0.0
		results_meta(17,ntype,j) = 0.0
		results_meta(18,ntype,j) = 0.0
		results_meta(19,ntype,j) = 0.0
		results_meta(20,ntype,j) = 0.0
		results_meta(21,ntype,j) = 0.0
		results_meta(22,ntype,j) = 0.0
		results_meta(24,ntype,j) = 0.0
		results_meta(28,ntype,j) = 0.0
		results_meta(29,ntype,j) = 0.0
		results_meta(30,ntype,j) = 0.0
		results_meta(31,ntype,j) = 0.0
		results_meta(32,ntype,j) = 0.0
		results_meta(33,ntype,j) = 0.0
		results_meta(34,ntype,j) = 0.0
		results_meta(35,ntype,j) = 0.0
		results_meta(36,ntype,j) = 0.0
		results_meta(37,ntype,j) = 0.0
		results_meta(38,ntype,j) = 0.0
		results_meta(39,ntype,j) = 0.0
		results_meta(40,ntype,j) = 0.0
		results_meta(42,ntype,j) = 0.0
		results_meta(46,ntype,j) = 0.0
		results_meta(47,ntype,j) = 0.0
		results_meta(48,ntype,j) = 0.0
		results_meta(49,ntype,j) = 0.0
		results_meta(50,ntype,j) = 0.0
		results_meta(54,ntype,j) = 0.0
		results_meta(57,ntype,j) = 0.0
		results_meta(58,ntype,j) = 0.0
	    endif
     	  end do
	  
10	continue

	do j=1,nlat
	do ntype=1,3
	do i=1,meta_nvar
	  results_mon(i,ntype,j) = results_mon (i,ntype,j)
     &				 + results_meta(i,ntype,j)
	end do
	end do
	end do
	
	nstep_meta = nstep_meta + 1

#endif
			
	return
	 end


!	===============================================
	subroutine chemmeta_mass(i,j,ktop,tmass1,
     &			         urban_area,rural_area,
     &			         ymass, yflux, x11)
!	===============================================

! --------------------------------------------------------
! A subroutine for recalculating zonal mean mixing ratrios
!	by incorperating GACM and META results
! --------------------------------------------------------

#include "chem_para"
#include "chem_com"

#include "chem_meta"
	
#include "BD2G04.COM"
	
	dimension x11(ktop)
	real urban_area, rural_area, xmass, xmix

#if ( defined CPL_CHEM ) && ( defined CPL_META )

        xmass    = 0.0		!rural tracer mass in 10^-9 kg
	do k=1,ktop
	  xmass = xmass
     &          + airmass(i,j,k)*x11(k)
	end do
	
	! === PM scheme
!        xmix = ( rural_area*xmass + 
!!     &		(urban_area*ymass +  yflux)*1.e9 )
!     &		(ymass +  yflux)*1.e9 )
!     &       * tmass1     				!ppbm
	! === FM scheme
        xmix = ( xmass + 
     &		(ymass +  yflux)*1.e9 )
     &       * tmass1     				!ppbm

	do k=1,ktop
		x11(k)  = xmix
	end do

#endif

	return
	 end

!	====================	 
	Block Data Meta_data
!	====================
	
#include "chem_para"
#include "chem_meta"

#if ( defined CPL_CHEM ) && ( defined CPL_META )

#if ( N_LAT == 24 )
	data rlati/90.0,82.2,74.3,66.5,58.7,50.9,
     &		   43.0,35.2,27.4,19.6,11.7, 3.9,
     &		    3.9,11.7,19.6,27.4,35.2,43.0,
     &		   50.9,58.7,66.5,74.3,82.2,90.0/


	data results_mon/4176*0.0/
	
     	data raqino /72*0.1/
	data raqinoo/72*0.1/
	data raqivoc/72*0.1/
	data raqiozo/72*0.1/
	data raqisox/72*0.1/	
#endif

#if ( N_LAT == 46 )
	data rlati
     &		 /
     &		 90.0,  86.0,  82.0,  78.0,  74.0, 
     &		 70.0,  66.0,  62.0,  58.0,  54.0, 
     &		 50.0,  46.0,  42.0,  38.0,  34.0, 
     &		 30.0,  26.0,  22.0,  18.0,  14.0, 
     &		 10.0,   6.0,   2.0,   2.0,   6.0,  
     &		 10.0,  14.0,  18.0,  22.0,  26.0,  
     &		 30.0,  34.0,  38.0,  42.0,  46.0,  
     &		 50.0,  54.0,  58.0,  62.0,  66.0,  
     &		 70.0,  74.0,  78.0,  82.0,  86.0,  
     &		 90.0
     &		/ 

	data results_mon/8004*0.0/
	
     	data raqino /138*0.1/
	data raqinoo/138*0.1/
	data raqivoc/138*0.1/
	data raqiozo/138*0.1/
	data raqisox/138*0.1/	
#endif

	data nstep_meta/0/
	     
	data xc_meta/  -86.4, -88.4,-117.6,
     &		      -111.8,  71.2,  81.1,
     &		       198.9,  82.7,-101.0,  
     &		        86.6,  89.7,  95.7,
     &		       104.6,-126.0,  47.9, 
     &		       -74.5,-114.3,-124.8,
     &		      -106.5,-121.6,-115.6,
     &		      -108.7,-102.3,  73.8,
     &		        85.7,  76.7,   0.0,
     &		       -84.4,-112.5,  81.9,
     &		       6*0.0,
     &			97.9,  99.4,20*0.0/    
     
     	data w_meta /183.6,184.3,184.8,
     &		     184.2,183.9,184.0,
     &		     183.6,183.1,184.9,
     &		     185.7,181.0,182.9,
     &		     183.9,184.1,189.3,
     &		     183.8,184.1,184.9,
     &		     183.0,184.6,182.9,
     &		     183.1,184.6,184.0,
     &		     182.1,184.2,  0.0,
     &		     183.3,183.4,183.3,
     &		    6* 0.0,
     &		     184.7,183.1,20*0.0/
     
     	data a_meta / 0.015,0.153,0.288,
     &		      0.466,8.040,9.880,
     &		      0.007,0.797,0.014,
     &		      0.886,1.016,0.154,
     &		      0.235,0.013,0.044,
     &		      0.323,0.188,0.007,
     &		      0.004,0.028,0.009,
     &		      0.042,0.052,0.053,
     &		      0.080,0.006,  0.0,
     &		      0.141,0.242,6.485,
     &		    6*0.0,
     &		      0.496,0.184,20*0.0/

	data c_meta / 11.4, -0.8,  4.0,
     &		       4.5,  3.6,  3.9,
     &		      10.5,  4.2,  1.1,
     &		      79.0, 87.0, 11.2,
     &		       7.6,  6.5, 10.3,
     &		       2.6,  5.3,  6.7,
     &		      12.4,  5.7, 12.8,
     &		      15.5, 10.1,  2.8,
     &		       2.0,  3.1,  0.0,
     &		      10.2,  8.8,  4.1,
     &		     6*0.0,
     &		     227.5,  7.6,20*0.0/
	
       	data urban_beta_1/1.0, 0.0,  0.0/
	data urban_beta_2/0.9, 1.3,  0.0/
	data urban_beta_3/0.8, 1.25, 3.0/

#endif
	   
	end
	


