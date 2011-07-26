
! ==================================================
!
! ctrparam.h
! ----------
!
!  Purpose:	A header file contains cpp control 
!			parameters for the model
!
!  -------------------------------------------------
!	
!  Usage:	1. (un)comment #define line; or 
!		2. #define/#undef x [number]
!		  to set given cpp parameters x
!		  to be true or false  
!
!  -------------------------------------------------
!
!  Author:	Chien Wang
!
!  Revision:
!	Date	By		Brief Description
!	----	--		-----------------
!	052400	Chien Wang	created
!	121800	Chien Wang	add NYR_CHEM
!	062304	Chien Wang	rev for igsm2
!
! ==================================================

!
! === number of latitudinal grid points
!
#define N_LAT 46

!
! === number of vertical pressure layers
!
#define N_LEV 11

#define N_LON0 72

! === number of vertical  layers in duffusive ocean model
!
! --- Current options include 12 and 11
!
#define N_LEVOCEAN 11

!
! === defines 3D ocean
!
!#define OCEAN_3D 1
#define ML_2D 1

!
! === coupling ocean co2 model
!
!
!#define CPL_OCEANCO2 1

! === coupling CLM
!
!#define CLM 1
!
!
! === coupling chemistry
!
 
!#define CPL_TEM 1

!
! === coupling natural emission model
!
!#define CPL_NEM 1
!
! === coupling chemistry
!
! --- IMPORTANT NOTE: In the current model version
!       if you choose to undefine CPL_CHEM, you have
!       to also close CPL_NEM, CPL_META, and CPL_OCEANCO2.
!       In addition, you should manage to set data reading
!       correctly if you want use PREDICTED_GASES and
!       PREDICTED_AEROSOLS.
!
!#define CPL_CHEM 1

!=== number of year of chemstry model integration
!
! --- 124 for period from 1977 to 2100
!
#define NYR_CHEM 1
!
! === coupling parameterized urban airshed model ( meta )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!               might not work properly
!
!#define CPL_META 1

!
! === coupling ocean co2 model
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!               might not work properly
!
!
! === use predicted concentrations of chemical species
!       by chemistry model to calculate radiative forcing
!       ( ifghgpredict )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!               might not work properly
!
!#define PREDICTED_GASES 1

!
! === use predicted aerosol concentrations by chemistry model
!       to calculate radiative forcing
!       ( ifaerpredict )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!               might not work properly
!
!#define PREDICTED_AEROSOL 1
!#define PREDICTED_BC 1

!
! === including "3 gases"
!
!#define INC_3GASES 1

!
! === use predicted ozone in radiation module
!
!#define O3_RAD 1

!
! === use forced ozone scenario ( ifo3forced )
!
! --- NOTE: currently CPL_CHEM true setting excludes
!		this parameter
!
!#define O3_FORCED 1

!
! === modify albedo to simulate sulfate aerosol forcing
!	( ifsulfalb )
!
! --- NOTE: currently CPL_CHEM true setting excludes
!		this parameter
!
!#define SVI_ALBEDO 1
!
!#define VOL_AER 1
