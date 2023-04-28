
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
! --- Current options include 24 and 46
!	DON'T DO ANYTHING SxxxxD!
!
!#define N_LAT 24
#define N_LAT 46

!
! === number of vertical pressure layers
!
! --- Current options include 9 and 11
!	DON'T DO ANYTHING SxxxxD!
!
!#define N_LEV 9
#define N_LEV 11

!#define N_LON0 36
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

!
! defines 2D ml ocean
!
#define ML_2D 1

!
! === coupling ocean co2 model
!
#define CPL_OCEANCO2 1

! === CLM
!
#define CLM 1
!
! === using CLM35
!
!#define CLM35 1
!
!

! === coupling TEM
!
#define CPL_TEM 1

! === data for offline run of TEM
!
!#define DATA4TEM 1
!#define DATA4CLM 1

!
! === coupling natural emission model
!
#define CPL_NEM 1

!
! === coupling chemistry
!
! --- IMPORTANT NOTE: In the current model version
!	if you choose to undefine CPL_CHEM, you have
!	to also close  CPL_META
!	In addition, you should manage to set data reading
!	correctly if you want use PREDICTED_GASES and
!	PREDICTED_AEROSOLS.
!
#define CPL_CHEM 1
!
!
!#define OLD_CHEM 1
!
!#define ACCRI 1
!
!
!
! === number of year of chemstry model integration
!
! --- 124 for period from 1977 to 2100
! --- 424 for period from 1977 to 2400
! --- 114 historical emissions for period from 1892 to 2005
!
!#define NYR_CHEM 124
#define NYR_CHEM 1

!
! === coupling parameterized urban airshed model ( meta )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!		might not work properly
!
#define CPL_META 1

!
! === use predicted concentrations of chemical species
!	by chemistry model to calculate radiative forcing
!	( ifghgpredict )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!		might not work properly
!
#define PREDICTED_GASES 1

!
! === use predicted aerosol concentrations by chemistry model
!	to calculate radiative forcing
!	( ifaerpredict )
!
! --- NOTE: with false setting of CPL_CHEM this parameter
!		might not work properly
!
#define PREDICTED_AEROSOL 1
#define PREDICTED_BC 1

!
! === including "3 gases"
!
#define INC_3GASES 1

!
! === use predicted ozone in radiation module
!
#define O3_RAD 1

!
#define VOL_AER 1
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
!#define ORBITAL_FOR 1
!#define TRACERS 1
#define IPCC_FORCING 1
