#include <misc.h>
#include <preproc.h>

module clmtype

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Define derived type hierarchy. Includes declaration of
! the clm derived type and 1d mapping arrays. 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar
!
! !PUBLIC TYPES:
  implicit none
!                              
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!*******************************************************************************
!----------------------------------------------------
! Begin definition of conservation check structures
!----------------------------------------------------
! energy balance structure
!----------------------------------------------------
type energy_balance_type
   real(r8):: errsoi        !soil/lake energy conservation error (W/m**2)
   real(r8):: errseb        !surface energy conservation error (W/m**2)
   real(r8):: errsol        !solar radiation conservation error (W/m**2)
   real(r8):: errlon        !longwave radiation conservation error (W/m**2)
   real(r8):: acc_errsoi    !accumulation of energy balance error
   real(r8):: acc_errseb    !accumulation of surface energy balance error
   real(r8):: acc_errsol    !accumulation of energy balance error
end type energy_balance_type

!----------------------------------------------------
! water balance structure
!----------------------------------------------------
type water_balance_type
   real(r8):: begwb         !water mass begining of the time step
   real(r8):: endwb         !water mass end of the time step
   real(r8):: errh2o        !water conservation error (mm H2O)
   real(r8):: acc_errh2o    !accumulation of water balance error
end type water_balance_type

!----------------------------------------------------
! carbon balance structure
!----------------------------------------------------
type carbon_balance_type
   real(r8):: dummy_entry
end type carbon_balance_type

!----------------------------------------------------
! nitrogen balance structure
!----------------------------------------------------
type nitrogen_balance_type
   real(r8):: dummy_entry
end type nitrogen_balance_type
!----------------------------------------------------
! End definition of conservation check structures
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the pft_type level
!----------------------------------------------------
! pft physical state variables structure
!----------------------------------------------------
type pft_pstate_type
   type(column_pstate_type), pointer:: cps !pointer to higher level ps struct 
   integer :: itype 		!vegetation type for this pft
   integer :: itypveg 	        !vegetation type for this pft (currently redundant to itype)
   integer :: frac_veg_nosno 	!fraction of vegetation not covered by snow (0 OR 1) [-] 
   integer :: frac_veg_nosno_alb!fraction of vegetation not covered by snow (0 OR 1) [-] 
   real(r8):: rsw 		!soil water content for root zone
   real(r8):: emv 		!vegetation emissivity
   real(r8):: z0mv 		!roughness length over vegetation, momentum [m]
   real(r8):: z0hv 		!roughness length over vegetation, sensible heat [m]
   real(r8):: z0qv 		!roughness length over vegetation, latent heat [m]
   real(r8):: rootfr(nlevsoi)	!fraction of roots in each soil layer
   real(r8):: rootr (nlevsoi)	!effective fraction of roots in each soil layer
   real(r8):: dewmx 		!Maximum allowed dew [mm]
   real(r8):: rssun 		!sunlit stomatal resistance (s/m)
   real(r8):: rssha 		!shaded stomatal resistance (s/m)
   real(r8):: laisun 		!sunlit leaf area
   real(r8):: laisha 		!shaded leaf area
   real(r8):: btran 		!transpiration wetness factor (0 to 1)
   real(r8):: fsun 		!sunlit fraction of canopy
   real(r8):: tlai 		!one-sided leaf area index, no burying by snow
   real(r8):: tsai 		!one-sided stem area index, no burying by snow
#if (defined PCP2PFT)
   real(r8):: pcp2pft           !Precipitation distribution factor
#endif
   real(r8):: elai 		!one-sided leaf area index with burying by snow
   real(r8):: esai 		!one-sided stem area index with burying by snow
   real(r8):: igs 		!growing season index (0=off, 1=on)
   real(r8):: stembio		!stem biomass (kg /m**2)
   real(r8):: rootbio		!root biomass (kg /m**2)
   real(r8):: fwet 		!fraction of canopy that is wet (0 to 1)
   real(r8):: fdry 		!fraction of foliage that is green and dry [-] (new)
   real(r8):: dt_veg 		!change in t_veg, last iteration (Kelvin)
   real(r8):: htop 		!canopy top (m)
   real(r8):: hbot 		!canopy bottom (m)
   real(r8):: z0m 		!momentum roughness length (m)
   real(r8):: displa 		!displacement height (m)
   real(r8):: albd(numrad)	!surface albedo (direct)
   real(r8):: albi(numrad)	!surface albedo (diffuse)
   real(r8):: fabd(numrad)	!flux absorbed by veg per unit direct flux
   real(r8):: fabi(numrad)	!flux absorbed by veg per unit diffuse flux
   real(r8):: ftdd(numrad)	!down direct flux below veg per unit dir flx
   real(r8):: ftid(numrad)	!down diffuse flux below veg per unit dir flx
   real(r8):: ftii(numrad)	!down diffuse flux below veg per unit dif flx
   real(r8):: ndvi 		!Normalized Difference Vegetation Index
   real(r8):: u10 		!10-m wind (m/s) (for dust model)
   real(r8):: ram1 		!aerodynamical resistance (s/m)
   real(r8):: fv 		!friction velocity (m/s) (for dust model)
   real(r8):: area 		!total land area for this pft (km^2)
   real(r8):: wt 		!weight (relative to column) for this pft (0-1)
   real(r8):: wtxy 		!weight (relative to gridcell) for this pft (0-1)
   integer :: ixy               !xy lon index
   integer :: jxy               !xy lat index
   integer :: mxy               !m index for laixy(i,j,m),etc.
   integer :: index1d           !index into 1d global pft array  
end type pft_pstate_type

!----------------------------------------------------
! pft ecophysiological constants structure
!----------------------------------------------------
type pft_epc_type
   integer :: ncorn 		!value for corn
   integer :: nwheat 		!value for wheat
   integer :: noveg 	        !value for not vegetated
   integer :: ntree 		!value for last type of tree
   real(r8):: smpmax            !Wilting point potential in mm
   real(r8):: foln 	        !foliage nitrogen (%)
   real(r8):: dleaf 		!characteristic leaf dimension (m)
   real(r8):: c3psn 		!photosynthetic pathway: 0. = c4, 1. = c3
   real(r8):: vcmx25 		!max rate of carboxylation at 25C (umol CO2/m**2/s)
   real(r8):: mp 		!slope of conductance-to-photosynthesis relationship
   real(r8):: qe25 		!quantum efficiency at 25C (umol CO2 / umol photon)
   real(r8):: xl 		!leaf/stem orientation index
   real(r8):: rhol(numrad)      !leaf reflectance: 1=vis, 2=nir
   real(r8):: rhos(numrad)      !stem reflectance: 1=vis, 2=nir
   real(r8):: taul(numrad)      !leaf transmittance: 1=vis, 2=nir
   real(r8):: taus(numrad)      !stem transmittance: 1=vis, 2=nir 
   real(r8):: z0mr 		!ratio of momentum roughness length to canopy top height (-)
   real(r8):: displar 		!ratio of displacement height to canopy top height (-)
   real(r8):: roota_par 	!CLM rooting distribution parameter [1/m]
   real(r8):: rootb_par	 	!CLM rooting distribution parameter [1/m]
   real(r8):: sla 		!sp. leaf area [m2 leaf g-1 carbon]
end type pft_epc_type

!----------------------------------------------------
! pft DGVM-specific ecophysiological constants structure
!----------------------------------------------------
type pft_dgvepc_type
   real(r8):: respcoeff       !maintenance respiration coefficient [-]
   real(r8):: flam            !flammability threshold [units?]
   real(r8):: resist          !fire resistance index [units?]
   real(r8):: l_turn          !leaf turnover period [years]
   real(r8):: l_long          !leaf longevity [years]
   real(r8):: s_turn          !sapwood turnover period [years]
   real(r8):: r_turn          !root turnover period [years]
   real(r8):: l_cton          !leaf C:N (mass ratio)
   real(r8):: s_cton          !sapwood C:N (mass ratio)
   real(r8):: r_cton          !root C:N (mass ratio)
   real(r8):: l_morph         !leaf morphology: 1=broad, 2=needle, 3=grass
   real(r8):: l_phen          !leaf phenology: 1=everg, 2=summerg, 3=raing, 4=any
   real(r8):: lmtorm          !leaf:root ratio under non-water stressed conditions
   real(r8):: crownarea_max   !tree maximum crown area [m2]
   real(r8):: init_lai        !sapling (or initial grass) LAI [-]
   real(r8):: x               !sapling: (heart+sapwood)/sapwood [-]
   real(r8):: tcmin           !minimum coldest monthly mean temperature [units?]
   real(r8):: tcmax           !maximum coldest monthly mean temperature [units?]
   real(r8):: gddmin          !minimum growing degree days (at or above 5 C)
   real(r8):: twmax           !upper limit of temperature of the warmest month [units?]
   real(r8):: lm_sapl 
   real(r8):: sm_sapl 
   real(r8):: hm_sapl 
   real(r8):: rm_sapl 
   logical :: tree
   logical :: summergreen
   logical :: raingreen
   real(r8):: reinickerp      !parameter in allometric equation
   real(r8):: wooddens	      !wood density (gC/m3)
   real(r8):: latosa	      !ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b)
   real(r8):: allom1	      !parameter in allometric
   real(r8):: allom2          !parameter in allometric				
   real(r8):: allom3	      !parameter in allometric
end type pft_dgvepc_type

!----------------------------------------------------
! pft energy state variables structure
!----------------------------------------------------
type pft_estate_type
   real(r8):: t_ref2m          !2 m height surface air temperature (Kelvin)
   real(r8):: t_af             !canopy air temperature (Kelvin)
   real(r8):: t_veg            !vegetation temperature (Kelvin)
   real(r8):: t_rad_pft        !radiative temperature (Kelvin)
end type pft_estate_type

!----------------------------------------------------
! pft water state variables structure
!----------------------------------------------------
type pft_wstate_type
   real(r8)::	h2ocan         !canopy water (mm H2O)
end type pft_wstate_type

!----------------------------------------------------
! pft carbon state variables structure
!----------------------------------------------------
type pft_cstate_type
   real(r8)::	dummy_entry
end type pft_cstate_type

!----------------------------------------------------
! pft nitrogen state variables structure
!----------------------------------------------------
type pft_nstate_type
   real(r8):: dummy_entry
end type pft_nstate_type

!----------------------------------------------------
! pft VOC state variables structure
!----------------------------------------------------
type pft_vstate_type
   real(r8):: dummy_entry
end type pft_vstate_type

!----------------------------------------------------
! pft DGVM state variables structure
!----------------------------------------------------
type pft_dgvstate_type
   real(r8):: agdd0               !accumulated growing degree days above 0 deg C
   real(r8):: agdd5               !accumulated growing degree days above -5
   real(r8):: agddtw              !accumulated growing degree days above twmax
   real(r8):: agdd                !accumulated growing degree days above 5
   real(r8):: t10                 !10-day running mean of the 2 m temperature (K)
   real(r8):: t_mo                !30-day average temperature (Kelvin)
   real(r8):: t_mo_min            !annual min of t_mo (Kelvin)
   real(r8):: fnpsn10             !10-day running mean net photosynthesis
   real(r8):: prec365             !365-day running mean of tot. precipitation
   real(r8):: agdd20              !20-yr running mean of agdd
   real(r8):: tmomin20            !20-yr running mean of tmomin
   real(r8):: t10min              !annual minimum of 10-day running mean (K)
   real(r8):: tsoi25              !soil temperature to 0.25 m (Kelvin)
   real(r8):: annpsn              !annual photosynthesis (umol CO2 /m**2)
   real(r8):: annpsnpot           !annual potential photosynthesis (same units)
   logical :: present             !whether PFT present in patch
   real(r8):: dphen               !phenology [0 to 1]
   real(r8):: leafon              !leafon days
   real(r8):: leafof              !leafoff days
   real(r8):: nind                !number of individuals (#/m**2)
   real(r8):: lm_ind              !individual leaf mass
   real(r8):: sm_ind              !individual sapwood mass
   real(r8):: hm_ind              !individual heartwood mass
   real(r8):: rm_ind              !individual root mass
   real(r8):: lai_ind             !LAI per individual
   real(r8):: fpcinc              !foliar projective cover increment (fraction) 
   real(r8):: fpcgrid             !foliar projective cover on gridcell (fraction)
   real(r8):: crownarea           !area that each individual tree takes up (m^2)
   real(r8):: bm_inc              !biomass increment
   real(r8):: afmicr              !annual microbial respiration
   real(r8):: firelength          !fire season in days
   real(r8):: litterag            !above ground litter
   real(r8):: litterbg            !below ground litter
   real(r8):: cpool_fast          !fast carbon pool 
   real(r8):: cpool_slow          !slow carbon pool
   real(r8):: k_fast_ave          !decomposition rate
   real(r8):: k_slow_ave          !decomposition rate 
   real(r8):: litter_decom_ave    !decomposition rate
   real(r8):: turnover_ind        !???
end type pft_dgvstate_type

!----------------------------------------------------
! pft energy flux variables structure
!----------------------------------------------------
type pft_eflux_type
   real(r8):: sabg 		!solar radiation absorbed by ground (W/m**2)
   real(r8):: sabv 		!solar radiation absorbed by vegetation (W/m**2)
   real(r8):: fsa 		!solar radiation absorbed (total) (W/m**2)
   real(r8):: fsr 		!solar radiation reflected (W/m**2)
   real(r8):: parsun            !average absorbed PAR for sunlit leaves (W/m**2)
   real(r8):: parsha            !average absorbed PAR for shaded leaves (W/m**2)
   real(r8):: dlrad             !downward longwave radiation below the canopy [W/m2]
   real(r8):: ulrad             !upward longwave radiation above the canopy [W/m2]
   real(r8):: eflx_lh_tot 	!total latent heat flux (W/m8*2)  [+ to atm]
   real(r8):: eflx_lh_grnd 	!ground evaporation heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_soil_grnd 	!soil heat flux (W/m**2) [+ = into soil]
   real(r8):: eflx_sh_tot       !total sensible heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_sh_grnd      !sensible heat flux from ground (W/m**2) [+ to atm]
   real(r8):: eflx_sh_veg       !sensible heat flux from leaves (W/m**2) [+ to atm]
   real(r8):: eflx_lh_vege      !veg evaporation heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_lh_vegt      !veg transpiration heat flux (W/m**2) [+ to atm]
   real(r8):: cgrnd             !deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8):: cgrndl            !deriv of soil latent heat flux wrt soil temp [w/m**2/k]
   real(r8):: cgrnds            !deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
   real(r8):: eflx_gnet         !net heat flux into ground (W/m**2)
   real(r8):: dgnetdT           !derivative of net ground heat flux wrt soil temp (W/m**2 K)
   real(r8):: eflx_lwrad_out 	!emitted infrared (longwave) radiation (W/m**2)
   real(r8):: eflx_lwrad_net 	!net infrared (longwave) rad (W/m**2) [+ = to atm]
end type pft_eflux_type

!----------------------------------------------------
! pft momentum flux variables structure
!----------------------------------------------------
type pft_mflux_type
   real(r8)::  taux             !wind (shear) stress: e-w (kg/m/s**2)
   real(r8)::  tauy             !wind (shear) stress: n-s (kg/m/s**2)
end type pft_mflux_type

!----------------------------------------------------
! pft water flux variables structure
!----------------------------------------------------
type pft_wflux_type
   real(r8):: qflx_prec_intr 	!interception of precipitation [mm/s]
   real(r8):: qflx_prec_grnd 	!water onto ground including canopy runoff [kg/(m2 s)]
   real(r8):: qflx_rain_grnd 	!rain on ground after interception (mm H2O/s) [+]
   real(r8):: qflx_snow_grnd 	!snow on ground after interception (mm H2O/s) [+]
   real(r8):: qflx_snowcap 	!excess precipitation due to snow capping (mm H2O /s) [+]
   real(r8):: qflx_efpot        !potential evaporation
   real(r8):: qflx_evap_veg 	!vegetation evaporation (mm H2O/s) (+ = to atm)
   real(r8):: qflx_tran_veg 	!vegetation transpiration (mm H2O/s) (+ = to atm)
   real(r8):: qflx_evap_can     !evaporation from leaves and stems 
   real(r8):: qflx_evap_soi 	!soil evaporation (mm H2O/s) (+ = to atm)
   real(r8):: qflx_evap_tot 	!qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
   real(r8):: qflx_evap_grnd 	!ground surface evaporation rate (mm H2O/s) [+]
   real(r8):: qflx_dew_grnd 	!ground surface dew formation (mm H2O /s) [+]
   real(r8):: qflx_sub_snow 	!sublimation rate from snow pack (mm H2O /s) [+]
   real(r8):: qflx_dew_snow 	!surface dew added to snow pack (mm H2O /s) [+]
end type pft_wflux_type

!----------------------------------------------------
! pft carbon flux variables structure
!----------------------------------------------------
type pft_cflux_type
   real(r8)::	psnsun 		!sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8)::	psnsha 		!shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8)::	fpsn 		!photosynthesis (umol CO2 /m**2 /s)
   real(r8)::	frm 		!total maintenance respiration (umol CO2 /m**2/s)
   real(r8)::	frmf 		!leaf maintenance respiration  (umol CO2 /m**2 /s)
   real(r8)::	frms 		!stem maintenance respiration  (umol CO2 /m**2 /s)
   real(r8)::	frmr 		!root maintenance respiration  (umol CO2 /m**2 /s)
   real(r8)::	frg 		!growth respiration (umol CO2 /m**2 /s)
   real(r8)::	dmi 		!total dry matter production (ug /m**2 /s)
   real(r8)::	fco2 		!net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
   real(r8)::	fmicr 		!microbial respiration (umol CO2 /m**2 /s)
end type pft_cflux_type

!----------------------------------------------------
! pft nitrogen flux variables structure
!----------------------------------------------------
type pft_nflux_type
   real(r8):: dummy_entry
end type pft_nflux_type

!----------------------------------------------------
! pft VOC flux variables structure
!----------------------------------------------------
type pft_vflux_type
   real(r8):: vocflx_tot    	!total VOC flux into atmosphere [ug C m-2 h-1]
   real(r8):: vocflx(nvoc) 	!VOC flux [ug C m-2 h-1]
   real(r8):: vocflx_1       	!vocflx(1) (for history output) [ug C m-2 h-1]
   real(r8):: vocflx_2       	!vocflx(2) (for history output) [ug C m-2 h-1]
   real(r8):: vocflx_3       	!vocflx(3) (for history output) [ug C m-2 h-1]
   real(r8):: vocflx_4       	!vocflx(4) (for history output) [ug C m-2 h-1]
   real(r8):: vocflx_5       	!vocflx(5) (for history output) [ug C m-2 h-1]
end type pft_vflux_type

!----------------------------------------------------
! pft dust flux variables structure
!----------------------------------------------------
type pft_dflux_type	
   real(r8):: flx_mss_vrt_dst(ndst)  !surface dust emission (kg/m**2/s) [ + = to atm]
   real(r8):: flx_mss_vrt_dst_tot    !total dust flux into atmosphere
   real(r8):: vlc_trb(ndst)	     !turbulent deposition velocity (m/s)
   real(r8):: vlc_trb_1 	     !turbulent deposition velocity 1(m/s)
   real(r8):: vlc_trb_2 	     !turbulent deposition velocity 2(m/s)
   real(r8):: vlc_trb_3 	     !turbulent deposition velocity 3(m/s)
   real(r8):: vlc_trb_4 	     !turbulent deposition velocity 4(m/s)
end type pft_dflux_type

!----------------------------------------------------
! End definition of structures defined at the pft_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the column_type level
!----------------------------------------------------
! column physical state variables structure
!----------------------------------------------------
type column_pstate_type
   type(landunit_pstate_type), pointer:: lps	!pointer to higher level ps struct
   type(pft_pstate_type):: pps_a     !pft-level pstate variables averaged to the column
   integer :: itype 		     !type for this column
   integer :: npfts 		     !number of pfts for this column
   integer :: snl  		     !number of snow layers
   real(r8):: gwc_thr 		     !threshold soil moisture based on clay content
   real(r8):: mss_frc_cly_vld 	     ![frc] Mass fraction clay limited to 0.20
   real(r8):: mbl_bsn_fct 	     !??
   logical :: do_capsnow             !true => do snow capping
   real(r8):: snowdp 		     !snow height (m)
   real(r8):: snowage		     !non dimensional snow age [-] (new)
   real(r8):: frac_sno 		     !fraction of ground covered by snow (0 to 1)
   real(r8):: zi(-nlevsno+0:nlevsoi) !interface level below a "z" level (m)
   real(r8):: dz(-nlevsno+1:nlevsoi) !layer thickness (m)
   real(r8):: z (-nlevsno+1:nlevsoi) !layer depth (m)
   real(r8):: frac_iceold(-nlevsno+1:nlevsoi)  !fraction of ice relative to the tot water (new)
   integer :: imelt(-nlevsno+1:nlevsoi)        !flag for melting (=1), freezing (=2), Not=0 (new)
   real(r8):: eff_porosity(nlevsoi)            !effective porosity = porosity - vol_ice
   real(r8):: sfact 		     !term for implicit correction to evaporation
   real(r8):: sfactmax 		     !maximim of "sfact"
   real(r8):: emg 		     !ground emissivity
   real(r8):: z0mg 		     !roughness length over ground, momentum [m]
   real(r8):: z0hg 		     !roughness length over ground, sensible heat [m]
   real(r8):: z0qg 		     !roughness length over ground, latent heat [m]
   real(r8):: htvp 		     !latent heat of vapor of water (or sublimation) [j/kg]
   real(r8):: beta 		     !coefficient of convective velocity [-]
   real(r8):: zii 		     !convective boundary height [m]
   real(r8):: ur 		     !wind speed at reference height [m/s]
   real(r8):: albgrd(numrad)         !ground albedo (direct)
   real(r8):: albgri(numrad)         !ground albedo (diffuse)
   real(r8):: rootr_column(nlevsoi)  !effective fraction of roots in each soil layer
   real(r8):: wf                     !soil water as frac. of whc for top 0.5 m
   real(r8):: area 		     !total land area for this column (km^2)
   real(r8):: wt 		     !weight (relative to landunit) for this column (0-1)
   real(r8):: wtxy 	    	     !weight (relative to gridcell) for this column (0-1)
   integer :: ixy                    !xy lon index
   integer :: jxy                    !xy lat index
   integer :: index1d                !index into 1d global column array  
end type column_pstate_type

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
type column_estate_type
   type(pft_estate_type):: pes_a	   !pft-level energy state variables averaged to the column
   real(r8):: t_grnd  			   !ground temperature (Kelvin)
   real(r8):: t_af  			   !canopy air temperature (Kelvin)
   real(r8):: dt_grnd	 		   !change in t_grnd, last iteration (Kelvin)
   real(r8):: t_rad_column 		   !radiative temperature (Kelvin)
   real(r8):: t_soisno(-nlevsno+1:nlevsoi) !soil temperature (Kelvin)
   real(r8):: t_lake(1:nlevlak)            !lake temperature (Kelvin)
   real(r8):: tssbef(-nlevsno:nlevsoi)     !soil/snow temperature before update
   real(r8):: t_snow 			   !vertically averaged snow temperature
   real(r8):: thv        		   !virtual potential temperature (kelvin)
   real(r8):: thm 			   !intermediate variable (forc_t+0.0098*forc_hgt_t)
end type column_estate_type

!----------------------------------------------------
! column water state variables structure
!----------------------------------------------------
type column_wstate_type
   type(pft_wstate_type):: pws_a		!pft-level water state variables averaged to the column
   real(r8):: h2osno 			        !snow water (mm H2O)
   real(r8):: h2osoi_liq(-nlevsno+1:nlevsoi)    !liquid water (kg/m2) (new)
   real(r8):: h2osoi_ice(-nlevsno+1:nlevsoi)    !ice lens (kg/m2) (new)
   real(r8):: h2osoi_vol(nlevsoi)		!volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8):: h2osno_old 			!snow mass for previous time step (kg/m2) (new)
   real(r8):: qg 				!ground specific humidity [kg/kg]
   real(r8):: dqgdT   			        !d(qg)/dT
   real(r8):: snowice 		                !average snow ice lens
   real(r8):: snowliq 		                !average snow liquid water
end type column_wstate_type

!----------------------------------------------------
! column carbon state variables structure
!----------------------------------------------------
type column_cstate_type
   type(pft_cstate_type):: pcs_a	!pft-level carbon state variables averaged to the column
   real(r8):: soilc 			!soil carbon (kg C /m**2)
end type column_cstate_type

!----------------------------------------------------
! column nitrogen state variables structure
!----------------------------------------------------
type column_nstate_type
   type(pft_nstate_type):: pns_a	!pft-level nitrogen state variables averaged to the column
end type column_nstate_type

!----------------------------------------------------
! column VOC state variables structure
!----------------------------------------------------
type column_vstate_type
   type(pft_vstate_type):: pvs_a	!pft-level VOC state variables averaged to the column
end type column_vstate_type

!----------------------------------------------------
! column DGVM state variables structure
!----------------------------------------------------
type column_dgvstate_type
   type(pft_dgvstate_type):: pdgvs_a
end type column_dgvstate_type

!----------------------------------------------------
! column dust state variables structure
!----------------------------------------------------
type column_dstate_type
   real(r8):: dummy_entry
end type column_dstate_type

!----------------------------------------------------
! column energy flux variables structure
!----------------------------------------------------
type column_eflux_type
   type(pft_eflux_type):: pef_a	!pft-level energy flux variables averaged to the column
   real(r8):: eflx_snomelt 	!snow melt heat flux (W/m**2)
   real(r8):: eflx_impsoil	!implicit evaporation for soil temperature equation
end type column_eflux_type

!----------------------------------------------------
! column momentum flux variables structure
!----------------------------------------------------
type column_mflux_type
   type(pft_mflux_type)::  pmf_a    !pft-level momentum flux variables averaged to the column
end type column_mflux_type 

!----------------------------------------------------
! column water flux variables structure
!----------------------------------------------------
type column_wflux_type
   type(pft_wflux_type):: pwf_a	!pft-level water flux variables averaged to the column
   real(r8):: qflx_infl		!infiltration (mm H2O /s)
   real(r8):: qflx_surf		!surface runoff (mm H2O /s)
   real(r8):: qflx_drain 	!sub-surface runoff (mm H2O /s)
   real(r8):: qflx_top_soil 	!net water input into soil from top (mm/s)
   real(r8):: qflx_snomelt 	!snow melt (mm H2O /s)
   real(r8):: qflx_qrgwl 	!qflx_surf at glaciers, wetlands, lakes
   real(r8):: qmelt 	    	!snow melt [mm/s]
   real(r8):: qchan2            !history file RTM river (channel) flow (m**3 H2O /s) 
   real(r8):: qchocn2           !history file RTM river (channel) flow into ocean (m**3/s)
end type column_wflux_type

!----------------------------------------------------
! column carbon flux variables structure
!----------------------------------------------------
type column_cflux_type
   type(pft_cflux_type):: pcf_a		!pft-level carbon flux variables averaged to the column
end type column_cflux_type

!----------------------------------------------------
! column nitrogen flux variables structure
!----------------------------------------------------
type column_nflux_type
   type(pft_nflux_type):: pnf_a	        !pft-level nitrogen flux variables averaged to the column
end type column_nflux_type

!----------------------------------------------------
! column VOC flux variables structure
!----------------------------------------------------
type column_vflux_type
   type(pft_vflux_type):: pvf_a		!pft-level VOC flux variables averaged to the column
end type column_vflux_type

!----------------------------------------------------
! column dust flux variables structure
!----------------------------------------------------
type column_dflux_type	
   type(pft_dflux_type):: pdf_a		!pft-level dust flux variables averaged to the column
end type column_dflux_type
!----------------------------------------------------
! End definition of structures defined at the column_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the landunit_type level
!----------------------------------------------------
! landunit physical state variables structure
!----------------------------------------------------
type landunit_pstate_type
   type(gridcell_pstate_type), pointer:: gps !pointer to higher level ps struct
   type(column_pstate_type):: cps_a  !column-level physical state variables averaged to landunit
   logical :: ifspecial              !true=>landunit is not vegetated
   integer :: itype 		     !type for this landunit
   integer :: ncolumns 		     !number of columns for this landunit
   integer :: itypwat 		     !water type
   integer :: isoicol 		     !soil color class
   logical :: lakpoi		     !BOOL: lake point
   real(r8):: bsw   (nlevsoi)        !Clapp and Hornberger "b"
   real(r8):: watsat(nlevsoi)        !volumetric soil water at saturation (porosity)
   real(r8):: hksat (nlevsoi)        !hydraulic conductivity at saturation (mm H2O /s)
   real(r8):: sucsat(nlevsoi)        !minimum soil suction (mm)
   real(r8):: csol  (nlevsoi)        !heat capacity, soil solids (J/m**3/Kelvin)
   real(r8):: tkmg  (nlevsoi)        !thermal conductivity, soil minerals  [W/m-K]  (new)
   real(r8):: tkdry (nlevsoi)        !thermal conductivity, dry soil (W/m/Kelvin)
   real(r8):: tksatu(nlevsoi)        !thermal conductivity, saturated soil [W/m-K]  (new)
   real(r8):: smpmin 		     !restriction for min of soil potential (mm) (new)
   real(r8):: area 		     !total land area for this landunit (km^2)
   real(r8):: wt 		     !weight (relative to gridcell) for this landunit (0-1)
   real(r8):: wtxy 		     !weight (relative to gridcell) for this landunit (0-1)
   integer :: ixy                    !xy lon index
   integer :: jxy                    !xy lat index
   integer :: index1d                !index into 1d global landunit array  
end type landunit_pstate_type

!----------------------------------------------------
! landunit energy state variables structure
!----------------------------------------------------
type landunit_estate_type
   type(column_estate_type)::	ces_a	     !column-level energy state variables averaged to landunit
end type landunit_estate_type

!----------------------------------------------------
! landunit water state variables structure
!----------------------------------------------------
type landunit_wstate_type
   type(column_wstate_type)::	cws_a	      !column-level water state variables averaged to landunit
end type landunit_wstate_type

!----------------------------------------------------
! landunit carbon state variables structure
!----------------------------------------------------
type landunit_cstate_type
   type(column_cstate_type)::	ccs_a	      !column-level carbon state variables averaged to landunit
end type landunit_cstate_type

!----------------------------------------------------
! landunit nitrogen state variables structure
!----------------------------------------------------
type landunit_nstate_type
   type(column_nstate_type)::	cns_a	      !column-level nitrogen state variables averaged to landunit
end type landunit_nstate_type

!----------------------------------------------------
! landunit VOC state variables structure
!----------------------------------------------------
type landunit_vstate_type
   type(column_vstate_type)::	cvs_a	      !column-level VOC state variables averaged to landunit
end type landunit_vstate_type

!----------------------------------------------------
! landunit DGVM state variables structure
!----------------------------------------------------
type landunit_dgvstate_type
   real(r8):: dummy_entry
end type landunit_dgvstate_type

!----------------------------------------------------
! landunit dust state variables structure
!----------------------------------------------------
type landunit_dstate_type
   type(column_dstate_type)::	cds_a		!column-level dust state variables averaged to landunit
end type landunit_dstate_type

!----------------------------------------------------
! landunit energy flux variables structure
!----------------------------------------------------
type landunit_eflux_type
   type(column_eflux_type)::	cef_a		!column-level energy flux variables averaged to landunit
end type landunit_eflux_type

!----------------------------------------------------
! landunit momentum flux variables structure
!----------------------------------------------------
type landunit_mflux_type
   type(pft_mflux_type)::	pmf_a		!pft-level momentum flux variables averaged to landunit
end type landunit_mflux_type

!----------------------------------------------------
! landunit water flux variables structure
!----------------------------------------------------
type landunit_wflux_type
   type(column_wflux_type)::	cwf_a		!column-level water flux variables averaged to landunit
end type landunit_wflux_type

!----------------------------------------------------
! landunit carbon flux variables structure
!----------------------------------------------------
type landunit_cflux_type
   type(column_cflux_type)::	ccf_a		!column-level carbon flux variables averaged to landunit
end type landunit_cflux_type

!----------------------------------------------------
! landunit nitrogen flux variables structure
!----------------------------------------------------
type landunit_nflux_type
   type(column_nflux_type)::	cnf_a		!column-level nitrogen flux variables averaged to landunit
end type landunit_nflux_type

!----------------------------------------------------
! landunit VOC flux variables structure
!----------------------------------------------------
type landunit_vflux_type
   type(pft_vflux_type)::	pvf_a		!pft-level VOC flux variables averaged to landunit
end type landunit_vflux_type

!----------------------------------------------------
! landunit dust flux variables structure
!----------------------------------------------------
type landunit_dflux_type
   type(pft_dflux_type)::	pdf_a		!pft-level dust flux variables averaged to landunit
end type landunit_dflux_type
!----------------------------------------------------
! End definition of structures defined at the landunit_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the gridcell_type level
!----------------------------------------------------
! gridcell physical state variables structure
!----------------------------------------------------
type gridcell_pstate_type
   type(model_pstate_type), pointer:: mps !pointer to higher-level ps struct
   type(column_pstate_type):: cps_a    !column-level physical state variables averaged to gridcell
   real(r8):: lat 	       !latitude (radians)
   real(r8):: lon 	       !longitude (radians)
   real(r8):: latdeg 	       !latitude (degrees)
   real(r8):: londeg 	       !longitude (degrees)
   real(r8):: wtfact 	       !Fraction of model area with high water table
   real(r8):: area 	       !total land area for this gridcell (km^2)
   real(r8):: wt 	       !weight for this gridcell (0-1)
   integer :: nlandunits       !number of landunits for this gridcell
   integer :: itype            !type for this gridcell
   integer :: ixy              !xy lon index
   integer :: jxy              !xy lat index
   integer :: index1d          !index into 1d global gridcell array  
end type gridcell_pstate_type

!----------------------------------------------------
! atmosphere -> land state variables structure
!----------------------------------------------------
type atm2lnd_state_type
   integer :: itypprc 		!precipitatin type
   real(r8):: forc_t 		!atmospheric temperature (Kelvin)
   real(r8):: forc_u 		!atmospheric wind speed in east direction (m/s)
   real(r8):: forc_v 		!atmospheric wind speed in north direction (m/s)
   real(r8):: forc_wind         !atmospheric wind speed   
   real(r8):: forc_q 		!atmospheric specific humidity (kg/kg)
   real(r8):: forc_hgt 		!atmospheric reference height (m)
   real(r8):: forc_hgt_u        !observational height of wind [m] (new)
   real(r8):: forc_hgt_t 	!observational height of temperature [m] (new)
   real(r8):: forc_hgt_q 	!observational height of humidity [m] (new)
   real(r8):: forc_pbot 	!atmospheric pressure (Pa)
   real(r8):: forc_th 		!atmospheric potential temperature (Kelvin)
   real(r8):: forc_vp 		!atmospheric vapor pressure (Pa)
   real(r8):: forc_rho 		!density (kg/m**3)
   real(r8):: forc_co2 		!atmospheric CO2 concentration (Pa)
   real(r8):: forc_o2 		!atmospheric O2 concentration (Pa)
   real(r8):: forc_psrf 	!surface pressure (Pa)
end type atm2lnd_state_type

!----------------------------------------------------
! land -> atmosphere state variables structure
!----------------------------------------------------
type lnd2atm_state_type
   real(r8):: t_veg 		!vegetation temperature (Kelvin)
   real(r8):: t_grnd  		!ground temperature (Kelvin)
   real(r8):: t_rad 		!radiative temperature (Kelvin)
   real(r8):: t_ref2m  		!2 m height surface air temperature (Kelvin)
   real(r8):: t_soisno(-nlevsno+1:nlevsoi)  !soil temperature (Kelvin)
   real(r8):: t_lake(1:nlevlak)	!lake temperature (Kelvin)
   real(r8):: t_snow 		!vertically averaged snow temperature
   real(r8):: dt_veg 		!change in t_veg, last iteration (Kelvin)
   real(r8):: dt_grnd 		!change in t_grnd, last iteration (Kelvin)
   real(r8):: albd(numrad)	!surface albedo (direct)
   real(r8):: albi(numrad)	!surface albedo (diffuse)
   real(r8):: albgrd(numrad)	!ground albedo (direct)
   real(r8):: albgri(numrad)	!ground albedo (diffuse)
end type lnd2atm_state_type 

!----------------------------------------------------
! gridcell energy state variables structure
!----------------------------------------------------
type gridcell_estate_type
   type(column_estate_type)::	ces_a	!column-level energy state variables averaged to gridcell
end type gridcell_estate_type

!----------------------------------------------------
! gridcell water state variables structure
!----------------------------------------------------
type gridcell_wstate_type
   type(column_wstate_type)::	cws_a	!column-level water state variables averaged to gridcell
end type gridcell_wstate_type

!----------------------------------------------------
! gridcell carbon state variables structure
!----------------------------------------------------
type gridcell_cstate_type
   type(column_cstate_type)::	ccs_a	!column-level carbon state variables averaged to gridcell
end type gridcell_cstate_type

!----------------------------------------------------
! gridcell nitrogen state variables structure
!----------------------------------------------------
type gridcell_nstate_type
   type(column_nstate_type)::	cns_a	!column-level nitrogen state variables averaged to gridcell
end type gridcell_nstate_type

!----------------------------------------------------
! gridcell VOC state variables structure
!----------------------------------------------------
type gridcell_vstate_type
   type(column_vstate_type)::	cvs_a	!column-level VOC state variables averaged to gridcell
end type gridcell_vstate_type

!----------------------------------------------------
! gridcell dust state variables structure
!----------------------------------------------------
type gridcell_dstate_type
   type(column_dstate_type)::	cds_a	!column-level dust state variables averaged to gridcell
end type gridcell_dstate_type

!----------------------------------------------------
! gridcell DGVM state variables structure
!----------------------------------------------------
type gridcell_dgvstate_type
   real(r8):: dummy_entry
end type gridcell_dgvstate_type

!----------------------------------------------------
! atmosphere -> land flux variables structure
!----------------------------------------------------
type atm2lnd_flux_type
   real(r8):: forc_lwrad 		!downward infrared (longwave) radiation (W/m**2)
   real(r8):: forc_solad(numrad) 	!direct beam radiation (vis=forc_sols , nir=forc_soll )
   real(r8):: forc_solai(numrad) 	!diffuse radiation     (vis=forc_solsd, nir=forc_solld)
   real(r8):: forc_solar                !incident solar radiation
   real(r8):: forc_rain 		!rain rate [mm/s]
   real(r8):: forc_snow 		!snow rate [mm/s]
!CAS, add Storm Statistics
   real(r8):: forc_strm_int    !intensity of storm event (mm H2O/s)
   real(r8):: forc_strm_dur    !duration of storm event (s)
   real(r8):: forc_strm_dry    !inter-storm arrival period (s)
   real(r8):: forc_strms    !# of events, inclusive
!CAS, add Storm Statistics
end type atm2lnd_flux_type

!----------------------------------------------------
! land -> atmosphere flux variables structure
!----------------------------------------------------
type lnd2atm_flux_type
   real(r8):: taux 		!wind stress: e-w (kg/m/s**2)
   real(r8):: tauy 		!wind stress: n-s (kg/m/s**2)
   real(r8):: eflx_lwrad_out 	!emitted infrared (longwave) radiation (W/m**2)
   real(r8):: eflx_lwrad_net 	!net infrared (longwave) rad (W/m**2) [+ = to atm]
   real(r8):: eflx_sh_tot 	!total sensible heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_sh_veg 	!sensible heat flux from leaves (W/m**2) [+ to atm]
   real(r8):: eflx_sh_grnd 	!sensible heat flux from ground (W/m**2) [+ to atm]
   real(r8):: eflx_lh_tot 	!total latent heat flux (W/m8*2)  [+ to atm]
   real(r8):: eflx_lh_vege 	!veg evaporation heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_lh_vegt 	!veg transpiration heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_lh_grnd 	!ground evaporation heat flux (W/m**2) [+ to atm]
   real(r8):: eflx_soil_grnd 	!soil heat flux (W/m**2) [+ = into soil]
   real(r8):: eflx_snomelt 	!snow melt heat flux (W/m**2)
   real(r8):: eflx_impsoil 	!implicit evaporation for soil temperature equation
end type lnd2atm_flux_type

!----------------------------------------------------
! gridcell energy flux variables structure
!----------------------------------------------------
type gridcell_eflux_type
   type(column_eflux_type)::	cef_a		!column-level energy flux variables averaged to gridcell
end type gridcell_eflux_type

!----------------------------------------------------
! gridcell momentum flux variables structure
!----------------------------------------------------
type gridcell_mflux_type
   type(pft_mflux_type)::	pmf_a		!pft-level momentum flux variables averaged to gridcell
end type gridcell_mflux_type

!----------------------------------------------------
! gridcell water flux variables structure
!----------------------------------------------------
type gridcell_wflux_type
   type(column_wflux_type)::	cwf_a		!column-level water flux variables averaged to gridcell
end type gridcell_wflux_type

!----------------------------------------------------
! gridcell carbon flux variables structure
!----------------------------------------------------
type gridcell_cflux_type
   type(column_cflux_type)::	ccf_a		!column-level carbon flux variables averaged to gridcell
end type gridcell_cflux_type

!----------------------------------------------------
! gridcell nitrogen flux variables structure
!----------------------------------------------------
type gridcell_nflux_type
   type(column_nflux_type)::	cnf_a		!column-level nitrogen flux variables averaged to gridcell
end type gridcell_nflux_type

!----------------------------------------------------
! gridcell VOC flux variables structure
!----------------------------------------------------
type gridcell_vflux_type
   type(pft_vflux_type)::	pvf_a		!pft-level VOC flux variables averaged to gridcell
end type gridcell_vflux_type

!----------------------------------------------------
! gridcell dust flux variables structure
!----------------------------------------------------
type gridcell_dflux_type
   type(pft_dflux_type)::	pdf_a		!pft-level dust flux variables averaged to gridcell
end type gridcell_dflux_type
!----------------------------------------------------
! End definition of structures defined at the gridcell_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the CLM level
!----------------------------------------------------
! CLM physical state variables structure
!----------------------------------------------------
type model_pstate_type
   type(column_pstate_type) ::	cps_a	!column-level physical state variables globally averaged
   real(r8):: area			!total land area for all gridcells (km^2)
   integer :: ngridcells          	!number of gridcells allocated for this process
end type model_pstate_type

!----------------------------------------------------
! CLM energy state variables structure
!----------------------------------------------------
type model_estate_type
   type(column_estate_type)::	ces_a		!column-level energy state variables globally averaged
end type model_estate_type

!----------------------------------------------------
! CLM water state variables structure
!----------------------------------------------------
type model_wstate_type
   type(column_wstate_type)::	cws_a		!column-level water state variables globally averaged
end type model_wstate_type

!----------------------------------------------------
! CLM carbon state variables structure
!----------------------------------------------------
type model_cstate_type
   type(column_cstate_type)::	ccs_a		!column-level carbon state variables globally averaged
end type model_cstate_type

!----------------------------------------------------
! CLM nitrogen state variables structure
!----------------------------------------------------
type model_nstate_type
   type(column_nstate_type)::	cns_a		!column-level nitrogen state variables globally averaged
end type model_nstate_type

!----------------------------------------------------
! CLM VOC state variables structure
!----------------------------------------------------
type model_vstate_type
   type(column_vstate_type)::	cvs_a		!column-level VOC state variables globally averaged
end type model_vstate_type

!----------------------------------------------------
! CLM dust state variables structure
!----------------------------------------------------
type model_dstate_type
   type(column_dstate_type)::	cds_a		!column-level dust state variables globally averaged
end type model_dstate_type

!----------------------------------------------------
! CLM energy flux variables structure
!----------------------------------------------------
type model_eflux_type
   type(column_eflux_type)::	cef_a		!column-level energy flux variables globally averaged
end type model_eflux_type

!----------------------------------------------------
! CLM momentum flux variables structure
!----------------------------------------------------
type model_mflux_type
   type(pft_mflux_type)::	pmf_a		!pft-level momentum flux variables globally averaged
end type model_mflux_type

!----------------------------------------------------
! CLM water flux variables structure
!----------------------------------------------------
type model_wflux_type
   type(column_wflux_type)::	cwf_a		!column-level water flux variables globally averaged
end type model_wflux_type

!----------------------------------------------------
! CLM carbon flux variables structure
!----------------------------------------------------
type model_cflux_type
   type(column_cflux_type)::	ccf_a		!column-level carbon flux variables globally averaged
end type model_cflux_type

!----------------------------------------------------
! CLM nitrogen flux variables structure
!----------------------------------------------------
type model_nflux_type
   type(column_nflux_type)::	cnf_a		!column-level nitrogen flux variables globally averaged
end type model_nflux_type

!----------------------------------------------------
! CLM VOC flux variables structure
!----------------------------------------------------
type model_vflux_type
   type(pft_vflux_type)::	pvf_a		!pft-level VOC flux variables globally averaged
end type model_vflux_type

!----------------------------------------------------
! CLM dust flux variables structure
!----------------------------------------------------
type model_dflux_type
   type(pft_dflux_type)::	pdf_a		!pft-level dust flux variables globally averaged
end type model_dflux_type

!----------------------------------------------------
! End definition of structures defined at the model_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of spatial scaling hierarchy
!----------------------------------------------------
! define the top-level (model) structure 
!----------------------------------------------------
! consists of a pointer to a 1d array of gridcells
! and data that apply to the entire model (i.e. shared by all gridcells)
type model_type
   ! pointers to next level in hierarchy
   type(gridcell_type),dimension(:),pointer:: 	g	!pointer to gridcells
   real(r8),dimension(:),pointer::		ga	!pointer to array of gridcell area (km^2)
   real(r8),dimension(:),pointer::		gw	!pointer to array of gridcell weight (0-1)
   integer,dimension(:),pointer::		gt	!pointer to array of gridcell type

   ! conservation check structures for the clm (global) level
   type(energy_balance_type)   :: mebal	!energy balance structure
   type(water_balance_type)    :: mwbal	!water balance structure
   type(carbon_balance_type)   :: mcbal	!carbon balnace structure
   type(nitrogen_balance_type) :: mnbal	!nitrogen balance structure
   
   ! globally average state variables 
   type(model_pstate_type) ::  mps       !clm physical state variables
   type(model_estate_type) ::  mes       !average of energy states over all gridcells
   type(model_wstate_type) ::  mws       !average of water states over all gridcells
   type(model_cstate_type) ::  mcs       !average of carbon states over all gridcells
   type(model_nstate_type) ::  mns       !average of nitrogen states over all gridcells
   type(model_vstate_type) ::  mvs       !average of VOC states over all gridcells
   type(model_dstate_type) ::  mds       !average of dust states over all gridcells
   
   ! globally averaged flux variables 
   type(model_eflux_type) ::   mef       !average of energy fluxes over all gridcells
   type(model_wflux_type) ::   mwf       !average of water fluxes over all gridcells
   type(model_cflux_type) ::   mcf       !average of carbon fluxes over all gridcells
   type(model_nflux_type) ::   mnf       !average of nitrogen fluxes over all gridcells
   type(model_vflux_type) ::   mvf       !average of VOC fluxes over all gridcells
   type(model_dflux_type) ::   mdf       !average of dust fluxes over all gridcells
end type model_type

!----------------------------------------------------
! define the gridcell structure
!----------------------------------------------------
! consists of a pointer to a 1d array of geomorphological land units occupying a cell
! and data that apply to the entire gridcell (i.e. shared by all landunits)
type gridcell_type
   ! pointer to next level in hierarchy
   type(landunit_type),dimension(:),pointer::l	!pointer to geomorphological landunits
   real(r8),dimension(:),pointer::	la	!pointer to array of landunit area (km^2)
   real(r8),dimension(:),pointer::	lw	!pointer to array of landunit weight (0-1)
   integer,dimension(:),pointer::	lt	!pointer to array of landunit type

   ! conservation check structures for the gridcell level
   type(energy_balance_type) ::	gebal		!energy balance structure
   type(water_balance_type) ::	gwbal		!water balance structure
   type(carbon_balance_type) ::	gcbal		!carbon balance structure
   type(nitrogen_balance_type) ::gnbal		!nitrogen balance structure
   
   ! state variables defined at the gridcell level
   type(gridcell_pstate_type) ::   gps      !gridcell physical state variables
   type(gridcell_estate_type) ::   ges      !average of energy states over all landunits
   type(gridcell_wstate_type) ::   gws      !average of water states over all landunits
   type(gridcell_cstate_type) ::   gcs      !average of carbon states over all landunits
   type(gridcell_nstate_type) ::   gns      !average of nitrogen states over all landunits
   type(gridcell_vstate_type) ::   gvs      !average of VOC states over all landunits
   type(gridcell_dstate_type) ::   gds      !average of dust states over all landunits
   type(atm2lnd_state_type) ::	a2ls		!atmospheric state variables required by the land
   type(lnd2atm_state_type) ::	l2as		!land state variables required by the atmosphere
   
   ! flux variables defined at the gridcell level
   type(gridcell_eflux_type) ::	gef		!average of energy fluxes over all landunits
   type(gridcell_wflux_type) ::	gwf		!average of water fluxes over all landunits
   type(gridcell_cflux_type) ::	gcf		!average of carbon fluxes over all landunits
   type(gridcell_nflux_type) ::	gnf		!average of nitrogen fluxes over all landunits
   type(gridcell_vflux_type) ::	gvf		!average of VOC fluxes over all landunits
   type(gridcell_dflux_type) ::	gdf		!average of dust fluxes over all landunits
   type(atm2lnd_flux_type) ::	a2lf		!atmospheric flux variables required by the land
   type(lnd2atm_flux_type) ::	l2af		!land flux variables required by the atmosphere
end type gridcell_type

!----------------------------------------------------
! define the geomorphological land unit structure
!----------------------------------------------------
! consists of a pointer to a 1d array of columns occupying a land unit
! and data that apply to the entire landunit (i.e. shared by all columns)
type landunit_type
   ! pointer to lower level in hierarchy
   type(column_type),dimension(:),pointer::c	!pointer to soil/snow/canopy columns
   real(r8),dimension(:),pointer::	ca	!pointer to array of column area (km^2)
   real(r8),dimension(:),pointer::	cw	!pointer to array of column weight (0-1)
   integer,dimension(:),pointer::	ct	!pointer to array of column type
   
   ! pointers to higher level in hierarchy
   type(atm2lnd_state_type),pointer::	a2ls	!atmospheric state variables required by the land
   type(atm2lnd_flux_type),pointer::	a2lf	!atmospheric flux variables required by the land
   
   ! conservation check structures for the landunit level
   type(energy_balance_type) ::	lebal		!energy balance structure
   type(water_balance_type) ::	lwbal		!water balance structure
   type(carbon_balance_type) ::	lcbal		!carbon balance structure
   type(nitrogen_balance_type) ::lnbal		!nitrogen balance structure
   
   ! state variables defined at the land unit level
   type(landunit_pstate_type) ::   lps   !land unit physical state variables
   type(landunit_estate_type) ::   les   !average of energy states over all columns
   type(landunit_wstate_type) ::   lws   !average of water states over all columns
   type(landunit_cstate_type) ::   lcs   !average of carbon states over all columns
   type(landunit_nstate_type) ::   lns   !average of nitrogen states over all columns
   type(landunit_vstate_type) ::   lvs   !average of VOC states over all columns
   type(landunit_dstate_type) ::   lds   !average of dust states over all columns
   
   ! flux variables defined at the landunit level
   type(landunit_eflux_type) :: lef      !average of energy fluxes over all columns
   type(landunit_mflux_type) :: lmf      !average of momentum fluxes over all columns
   type(landunit_wflux_type) :: lwf      !average of water fluxes over all columns
   type(landunit_cflux_type) :: lcf      !average of carbon fluxes over all columns
   type(landunit_nflux_type) :: lnf      !average of nitrogen fluxes over all columns
   type(landunit_vflux_type) :: lvf      !average of VOC fluxes over all columns
   type(landunit_dflux_type) :: ldf      !average of dust fluxes over all columns
end type landunit_type

!----------------------------------------------------
! define the column structure
!----------------------------------------------------
! consists of a pointer to a 1d array of plant types occupying a column
! and data that appy to the entire column (i.e. shared by all plant types)
type column_type
   ! pointer to the next level in hierarcy
   type(pft_type), dimension(:),pointer:: p !pointer to plant functional types (pfts)
   real(r8), dimension(:), pointer :: pa    !pointer to array of pft area (km^2)
   real(r8), dimension(:), pointer :: pw    !pointer to array of pft weight (0-1)
   integer , dimension(:), pointer :: pt    !pointer to array of pft type
   
   ! pointers to higher level in hierarchy
   type(atm2lnd_state_type), pointer ::	a2ls !atmospheric state variables required by the land
   type(atm2lnd_flux_type) , pointer ::	a2lf !atmospheric flux variables required by the land
   
   ! conservation check structures for the column level
   type(energy_balance_type) ::	cebal	!energy balance structure
   type(water_balance_type) :: 	cwbal	!water balance structure
   type(carbon_balance_type) ::	ccbal	!carbon balance structure
   type(nitrogen_balance_type) :: cnbal	!nitrogen balance structure
   
   ! state variables defined at the column level
   type(column_pstate_type) :: cps   !column physical state variables
   type(column_estate_type) :: ces   !column energy state
   type(column_wstate_type) :: cws   !column water state
   type(column_cstate_type) :: ccs   !column carbon state
   type(column_nstate_type) :: cns   !column nitrogen state
   type(column_vstate_type) :: cvs   !column VOC state
   type(column_dstate_type) :: cds   !column dust state
   
   ! flux variables defined at the column level
   type(column_eflux_type) :: cef    !column energy flux
   type(column_mflux_type) :: cmf    !column momentum flux
   type(column_wflux_type) :: cwf    !column water flux
   type(column_cflux_type) :: ccf    !column carbon flux
   type(column_nflux_type) :: cnf    !column nitrogen flux
   type(column_vflux_type) :: cvf    !column VOC flux
   type(column_dflux_type) :: cdf    !column dust flux

   ! dgvm variables defined at the column level
   type (column_dgvstate_type) :: cdgvs !column DGVM structure
end type column_type

!----------------------------------------------------
! define the pft structure
!----------------------------------------------------
! no pointer to array of lower types - this is the bottom of the 
! type hierarchy. Only data common to the pft, and pointers to higher levels
type pft_type
   ! pointer to higher level in hierarchy
   type(atm2lnd_state_type), pointer ::	a2ls !atmospheric state variables required by the land
   type(atm2lnd_flux_type) , pointer ::	a2lf !atmospheric flux variables required by the land
 
   ! conservation check structures for the pft level
   type(energy_balance_type)   :: pebal	     !energy balance structure
   type(water_balance_type)    :: pwbal	     !water balance structure
   type(carbon_balance_type)   :: pcbal	     !carbon balance structure
   type(nitrogen_balance_type) :: pnbal	     !nitrogen balance structure
   
   ! pft parameters and ecophysiological variables 
   type(pft_epc_type)   , pointer:: pepc     !pointer to array of ecophysiological constants 
   type(pft_dgvepc_type), pointer:: pdgvepc  !pointer to array of DGVM ecophysiological constants 

   ! DGVM state variables
   type(pft_dgvstate_type) :: pdgvs	     !pft DGVM state variables
   
   ! state variables defined at the pft level
   type(pft_pstate_type) :: pps              !physical state variables
   type(pft_estate_type) :: pes              !pft energy state
   type(pft_wstate_type) :: pws              !pft water state
   type(pft_cstate_type) :: pcs              !pft carbon state
   type(pft_nstate_type) :: pns              !pft nitrogen state
   type(pft_vstate_type) :: pvs              !pft VOC state

   ! flux variables defined at the pft level
   type(pft_eflux_type)  :: pef              !pft energy flux
   type(pft_mflux_type)  :: pmf              !pft momentum flux
   type(pft_wflux_type)  :: pwf              !pft water flux
   type(pft_cflux_type)  :: pcf              !pft carbon flux
   type(pft_nflux_type)  :: pnf              !pft nitrogen flux
   type(pft_vflux_type)  :: pvf              !pft VOC flux
   type(pft_dflux_type)  :: pdf              !pft dust flux
end type pft_type
!----------------------------------------------------
! End definition of spatial scaling hierarchy
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Declare single instance of clmtype
!----------------------------------------------------
type(model_type), target, save :: clm

!----------------------------------------------------
! Declare single instance of array of ecophysiological constant types
!----------------------------------------------------
type(pft_epc_type), target, save :: pftcon(0:numpft)

!----------------------------------------------------
! Declare single instance of array of dgvm ecophysiological constant types
!----------------------------------------------------
type(pft_dgvepc_type), target, save :: dgv_pftcon(0:numpft)

!----------------------------------------------------
! Declare data types for 1d mapping
!----------------------------------------------------
type gridcell_1d
   character(len=16) :: name     !name of level type
   integer :: beg                !per-proc beginning land index (min val is 1)
   integer :: end                !per-proc ending land index (max val is total gridcells)
   integer :: num                !total number of land grid cells (over all procs)
   real(r8),pointer :: wtxy(:)   !weight relative to grid cell
   integer, pointer :: ixy(:)    !longitude indices 
   integer, pointer :: jxy(:)    !latitude indices 
   real(r8),pointer :: latdeg(:) !latitudes in degrees north
   real(r8),pointer :: londeg(:) !longitudes in degrees east
end type gridcell_1d

type landunit_1d
   character(len=16) :: name     !name of level type
   integer :: beg                !per-proc beginning land index (min val is 1)
   integer :: end                !per-proc ending land index (max val is total landunits)
   integer :: num                !total number of land units (over all procs)
   real(r8),pointer :: wtxy(:)   !weight relative to grid cell
   integer, pointer :: ixy(:)    !longitude indices 
   integer, pointer :: jxy(:)    !latitude indices 
   integer, pointer :: gindex(:) !1d grid indices
   real(r8),pointer :: latdeg(:) !latitudes in degrees north
   real(r8),pointer :: londeg(:) !longitudes in degrees east
   integer, pointer :: itypwat(:)!water type (soil,lake, wetland, glackier)
end type landunit_1d

type column_1d
   character(len=16) :: name     !name of level type
   integer :: beg                !per-proc beginning column index (min val is 1)
   integer :: end                !per-proc ending column index (max val is total columns)
   integer :: num                !total number of columns 
   real(r8),pointer :: wtxy(:)   !weight relative to gridcell
   real(r8),pointer :: wtlnd(:)  !weight relative to landunit
   integer, pointer :: ixy(:)    !2d longitude indices 
   integer, pointer :: jxy(:)    !2d latitude indices 
   integer, pointer :: gindex(:) !1d grid indices
   integer, pointer :: lindex(:) !1d landunit indices
   real(r8),pointer :: latdeg(:) !latitudes in degrees north
   real(r8),pointer :: londeg(:) !longitudes in degrees east
   integer, pointer :: itypwat(:)!water type (soil,lake, wetland, glackier)
end type column_1d

type pft_1d
   character(len=16) :: name     !name of level type
   integer :: beg                !per-proc beginning pft index (min value is 1)
   integer :: end                !per-proc ending pft index (max val is total pfts)
   integer :: num                !total number of vegetated patches
   real(r8),pointer :: wtxy(:)   !weight relative to grid cell
   real(r8),pointer :: wtlnd(:)  !weight relative to landunit
   real(r8),pointer :: wtcol(:)  !weight relative to column
   integer, pointer :: ixy(:)    !grid longitude indices
   integer, pointer :: jxy(:)    !grid latitude indices
   integer, pointer :: mxy(:)    !grid patch index
   integer, pointer :: gindex(:) !1d grid indices
   integer, pointer :: lindex(:) !1d landunit indices
   integer, pointer :: cindex(:) !1d column indices
   real(r8),pointer :: latdeg(:) !latitudes in degrees north
   real(r8),pointer :: londeg(:) !longitudes in degrees east
   integer, pointer :: itypwat(:)!water type (soil,lake, wetland, glackier)
   integer, pointer :: itypveg(:)!vegetation type
end type pft_1d

!----------------------------------------------------
! Declare single instance of 1d mapping types
!----------------------------------------------------
type(gridcell_1d), save :: grid1d
type(landunit_1d), save :: land1d
type(column_1d)  , save :: cols1d
type(pft_1d)     , save :: pfts1d
!
!EOP
!----------------------------------------------------------------------- 

end module clmtype
