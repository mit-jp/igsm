#include <misc.h>
#include <preproc.h>

module Biogeophysics2Mod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: Biogeophysics2Mod
! 
! !DESCRIPTION: 
! Performs the calculation of soil/snow and ground temperatures 
! and updates surface fluxes based on the new ground temperature. 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics2   ! Calculate soil/snow and ground temperatures
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Biogeophysics2
!
! !INTERFACE:
  subroutine Biogeophysics2 (c)
!
! !DESCRIPTION: 
! This is the main subroutine to execute the calculation of soil/snow and
! ground temperatures and update surface fluxes based on the new ground
! temperature 
! Calling sequence is:
! Biogeophysics2:                    surface biogeophysics driver
!    -> SoilTemperature:             soil/snow and ground temperatures      
!          -> SoilTermProp           thermal conductivities and heat 
!                                     capacities        
!          -> Tridiagonal            tridiagonal matrix solution            
!          -> PhaseChange            phase change of liquid/ice contents        
!
! (1) Snow and soil temperatures
!     o The volumetric heat capacity is calculated as a linear combination 
!       in terms of the volumetric fraction of the constituent phases. 
!     o The thermal conductivity of soil is computed from 
!       the algorithm of Johansen (as reported by Farouki 1981), and the 
!       conductivity of snow is from the formulation used in
!       SNTHERM (Jordan 1991).
!     o Boundary conditions:  
!       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
!     o Soil / snow temperature is predicted from heat conduction 
!       in 10 soil layers and up to 5 snow layers. 
!       The thermal conductivities at the interfaces between two 
!       neighboring layers (j, j+1) are derived from an assumption that 
!       the flux across the interface is equal to that from the node j 
!       to the interface and the flux from the interface to the node j+1. 
!       The equation is solved using the Crank-Nicholson method and 
!       results in a tridiagonal system equation.
!
! (2) Phase change (see PhaseChange.F90)
!
! !USES:
    use clmtype
    use globals
    use clm_varcon, only : hvap, cpair, grav, vkc, tfrz, sb
    use clm_varpar, only : nlevsoi
    use SoilTemperatureMod
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! Migrated to clm2.0 by Keith Oleson and Mariana Vertenstein
! Migrated to clm2.1 new data structures by Peter Thornton and M. Vertenstein 
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: snl             !number of snow layers
    logical , pointer :: do_capsnow      !true => do snow capping
    real(r8), pointer :: forc_lwrad      !downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: emg             !ground emissivity
    real(r8), pointer :: htvp            !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: t_grnd          !ground temperature (Kelvin)
    integer , pointer :: frac_veg_nosno  !fraction of vegetation not covered by snow (0 OR 1 now) [-]
    real(r8), pointer :: cgrnds          !deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl          !deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: sabg            !solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: dlrad           !downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad           !upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: eflx_sh_veg     !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_veg   !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg   !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can   !evaporation from leaves and stems (mm H2O/s) (+ = to atm)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: eflx_sh_grnd    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi   !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_snowcap    !excess precipitation due to snow capping (mm H2O /s)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: dt_grnd         !change in t_grnd, last iteration (Kelvin)
    real(r8), pointer :: eflx_soil_grnd  !soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_sh_tot     !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_tot   !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: eflx_lh_tot     !total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: qflx_evap_grnd  !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow   !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow   !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd   !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: eflx_lwrad_out  !emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net  !net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lh_vege    !veg evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt    !veg transpiration heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd    !ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: t_rad_pft       !pft-level radiative temperature (Kelvin)
    real(r8), pointer :: t_rad_column    !column-level radiative temperature (Kelvin)
    real(r8), pointer :: errsoi          !soil/lake energy conservation error (W/m**2)
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: wt              !pft weight relative to column
!
! local pointers to implicit in arrays
!
    real(r8), dimension(:), pointer:: tssbef     !soil/snow temperature before update
    real(r8), dimension(:), pointer:: t_soisno   !soil temperature (Kelvin)
    real(r8), dimension(:), pointer:: h2osoi_ice !ice lens (kg/m2) (new)
    real(r8), dimension(:), pointer:: h2osoi_liq !liquid water (kg/m2) (new)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j                           !do loop index
    real(r8) :: fact(c%cps%snl+1 : nlevsoi) !used in computing tridiagonal matrix
    real(r8) :: egsmax  !max. evaporation which soil can provide at one time step
    real(r8) :: egirat  !ratio of topsoil_evap_tot : egsmax
    real(r8) :: xmf     !total latent heat of phase change of ground water
    real(r8) :: tinc    !temperature difference of two time step
    integer  :: pi                     !pft index
    real(r8) :: column_eflx_lwrad_out  !column-average of elfx_lwrad_out
    real(r8) :: topsoil_evap_tot       !column-level total evaporation from top soil layer
    real(r8) :: save_qflx_evap_soi     !temporary storage for qflx_evap_soi
    real(r8) :: pftsum                 !temporary used for pft averaging for columns
    real(r8) :: sumwt                  !temporary 
    real(r8) :: evaprat                !ratio of qflx_evap_soi/topsoil_evap_tot 
    type(atm2lnd_flux_type) ,pointer :: a2lf ! local pointers to derived subtypes
    type(column_pstate_type),pointer :: cps  ! local pointers to derived subtypes
    type(column_estate_type),pointer :: ces  ! local pointers to derived subtypes
    type(column_wstate_type),pointer :: cws  ! local pointers to derived subtypes
    type(pft_type)          ,pointer :: p    ! local pointers to derived subtypes
    type(pft_pstate_type)   ,pointer :: pps  ! local pointers to derived subtypes
    type(pft_estate_type)   ,pointer :: pes  ! local pointers to derived subtypes
    type(pft_eflux_type)    ,pointer :: pef  ! local pointers to derived subtypes
    type(pft_wflux_type)    ,pointer :: pwf  ! local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    a2lf => c%a2lf
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! Assign local pointers to derived subtypes components (column-level)
    forc_lwrad => a2lf%forc_lwrad
    snl => cps%snl
    do_capsnow => cps%do_capsnow
    htvp => cps%htvp
    emg => cps%emg
    t_grnd => ces%t_grnd
    dt_grnd => ces%dt_grnd
    t_rad_column => ces%t_rad_column  

    ! Assign local pointers to implicit in arrays (column-level)
    t_soisno => ces%t_soisno
    tssbef => ces%tssbef
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! Determine soil temperatures including surface soil temperature
    call SoilTemperature(c, xmf , fact)

    ! calculate difference in soil temperature from last time step, for
    ! flux corrections
    tinc = t_soisno(snl+1) - tssbef(snl+1)

    ! a preliminary pft loop to determine if corrections are required for
    ! excess evaporation from the top soil layer... Includes new logic
    ! to distribute the corrections between pfts on the basis of their
    ! evaporative demands. 
    ! egirat holds the ratio of demand to availability if demand is
    ! greater than availability, or 1.0 otherwise.

    egsmax = (h2osoi_ice(snl+1)+h2osoi_liq(snl+1)) / dtime
    ! added to trap very small negative soil water,ice 
    if (egsmax < 0._r8) then
       egsmax = 0._r8
    endif

    topsoil_evap_tot = 0._r8
    sumwt = 0._r8
    do pi=1,cps%npfts

       ! Assign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pef => p%pef
       pwf => p%pwf

       ! Assign local pointers to derived subtypes components (pft-level)
       cgrnds => pef%cgrnds
       cgrndl => pef%cgrndl
       eflx_sh_grnd => pef%eflx_sh_grnd
       qflx_evap_soi => pwf%qflx_evap_soi
       wt => pps%wt

       ! Correct fluxes to present soil temperature
       eflx_sh_grnd =  eflx_sh_grnd + tinc*cgrnds
       qflx_evap_soi =  qflx_evap_soi + tinc*cgrndl

       ! Set the column-average qflx_evap_soi as the weighted average over all pfts
       ! but only count the pfts that are evaporating
       if (qflx_evap_soi > 0._r8) then
          topsoil_evap_tot = topsoil_evap_tot + qflx_evap_soi * wt
          sumwt = sumwt + wt
       endif

    end do ! end preliminary pfts loop

    if (sumwt > 0._r8) then
       topsoil_evap_tot = topsoil_evap_tot / sumwt
    endif

    ! calculate ratio for rescaling pft-level fluxes to meet availability
    if (topsoil_evap_tot > egsmax) then
       egirat = (egsmax/topsoil_evap_tot)
    else
       egirat = 1.0
    end if

    ! begin second pfts loop
    ! also using this loop to calculate the column-average eflx_lwrad_out
    ! for calculating column-level t_rad
    column_eflx_lwrad_out = 0._r8
    do pi=1,cps%npfts

       ! Assign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pes => p%pes
       pef => p%pef
       pwf => p%pwf

       ! Assign local pointers to derived subtypes components (pft-level)
       wt => pps%wt
       frac_veg_nosno => pps%frac_veg_nosno
       sabg => pef%sabg
       dlrad => pef%dlrad
       ulrad => pef%ulrad
       eflx_sh_grnd => pef%eflx_sh_grnd
       eflx_sh_veg => pef%eflx_sh_veg
       qflx_evap_soi => pwf%qflx_evap_soi
       qflx_evap_veg => pwf%qflx_evap_veg
       qflx_tran_veg => pwf%qflx_tran_veg
       qflx_evap_can => pwf%qflx_evap_can
       qflx_snowcap => pwf%qflx_snowcap
       qflx_evap_tot => pwf%qflx_evap_tot  
       qflx_evap_grnd => pwf%qflx_evap_grnd  
       qflx_sub_snow => pwf%qflx_sub_snow  
       qflx_dew_snow => pwf%qflx_dew_snow  
       qflx_dew_grnd => pwf%qflx_dew_grnd  
       eflx_soil_grnd => pef%eflx_soil_grnd  
       eflx_sh_tot => pef%eflx_sh_tot  
       eflx_lh_tot => pef%eflx_lh_tot  
       eflx_lwrad_out => pef%eflx_lwrad_out  
       eflx_lwrad_net => pef%eflx_lwrad_net  
       eflx_lh_vege => pef%eflx_lh_vege    
       eflx_lh_vegt => pef%eflx_lh_vegt    
       eflx_lh_grnd => pef%eflx_lh_grnd    
       errsoi => p%pebal%errsoi  
       t_rad_pft => pes%t_rad_pft  

       ! calculate the ratio of soil evaporation for this pft relative 
       ! to the total of all evaporating PFTs
       if (qflx_evap_soi > 0._r8) then
          evaprat = qflx_evap_soi/topsoil_evap_tot
       else
          evaprat = 0._r8
       endif

       ! correct soil fluxes for possible evap in excess of top layer water 
       ! excess energy is added to the sensible heat flux from soil
       if (egirat < 1.0 .AND. qflx_evap_soi > 0._r8) then
          save_qflx_evap_soi = qflx_evap_soi
          qflx_evap_soi = qflx_evap_soi * egirat
          eflx_sh_grnd = eflx_sh_grnd + (save_qflx_evap_soi - qflx_evap_soi)*htvp
       end if

       ! Ground heat flux
       eflx_soil_grnd = sabg + dlrad + (1-frac_veg_nosno)*emg*forc_lwrad &
            - emg*sb*tssbef(snl+1)**3*(tssbef(snl+1) + 4.*tinc) &
            - (eflx_sh_grnd+qflx_evap_soi*htvp)

       ! Total fluxes (vegetation + ground)
       eflx_sh_tot = eflx_sh_veg + eflx_sh_grnd
       qflx_evap_tot = qflx_evap_veg + qflx_evap_soi
       eflx_lh_tot= hvap*qflx_evap_veg + htvp*qflx_evap_soi   

       ! Assign ground evaporation to sublimation from soil ice or to dew
       ! on snow or ground 
       qflx_evap_grnd = 0.
       qflx_sub_snow = 0.
       qflx_dew_snow = 0.
       qflx_dew_grnd = 0.

       if (qflx_evap_soi >= 0.) then
          qflx_evap_grnd = min(evaprat*h2osoi_liq(snl+1)/dtime, qflx_evap_soi)
          qflx_sub_snow = qflx_evap_soi - qflx_evap_grnd
       else
          if (t_grnd < tfrz) then
             qflx_dew_snow = abs(qflx_evap_soi)
          else
             qflx_dew_grnd = abs(qflx_evap_soi)
          endif
       endif

       ! Update the pft-level qflx_snowcap
       ! This was moved in from Hydrology2 to keep all pft-level
       ! calculations out of Hydrology2
       if (snl < 0 .and. do_capsnow) then
          qflx_snowcap = qflx_snowcap + qflx_dew_snow + qflx_dew_grnd
       endif

       ! Outgoing long-wave radiation from vegetation + ground
       ! For conservation we put the increase of ground longwave to outgoing
       eflx_lwrad_out = ulrad &
            + (1-frac_veg_nosno)*(1.-emg)*forc_lwrad &
            + (1-frac_veg_nosno)*emg*sb * tssbef(snl+1)**4 &
            + 4.*emg*sb*tssbef(snl+1)**3*tinc

       ! pft-level Radiative temperature
       t_rad_pft = (eflx_lwrad_out/sb)**0.25

       ! weighted average of pft-level lwrad_out for column-level t_rad
       column_eflx_lwrad_out = column_eflx_lwrad_out + eflx_lwrad_out * wt

       ! Soil Energy balance check
       errsoi = 0.
       do j = snl+1, nlevsoi
          errsoi = errsoi - (t_soisno(j)-tssbef(j))/fact(j)
       enddo
       errsoi = errsoi + eflx_soil_grnd - xmf

       ! Variables needed by history tape
       qflx_evap_can  = qflx_evap_veg - qflx_tran_veg
       eflx_lh_vege   = (qflx_evap_veg - qflx_tran_veg) * hvap
       eflx_lh_vegt   = qflx_tran_veg * hvap
       eflx_lh_grnd   = qflx_evap_soi * htvp
       eflx_lwrad_net = eflx_lwrad_out - forc_lwrad

    end do ! end of second pfts loop

    ! column-level Radiative temperature

    t_rad_column = (column_eflx_lwrad_out/sb)**0.25

    ! averaging for pft energy balance variables
    ! lake balance for errsoi is not over pft - therefore
    ! this calculation is done here rather than in pft_to_column 

    pftsum = 0.0
    do pi=1,cps%npfts
       pftsum = pftsum + c%p(pi)%pebal%errsoi * c%pw(pi)   
    end do
    c%cebal%errsoi = pftsum

  end subroutine Biogeophysics2

end module Biogeophysics2Mod
