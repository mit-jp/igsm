#include <misc.h>
#include <preproc.h>

module Biogeophysics1Mod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: Biogeophysics1Mod
! 
! !DESCRIPTION: 
! Performs calculation of leaf temperature and surface fluxes. 
! Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics1   ! Calculate leaf temperature and surface fluxes
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
! !IROUTINE: Biogeophysics1
!
! !INTERFACE:
  subroutine Biogeophysics1 (c)
!
! !DESCRIPTION: 
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
!  Biogeophysics1:                   surface biogeophysics driver
!    -> QSat:                        saturated vapor pressure, specific humidity,
!                                     and derivatives at ground surface
!    -> SurfaceRadiation:            surface solar radiation
!    -> BareGroundFluxes:            surface fluxes for bare soil or
!                                     snow-covered vegetation patches
!          -> MoninObukIni:          first-guess Monin-Obukhov length and
!                                     wind speed
!          -> FrictionVelocity:      friction velocity and potential 
!                                     temperature and humidity profiles
!    -> CanopyFluxes:                leaf temperature and surface fluxes
!                                     for vegetated patches 
!          -> QSat                   saturated vapor pressure, specific humidity,
!                                     and derivatives at leaf surface
!          -> MoninObukIni           first-guess Monin-Obukhov length and 
!                                     wind speed
!          -> FrictionVelocity       friction velocity and potential
!                                     temperature and humidity profiles
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for sunlit leaves
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for shaded leaves
!          -> QSat                   recalculation of saturated vapor pressure,
!                                     specific humidity, and derivatives at
!                                     leaf surface using updated leaf temperature
!  Leaf temperature
!  Foliage energy conservation is given by the foliage energy budget 
!  equation:
!                 Rnet - Hf - LEf = 0 
!  The equation is solved by Newton-Raphson iteration, in which this 
!  iteration includes the calculation of the photosynthesis and 
!  stomatal resistance, and the integration of turbulent flux profiles. 
!  The sensible and latent heat transfer between foliage and atmosphere 
!  and ground is linked by the equations:  
!                 Ha = Hf + Hg and Ea = Ef + Eg
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon, only : denh2o, denice, roverg, hvap, hsub, istice, istwet, &
         zlnd, zsno
    use clm_varpar, only : nlevsoi
    use QSatMod, only : QSat
    use SurfaceRadiationMod, only : SurfaceRadiation
    use BareGroundFluxesMod, only : BareGroundFluxes
    use CanopyFluxesMod, only : CanopyFluxes
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
! local pointers to implicit in scalars
!
    real(r8),pointer:: forc_pbot    !atmospheric pressure (Pa)
    real(r8),pointer:: forc_q       !atmospheric specific humidity (kg/kg)
    real(r8),pointer:: forc_t       !atmospheric temperature (Kelvin)
    real(r8),pointer:: forc_hgt_t   !observational height of temperature [m]
    real(r8),pointer:: forc_th      !atmospheric potential temperature (Kelvin)
    real(r8),pointer:: forc_u       !atmospheric wind speed in east direction (m/s)
    real(r8),pointer:: forc_v       !atmospheric wind speed in north direction (m/s)
    integer,pointer:: itypwat       !water type
    real(r8),pointer:: smpmin       !restriction for min of soil potential (mm)
    integer,pointer :: snl          !number of snow layers
    real(r8),pointer:: frac_sno     !fraction of ground covered by snow (0 to 1)
    real(r8),pointer:: h2osno       !snow water (mm H2O)
    real(r8),pointer:: elai         !one-sided leaf area index with burying by snow
    real(r8),pointer:: esai         !one-sided stem area index with burying by snow
    real(r8),pointer:: z0mr         !ratio of momentum roughness length to canopy top height (-)
    real(r8),pointer:: htop         !canopy top (m)
    real(r8),pointer:: displar      !ratio of displacement height to canopy top height (-)
    integer,pointer :: frac_veg_nosno  !fraction of vegetation not covered by snow (0 OR 1 now) [-]
!
! local pointers to implicit out scalars
!
    real(r8),pointer:: t_grnd       !ground temperature (Kelvin)
    real(r8),pointer:: qg           !ground specific humidity [kg/kg]
    real(r8),pointer:: dqgdT        !d(qg)/dT
    real(r8),pointer:: emg          !ground emissivity
    real(r8),pointer:: htvp         !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8),pointer:: beta         !coefficient of convective velocity [-]
    real(r8),pointer:: zii          !convective boundary height [m]
    real(r8),pointer:: thm          !intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8),pointer:: thv          !virtual potential temperature (kelvin)
    real(r8),pointer:: z0mg         !roughness length over ground, momentum [m]
    real(r8),pointer:: z0hg         !roughness length over ground, sensible heat [m]
    real(r8),pointer:: z0qg         !roughness length over ground, latent heat [m]
    real(r8),pointer:: emv          !vegetation emissivity
    real(r8),pointer:: z0m          !momentum roughness length (m)
    real(r8),pointer:: displa       !displacement height (m)
    real(r8),pointer:: z0mv         !roughness length over vegetation, momentum [m]
    real(r8),pointer:: z0hv         !roughness length over vegetation, sensible heat [m]
    real(r8),pointer:: z0qv         !roughness length over vegetation, latent heat [m]
    real(r8),pointer:: ur           !wind speed at reference height [m/s]
    real(r8),pointer:: psnsun       !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8),pointer:: psnsha       !shaded leaf photosynthesis (umol CO2 /m**2/ s)
!
! local pointers to implicit in arrays
!
    real(r8),dimension(:),pointer:: dz           !layer depth (m)
    real(r8),dimension(:),pointer:: t_soisno     !soil temperature (Kelvin)
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2)
    real(r8),dimension(:),pointer:: watsat       !volumetric soil water at saturation (porosity)
    real(r8),dimension(:),pointer:: sucsat       !minimum soil suction (mm)
    real(r8),dimension(:),pointer:: bsw          !Clapp and Hornberger "b"
!
! local pointers to implicit out arrays
!
    real(r8),dimension(:),pointer:: tssbef       !soil/snow temperature before update
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer:: i,j      ! loop indices
    real(r8):: qred    ! soil surface relative humidity
    real(r8):: avmuir  ! ir inverse optical depth per unit leaf area
    real(r8):: eg      ! water vapor pressure at temperature T [pa]
    real(r8):: qsatg   ! saturated humidity [kg/kg]
    real(r8):: degdT   ! d(eg)/dT
    real(r8):: qsatgdT ! d(qsatg)/dT
    real(r8):: fac     ! soil wetness of surface layer
    real(r8):: psit    ! negative potential of soil
    real(r8):: hr      ! relative humidity
    real(r8):: wx      ! partial volume of ice and water of surface layer
    integer :: pi           ! pft index 	
    type(atm2lnd_state_type)  ,pointer :: a2ls ! local pointers to derived subtypes
    type(landunit_pstate_type),pointer :: lps  ! local pointers to derived subtypes
    type(column_pstate_type),pointer :: cps    ! local pointers to derived subtypes
    type(column_estate_type),pointer :: ces    ! local pointers to derived subtypes
    type(column_wstate_type),pointer :: cws    ! local pointers to derived subtypes
    type(column_eflux_type) ,pointer :: cef    ! local pointers to derived subtypes
    type(column_wflux_type) ,pointer :: cwf    ! local pointers to derived subtypes
    type(pft_type)       , pointer :: p        ! local pointers to derived subtypes
    type(pft_pstate_type), pointer :: pps      ! local pointers to derived subtypes
    type(pft_eflux_type) , pointer :: pef      ! local pointers to derived subtypes
    type(pft_wflux_type) , pointer :: pwf      ! local pointers to derived subtypes
    type(pft_epc_type)   , pointer :: pepc     ! local pointers to derived subtypes
    type(pft_cflux_type) , pointer :: pcf      ! local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes (column-level)
    a2ls => c%a2ls
    cps => c%cps
    ces => c%ces
    cws => c%cws
    cef => c%cef
    cwf => c%cwf
    lps => cps%lps

    ! assign local pointers to derived type scalar members (column-level)
    forc_pbot => a2ls%forc_pbot 
    forc_q => a2ls%forc_q   
    forc_t => a2ls%forc_t    
    forc_hgt_t => a2ls%forc_hgt_t 
    forc_th => a2ls%forc_th    
    forc_u => a2ls%forc_u     
    forc_v => a2ls%forc_v     
    itypwat => lps%itypwat
    smpmin => lps%smpmin
    snl => cps%snl
    frac_sno => cps%frac_sno
    h2osno => cws%h2osno
    t_grnd => ces%t_grnd
    qg => cws%qg
    dqgdT => cws%dqgdT
    emg => cps%emg
    htvp => cps%htvp
    z0mg => cps%z0mg
    z0hg => cps%z0hg
    z0qg => cps%z0qg
    beta => cps%beta
    zii => cps%zii
    ur => cps%ur
    thm => ces%thm
    thv => ces%thv

    ! assign local pointers to derived type array members (column-level)
    watsat => lps%watsat
    sucsat => lps%sucsat
    bsw => lps%bsw
    dz => cps%dz
    t_soisno => ces%t_soisno
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice
    tssbef => ces%tssbef 

    ! begin calculations that relate only to the column level

    ! Ground and soil temperatures from previous time step
    t_grnd = t_soisno(snl+1)
    do i = snl+1, nlevsoi
       tssbef(i) = t_soisno(i)
    enddo

    ! Saturated vapor pressure, specific humidity and their derivatives
    ! at ground surface
    qred = 1.
    if (itypwat/=istwet .AND. itypwat/=istice) then
       wx   = (h2osoi_liq(1)/denh2o+h2osoi_ice(1)/denice)/dz(1)
       fac  = min(1._r8, wx/watsat(1))
       fac  = max( fac, 0.01_r8 )
       psit = -sucsat(1) * fac ** (-bsw(1))
       psit = max(smpmin, psit)
       hr   = exp(psit/roverg/t_grnd)
       qred = (1.-frac_sno)*hr + frac_sno
    endif

    call QSat(t_grnd, forc_pbot, eg, degdT, qsatg, qsatgdT)

    qg = qred*qsatg
    dqgdT = qred*qsatgdT

    if (qsatg > forc_q .AND. forc_q > qred*qsatg) then
       qg = forc_q
       dqgdT = 0.
    endif

    ! Ground Emissivity
    if (h2osno>0. .OR.itypwat==istice) then
       emg = 0.97
    else
       emg = 0.96
    endif

    ! Latent heat. We arbitrarily assume that the sublimation occurs 
    ! only as h2osoi_liq = 0
    htvp = hvap
    if (h2osoi_liq(snl+1) <= 0. .AND. h2osoi_ice(snl+1) > 0.) htvp = hsub

    ! Switch between vaporization and sublimation causes rapid solution
    ! separation in perturbation growth test
#if (defined PERGRO)
    htvp = hvap
#endif

    ! Roughness lengths over bare ground
    if (frac_sno > 0.) then
       z0mg = zsno
       z0hg = z0mg            ! initial set only
       z0qg = z0mg            ! initial set only
    else
       z0mg = zlnd
       z0hg = z0mg            ! initial set only
       z0qg = z0mg            ! initial set only
    endif

    ! Potential, virtual potential temperature, and wind speed at the 
    ! reference height
    beta=1.
    zii = 1000.
    thm = forc_t + 0.0098*forc_hgt_t
    thv = forc_th*(1.+0.61*forc_q)
    ur = max(1.0_r8,sqrt(forc_u*forc_u+forc_v*forc_v))


    ! begin calculations that relate to both column and pft

    ! Surface Radiation
    ! (this routine includes an internal pft loop)
    call SurfaceRadiation (c)

    ! begin loop through pfts
    do pi=1,cps%npfts
       ! assign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pwf => p%pwf
       pef => p%pef
       pepc => p%pepc
       pcf => p%pcf

       ! assign local pointers to derived type scalar members (fpt-level)
       elai => pps%elai
       esai => pps%esai
       htop => pps%htop
       frac_veg_nosno => pps%frac_veg_nosno
       z0mr => pepc%z0mr
       displar => pepc%displar
       emv => pps%emv
       z0m => pps%z0m
       displa => pps%displa
       z0mv => pps%z0mv
       z0hv => pps%z0hv
       z0qv => pps%z0qv
       psnsun => pcf%psnsun
       psnsha => pcf%psnsha

       ! Initial set (needed for history tape fields)
       pef%eflx_sh_tot = 0._r8
       pef%eflx_lh_tot = 0._r8
       pef%eflx_sh_veg = 0._r8  
       pwf%qflx_evap_tot = 0._r8
       pwf%qflx_evap_veg = 0._r8  
       pwf%qflx_tran_veg = 0._r8  

       ! Initial set for calculation 
       pef%cgrnd  = 0._r8
       pef%cgrnds = 0._r8
       pef%cgrndl = 0._r8

       ! Vegetation Emissivity   
       avmuir=1.
       emv=1.-exp(-(elai+esai)/avmuir)

       ! Roughness lengths over vegetation
       z0m = z0mr*htop
       displa = displar*htop
       z0mv = z0m
       z0hv = z0mv
       z0qv = z0mv

       ! Surface Temperature and Fluxes

       if (frac_veg_nosno == 0) then
          ! BARE SOIL OR SNOW-COVERED VEGETATION
          ! Ground fluxes
          ! NOTE: in the current scheme frac_veg_nosno is EITHER 1 or 0

          call BareGroundFluxes (c, p)
          psnsun = 0.
          psnsha = 0. !put these lines here to avoid psn = NaN
       else
          ! VEGETATION
          ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
          ! and leaf water change by evapotranspiration
          call CanopyFluxes (c, p)
       endif

    end do ! (end of pfts loop)

  end subroutine Biogeophysics1

end module Biogeophysics1Mod
