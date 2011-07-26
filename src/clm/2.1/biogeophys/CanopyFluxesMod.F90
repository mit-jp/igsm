#include <misc.h>
#include <preproc.h>

module CanopyFluxesMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CanopyFluxesMod
! 
! !DESCRIPTION: 
! Calculates the leaf temperature and the leaf fluxes, 
! transpiration, photosynthesis and  updates the dew
! accumulation due to evaporation.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxes  !Calculates the leaf temperature and leaf fluxes 
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: Stomata      !Leaf stomatal resistance and leaf photosynthesis
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
! !IROUTINE: CanopyFluxes 
!
! !INTERFACE:
  subroutine CanopyFluxes (c, p)
!
! !DESCRIPTION: 
! 1. Calculates the leaf temperature: 
! 2. Calculates the leaf fluxes, transpiration, photosynthesis and 
!    updates the dew accumulation due to evaporation.
!
! Method:
! Use the Newton-Raphson iteration to solve for the foliage 
! temperature that balances the surface energy budget:
!
! f(t_veg) = Net radiation - Sensible - Latent = 0
! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
!
! Note:
! (1) In solving for t_veg, t_grnd is given from the previous timestep.
! (2) The partial derivatives of aerodynamical resistances, which cannot 
!     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
! (3) The weighted stomatal resistance of sunlit and shaded foliage is used 
! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
!                                                          => Ec + Eg = Ea
! (5) Energy loss is due to: numerical truncation of energy budget equation
!     (*); and "ecidif" (see the code) which is dropped into the sensible 
!     heat 
! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n) and 
!     del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference of 
!     water flux from the leaf between the iteration step (n+1) and (n) 
!     less than 0.1 W/m2; or the iterative steps over 40.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals
    use clm_varcon, only : sb, cpair, hvap, vkc, grav, denice, denh2o, tfrz, &
         csoilc
    use QSatMod, only : QSat
    use FrictionVelocityMod, only : MoninObukIni, FrictionVelocity
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout)   :: c     !column-level data
    type (pft_type), target, intent(inout)      :: p     !pft-level data
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 12/19/01, Peter Thornton
! Changed tg to t_grnd for consistency with other routines
! 1/29/02, Peter Thornton
! Migrate to new data structures, new calling protocol. For now co2 and
! o2 partial pressures are hardwired, but they should be coming in from
! forc_co2 and forc_o2. Keeping the same hardwired values as in CLM2 to
! assure bit-for-bit results in the first comparisons.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: th          !atmospheric potential temperature (Kelvin)
    real(r8), pointer :: t_grnd      !ground surface temperature [K]
    real(r8), pointer :: thm         !intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8), pointer :: qg          !specific humidity at ground surface [kg/kg]
    real(r8), pointer :: thv         !virtual potential temperature (kelvin)
    real(r8), pointer :: z0mv        !roughness length over vegetation, momentum [m]  
    real(r8), pointer :: z0hv        !roughness length over vegetation, sensible heat [m]
    real(r8), pointer :: z0qv        !roughness length over vegetation, latent heat [m]
    real(r8), pointer :: dqgdT       !temperature derivative of "qg"
    real(r8), pointer :: htvp        !latent heat of evaporation (/sublimation) [J/kg]
    real(r8), pointer :: emv         !ground emissivity
    real(r8), pointer :: emg         !vegetation emissivity
    real(r8), pointer :: forc_pbot   !atmospheric pressure (Pa)
    real(r8), pointer :: forc_q      !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_u      !atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v      !atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_hgt    !atmospheric reference height (m)
    real(r8), pointer :: forc_hgt_u  !observational height of wind [m]
    real(r8), pointer :: forc_hgt_t  !observational height of temperature [m]
    real(r8), pointer :: forc_hgt_q  !observational height of humidity [m]
    real(r8), pointer :: forc_rho    !density (kg/m**3)
    real(r8), pointer :: forc_co2    !atmospheric CO2 concentration (Pa)
    real(r8), pointer :: forc_o2     !atmospheric O2 concentration (Pa)
    real(r8), pointer :: forc_lwrad  !downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: smpmax      !wilting point potential in mm
    real(r8), pointer :: displa      !displacement height (m)
    real(r8), pointer :: dleaf       !characteristic leaf dimension (m)
    real(r8), pointer :: parsun      !average absorbed PAR for sunlit leaves (W/m**2)
    real(r8), pointer :: parsha      !average absorbed PAR for shaded leaves (W/m**2)
    real(r8), pointer :: qe25        !quantum efficiency at 25C (umol CO2 / umol photon)
    real(r8), pointer :: vcmx25      !max rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8), pointer :: mp          !slope of conductance-to-photosynthesis relationship
    real(r8), pointer :: c3psn       !photosynthetic pathway: 0. = c4, 1. = c3
    real(r8), pointer :: elai        !one-sided leaf area index with burying by snow  
    real(r8), pointer :: esai        !one-sided stem area index with burying by snow
    real(r8), pointer :: fdry        !fraction of foliage that is green and dry [-]
    real(r8), pointer :: fwet        !fraction of canopy that is wet (0 to 1)
    real(r8), pointer :: laisun      !sunlit leaf area
    real(r8), pointer :: laisha      !shaded leaf area
    integer , pointer :: frac_veg_nosno  !frac of veg not covered by snow (0 OR 1 now) [-]
    real(r8), pointer :: sabv        !solar radiation absorbed by vegetation (W/m**2)
!
! local pointers to implicit inout scalars
!
    real(r8), pointer :: cgrnds      !deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl      !deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: t_veg       !vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m     !2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_af        !canopy air temperature (Kelvin)
    real(r8), pointer :: h2ocan      !canopy water (mm H2O)
    real(r8), pointer :: annpsnpot   !annual potential photosynthesis (umol CO2 /m**2)
    real(r8), pointer :: annpsn      !annual photosynthesis (umol CO2 /m**2)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: cgrnd       !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: dlrad       !downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad       !upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: u10         !10-m wind (m/s) (for dust model)
    real(r8), pointer :: fv          !friction velocity (m/s) (for dust model)
    real(r8), pointer :: ram1        !aerodynamical resistance (s/m)
    real(r8), pointer :: btran       !transpiration wetness factor (0 to 1)
    real(r8), pointer :: rssun       !sunlit stomatal resistance (s/m)
    real(r8), pointer :: rssha       !shaded stomatal resistance (s/m)
    real(r8), pointer :: psnsun      !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha      !shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: qflx_tran_veg   !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: dt_veg      !change in t_veg, last iteration (Kelvin)
    real(r8), pointer :: qflx_evap_veg   !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: eflx_sh_veg !sensible heat flux from leaves (W/m**2) [+ to atm]        
    real(r8), pointer :: taux        !wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy        !wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi   !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: efpot       !potential latent energy flux [kg/m2/s]
!
! local pointers to implicit in arrays
!
    real(r8), dimension(:), pointer :: watsat     !volumetric soil water at saturation (porosity)
    real(r8), dimension(:), pointer :: h2osoi_ice !ice lens (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_liq !liquid water (kg/m2)
    real(r8), dimension(:), pointer :: dz         !layer depth (m)
    real(r8), dimension(:), pointer :: t_soisno   !soil temperature (Kelvin)
    real(r8), dimension(:), pointer :: sucsat     !minimum soil suction (mm)
    real(r8), dimension(:), pointer :: bsw        !Clapp and Hornberger "b"
    real(r8), dimension(:), pointer :: rootfr     !fraction of roots in each soil layer
!
! local pointers to implicit out arrays
!
    real(r8), dimension(:), pointer :: rootr      !effective fraction of roots in each soil layer
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    real(r8):: zldis       !reference height "minus" zero displacement height [m]
    real(r8):: zii         !convective boundary layer height [m]
    real(r8):: zeta        !dimensionless height used in Monin-Obukhov theory
    real(r8):: beta        !coefficient of conective velocity [-]
    real(r8):: wc          !convective velocity [m/s]
    real(r8):: dth         !diff of virtual temp. between ref. height and surface 
    real(r8):: dthv        !diff of vir. poten. temp. between ref. height and surface
    real(r8):: dqh         !diff of humidity between ref. height and surface
    real(r8):: obu         !Monin-Obukhov length (m)
    real(r8):: um          !wind speed including the stablity effect [m/s]
    real(r8):: ur          !wind speed at reference height [m/s]
    real(r8):: uaf         !velocity of air within foliage [m/s]
    real(r8):: temp1       !relation for potential temperature profile
    real(r8):: temp2       !relation for specific humidity profile
    real(r8):: ustar       !friction velocity [m/s]
    real(r8):: tstar       !temperature scaling parameter
    real(r8):: qstar       !moisture scaling parameter
    real(r8):: thvstar     !virtual potential temperature scaling parameter
    real(r8):: taf         !air temperature within canopy space [K]
    real(r8):: qaf         !humidity of canopy air [kg/kg]
    real(r8):: rpp         !fraction of potential evaporation from leaf [-]
    real(r8):: rppdry      !fraction of potential evaporation through transp [-]
    real(r8):: cf          !heat transfer coefficient from leaves [-]
    real(r8):: rb          !leaf boundary layer resistance [s/m]
    real(r8):: ram(2)      !aerodynamical resistance [s/m]
    real(r8):: rah(2)      !thermal resistance [s/m]
    real(r8):: raw(2)      !moisture resistance [s/m]
    real(r8):: wta         !heat conductance for air [m/s]
    real(r8):: wtg         !heat conductance for ground [m/s]
    real(r8):: wtl         !heat conductance for leaf [m/s]
    real(r8):: wta0        !normalized heat conductance for air [-]
    real(r8):: wtl0        !normalized heat conductance for leaf [-]
    real(r8):: wtg0        !normalized heat conductance for ground [-]
    real(r8):: wtal        !normalized heat conductance for air and leaf [-]
    real(r8):: wtgl        !normalized heat conductance for leaf and ground [-]
    real(r8):: wtga        !normalized heat cond. for air and ground  [-]
    real(r8):: wtaq        !latent heat conductance for air [m/s]
    real(r8):: wtlq        !latent heat conductance for leaf [m/s]
    real(r8):: wtgq        !latent heat conductance for ground [m/s]
    real(r8):: wtaq0       !normalized latent heat conductance for air [-]
    real(r8):: wtlq0       !normalized latent heat conductance for leaf [-]
    real(r8):: wtgq0       !normalized heat conductance for ground [-]
    real(r8):: wtalq       !normalized latent heat cond. for air and leaf [-]
    real(r8):: wtglq       !normalized latent heat cond. for leaf and ground [-]
    real(r8):: wtgaq       !normalized latent heat cond. for air and ground [-]
    real(r8):: el          !vapor pressure on leaf surface [pa]
    real(r8):: deldT       !derivative of "el" on "t_veg" [pa/K]
    real(r8):: qsatl       !leaf specific humidity [kg/kg]
    real(r8):: qsatldT     !derivative of "qsatl" on "t_veg"
    real(r8):: air,bir,cir !atmos. radiation temporay set
    real(r8):: dc1,dc2     !derivative of energy flux [W/m2/K]
    real(r8):: delt        !temporary
    real(r8):: delq        !temporary
    real(r8):: del         !absolute change in leaf temp in current iteration [K]
    real(r8):: del2        !change in leaf temperature in previous iteration [K]
    real(r8):: dele        !change in latent heat flux from leaf [K]
    real(r8):: delmax      !maximum change in  leaf temperature [K]
    real(r8):: dels        !change in leaf temperature in current iteration [K]
    real(r8):: det         !maximum leaf temp. change in two consecutive iter [K]
    real(r8):: dlemin      !max limit for energy flux convergence [w/m2]
    real(r8):: dtmin       !max limit for temperature convergence [K]
    real(r8):: efeb        !latent heat flux from leaf (previous iter) [mm/s]
    real(r8):: efeold      !latent heat flux from leaf (previous iter) [mm/s]
!CAS    real(r8):: efpot       !potential latent energy flux [kg/m2/s]
    real(r8):: efe         !water flux from leaf [mm/s]
    real(r8):: efsh        !sensible heat from leaf [mm/s]
    real(r8):: obuold      !monin-obukhov length from previous iteration
    real(r8):: tlbef       !leaf temperature from previous iteration [K]
    real(r8):: ecidif      !excess energies [W/m2]
    real(r8):: err         !balance error
    real(r8):: erre        !balance error
    real(r8):: co2         !atmospheric co2 concentration (pa)
    real(r8):: o2          !atmospheric o2 concentration (pa)
    real(r8):: svpts       !saturation vapor pressure at t_veg (pa)
    real(r8):: eah         !canopy air vapor pressure (pa)
    real(r8):: s_node      !vol_liq/eff_porosity
    real(r8):: smp_node    !matrix potential
    real(r8):: vol_ice(1:nlevsoi)      !partial volume of ice lens in layer
    real(r8):: eff_porosity(1:nlevsoi) !effective porosity in layer
    real(r8):: vol_liq(1:nlevsoi)      !partial volume of liquid water in layer
    real(r8):: rresis(1:nlevsoi)       !soil water contribution to root resistance
    real(r8):: po2         !partial pressure  o2 (mol/mol)
    real(r8):: pco2        !partial pressure co2 (mol/mol)
    data po2,pco2 /0.209,355.e-06/
    real(r8):: mpe = 1.e-6 ! prevents overflow error if division by zero
    integer :: i            !loop index
    integer :: itlef        !counter for leaf temperature iteration [-]
    integer :: itmax        !maximum number of iteration [-]
    integer :: itmin        !minimum number of iteration [-]
    integer :: nmozsgn      !number of times stability changes sign
    real(r8):: fm           !needed for BGC only to diagnose 10m wind speed
    real(r8):: wtshi        !sensible heat resistance for air, grnd and leaf [-]
    real(r8):: wtsqi        !latent heat resistance for air, grnd and leaf [-]
    real(r8), target :: annpsnpot_junk = 0. ! when DGVM not used
    real(r8), target :: annpsn_junk = 0.    ! when DGVM not used
    type(atm2lnd_state_type), pointer :: a2ls  !local pointers to derived subtypes
    type(atm2lnd_flux_type) , pointer :: a2lf  !local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps   !local pointers to derived subtypes
    type(column_estate_type), pointer :: ces   !local pointers to derived subtypes
    type(column_wstate_type), pointer :: cws   !local pointers to derived subtypes
    type(pft_pstate_type)   , pointer :: pps   !local pointers to derived subtypes
    type(pft_dgvstate_type) , pointer :: pdgvs !local pointers to derived subtypes
    type(pft_epc_type)      , pointer :: pepc  !local pointers to derived subtypes
    type(pft_estate_type)   , pointer :: pes   !local pointers to derived subtypes
    type(pft_wstate_type)   , pointer :: pws   !local pointers to derived subtypes
    type(pft_mflux_type)    , pointer :: pmf   !local pointers to derived subtypes
    type(pft_eflux_type)    , pointer :: pef   !local pointers to derived subtypes 
    type(pft_wflux_type)    , pointer :: pwf   !local pointers to derived subtypes
    type(pft_cflux_type)    , pointer :: pcf   !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes
    a2ls  => c%a2ls
    a2lf  => c%a2lf
    cps   => c%cps
    ces   => c%ces
    cws   => c%cws
    pps   => p%pps
    pepc  => p%pepc
    pes   => p%pes
    pws   => p%pws
    pmf   => p%pmf
    pef   => p%pef
    pwf   => p%pwf
    pcf   => p%pcf

    ! assign local pointers to derived type array members
    watsat      => cps%lps%watsat
    h2osoi_ice  => cws%h2osoi_ice
    h2osoi_liq  => cws%h2osoi_liq
    dz          => cps%dz
    t_soisno    => ces%t_soisno
    sucsat      => cps%lps%sucsat
    bsw         => cps%lps%bsw
    rootfr      => pps%rootfr
    rootr       => pps%rootr

    ! assign local pointers to derived type scalar members
    th     => a2ls%forc_th
    forc_pbot   => a2ls%forc_pbot
    forc_q      => a2ls%forc_q
    forc_u      => a2ls%forc_u
    forc_v      => a2ls%forc_v
    forc_hgt    => a2ls%forc_hgt
    forc_hgt_u  => a2ls%forc_hgt_u
    forc_hgt_t  => a2ls%forc_hgt_t
    forc_hgt_q  => a2ls%forc_hgt_q
    forc_rho    => a2ls%forc_rho
    forc_co2    => a2ls%forc_co2
    forc_o2     => a2ls%forc_o2
    forc_lwrad  => a2lf%forc_lwrad
    thm         => ces%thm
    htvp        => cps%htvp
    emg         => cps%emg
    t_grnd      => ces%t_grnd
    thv         => ces%thv
    qg          => cws%qg
    dqgdT       => cws%dqgdT
    z0mv        => pps%z0mv
    z0hv        => pps%z0hv
    z0qv        => pps%z0qv
    emv         => pps%emv
    displa      => pps%displa
    elai        => pps%elai
    esai        => pps%esai
    fdry        => pps%fdry
    fwet        => pps%fwet
    laisun      => pps%laisun
    laisha      => pps%laisha
    frac_veg_nosno => pps%frac_veg_nosno
    smpmax      => pepc%smpmax
    dleaf       => pepc%dleaf
    qe25        => pepc%qe25
    vcmx25      => pepc%vcmx25
    mp          => pepc%mp
    c3psn       => pepc%c3psn
    parsun      => pef%parsun
    parsha      => pef%parsha
    sabv        => pef%sabv
    eflx_sh_veg => pef%eflx_sh_veg
    t_veg       => pes%t_veg
    t_ref2m     => pes%t_ref2m
    t_af        => pes%t_af
    h2ocan      => pws%h2ocan
    cgrnds      => pef%cgrnds
    cgrndl      => pef%cgrndl
    annpsnpot => annpsnpot_junk
    annpsn => annpsn_junk
    u10 => pps%u10        
    fv => pps%fv         
    ram1 => pps%ram1
    btran => pps%btran   
    rssun => pps%rssun   
    rssha => pps%rssha   
    dt_veg => pps%dt_veg 
    cgrnd => pef%cgrnd   
    eflx_sh_veg => pef%eflx_sh_veg
    eflx_sh_grnd => pef%eflx_sh_grnd
    dlrad => pef%dlrad 
    ulrad => pef%ulrad 
    taux => pmf%taux   
    tauy => pmf%tauy   
    qflx_tran_veg => pwf%qflx_tran_veg
    qflx_evap_veg => pwf%qflx_evap_veg
    qflx_evap_soi => pwf%qflx_evap_soi
    efpot => pwf%qflx_efpot
    psnsun => pcf%psnsun
    psnsha => pcf%psnsha

    ! Initialize
    del   = 0.0  ! change in leaf temperature from previous iteration
    itlef = 0    ! counter for leaf temperature iteration
    efeb  = 0.0  ! latent head flux from leaf for previous iteration
    wtlq = 0.0
    wtlq0 = 0.0
    wtgq0 = 0.0
    wtalq = 0.0
    wtgaq = 0.0
    wtglq = 0.0
    wtaq = 0.0
    wtgq = 0.0
    wtaq0 = 0.0
    wtlq0 = 0.0
    wtgq0 = 0.0
    wtalq = 0.0
    wtgaq = 0.0
    wtglq = 0.0


    ! Assign iteration parameters
    delmax = 1.0  ! maximum change in  leaf temperature
    itmax  = 40   ! maximum number of iteration
    itmin  = 2    ! minimum number of iteration
    dtmin  = 0.01 ! max limit for temperature convergence
    dlemin = 0.1  ! max limit for energy flux convergence

    ! Effective porosity of soil, partial volume of ice and liquid (needed
    ! for btran)
    do i = 1,nlevsoi
       vol_ice(i) = min(watsat(i), h2osoi_ice(i)/(dz(i)*denice))
       eff_porosity(i) = watsat(i)-vol_ice(i)
       vol_liq(i) = min(eff_porosity(i), h2osoi_liq(i)/(dz(i)*denh2o))
    enddo

    ! Root resistance factors
    btran = 1.e-10
    do i = 1,nlevsoi
       if (t_soisno(i) > tfrz) then
          s_node = max(vol_liq(i)/eff_porosity(i),0.01_r8)
          smp_node = max(smpmax, -sucsat(i)*s_node**(-bsw(i)))
          rresis(i) = (1.-smp_node/smpmax)/(1.+sucsat(i)/smpmax)
          rootr(i) = rootfr(i)*rresis(i)
          btran = btran + rootr(i)
       else
          rootr(i) = 0.
       endif
    enddo

    ! Normalize root resistances to get layer contribution to ET
    do i = 1,nlevsoi
       rootr(i)  = rootr(i)/btran
    enddo

    ! Net absorbed longwave radiation by canopy and ground
    ! =air+bir*t_veg**4+cir*t_grnd**4
    air =   emv * (1.+(1.-emv)*(1.-emg)) * forc_lwrad
    bir = - (2.-emv*(1.-emg)) * emv * sb
    cir =   emv*emg*sb

    ! Saturated vapor pressure, specific humidity, and their derivatives
    ! at the leaf surface
    call QSat (t_veg, forc_pbot, el, deldT, qsatl, qsatldT)

    ! Determine atmospheric co2 and o2
    ! for bit-for-bit results with clm2, using the hardwired local co2 and o2
    !co2 = forc_co2
    !o2 = forc_o2
    co2 = pco2*forc_pbot
    o2  = po2*forc_pbot

    ! Initialize flux profile
    nmozsgn = 0
    obuold = 0.
    zii=1000.         ! m  (pbl height)
    beta=1.           ! -  (in computing W_*)

    taf = (t_grnd + thm)/2.
    qaf = (forc_q+qg)/2.

    ur = max(1.0_r8,sqrt(forc_u*forc_u+forc_v*forc_v))
    dth = thm-taf
    dqh = forc_q-qaf
    dthv = dth*(1.+0.61*forc_q)+0.61*th*dqh
    zldis = forc_hgt_u - displa

    ! Initialize Monin-Obukhov length and wind speed
    call MoninObukIni (ur, thv, dthv, zldis, z0mv, um, obu  )

    ! Begin stability iteration
    ITERATION : do while (itlef <= itmax) 

       tlbef = t_veg
       del2 = del

       ! Determine friction velocity, and potential temperature and humidity
       ! profiles of the surface boundary layer
       call FrictionVelocity (displa, z0mv,  z0hv,  z0qv,  obu, &
            itlef+1, ur, um, ustar, temp1, &
            temp2, forc_hgt, forc_hgt_u, forc_hgt_t, forc_hgt_q, &
            u10, fm, fv)

       ! Determine aerodynamic resistances
       ram(1)=1./(ustar*ustar/um)
       rah(1)=1./(temp1*ustar) 
       raw(1)=1./(temp2*ustar) 
       ram1 = ram(1) !pass value to global variable

       ! Bulk boundary layer resistance of leaves
       uaf = um*sqrt( 1./(ram(1)*um) )
       cf = 0.01/(sqrt(uaf)*sqrt(dleaf))
       rb = 1./(cf*uaf)

       ! Aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
       ! if no vegetation, rah(2)=0 because zpd+z0h = z0hg.
       ! (Dickinson et al., 1993, pp.54)
       ram(2) = 0.               ! not used
       rah(2) = 1./(csoilc*uaf)
       raw(2) = rah(2) 

       ! Stomatal resistances for sunlit and shaded fractions of canopy.
       ! Done each iteration to account for differences in eah, tv.
       svpts = el                        ! pa
       eah = forc_pbot * qaf / 0.622     ! pa

       call Stomata (mpe, parsun, svpts, eah, thm, o2, co2, btran, rb, &
            rssun, psnsun, qe25, vcmx25, mp, c3psn, forc_pbot, t_veg, &
            annpsnpot, annpsn)

       call Stomata (mpe, parsha, svpts, eah, thm, o2, co2, btran, rb, &
            rssha, psnsha, qe25, vcmx25, mp, c3psn, forc_pbot, t_veg, &
            annpsnpot, annpsn)

       ! Sensible heat conductance for air, leaf and ground 
       ! Moved the original subroutine in-line...
       wta   = 1./rah(1)             ! air
       wtl   = (elai+esai)/rb        ! leaf
       wtg   = 1./rah(2)             ! ground
       wtshi = 1./(wta+wtl+wtg)

       wtl0  = wtl*wtshi         ! leaf
       wtg0  = wtg*wtshi         ! ground
       wta0  = wta*wtshi         ! air

       wtgl  = wtl0+wtg0         ! ground + leaf
       wtga  = wta0+wtg0         ! ground + air
       wtal  = wta0+wtl0         ! air + leaf

       ! Fraction of potential evaporation from leaf
       if (fdry .gt. 0.0) then
          rppdry  = fdry*rb*(laisun/(rb+rssun) + laisha/(rb+rssha))/elai
       else
          rppdry = 0.0
       endif
       efpot = forc_rho*wtl*(qsatl-qaf)

       if (efpot > 0.) then
          if (btran > 1.e-10) then
             qflx_tran_veg = efpot*rppdry
             rpp = rppdry + fwet
          else
             !No transpiration if btran below 1.e-10
             rpp = fwet
             qflx_tran_veg = 0.
          endif
          !Check total evapotranspiration from leaves
          rpp = min(rpp, (qflx_tran_veg+h2ocan/dtime)/efpot)
       else
          !No transpiration if potential evaporation less than zero
          rpp = 1.
          qflx_tran_veg = 0.
       endif

       ! Update conductances for changes in rpp 
       ! Latent heat conductances for ground and leaf.
       ! Air has same conductance for both sensible and latent heat.
       ! Moved the original subroutine in-line...
       wtaq  = frac_veg_nosno/raw(1)               ! air
       wtlq  = frac_veg_nosno*(elai+esai)/rb * rpp ! leaf
       wtgq  = frac_veg_nosno/raw(2)               ! ground
       wtsqi = 1./(wtaq+wtlq+wtgq)

       wtgq0 = wtgq*wtsqi   ! ground
       wtlq0 = wtlq*wtsqi   ! leaf
       wtaq0 = wtaq*wtsqi   ! air

       wtglq = wtgq0+wtlq0  ! ground + leaf
       wtgaq = wtaq0+wtgq0  ! air + ground
       wtalq = wtaq0+wtlq0  ! air + leaf

       dc1 = forc_rho*cpair*wtl
       dc2 = hvap*forc_rho*wtlq

       efsh = dc1*(wtga*t_veg-wtg0*t_grnd-wta0*thm)
       efe = dc2*(wtgaq*qsatl-wtgq0*qg-wtaq0*forc_q)

       ! Evaporation flux from foliage
       erre = 0.
       if (efe*efeb < 0.0) then
          efeold = efe
          efe  = 0.1*efeold
          erre = efe - efeold
       endif
       dt_veg = (sabv + air + bir*t_veg**4 + cir*t_grnd**4 - efsh - efe) &
            / (- 4.*bir*t_veg**3 +dc1*wtga +dc2*wtgaq*qsatldT)
       t_veg = tlbef + dt_veg
       dels = dt_veg
       del  = abs(dels)
       err = 0.
       if (del > delmax) then
          dt_veg = delmax*dels/del
          t_veg = tlbef + dt_veg
          err = sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
               + cir*t_grnd**4 - (efsh + dc1*wtga*dt_veg) &
               - (efe + dc2*wtgaq*qsatldT*dt_veg)
       endif

       ! Fluxes from leaves to canopy space
       ! "efe" was limited as its sign changes frequently.  This limit may
       ! result in an imbalance in "hvap*qflx_evap_veg" and "efe + dc2*wtgaq*qsatldT*dt_veg" 
       efpot = forc_rho*wtl*(wtgaq*(qsatl+qsatldT*dt_veg) &
            -wtgq0*qg-wtaq0*forc_q)
       qflx_evap_veg = rpp*efpot

       ! Calculation of evaporative potentials (efpot) and
       ! interception losses; flux in kg m**-2 s-1.  ecidif 
       ! holds the excess energy if all intercepted water is evaporated
       ! during the timestep.  This energy is later added to the
       ! sensible heat flux.
       ecidif = 0.
       if (efpot > 0. .AND. btran > 1.e-10) then
          qflx_tran_veg = efpot*rppdry
       else
          qflx_tran_veg = 0.
       endif
       ecidif = max(0._r8, qflx_evap_veg-qflx_tran_veg-h2ocan/dtime)
       qflx_evap_veg = min(qflx_evap_veg,qflx_tran_veg+h2ocan/dtime)

       ! The energy loss due to above two limits is added to 
       ! the sensible heat flux.
       eflx_sh_veg = efsh + dc1*wtga*dt_veg + err + erre + hvap*ecidif

       ! Re-calculate saturated vapor pressure, specific humidity, and their
       ! derivatives at the leaf surface
       call QSat (t_veg, forc_pbot, el, deldT, qsatl, qsatldT)

       ! Update vegetation/ground surface temperature, canopy air temperature, 
       ! canopy vapor pressure, aerodynamic temperature, and
       ! Monin-Obukhov stability parameter for next iteration. 
       taf = wtg0*t_grnd + wta0*thm + wtl0*t_veg
       qaf = wtlq0*qsatl + wtgq0*qg + forc_q*wtaq0

       ! Update Monin-Obukhov length and wind speed including the stability effect
       dth = thm-taf       
       dqh = forc_q-qaf

       tstar=temp1*dth
       qstar=temp2*dqh

       thvstar=tstar*(1.+0.61*forc_q) + 0.61*th*qstar
       zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)

       if (zeta >= 0.) then     !stable
          zeta = min(2._r8,max(zeta,0.01_r8))
          um = max(ur,0.1_r8)
       else                     !unstable
          zeta = max(-100._r8,min(zeta,-0.01_r8))
          wc = beta*(-grav*ustar*thvstar*zii/thv)**0.333
          um = sqrt(ur*ur+wc*wc)
       endif
       obu = zldis/zeta

       if (obuold*obu < 0.) nmozsgn = nmozsgn+1
       if (nmozsgn >= 4) then 
          obu = zldis/(-0.01)
       endif

       obuold = obu

       ! Test for convergence
       itlef = itlef+1
       if (itlef > itmin) then
          dele = abs(efe-efeb)
          efeb = efe
          det  = max(del,del2)
          if (det < dtmin .AND. dele < dlemin) exit 
       endif

    enddo ITERATION     ! End stability iteration

    ! Energy balance check in canopy
    err = sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
         + cir*t_grnd**4 - eflx_sh_veg - hvap*qflx_evap_veg
    if (abs(err) > 0.1) then
       write(6,*) 'energy balance in canopy X',err
    endif

    ! Fluxes from ground to canopy space 
    delt  = wtal*t_grnd-wtl0*t_veg-wta0*thm
    delq  = wtalq*qg-wtlq0*qsatl-wtaq0*forc_q
    taux  = -forc_rho*forc_u/ram(1)
    tauy  = -forc_rho*forc_v/ram(1)
    eflx_sh_grnd = cpair*forc_rho*wtg*delt
    qflx_evap_soi = forc_rho*wtgq*delq

!CAS
!     Canopy air temperature
!CAS
    t_af = taf

    ! 2 m height air temperature
    t_ref2m = (taf + temp1*dth* 1./vkc *log((2.+z0hv)/z0hv))

    ! Downward longwave radiation below the canopy    
    dlrad = (1.-emv)*emg*forc_lwrad + emv*emg*sb*tlbef**3*(tlbef + 4.*dt_veg)

    ! Upward longwave radiation above the canopy    
    ulrad = ((1.-emg)*(1.-emv)*(1.-emv)*forc_lwrad &
         + emv*(1.+(1.-emg)*(1.-emv))*sb*tlbef**3*(tlbef + 4.*dt_veg) &
         + emg*(1.-emv)*sb*t_grnd**4)

    ! Derivative of soil energy flux with respect to soil temperature 
    cgrnds = cgrnds + cpair*forc_rho*wtg*wtal
    cgrndl = cgrndl + forc_rho*wtgq*wtalq*dqgdT
    cgrnd  = cgrnds  + cgrndl*htvp

    ! Update dew accumulation (kg/m2) 
    h2ocan = max(0._r8,h2ocan + (qflx_tran_veg-qflx_evap_veg)*dtime)

  end subroutine CanopyFluxes

  !=======================================================================

  subroutine Stomata (mpe,  apar,   ei,    ea,   tgcm,   &
       o2,   co2,    btran, rb,   rs,     &
       psn,  qe25,   vcmx25,mp,   c3psn,  &
       forc_pbot, tl, annpsnpot, annpsn) 

    !-----------------------------------------------------------------------
    ! Purpose:
    ! Leaf stomatal resistance and leaf photosynthesis.
    !
    ! Method:
    !
    ! Author:
    ! author:            Gordon Bonan
    ! standardized:      J. Truesdale, Feb. 1996
    ! reviewed:          G. Bonan, Feb. 1996
    ! 15 September 1999: Yongjiu Dai; Initial code
    ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
    !
    ! Revisions:
    ! 1/30/02, PET
    ! Made all input and output parameters explicit (removed reference to clm)
    !
    !-----------------------------------------------------------------------
    ! $Id$
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_RGAS
    implicit none

    !----Arguments----------------------------------------------------------

    real(r8), intent(in) :: mpe    ! prevents division by zero errors
    real(r8), intent(in) :: ei     ! vapor pressure inside leaf (sat vapor press at tl) (pa)
    real(r8), intent(in) :: ea     ! vapor pressure of canopy air (pa)
    real(r8), intent(in) :: apar   ! par absorbed per unit lai (w/m**2)
    real(r8), intent(in) :: o2     ! atmospheric o2 concentration (pa)
    real(r8), intent(in) :: co2    ! atmospheric co2 concentration (pa)
    real(r8), intent(in) :: tgcm   ! air temperature at agcm reference height (kelvin)
    real(r8), intent(in) :: btran  ! soil water transpiration factor (0 to 1)
    real(r8), intent(in) :: qe25   ! quantum efficiency at 25c (umol co2 / umol photon)
    real(r8), intent(in) :: vcmx25 ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real(r8), intent(in) :: mp     ! slope for conductance-to-photosynthesis relationship 
    real(r8), intent(in) :: c3psn  ! photosynthetic pathway: 0. = c4, 1. = c3
    real(r8), intent(in) :: forc_pbot !atmospheric pressure (Pa)
    real(r8), intent(in) :: tl     !leaf temperature (Kelvin)

    real(r8), intent(inout) :: rb  ! boundary layer resistance (s/m)
    real(r8), intent(inout) :: annpsnpot !annual potential photosynthesis (umol CO2 /m**2)
    real(r8), intent(inout) :: annpsn    !annual photosynthesis (umol CO2 /m**2)

    real(r8), intent(out)   :: rs  ! leaf stomatal resistance (s/m)
    real(r8), intent(out)   :: psn ! foliage photosynthesis (umol co2 /m**2/ s) [always +]

    !----Local Variables----------------------------------------------------

    integer, parameter :: niter = 3  ! number of iterations
    integer  iter                    ! iteration index

    real(r8) ab      ! used in statement functions
    real(r8) bc      ! used in statement functions
    real(r8) f1      ! generic temperature response (statement function)
    real(r8) f2      ! generic temperature inhibition (statement function)
    real(r8) tc      ! leaf temperature (degree celsius)
    real(r8) cs      ! co2 concentration at leaf surface (pa)
    real(r8) kc      ! co2 michaelis-menten constant (pa)
    real(r8) ko      ! o2 michaelis-menten constant (pa)
    real(r8) a,b,c,q ! intermediate calculations for rs
    real(r8) r1,r2   ! roots for rs
    real(r8) ppf     ! absorb photosynthetic photon flux (umol photons/m**2/s)
    real(r8) wc      ! rubisco limited photosynthesis (umol co2/m**2/s)
    real(r8) wj      ! light limited photosynthesis (umol co2/m**2/s)
    real(r8) we      ! export limited photosynthesis (umol co2/m**2/s)
    real(r8) cp      ! co2 compensation point (pa)
    real(r8) ci      ! internal co2 (pa)
    real(r8) awc     ! intermediate calcuation for wc
    real(r8) vcmx    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8) j       ! electron transport (umol co2/m**2/s)
    real(r8) cea     ! constrain ea or else model blows up
    real(r8) cf      ! s m**2/umol -> s/m
    real(r8) rsmax0  ! maximum stomatal resistance [s/m]

    real(r8) kc25    ! co2 michaelis-menten constant at 25c (pa)
    real(r8) akc     ! q10 for kc25
    real(r8) ko25    ! o2 michaelis-menten constant at 25c (pa)
    real(r8) ako     ! q10 for ko25
    real(r8) avcmx   ! q10 for vcmx25
    real(r8) bp      ! minimum leaf conductance (umol/m**2/s)

    !----End Variable List--------------------------------------------------

    f1(ab,bc) = ab**((bc-25.)/10.)
    f2(ab) = 1. + exp((-2.2e05+710.*(ab+SHR_CONST_TKFRZ))/(SHR_CONST_RGAS*0.001*(ab+SHR_CONST_TKFRZ)))

    kc25 = 30.
    akc = 2.1
    ko25 = 30000.
    ako = 1.2
    avcmx = 2.4
    bp = 2000.

    !
    ! Initialize rs=rsmax and psn=0 because calculations are performed only
    ! when apar > 0, in which case rs <= rsmax and psn >= 0
    ! Set constants
    !

    rsmax0 = 2.e4
    cf = forc_pbot/(SHR_CONST_RGAS*0.001*tgcm)*1.e06 

    if (apar <= 0.) then          ! night time
       rs = min(rsmax0, 1./bp * cf)
       psn = 0.
       return
    else                          ! day time
       tc = tl-SHR_CONST_TKFRZ
       ppf = 4.6*apar                  
       j = ppf*qe25
       kc = kc25 * f1(akc,tc)       
       ko = ko25 * f1(ako,tc)
       awc = kc * (1.+o2/ko)
       cp = 0.5*kc/ko*o2*0.21
       vcmx = vcmx25 * f1(avcmx,tc) / f2(tc) * btran

       !
       ! First guess ci
       !

       ci = 0.7*co2*c3psn + 0.4*co2*(1.-c3psn)  

       !
       ! rb: s/m -> s m**2 / umol
       !

       rb = rb/cf 

       !
       ! Constrain ea
       !

       cea = max(0.25*ei*c3psn+0.40*ei*(1.-c3psn), min(ea,ei) ) 

       !
       ! ci iteration for 'actual' photosynthesis
       !

       do iter = 1, niter
          wj = max(ci-cp,0._r8)*j/(ci+2.*cp)*c3psn + j*(1.-c3psn)
          wc = max(ci-cp,0._r8)*vcmx/(ci+awc)*c3psn + vcmx*(1.-c3psn)
          we = 0.5*vcmx*c3psn + 4000.*vcmx*ci/forc_pbot*(1.-c3psn) 
          psn = min(wj,wc,we) 
          cs = max( co2-1.37*rb*forc_pbot*psn, mpe )
          a = mp*psn*forc_pbot*cea / (cs*ei) + bp
          b = ( mp*psn*forc_pbot/cs + bp ) * rb - 1.
          c = -rb
          if (b >= 0.) then
             q = -0.5*( b + sqrt(b*b-4.*a*c) )
          else
             q = -0.5*( b - sqrt(b*b-4.*a*c) )
          endif
          r1 = q/a
          r2 = c/q
          rs = max(r1,r2)
          ci = max( cs-psn*forc_pbot*1.65*rs, 0._r8 )
       enddo

       ! rs, rb:  s m**2 / umol -> s/m 
       rs = min(rsmax0, rs*cf)
       rb = rb*cf 

    endif

  end subroutine Stomata

end module CanopyFluxesMod
