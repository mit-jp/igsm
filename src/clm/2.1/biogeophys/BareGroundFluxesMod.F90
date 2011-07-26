#include <misc.h>
#include <preproc.h>

module BareGroundFluxesMod 

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: 
! 
! !DESCRIPTION: 
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
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
! !IROUTINE: BareGroundFluxes
!
! !INTERFACE:
  subroutine BareGroundFluxes (c, p)
!
! !DESCRIPTION: 
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
    use clmtype
    use clm_varcon, only : cpair, vkc, grav
    use shr_const_mod, only : SHR_CONST_RGAS
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
! This routine originally had a long list of parameters, and also a reference to
! the entire clm derived type.  For consistency, only the derived type reference
! is passed (now pointing to the current column and pft), and the other original
! parameters are initialized locally. Using t_grnd instead of tg (tg eliminated 
! as redundant).
! 1/23/02, PET: Added pft reference as parameter. All outputs will be written
! to the pft data structures, and averaged to the column level outside of
! this routine.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: t_grnd       !ground surface temperature [K]
    real(r8), pointer :: thm          !intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8), pointer :: qg	          !specific humidity at ground surface [kg/kg]
    real(r8), pointer :: thv          !virtual potential temperature (kelvin)
    real(r8), pointer :: z0mg         !roughness length, momentum [m]
    real(r8), pointer :: dqgdT        !temperature derivative of "qg"
    real(r8), pointer :: htvp         !latent heat of evaporation (/sublimation) [J/kg]
    real(r8), pointer :: beta         !coefficient of conective velocity [-]
    real(r8), pointer :: zii          !convective boundary height [m]
    real(r8), pointer :: ur	          !wind speed at reference height [m/s]
    real(r8), pointer :: forc_u       !atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v       !atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_t       !atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_th      !atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q       !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rho     !density (kg/m**3)
    real(r8), pointer :: forc_pbot    !atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt     !atmospheric reference height (m)
    real(r8), pointer :: forc_hgt_u   !observational height of wind [m]
    real(r8), pointer :: forc_hgt_t   !observational height of temperature [m] 
    real(r8), pointer :: forc_hgt_q   !observational height of humidity [m]
!
! local pointers to implicit inout scalars
!
    real(r8), pointer :: z0hg         !roughness length, sensible heat [m]
    real(r8), pointer :: z0qg         !roughness length, latent heat [m]
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: dlrad        !downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad        !upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: cgrnds       !deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl       !deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: cgrnd        !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: taux         !wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy         !wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd !sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot  !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi!soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot!qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: t_ref2m      !2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_veg        !vegetation temperature (Kelvin)
    real(r8), pointer :: btran        !transpiration wetness factor (0 to 1)
    real(r8), pointer :: rssun        !sunlit stomatal resistance (s/m)
    real(r8), pointer :: rssha        !shaded stomatal resistance (s/m)
    real(r8), pointer :: u10          !10-m wind (m/s) (for dust model)
    real(r8), pointer :: fv           !friction velocity (m/s) (for dust model)
    real(r8), pointer :: ram1         !aerodynamical resistance (s/m)
!
! local pointers to implicit out arrays
!
    real(r8), dimension(:), pointer :: rootr !effective fraction of roots in each soil layer
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: nmozsgn !number of times moz changes sign
    integer  :: niters  !maximum number of iterations for surface temperature
    integer  :: iter    !iteration index
    real(r8) :: zldis   !reference height "minus" zero displacement height [m]
    real(r8) :: displa  !displacement height [m]
    real(r8) :: zeta    !dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc      !convective velocity [m/s]
    real(r8) :: dth     !diff of virtual temp. between ref. height and surface
    real(r8) :: dthv    !diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh     !diff of humidity between ref. height and surface
    real(r8) :: obu     !Monin-Obukhov length (m)
    real(r8) :: um      !wind speed including the stablity effect [m/s]
    real(r8) :: temp1   !relation for potential temperature profile
    real(r8) :: temp2   !relation for specific humidity profile
    real(r8) :: ustar   !friction velocity [m/s]
    real(r8) :: tstar   !temperature scaling parameter
    real(r8) :: qstar   !moisture scaling parameter
    real(r8) :: thvstar !virtual potential temperature scaling parameter
    real(r8) :: cf      !heat transfer coefficient from leaves [-]
    real(r8) :: ram     !aerodynamical resistance [s/m]
    real(r8) :: rah     !thermal resistance [s/m]
    real(r8) :: raw     !moisture resistance [s/m]
    real(r8) :: raih    !temporary variable [kg/m2/s]
    real(r8) :: raiw    !temporary variable [kg/m2/s]
    real(r8) :: obuold  !monin-obukhov length from previous iteration
    real(r8) :: fm     !needed for BGC only to diagnose 10m wind speed
    type(atm2lnd_state_type), pointer :: a2ls !local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps  !local pointers to derived subtypes 
    type(column_estate_type), pointer :: ces  !local pointers to derived subtypes
    type(column_wstate_type), pointer :: cws  !local pointers to derived subtypes
    type(pft_pstate_type)   , pointer :: pps  !local pointers to derived subtypes
    type(pft_estate_type)   , pointer :: pes  !local pointers to derived subtypes
    type(pft_mflux_type)    , pointer :: pmf  !local pointers to derived subtypes
    type(pft_eflux_type)    , pointer :: pef  !local pointers to derived subtypes
    type(pft_wflux_type)    , pointer :: pwf  !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    a2ls => c%a2ls
    cps => c%cps
    ces => c%ces
    cws => c%cws
    pps => p%pps
    pes => p%pes
    pmf => p%pmf
    pef => p%pef
    pwf => p%pwf

    ! assign local pointers to derived type scalar members
    forc_u => a2ls%forc_u
    forc_v => a2ls%forc_v
    forc_t => a2ls%forc_t
    forc_th => a2ls%forc_th
    forc_q => a2ls%forc_q
    forc_rho => a2ls%forc_rho
    forc_pbot => a2ls%forc_pbot
    forc_hgt => a2ls%forc_hgt
    forc_hgt_u => a2ls%forc_hgt_u
    forc_hgt_t => a2ls%forc_hgt_t
    forc_hgt_q => a2ls%forc_hgt_q
    z0mg => cps%z0mg
    z0hg => cps%z0hg
    z0qg => cps%z0qg
    htvp => cps%htvp
    beta => cps%beta
    zii => cps%zii
    ur => cps%ur
    t_grnd => ces%t_grnd
    thm => ces%thm
    thv => ces%thv
    qg => cws%qg
    dqgdT => cws%dqgdT
    btran => pps%btran
    rssun => pps%rssun
    rssha => pps%rssha
    u10 => pps%u10
    fv => pps%fv
    ram1 => pps%ram1
    t_ref2m => pes%t_ref2m
    t_veg => pes%t_veg
    dlrad => pef%dlrad
    ulrad => pef%ulrad
    cgrnds => pef%cgrnds
    cgrndl => pef%cgrndl
    cgrnd => pef%cgrnd
    eflx_sh_grnd => pef%eflx_sh_grnd
    eflx_sh_tot => pef%eflx_sh_tot
    taux => pmf%taux
    tauy => pmf%tauy
    qflx_evap_soi => pwf%qflx_evap_soi
    qflx_evap_tot => pwf%qflx_evap_tot

    ! assign local pointers to derived type array members
    rootr => pps%rootr

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.

    ! Initialization variables
    dlrad  = 0.
    ulrad  = 0.
    nmozsgn = 0
    obuold = 0.
    dth   = thm-t_grnd
    dqh   = forc_q-qg
    dthv  = dth*(1.+0.61*forc_q)+0.61*forc_th*dqh
    zldis = forc_hgt_u-0.

    ! Initialize Monin-Obukhov length and wind speed
    call MoninObukIni(ur, thv, dthv, zldis, z0mg, um, obu  )

    ! Begin stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer
    niters=3
    do iter = 1, niters

       displa = 0.0_r8
       call FrictionVelocity(displa, z0mg, z0hg, z0qg, obu, &
            iter, ur, um, ustar, temp1, &
            temp2, forc_hgt, forc_hgt_u, forc_hgt_t, forc_hgt_q,  &
            u10, fm, fv)

       tstar = temp1*dth
       qstar = temp2*dqh
       z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
       z0qg = z0hg

       thvstar=tstar*(1.+0.61*forc_q) + 0.61*forc_th*qstar
       zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)
       if (zeta >= 0.) then	 !stable
          zeta = min(2._r8,max(zeta,0.01_r8))
          um = max(ur,0.1_r8)
       else					 !unstable
          zeta = max(-100._r8,min(zeta,-0.01_r8))
          wc = beta*(-grav*ustar*thvstar*zii/thv)**0.333
          um = sqrt(ur*ur+wc*wc)
       endif
       obu = zldis/zeta

       if (obuold*obu < 0.) nmozsgn = nmozsgn+1
       if (nmozsgn >= 4) EXIT

       obuold = obu

    enddo						! end stability iteration

    ! Determine aerodynamic resistances
    ram    = 1./(ustar*ustar/um)
    rah    = 1./(temp1*ustar)
    raw    = 1./(temp2*ustar)
    raih   = forc_rho*cpair/rah
    raiw   = forc_rho/raw
    ram1   = ram  !pass value to global variable

    ! Output to pft-level data structures
    ! Derivative of fluxes with respect to ground temperature
    cgrnds = raih
    cgrndl = raiw*dqgdT
    cgrnd  = cgrnds + htvp*cgrndl

    ! Surface fluxes of momentum, sensible and latent heat
    ! using ground temperatures from previous time step
    taux           = -forc_rho*forc_u/ram
    tauy           = -forc_rho*forc_v/ram
    eflx_sh_grnd   = -raih*dth
    eflx_sh_tot    = eflx_sh_grnd
    qflx_evap_soi  = -raiw*dqh
    qflx_evap_tot  = qflx_evap_soi

    ! 2 m height air temperature
    t_ref2m = (t_grnd+temp1*dth * 1./vkc *log((2.+z0hg)/z0hg))

    ! Variables needed by history tape
    t_veg = forc_t
    btran = 0.
    rootr(:) = 0.
    cf = forc_pbot/(SHR_CONST_RGAS*0.001*thm)*1.e06
    rssun = 1./1.e15 * cf
    rssha = 1./1.e15 * cf

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
