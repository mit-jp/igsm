#include <misc.h>
#include <preproc.h>

module FrictionVelocityMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: 
! 
! !DESCRIPTION: 
! Calculation of the friction velocity, relation for potential 
! temperature and humidity profiles of surface boundary layer. 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FrictionVelocity      ! Calculate friction velocity
  public :: MoninObukIni          ! Initialization of the Monin-Obukhov length
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: StabilityFunc        ! Stability function for rib < 0.
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
! !IROUTINE: FrictionVelocity
!
! !INTERFACE:
  subroutine FrictionVelocity (displa, z0m,   z0h,   z0q,   obu, &
                               iter, ur, um, ustar, temp1,  &
                               temp2, forc_hgt, forc_hgt_u, forc_hgt_t, forc_hgt_q, &
                               u10, fm, fv)
!
! !DESCRIPTION: 
! Calculation of the friction velocity, relation for potential 
! temperature and humidity profiles of surface boundary layer. 
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! !USES:
    use clm_varcon, only : vkc
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: displa     !displacement height [m]
    real(r8), intent(in) :: z0m        !roughness length, momentum [m]
    real(r8), intent(in) :: z0h        !roughness length, sensible heat [m]
    real(r8), intent(in) :: z0q        !roughness length, latent heat [m]
    real(r8), intent(in) :: obu        !monin-obukhov length (m)
    real(r8), intent(in) :: um         !wind speed including the stablity effect [m/s]
    real(r8), intent(in) :: ur   
    integer , intent(in) :: iter       !iteration number
    real(r8), intent(in) :: forc_hgt   !atmospheric reference height (m)
    real(r8), intent(in) :: forc_hgt_u !observational height of wind [m]
    real(r8), intent(in) :: forc_hgt_t !observational height of temperature [m]
    real(r8), intent(in) :: forc_hgt_q !observational height of humidity [m]
    real(r8), intent(out) :: ustar     !friction velocity [m/s]
    real(r8), intent(out) :: temp1     !relation for potential temperature profile
    real(r8), intent(out) :: temp2     !relation for specific humidity profile
    real(r8), intent(out) :: u10       !10-m wind (m/s) (for dust model)
    real(r8), intent(out) :: fv	       !friction velocity (m/s) (for dust model)
    real(r8), intent(inout) :: fm      !needed for BGC only to diagnose 10m wind speed
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 12/19/01, Peter Thornton
! Added arguments to eliminate passing clm derived type into this function.
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) zldis   ! reference height "minus" zero displacement heght [m]
    real(r8) zetam   ! transition point of flux-gradient relation (wind profile)
    real(r8) zetat   ! transition point of flux-gradient relation (temp. profile)
    real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory
!-----------------------------------------------------------------------

    ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.
    ! Wind profile

    zldis=forc_hgt_u-displa
    zeta=zldis/obu
    zetam=1.574

    if (zeta < -zetam) then           ! zeta < -1
       ustar=vkc*um/(log(-zetam*obu/z0m)- &
            StabilityFunc(1,-zetam) +StabilityFunc(1,z0m/obu) &
            +1.14*((-zeta)**0.333-(zetam)**0.333))
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       ustar=vkc*um/(log(zldis/z0m)- &
            StabilityFunc(1,zeta)+StabilityFunc(1,z0m/obu))
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       ustar=vkc*um/(log(zldis/z0m) + &
            5.*zeta -5.*z0m/obu)
    else                             !  1 < zeta, phi=5+zeta
       ustar=vkc*um/(log(obu/z0m)+5.-5.*z0m/obu &
            +(5.*log(zeta)+zeta-1.))
    endif

#if (defined PERGRO)
    if (zeta < -zetam) then           ! zeta < -1
       ustar=vkc*um/log(-zetam*obu/z0m)
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       ustar=vkc*um/log(zldis/z0m)
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       ustar=vkc*um/log(zldis/z0m)
    else                             !  1 < zeta, phi=5+zeta
       ustar=vkc*um/log(obu/z0m)
    endif
#endif

    ! Temperature profile

    zldis=forc_hgt_t-displa
    zeta=zldis/obu
    zetat=0.465
    if (zeta < -zetat) then           ! zeta < -1
       temp1=vkc/(log(-zetat*obu/z0h)-StabilityFunc(2,-zetat) &
            + StabilityFunc(2,z0h/obu) &
            + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       temp1=vkc/(log(zldis/z0h) - StabilityFunc(2,zeta) + &
            StabilityFunc(2,z0h/obu))
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       temp1=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
    else                             !  1 < zeta, phi=5+zeta
       temp1=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
            + (5.*log(zeta)+zeta-1.))
    endif

#if (defined PERGRO)
    if (zeta < -zetat) then           ! zeta < -1
       temp1=vkc/log(-zetat*obu/z0h)
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       temp1=vkc/log(zldis/z0h)
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       temp1=vkc/log(zldis/z0h)
    else                             !  1 < zeta, phi=5+zeta
       temp1=vkc/log(obu/z0h)
    endif
#endif

    ! Humidity profile

    zldis=forc_hgt_q-displa
    zeta=zldis/obu
    zetat=0.465
    if (zeta < -zetat) then          ! zeta < -1
       temp2=vkc/(log(-zetat*obu/z0q) - &
            StabilityFunc(2,-zetat) + StabilityFunc(2,z0q/obu) &
            + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       temp2=vkc/(log(zldis/z0q) - &
            StabilityFunc(2,zeta)+StabilityFunc(2,z0q/obu))
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       temp2=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
    else                             !  1 < zeta, phi=5+zeta
       temp2=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
            + (5.*log(zeta)+zeta-1.))
    endif

#if (defined PERGRO)
    if (zeta < -zetat) then          ! zeta < -1
       temp2=vkc/log(-zetat*obu/z0q)
    else if (zeta < 0.) then         ! -1 <= zeta < 0
       temp2=vkc/log(zldis/z0q)
    else if (zeta <= 1.) then        !  0 <= ztea <= 1
       temp2=vkc/log(zldis/z0q)
    else                             !  1 < zeta, phi=5+zeta
       temp2=vkc/log(obu/z0q)
    endif
#endif

  end subroutine FrictionVelocity

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc
!
! !INTERFACE:
  real(r8) function StabilityFunc(k, zeta)
!
! !DESCRIPTION: 
! Stability function for rib < 0.
!
! !USES:
    use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: k 
    real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! subroutine FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) chik     
!-----------------------------------------------------------------------

    chik = (1.-16.*zeta)**0.25
    if (k == 1) then
       StabilityFunc = 2.*log((1.+chik)*0.5) &
            + log((1.+chik*chik)*0.5)-2.*atan(chik)+SHR_CONST_PI*0.5
    else
       StabilityFunc = 2.*log((1.+chik*chik)*0.5)
    endif

  end function StabilityFunc

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MoninObukIni
!
! !INTERFACE:
  subroutine MoninObukIni (ur, thv, dthv, zldis, z0m, &
       um, obu  )
!
! !DESCRIPTION: 
! Initialization of the Monin-Obukhov length.
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! !USES:
    use clm_varcon, only : grav
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in) :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in) :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in) :: z0m   ! roughness length, momentum [m]

    real(r8), intent(out) :: um   ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu  ! monin-obukhov length (m)
!
! !CALLED FROM:
! subroutine BareGroundFluxes in module BareGroundFluxesMod.F90
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod.F90
! subroutine CanopyFluxes in module CanopyFluxesMod.F90
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8)  wc    ! convective velocity [m/s]
    real(r8)  rib   ! bulk Richardson number
    real(r8)  zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8)  ustar ! friction velocity [m/s]     
!-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06
    wc=0.5
    if (dthv >= 0.) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)
#if (defined PERGRO)
    rib = 0.
#endif

    if (rib >= 0.) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19_r8))
       zeta = min(2._r8,max(zeta,0.01_r8 ))
    else                    !unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityMod
