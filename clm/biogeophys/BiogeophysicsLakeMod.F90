#include <misc.h>
#include <preproc.h>

module BiogeophysicsLakeMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: BiogeophysicsLakeMod
! 
! !DESCRIPTION: 
! Calculates lake temperatures  and surface fluxes.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysicsLake 

! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BiogeophysicsLake 
!
! !INTERFACE:
  subroutine BiogeophysicsLake (c) 
!
! !DESCRIPTION: 
! Calculates lake temperatures  and surface fluxes.
! Lake temperatures are determined from a one-dimensional thermal
! stratification model based on eddy diffusion concepts to 
! represent vertical mixing of heat.
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz   
!
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
! There are two types of lakes: 
!   Deep lakes are 50 m. 
!   Shallow lakes are 10 m deep.
!
!   For unfrozen deep lakes:    ke > 0 and    convective mixing
!   For unfrozen shallow lakes: ke = 0 and no convective mixing
!
! Use the Crank-Nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
!
! The solution conserves energy as:
!
! cw*([ts(      1)] n+1 - [ts(      1)] n)*dz(      1)/dt + ... +
! cw*([ts(nlevlak)] n+1 - [ts(nlevlak)] n)*dz(nlevlak)/dt = fin
!
! where:
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot 
!            - hm + phi(1) + ... + phi(nlevlak) 
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals
    use clm_varpar, only : nlevlak
    use clm_varcon, only : hvap, hsub, hfus, cpair, cpliq, tkwat, tkice, &
         sb, vkc, grav, denh2o, tfrz, spval
    use SurfaceRadiationMod, only : SurfaceRadiation
    use QSatMod, only : QSat
    use FrictionVelocityMod, only : MoninObukIni, FrictionVelocity
    use TridiagonalMod, only : Tridiagonal
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout) :: c  !column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! Migrated to clm2.1 new data structures by Peter Thornton and M. Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: forc_t          !atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_pbot       !atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt        !atmospheric reference height (m)
    real(r8), pointer :: forc_hgt_q      !observational height of humidity [m]
    real(r8), pointer :: forc_hgt_t      !observational height of temperature [m]
    real(r8), pointer :: forc_hgt_u      !observational height of wind [m]
    real(r8), pointer :: forc_th         !atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q          !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_u          !atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v          !atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_lwrad      !downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: forc_rho        !density (kg/m**3)
    real(r8), pointer :: forc_snow       !snow rate [mm/s]
    real(r8), pointer :: forc_rain       !rain rate [mm/s]
    real(r8), pointer :: t_grnd          !ground temperature (Kelvin)
    real(r8), pointer :: h2osno          !snow water (mm H2O)
    real(r8), pointer :: snowdp          !snow height (m)
    real(r8), pointer :: sabg            !solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: lat             !latitude (radians)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: begwb           !water mass begining of the time step
    real(r8), pointer :: qflx_prec_grnd  !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_evap_soi   !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot   !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: eflx_sh_grnd    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lwrad_out  !emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net  !net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_soil_grnd  !soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_sh_tot     !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot     !total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_lh_grnd    !ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: t_rad_column    !radiative temperature (Kelvin)
    real(r8), pointer :: t_rad_pft       !radiative temperature (Kelvin)
    real(r8), pointer :: t_veg           !vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m         !2 m height surface air temperature (Kelvin)
    real(r8), pointer :: taux            !wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy            !wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: qmelt           !snow melt [mm/s]
    real(r8), pointer :: u10             !10-m wind (m/s) (for dust model)
    real(r8), pointer :: fv              !friction velocity (m/s) (for dust model)
    real(r8), pointer :: ram1            !aerodynamical resistance (s/m)
    real(r8), pointer :: errsoi          !soil/lake energy conservation error (W/m**2)
!
! local pointers to implicit in arrays
!
    real(r8), dimension(:), pointer :: dz !layer thickness (m)
    real(r8), dimension(:), pointer :: z  !layer depth (m)
!
! local pointers to implicit out arrays
!
    real(r8), dimension(:), pointer :: t_lake !lake temperature (Kelvin)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: i,j       !do loop or array index
    integer  :: idlak = 1 !index of lake, 1 = deep lake, 2 = shallow lake
    integer  :: niters    !maximum number of iterations for surface temperature
    integer  :: iter      !iteration index
    integer  :: nmozsgn   !number of times moz changes sign
    real(r8) :: ax       !
    real(r8) :: bx       !
    real(r8) :: beta1    !coefficient of conective velocity [-]
    real(r8) :: degdT    !d(eg)/dT
    real(r8) :: dqh      !diff of humidity between ref. height and surface
    real(r8) :: dth      !diff of virtual temp. between ref. height and surface
    real(r8) :: dthv     !diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dzsur    !
    real(r8) :: eg       !water vapor pressure at temperature T [pa]
    real(r8) :: emg      !ground emissivity (0.97 for snow,
    real(r8) :: hm       !energy residual [W/m2]
    real(r8) :: htvp     !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8) :: obu      !monin-obukhov length (m)
    real(r8) :: obuold   !monin-obukhov length of previous iteration
    real(r8) :: qsatg    !saturated humidity [kg/kg]
    real(r8) :: qsatgdT  !d(qsatg)/dT
    real(r8) :: qstar    !moisture scaling parameter
    real(r8) :: ram      !aerodynamical resistance [s/m]
    real(r8) :: rah      !thermal resistance [s/m]
    real(r8) :: raw      !moisture resistance [s/m]
    real(r8) :: stftg3   !
    real(r8) :: temp1    !relation for potential temperature profile
    real(r8) :: temp2    !relation for specific humidity profile
    real(r8) :: tgbef    !
    real(r8) :: thm      !intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv      !virtual potential temperature (kelvin)
    real(r8) :: thvstar  !virtual potential temperature scaling parameter
    real(r8) :: tksur    !thermal conductivity of snow/soil (w/m/kelvin)
    real(r8) :: tstar    !temperature scaling parameter
    real(r8) :: um       !wind speed including the stablity effect [m/s]
    real(r8) :: ur       !wind speed at reference height [m/s]
    real(r8) :: ustar    !friction velocity [m/s]
    real(r8) :: wc       !convective velocity [m/s]
    real(r8) :: zeta     !dimensionless height used in Monin-Obukhov theory
    real(r8) :: zii      !convective boundary height [m]
    real(r8) :: zldis    !reference height "minus" zero displacement height [m]
    real(r8) :: displa   !displacement height [m]
    real(r8) :: z0mg     !roughness length over ground, momentum [m]
    real(r8) :: z0hg     !roughness length over ground, sensible heat [m]
    real(r8) :: z0qg     !roughness length over ground, latent heat [m]
    real(r8) :: beta(2)  !fraction solar rad absorbed at surface: depends on lake type
    real(r8) :: za(2)    !base of surface absorption layer (m): depends on lake type
    real(r8) :: eta(2)   !light extinction coefficient (/m): depends on lake type
    real(r8) :: p0       !neutral value of turbulent prandtl number
    real(r8) :: a(nlevlak)    !"a" vector for tridiagonal matrix
    real(r8) :: b(nlevlak)    !"b" vector for tridiagonal matrix
    real(r8) :: c1(nlevlak)   !"c" vector for tridiagonal matrix
    real(r8) :: r(nlevlak)    !"r" vector for tridiagonal solution
    real(r8) :: rhow(nlevlak) !density of water (kg/m**3)
    real(r8) :: phi(nlevlak)  !solar radiation absorbed by layer (w/m**2)
    real(r8) :: kme(nlevlak)  !molecular + eddy diffusion coefficient (m**2/s)
    real(r8) :: cwat     !specific heat capacity of water (j/m**3/kelvin)
    real(r8) :: ws       !surface friction velocity (m/s)
    real(r8) :: ks       !coefficient
    real(r8) :: in       !relative flux of solar radiation into layer
    real(r8) :: out      !relative flux of solar radiation out of layer
    real(r8) :: ri       !richardson number
    real(r8) :: fin      !heat flux into lake - flux out of lake (w/m**2)
    real(r8) :: ocvts    !(cwat*(t_lake[n  ])*dz
    real(r8) :: ncvts    !(cwat*(t_lake[n+1])*dz
    real(r8) :: m1       !intermediate variable for calculating r, a, b, c
    real(r8) :: m2       !intermediate variable for calculating r, a, b, c
    real(r8) :: m3       !intermediate variable for calculating r, a, b, c
    real(r8) :: ke       !eddy diffusion coefficient (m**2/s)
    real(r8) :: km       !molecular diffusion coefficient (m**2/s)
    real(r8) :: zin      !depth at top of layer (m)
    real(r8) :: zout     !depth at bottom of layer (m)
    real(r8) :: drhodz   !d [rhow] /dz (kg/m**4)
    real(r8) :: n2       !brunt-vaisala frequency (/s**2)
    real(r8) :: num      !used in calculating ri
    real(r8) :: den      !used in calculating ri
    real(r8) :: tav      !used in aver temp for convectively mixed layers
    real(r8) :: nav      !used in aver temp for convectively mixed layers
    real(r8) :: phidum   !temporary value of phi
    real(r8) :: u2m      !2 m wind speed (m/s)
    real(r8) :: cf       !s m**2/umol -> s/m
    real(r8) :: fm       !needed for BGC only to diagnose 10m wind speed
    real(r8) :: t_rad 
    type(atm2lnd_state_type), pointer :: a2ls !local pointers to derived subtypes
    type(atm2lnd_flux_type) , pointer :: a2lf !local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps  !local pointers to derived subtypes
    type(column_estate_type), pointer :: ces  !local pointers to derived subtypes
    type(column_wstate_type), pointer :: cws  !local pointers to derived subtypes
    type(column_eflux_type) , pointer :: cef  !local pointers to derived subtypes
    type(column_wflux_type) , pointer :: cwf  !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    a2ls => c%a2ls
    a2lf => c%a2lf
    cps => c%cps
    ces => c%ces
    cws => c%cws
    cef => c%cef
    cwf => c%cwf

    ! assign local pointers to derived type array members
    dz => cps%dz
    z => cps%z
    t_lake => ces%t_lake

    ! assign local pointers to derived type scalar members
    forc_t => a2ls%forc_t
    forc_pbot => a2ls%forc_pbot
    forc_hgt => a2ls%forc_hgt
    forc_hgt_u => a2ls%forc_hgt_u
    forc_hgt_t => a2ls%forc_hgt_t
    forc_hgt_q => a2ls%forc_hgt_q
    forc_th => a2ls%forc_th
    forc_q => a2ls%forc_q
    forc_u => a2ls%forc_u
    forc_v => a2ls%forc_v
    forc_rho => a2ls%forc_rho
    forc_lwrad  => a2lf%forc_lwrad
    forc_snow => a2lf%forc_snow
    forc_rain => a2lf%forc_rain
    h2osno => cws%h2osno
    snowdp => cps%snowdp
    lat => cps%lps%gps%lat
    sabg => c%p(1)%pef%sabg    
    t_grnd => ces%t_grnd
    t_ref2m => c%p(1)%pes%t_ref2m
    t_veg => c%p(1)%pes%t_veg
    t_rad_column => ces%t_rad_column
    t_rad_pft => c%p(1)%pes%t_rad_pft
    begwb => c%cwbal%begwb
    errsoi => c%cebal%errsoi
    eflx_lwrad_out => c%p(1)%pef%eflx_lwrad_out
    eflx_lwrad_net => c%p(1)%pef%eflx_lwrad_net
    eflx_soil_grnd => c%p(1)%pef%eflx_soil_grnd
    eflx_lh_tot => c%p(1)%pef%eflx_lh_tot
    eflx_lh_grnd => c%p(1)%pef%eflx_lh_grnd
    eflx_sh_grnd => c%p(1)%pef%eflx_sh_grnd
    eflx_sh_tot => c%p(1)%pef%eflx_sh_tot
    qmelt => cwf%qmelt
    u10 => c%p(1)%pps%u10
    fv => c%p(1)%pps%fv
    ram1 => c%p(1)%pps%ram1
    taux => c%p(1)%pmf%taux
    tauy => c%p(1)%pmf%tauy
    qflx_prec_grnd => c%p(1)%pwf%qflx_prec_grnd
    qflx_evap_soi => c%p(1)%pwf%qflx_evap_soi
    qflx_evap_tot => c%p(1)%pwf%qflx_evap_tot

    ! Surface Radiation 

    call SurfaceRadiation(c)

    ! Determine beginning water balance

    begwb = h2osno  

    ! ------------------------------------------------------- 
    ! Constants and model parameters
    ! ------------------------------------------------------- 

    ! Constants for lake temperature model

    beta = (/0.4, 0.4/)                            ! (deep lake, shallow lake)
    za   = (/0.6, 0.5/)
    eta  = (/0.1, 0.5/)
    p0   = 1.

    ! Roughness lengths

    if (t_grnd >= tfrz)then                    ! for unfrozen lake
       z0mg = 0.01
    else                                           ! for frozen lake
       z0mg = 0.04
    endif
    z0hg = z0mg
    z0qg = z0mg

    ! Latent heat 

    if (forc_t > tfrz) then
       htvp = hvap
    else
       htvp = hsub
    endif

#if (defined PERGRO)
    htvp = hvap
#endif

    ! Emissivity

    emg = 0.97

    ! ------------------------------------------------------- 
    ! Surface temperature and fluxes
    ! ------------------------------------------------------- 

    dzsur = dz(1) + snowdp

    ! Saturated vapor pressure, specific humidity and their derivatives
    ! at lake surface

    call QSat(t_grnd, forc_pbot, eg, degdT, qsatg, qsatgdT)

    ! Potential, virtual potential temperature, and wind speed at the
    ! reference height

    beta1=1.       ! -  (in computing W_*)
    zii = 1000.    ! m  (pbl height)
    thm = forc_t + 0.0098*forc_hgt_t         ! intermediate variable
    thv = forc_th*(1.+0.61*forc_q)           ! virtual potential T
    ur = max(1.0_r8,sqrt(forc_u*forc_u+forc_v*forc_v))

    ! Initialize stability variables

    nmozsgn = 0
    obuold = 0.

    dth   = thm-t_grnd
    dqh   = forc_q-qsatg
    dthv  = dth*(1.+0.61*forc_q)+0.61*forc_th*dqh
    zldis = forc_hgt_u-0.

    ! Initialize Monin-Obukhov length and wind speed

    call MoninObukIni(ur, thv, dthv, zldis, z0mg, um, obu)

    ! Begin stability iteration

    niters = 3
    do iter = 1, niters
       tgbef = t_grnd
       if (t_grnd > tfrz) then
          tksur = tkwat
       else
          tksur = tkice
       endif

       ! Determine friction velocity, and potential temperature and humidity
       ! profiles of the surface boundary layer

       displa = 0.0_r8
       call FrictionVelocity (displa, z0mg,  z0hg,  z0qg,  obu, &
            iter, ur, um, ustar, temp1, &
            temp2, forc_hgt, forc_hgt_u, forc_hgt_t, forc_hgt_q, &
            u10, fm, fv)
       obuold = obu

       ! Determine aerodynamic resistances

       ram    = 1./(ustar*ustar/um)
       rah    = 1./(temp1*ustar)
       raw    = 1./(temp2*ustar)
       ram1   = ram !pass value to global variable

       ! Get derivative of fluxes with respect to ground temperature

       stftg3 = emg*sb*tgbef*tgbef*tgbef

       ax  = sabg + emg*forc_lwrad + 3.*stftg3*tgbef &
            + forc_rho*cpair/rah*thm &
            - htvp*forc_rho/raw*(qsatg-qsatgdT*tgbef - forc_q) &
            + tksur*t_lake(1)/dzsur

       bx  = 4.*stftg3 + forc_rho*cpair/rah &
            + htvp*forc_rho/raw*qsatgdT + tksur/dzsur

       t_grnd = ax/bx

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       eflx_sh_grnd = forc_rho*cpair*(t_grnd-thm)/rah
       qflx_evap_soi = forc_rho*(qsatg+qsatgdT*(t_grnd-tgbef)-forc_q)/raw

       ! Re-calculate saturated vapor pressure, specific humidity and their
       ! derivatives at lake surface

       call QSat(t_grnd, forc_pbot, eg, degdT, qsatg, qsatgdT)

       dth=thm-t_grnd
       dqh=forc_q-qsatg

       tstar = temp1*dth
       qstar = temp2*dqh

       dthv=dth*(1.+0.61*forc_q)+0.61*forc_th*dqh
       thvstar=tstar*(1.+0.61*forc_q) + 0.61*forc_th*qstar
       zeta=zldis*vkc * grav*thvstar/(ustar**2*thv)

       if (zeta >= 0.) then     !stable
          zeta = min(2._r8,max(zeta,0.01_r8))
          um = max(ur,0.1_r8)
       else                     !unstable
          zeta = max(-100._r8,min(zeta,-0.01_r8))
          wc = beta1*(-grav*ustar*thvstar*zii/thv)**0.333
          um = sqrt(ur*ur+wc*wc)
       endif
       obu = zldis/zeta

       if (obuold*obu < 0.) nmozsgn = nmozsgn+1
       if (nmozsgn >= 4) EXIT

    enddo

    ! If there is snow on the ground and t_grnd > tfrz: reset t_grnd = tfrz.
    ! Re-evaluate ground fluxes. Energy imbalance used to melt snow.  
    ! h2osno > 0.5 prevents spurious fluxes.
    ! note that qsatg and qsatgdT should be f(tgbef) (PET: not sure what this
    ! comment means)

    if (h2osno > 0.5 .AND. t_grnd > tfrz) then
       t_grnd = tfrz
       eflx_sh_grnd = forc_rho*cpair*(t_grnd-thm)/rah
       qflx_evap_soi = forc_rho*(qsatg+qsatgdT*(t_grnd-tgbef) - forc_q)/raw  
    endif

    ! Net longwave from ground to atmosphere

    eflx_lwrad_out = (1.-emg)*forc_lwrad + stftg3*(-3.*tgbef+4.*t_grnd)

    ! Radiative temperature

    t_rad = (eflx_lwrad_out/sb)**0.25
    t_rad_column = t_rad
    t_rad_pft = t_rad

    ! Ground heat flux

    eflx_soil_grnd = sabg + forc_lwrad - eflx_lwrad_out - &
         eflx_sh_grnd - htvp*qflx_evap_soi

    taux   = -forc_rho*forc_u/ram
    tauy   = -forc_rho*forc_v/ram

    eflx_sh_tot   = eflx_sh_grnd
    qflx_evap_tot = qflx_evap_soi
    eflx_lh_tot   = htvp*qflx_evap_soi
    eflx_lh_grnd  = htvp*qflx_evap_soi

    ! 2 m height air temperature
    t_ref2m   = (t_grnd + temp1*dth * 1./vkc *log((2.+z0hg)/z0hg))

    ! Energy residual used for melting snow
    if (h2osno > 0. .AND. t_grnd >= tfrz) then
       hm = min(h2osno*hfus/dtime, max(eflx_soil_grnd,0._r8))
    else
       hm = 0.
    endif
    qmelt = hm/hfus             ! snow melt (mm/s)

    ! Lake layer temperature

    ! Lake density

    do j = 1, nlevlak
       rhow(j) = 1000.*( 1.0 - 1.9549e-05*(abs(t_lake(j)-277.))**1.68 )
    enddo

    ! Eddy diffusion +  molecular diffusion coefficient:
    ! eddy diffusion coefficient used for unfrozen deep lakes only

    cwat = cpliq*denh2o
    km = tkwat/cwat

    fin = beta(idlak) * sabg + forc_lwrad - (eflx_lwrad_out + &
         eflx_sh_tot + eflx_lh_tot + hm)
    u2m = max(1.0_r8,ustar/vkc*log(2./z0mg))

    ws = 1.2e-03 * u2m
    ks = 6.6*sqrt(abs(sin(lat)))*(u2m**(-1.84))

    do j = 1, nlevlak-1
       drhodz = (rhow(j+1)-rhow(j)) / (z(j+1)-z(j))
       n2 = -grav / rhow(j) * drhodz
       num = 40. * n2 * (vkc*z(j))**2
       den = max( (ws**2) * exp(-2.*ks*z(j)), 1.e-10_r8 )
       ri = ( -1. + sqrt( max(1.+num/den, 0._r8) ) ) / 20.
       if (idlak == 1 .AND. t_grnd > tfrz) then
          ke = vkc*ws*z(j)/p0 * exp(-ks*z(j)) / (1.+37.*ri*ri)
       else
          ke = 0.
       endif
       kme(j) = km + ke
    enddo

    kme(nlevlak) = kme(nlevlak-1)

    ! Heat source term: unfrozen lakes only

    do j = 1, nlevlak
       zin  = z(j) - 0.5*dz(j)
       zout = z(j) + 0.5*dz(j)
       in  = exp( -eta(idlak)*max(  zin-za(idlak),0._r8 ) )
       out = exp( -eta(idlak)*max( zout-za(idlak),0._r8 ) )

       ! Assume solar absorption is only in the considered depth
       if (j == nlevlak) out = 0.  
       if (t_grnd > tfrz) then
          phidum = (in-out) * sabg * (1.-beta(idlak))
       else if (j == 1) then
          phidum= sabg * (1.-beta(idlak))
       else
          phidum = 0.
       endif
       phi(j) = phidum
    enddo

    ! Sum cwat*t_lake*dz for energy check

    ocvts = 0.
    do j = 1, nlevlak
       ocvts = ocvts + cwat*t_lake(j)*dz(j)
    enddo

    ! Set up vector r and vectors a, b, c that define tridiagonal matrix

    j = 1
    m2 = dz(j)/kme(j) + dz(j+1)/kme(j+1)
    m3 = dtime/dz(j)
    r(j) = t_lake(j) + (fin+phi(j))*m3/cwat - (t_lake(j)-t_lake(j+1))*m3/m2
    a(j) = 0.
    b(j) = 1. + m3/m2
    c1(j) = -m3/m2

    j = nlevlak
    m1 = dz(j-1)/kme(j-1) + dz(j)/kme(j)
    m3 = dtime/dz(j)
    r(j) = t_lake(j) + phi(j)*m3/cwat + (t_lake(j-1)-t_lake(j))*m3/m1
    a(j) = -m3/m1
    b(j) = 1. + m3/m1
    c1(j) = 0.

    do j = 2, nlevlak-1
       m1 = dz(j-1)/kme(j-1) + dz(j  )/kme(j  )
       m2 = dz(j  )/kme(j  ) + dz(j+1)/kme(j+1)
       m3 = dtime/dz(j)
       r(j) = t_lake(j) + phi(j)*m3/cwat + &
            (t_lake(j-1) - t_lake(j))*m3/m1 - &
            (t_lake(j)-t_lake(j+1))*m3/m2

       a(j) = -m3/m1
       b(j) = 1. + m3/m1 + m3/m2
       c1(j) = -m3/m2
    enddo

    ! Solve for t_lake: a, b, c, r, u 

    call Tridiagonal (nlevlak, a, b, c1, r, t_lake(1:nlevlak)) 

    ! Convective mixing: make sure cwat*dz*ts is conserved.

    if (idlak == 1 .AND. t_grnd > tfrz) then
       do j = 1, nlevlak-1
          if (rhow(j) > rhow(j+1)) then
             tav = 0.
             nav = 0.
             do i = 1, j+1
                tav = tav + t_lake(i)*dz(i)
                nav = nav + dz(i)
             enddo
             tav = tav/nav
             do i = 1, j+1
                t_lake(i) = tav
                rhow(i) = 1000.*( 1.0 - 1.9549e-05*(abs(t_lake(i)-277.))**1.68 )
             enddo
          endif
       enddo
    endif

    ! Sum cwat*t_lake*dz and total energy into lake for energy check

    ncvts = 0.
    do j = 1, nlevlak
       ncvts = ncvts + cwat*t_lake(j)*dz(j)
       fin = fin + phi(j)
    enddo

    errsoi = (ncvts-ocvts) / dtime - fin

    ! The following are needed for global average on history tape.

    t_veg = forc_t  
    eflx_lwrad_net  = eflx_lwrad_out - forc_lwrad
    qflx_prec_grnd = forc_rain + forc_snow

    return
  end subroutine BiogeophysicsLake

end module BiogeophysicsLakeMod
