#include <misc.h>
#include <preproc.h>

module Hydrology1Mod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE:  Hydrology1Mod
! 
! !DESCRIPTION: 
! Calculation of 
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring. 
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydrology1
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
! !IROUTINE: Hydrology1 
!
! !INTERFACE:
  subroutine Hydrology1 (c)
!
! !DESCRIPTION: 
! Calculation of 
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring. 
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
! Note:  The evaporation loss is taken off after the calculation of leaf 
! temperature in the subroutine clm_leaftem.f90, not in this subroutine.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : tfrz, istice, istwet, istsoil
    use FracWetMod, only : FracWet
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
! 2/15/02, Peter Thornton: Migrated to new data structures. Required
! adding a PFT loop.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: itypwat      !water type
    integer , pointer :: itypprc      !precipitatin type
    real(r8), pointer :: forc_rain    !rain rate [mm/s]
    real(r8), pointer :: forc_snow    !snow rate [mm/s]
    real(r8), pointer :: forc_t       !atmospheric temperature (Kelvin)
    logical , pointer :: do_capsnow   !true => do snow capping
    real(r8), pointer :: t_grnd       !ground temperature (Kelvin)
    real(r8), pointer :: dewmx        !Maximum allowed dew [mm]
    integer , pointer :: frac_veg_nosno !fraction of veg not covered by snow (0/1 now) [-]
    real(r8), pointer :: elai         !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai         !one-sided stem area index with burying by snow
#if (defined PCP2PFT)
!CAS  For redistribution of precipitation across PFT types
!CAS  
!CAS  NOTE: This modification has been employed for use with the zonal configuration
!CAS        of CLM coupled to the MIT Integrated Global System Model (IGSM)
!CAS
    real(r8), pointer :: pcp2pft      !PCP Distribution factor across PFTs
#endif
!
! local pointers to original implicit inout scalars
!
    integer , pointer :: snl          !number of snow layers
    real(r8), pointer :: snowage      !non dimensional snow age [-]
    real(r8), pointer :: snowdp       !snow height (m)
    real(r8), pointer :: h2osno       !snow water (mm H2O)
    real(r8), pointer :: h2ocan       !canopy water (mm H2O)
!
! local pointers to original implicit out scalars
!
    real(r8), pointer :: qflx_prec_intr  !interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd  !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snowcap    !excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snow_grnd  !snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_rain_grnd  !rain on ground after interception (mm H2O/s) [+]
!
! local pointers to original implicit out arrays
!
    real(r8), dimension(:), pointer :: zi           !interface level below a "z" level (m)
    real(r8), dimension(:), pointer :: dz           !layer depth (m)
    real(r8), dimension(:), pointer :: z            !layer thickness (m)
    real(r8), dimension(:), pointer :: t_soisno     !soil temperature (Kelvin)
    real(r8), dimension(:), pointer :: h2osoi_ice   !ice lens (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_liq   !liquid water (kg/m2)
    real(r8), dimension(:), pointer :: frac_iceold  !fraction of ice relative to the tot water
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: newnode      !flag when new snow node is set, (1=yes, 0=no)
    real(r8):: prcp         !precipitation rate [mm/s]
    real(r8):: h2ocanmx     !maximum allowed water on canopy [mm]
    real(r8):: fpi          !coefficient of interception
    real(r8):: xrun         !excess water that exceeds the leaf capacity [mm/s]
    real(r8):: qflx_candrip !rate of canopy runoff and snow falling off canopy [mm/s]
    real(r8):: qflx_through !direct throughfall [mm/s]
    real(r8):: dz_snowf     !layer thickness rate change due to precipitation [mm/s]
    real(r8):: flfall       !fraction of liquid water within falling precip.
    real(r8):: bifall       !bulk density of newly fallen dry snow [kg/m3]
    real(r8):: coefb        !Slope of "Alta" expression for dependence of flfall on temp
    real(r8):: coefa        !Offset of  of "Alta" expression for dependence of flfall on temp
    integer :: pi	    !pft index
    real(r8):: pftsum       !temporary used for pft averaging for columns
    type(atm2lnd_flux_type)    , pointer :: a2lf ! local pointers to derived subtypes
    type(atm2lnd_state_type)   , pointer :: a2ls ! local pointers to derived subtypes 
    type(landunit_pstate_type) , pointer :: lps  ! local pointers to derived subtypes
    type(column_pstate_type)   , pointer :: cps  ! local pointers to derived subtypes
    type(column_wflux_type)    , pointer :: cwf  ! local pointers to derived subtypes
    type(column_wstate_type)   , pointer :: cws  ! local pointers to derived subtypes
    type(column_estate_type)   , pointer :: ces  ! local pointers to derived subtypes
    type(pft_type)             , pointer :: p    ! local pointers to derived subtypes
    type(pft_pstate_type)      , pointer :: pps  ! local pointers to derived subtypes
    type(pft_wstate_type)      , pointer :: pws  ! local pointers to derived subtypes
    type(pft_wflux_type)       , pointer :: pwf	 ! local pointers to derived subtypes 
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    cps => c%cps
    lps => cps%lps
    a2lf => c%a2lf
    a2ls => c%a2ls
    cwf => c%cwf
    cws => c%cws
    ces => c%ces

    ! Assign local pointers to derived type scalar members (column-level)
    itypwat => lps%itypwat
    forc_rain => a2lf%forc_rain
    forc_snow => a2lf%forc_snow
    itypprc => a2ls%itypprc
    forc_t => a2ls%forc_t
    do_capsnow => cps%do_capsnow
    t_grnd => ces%t_grnd
    snl => cps%snl
    snowdp => cps%snowdp
    snowage => cps%snowage
    h2osno => cws%h2osno

    ! Assign local pointers to derived type array members (column-level)
    zi => cps%zi
    dz => cps%dz
    z => cps%z
    frac_iceold => cps%frac_iceold
    t_soisno => ces%t_soisno
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! Start pft loop
    do pi=1,cps%npfts

       ! Anssign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pws => p%pws
       pwf => p%pwf

       ! Assign local pointers to derived type scalar members (pft-level)
       dewmx => pps%dewmx
       frac_veg_nosno => pps%frac_veg_nosno
       elai => pps%elai
       esai => pps%esai
       h2ocan => pws%h2ocan
       qflx_prec_intr => pwf%qflx_prec_intr
       qflx_prec_grnd => pwf%qflx_prec_grnd
       qflx_snowcap => pwf%qflx_snowcap
       qflx_snow_grnd => pwf%qflx_snow_grnd
       qflx_rain_grnd => pwf%qflx_rain_grnd
#if (defined PCP2PFT)
       !CAS  Redistribute precipitation across PFT types
       !CAS  
       !CAS  NOTE: This modification has been employed for use with the zonal configuration
       !CAS        of CLM coupled to the MIT Integrated Global System Model (IGSM)
       !CAS
       pcp2pft => pps%pcp2pft
       !write (6,*) 'Hydrology1: PCP2PFT @ column: ',pi,': ',pcp2pft
       !write (6,*) 'at i = ',pps%ixy,' j = ',pps%jxy,' and PFT = ',pps%itype
       forc_rain = forc_rain*pcp2pft
       forc_snow = forc_snow*pcp2pft
#endif

       ! Canopy interception and precipitation onto ground surface
       ! Add precipitation to leaf water
       if (itypwat==istsoil .OR. itypwat==istwet) then

          qflx_candrip = 0.                 ! rate of canopy runoff
          qflx_through = 0.                 ! precipitation direct through canopy
          qflx_prec_intr = 0.           ! intercepted precipitation
          prcp = forc_rain + forc_snow  ! total precipitation

          if (frac_veg_nosno == 1 .AND. prcp > 0.) then

             ! The leaf water capacities for solid and liquid are different, 
             ! generally double for snow, but these are of somewhat less significance
             ! for the water budget because of lower evap. rate at lower temperature.
             ! Hence, it is reasonable to assume that vegetation storage of solid water 
             ! is the same as liquid water.
             h2ocanmx = dewmx * (elai + esai)

             ! Direct throughfall
             fpi = 1. - exp(-0.5*(elai + esai))
             qflx_through  = prcp*(1.-fpi)

             ! Water storage of intercepted precipitation and dew
             qflx_prec_intr = prcp*fpi
             h2ocan = max(0._r8, h2ocan + dtime*qflx_prec_intr)

             ! Initialize rate of canopy runoff and snow falling off canopy
             qflx_candrip = 0.0

             ! Excess water that exceeds the leaf capacity
             xrun = (h2ocan - h2ocanmx)/dtime

             ! Test on maximum dew on leaf
             if (xrun > 0.) then
                qflx_candrip = xrun
                h2ocan = h2ocanmx
             endif

          endif

       else if (itypwat == istice) then

          qflx_prec_intr = 0.
          h2ocan = 0.
          qflx_candrip = 0.
          qflx_through = 0.

       endif

       ! Precipitation onto ground (kg/(m2 s))
       if (frac_veg_nosno == 0) then
          qflx_prec_grnd = forc_rain + forc_snow
       else
          qflx_prec_grnd = qflx_through + qflx_candrip
       endif

       if (itypprc <= 1) then

          ! The percentage of liquid water by mass, which is arbitrarily set to 
          ! vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.
          ! rain (itypprc = 1) or no precip (itypprc = 0)
          flfall = 1.                        ! fraction of liquid water within falling precip.
          if (do_capsnow) then
             qflx_snowcap = qflx_prec_grnd
             qflx_snow_grnd = 0.
             qflx_rain_grnd = 0.
          else
             qflx_snowcap = 0.
             qflx_snow_grnd = 0.              ! ice onto ground (mm/s)
             qflx_rain_grnd = qflx_prec_grnd  ! liquid water onto ground (mm/s)
          endif

       else

          ! snow (itypprc = 2)
          ! calculate the fraction of liquid water in snow
#if (defined PERGRO)
          ! Note for cleanup: this PERGRO block has the same functional form as the
          ! non-PERGRO block - not sure why.
          coefb = 0.4/2.0
          coefa = -coefb*tfrz
          if (forc_t <= tfrz) then
             flfall = 0.0
          else if (forc_t <= tfrz+2.) then
             flfall = coefa + coefb*forc_t
          else
             flfall = coefa + coefb*(tfrz+2.)
          endif
#else
          if (forc_t <= tfrz) then
             flfall = 0.
          else if (forc_t <= tfrz+2.) then
             flfall = -54.632 + 0.2*forc_t
          else
             flfall = 0.4
          endif
#endif
          if (do_capsnow) then
             qflx_snowcap = qflx_prec_grnd
             qflx_snow_grnd = 0.
             qflx_rain_grnd = 0.
          else
             qflx_snowcap = 0.
             qflx_snow_grnd = qflx_prec_grnd*(1.-flfall)
             qflx_rain_grnd = qflx_prec_grnd*flfall
          endif

       endif

       ! Determine the fraction of foliage covered by water and the 
       ! fraction of foliage that is dry and transpiring.

       call FracWet(p)

#if (defined PCP2PFT)
          !CAS  Reverse the Redistribution of precipitation across PFT types
          !CAS  
          !CAS  NOTE: THIS IS DONE SO THAT CLM WILL PASS WATER/MASS CONSERVATION
          !CAS        AS WELL AS TO ALLOW FOR SUBSEQUENT PFT ADJUSTMENT OF PCP
          !CAS  
          !CAS  NOTE: This modification has been employed for use with the zonal configuration
          !CAS        of CLM coupled to the MIT Integrated Global System Model (IGSM)
          !CAS
!CAS:DEBUG
!          write (6,*) 'PCP2PFT  = ',pcp2pft
!CAS:DEBUG
          forc_rain = forc_rain/pcp2pft
          forc_snow = forc_snow/pcp2pft
          prcp = forc_rain + forc_snow  ! total precipitation
#endif
    end do ! (end pfts loop)

    ! This routine includes an update of some of the column level
    ! state variables for snow.  This requires knowing the column-level
    ! average of qflx_snow_grnd (taken over all pfts).
    ! This averaging normally happens outside the Hydrology1
    ! subroutine, but for the first pass at conversion to new data
    ! structures, I am introducing a special averaging loop here
    ! for only this variable

    pftsum = 0.0
    do pi=1,cps%npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_snow_grnd * c%pw(pi)   
    end do
    cwf%pwf_a%qflx_snow_grnd = pftsum
    qflx_snow_grnd => cwf%pwf_a%qflx_snow_grnd

    ! Use Alta relationship, Anderson(1976); LaChapelle(1961), 
    ! U.S.Department of Agriculture Forest Service, Project F, 
    ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.
    if (itypprc <= 1) then

       ! rain (itypprc = 1) or no precip (itypprc = 0)
       dz_snowf = 0.     ! rate of snowfall, snow depth/s (m/s)

    else

       ! snow (itypprc = 2)
       if (do_capsnow) then
          dz_snowf = 0.
       else
#if (defined PERGRO)
          if (forc_t > tfrz + 2.) then
             bifall=50. + 1.7*(17.0)**1.5
          else if (forc_t > tfrz - 15.) then
             bifall=50. + 1.7*(forc_t - tfrz + 15.)**1.5
          else
             bifall=50.
          endif
#else
          if (forc_t > tfrz + 2.) then
             bifall =189.
          else if (forc_t > tfrz - 15.) then
             bifall=50. + 1.7*(forc_t - tfrz + 15.)**1.5
          else
             bifall=50.
          endif
#endif
          dz_snowf = qflx_snow_grnd/bifall
          snowdp = snowdp + dz_snowf*dtime
          h2osno = h2osno + qflx_snow_grnd*dtime  ! snow water equivalent (mm)
       endif

       if (itypwat==istwet .AND. t_grnd>tfrz) then
          h2osno=0.
          snowdp=0.
          snowage=0.
       endif

    endif

    ! When the snow accumulation exceeds 10 mm, initialize snow layer
    ! Currently, the water temperature for the precipitation is simply set 
    ! as the surface air temperature

    newnode = 0    ! flag for when snow node will be initialized
    if (snl == 0 .AND. qflx_snow_grnd > 0.0 .AND. snowdp >= 0.01) then
       newnode = 1
       snl = -1
       dz(0) = snowdp                       ! meter
       z(0) = -0.5*dz(0)
       zi(-1) = -dz(0)
       snowage = 0.                         ! snow age
       t_soisno (0) = min(tfrz, forc_t)     ! K
       h2osoi_ice(0) = h2osno               ! kg/m2
       h2osoi_liq(0) = 0.                   ! kg/m2
       frac_iceold(0) = 1.
    endif

    ! The change of ice partial density of surface node due to precipitation.
    ! Only ice part of snowfall is added here, the liquid part will be added later

    if (snl < 0 .AND. newnode == 0) then
       h2osoi_ice(snl+1) = h2osoi_ice(snl+1)+dtime*qflx_snow_grnd
       dz(snl+1) = dz(snl+1)+dz_snowf*dtime
    endif


  end subroutine Hydrology1

end module Hydrology1Mod




