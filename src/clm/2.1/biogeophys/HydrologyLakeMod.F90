#include <misc.h>
#include <preproc.h>

module HydrologyLakeMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: HydrologyLakeMod
! 
! !DESCRIPTION: 
! Calculate lake hydrology
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: HydrologyLake
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
! !IROUTINE: HydrologyLake
!
! !INTERFACE:
  subroutine HydrologyLake (c) 
!
! !DESCRIPTION: 
! Calculate lake hydrology
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : hfus, tfrz, spval
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout) :: c !column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 3/4/02: Peter Thornton; Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: begwb         !water mass begining of the time step
    real(r8), pointer :: forc_snow     !snow rate [mm/s]
    real(r8), pointer :: forc_rain     !rain rate [mm/s]
    logical , pointer :: do_capsnow    !true => do snow capping
    real(r8), pointer :: t_grnd        !ground temperature (Kelvin)
    real(r8), pointer :: qmelt         !snow melt [mm/s]
    real(r8), pointer :: qflx_evap_soi !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
!
! local pointers to implicit inout scalars
!
    real(r8), pointer :: h2osno        !snow water (mm H2O)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: endwb         !water mass end of the time step
    real(r8), pointer :: snowdp        !snow height (m)
    real(r8), pointer :: snowice       !average snow ice lens
    real(r8), pointer :: snowliq       !average snow liquid water
    real(r8), pointer :: eflx_snomelt  !snow melt heat flux (W/m**2)
    real(r8), pointer :: qflx_infl     !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_snomelt  !snow melt (mm H2O /s)
    real(r8), pointer :: qflx_surf     !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_drain    !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl    !qflx_surf at glaciers, wetlands, lakes
!
! local pointers to implicit out arrays
!
    real(r8), dimension(:), pointer :: rootr_column !effective fraction of roots in each soil layer
    real(r8), dimension(:), pointer :: h2osoi_vol   !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), dimension(:), pointer :: h2osoi_ice   !ice lens (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_liq   !liquid water (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: j              !do loop index
    real(r8):: qflx_evap_grnd !ground surface evaporation rate (mm h2o/s)
    real(r8):: qflx_dew_grnd  !ground surface dew formation (mm h2o /s) [+]
    real(r8):: qflx_sub_snow  !sublimation rate from snow pack (mm h2o /s) [+]
    real(r8):: qflx_dew_snow  !surface dew added to snow pack (mm h2o /s) [+]
    type(water_balance_type),pointer :: cwbal !local pointers to derived subtypes
    type(atm2lnd_flux_type) ,pointer :: a2lf  !local pointers to derived subtypes 
    type(column_pstate_type),pointer :: cps   !local pointers to derived subtypes
    type(column_estate_type),pointer :: ces   !local pointers to derived subtypes
    type(column_wstate_type),pointer :: cws   !local pointers to derived subtypes
    type(column_eflux_type) ,pointer :: cef   !local pointers to derived subtypes
    type(column_wflux_type) ,pointer :: cwf   !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! Asign local pointers to derived subtypes

    cwbal => c%cwbal
    a2lf => c%a2lf
    cps => c%cps
    ces => c%ces
    cws => c%cws
    cef => c%cef
    cwf => c%cwf

    ! Assign local pointers to derived type scalar members

    begwb => cwbal%begwb
    endwb => cwbal%endwb
    forc_snow => a2lf%forc_snow
    forc_rain => a2lf%forc_rain
    do_capsnow => cps%do_capsnow
    snowdp => cps%snowdp
    t_grnd => ces%t_grnd
    h2osno => cws%h2osno
    snowice => cws%snowice
    snowliq => cws%snowliq
    eflx_snomelt => cef%eflx_snomelt
    qmelt => cwf%qmelt
    qflx_snomelt => cwf%qflx_snomelt
    qflx_surf => cwf%qflx_surf
    qflx_qrgwl => cwf%qflx_qrgwl
    qflx_drain => cwf%qflx_drain
    qflx_infl => cwf%qflx_infl
    qflx_evap_soi => c%p(1)%pwf%qflx_evap_soi
    qflx_evap_tot => c%p(1)%pwf%qflx_evap_tot

    ! Assign local pointers to derived type array members

    rootr_column => cps%rootr_column
    h2osoi_vol => cws%h2osoi_vol
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! Snow on the lake ice 
    ! Note that these are only local variables, as per the original
    ! Hydrology_Lake code. So even though these names correspond to
    ! variables in clmtype, this routine is not updating the
    ! values of the clmtype variables. (PET, 3/4/02)

    qflx_evap_grnd = 0.
    qflx_sub_snow = 0.
    qflx_dew_snow = 0.
    qflx_dew_grnd = 0.

    if (qflx_evap_soi >= 0.) then

       ! Sublimation: do not allow for more sublimation than there is snow
       ! after melt.  Remaining surface evaporation used for infiltration.

       qflx_sub_snow = min( qflx_evap_soi, h2osno/dtime-qmelt )
       qflx_evap_grnd = qflx_evap_soi - qflx_sub_snow

    else

       if (t_grnd < tfrz-0.1) then
          qflx_dew_snow = abs(qflx_evap_soi)
       else
          qflx_dew_grnd = abs(qflx_evap_soi)
       endif

    endif

    ! Update snow pack

    if (do_capsnow) then
       h2osno = h2osno - (qmelt + qflx_sub_snow)*dtime
    else
       h2osno = h2osno + (forc_snow-qmelt-qflx_sub_snow+qflx_dew_snow)*dtime
    endif
    h2osno = max( h2osno, 0._r8 )

    ! No snow if lake unfrozen

    if (t_grnd > tfrz) h2osno = 0.

    ! Snow depth

    snowdp = h2osno/250.    !Assume a constant snow bulk density = 250.

    ! Determine ending water balance

    endwb = h2osno

    ! The following are needed for global average on history tape. 
    ! Note that components that are not displayed over lake on history tape
    ! must be set to spval here

    eflx_snomelt = qmelt*hfus
    qflx_infl = 0.
    qflx_snomelt = qmelt
    qflx_surf = 0.
    qflx_drain = 0.
    rootr_column(:) = spval
    snowice = spval
    snowliq = spval
    h2osoi_vol(:) = spval
    h2osoi_ice(:) = spval
    h2osoi_liq(:) = spval
    qflx_qrgwl = forc_rain + forc_snow - qflx_evap_tot - (endwb-begwb)/dtime

    ! The pft average must be done here for output to history tape

    c%cwf%pwf_a%qflx_evap_tot = qflx_evap_tot

  end subroutine HydrologyLake

end module HydrologyLakeMod
