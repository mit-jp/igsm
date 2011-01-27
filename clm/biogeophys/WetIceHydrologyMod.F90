#include <misc.h>
#include <preproc.h>

module WetIceHydrologyMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: WetIceHydrologyMod
! 
! !DESCRIPTION: 
! Calculate hydrology for ice and wetland!
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: WetIceHydrology  ! Calculate hydrology for ice and wetland
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
! !IROUTINE: WetIceHydrology
!
! !INTERFACE:
  subroutine WetIceHydrology (c)
!
! !DESCRIPTION: 
! Calculate hydrology for ice and wetland
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : istwet, istice
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 7 November 2000: Mariana Vertenstein; Initial code
! 2/28/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    real(r8),pointer:: forc_rain     !rain rate [mm/s]
    real(r8),pointer:: forc_snow     !snow rate [mm/s]
    real(r8),pointer:: qflx_evap_tot !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8),pointer:: endwb         !water mass end of the time step
    real(r8),pointer:: begwb         !water mass begining of the time step
!
! local pointers to original implicit out scalars
!
    real(r8),pointer:: qflx_drain    !sub-surface runoff (mm H2O /s)
    real(r8),pointer:: qflx_surf     !surface runoff (mm H2O /s)
    real(r8),pointer:: qflx_infl     !infiltration (mm H2O /s)
    real(r8),pointer:: qflx_qrgwl    !qflx_surf at glaciers, wetlands, lakes
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
! local pointers to derived subtypes
!
    type(atm2lnd_flux_type),pointer::   a2lf
    type(water_balance_type),pointer::  cwbal
    type(column_wflux_type),pointer::   cwf
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    a2lf => c%a2lf
    cwbal => c%cwbal
    cwf => c%cwf

    ! assign local pointers to derived type scalar members
    forc_rain => a2lf%forc_rain
    forc_snow => a2lf%forc_snow
    endwb => cwbal%endwb
    begwb => cwbal%begwb
    qflx_evap_tot => cwf%pwf_a%qflx_evap_tot
    qflx_drain => cwf%qflx_drain
    qflx_surf => cwf%qflx_surf
    qflx_infl => cwf%qflx_infl
    qflx_qrgwl => cwf%qflx_qrgwl

    ! Wetland and land ice runoff

    qflx_drain  = 0.
    qflx_surf   = 0.
    qflx_infl   = 0.
    qflx_qrgwl  = forc_rain + forc_snow - qflx_evap_tot - &
         (endwb - begwb)/dtime
  end subroutine WetIceHydrology

end module WetIceHydrologyMod
