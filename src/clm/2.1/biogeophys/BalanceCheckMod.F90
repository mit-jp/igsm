#include <misc.h>
#include <preproc.h>

module BalanceCheckMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: 
! 
! !DESCRIPTION: 
! Water and energy balance check  
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BalanceCheck ! Water and energy balance check  
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
! !IROUTINE: BalanceCheck
!
! !INTERFACE:
  subroutine BalanceCheck (c) 
!
! !DESCRIPTION: 
! Water and energy balance check  
! This subroutine accumulates the numerical truncation errors of the water
! and energy balance calculation. It is helpful to see the performance of 
! the process of integration.
!
! The error for energy balance: 
! error = abs(Net radiation - the change of internal energy - Sensible heat
!             - Latent heat) 
! The error should be less than 0.02 W/m2 in each time integration interval;
!
! The error for water balance:
! error = abs(precipitation - change of water storage - evaporation - runoff)
! The error should be less than 0.001 mm in  each time integration interval.
!
! !USES:
    use clmtype
    use globals, only: dtime, nstep
#if (defined PCP2PFT || defined COUP_MIT2D)
    use clmpoint
#endif
#if (defined COUP_MIT2D)
    use clm_varpar, only : maxpatch_pft
#endif
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout) :: c	!column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 10 November 2000: Mariana Vertenstein
! Migrated to new data structures by Mariana Vertenstein and 
! Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
    real(r8), pointer :: forc_rain     !rain rate [mm/s]
    real(r8), pointer :: forc_snow     !snow rate [mm/s]
    real(r8), pointer :: forc_lwrad    !downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: endwb         !water mass end of the time step
    real(r8), pointer :: begwb         !water mass begining of the time step
    real(r8), pointer :: fsa           !solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsr           !solar radiation reflected (W/m**2)
    real(r8), pointer :: eflx_lwrad_out!emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net!net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: sabv          !solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: sabg          !solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: eflx_sh_tot   !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot   !total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_soil_grnd!soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: qflx_evap_tot !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_surf     !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl    !qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_drain    !sub-surface runoff (mm H2O /s)
!
! local pointers to original implicit inout scalars
!
    real(r8), pointer :: acc_errh2o    !accumulation of water balance error
    real(r8), pointer :: acc_errseb    !accumulation of surface energy balance error
!
! local pointers to original implicit out scalars
!
    real(r8), pointer :: errh2o        !water conservation error (mm H2O)
    real(r8), pointer :: errsol        !solar radiation conservation error (W/m**2)
    real(r8), pointer :: errlon        !longwave radiation conservation error (W/m**2)
    real(r8), pointer :: errseb        !surface energy conservation error (W/m**2)
!
! local pointers to original implicit in scalars
!
    real(r8), pointer :: latdeg    !latitude (radians)
!
! local pointers to original implicit in arrays
!
    real(r8),dimension(:), pointer :: forc_solad   !direct beam radiation (vis=forc_sols , nir=forc_soll )
    real(r8),dimension(:), pointer :: forc_solai   !diffuse radiation     (vis=forc_solsd, nir=forc_solld)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: j                   !do loop index
    integer :: m                   !indices
    logical :: constop = .false.   !true => stop if energy balance err too great
    type(atm2lnd_flux_type)  , pointer :: a2lf   !local pointers to derived subtypes
    type(water_balance_type) , pointer :: cwbal  !local pointers to derived subtypes
    type(energy_balance_type), pointer :: cebal  !local pointers to derived subtypes
    type(column_pstate_type) , pointer :: cps    !local pointers to derived subtypes
    type(column_eflux_type)  , pointer :: cef    !local pointers to derived subtypes
    type(column_wflux_type)  , pointer :: cwf    !local pointers to derived subtypes
#if (defined PCP2PFT || defined COUP_MIT2D)
!CAS  For redistribution of precipitation across PFT types
!CAS  
!CAS  NOTE: This modification has been employed for use with the zonal configuration
!CAS        of CLM coupled to the MIT Integrated Global System Model (IGSM)
!CAS
    type(pft_type)             , pointer :: p    ! local pointers to derived subtypes
!    type(pft_pstate_type)      , pointer :: pps  ! local pointers to derived subtypes
!    real(r8), pointer :: pcp2pft                 !PCP Distribution factor across PFTs
    integer :: pi                                !Index to PFT 
#endif

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes
    a2lf => c%a2lf
    cwbal => c%cwbal
    cebal => c%cebal
    cps => c%cps
    cef => c%cef
    cwf => c%cwf
#if (defined PCP2PFT || defined COUP_MIT2D)
    pi = cps%index1d
#endif
#if (defined PCP2PFT)
    p => ppoint%pft(pi)%p
#endif
    latdeg => cps%lps%gps%latdeg
    ! Assign local pointers to derived type scalar members
    forc_rain => a2lf%forc_rain
    forc_snow => a2lf%forc_snow
    forc_lwrad => a2lf%forc_lwrad
    endwb => cwbal%endwb
    begwb => cwbal%begwb
    fsa => cef%pef_a%fsa
    fsr => cef%pef_a%fsr
    eflx_lwrad_out => cef%pef_a%eflx_lwrad_out
    eflx_lwrad_net => cef%pef_a%eflx_lwrad_net
    sabv => cef%pef_a%sabv
    sabg => cef%pef_a%sabg
    eflx_sh_tot => cef%pef_a%eflx_sh_tot
    eflx_lh_tot => cef%pef_a%eflx_lh_tot
    eflx_soil_grnd => cef%pef_a%eflx_soil_grnd
    qflx_evap_tot => cwf%pwf_a%qflx_evap_tot
    qflx_surf => cwf%qflx_surf
    qflx_qrgwl => cwf%qflx_qrgwl
    qflx_drain => cwf%qflx_drain
    errh2o => cwbal%errh2o
    errsol => cebal%errsol
    errseb => cebal%errseb
    errlon => cebal%errlon
    acc_errh2o => cwbal%acc_errh2o
    acc_errseb => cebal%acc_errseb

    ! Assign local pointers to derived type array members
    forc_solad => a2lf%forc_solad
    forc_solai => a2lf%forc_solai

    ! Water balance 
#if (defined PCP2PFT)
    !write (6,*) 'WaterBalance: PCP2PFT = ',p%pps%pcp2pft
    forc_rain = forc_rain*p%pps%pcp2pft
    forc_snow = forc_snow*p%pps%pcp2pft
#endif
    errh2o = endwb - begwb - &
         ( forc_rain  + forc_snow - qflx_evap_tot - qflx_surf &
         - qflx_qrgwl - qflx_drain ) * dtime

#if (defined PCP2PFT)
    !write (6,*) 'WaterBalance: PCP2PFT = ',p%pps%pcp2pft
    forc_rain = forc_rain/p%pps%pcp2pft
    forc_snow = forc_snow/p%pps%pcp2pft
#endif
#if (defined COUP_MIT2D)
    !m = pfts1d%itypwat(pi)
    !if ( m > 1 ) then
    !if (forc_rain+forc_snow == 0.0 ) then
    !  qflx_qrgwl=0.0+((begwb-endwb)/dtime)
    !else
    !  qflx_qrgwl=max(qflx_qrgwl,0.0)
    !endif
    !endif
#endif
    if (abs(errh2o) > .10) then
       write(6,200)'water balance error',nstep,cps%index1d,errh2o
       write(6,*)'clm model is stopping'
       call endrun
    endif

    ! Solar radiation energy balance
    errsol = fsa + fsr - (forc_solad(1) + forc_solad(2) &
         + forc_solai(1) + forc_solai(2))

    if (abs(errsol) > .10 ) then
       write(6,100)'solar radiation balance error',nstep,cps%index1d,errsol
       write(6,*)'clm model is stopping'
       call endrun
    endif

    ! Longwave radiation energy balance
    errlon = eflx_lwrad_out - eflx_lwrad_net - forc_lwrad

    if (abs(errlon) > .10 ) then
       write(6,100)'longwave enery balance error',nstep,cps%index1d,errlon
       write(6,*)'at latitude: ',latdeg
       write(6,*)'clm model is stopping'
       call endrun
    endif

    ! Surface energy balance
    errseb = sabv + sabg  &
         + forc_lwrad - eflx_lwrad_out &
         - eflx_sh_tot &
         - eflx_lh_tot &
         - eflx_soil_grnd

    if (abs(errseb) > .10 ) then
       write(6,100)'surface flux energy balance error',nstep,cps%index1d,errseb
       write(6,*)'clm model is stopping'
       call endrun
    endif

    ! Accumulation of water and surface energy balance error
    acc_errh2o = acc_errh2o + errh2o
    acc_errseb = acc_errseb + errseb

100 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' W/m2') 
200 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' mm') 

  end subroutine BalanceCheck

end module BalanceCheckMod
