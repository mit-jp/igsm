#include <misc.h>
#include <preproc.h>

module EcosystemDynMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: EcosystemDynMod
! 
! !DESCRIPTION: 
! Ecosystem dynamics: phenology, vegetation
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: EcosystemDyn
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
! !IROUTINE: EcosystemDyn
!
! !INTERFACE:
  subroutine EcosystemDyn (c, doalb, endofyr)
!
! !DESCRIPTION: 
! Ecosystem dynamics: phenology, vegetation
! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
! height (htop)
!
! !USES:
    use clmtype
    use globals
    use mvegFileMod , only : mlai1, msai1, mlai2, msai2, &
         mhvt1, mhvt2, mhvb1, mhvb2, timwt
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c	!column derived type
    logical, intent(in):: doalb   !true = surface albedo calculation time step
    logical, intent(in):: endofyr
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrated to new data structure.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: snowdp   !snow height (m)
!
! local pointers to implicit in/out scalars
!
    real(r8), pointer :: htop     !canopy top (m)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: tlai     !one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai     !one-sided stem area index, no burying by snow
    real(r8), pointer :: hbot     !canopy bottom (m)
    real(r8), pointer :: elai     !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai     !one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb !frac of vegetation not covered by snow [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer ::  indexp            !index into global 1d pft array
    real(r8) :: ol                !thickness of canopy layer covered by snow (m)
    real(r8) :: fb                !fraction of canopy layer covered by snow
    integer  :: pi                !pft index

    ! local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
!-----------------------------------------------------------------------

    if (doalb) then

       ! Assign local pointers to derived subtypes (column-level)
       cps => c%cps

       ! Assign local pointers to derived type scalar members (column-level)
       snowdp => cps%snowdp

       ! Loop through pfts
       do pi=1, cps%npfts

          ! Assign local pointers to pft-level derived subtypes
          p => c%p(pi)
          pps => p%pps

          ! Assign local pointers to derived type scalar members (pft-level)
          tlai => pps%tlai 
          tsai => pps%tsai 
          elai => pps%elai 
          esai => pps%esai 
          htop => pps%htop 
          hbot => pps%hbot 
          frac_veg_nosno_alb => pps%frac_veg_nosno_alb 

          ! Determine 1d index
          indexp = pps%index1d

          ! need to update elai and esai only every albedo time step so do not 
          ! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e., 
          ! if albedos are not done every time step).
          ! leaf phenology
          ! Set leaf and stem areas based on day of year
          ! Interpolate leaf area index, stem area index, and vegetation heights
          ! between two monthly 
          ! The weights below (timwt(1) and timwt(2)) were obtained by a call to 
          ! routine InterpMonthlyVeg in subroutine NCARlsm. 
          !                 Field   Monthly Values
          !                -------------------------
          ! leaf area index LAI  <- mlai1 and mlai2
          ! leaf area index SAI  <- msai1 and msai2
          ! top height      HTOP <- mhvt1 and mhvt2
          ! bottom height   HBOT <- mhvb1 and mhvb2

          tlai = timwt(1)*mlai1(indexp) + timwt(2)*mlai2(indexp)
          tsai = timwt(1)*msai1(indexp) + timwt(2)*msai2(indexp)
          htop = timwt(1)*mhvt1(indexp) + timwt(2)*mhvt2(indexp)
          hbot = timwt(1)*mhvb1(indexp) + timwt(2)*mhvb2(indexp)

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05, set equal to zero to prevent numerical 
          ! problems associated with very small lai and sai.

          ol = min( max(snowdp-hbot,0._r8), htop-hbot)
          fb = 1. - ol / max(1.e-06, htop-hbot)
          elai = max(tlai*fb,0.0_r8)
          esai = max(tsai*fb,0.0_r8)
          if (elai < 0.05) elai = 0._r8
          if (esai < 0.05) esai = 0._r8

          ! Fraction of vegetation free of snow

          if ((elai + esai) >= 0.05) then
             frac_veg_nosno_alb = 1
          else
             frac_veg_nosno_alb = 0
          endif

       end do ! end of pft loop

    endif  !end of if-doalb block

    return
  end subroutine EcosystemDyn

end  module EcosystemDynMod





