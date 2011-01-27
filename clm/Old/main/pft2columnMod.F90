#include <misc.h>
#include <preproc.h>

module pft2columnMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: pft2columnMod
! 
! !DESCRIPTION: 
! Contains methods to perfom averages over all pfts on a column to provide
! column level values 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pft_to_col_glob  !averages over column pfts for soil/lake variables
  public :: pft_to_col_soil  !averages over column pfts for soil variables
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
! !IROUTINE: pft_to_col_glob
!
! !INTERFACE: subroutine pft_to_col_glob(c)
  subroutine pft_to_col_glob(c)
!
! !DESCRIPTION: 
! Averages over all pfts for variables defined over both soil and lake 
! to provide the column-level averages of state and flux variables 
! defined at the PFT level.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout) :: c !column derived type
!
! !REVISION HISTORY:
! Peter Thornton, 02.03.05
! Mariana Vertenstein, 02.09.17 Modified to split lake and non-lake
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), dimension(:), pointer ::  wt !pointer to array of pft weights
    integer :: pi      !pft index
    integer :: npfts   !number of pfts
    real(r8):: pftsum  !temporary used for pft averaging for columns
! -----------------------------------------------------------------

    ! Set initial pointers

    wt => c%pw
    npfts = c%cps%npfts

    ! Averaging for pft physical state variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%ndvi * wt(pi)   
    end do
    c%cps%pps_a%ndvi = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%elai * wt(pi)   
    end do
    c%cps%pps_a%elai = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%esai * wt(pi)   
    end do
    c%cps%pps_a%esai = pftsum

    ! Averaging for pft energy state variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pes%t_ref2m * wt(pi)   
    end do
    c%ces%pes_a%t_ref2m = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pes%t_af * wt(pi)   
    end do
    c%ces%pes_a%t_af = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pes%t_veg * wt(pi)   
    end do
    c%ces%pes_a%t_veg = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pes%t_rad_pft * wt(pi)   
    end do
    c%ces%pes_a%t_rad_pft = pftsum

    ! Averaging for pft energy flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lh_tot * wt(pi)   
    end do
    c%cef%pef_a%eflx_lh_tot = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lh_grnd * wt(pi)   
    end do
    c%cef%pef_a%eflx_lh_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_soil_grnd * wt(pi)   
    end do
    c%cef%pef_a%eflx_soil_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_sh_tot * wt(pi)   
    end do
    c%cef%pef_a%eflx_sh_tot = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_sh_grnd * wt(pi)   
    end do
    c%cef%pef_a%eflx_sh_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lwrad_out * wt(pi)   
    end do
    c%cef%pef_a%eflx_lwrad_out = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lwrad_net * wt(pi)   
    end do
    c%cef%pef_a%eflx_lwrad_net = pftsum

    ! Averaging for pft momentum flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pmf%taux * wt(pi)   
    end do
    c%cmf%pmf_a%taux = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pmf%tauy * wt(pi)   
    end do
    c%cmf%pmf_a%tauy = pftsum

    ! Averaging for pft energy flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%sabg * wt(pi)   
    end do
    c%cef%pef_a%sabg = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%sabv * wt(pi)   
    end do
    c%cef%pef_a%sabv = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%fsa * wt(pi)   
    end do
    c%cef%pef_a%fsa = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%fsr * wt(pi)   
    end do
    c%cef%pef_a%fsr = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%parsun * wt(pi)   
    end do
    c%cef%pef_a%parsun = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%parsha * wt(pi)   
    end do
    c%cef%pef_a%parsha = pftsum

    ! Averaging for pft water flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_prec_grnd * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_prec_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_evap_soi * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_evap_soi = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_evap_tot * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_evap_tot = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_efpot * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_efpot = pftsum

    return
  end subroutine pft_to_col_glob

!-----------------------------------------------------------------------
!!BOP
!
! !IROUTINE: pft_to_col_soil
!
! !INTERFACE: subroutine pft_to_col_soil(c)
  subroutine pft_to_col_soil(c)
!
! !DESCRIPTION: 
! Averages over all pfts for variables defined only over soil 
! (i.e. non-lake) to provide the column-level averages of state and 
! flux variables defined at the PFT level.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
!
! !ARGUMENTS:
  implicit none
  type (column_type), target, intent(inout) :: c !column derived type
!   	
! !REVISION HISTORY:
! Peter Thornton, 02.03.05
! Mariana Vertenstein, 02.09.17 Modified to split lake and non-lake
!
!!EOP
!
! !LOCAL VARIABLES:
    real(r8), dimension(:), pointer ::  wt !pointer to array of pft weights
    integer :: pi      !pft index
    integer :: npfts   !number of pfts
    real(r8):: pftsum  !temporary used for pft averaging for columns
! -----------------------------------------------------------------

    ! Set initial pointers

    wt => c%pw
    npfts = c%cps%npfts

    ! Averaging for pft physical state variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%rssun * wt(pi)   
    end do
    c%cps%pps_a%rssun = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%rssha * wt(pi)   
    end do
    c%cps%pps_a%rssha = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pps%btran * wt(pi)   
    end do
    c%cps%pps_a%btran = pftsum

    ! Averaging for pft water state variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pws%h2ocan * wt(pi)   
    end do
    c%cws%pws_a%h2ocan = pftsum

    ! Averaging for pft energy flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%dlrad * wt(pi)   
    end do
    c%cef%pef_a%dlrad = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%ulrad * wt(pi)   
    end do
    c%cef%pef_a%ulrad = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_sh_veg * wt(pi)   
    end do
    c%cef%pef_a%eflx_sh_veg = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lh_vege * wt(pi)   
    end do
    c%cef%pef_a%eflx_lh_vege = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_lh_vegt * wt(pi)   
    end do
    c%cef%pef_a%eflx_lh_vegt = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%cgrnd * wt(pi)   
    end do
    c%cef%pef_a%cgrnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%cgrndl * wt(pi)   
    end do
    c%cef%pef_a%cgrndl = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%cgrnds * wt(pi)   
    end do
    c%cef%pef_a%cgrnds = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%eflx_gnet * wt(pi)   
    end do
    c%cef%pef_a%eflx_gnet = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pef%dgnetdT * wt(pi)   
    end do
    c%cef%pef_a%dgnetdT = pftsum

    ! Averaging for pft water flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_prec_intr * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_prec_intr = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_rain_grnd * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_rain_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_snow_grnd * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_snow_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_snowcap * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_snowcap = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_evap_veg * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_evap_veg = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_tran_veg * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_tran_veg = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_evap_can * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_evap_can = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_evap_grnd * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_evap_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_dew_grnd * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_dew_grnd = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_sub_snow * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_sub_snow = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pwf%qflx_dew_snow * wt(pi)   
    end do
    c%cwf%pwf_a%qflx_dew_snow = pftsum


    ! Averaging for pft carbon flux variables

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pcf%psnsun * wt(pi)   
    end do
    c%ccf%pcf_a%psnsun = pftsum

    pftsum = 0.0
    do pi=1,npfts 
       pftsum = pftsum + c%p(pi)%pcf%psnsha * wt(pi)   
    end do
    c%ccf%pcf_a%psnsha = pftsum

    return
  end subroutine pft_to_col_soil

end module pft2columnMod
