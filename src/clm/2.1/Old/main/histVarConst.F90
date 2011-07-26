#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: subourinte histvarconst
!
! !INTERFACE: 
subroutine histVarConst()
!
! !DESCRIPTION: 
! Sets variables values over lake points that are needed for
! global 2d grid output on history file
!
! !USES:
  use clmtype
  use clmpoint, only : cpoint
  use shr_kind_mod, only: r8 => shr_kind_r8
  use histFileMod, only : spval
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine iniTimeConst in module iniTimeConstMod
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: ci,pi                  !indices
  type(column_type), pointer :: c   !pointer to derived subtype
! -----------------------------------------------------------------

  do ci = cols1d%beg,cols1d%end
     c => cpoint%col(ci)%c

     if (c%cps%lps%lakpoi) then  

        c%ces%t_snow = spval
        do pi = 1, c%cps%npfts

           c%p(pi)%pps%rssun = spval
           c%cps%pps_a%rssun = spval

           c%p(pi)%pps%rssha = spval
           c%cps%pps_a%rssha = spval

           c%p(pi)%pps%btran = spval
           c%cps%pps_a%btran = spval

           c%p(pi)%pwf%qflx_tran_veg = 0.
           c%cwf%pwf_a%qflx_tran_veg = 0.

           c%p(pi)%pwf%qflx_efpot = 0.
           c%cwf%pwf_a%qflx_efpot = 0.

           c%p(pi)%pwf%qflx_evap_can = 0.
           c%cwf%pwf_a%qflx_evap_can = 0.

           c%p(pi)%pwf%qflx_prec_intr = 0.
           c%cwf%pwf_a%qflx_prec_intr = 0.

           c%p(pi)%pws%h2ocan = 0.
           c%cws%pws_a%h2ocan = 0.

           c%p(pi)%pef%eflx_lh_vegt = 0.
           c%cef%pef_a%eflx_lh_vegt = 0.

           c%p(pi)%pef%eflx_lh_vege = 0.
           c%cef%pef_a%eflx_lh_vege = 0.

        end do

     endif

  end do

end subroutine histVarConst
