#include <misc.h>
#include <preproc.h>

module FracWetMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: FracWetMod
! 
! !DESCRIPTION: 
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry 
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FracWet
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
! !IROUTINE: FracWet
!
! !INTERFACE:
  subroutine FracWet (p)
!
! !DESCRIPTION: 
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry. The variable "fwet" is the  
! fraction of all vegetation surfaces which are wet including 
! stem area which contribute to evaporation. The variable "fdry"
! is the fraction of elai which is dry because only leaves
! can transpire.  Adjusted for stem area which does not transpire.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type (pft_type),target,intent(inout):: p		!pft derived type
!
! !CALLED FROM:
! subroutine Hydrology in module Hydrology1Mod     
!
! !REVISION HISTORY:
! Created by Keith Oleson and M. Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer,pointer :: frac_veg_nosno
    real(r8),pointer:: dewmx
    real(r8),pointer:: elai
    real(r8),pointer:: esai
    real(r8),pointer:: h2ocan
!
! local pointers to original implicit out scalars
!
    real(r8),pointer:: fwet
    real(r8),pointer:: fdry
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    real(r8):: vegt             ! frac_veg_nosno*lsai
    real(r8):: dewmxi           ! inverse of maximum allowed dew [1/mm]
    type(pft_pstate_type),pointer:: pps !local pointers to derived subtypes
    type(pft_wstate_type),pointer:: pws !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    pps => p%pps
    pws => p%pws

    ! assign local pointers to derived type scalar members
    frac_veg_nosno => pps%frac_veg_nosno
    dewmx => pps%dewmx
    elai => pps%elai
    esai => pps%esai
    h2ocan => pws%h2ocan
    fwet => pps%fwet
    fdry => pps%fdry

    ! begin calculation
    if (frac_veg_nosno == 1) then
       if (h2ocan > 0.) then
          vegt     = frac_veg_nosno*(elai + esai)
          dewmxi   = 1.0/dewmx
          fwet = ((dewmxi/vegt)*h2ocan)**.666666666666
          fwet = min (fwet,1.0_r8)     ! Check for maximum limit of fwet
       else
          fwet = 0.
       endif
       fdry = (1.-fwet)*elai/(elai+esai)
#if (defined PERGRO)
       fwet = 0.
       fdry = elai/(elai+esai)
#endif
    else
       fwet = 0.
       fdry = 0.
    endif

  end subroutine FracWet

end module FracWetMod
