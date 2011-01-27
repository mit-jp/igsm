#include <misc.h>
#include <preproc.h>

module TridiagonalMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: TridiagonalMod
! 
! !DESCRIPTION: 
! Tridiagonal matrix solution
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Tridiagonal
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
! !IROUTINE: Tridiagonal
!
! !INTERFACE:
  subroutine Tridiagonal (n, a, b, c, r, u)
!
! !DESCRIPTION: 
! Tridiagonal matrix solution
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: n
    real(r8), intent(in)  :: a(1:n), b(1:n), c(1:n), r(1:n)
    real(r8), intent(out) :: u(1:n)
!
! !CALLED FROM:
!subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!subroutine SoilTemperature in module SoilTemperatureMod
!subroutine SoilWater in module HydrologyMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!EOP
!
! !LOCAL VARIABLES:
    integer j
    real(r8) gam(1:n)
    real(r8) bet
!-----------------------------------------------------------------------

    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j) * gam(j)
       u(j) = (r(j) - a(j)*u(j-1)) / bet
    enddo

    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1) * u(j+1)
    enddo

  end subroutine Tridiagonal

end module TridiagonalMod
