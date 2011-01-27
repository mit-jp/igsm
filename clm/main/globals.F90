module globals

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: globals
!
! !DESCRIPTION: 
! Module of global time-related control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
  real(r8):: dtime      !land model time step (sec)
  integer :: nstep      !time step number
  type control_type
     integer :: year    !current year (0 -> ...)
     integer :: month   !current month (1 -> 12)
     integer :: day     !current day (1 -> 31)
     integer :: secs    !seconds of current date
  end type control_type
  type(control_type) :: ctl !control variables
!                              
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

end module globals
