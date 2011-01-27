!===============================================================================
! CVS: $Id$
! CVS: $Source$
! CVS: $Name$
!===============================================================================

module shr_timer_mod

   !----------------------------------------------------------------------------
   !
   ! routines that support multiple CPU timers via F90 intrisics
   !
   ! Note: 
   ! o if   an operation is requested on an invalid timer number n
   !   then nothing is done in a routine
   ! o if   more than max_timers are requested, 
   !   then timer n=max_timers is "overloaded" and becomes invalid/undefined
   !----------------------------------------------------------------------------

   use shr_kind_mod

   implicit none

   private  ! resticted access
   public  :: shr_timer_init , shr_timer_get      , &
   &          shr_timer_start, shr_timer_stop     , &
   &          shr_timer_print, shr_timer_print_all, &
   &          shr_timer_check, shr_timer_check_all, &
   &          shr_timer_zero , shr_timer_zero_all , &
   &          shr_timer_free , shr_timer_free_all , &
   &          shr_timer_sleep

   integer(shr_kind_in),parameter :: stat_free    = 0  ! timer status constants
   integer(shr_kind_in),parameter :: stat_inuse   = 1
   integer(shr_kind_in),parameter :: stat_started = 2
   integer(shr_kind_in),parameter :: stat_stopped = 3
   integer(shr_kind_in),parameter :: max_timers   = 100 ! max number of timers

   integer(shr_kind_in) :: status (max_timers) ! status of each timer
   integer(shr_kind_in) :: cycles1(max_timers) ! cycle number at timer start 
   integer(shr_kind_in) :: cycles2(max_timers) ! cycle number at timer stop  
   character   (len=80) :: name   (max_timers) ! name assigned to each timer
   real   (shr_kind_r8) :: dt     (max_timers) ! accumulated time
   integer(shr_kind_in) :: calls  (max_timers) ! # of samples in accumulation
   integer(shr_kind_in) :: cycles_max = -1     ! max cycles before wrapping
   real   (shr_kind_r8) :: clock_rate          ! clock_rate: seconds per cycle

   save

!===============================================================================

   contains

!===============================================================================

subroutine shr_timer_init

   !----- local -----
   integer(shr_kind_in) :: cycles ! count rate return by system clock

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_init) ',a,i5)"

!-------------------------------------------------------------------------------
!
! This routine initializes:
! 1) values in all timer array locations
! 2) machine parameters necessary for computing cpu time from F90 intrinsics.
!    F90 intrinsic: system_clock(count_rate=cycles, count_max=cycles_max)
!-------------------------------------------------------------------------------

   call shr_timer_free_all

   call system_clock(count_rate=cycles, count_max=cycles_max)

   if (cycles /= 0) then
     clock_rate = 1.0/real(cycles)
   else
     clock_rate = 0
     write(6,F00) 'ERROR: no system clock available' 
   endif

end subroutine shr_timer_init

!===============================================================================

subroutine shr_timer_get(n, str)

   !----- arguments -----
   integer(shr_kind_in),intent(out) :: n    ! timer number 
   character (*)       ,intent( in) :: str  ! text string with timer name

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_get) ',a,i5)"

!-----------------------------------------------------------------------
!
!  search for next free timer
!
!-----------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_free) then
       status(n) = stat_inuse
       name  (n) = str
       calls (n) = 0
       return
     endif
   end do

   n=max_timers
   name  (n) = "<invalid - undefined - overloaded>"
   write(6,F00) 'ERROR: exceeded maximum number of timers'

end subroutine shr_timer_get

!===============================================================================

subroutine shr_timer_start(n)

   !----- arguments -----
   integer(shr_kind_in), intent(in) :: n      ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_start) ',a,i5)"

!-----------------------------------------------------------------------
!
!  This routine starts a given timer.
!
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) call shr_timer_stop(n)

     status(n) = stat_started
     call system_clock(count=cycles1(n))
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_start
 
!===============================================================================

subroutine shr_timer_stop(n)

   !----- arguments -----
   integer(shr_kind_in), intent(in) :: n  ! timer number

   !----- local -----
   real (shr_kind_r8) :: elapse      ! elapsed time returned by system counter

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_stop) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine stops a given timer, checks for cycle wrapping, computes the 
!  elpased time, and accumulates the elpased time in the dt(n) array
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if ( status(n) == stat_started) then
       call system_clock(count=cycles2(n))
       if (cycles2(n) >= cycles1(n)) then
         dt(n) = dt(n) + clock_rate*(cycles2(n) - cycles1(n))
       else
         dt(n) = dt(n) + clock_rate*(cycles_max + cycles2(n) - cycles1(n))
       endif
       calls (n) = calls (n) + 1
       status(n) = stat_stopped
     end if
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_stop
 
!===============================================================================

subroutine shr_timer_print(n)

   !----- arguments -----
   integer(shr_kind_in), intent(in) :: n     ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_print) ',a,i5)"
   character(len=*),parameter :: F01 = "('(shr_timer_print) timer',i3,&
   &                                     ':',i8,' calls,',f10.3,'s, id: ',a)"
!-------------------------------------------------------------------------------
!
!  prints the accumulated time for a given timer
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call shr_timer_stop(n)
       write (6,F01) n,dt(n),trim(name(n))
       call shr_timer_start(n)
     else
       write (6,F01) n,calls(n),dt(n),trim(name(n))
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_print

!===============================================================================

subroutine shr_timer_print_all

   !----- local -----
   integer(shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_print_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  prints accumulated time for all timers in use
!
!-------------------------------------------------------------------------------

   write(6,F00) 'print all timing info:'

   do n=1,max_timers
     if (status(n) /= stat_free) call shr_timer_print(n)
   end do

end subroutine shr_timer_print_all

!===============================================================================

subroutine shr_timer_zero(n)

   !----- arguments -----
   integer(shr_kind_in), intent(in) :: n       ! timer number
   
   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_zero) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine resets a given timer.
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     dt(n) = 0.0
     calls(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_zero

!===============================================================================

subroutine shr_timer_zero_all

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_zero_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine resets all timers.
!
!-------------------------------------------------------------------------------

   dt = 0.0
   calls = 0

end subroutine shr_timer_zero_all

!===============================================================================

subroutine shr_timer_check(n)

   !----- arguments -----
   integer(shr_kind_in), intent(in) ::  n   ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_check) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  This routine checks a given timer.  This is primarily used to
!  periodically accumulate time in the timer to prevent timer cycles
!  from wrapping around max_cycles.
!
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call shr_timer_stop (n)
       call shr_timer_start(n)
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_check

!===============================================================================

subroutine shr_timer_check_all

   !----- local -----
   integer(shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_check_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  Call shr_timer_check for all timers in use
!
!-------------------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_started) then
       call shr_timer_stop (n)
       call shr_timer_start(n)
     endif
   end do

end subroutine shr_timer_check_all

!===============================================================================

subroutine shr_timer_free(n)

   !----- arguments -----
   integer(shr_kind_in),intent(in) :: n    ! timer number 

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_free) ',a,i5)"

!-----------------------------------------------------------------------
!
!  initialize/free all timer array values
!
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     status (n) = stat_free
     name   (n) = "<invalid - undefined>"
     dt     (n) = 0.0
     cycles1(n) = 0
     cycles2(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_free

!===============================================================================

subroutine shr_timer_free_all

   !----- local -----
   integer(shr_kind_in) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_free_all) ',a,i5)"

!-------------------------------------------------------------------------------
!
!  initialize/free all timer array values
!
!-------------------------------------------------------------------------------

   do n=1,max_timers
     call shr_timer_free(n)
   end do

end subroutine shr_timer_free_all

!===============================================================================

subroutine shr_timer_sleep(sec)

   use shr_sys_mod     ! share system calls (namely, shr_sys_sleep)

   !----- local -----
   real   (shr_kind_r8),intent(in) :: sec  ! number of seconds to sleep

!-------------------------------------------------------------------------------
! Sleep for approximately sec seconds
!
! Note: sleep is typically a system call, hence it is implemented in 
!       shr_sys_mod, although it probably would only be used in a timing 
!       context, which is why there is a shr_timer_* wrapper provided here.
!-------------------------------------------------------------------------------

   call shr_sys_sleep(sec)

end subroutine shr_timer_sleep

!===============================================================================
end module shr_timer_mod
!===============================================================================
