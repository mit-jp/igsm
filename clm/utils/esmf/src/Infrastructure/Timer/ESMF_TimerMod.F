! $Id$
	module ESMF_TimerMod
!===============================================================================
!BOP
! !MODULE: ESMF_TimerMod
!
! !USES:
	Use ESMF_LogMod
!
! !PUBLIC TYPES:
	implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!     ESMF_TimerInit
!     ESMF_TimerStart
!     ESMF_TimerStop
!     ESMF_TimerStamp
!     ESMF_TimerPrint
!
! !PUBLIC DATA MEMBERS:
	integer, parameter ::  ! This list is documented under ESMF_TimerInit
     &	ESMF_RESERVED = 0,
     &	ESMF_USRSYS = 1,
     &	ESMF_WALL = 2,
     &  ESMF_RESERVED1 = 3, 
     &	ESMF_LD1CACHE_MISS = 4,
     &	ESMF_LD2CACHE_MISS = 5,
     &	ESMF_CYCLES = 6,
     &	ESMF_ELAPSED_CYCLES = 7,
     &	ESMF_FP_INSTR = 8,
     &	ESMF_LOADSTORE_INSTR = 9,
     &	ESMF_INSTR = 10,
     &	ESMF_STALL = 11,
     &  ESMF_RESERVED2 = 12, 
     &	ESMF_STUB_TIMERS = 13
	character*(*), parameter ::
     &	ESMF_ALL_TIMERS = "*timers*"
!
! !DESCRIPTION:
!     Outputs diagnostic information.  Sorts the information based on process
!     and provides a more robust solution that simple writes.
!EOP
!===============================================================================
	contains
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_TimerStart
!
! !INTERFACE:
        subroutine ESMF_TimerStart(name, rc)
!
! !PARAMETERS:
        character *(*) :: name                          ! Timer name
        integer, intent(out), optional :: rc            ! Return code.
!
! !DESCRIPTION:
!     Starts the named timer.  If the timer does not yet exist, it is created.
!
!EOP
!-------------------------------------------------------------------------------
	integer stub

        call ESMC_TimerStart(name, stub)
	if (present(rc)) rc = stub

        end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_TimerStop
!
! !INTERFACE:
        subroutine ESMF_TimerStop(name, rc)
!
! !PARAMETERS:
        character *(*) :: name                          ! Timer name.
        integer, intent(out), optional :: rc            ! Return code.
!
! !DESCRIPTION:
!     Stops the named timer.
!
!EOP
!-------------------------------------------------------------------------------
	integer stub

        call ESMC_TimerStop(name, stub)
	if (present(rc)) rc = stub

        end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_TimerStamp
!
! !INTERFACE:
        subroutine ESMF_TimerStamp(wall, user, sys, rc)
!
! !PARAMETERS:
        real(8), intent(out) :: wall                    ! Wall clock time.
        real(8), intent(out) :: user                    ! User time.
        real(8), intent(out) :: sys                     ! System clock time.
        integer, intent(out), optional :: rc            ! Return code.
!
! !DESCRIPTION:
!     Retrieve the clock times at the time of calling.
!
!EOP
!-------------------------------------------------------------------------------
	integer stub

        call ESMC_TimerStamp(wall, user, sys, stub)
	if (present(rc)) rc = stub

        end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_TimerPrint
!
! !INTERFACE:
        subroutine ESMF_TimerPrint(name, log, rc)
        integer, intent(out), optional :: rc            ! Return code.
!
! !PARAMETERS:
        character *(*), optional :: name                ! Timer name.
        type(ESMF_Log), optional :: log                 ! Log to use for printing.
!
! !DESCRIPTION:
!     Print timing results.  If {\tt ESMF\_ALL\_TIMERS} is sent as the parameter,
!     the results of all timers are sent to the log.  If another timer name is
!     used, only that timer is printed.  
!     If the {\tt log} argument is not present, the standard application log
!     is used.
!
!EOP
!-------------------------------------------------------------------------------
	type(ESMF_Log) :: log1
	integer :: stub
	integer*8 :: nolog

	if (present(name)) then
	else
	    print *, "Name not yet optional"
	    call exit(1)
	endif

	if (present(log)) then 
	    call ESMC_TimerPrint(name, log, stub)
	else
	    nolog = 0
            call ESMC_TimerPrint(name, nolog, stub)
	endif

	if (present(rc)) rc = stub

	end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_TimerInit
!
! !INTERFACE:
!       subroutine ESMF_TimerInit(name, option1, value1, ....)
!
! !PARAMETERS:
!       character *(*) :: name               ! Timer name
!       enum :: option1                      ! First option
!       integer :: value1                    ! First value
!       ... (repeat option, value)
!
! !DESCRIPTION:
!     Initialize the timer.  As many option,value pairs may be passed in as
!     needed.  The list must be terminated with a zero to signify the
!     end of the arguments.
!     This table summarizes the options that may be passed to the
!     {\tt ESMF\_TimerInit} call.  Each option takes an associated value of
!     either 1 (true/enable) or 0 (false/disable).
!\smallskip
!
!\begin{tabular}{lll}
!{\tt \bf Timer Opt.} & {\tt ESMF\_RESERVED[\#]} & Not to be used.  Reserve values for library. \\
!                 & {\tt ESMF\_USRSYS} & (Tracks/reports) Unix user and system time.  \\
!                 & {\tt ESMF\_WALL} & Wall clock time (absolute run-time).  \\
!{\tt \bf Profiling Opt.} & {\tt ESMF\_LD1CACHE\_MISS} & Level 1 cache misses.  \\
!                 & {\tt ESMF\_LD2CACHE\_MISS} & Level 2 cache misses.  \\
!                 & {\tt ESMF\_CYCLES} & Cycles spent in process/(thread).  \\
!                 & {\tt ESMF\_ELAPSED\_CYCLES} & Total cycles spent (all processes).  \\
!                 & {\tt ESMF\_FP\_INSTR} & Floating point operations. \\
!                 & {\tt ESMF\_LOADSTORE\_INSTR} & Load and store operations. \\
!                 & {\tt ESMF\_INSTR} & Sum of all completed instructions. \\
!                 & {\tt ESMF\_STALL} & Sum of all cycles the process/(thread) is stalled. \\
!                 & {\tt ESMF\_STUB\_TIMERS} & Causes the library to return immediately from calls. \\
!\end{tabular}
!
!EOP
!-------------------------------------------------------------------------------

	end module
