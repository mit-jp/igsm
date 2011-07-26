! $Id$
        module ESMF_LogMod
!===============================================================================
!BOP
! !MODULE: ESMF_LogMod
!
! !USES:
!
! !PUBLIC TYPES:
#include "ESMF_conf.h"
        implicit none

        type ESMF_Log
        private
        character(ESMC_POINTER_SIZE) :: this
        end type
!
! !PUBLIC MEMBER FUNCTIONS:
!     ESMF_LogPrint
!     ESMF_LogNew
!     ESMF_LogDelete
!     ESMF_LogSetState
!     ESMF_LogFlush
!
! !PUBLIC DATA MEMBERS:
! Log Levels
        integer, parameter :: 
     &  ESMF_LOGLEVEL_INFO = 0,
     &  ESMF_LOGLEVEL_ERROR = 1,
     &  ESMF_LOGLEVEL_TIMING = 2

! Log states.
        integer, parameter ::
     &	ESMF_LOGSTATE_QUIET = 0,
     &	ESMF_LOGSTATE_NORMAL = 1,
     &	ESMF_LOGSTATE_TIMER = 2,
     &	ESMF_LOGSTATE_VERBOSE = 3

! Log configurables
        integer, parameter ::
     &	ESMF_LOGCONFIG_PRINTPID = 1,
     &	ESMF_LOGCONFIG_PRINTTID = 2,
     &	ESMF_LOGCONFIG_PRINTTIME = 3,
     &	ESMF_LOGCONFIG_LOGSTATE = 4
!
! !DESCRIPTION:
!     Outputs diagnostic information.  Sorts the information based on process
!     and provides a more robust solution that simple writes.
!
!EOP
!===============================================================================
        contains
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogPrint
!
! !INTERFACE:
!      subroutine ESMF_LogPrint(log, priority, format, ...)

! !PARAMETERS:
!      type(ESMF_Log) :: log          ! Log to use.
!      integer :: priority            ! The log level
!      character(*) :: format         ! A format string
!      ...                            ! Variables

! !DESCRIPTION:
!     Interface to the logging mechanism.  Works as a printf for Fortran.  
!     This function is not actually in this module, but is in the
!     c-file interface.  This is the case since F90 Modules do not support
!     a varargs interface.  This is one of the few functions in the
!     ESMF library that does not take an optional error return code.  If
!     one is passed in, it will be ignored.
!
!     Each message that is output via the log will have an associated logging 
!     level.  This level is used to categorize what type of message is being output
!     and what the priority/category of the message is.
!     The logging levels will be defined as a list of integer parameters:
!
!\begin{itemize}
!\item{\bf ESMF\_LOGLEVEL\_INFO}
!
!     Marks a message that would be useful for someone debugging the code to have.
!
!\item{\bf ESMF\_LOGLEVEL\_ERROR}
!
!     Reports an error condition.
!
!\item{\bf ESMF\_LOGLEVEL\_TIMING}
!
!     Used to write information related to code profiling or timing. 
!
!\end{itemize}
!
!EOP
!-------------------------------------------------------------------------------
!
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogNew
!
! !INTERFACE:
        function ESMF_LogNew(logname, logstate, labelio, rc)
!
! !RETURN VALUE:
        type(ESMF_Log) :: ESMF_LogNew           ! Returns a handle to the log.
!
! !PARAMETERS:
         character*(*), intent(in) :: logname    ! Log filename.
         integer, intent(in) :: logstate         ! Initial log state.    
         integer, intent(in) :: labelio          ! True == separate files/process.
         integer, intent(out), optional :: rc    ! Return code.
!
! !DESCRIPTION:
!     Creates and returns the log. The log will output to a file or file
!     formed from {\tt logname}.  The {\tt logstate} argument sets the default
!     log state.  {\tt Labelio} set to true creates a separate file for each
!     process.
!
!EOP
!------------------------------------------------------------------------------- 
        integer stub

        call ESMC_LogNew(ESMF_LogNew, logname, logstate, labelio, stub)
        if (present(rc)) rc = stub
                                                                                 
        end function
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogDelete
!
! !INTERFACE:
        subroutine ESMF_LogDelete(log, rc)
!
! !PARAMETERS:
        type(ESMF_Log) :: log                   ! The log to destroy
        integer, intent(out), optional :: rc    ! Return code.
!
! !DESCRIPTION:
!     Closes all log files and frees space created by log.
!
!EOP
!------------------------------------------------------------------------------- 
        integer stub

        call ESMC_LogDelete(log, stub)
        if (present(rc)) rc = stub
                                                                                 
        end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogSetConfig
!
! !INTERFACE:
!	subroutine ESMF_LogSetConfig(log, option1, value1, ...)
!
! !PARAMETERS:
!       type(ESMF_Log) :: log                   ! The log to configure.
!       integer        :: option1               ! First option to set
!       ????           :: value1                ! First option value
!       ...                                     ! Repeat
!
! !DESCRIPTION:
!     Configures the log object.  List must end in {\tt ESMF\_MULL} to signal the 
!     end of the option list.  This function is a stdarg function and is implemented
!     in the C interface.  This function does not take an optional return code.
!
!EOP
!------------------------------------------------------------------------------- 
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogGetConfig
!
! !INTERFACE:
!	subroutine ESMF_LogGetConfig(log, option1, value1, ...)
!
! !PARAMETERS:
!       type(ESMF_Log) :: log                   ! The log to configure.
!       integer        :: option1               ! First option to get
!       ????           :: value1                ! Variable for first option value
!       ...                                     ! Repeat
!
! !DESCRIPTION:
!     Retrieves configuration of a log object.  List must end in {\tt ESMF\_MULL} to signal the 
!     end of the option list.  This function is a stdarg function and is implemented
!     in the C interface.  This function does not take an optional return code.
!
!EOP
!------------------------------------------------------------------------------- 
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogSetState
!
! !INTERFACE:
        subroutine ESMF_LogSetState(log, logstate, rc)
!
! !PARAMETERS:
        type(ESMF_Log) :: log                   ! Log to set.
        integer, intent(in) :: logstate         ! Initial log state.    
        integer, intent(out), optional :: rc    ! Return code.
!
! !DESCRIPTION:
!     Set a new log state.
!     The program may be run at several different logging states.  Each state will
!     allow the messages at various logging levels to be output. 
!
!\begin{itemize}
!\item{\bf ESMF\_LOGSTATE\_QUIET}
!
!     Suppresses the output of all message types.  All calls to the log will
!     return immediately after verifying this flag.  This will be the fastest mode of 
!     operation, the only faster being to turn the library off via the compile flags.
!
!\item{\bf ESMF\_LOGSTATE\_NORMAL}
!
!     This option will print only the types of messages that you would want to have 
!     under normal operations.
!
!\item{\bf ESMF\_LOGSTATE\_TIMER}
!
!     Prints out only profiling/timing information.
!
!\item{\bf ESMF\_LOGSTATE\_VERBOSE}
!
!     Outputs every type of message.
!
!\end{itemize}
!
!
!EOP
!------------------------------------------------------------------------------- 
        integer stub

        call ESMC_LogSetState(log, logstate, stub)
        if (present(rc)) rc = stub
                                                                                 
        end subroutine
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_LogFlush()
!
! !INTERFACE:
        subroutine ESMF_LogFlush(log, rc)
!
! !PARAMETERS:
        type(ESMF_Log), intent(in) :: log       ! The log to Flush.
        integer, intent(out), optional :: rc    ! Return code.
!
! !DESCRIPTION:
!     Flushes the log's output file.
!
!EOP
!------------------------------------------------------------------------------- 
        integer stub

        call ESMC_LogFlush(log, stub)
        if (present(rc)) rc = stub
                                                                                 
        end subroutine

        end module
