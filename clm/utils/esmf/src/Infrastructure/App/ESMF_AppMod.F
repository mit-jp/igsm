! $Id$
	module ESMF_AppMod
!===============================================================================
!BOP
! !MODULE: ESMF_AppMod
!
! !USES:
	use ESMF_BasicUtilMod
	use ESMF_ErrorMod
	use ESMF_LogMod
!
! !PUBLIC TYPES:
	implicit none
#include "ESMF_conf.h"
	type ESMF_App
	private
	character(ESMC_POINTER_SIZE) :: this
	end type
!
! !PUBLIC MEMBER FUNCTIONS:
!      ESMF_AppNew
!      ESMF_AppDelete
!
! !DESCRIPTION:
! The object contains the underlying structure for the entire ESMF library.  
! This object must be created before any other part of the library may be used.
!
!EOP
!===============================================================================

	contains

!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_AppNew
!
! !INTERFACE:
	function ESMF_AppNew(rc)
!
! !RETURN VALUE:
	type(ESMF_App) :: ESMF_AppNew
!
! !PARAMETERS:
	integer, intent(out), optional :: rc           ! Return code.
!
! !DESCRIPTION:
!     Initialize the Application.
!
!EOP
!------------------------------------------------------------------------------
	integer stub

	call ESMC_AppNew(ESMF_AppNew, stub)
	if (present(rc)) rc = stub
                                                                                 
	end function
!===============================================================================
!BOP
!
! !IROUTINE:  ESMF_AppDelete
!
! !INTERFACE:
	subroutine ESMF_AppDelete(app, rc)
!
! !PARAMETERS:
	type(ESMF_App) :: app
	integer, intent(out), optional :: rc          ! Return code.
!
! !DESCRIPTION:
!     Finalizes the application.
!
!EOP
!------------------------------------------------------------------------------- 
	integer stub

	call ESMC_AppDelete(app, stub)
	if (present(rc)) rc = stub

	end subroutine

	end module
