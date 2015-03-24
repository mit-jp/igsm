!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id: m_clmmpif.F90 18 2005-12-12 17:49:42Z mvr $
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_clmmpif - a portable interface to the CLMMPI "mpif.h" COMMONs.
!
! !DESCRIPTION:
!
!   The purpose of \verb"m_clmmpif" module is to provide a portable
!   interface of \verb"mpif.h" with different CLMMPI implementation.
!   By combining module \verb"m_clmmpif" and \verb"m_clmmpif90", it may be
!   possible to build a Fortran 90 CLMMPI binding module graduately.
!
!   Although it is possible to use \verb'include "mpif.h"' directly
!   in individual modules, it has several problems:
!   \begin{itemize}
!   \item It may conflict with either the source code of a {\sl fixed}
!	format or the code of a {\sl free} format;
!   \item It does not provide the protection and the safety of using
!	these variables as what a \verb"MODULE" would provide.
!   \end{itemize}
!
!   More information may be found in the module \verb"m_clmmpif90".
!
! !INTERFACE:

	module m_clmmpif
	  implicit none
	  private	! except

	  public :: CLMMPI_INTEGER
	  public :: CLMMPI_REAL
	  public :: CLMMPI_DOUBLE_PRECISION
	  public :: CLMMPI_LOGICAL
	  public :: CLMMPI_CHARACTER

	  public :: CLMMPI_REAL4
	  public :: CLMMPI_REAL8

	  public :: CLMMPI_COMM_WORLD
	  public :: CLMMPI_COMM_NULL

	  public :: CLMMPI_SUM
	  public :: CLMMPI_PROD
	  public :: CLMMPI_MIN
  	  public :: CLMMPI_MAX

	  public :: CLMMPI_MAX_ERROR_STRING
	  public :: CLMMPI_STATUS_SIZE
	  public :: CLMMPI_ANY_SOURCE

#ifdef CLMMPICH_
	  public :: CLMMPIPRIV	! the common block name
#endif

	  include "mpif.h"

! !REVISION HISTORY:
! 	01Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP
!_______________________________________________________________________
	character(len=*),parameter :: myname='MCT(MPEU)::m_clmmpif'

	end module m_clmmpif
!.
