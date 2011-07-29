!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id: m_clmmpif90.F90 18 2005-12-12 17:49:42Z mvr $
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_clmmpif90 - a Fortran 90 style CLMMPI module interface.
!
! !DESCRIPTION:
!
!   By wrapping \verb'include "mpif.h"' into a module, \verb"m_clmmpif()"
!   provides an easy way to
!\begin{itemize}
!  \item avoid the problem with {\sl fixed} or {\sl free} formatted
!	Fortran 90 files;
!  \item provide protections with only a limited set of \verb"PUBLIC"
!	variables; and
!  \item be extended to a CLMMPI Fortran 90 binding.
!\end{itemize}
!
! !INTERFACE:

    module m_clmmpif90
      use m_clmmpif, only : MP_INTEGER	=> CLMMPI_INTEGER
      use m_clmmpif, only : MP_REAL	=> CLMMPI_REAL
      use m_clmmpif, only : MP_DOUBLE_PRECISION	&
					=> CLMMPI_DOUBLE_PRECISION
      use m_clmmpif, only : MP_LOGICAL	=> CLMMPI_LOGICAL
      use m_clmmpif, only : MP_CHARACTER	=> CLMMPI_CHARACTER

      use m_clmmpif, only : MP_REAL4	=> CLMMPI_REAL4
      use m_clmmpif, only : MP_REAL8	=> CLMMPI_REAL8

      use m_clmmpif, only : MP_COMM_WORLD	=> CLMMPI_COMM_WORLD
      use m_clmmpif, only : MP_COMM_NULL	=> CLMMPI_COMM_NULL
      use m_clmmpif, only : MP_SUM		=> CLMMPI_SUM
      use m_clmmpif, only : MP_PROD	=> CLMMPI_PROD
      use m_clmmpif, only : MP_MIN 	=> CLMMPI_MIN
      use m_clmmpif, only : MP_MAX 	=> CLMMPI_MAX
      use m_clmmpif, only : MP_MAX_ERROR_STRING	&
					=> CLMMPI_MAX_ERROR_STRING
      use m_clmmpif, only : MP_STATUS_SIZE => CLMMPI_STATUS_SIZE
      use m_clmmpif, only : MP_ANY_SOURCE	=> CLMMPI_ANY_SOURCE

      implicit none
      private

      public :: MP_type

      public :: MP_INTEGER
      public :: MP_REAL
      public :: MP_DOUBLE_PRECISION
      public :: MP_LOGICAL
      public :: MP_CHARACTER

      public :: MP_REAL4
      public :: MP_REAL8

      public :: MP_COMM_WORLD
      public :: MP_COMM_NULL

      public :: MP_SUM
      public :: MP_PROD
      public :: MP_MIN
      public :: MP_MAX

      public :: MP_ANY_SOURCE

      public :: MP_MAX_ERROR_STRING

      public :: MP_init
      public :: MP_initialized
      public :: MP_finalize
      public :: MP_abort

      public :: MP_wtime
      public :: MP_wtick

      public :: MP_comm_size
      public :: MP_comm_rank
      public :: MP_comm_dup
      public :: MP_comm_free

      public :: MP_cart_create
      public :: MP_dims_create
      public :: MP_cart_coords
      public :: MP_cart_rank

      public :: MP_error_string

      public :: MP_perr

      public :: MP_STATUS_SIZE
      public :: MP_status

      public :: MP_log2

! !REVISION HISTORY:
! 	09Dec97 - Jing Guo <guo@thunder> - initial prototyping/coding.
!		. started with everything public, without any interface
!		  declaration.
!		. Then limited to only variables current expected to
!		  be used.
!	
!EOP
!_______________________________________________________________________

integer,dimension(MP_STATUS_SIZE) :: MP_status

	!----------------------------------------

interface MP_init
  subroutine CLMMPI_init(ier)
    integer,intent(out) :: ier
  end subroutine CLMMPI_init
end interface

interface MP_initialized
  subroutine CLMMPI_initialized(flag,ier)
    logical,intent(out) :: flag
    integer,intent(out) :: ier
  end subroutine CLMMPI_initialized
end interface

interface MP_finalize
  subroutine CLMMPI_finalize(ier)
    integer,intent(out) :: ier
  end subroutine CLMMPI_finalize
end interface

interface MP_error_string
  subroutine CLMMPI_error_string(ierror,cerror,ln,ier)
    integer,intent(in) :: ierror
    character(len=*),intent(out) :: cerror
    integer,intent(out) :: ln
    integer,intent(out) :: ier
  end subroutine CLMMPI_error_string
end interface

interface MP_type; module procedure	&
  typeI_,	& ! CLMMPI_INTEGER
  typeL_,	& ! CLMMPI_LOGICAL
  typeC_,	& ! CLMMPI_CHARACTER
  typeSP_,	& ! CLMMPI_REAL
  typeDP_,	& ! CLMMPI_DOUBLE_PRECISION
  typeI1_,	& ! CLMMPI_INTEGER
  typeL1_,	& ! CLMMPI_LOGICAL
  typeC1_,	& ! CLMMPI_CHARACTER
  typeSP1_,	& ! CLMMPI_REAL
  typeDP1_,	& ! CLMMPI_DOUBLE_PRECISION
  typeI2_,	& ! CLMMPI_INTEGER
  typeL2_,	& ! CLMMPI_LOGICAL
  typeC2_,	& ! CLMMPI_CHARACTER
  typeSP2_,	& ! CLMMPI_REAL
  typeDP2_	  ! CLMMPI_DOUBLE_PRECISION
end interface

interface MP_perr; module procedure perr_; end interface

interface MP_abort
  subroutine CLMMPI_abort(comm,errorcode,ier)
    integer,intent(in) :: comm
    integer,intent(in) :: errorcode
    integer,intent(out) :: ier
  end subroutine CLMMPI_abort
end interface

	!----------------------------------------
interface MP_wtime
  function CLMMPI_wtime()
    double precision :: CLMMPI_wtime
  end function CLMMPI_wtime
end interface

interface MP_wtick
  function CLMMPI_wtick()
    double precision :: CLMMPI_wtick
  end function CLMMPI_wtick
end interface

	!----------------------------------------
interface MP_comm_size
  subroutine CLMMPI_comm_size(comm,size,ier)
    integer,intent(in) :: comm
    integer,intent(out) :: size
    integer,intent(out) :: ier
  end subroutine CLMMPI_comm_size
end interface

interface MP_comm_rank
  subroutine CLMMPI_comm_rank(comm,rank,ier)
    integer,intent(in) :: comm
    integer,intent(out) :: rank
    integer,intent(out) :: ier
  end subroutine CLMMPI_comm_rank
end interface

interface MP_comm_dup
  subroutine CLMMPI_comm_dup(comm,newcomm,ier)
    integer,intent(in) :: comm
    integer,intent(out) :: newcomm
    integer,intent(out) :: ier
  end subroutine CLMMPI_comm_dup
end interface

interface MP_comm_free
  subroutine CLMMPI_comm_free(comm,ier)
    integer,intent(inout) :: comm
    integer,intent(out) :: ier
  end subroutine CLMMPI_comm_free
end interface

	!----------------------------------------
interface MP_cart_create
  subroutine CLMMPI_cart_create(comm_old,ndims,dims,periods,	&
  	reorder,comm_cart,ier)
    integer,intent(in) :: comm_old
    integer,intent(in) :: ndims
    integer,dimension(*),intent(in) :: dims
    logical,dimension(*),intent(in) :: periods
    logical,             intent(in) :: reorder
    integer,intent(out) :: comm_cart
    integer,intent(out) :: ier
  end subroutine CLMMPI_cart_create
end interface

interface MP_dims_create
  subroutine CLMMPI_dims_create(nnodes,ndims,dims,ier)
    integer,intent(in) :: nnodes
    integer,intent(in) :: ndims
    integer,dimension(*),intent(inout) :: dims
    integer,intent(out) :: ier
  end subroutine CLMMPI_dims_create
end interface

interface MP_cart_coords
  subroutine CLMMPI_cart_coords(comm,rank,maxdims,coords,ier)
    integer,intent(in) :: comm
    integer,intent(in) :: rank
    integer,intent(in) :: maxdims
    integer,dimension(*),intent(out) :: coords
    integer,intent(out) :: ier
  end subroutine CLMMPI_cart_coords
end interface

interface MP_cart_rank
  subroutine CLMMPI_cart_rank(comm,coords,rank,ier)
    integer,intent(in) :: comm
    integer,dimension(*),intent(in) :: coords
    integer,intent(out) :: rank
    integer,intent(out) :: ier
  end subroutine CLMMPI_cart_rank
end interface
	!----------------------------------------

  character(len=*),parameter :: myname='m_clmmpif90'
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI_ - return CLMMPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeI_(ival)
      implicit none
      integer,intent(in) :: ival
      integer :: typeI_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeI_'

  typeI_=MP_INTEGER

end function typeI_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL_ - return CLMMPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeL_(lval)
      implicit none
      logical,intent(in) :: lval
      integer :: typeL_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeL_'

  typeL_=MP_LOGICAL

end function typeL_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC_ - return CLMMPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeC_(cval)
      implicit none
      character(len=*),intent(in) :: cval
      integer :: typeC_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeC_'

  typeC_=MP_CHARACTER

end function typeC_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP_ - return CLMMPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeSP_(rval)
      use m_realkinds,only : SP
      implicit none
      real(SP),intent(in) :: rval
      integer :: typeSP_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeSP_'

  typeSP_=MP_REAL

end function typeSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP_ - return CLMMPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeDP_(rval)
      use m_realkinds,only : DP
      implicit none
      real(DP),intent(in) :: rval
      integer :: typeDP_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeDP_'

  typeDP_=MP_DOUBLE_PRECISION

end function typeDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI1_ - return CLMMPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeI1_(ival)
      implicit none
      integer,dimension(:),intent(in) :: ival
      integer :: typeI1_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeI1_'

  typeI1_=MP_INTEGER

end function typeI1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL1_ - return CLMMPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeL1_(lval)
      implicit none
      logical,dimension(:),intent(in) :: lval
      integer :: typeL1_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeL1_'

  typeL1_=MP_LOGICAL

end function typeL1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC1_ - return CLMMPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeC1_(cval)
      implicit none
      character(len=*),dimension(:),intent(in) :: cval
      integer :: typeC1_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeC1_'

  typeC1_=MP_CHARACTER

end function typeC1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP1_ - return CLMMPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeSP1_(rval)
      use m_realkinds,only : SP
      implicit none
      real(SP),dimension(:),intent(in) :: rval
      integer :: typeSP1_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeSP1_'

  typeSP1_=MP_REAL

end function typeSP1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP1_ - return CLMMPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeDP1_(rval)
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:),intent(in) :: rval
      integer :: typeDP1_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeDP1_'

  typeDP1_=MP_DOUBLE_PRECISION

end function typeDP1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI2_ - return CLMMPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeI2_(ival)
      implicit none
      integer,dimension(:,:),intent(in) :: ival
      integer :: typeI2_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeI2_'

  typeI2_=MP_INTEGER

end function typeI2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL2_ - return CLMMPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeL2_(lval)
      implicit none
      logical,dimension(:,:),intent(in) :: lval
      integer :: typeL2_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeL2_'

  typeL2_=MP_LOGICAL

end function typeL2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC2_ - return CLMMPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeC2_(cval)
      implicit none
      character(len=*),dimension(:,:),intent(in) :: cval
      integer :: typeC2_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeC2_'

  typeC2_=MP_CHARACTER

end function typeC2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP2_ - return CLMMPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeSP2_(rval)
      use m_realkinds,only : SP
      implicit none
      real(SP),dimension(:,:),intent(in) :: rval
      integer :: typeSP2_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeSP2_'

  typeSP2_=MP_REAL

end function typeSP2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP2_ - return CLMMPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeDP2_(rval)
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:,:),intent(in) :: rval
      integer :: typeDP2_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeDP2_'

  typeDP2_=MP_DOUBLE_PRECISION

end function typeDP2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: perr_ - CLMMPI error information hanlder
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine perr_(proc,MP_proc,ierror)
      use m_stdio, only : stderr
      implicit none
      character(len=*),intent(in) :: proc
      character(len=*),intent(in) :: MP_proc
      integer,intent(in) :: ierror

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::perr_'

  character(len=MP_MAX_ERROR_STRING) :: estr
  integer :: ln,ier

  call MP_error_string(ierror,estr,ln,ier)
  if(ier /= 0 .or. ln<=0) then
    write(stderr,'(4a,i4)') proc,': ',	&
	MP_proc,' error, ierror =',ierror
  else
    write(stderr,'(6a)') proc,': ',	&
	MP_proc,' error, "',estr(1:ln),'"'
  endif

end subroutine perr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MP_log2 - The smallest integer its power of 2 is >= nPE
!
! !DESCRIPTION:
!
! !INTERFACE:

    function MP_log2(nPE)
      implicit none
      integer,intent(in) :: nPE
      integer :: MP_log2

! !REVISION HISTORY:
! 	01Feb00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MP_log2'

  integer :: n2

  MP_log2=0
  n2=1
  do while(n2<nPE)
    MP_log2 = MP_log2+1
    n2 = n2+n2
  end do

end function MP_log2

end module m_clmmpif90
!.
