!===============================================================================
! SVN $Id: shr_clmmpi_mod.F90 3249 2007-02-21 00:45:09Z erik $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/shr/shr_clmmpi_mod.F90 $
!===============================================================================

Module shr_clmmpi_mod

!-------------------------------------------------------------------------------
! PURPOSE: general layer on CLMMPI functions
!-------------------------------------------------------------------------------

   use shr_kind_mod

   implicit none
   private

! PUBLIC: Public interfaces

   public :: shr_clmmpi_chkerr
   public :: shr_clmmpi_send
   public :: shr_clmmpi_recv
   public :: shr_clmmpi_bcast
   public :: shr_clmmpi_gathScatVInit
   public :: shr_clmmpi_gatherV
   public :: shr_clmmpi_scatterV
   public :: shr_clmmpi_sum
   public :: shr_clmmpi_min
   public :: shr_clmmpi_max
   public :: shr_clmmpi_commsize
   public :: shr_clmmpi_commrank
   public :: shr_clmmpi_initialized
   public :: shr_clmmpi_abort
   public :: shr_clmmpi_barrier
   public :: shr_clmmpi_init
   public :: shr_clmmpi_finalize

   interface shr_clmmpi_send ; module procedure &
     shr_clmmpi_sendi0, &
     shr_clmmpi_sendi1, &
     shr_clmmpi_sendr0, &
     shr_clmmpi_sendr1, &
     shr_clmmpi_sendr3
   end interface
   interface shr_clmmpi_recv ; module procedure &
     shr_clmmpi_recvi0, &
     shr_clmmpi_recvi1, &
     shr_clmmpi_recvr0, &
     shr_clmmpi_recvr1, &
     shr_clmmpi_recvr3
   end interface
   interface shr_clmmpi_bcast ; module procedure &
     shr_clmmpi_bcastc0, &
     shr_clmmpi_bcastl0, &
     shr_clmmpi_bcasti0, &
     shr_clmmpi_bcasti1, &
     shr_clmmpi_bcasti2, &
     shr_clmmpi_bcastr0, &
     shr_clmmpi_bcastr1, &
     shr_clmmpi_bcastr2, &
     shr_clmmpi_bcastr3
   end interface
   interface shr_clmmpi_gathScatVInit ; module procedure &
     shr_clmmpi_gathScatVInitr1
   end interface
   interface shr_clmmpi_gatherv ; module procedure &
     shr_clmmpi_gatherVr1
   end interface
   interface shr_clmmpi_scatterv ; module procedure &
     shr_clmmpi_scatterVr1
   end interface
   interface shr_clmmpi_sum ; module procedure &
     shr_clmmpi_sumi0, &
     shr_clmmpi_sumi1, &
     shr_clmmpi_sumr0, &
     shr_clmmpi_sumr1, &
     shr_clmmpi_sumr2, &
     shr_clmmpi_sumr3
   end interface
   interface shr_clmmpi_min ; module procedure &
     shr_clmmpi_mini0, &
     shr_clmmpi_mini1, &
     shr_clmmpi_minr0, &
     shr_clmmpi_minr1
   end interface
   interface shr_clmmpi_max ; module procedure &
     shr_clmmpi_maxi0, &
     shr_clmmpi_maxi1, &
     shr_clmmpi_maxr0, &
     shr_clmmpi_maxr1
   end interface

#include <mpif.h>         ! clmmpi library include file

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE shr_clmmpi_chkerr(rcode,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: rcode  ! input CLMMPI error code
   character(*),         intent(in) :: string ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_chkerr) '
   character(CLMMPI_MAX_ERROR_STRING)  :: lstring
   integer(SHR_KIND_IN)             :: len
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: layer on CLMMPI error checking
!-------------------------------------------------------------------------------

   if (rcode /= CLMMPI_SUCCESS) then
     call CLMMPI_ERROR_STRING(rcode,lstring,len,ierr)
     write(6,*) trim(subName),":",lstring(1:len)
     call shr_clmmpi_abort(string,rcode)
   endif

END SUBROUTINE shr_clmmpi_chkerr

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sendi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! send value
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sendi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a single integer
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_SEND(lvec,lsize,CLMMPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_sendi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sendi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sendi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call CLMMPI_SEND(lvec,lsize,CLMMPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_sendi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sendr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sendr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a real scalar
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_SEND(lvec,lsize,CLMMPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_sendr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sendr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sendr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call CLMMPI_SEND(lvec,lsize,CLMMPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_sendr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sendr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(in) :: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! clmmpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sendr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call CLMMPI_SEND(array,lsize,CLMMPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_sendr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_recvi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_recvi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(CLMMPI_STATUS_SIZE)  ! clmmpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_RECV(lvec,lsize,CLMMPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_recvi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_recvi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_recvi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(CLMMPI_STATUS_SIZE)  ! clmmpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call CLMMPI_RECV(lvec,lsize,CLMMPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_recvi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_recvr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_recvr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(CLMMPI_STATUS_SIZE)  ! clmmpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_RECV(lvec,lsize,CLMMPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_recvr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_recvr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_recvr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(CLMMPI_STATUS_SIZE)  ! clmmpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call CLMMPI_RECV(lvec,lsize,CLMMPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_recvr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_recvr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(out):: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! clmmpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_recvr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(CLMMPI_STATUS_SIZE)  ! clmmpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call CLMMPI_RECV(array,lsize,CLMMPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_recvr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcasti0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcasti0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_BCAST(vec,lsize,CLMMPI_INTEGER,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcasti0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastl0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcastl0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_BCAST(vec,lsize,CLMMPI_LOGICAL,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastl0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastc0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   character(len=*), intent(inout)    :: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcastc0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a character string
!-------------------------------------------------------------------------------

   lsize = len(vec)

   call CLMMPI_BCAST(vec,lsize,CLMMPI_CHARACTER,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastc0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastr0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcastr0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a real
!-------------------------------------------------------------------------------

   lsize = 1

   call CLMMPI_BCAST(vec,lsize,CLMMPI_REAL8,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcasti1(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcasti1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(vec)

   call CLMMPI_BCAST(vec,lsize,CLMMPI_INTEGER,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcasti1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastr1(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_bcastr1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(vec)

   call CLMMPI_BCAST(vec,lsize,CLMMPI_REAL8,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastr2(arr,comm,string)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_bcastr2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)

   call CLMMPI_BCAST(arr,lsize,CLMMPI_REAL8,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcasti2(arr,comm,string)

   IMPLICIT none

   !----- arguments -----
   integer,              intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! clmmpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_bcasti2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of integers
!-------------------------------------------------------------------------------

   lsize = size(arr)

   call CLMMPI_BCAST(arr,lsize,CLMMPI_INTEGER,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcasti2

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_bcastr3(arr,comm,string)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:,:) ! array, 3d 
   integer(SHR_KIND_IN), intent(in)   :: comm       ! clmmpi communicator
   character(*),optional,intent(in)   :: string     ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_bcastr3) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 3d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)

   call CLMMPI_BCAST(arr,lsize,CLMMPI_REAL8,0,comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_bcastr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_gathScatvInitr1(comm, rootid, locArr, glob1DArr, globSize, &
                                   displs, string )

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in)   :: comm          ! clmmpi communicator
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! CLMMPI task to gather/scatter on
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array of distributed data
   real(SHR_KIND_R8),    pointer      :: glob1DArr(:)  ! Global 1D array of gathered data
   integer(SHR_KIND_IN), pointer      :: globSize(:)   ! Size of each distributed piece
   integer(SHR_KIND_IN), pointer      :: displs(:)     ! Displacements for receive
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: npes          ! Number of CLMMPI tasks
   integer(SHR_KIND_IN)               :: locSize       ! Size of local distributed data
   integer(SHR_KIND_IN), pointer      :: sendSize(:)   ! Size to send for initial gather
   integer(SHR_KIND_IN)               :: i             ! Index
   integer(SHR_KIND_IN)               :: rank          ! Rank of this CLMMPI task
   integer(SHR_KIND_IN)               :: nSize         ! Maximum size to send
   integer(SHR_KIND_IN)               :: ierr          ! Error code
   integer(SHR_KIND_IN)               :: nSiz1D        ! Size of 1D global array
   integer(SHR_KIND_IN)               :: maxSize       ! Maximum size

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_gathScatvInitr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Setup arrays for a gatherv/scatterv operation
!-------------------------------------------------------------------------------

   locSize = size(locarr)
   call shr_clmmpi_commsize( comm, npes )
   call shr_clmmpi_commrank( comm, rank )
   allocate( globSize(npes) )
   !
   ! --- Gather the send global sizes from each CLMMPI task -----------------------
   !
   allocate( sendSize(npes) )
   sendSize(:) = 1
   call CLMMPI_GATHER( locSize, 1, CLMMPI_INTEGER, globSize, sendSize, &
                    CLMMPI_INTEGER, rootid, comm, ierr )
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif
   deallocate( sendSize )
   !
   ! --- Prepare the displacement and allocate arrays -------------------------
   !
   allocate( displs(npes) )
   displs(1) = 0
   if ( rootid /= rank )then
      maxSize = 1
   else
      maxSize = maxval(globSize)
   end if
   nsiz1D    = min(maxSize,globSize(1))
   do i = 2, npes
      nSize = min(maxSize,globSize(i-1))
      displs(i) = displs(i-1) + nSize
      nsiz1D = nsiz1D + min(maxSize,globSize(i))
   end do
   allocate( glob1DArr(nsiz1D) )
   !----- Do some error checking for the root task arrays computed ----
   if ( rootid == rank )then
      if ( nsiz1D /= sum(globSize) ) &
         call shr_clmmpi_abort( subName//" : Error, size of global array not right" )
      if ( any(displs < 0) .or. any(displs >= nsiz1D) ) &
         call shr_clmmpi_abort( subName//" : Error, displacement array not right" )
      if ( (displs(npes)+globSize(npes)) /= nsiz1D ) &
         call shr_clmmpi_abort( subName//" : Error, displacement array values too big" )
   end if

END SUBROUTINE shr_clmmpi_gathScatvInitr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_gathervr1(locarr, locSize, glob1DArr, globSize, displs, rootid, &
                             comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array
   real(SHR_KIND_R8),    intent(inout):: glob1DArr(:)  ! Global 1D array to receive in on
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to send this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to receive each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for receive
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! CLMMPI task to gather on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! clmmpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_gathervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Gather a 1D array of reals
!-------------------------------------------------------------------------------

   call CLMMPI_GATHERV( locarr, locSize, CLMMPI_REAL8, glob1Darr, globSize, displs, &
                     CLMMPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_gathervr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_scattervr1(locarr, locSize, glob1Darr, globSize, displs, rootid, &
                              comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(out)  :: locarr(:)     ! Local array
   real(SHR_KIND_R8),    intent(in)   :: glob1Darr(:)  ! Global 1D array to send from
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to receive this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to send to each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for send
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! CLMMPI task to scatter on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! clmmpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_clmmpi_scattervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Scatter a 1D array of reals
!-------------------------------------------------------------------------------


   call CLMMPI_SCATTERV( glob1Darr, globSize, displs, CLMMPI_REAL8, locarr, locSize, &
                      CLMMPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_scattervr1


!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumr2(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:)! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:)! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumr2) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_sumr3(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:,:) ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:,:) ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_sumr3) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_sumr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_mini0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_mini0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_mini0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_mini1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_mini1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_mini1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_minr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_minr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_minr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_minr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_minr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_minr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_maxi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_maxi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_maxi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_maxi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_maxi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_maxi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_maxr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_maxr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_maxr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_maxr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! clmmpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_clmmpi_maxr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! clmmpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = CLMMPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_clmmpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call CLMMPI_ALLREDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_ALLREDUCE")
   else
     call CLMMPI_REDUCE(lvec,gvec,gsize,CLMMPI_REAL8,reduce_type,0,comm,ierr)
     call shr_clmmpi_chkerr(ierr,trim(lstring)//" CLMMPI_REDUCE")
   endif

END SUBROUTINE shr_clmmpi_maxr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_commsize(comm,size,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: size
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_commsize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI commsize
!-------------------------------------------------------------------------------

   call CLMMPI_COMM_SIZE(comm,size,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_commsize

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_commrank(comm,rank,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: rank
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_commrank) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI commrank
!-------------------------------------------------------------------------------

   call CLMMPI_COMM_RANK(comm,rank,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_commrank

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_initialized(flag,string)

   IMPLICIT none

   !----- arguments ---
   logical,intent(out)                :: flag
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_initialized) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI initialized
!-------------------------------------------------------------------------------

   call CLMMPI_INITIALIZED(flag,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_initialized

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_abort(string,rcode)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message
   integer,optional,intent(in)        :: rcode    ! optional code

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_abort) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI abort
!-------------------------------------------------------------------------------

   if ( present(string) .and. present(rcode) ) write(6,*) trim(subName),":",trim(string),rcode
   call CLMMPI_ABORT(CLMMPI_COMM_WORLD,rcode,ierr)

END SUBROUTINE shr_clmmpi_abort

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_barrier(comm,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_barrier) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI barrier
!-------------------------------------------------------------------------------

   call CLMMPI_BARRIER(comm,ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_barrier

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_init(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_init) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI init
!-------------------------------------------------------------------------------

   call CLMMPI_INIT(ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_init

!===============================================================================
!===============================================================================

SUBROUTINE shr_clmmpi_finalize(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_clmmpi_finalize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: CLMMPI finalize
!-------------------------------------------------------------------------------

   call CLMMPI_BARRIER(CLMMPI_COMM_WORLD,ierr)
   call CLMMPI_FINALIZE(ierr)
   if (present(string)) then
     call shr_clmmpi_chkerr(ierr,subName//trim(string))
   else
     call shr_clmmpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_clmmpi_finalize

!===============================================================================
!===============================================================================

END MODULE shr_clmmpi_mod
