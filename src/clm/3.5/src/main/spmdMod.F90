#include <misc.h>
#include <preproc.h>

module spmdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdMod
!
! !DESCRIPTION:
! SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private

  ! Default settings valid even if there is no spmd 

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors for clm
  integer, public :: clmmpicom          ! communicator group for clm
  integer, public :: comp_id         ! component id

  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

  !
  ! Values from clmmpif.h that can be used
  !
  public :: CLMMPI_INTEGER
  public :: CLMMPI_REAL8
  public :: CLMMPI_LOGICAL
  public :: CLMMPI_SUM
  public :: CLMMPI_MIN
  public :: CLMMPI_MAX
  public :: CLMMPI_STATUS_SIZE
  public :: CLMMPI_ANY_SOURCE
  public :: CLMMPI_CHARACTER
  public :: CLMMPI_COMM_WORLD
  public :: CLMMPI_MAX_PROCESSOR_NAME

#include <mpif.h>  

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init( clm_clmmpicom )
!
! !INTERFACE:
  subroutine spmd_init( clm_clmmpicom )
!
! !DESCRIPTION:
! CLMMPI initialization (number of cpus, processes, tids, etc)
!
! !USES
#if (defined COUP_CSM)
    use cpl_comm_mod, only : cpl_comm_mph_cid
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: clm_clmmpicom
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         ! indices
    integer :: ier         ! return error status
    logical :: clmmpi_running ! temporary
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(CLMMPI_MAX_PROCESSOR_NAME), allocatable :: procname(:)
!-----------------------------------------------------------------------

    ! Initialize clmmpi communicator group

    clmmpicom = clm_clmmpicom

    comp_id = 1
#if (defined COUP_CSM)
    comp_id = cpl_comm_mph_cid
#endif

    ! Initialize clmmpi

#ifdef OFFLINE
    call clmmpi_initialized (clmmpi_running, ier)
    if (.not. clmmpi_running) call clmmpi_init (ier)
#endif

    ! Get my processor id

    call clmmpi_comm_rank(clmmpicom, iam, ier)
    if (iam==0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors

    call clmmpi_comm_size(clmmpicom, npes, ier)

    ! Get my processor names

    allocate (length(0:npes-1), displ(0:npes-1), procname(0:npes-1))

    call clmmpi_get_processor_name (procname(iam), length(iam), ier)
    call clmmpi_allgather (length(iam),1,CLMMPI_INTEGER,length,1,CLMMPI_INTEGER,clmmpicom,ier)
    do i = 0,npes-1
       displ(i)=i*CLMMPI_MAX_PROCESSOR_NAME
    end do
    call clmmpi_gatherv (procname(iam),length(iam),CLMMPI_CHARACTER, &
                      procname,length,displ,CLMMPI_CHARACTER,0,clmmpicom,ier)
    if (masterproc) then
       write(6,100)npes
       write(6,200)
       write(6,220)
       do i=0,npes-1
          write(6,250)i,(procname((i))(j:j),j=1,length(i))
       end do
    endif

    deallocate (length, displ, procname)

100 format(i3," pes participating in computation")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i3,")",2x,100a1)

  end subroutine spmd_init

end module spmdMod
