#include <misc.h>
#include <preproc.h>

module spmdMod

!----------------------------------------------------------------------- 
! 
! Purpose: 
! MPI routines for initialization and computing arguments for
! gatherv and scatterv operations
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined COUP_CSM)
  use shr_msg_mod
#endif

#if (!defined SPMD)
! ----------------------- SPMD OFF --------------------------------------

  logical :: masterproc = .true. ! proc 0 logical for printing msgs
  integer :: iam = 0
  integer :: npes = 1
#endif

#if (defined SPMD)
! ----------------------- SPMD ON ---------------------------------------

#if (defined COUP_CAM)
  use mpishorthand
  use spmd_dyn, only: npes
  use pmgrid  , only: masterproc, iam 
#endif 

#if (defined OFFLINE)
  use mpiinc
#endif

#if (defined OFFLINE) || (defined COUP_CSM)
  integer, public :: npes        !number of processors
  integer, public :: iam         !proc number
  logical, public :: masterproc  !proc 0 logical for printing msgs
  integer, public :: mpicom 
#endif

  integer, public, allocatable :: proc_gridtot(:)
  integer, public, allocatable :: proc_landtot(:)
  integer, public, allocatable :: proc_coltot(:)
  integer, public, allocatable :: proc_pfttot(:)

  public :: scatter_data_from_master, gather_data_to_master

  INTERFACE scatter_data_from_master
     MODULE procedure scatter_1darray_int
     MODULE procedure scatter_1darray_real
     MODULE procedure scatter_2darray_int
     MODULE procedure scatter_2darray_real
  END INTERFACE

  INTERFACE gather_data_to_master
     MODULE procedure gather_1darray_int
     MODULE procedure gather_1darray_real
     MODULE procedure gather_2darray_int
     MODULE procedure gather_2darray_real
  END INTERFACE

  SAVE

!===============================================================================
CONTAINS
!===============================================================================

#if (defined OFFLINE) || (defined COUP_CSM)

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init
!
! !INTERFACE:
  subroutine spmd_init
!
! !DESCRIPTION: 
! MPI initialization (number of cpus, processes, tids, etc)
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i,j        ! indices
    integer ier        ! return error status      
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(MPI_MAX_PROCESSOR_NAME), allocatable :: proc_name(:)
#if (defined OFFLINE)
    logical mpi_running
#endif
!-----------------------------------------------------------------------

#if (defined OFFLINE)

    ! Initialize mpi

    call mpi_initialized (mpi_running, ier)
    if (.not. mpi_running) call mpi_init (ier)       

#endif

    ! Set communication group 

#if (defined OFFLINE)
    mpicom  = MPI_COMM_WORLD
#elif (defined COUP_CSM)
    mpicom  = SHR_MSG_COMM_LND
#endif

    ! Get my processor id  

    call mpi_comm_rank(mpicom, iam, ier)  
    if (iam==0) then 
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors

    call mpi_comm_size(mpicom, npes, ier) 

    ! Get my processor names

    allocate (length(0:npes-1))
    allocate (displ(0:npes-1))
    allocate (proc_name(0:npes-1))
    call mpi_get_processor_name (proc_name(iam),length(iam),ier)
    call mpi_allgather (length(iam),1,MPI_INTEGER,length,1,MPI_INTEGER,mpicom,ier)
    do i =0,npes-1
       displ(i)=i*MPI_MAX_PROCESSOR_NAME
    end do
    call mpi_gatherv (proc_name(iam),length(iam),MPI_CHARACTER, &
                      proc_name,length,displ,MPI_CHARACTER,0,mpicom,ier)
    if (masterproc) then
       write(6,100)npes
       write(6,200)
       write(6,220)
       do i=0,npes-1
          write(6,250)i,(proc_name((i))(j:j),j=1,length(i))
       end do
    endif
    deallocate (length)
    deallocate (displ)
    deallocate (proc_name)

100 format(i3," pes participating in computation")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i3,")",2x,100a1)

    return
  end subroutine spmd_init

#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init_arrays
!
! !INTERFACE:
  subroutine spmd_init_arrays
!
! !DESCRIPTION: 
! Initialize arrays for number of land/patch points per proc
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    allocate (proc_gridtot(0:npes-1))
    allocate (proc_landtot(0:npes-1))
    allocate (proc_coltot(0:npes-1))
    allocate (proc_pfttot(0:npes-1))
    
    return
  end subroutine spmd_init_arrays

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_compute_mpigs
!
! !INTERFACE:
  subroutine spmd_compute_mpigs (clmlevel, nfact, numtot, numperproc, displs, indexi)
!
! !DESCRIPTION: 
! Compute arguments for gatherv, scatterv for vectors
!
! !USES:
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: clmlevel     !type of input data 
    integer, intent(in ) :: nfact                !multiplicative factor for patches
    integer, intent(out) :: numtot               !total number of elements (to send or recv)
    integer, intent(out) :: numperproc(0:npes-1) !per-PE number of items to receive
    integer, intent(out) :: displs(0:npes-1)     !per-PE displacements
    integer, intent(out) :: indexi               !beginning array index (grid,land,col or pft)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p                                 ! index
!-----------------------------------------------------------------------

    select case (clmlevel)
    case('gridcell')
       numtot = (proc_gridtot(iam))*nfact
       do p=0,npes-1
          numperproc(p) = proc_gridtot(p)*nfact
       end do
       indexi = grid1d%beg
    case('landunit')
       numtot = (proc_landtot(iam))*nfact
       do p=0,npes-1
          numperproc(p) = proc_landtot(p)*nfact
       end do
       indexi = land1d%beg
    case('column')
       numtot = (proc_coltot(iam))*nfact
       do p=0,npes-1
          numperproc(p) = proc_coltot(p)*nfact
       end do
       indexi = cols1d%beg
    case('pft')
       numtot = (proc_pfttot(iam))*nfact
       do p=0,npes-1
          numperproc(p) = proc_pfttot(p)*nfact
       end do
       indexi = pfts1d%beg
    case default
       write(6,*) 'COMPUTE_MPIGS: Invalid expansion character: ',trim(clmlevel)
       call endrun
    end select

    displs(0) = 0
    do p=1,npes-1
       displs(p) = displs(p-1) + numperproc(p-1)
    end do

    return
  end subroutine spmd_compute_mpigs

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_int
!
! !INTERFACE:
  subroutine scatter_1darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to scatter integer 1d array 
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:) :: ilocal(:)   !local read data
    integer, pointer, dimension(:) :: iglobal(:)  !global read data
    character(len=*), intent(in) :: clmlevel      !type of input data 
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: beg                    !temporary
    integer :: numsendv(0:npes-1)     !vector of items to be sent
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numrecv                !number of items to be received
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (iglobal, numsendv, displsv, MPI_INTEGER, &
            ilocal(beg), numrecv , MPI_INTEGER , 0, mpicom, ier)
    else
       call mpi_scatterv (0, numsendv, displsv, MPI_INTEGER, &
            ilocal(beg), numrecv , MPI_INTEGER , 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_1darray_int error: ',ier
       call endrun
    endif
    return
  end subroutine scatter_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_real
!
! !INTERFACE:
  subroutine scatter_1darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to scatter 1d real array from master processor 
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:) :: rlocal  !local read data
    real(r8), pointer, dimension(:) :: rglobal !global read data
    character(len=*), intent(in) :: clmlevel   !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                             !return code
    integer :: beg                             !temporaries
    integer :: numsendv(0:npes-1)              !vector of items to be sent
    integer :: displsv(0:npes-1)               !displacement vector
    integer :: numrecv                         !number of items to be received
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (rglobal, numsendv, displsv, MPI_REAL8, &
            rlocal(beg), numrecv , MPI_REAL8 , 0, mpicom, ier)
    else
       call mpi_scatterv (0._r8, numsendv, displsv, MPI_REAL8, &
            rlocal(beg), numrecv , MPI_REAL8 , 0, mpicom, ier)
    endif
    if (ier/=0 ) then
       write(6,*)'scatter_1darray_real error: ',ier
       call endrun
    endif
    return
  end subroutine scatter_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_2darray_int
!
! !INTERFACE:
  subroutine scatter_2darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to scatter 2d integer array from master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:,:) :: ilocal  !local read data
    integer, pointer, dimension(:,:) :: iglobal !global read data
    character(len=*), intent(in) :: clmlevel    !type of input data 
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                   !error status         
    integer :: ndim1                 !size of first dimension
    integer :: l1                    !lower bound of first dimension
    integer :: beg                   !temporaries                 
    integer :: numsendv(0:npes-1)    !vector of items to be sent
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(ilocal,dim=1) /= lbound(iglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(ilocal,dim=1)
	  write(6,*)' l1 global= ',lbound(iglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(ilocal,dim=1)
    ndim1 = size(ilocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (iglobal(l1,1), numsendv, displsv, MPI_INTEGER, &
            ilocal(l1,beg), numrecv ,  MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_scatterv (0, numsendv, displsv, MPI_INTEGER, &
            ilocal(l1,beg), numrecv ,  MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_2darray_int error: ',ier
       call endrun
    endif
    return
  end subroutine scatter_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_2darray_real
!
! !INTERFACE:
  subroutine scatter_2darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to scatter 2d integer array from master processor
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:,:) :: rlocal  !local read data
    real(r8), pointer, dimension(:,:) :: rglobal !global read data
    character(len=*), intent(in) :: clmlevel     !type of input data 
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            !return code 
    integer :: ndim1                          !size of second dimension         
    integer :: l1                             !lower,upper bounds of input arrays
    integer :: beg                            !temporaries                 
    integer :: numsendv(0:npes-1)             !vector of items to be sent
    integer :: displsv(0:npes-1)              !displacement vector
    integer :: numrecv                        !number of items to be received
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(rlocal,dim=1) /= lbound(rglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(rlocal,dim=1)
	  write(6,*)' l1 global= ',lbound(rglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(rlocal,dim=1)
    ndim1 = size(rlocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (rglobal(l1,1), numsendv, displsv, MPI_REAL8, &
            rlocal(l1,beg), numrecv ,  MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_scatterv (0._r8, numsendv, displsv, MPI_REAL8, &
            rlocal(l1,beg), numrecv ,  MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_2darray_real error: ',ier
       call endrun
    endif
    return
  end subroutine scatter_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_int
!
! !INTERFACE:
  subroutine gather_1darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to gather 1d integer array on master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:) :: ilocal     !output data
    integer, pointer, dimension(:) :: iglobal    !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !errorcode
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received  
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (ilocal(beg), numsend , MPI_INTEGER, &
            iglobal, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_gatherv (ilocal(beg), numsend , MPI_INTEGER, &
            0, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier/=0 ) then
       write(6,*)'gather_1darray_int error: ',ier
       call endrun
    endif
    return
  end subroutine gather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_real
!
! !INTERFACE:
  subroutine gather_1darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to gather 1d real array on master processor
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:) :: rlocal    !output data
    real(r8), pointer, dimension(:) :: rglobal   !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received  
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (rlocal(beg), numsend , MPI_REAL8, &
            rglobal, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_gatherv (rlocal(beg), numsend , MPI_REAL8, &
            0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_1darray_real error: ',ier
       call endrun
    endif
    return
  end subroutine gather_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_int
!
! !INTERFACE:
  subroutine gather_2darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to gather 2d integer array on master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:,:) :: ilocal   !read data
    integer, pointer, dimension(:,:) :: iglobal  !global data
    character(len=*), intent(in) :: clmlevel     !type of input data 
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                   !return code
    integer :: ndim1                 !size of second dimension         
    integer :: l1                    !lower bounds of input arrays
    integer :: beg                   !temporaries                 
    integer :: numrecvv(0:npes-1)    !vector of items to be received  
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(ilocal,dim=1) /= lbound(iglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(ilocal,dim=1)
	  write(6,*)' l1 global= ',lbound(iglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(ilocal,dim=1)
    ndim1 = size(ilocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (ilocal(l1,beg), numsend , MPI_INTEGER, &
            iglobal(l1,1), numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_gatherv (ilocal(l1,beg), numsend , MPI_INTEGER, &
            0, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_2darray_int error: ',ier
       call endrun
    endif
    return
  end subroutine gather_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_real
!
! !INTERFACE:
  subroutine gather_2darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to gather 2d real array
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:,:) :: rlocal   !local data
    real(r8), pointer, dimension(:,:) :: rglobal  !global data
    character(len=*), intent(in) :: clmlevel      !type of input data 
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: ndim1                  !size of second dimension         
    integer :: l1                     !lower bounds of input arrays
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received  
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(rlocal,dim=1) /= lbound(rglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(rlocal,dim=1)
	  write(6,*)' l1 global= ',lbound(rglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(rlocal,dim=1)
    ndim1 = size(rlocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (rlocal(l1,beg), numsend , MPI_REAL8, &
            rglobal(l1,1), numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_gatherv (rlocal(l1,beg), numsend , MPI_REAL8, &
            0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_2darray_real error: ',ier
       call endrun
    endif
    return
  end subroutine gather_2darray_real

#endif

end module spmdMod


