#include <misc.h>
#include <preproc.h>

module iobinary

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iobinary
! 
! !DESCRIPTION: 
! Set of wrappers to write binary I/O
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined SPMD)
  use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master, &
       mpicom, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL 
#else
  use spmdMod, only : masterproc
#endif
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: readin, wrtout  ! read and write bindary I/O
  interface readin
     module procedure readin_1darray_int
     module procedure readin_2darray_int
     module procedure readin_1darray_real
     module procedure readin_2darray_real
  end interface
  interface wrtout
     module procedure wrtout_1darray_int
     module procedure wrtout_2darray_int
     module procedure wrtout_1darray_real
     module procedure wrtout_2darray_real
  end interface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: get_num1d  !get 1d type
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_1d_array_int
!
! !INTERFACE:
  subroutine readin_1darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to read integer 1d array from restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    integer, pointer, dimension(:) :: iarr     !input data
    character(len=*), intent(in) :: clmlevel   !type of input data 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                             !return code
#if (defined SPMD)
    integer :: nsize                           !size of first dimension 
    integer, pointer, dimension(:) :: iglobal  !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (iglobal(nsize))
       read (iu,iostat=ier) iglobal
       if (ier /= 0 ) then
          write (6,*) 'READIN_1DARRAY_INT error ',ier,' on i/o unit = ',iu
          call endrun
       endif
    endif
    call scatter_data_from_master(iarr, iglobal, clmlevel=clmlevel)
    if (masterproc) then
       deallocate(iglobal)
    endif
#else
    read (iu,iostat=ier) iarr
    if (ier /= 0 ) then
       write (6,*) 'READIN_1DARRAY_INT error ',ier,' on i/o unit = ',iu
       call endrun
    endif
#endif
    return
  end subroutine readin_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_1darray_real
!
! !INTERFACE:
  subroutine readin_1darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to read real 1d array from restart binary file 

! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    real(r8), pointer, dimension(:) :: rarr    !input data
    character(len=*), intent(in) :: clmlevel   !input data type
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                              !return code
#if (defined SPMD)
    integer :: nsize                            !size of first dimension 
    real(r8), pointer, dimension(:) :: rglobal  !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (rglobal(nsize))
       read (iu,iostat=ier) rglobal
       if (ier /= 0 ) then
          write (6,*) 'READIN_1DARRAY_REAL error ',ier,' on i/o unit = ',iu
          call endrun
       endif
    endif
    call scatter_data_from_master(rarr, rglobal, clmlevel=clmlevel)
    if (masterproc) then
       deallocate(rglobal)
    endif
#else
    read (iu,iostat=ier) rarr
    if (ier /= 0 ) then
       write (6,*) 'READIN_1DARRAY_REAL error ',ier,' on i/o unit = ',iu
       call endrun
    endif
#endif
    return
  end subroutine readin_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_2d_arrayint
!
! !INTERFACE:
  subroutine readin_2darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to read integer 2d array from restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    integer, pointer, dimension(:,:) :: iarr   !input data
    character(len=*), intent(in) :: clmlevel   !type of input data 
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                               !error status         
#if (defined SPMD)
    integer :: lb,ub                             !bounds of first dimension 
    integer :: nsize                             !size of second dimension         
    integer, pointer, dimension(:,:) :: iglobal  !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       lb = lbound(iarr, dim=1)
       ub = ubound(iarr, dim=1)
       allocate (iglobal(lb:ub,nsize))
       read (iu,iostat=ier) iglobal
       if (ier /= 0 ) then
          write (6,*) 'READIN_2DARRAY_INT error ',ier,' on i/o unit = ',iu
          call endrun
       endif
    endif
    call scatter_data_from_master(iarr, iglobal, clmlevel=clmlevel )
    if (masterproc) then
       deallocate(iglobal)
    endif
#else
    read (iu,iostat=ier) iarr
    if (ier /= 0 ) then
       write (6,*) 'READIN_2DARRAY_INT error ',ier,' on i/o unit = ',iu
       call endrun
    endif
#endif
    return
  end subroutine readin_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_2darray_real
!
! !INTERFACE:
  subroutine readin_2darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to read real 2d array from restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !input unit
    real(r8), pointer, dimension(:,:) :: rarr !input data
    character(len=*), intent(in) :: clmlevel  !type of input data 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein  
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                               !return code 
#if (defined SPMD)
    integer :: lb,ub                             !bounds of first dimension 
    integer :: nsize                             !size of second dimension         
    real(r8), pointer, dimension(:,:) :: rglobal !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       lb = lbound(rarr, dim=1)
       ub = ubound(rarr, dim=1)
       allocate (rglobal(lb:ub,nsize))
       read (iu,iostat=ier) rglobal
       if (ier /= 0 ) then
          write (6,*) 'READIN_2DARRAY_REAL error ',ier,' on i/o unit = ',iu
          call endrun
       endif
    endif
    call scatter_data_from_master(rarr, rglobal, clmlevel=clmlevel)
    if (masterproc) then
       deallocate(rglobal)
    endif
#else
    read (iu,iostat=ier) rarr
    if (ier /= 0 ) then
       write (6,*) 'READIN_2DARRAY_REAL error ',ier,' on i/o unit = ',iu
       call endrun
    endif
#endif
    return
  end subroutine readin_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_1d_array_int
!
! !INTERFACE:
  subroutine wrtout_1darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to write integer array to restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !output unit
    integer, pointer, dimension(:) :: iarr    !output data
    character(len=*), intent(in) :: clmlevel  !output 1d type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                              !errorcode
#if (defined SPMD)
    integer :: nsize                            !size of first dimension 
    integer, pointer, dimension(:) :: iglobal   !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (iglobal(nsize))
    endif
    call gather_data_to_master(iarr, iglobal, clmlevel=clmlevel)
    if (masterproc) then
       write (iu, iostat=ier) iglobal
       if (ier /= 0 ) then
          write(6,*) 'WRTOUT_1DARRAY_INT error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(iglobal)
    endif
#else
    write (iu, iostat=ier) iarr
    if (ier /= 0 ) then
       write(6,*) 'WRTOUT_1DARRAY_INT error ',ier,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_1d_array_real
!
! !INTERFACE:
  subroutine wrtout_1darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to write real array to restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !output unit
    real(r8), pointer, dimension(:) :: rarr   !output data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
#if (defined SPMD)
    integer :: nsize                            !size of first dimension 
    real(r8), pointer, dimension(:) :: rglobal  !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (rglobal(nsize))
    endif
    call gather_data_to_master(rarr, rglobal, clmlevel=clmlevel)
    if (masterproc) then
       write (iu, iostat=ier) rglobal
       if (ier /= 0 ) then
          write(6,*) 'WRTOUT_1DARRAY_REAL error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(rglobal)
    endif
#else
    write (iu, iostat=ier) rarr
    if (ier /= 0 ) then
       write(6,*) 'WRTOUT_1DARRAY_REAL error ',ier,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_2d_array_int
!
! !INTERFACE:
  subroutine wrtout_2darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to write integer array to restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                !output unit
    integer, pointer, dimension(:,:) :: iarr !output data
    character(len=*), intent(in) :: clmlevel !output data 1d type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                              !return code
#if (defined SPMD)
    integer :: lb,ub                             !bounds of first dimension 
    integer :: nsize                             !size of second dimension         
    integer, pointer, dimension(:,:) :: iglobal !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       lb = lbound(iarr, dim=1)
       ub = ubound(iarr, dim=1)
       allocate (iglobal(lb:ub,nsize))
    endif
    call gather_data_to_master(iarr, iglobal, clmlevel=clmlevel)
    if (masterproc) then
       write (iu, iostat=ier) iglobal
       if (ier /= 0 ) then
          write(6,*) 'WRTOUT_2DARRAY_INT error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(iglobal)
    endif
#else
    write (iu, iostat=ier) iarr
    if (ier /= 0 ) then
       write(6,*) 'WRTOUT_2DARRAY_INT error ',ier,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout_2darray_int

!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_2darray_real
!
! !INTERFACE:
  subroutine wrtout_2darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION: 
! Wrapper routine to write real array to restart binary file 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !input unit
    character(len=*), intent(in) :: clmlevel  !type of input data 
    real(r8), pointer, dimension(:,:) :: rarr !output data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                               !return code
#if (defined SPMD)
    integer :: lb,ub                             !bounds of first dimension 
    integer :: nsize                             !size of second dimension         
    real(r8), pointer, dimension(:,:) :: rglobal !temporary 
#endif
!-----------------------------------------------------------------------
#if (defined SPMD)
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       lb = lbound(rarr, dim=1)
       ub = ubound(rarr, dim=1)
       allocate (rglobal(lb:ub,nsize))
    endif
    call gather_data_to_master(rarr, rglobal, clmlevel=clmlevel)
    if (masterproc) then
       write (iu, iostat=ier) rglobal
       if (ier /= 0 ) then
          write(6,*) 'WRTOUT_2DARRAY_REAL error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(rglobal)
    endif
#else
    write (iu, iostat=ier) rarr
    if (ier /= 0 ) then
       write(6,*) 'WRTOUT_2DARRAY_REAL error ',ier,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_num1d
!
! !INTERFACE:
  integer function get_num1d (clmlevel)
!
! !DESCRIPTION: 
! Determines 1d size from clmlevel type
!
! !USES: 
  use clmtype
!
! !ARGUMENTS:
  implicit none
  character(len=*), intent(in) :: clmlevel      !type of clm 1d array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    select case (clmlevel)
    case('gridcell')
       get_num1d = grid1d%num
    case('landunit')
       get_num1d = land1d%num
    case('column')
       get_num1d = cols1d%num
    case('pft')
       get_num1d = pfts1d%num
    case default
       write(6,*) 'GET1DSIZE does not match level: ',trim(clmlevel)
       call endrun
    end select
    return
  end function get_num1d

!=======================================================================

end module iobinary
