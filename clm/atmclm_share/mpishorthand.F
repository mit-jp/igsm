#include <misc.h>
!-----------------------------------------------------------------------------------
!
! Purpose: 
!
! 	Data and parameters used for MPI. Some shorthand variables with shorter
!	names than the standard MPI parameters. Also some variables used for heap
!	management.
!
! Note: The #include of "mpif.h" which is typically in f77 fixed format means that
! 	 this module MUST be in fixed format.
!
! Author: Jim Rosinski
!
! $Id$
!
!-----------------------------------------------------------------------------------
      module mpishorthand
      public
#if (defined SPMD || defined COUP_CSM)
#include <mpif.h>

!
! Need to set as variables rather than parameters since some MPI implementations 
! set values for MPI tags at run time
!
      save
      integer, public :: mpiint    ! MPI integers
      integer, public :: mpichar   ! MPI character data
      integer, public :: mpilog    ! MPI logical data
      integer, public :: mpir4     ! MPI real data for r4
      integer, public :: mpir8     ! MPI real data
      integer, public :: mpicom    ! MPI communication
      integer, public :: mpipk     ! MPI packed data
!
! Common info for heap manager
!
      integer, public::  nsend = 0  ! Number of MPI messages sent
      integer, public::  nrecv = 0  ! Number of MPI messages received
      integer, public::  nwsend = 0 ! Number of MPI words sent
      integer, public::  nwrecv = 0 ! Number of MPI words received

#endif
      end module mpishorthand
