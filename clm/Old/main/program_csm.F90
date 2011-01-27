#include <misc.h>
#include <preproc.h>

#if (defined COUP_CSM)

!================================================================================
! CCSM was developed in cooperation with the National Science Foundation,
! the Department of Energy Los Alamos National Laboratory (LANL), the
! National Aeronautics and Space Administration Data Assimilation Office,
! and the University Corporation for Atmospheric Research's National
! Center for Atmospheric Research.*  Except for POP, SCRIP, and CICE,
! which are segregable components of the model, CCSM is public domain
! software.
! 
! POP, SCRIP, and CICE were developed at LANL and are protected by the
! following copyright notice.
! 
! Copyright ©2002 The Regents of the University of California. All Rights
! Reserved.
! 
! Permission to use, copy, modify, and distribute this software and its
! documentation for educational, research and non-profit purposes, without
! fee, and without a written agreement is hereby granted, provided that
! the above copyright notice, this paragraph and the following three
! paragraphs appear in all copies.
! 
! Permission to incorporate this software into commercial products may be
! obtained by contacting the University of California, Charles Rzeszutko,
! Campus Liaison Officer Office of Technology Transfer, 1111 Franklin
! Street, 5th Floor, Oakland, CA 94607, (510) 587-6063,
! charles.rzeszutko@ucop.edu.
! 
! THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE COPYRIGHTED BY
! THE REGENTS OF THE UNIVERSITY OF CALIFORNIA. THE SOFTWARE
! PROGRAM AND DOCUMENTATION ARE SUPPLIED "AS IS", WITHOUT ANY
! ACCOMPANYING SERVICES FROM THE REGENTS. THE REGENTS DOES NOT
! WARRANT THAT THE OPERATION OF THE PROGRAM WILL BE
! UNINTERRUPTED OR ERROR-FREE. THE END-USER UNDERSTANDS THAT
! THE PROGRAM WAS DEVELOPED FOR RESEARCH PURPOSES AND IS
! ADVISED NOT TO RELY EXCLUSIVELY ON THE PROGRAM FOR ANY
! REASON.
! 
! IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
! PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
! CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
! THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
! UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
! SUCH DAMAGE.  THE UNIVERSITY OF CALIFORNIA SPECIFICALLY
! DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN
! "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO
! OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
! ENHANCEMENTS, OR MODIFICATIONS.
! 
! *The National Center for Atmospheric Research is sponsored by the
! National Science Foundation.  Any opinions, findings and conclusions or
! recommendations expressed in this publication are those of the author(s)
! and do not necessarily reflect the views of the National Science
! Foundation.
!================================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: program_csm
!
! !INTERFACE:
PROGRAM program_csm
!
! !DESCRIPTION: 
! Driver for CLM as the land component of CCSM.
! This program is the driver for CLM to work as the land component of
! CCSM.  The flux coupler will provide all the appropriate atmospheric
! forcing for the land model to run.
! o the land surface model returns to the CCSM flux coupler surface
!   fluxes, temperatures, and albedos for only the land points on the
!   [lsmlon x lsmlat] grid 
! o the land surface model uses its own surface type data set. because
!   it processes only land points, this data set must correctly define
!   what points are land and what are not land
! o the land surface model uses its own grid dimensions (lsmlon and
!   lsmlat). currently these must equal lsmlon and lsmlat so that there
!   is a direct correspondence between the atmosphere and land grids
! o the zenith angle calculation is calculated for the
!   NEXT time step rather than the current time step. make sure
!   the calendar day is for the NEXT time step. make sure the
!   solar declination calculation is the same as in the 
!   atmospheric model but for the NEXT time step. make sure the
!   calendar day is for greenwich time (see next comment).
! o subroutine calendr: this generates a julian day (with fraction)
!   based on the time step, which is used to calculate the solar
!   zenith angle. this time must be at the greenwich meridian to
!   get the correct zenith angle. also, output from this subroutine 
!   is used to calculate the month (1, ..., 12), day (1, ..., 31), 
!   and year (00, ...) of the simulation. 
! o the land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. the net longwave radiation
!   at the surface will differ somewhat from that calculated from the
!   CCSM flux coupler because the cpl model will use the upward 
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. the net solar radiation should equal that calculated
!   from the flux coupler. if not, there is a problem.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar         !parameters
  use clm_varctl         !run control variables    
  use shr_orb_mod        !orbital parameters and routines 
  use shr_msg_mod        !csm message passing routines and variables    
#if (defined SPMD)
  use spmdMod      , only : masterproc, iam, spmd_init
#else
  use spmdMod      , only : masterproc, iam
#endif
  use initializeMod, only : initialize
  use clm_csmMod   , only : csmstop_now
  use time_manager , only : advance_timestep, get_nstep 
!
! !ARGUMENTS:
    implicit none
#include "gpt.inc"
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: i,j        !loop indices 	
  integer :: nstep      !time step 
  logical :: doalb      !true if surface albedo calculation time step
!
! Earth's orbital characteristics
!
  integer  :: iyear_AD  !Year (AD) to simulate above earth's orbital parameters for
  real(r8) :: eccen     !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) :: obliq     !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) :: mvelp     !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)
!
! Orbital information after call to routine shr_orbit_params
!
  real(r8) :: obliqr    !Earth's obliquity in radians
  real(r8) :: lambm0    !Mean longitude (radians) of perihelion at the vernal equinox 
  real(r8) :: mvelpp    !Earth's moving vernal equinox longitude
  logical  :: log_print !true=> print diagnostics  
!-----------------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! Initialize model
  ! -----------------------------------------------------------------
  
  ! Initialize timing library.  2nd arg 0 means disable, 1 means enable
  
  call t_setoptionf (usrsys, 0)
  call t_initializef ()

  ! Determine input/output units 

  call shr_msg_stdio ('lnd')

  ! Initialize MPI communication groups for flux coupler

  call shr_msg_init ('lnd')
  call shr_msg_groups ('lnd')

  ! Initialize intra-MPI communication stuff or 
  ! set masterproc if not running in SPMD mode

#if (defined SPMD)
  call spmd_init()
#endif

  ! Initialize land model - initialize communication with flux coupler

  call initialize (eccen, obliqr, lambm0, mvelpp)

! -----------------------------------------------------------------
! Time stepping loop
! -----------------------------------------------------------------

  call t_startf('lnd_timeloop')

  ! begin time stepping loop

  do 

     ! doalb is true when the next time step is a radiation time step
     ! this allows for the fact that an atmospheric model may not do 
     ! the radiative calculations every time step. for example:
     !      nstep dorad doalb
     !        1     F     F
     !        2     F     T
     !        3     T     F
     ! The following expression for doalb is specific to CAM
     
     nstep = get_nstep() 
     doalb = ((irad==1 .and. nstep+1/=1) .or. (mod(nstep,irad)==0 .and. nstep+1/=1))
     
     ! Call land surface model driver 
     
     call driver (doalb, eccen, obliqr, lambm0 ,mvelpp)
     
     ! determine if time to stop
     
     if (csmstop_now) exit 
     
     ! increment time step 
     
     call advance_timestep()

  end do

  call t_stopf ('lnd_timeloop')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
  call shr_msg_finalize

  stop
end program program_csm

#else

!The following is only here since empty file won't compile
subroutine program_csm_stub
  write(6,*) 'PROGRAM_CSM: this routine should not be called'
  stop 99
end subroutine program_csm_stub

#endif





