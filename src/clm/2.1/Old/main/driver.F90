#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: driver
!
! !INTERFACE: 
subroutine driver (doalb, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION: 
! Main CLM2 driver calling sequence:
! -> loop over subgrid hierarchy, calling for each column:
!    -> Hydrology1          canopy interception and precip on ground
!    -> Biogeophysics1      leaf temperature and surface fluxes
!    -> Biogeophysics_Lake  lake temperature and surface fluxes
!    -> Biogeophysics2      soil/snow and ground temp and update surface fluxes
!    -> Hydrology2          surface and soil hydrology
!    -> Hydrology_Lake      lake hydrology
!    -> Biogeochemistry     surface biogeochemical fluxes (LSM)
!    -> EcosystemDyn:       ecosystem dynamics: phenology, vegetation, soil carbon
!    -> SurfaceAlbedo:      albedos for next time step
!    -> BalanceCheck        check for errors in energy and water balances
!  -> write_diagnostic      output diagnostic if appropriate
!  -> Rtmriverflux          calls RTM river routing model
!  -> update_hbuf           accumulate history fields over history time interval
!  -> htapes_wrapup         write history tapes if appropriate
!  -> restart               write restart file if appropriate
!  -> inicwrt               write initial file if appropriate 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clmpoint
  use globals
  use spmdMod, only : masterproc 
  use clm_varcon, only : zlnd		
  use time_manager, only : get_step_size, get_curr_calday, get_curr_date, get_nstep
  use histFileMod, only : update_hbuf, htapes_wrapup 
  use restFileMod, only : restart
  use inicFileMod, only : inicwrt, do_inicwrite 
  use BalanceCheckMod, only : BalanceCheck
  use Hydrology1Mod, only : Hydrology1
  use Hydrology2Mod, only : Hydrology2
  use HydrologyLakeMod, only : HydrologyLake
  use Biogeophysics1Mod, only : Biogeophysics1
  use Biogeophysics2Mod, only : Biogeophysics2
  use BiogeophysicsLakeMod, only : BiogeophysicsLake
  use SurfaceAlbedoMod, only : SurfaceAlbedo, Snowage
  use pft2columnMod, only : pft_to_col_glob, pft_to_col_soil
#if (defined COUP_CSM)
  use clm_varctl, only : wrtdia, fsurdat, csm_doflxave
#endif
!CAS  3/2005
!CAS  Added to incorporate PFT distr. of PCP into zonal-CLM framework
!CAS  3/2005
#if (defined PCP2PFT)
  use clm_varctl, only : wrtdia, fsurdat, fpcp2pft
#else
  use clm_varctl, only : wrtdia, fsurdat
#endif
  use EcosystemDynMod, only : EcosystemDyn
  use mvegFileMod, only : interpmonthlyveg
!CAS  3/2005
!CAS  Added to incorporate PFT distr. of PCP into zonal-CLM framework
!CAS  3/2005
#if (defined PCP2PFT)
  use pcp2pftMod
#endif
!CAS  3/2005
!CAS  Added to incorporate PFT distr. of PCP into zonal-CLM framework
!CAS  3/2005
#if (defined RTM)
  use RtmMod, only : Rtmriverflux
#endif
#if (defined COUP_CSM)
  use clm_csmMod, only : csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
       dorecv, dosend, csmstop_now
#endif
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: doalb  !true if time for surface albedo calculation
  real(r8), intent(in) :: eccen  !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0 !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp !Earth's moving vernal equinox long. of perihelion + pi (radians)
!
! !CALLED FROM:
! program program_off (if COUP_OFFLINE cpp variable is defined)
! program program_csm (if COUP_CSM cpp variable is defined)
! subroutine atm_lnddrv in module atm_lndMod (if COUP_CAM cpp variable is defined)
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: ci,pi,j               !indices
  integer  :: yrp1                  !year (0, ...) for nstep+1
  integer  :: monp1                 !month (1, ..., 12) for nstep+1
  integer  :: dayp1                 !day of month (1, ..., 31) for nstep+1
  integer  :: secp1                 !seconds into current date for nstep+1   
#if (defined PCP2PFT)
  integer  :: yrp2p                  !year (0, ...) for nstep+1
  integer  :: monp2p                 !month (1, ..., 12) for nstep+1
  integer  :: dayp2p                 !day of month (1, ..., 31) for nstep+1
  integer  :: secp2p                 !seconds into current date for nstep+1   
#endif
  real(r8) :: caldayp1              !calendar day for nstep+1
  logical, external :: do_restwrite !determine if time to write restart   
  type(column_type), pointer :: c   !local pointer to derived subtype
  integer  :: n                     !temporary    
!-----------------------------------------------------------------------

  call t_startf('clm_driver')

  ! ----------------------------------------------------------------------
  ! Calendar information for current time step
  ! ----------------------------------------------------------------------

  nstep = get_nstep() 
  call get_curr_date (ctl%year, ctl%month, ctl%day, ctl%secs)

  ! ----------------------------------------------------------------------
  ! Coupler receive
  ! ----------------------------------------------------------------------

#if (defined COUP_CSM)
  ! Determine if information should be sent/received to/from flux coupler
  call csm_dosndrcv(doalb)

  ! Get atmospheric state and fluxes from flux coupler
  if (dorecv) then
     call csm_recv()
     if (csmstop_now) then
        call t_stopf('clm_driver')
        RETURN
     endif
  endif
#endif

  ! ----------------------------------------------------------------------
  ! Calendar information for next time step
  ! o caldayp1 = calendar day (1.00 -> 365.99) for cosine solar zenith angle
  !              calday is based on Greenwich time
  ! o monp1    = month (1 -> 12) for leaf area index and stem area index
  ! o dayp1    = day (1 -> 31)   for leaf area index and stem area index
  ! ----------------------------------------------------------------------

  dtime = get_step_size()
  caldayp1 = get_curr_calday(offset=int(dtime))
  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))

#if (defined PCP2PFT)

  ! ----------------------------------------------------------------------
  ! Determine weights for monthly PFT distribution of Precpitation.
  ! This also determines whether it is time to read new monthly data
  ! (variable name 'pcp2pft').
  ! ----------------------------------------------------------------------

  call get_curr_date(yrp2p, monp2p, dayp2p, secp2p)
  if (monp2p /= monp1) call readMonthlyPcp2Pft(fpcp2pft,monp1)

#endif

  ! ----------------------------------------------------------------------
  ! Determine weights for time interpolation of monthly vegetation data.
  ! This also determines whether it is time to read new monthly vegetation and
  ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
  ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
  ! weights obtained here are used in subroutine ecosystemdyn to obtain time
  ! interpolated values.
  ! ----------------------------------------------------------------------

  if (doalb) call interpMonthlyVeg (fsurdat, monp1, dayp1)

  ! ----------------------------------------------------------------------
  ! LOOP 1
  ! ----------------------------------------------------------------------

  ! Initialize driver loops

  call t_startf('clm_loop1')
!$OMP PARALLEL DO PRIVATE (ci,c,n)
  do ci = cols1d%beg,cols1d%end
     c => cpoint%col(ci)%c

     ! Save snow mass at previous time step
     c%cws%h2osno_old = c%cws%h2osno

     ! Initialize fraction of vegetatino not covered by snow
     c%cps%pps_a%frac_veg_nosno = c%cps%pps_a%frac_veg_nosno_alb
     do pi = 1, c%cps%npfts
        c%p(pi)%pps%frac_veg_nosno = c%p(pi)%pps%frac_veg_nosno_alb
     end do

     ! Decide whether to cap snow
     if (c%cws%h2osno > 1000.) then
        c%cps%do_capsnow = .true.
     else
        c%cps%do_capsnow = .false.
     end if

     ! Non-lake points
     if (.not. c%cps%lps%lakpoi) then

        ! Initial set of previous time-step variables
        ! ice fraction of snow at previous time step
        n = c%cps%snl+1
        do j=n,0
           c%cps%frac_iceold(j) = c%cws%h2osoi_ice(j)/(c%cws%h2osoi_liq(j)+c%cws%h2osoi_ice(j))
        end do

        ! Determine beginning water balance (at previous time step)
        c%cwbal%begwb = c%cws%pws_a%h2ocan + c%cws%h2osno
        do j=1, nlevsoi
           c%cwbal%begwb = c%cwbal%begwb + c%cws%h2osoi_ice(j) + c%cws%h2osoi_liq(j)
        end do

     endif
  end do ! end column loop

!$OMP PARALLEL DO PRIVATE (ci,c)
  do ci = cols1d%beg,cols1d%end
     c => cpoint%col(ci)%c

     if (.not. c%cps%lps%lakpoi) then

        ! Determine canopy interception and precipitation onto ground surface.
        ! Determine the fraction of foliage covered by water and the fraction
        ! of foliage that is dry and transpiring. Initialize snow layer if the
        ! snow accumulation exceeds 10 mm.
        ! Hydrology1() includes a loop through all the pfts for this column
        call Hydrology1(c)

        ! Determine leaf temperature and surface fluxes based on ground
        ! temperature from previous time step.
        call Biogeophysics1(c)

     else if (c%cps%lps%lakpoi) then

        ! Determine lake temperature and surface fluxes
        call BiogeophysicsLake(c)

     endif

     if (.not. c%cps%lps%lakpoi) then

        ! Ecosystem dynamics: 
        ! phenology, vegetation, soil carbon, snow fraction
        call EcosystemDyn(c, doalb, .false.)

     end if

     if (doalb) then

        ! Determine albedos for next time step
        call SurfaceAlbedo(c, caldayp1, eccen, obliqr, lambm0, mvelpp)

     endif

     if (.not. c%cps%lps%lakpoi) then

        ! Determine soil/snow temperatures including ground temperature and
        ! update surface fluxes for new ground temperature.
        call Biogeophysics2(c)

        ! Perform averaging from PFT level to column level - non lake points
        call pft_to_col_soil(c)

     endif

     ! Perform averaging from PFT level to column level - lake and non-lake
     call pft_to_col_glob(c)

  end do ! end column loop
  call t_stopf('clm_loop1')

  ! ----------------------------------------------------------------------
  ! Coupler send
  ! ----------------------------------------------------------------------

#if (defined COUP_CSM)
  ! Average fluxes over interval if appropriate
  ! Surface states sent to the flux coupler states are not time averaged
  if (csm_doflxave) call csm_flxave()

  ! Send fields to flux coupler
  ! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
  ! albedos[n+1], and ocnrof_vec[n]
  if (dosend) call csm_send()
#endif

  ! ----------------------------------------------------------------------
  ! LOOP 2
  ! ----------------------------------------------------------------------

  call t_startf('clm_loop2')
!$OMP PARALLEL DO PRIVATE (ci,c)
  do ci = cols1d%beg,cols1d%end
     c => cpoint%col(ci)%c

     ! Vertical (column) soil and surface hydrology
     if (.not. c%cps%lps%lakpoi) call Hydrology2 (c)

     ! Lake hydrology
     if (c%cps%lps%lakpoi) call HydrologyLake (c)

     ! Update Snow Age (needed for surface albedo calculation - but is
     ! really a column type property
     call SnowAge (c)

     ! Fraction of soil covered by snow - really a column property
     c%cps%frac_sno = c%cps%snowdp/(10.*zlnd + c%cps%snowdp)

     ! Check the energy and water balance
     call BalanceCheck (c)

  end do   ! end column loop
  call t_stopf('clm_loop2')

  ! ----------------------------------------------------------------------
  ! Write global average diagnostics to standard output
  ! ----------------------------------------------------------------------

  call write_diagnostic (wrtdia, nstep)

#if (defined RTM)
  ! ----------------------------------------------------------------------
  ! Route surface and subsurface runoff into rivers
  ! ----------------------------------------------------------------------

  call t_startf('clm_rtm')
  call Rtmriverflux ()
  call t_stopf('clm_rtm')
#endif

  ! Update history buffer

  call t_startf('updat_hbuf')
  call update_hbuf()
  call t_stopf('updat_hbuf')

  call t_startf('clm_output')

  ! Create history tapes if  appropriate and  write to history tapes if appropriate

  call htapes_wrapup()

  ! Write restart files if appropriate

  if (do_restwrite()) call restart('write')

  ! Write intial files if appropriate

  if (do_inicwrite()) call inicwrt()

  call t_stopf('clm_output')

  call t_stopf('clm_driver')

  return
end subroutine driver

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
subroutine write_diagnostic (wrtdia, nstep)
!
! !DESCRIPTION: 
! Write diagnostic output.
!
! !USES:
  use clmtype
  use clmpoint
#if (defined SPMD)
  use spmdMod, only : masterproc, gather_data_to_master, mpicom
#else
  use spmdMod, only : masterproc
#endif
  use shr_sys_mod, only : shr_sys_flush
  use system_messages
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia                !true => write diagnostic   
  integer, intent(in) :: nstep                 !model time step
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: ci                               !index
  integer  :: ier                              !error status 
  integer  :: begc,endc,numc                   !column 1d indices
  real(r8) :: tsxyav                           !average ts for diagnostic output
  real(r8), pointer :: rloc(:)                 !temporary buffer 
  real(r8), pointer :: rglob(:)                !temporary buffer
!------------------------------------------------------------------------

  begc = cols1d%beg
  endc = cols1d%end
  numc = cols1d%num

  if (wrtdia) then

#if (defined TIMING_BARRIERS)
     call t_startf ('sync_write_diag')
     call mpibarrier (mpicom)
     call t_stopf ('sync_write_diag')
#endif
     allocate (rloc(begc:endc), stat=ier)
     call allocation_err(ier, 'write_diagnostic', 'rloc', endc-begc+1)
!$OMP PARALLEL DO PRIVATE (ci)
     do ci = begc,endc
        rloc(ci) = clmpointer(ic_ces_t_rad_column)%val(ci)%rp
     end do
#if (defined SPMD)
     if (masterproc) then
        allocate (rglob(numc), stat=ier)
        call allocation_err(ier, 'write_diagnostic', 'rglob', numc)
     endif
     call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
#else
     rglob => rloc
#endif
     if (masterproc) then
        tsxyav = 0._r8
        do ci = 1,numc
           tsxyav = tsxyav + rglob(ci)*cols1d%wtxy(ci)
        end do
        tsxyav = tsxyav / grid1d%num
!        write (6,1000) nstep, tsxyav
        call shr_sys_flush(6)
     end if
#if (defined SPMD)
     if (masterproc) deallocate(rglob)
#endif
     deallocate (rloc)

  else

#if (!defined COUP_CAM)
     if (masterproc) then
        write(6,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(6)
     endif
#endif

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',e21.15)
  return
end subroutine write_diagnostic
