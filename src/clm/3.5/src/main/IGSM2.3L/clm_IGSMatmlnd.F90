!=======================================================================

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: igsmdrv
!
! !INTERFACE:
  subroutine igsmdrv(nstep)
!
! !DESCRIPTION:  CAS 05/2009
! This code reads in atmospheric fields from an input file or from a
! memory block for the IGSM and generates the required atmospheric forcing.
! The atmospheric forcing to the Global Land System (GLS) of the IGSM can either
! be in the zonal, mosaic structure (IGSM2.2) or the zonal atmosphere variables
! can be projected onto a 3D configure for the land grid (and projection coefficients
! are based on an observed climatology with trends (if selected) based on an
! AR4 GCM.
!
! ============================
! Possible atmospheric fields:
! ============================
! Name     Description                              Required/Optional
! -----------------------------------------------------------------------------
! TBOT     temperature (K)                          Required
! WIND     wind:sqrt(u**2+v**2) (m/s)               Required
! QBOT     specific humidity (kg/kg)                Required
! Tdew     dewpoint temperature (K)                 Alternative to Q
! RH       relative humidity (percent)              Alternative to Q
! ZBOT     reference height (m)                     optional
! PSRF     surface pressure (Pa)                    optional
! FSDS     total incident solar radiation (W/m**2)  Required
! FSDSdir  direct incident solar radiation (W/m**2) optional (replaces FSDS)
! FSDSdif  diffuse incident solar rad (W/m**2)      optional (replaces FSDS)
! FLDS     incident longwave radiation (W/m**2)     optional
! PRECTmms total precipitation (mm H2O / sec)       Required
! PRECCmms convective precipitation (mm H2O / sec)  optional (replaces PRECT)
! PRECLmms large-scale precipitation (mm H2O / sec) optional (replaces PRECT)
!
! ===============
! Namelist input:
! ===============
! character*256 offline_atmdir = directory for input atm data files (can be Mass Store)
!
! !USES:
    use nanMod
    use decompMod   , only : adecomp, get_proc_bounds_atm
    use clm_atmlnd  , only : clm_mapa2l, atm_a2l, clm_a2l
    use clm_varctl  , only : offline_atmdir, pertlim
#if (defined STOCHASTIC)
    use clm_varcon  , only : rair, cpair, co2_ppmv_const, o2_molar_const, tcrit, c13ratio, spval
#else
    use clm_varcon  , only : rair, cpair, co2_ppmv_const, o2_molar_const, tcrit, c13ratio
#endif
    use clm_time_manager, only : get_step_size, get_curr_calday, get_curr_date
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'IGSM2.inc'
    integer, intent(in) :: nstep    !current time step
!
! !REVISION HISTORY:
! Created by Adam Schlosser
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,n,k,g,g1           !indices
    integer :: j1                     !indices
    integer :: itimlast               !last time index used in atmrd
    real(r8):: calday                 !calendar day at Greenwich (1.00 -> 365.99)
    integer :: kda                    !day (1 -> 31)
    integer :: kmo                    !month (1 -> 12)
    integer :: kyr                    !year (0 -> ...)
    integer :: ksec                   !current seconds of current date (0 -> 86400)
    integer :: mcdate                 !current date in integer format [yyyymmdd]
    integer :: mcsec                  !current time of day [seconds]
    integer :: dtime                  !time step size
    integer :: minpday = 1440         !minutes per day
    integer :: secpmin = 60           !seconds per minute
    integer, SAVE :: itim             !time index used in atmrd
    integer, SAVE :: atmmin           !temporal resolution of atm data (in minutes)
    character(len=256), SAVE :: locfn !full file name in case atmdir is in MSS
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
!------------------------------------------------------------------------
    j1 = 0
    !write (6,*) 'ATMDRV: datlon = ',datlon
    !write (6,*) 'ATMDRV: datlat = ',datlat
#if (!defined 2dato3dl)
    cpxy(:,:) = 1.
    ctxy(:,:) = 1.
#endif
    do j = 1, datlat
        do i = 1, datlon
            n = (j1-1)*datlon + i
            aV_atm_d2a%rAttr(if_txy,n) = ctxy(i,j)*tslmit(1,j)
            !write (6,*) 'ATMDRV: i = ',i,' j = ',j,' n = ',n
            !write (6,*) 'ATMDRV:  tslmit = ',tslmit(i,j)
            !write (6,*) 'ATMDRV:  aV_atm_txy = ',aV_atm_d2a%rAttr(if_txy,n)
            aV_atm_d2a%rAttr(if_uxy,n) = usmit(1,j)
            aV_atm_d2a%rAttr(if_vxy,n) = vsmit(1,j)
            aV_atm_d2a%rAttr(if_psrfxy,n) = psmit(1,j)*100
            aV_atm_d2a%rAttr(if_pbotxy,n) = psmit(1,j)*100
            aV_atm_d2a%rAttr(if_qxy,n) = qsmit(1,j)
            aV_atm_d2a%rAttr(izgcmxy,n) = 100.0
            aV_atm_d2a%rAttr(if_sols,n) = 0.7*swparmit(1,j)
            aV_atm_d2a%rAttr(if_soll,n) = 0.7*swnirmit(1,j)
            aV_atm_d2a%rAttr(if_solld,n) = 0.3*swnirmit(1,j)
            aV_atm_d2a%rAttr(if_solsd,n) = 0.3*swparmit(1,j)
            aV_atm_d2a%rAttr(iprcxy,n) = cpxy(i,j)*pcpcmit(1,j)/3600.
            aV_atm_d2a%rAttr(iprlxy,n) = cpxy(i,j)*pcplmit(1,j)/3600.

#if (defined STOCHASTIC)
!CAS
!CAS    Stochastic treatment of precipitation added for MIT2D model implementation/coupling
!CAS
!CAS    C. Adam Schlosser, June, 2004
!CAS   
!CAS   
!      
!       Check for beginning of month, if so, reset storm count
!      
       
            if (kda == 1 .and. mcsec == 0) storms = 0
            call stoch(aV_atm_d2a%rAttr(iprcxy,n),aV_atm_d2a%rAttr(iprlxy,n), &
                     prc_poiss(i,j),prl_poiss(i,j), &
                     storm_dur(i,j),exp_tstorm(i,j,kmo),exp_tdry(i,j,kmo), &
                     t_storm(i,j),t_dry(i,j),dtcumu(i,j),storms(i,j),pcpc_resid(i,j),pcpl_resid(i,j))
#endif
            aV_atm_d2a%rAttr(iflwdsxy,n) = dlwmit(1,j)
        enddo
     enddo
