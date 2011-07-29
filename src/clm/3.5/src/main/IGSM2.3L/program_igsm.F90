#include <misc.h>
#include <preproc.h>

#if (defined COUP_MIT2D)

!================================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: program_igsm
!
PROGRAM program_igsm
!
! !DESCRIPTION: 
!   "off-line" code to mimic coupling to the IGSM atmospheric model.
!   This program is an "off-line" driver for clm3.
!   The appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
! o The zenith angle calculation is for the NEXT time step rather 
!   than the current time step. Make sure the calendar day is for 
!   the NEXT time step. Make sure the calendar day is for Greenwich 
!   time (see next comment).
! o The land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. The net longwave radiation
!   at the surface will differ somewhat from that calculated in the
!   atmospheric model because the atm model will use the upward 
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. The net solar radiation should equal that calculated
!   in the atmospheric model. If not, there is a problem in how
!   the models are coupled.
!
! !REVISION HISTORY:
! Author: Adam Schlosser
!
!EOP
!
!   !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
!    use controlMod, only : control_init, control_print
!    use clm_varpar
!    use clm_varctl
!    use clm_varsur
!    use initializeMod
!    use atmdrvMod
!    use shr_orb_mod
    use clm_time_manager

    implicit none
    include 'IGSM2.inc'
#if (defined COUP_TEM)
    include 'TEM.inc'
#endif
#if (defined COUP_URBAN)
    real :: tscreen(lsmlon,lsmlat,mxmsaics)      ! Screen temperature (for each PFT) to Urban airshed model
    real :: dlypcp(lsmlon,lsmlat,mxmsaics)       ! Daily precipitation (for each PFT) to Urbran airshed model
    common/chem_meta_gls/tscreen,dlypcp
#endif
!
! !LOCAL VARIABLES:
  integer ier       !error code

! Earth's orbital characteristics
! Orbital information after call to routine shr_orbit_params

    integer :: igsmloop
    integer :: JMONTH,JDATE,TOFDAY,njyear
    integer :: jyear
    integer, save :: nstep
!  integer, save :: iyear_AD, atmlon, atmlat
!  logical, save :: log_print
    logical, save :: addSAtip
    logical, save :: moredata
    data addSAtip/.false./
    data moredata/.true./
!  real(r8), save :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
!-----------------------------------------------------------------------
    call get_start_date(jyear,jmonth,jday,jtod)
    write (flyr,'i4.4') jyear
    fl4clm = trim(igsmdir) // 'data4clm' // flyr
    !fl4tem = trim(igsmdir) // 'data4tem' // flyr
    print *,'Opening data files: '
    print *, fl4clm
    print *, 'Opening CLM forcing:',fl4clm
    call opnfil(fl4clm,20,'u')
    !print *, fl4tem
    !print *, 'Opening TEM forcing:',fl4clm
    !call opnfil(fl4tem,21,'u')
!
! initialization call (done before igsmloop)
!
!     call t_setoptionf (usrsys, 0)
!     call t_initializef ()
!     log_print = .true.
!     iyear_AD= 1950
!     eccen = SHR_ORB_UNDEF_REAL
!     obliq = SHR_ORB_UNDEF_REAL
!     mvelp = SHR_ORB_UNDEF_REAL
!     call shr_orb_params (iyear_AD, eccen, obliq, mvelp, obliqr,lambm0, mvelpp, log_print)
!     call initialize (eccen, obliqr, lambm0, mvelpp)
!     call t_startf('total')
!     call clm4mit2d
1000  rewind 20
      do while (moredata)
#if (defined PCPTRND)
        read(20,end=1000) psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit,tglb
#else
        read(20,end=1000) psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit
#endif
        !nstep = get_nstep()
        call clm4mit2d
        call get_curr_date(njyear,jmonth,jday,jtod)
        if (is_last_step) moredata=.false.
        if ( njyear > jyear ) then
          jyear = njyear
#if (defined PCPTRND)
          if ( jyear < nyt ) then
            tglbtrnd=ttrnd(jyear-year1+1)
            print *,'Temperature change for ',jyear,' is : ',tglbtrnd
          else
            print *,'Ran out of temperature trend data in: ',fl4ttrnd
            print *,'Continuing with no further trend to affect pcp frequency'
            tglbtrnd=ttrnd(nyt)
          endif
#endif
          write (flyr,'i4.4') jyear
          fl4clm = trim(igsmdir) // 'data4clm' // flyr
          fl4tem = trim(igsmdir) // 'data4tem' // flyr
          call getfil(fl4clm,locfn, 1)
          inquire (file = locfn, exist = lexist)
          if (lexist) then
            close (20)
            close (21)
            print *,'Opening data files: '
            print *, fl4clm
            print *, fl4tem
            print *, 'for CLM and TEM forcing'
            call opnfil(fl4clm,20,'u')
            call opnfil(fl4tem,21,'u')
          else
            print *,'Reached end of atmos. forcing timeseries'
            print *,'========================================'
            print *,'Will repeat last year'
            print *,'====================='
            rewind (20)
            rewind (21)
          endif
        endif
      end do
      end program program_igsm

#else

!The following is only here since empty file won't compile
subroutine program_igsm_sub
  write(6,*) 'PROGRAM_IGSM: this routine should not be called'
  return
end subroutine program_igsm_sub

#endif
