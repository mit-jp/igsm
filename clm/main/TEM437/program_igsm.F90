#include <misc.h>
#include <preproc.h>

#if (defined COUP_MIT2D)

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
! !IROUTINE: program_igsm
!
PROGRAM program_igsm
!
! !DESCRIPTION: 
! "off-line" code to mimic coupling to an atmospheric model.
! This program is an "off-line" driver for clm2.
! This code can be used to run the clm2 uncoupled from any atmospheric model. 
! The appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
! o If running as an offline driver, the land surface model may use
!   a different grid than the input atmospheric data. The atmospheric
!   data is then interpolated to the land model grid inside the
!   atmospheric driver module [atmdrvMod.F90].
! o If running as part of cam, the land surface model must use the
!   same grid as the cam.
! o If running through the flux coupler, the land surface model grid
!   is interpolated to the atmospheric grid inside the flux coupler
! o To map from the atmospheric grid to the land grid, the atmospheric
!   model must provide latitudes and longitudes (degrees) for each grid
!   point and the North, East, South, and West edges of atmospheric grid.
!   Comparable data for the land grid are provided by the land model. 
!   When mapping from land to atm grid, an atm grid cell that is part
!   land and part ocean (as defined by the land surface grid) will have
!   fluxes only based on the land portion.
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
! Author: Gordon Bonan and Mariana Vertenstein
! 11/30/01 Peter Thornton : Added use globals, removed use clm_varctl
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
    use time_manager

    implicit none
 
!
! atmospheric gridded data
!
    integer,parameter :: lsmlon=1
    integer,parameter :: lsmlat=46
    real(r8) :: psmit (lsmlon,lsmlat)
    real(r8) :: pcptmit (lsmlon,lsmlat)
    real(r8) :: pcplmit (lsmlon,lsmlat)
    real(r8) :: pcpcmit (lsmlon,lsmlat)
    real(r8) :: tprmit (lsmlon,lsmlat)
    real(r8) :: tslmit (lsmlon,lsmlat)
    real(r8) :: qsmit (lsmlon,lsmlat)
    real(r8) :: wsmit (lsmlon,lsmlat)
    real(r8) :: usmit (lsmlon,lsmlat)
    real(r8) :: vsmit (lsmlon,lsmlat)
    real(r8) :: dswmit (lsmlon,lsmlat)
    real(r8) :: dlwmit (lsmlon,lsmlat)
    real(r8) :: pco2mit (lsmlon,lsmlat)
!    real(r8) :: po3mit (lsmlon,lsmlat)
    real(r8) :: swnirmit (lsmlon,lsmlat)
    real(r8) :: swparmit (lsmlon,lsmlat)
!  common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,po3mit,swnirmit,swparmit
    common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit

!	TEM variables
!
    integer, parameter :: mxmsaics = 35
    integer, parameter :: mxmthdys = 31
    integer, parameter :: mxdayhrs = 24
    integer, parameter :: max3hrs = 8
    integer, parameter :: mxnlayers = 6
    integer, parameter :: totlayers = 10
    real(r8) :: mitco2(lsmlat)
    real(r8) :: mito3(max3hrs,lsmlat)
    real(r8) :: mittemp(lsmlat)
    real(r8) :: mitdaytemp(mxmthdys,lsmlat)
    real(r8) :: mitswrs(lsmlat)
    real(r8) :: mitpre(lsmlat)
    real(r8) :: mitstrmdur(mxmthdys,lsmlat)
    real(r8) :: mitqstrm(mxmthdys,lsmlat)
    real(r8) :: mitaet(mxmsaics,lsmlat)
    real(r8) :: mitsh2o1m(mxmsaics,lsmlat)
    real(r8) :: mitsh2o2m(mxmsaics,lsmlat)
    real(r8) :: mitswe(mxmsaics,lsmlat)
    real(r8) :: mitsfr(mxmsaics,lsmlat)
    real(r8) :: mitdrn(mxmsaics,lsmlat)
    real(r8) :: mitdaytsoil(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mitdaysh2o(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mithrsh2o(mxdayhrs,mxmthdys,mxmsaics,lsmlat,mxnlayers)
    common/climate4tem/mitco2,mito3,mittemp,mitdaytemp,mitswrs,mitpre,mitstrmdur,&
                       mitqstrm,mitaet,mitsh2o1m,mitsh2o2m,mitswe,mitsfr,mitdrn,&
                       mitdaytsoil,mitdaysh2o,&
                       mithrsh2o
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
    integer :: JYEAR,JMONTH,JDATE,TOFDAY
    integer, save :: nstep
!  integer, save :: iyear_AD, atmlon, atmlat
!  logical, save :: log_print
    logical, save :: addSAtip
    logical, save :: moredata
    data addSAtip/.false./
    data moredata/.true./
!  real(r8), save :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
!-----------------------------------------------------------------------
    open (20,file='data4clm',form='unformatted')
    open (21,file='data4tem',form='unformatted')
!    open (20,file='/d1/adam/CLM/inputdata/MIT2d/data4clmsp11',form='unformatted')
!    open (21,file='/d1/adam/CLM/inputdata/MIT2d/data4tem11',form='unformatted')
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
        read(20,end=1000) psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit
        nstep = get_nstep()
        call clm4mit2d
#if (defined COUP_URBAN)
        print *,'TSCREEN to Urban Airshed:'
        print *,'========================'
        print *,tscreen(1,:,:)
        print *,'DALYPCP to Urban Airshed:'
        print *,'========================'
        print *,dlypcp(1,:,:)
#endif
        if (is_last_step) moredata=.false.
      end do
      end program program_igsm

#else

!The following is only here since empty file won't compile
subroutine program_igsm_sub
  write(6,*) 'PROGRAM_OFF: this routine should not be called'
  return
end subroutine program_igsm_sub

#endif
