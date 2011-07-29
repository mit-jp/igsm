#include <misc.h>
#include <preproc.h>

module stochrdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stochrdMod
!
! !DESCRIPTION:
!	Reading in stochastic precipitation information for zonal CLM/IGSM2 configuration
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : lsmlon, lsmlat
  use clm_varpar  , only : nlevsoi, numpft, &
                           maxpatch_pft, maxpatch_cft, maxpatch, &
                           npatch_urban, npatch_lake, npatch_wet, npatch_glacier
  use clm_varsur  , only : wtxy, vegxy
  use clm_varsur  , only : pctspec
  use ncdio
  use clmtype
  use spmdMod                         
  use clm_varctl,   only : scmlat, scmlon, single_column
  use domainMod   , only : domain_type
  use decompMod   , only : get_proc_bounds,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo

!
! !PUBLIC TYPES:
  implicit none
  save

  integer, parameter :: mthsperyr=12
  real,allocatable :: exp_tstorm(:,:,:)
  real,allocatable :: exp_tdry(:,:,:)
  real,allocatable :: trnd_tdry(:,:,:)
  real,allocatable :: t_dry(:,:)
  real,allocatable :: t_storm(:,:)
  real,allocatable :: pcpc_resid(:,:)
  real,allocatable :: pcpl_resid(:,:)
  real,allocatable :: prc_poiss(:,:)
  real,allocatable :: prl_poiss(:,:)
  integer,allocatable :: storms(:,:)
  integer,allocatable :: storm_dur(:,:)
  integer,allocatable :: dtcumu(:,:)
  integer :: ist,jst

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd_get_stoch  ! Read surface dataset and determine subgrid weights

!
! !REVISION HISTORY:
! Created by A. Schlosser
!
!EOP
!
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_stoch
!
! !INTERFACE:
  subroutine surfrd_get_stoch(domain,filename)
!
! !DESCRIPTION:
! Read the topo dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by A. Schlosser
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: beg3d(3),len3d(3)   ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    real(r8):: ival = 0.0_r8       ! initial value
    real(r8):: ival1 = 3600_r8    ! initial value
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    character(len=32) :: subname = 'surfrd_get_stoch'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       allocate (exp_tdry(lsmlon,lsmlat,mthsperyr))
       allocate (exp_tstorm(lsmlon,lsmlat,mthsperyr))
       allocate (trnd_tdry(lsmlon,lsmlat,mthsperyr))
       allocate (t_dry(lsmlon,lsmlat))
       allocate (t_storm(lsmlon,lsmlat))
       allocate (pcpc_resid(lsmlon,lsmlat))
       allocate (pcpl_resid(lsmlon,lsmlat))
       allocate (prc_poiss(lsmlon,lsmlat))
       allocate (prl_poiss(lsmlon,lsmlat))
       allocate (storms(lsmlon,lsmlat))
       allocate (storm_dur(lsmlon,lsmlat))
       allocate (dtcumu(lsmlon,lsmlat))
       t_dry = ival1
       t_storm = ival1
       pcpc_resid = ival
       pcpl_resid = ival
       prc_poiss = ival
       prl_poiss = ival
       storms = int(ival)
       storm_dur = ival
       dtcumu = ival

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)

       if (lsmlon /= ni .or. lsmlat /= nj) then
          write(6,*) trim(subname),' ERROR: Stochastic file mismatch ni,nj',domain%ni,ni,domain%nj,nj
          call endrun()
       endif

       beg3d(1) = 1     ;  len3d(1) = lsmlon
       beg3d(2) = 1     ;  len3d(2) = lsmlat
       beg3d(3) = 1     ;  len3d(3) = mthsperyr

       call check_ret(nf_inq_varid(ncid, 'EXP_TSTORM', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, exp_tstorm), subname)

       call check_ret(nf_inq_varid(ncid, 'EXP_TDRY', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, exp_tdry), subname)

#if (defined PCPTRND)
       call check_ret(nf_inq_varid(ncid, 'TRND_TDRY', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, trnd_tdry, beg3d, len3d), subname)
#endif

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

!CAS: comment CLMMPI stuff out for now
    !call clmmpi_bcast (domain%topo , size(domain%topo) , CLMMPI_REAL8  , 0, clmmpicom, ier)

  end subroutine surfrd_get_stoch
end module stochrdMod
