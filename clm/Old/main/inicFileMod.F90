#include <misc.h>
#include <preproc.h>

module inicFileMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: inicFileMod
!
! !DESCRIPTION: 
!
! Read and writes CLM2 initial data netCDF files
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clmpoint
  use clm_varpar, only : nlevsno, nlevsoi, nlevlak, rtmlon, rtmlat
#if (defined RTM)
  use RtmMod
#endif
#if (defined SPMD)
  use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master, iam
#else
  use spmdMod, only : masterproc
#endif
  use shr_sys_mod, only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: do_inicwrite ! true=> write initial data file
  public  :: inicrd       ! read in initial data
  public  :: inicwrt      ! write out initial data  
!                              
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!  
! !PRIVATE METHODS:
  private :: set_init_filename
  interface write_ncd
     module procedure write_real_2d
     module procedure write_real_1d
     module procedure write_int_1d
  end interface
  interface read_ncd
     module procedure read_real_2d
     module procedure read_real_1d
     module procedure read_int_1d
  end interface
!----------------------------------------------------------------------- 

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine inicrd
!
! !INTERFACE:
  subroutine inicrd ()
!
! !DESCRIPTION: 
! Read initial data from netCDF instantaneous initial data file 
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : finidat
    use fileutils, only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 
!
!EOP
!
! !LOCAL VARIABLES
    character(len=256) :: locfn         !local file name
    integer :: ncid                     !netCDF dataset id
    integer :: dimid                    !netCDF dimension id 
    integer :: varid                    !netCDF variable id
    integer :: dimlen                   !input dimension length
    integer :: gi,li,ci,pi,j,index      !indices 
    integer :: ret                      !netcdf return code
    real(r8):: pftsum                   !temporary for pft averaging for columns
    type(gridcell_type), pointer :: g   !local pointer to derived subtype
    type(landunit_type), pointer :: l   !local pointer to derived subtype
    type(column_type)  , pointer :: c   !local pointer to derived subtype
    type(pft_type)     , pointer :: p   !local pointer to derived subtype
    integer , pointer :: ibufslc(:)     !pointer to memory to be allocated
    real(r8), pointer :: rbufslc(:)     !pointer to memory to be allocated
    real(r8), pointer :: rbufmlc(:,:)   !pointer to memory to be allocated
    integer , pointer :: ibufslp(:)     !pointer to memory to be allocated
    real(r8), pointer :: rbufslp(:)     !pointer to memory to be allocated
    real(r8), pointer :: rbufmlp(:,:)   !pointer to memory to be allocated
#if (defined SPMD)
    real(r8), pointer :: rloc(:)        !temporary 
    real(r8), pointer :: rglob(:)       !temporary 
#endif                                  
    integer begc,endc,numc              !column 1d indices
    integer begp,endp,nump              !pft 1d indices
    character(len=16) namec             !column 1d name 
    character(len=16) namep             !pft 1d name 
!------------------------------------------------------------------------

    ! Open netCDF data file and read data

    if (masterproc) then

       call getfil (finidat, locfn, 0)
       call wrap_open (locfn, nf_nowrite, ncid)

       ! Check consistency input dimensions

       call wrap_inq_dimid (ncid, 'gridcell', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= grid1d%num) then
          write (6,*) 'INICRD error: numgrid grid cell values disagree'
          write (6,*) 'finidat numgrid = ',dimlen,' model numgrid = ',grid1d%num
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'landunit', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= land1d%num) then
          write (6,*) 'INICRD error: '
          write (6,*) 'finidat numland = ',dimlen,' model numland = ',land1d%num
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'column', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= cols1d%num) then
          write (6,*) 'INICRD error: '
          write (6,*) 'finidat numcols = ',dimlen,' model numcols = ',cols1d%num
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'pft', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= pfts1d%num) then
          write (6,*) 'INICRD error: '
          write (6,*) 'finidat numpfts = ',dimlen,' model numpfts = ',pfts1d%num
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'levsno', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= nlevsno) then
          write (6,*) 'INICRD error: '
          write (6,*) 'finidat levsno = ',dimlen,' model nlevsno = ',nlevsno
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'levsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, dimlen)
       if (dimlen /= nlevsoi) then
          write (6,*) 'INICRD error: '
          write (6,*) 'finidat nlevsoi = ',dimlen,' model nlevsoi = ',nlevsoi
          call endrun
       end if

       ret = nf_inq_dimid (ncid, 'levlak', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, dimlen)
          if (dimlen /= nlevlak) then
             write (6,*) 'INICRD error: '
             write (6,*) 'finidat nlevlak = ',dimlen,' model nlevlak = ',nlevlak
             call endrun
          end if
       endif

#if (defined RTM)
       ret = nf_inq_dimid (ncid, 'rtmlon', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, dimlen)
          if (dimlen /= rtmlon) then
             write (6,*) 'INICRD error: '
             write (6,*) 'finidat rtmlon = ',dimlen,' model rtmlon = ',rtmlon
             call endrun
          end if
       endif
       ret = nf_inq_dimid (ncid, 'rtmlat', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, dimlen)
          if (dimlen /= rtmlat) then
             write (6,*) 'INICRD error: '
             write (6,*) 'finidat rtmlat = ',dimlen,' model rtmlat = ',rtmlat
             call endrun
          end if
       endif
#endif

    endif ! if-masterproc block

    ! Obtain data - for the snow interfaces, are only examing the snow 
    ! interfaces above zi=0 which is why zisno and zsno have the same 
    ! level dimension below

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num
    namec = cols1d%name

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num
    namep = pfts1d%name

    allocate (rbufslc(begc:endc))
    allocate (ibufslc(begc:endc))
    allocate (rbufslp(begp:endp))
    allocate (ibufslp(begp:endp))
    
    ! Read in zisno - column
    ! NOTE: zi(0) is set to 0 in routine iniTimeConst

    allocate(rbufmlc(-nlevsno:-1,begc:endc))
    call read_ncd (ncid, 'ZISNO', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%zi(-nlevsno:-1) = rbufmlc(-nlevsno:-1,ci)
    end do
    deallocate (rbufmlc)

    ! Read in zsno - column
    allocate (rbufmlc(-nlevsno+1:0,begc:endc))
    call read_ncd (ncid, 'ZSNO', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%z(-nlevsno+1:0) = rbufmlc(-nlevsno+1:0,ci)
    end do
    deallocate (rbufmlc)
    
    ! Read in dzsno - column
    allocate (rbufmlc(-nlevsno+1:0,begc:endc))
    call read_ncd (ncid, 'DZSNO', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%dz(-nlevsno+1:0) = rbufmlc(-nlevsno+1:0,ci)
    end do
    deallocate (rbufmlc)

    ! Read in h2osoi_liq - column
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    call read_ncd (ncid, 'H2OSOI_LIQ', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cws%h2osoi_liq(-nlevsno+1:nlevsoi) = rbufmlc(-nlevsno+1:nlevsoi,ci)
    end do
    deallocate (rbufmlc)
    
    ! Read in h2osoi_ice - column
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    call read_ncd (ncid, 'H2OSOI_ICE', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cws%h2osoi_ice(-nlevsno+1:nlevsoi) = rbufmlc(-nlevsno+1:nlevsoi,ci)
    end do
    deallocate (rbufmlc)
    
    ! Read in t_soisno - column
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    call read_ncd (ncid, 'T_SOISNO', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%ces%t_soisno(-nlevsno+1:nlevsoi) = rbufmlc(-nlevsno+1:nlevsoi,ci)
    end do
    deallocate (rbufmlc)
    
    ! Read in t_lake - column
    allocate (rbufmlc(1:nlevlak,begc:endc))
    call read_ncd (ncid, 'T_LAKE', rbufmlc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%ces%t_lake(1:nlevlak) = rbufmlc(1:nlevlak,ci)
    end do
    deallocate (rbufmlc)
    
    ! Read in t_veg - pft
    call read_ncd (ncid, 'T_VEG', rbufslp, clmlevel=namep)
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       p%pes%t_veg = rbufslp(pi)
    end do
    
    ! Read in t_grnd - column
    call read_ncd (ncid, 'T_GRND', rbufslc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%ces%t_grnd = rbufslc(ci)
    end do
    
    ! Read in h2ocan - pft 
    call read_ncd (ncid, 'H2OCAN', rbufslp, clmlevel=namep)
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       p%pws%h2ocan = rbufslp(pi)
    end do

    ! Obtain column average of h2ocan
    ! driver.F90 needs this upon input for an initial run
    do ci = begc,endc
       c => cpoint%col(ci)%c
       pftsum = 0.0
       do pi = 1,c%cps%npfts
          p => c%p(pi)
          pftsum = pftsum + p%pws%h2ocan * c%pw(pi)   
       end do
       c%cws%pws_a%h2ocan = pftsum
    end do

    ! Read in h2osno - column 
    call read_ncd (ncid, 'H2OSNO', rbufslc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cws%h2osno = rbufslc(ci)
    end do
    
    ! Read in snowdp - column
    call read_ncd (ncid, 'SNOWDP', rbufslc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%snowdp = rbufslc(ci)
    end do
    
    ! Read in snowage - column
    call read_ncd (ncid, 'SNOWAGE', rbufslc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%snowage = rbufslc(ci)
    end do
    
    ! Read in snlsno - column
    call read_ncd (ncid, 'SNLSNO', ibufslc, clmlevel=namec)
    do ci = begc,endc
       c => cpoint%col(ci)%c
       c%cps%snl = ibufslc(ci)
    end do

#if (defined RTM)
    ! Read in RTM volr
    if (masterproc) then
       ret = nf_inq_varid (ncid, 'RTMVOLR', varid)
       if (ret == NF_NOERR) then
          write(6,*)'INICFILE: reading in river volr'
          call wrap_get_var_realx(ncid, varid, volr)
       endif
    endif

    ! Determine initial fluxout from volr
    call Rtmfluxini()
#endif

    deallocate (ibufslc)
    deallocate (rbufslc)
    deallocate (ibufslp)
    deallocate (rbufslp)

    return
  end subroutine inicrd

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: inicwrt
!
! !INTERFACE:
  subroutine inicwrt ()
!
! !DESCRIPTION: 
! Write instantaneous initial data to netCDF initial data file
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varctl, only : caseid, ctitle, version, fsurdat, archive_dir, &
         mss_wpass, mss_irt
    use time_manager, only : get_nstep
    use fileutils, only : set_filename, putfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ncid                   !netCDF dataset id
    integer :: varid                  !netCDF variable id
    integer :: gi,li,ci,pi,j          !indices 
    integer :: dim1_id(1)             !netCDF dimension id for 1-d variables   
    integer :: dim2_id(2)             !netCDF dimension id for 2-d variables
    integer :: omode                  !netCDF dummy variable
    integer :: status                 !netCDF error status
    character(len=256) :: loc_fn      !local 
    character(len=256) :: rem_dir     !remote (archive) directory
    character(len=256) :: rem_fn      !remote (archive) filename
    character(len=256) :: str         !global attribute string 
    character(len=  8) :: curdate     !current date
    character(len=  8) :: curtime     !current time 
    type(gridcell_type), pointer :: g !local pointer to derived subtype
    type(landunit_type), pointer :: l !local pointer to derived subtype
    type(column_type)  , pointer :: c !local pointer to derived subtype
    type(pft_type)     , pointer :: p !local pointer to derived subtype
!
! Netcdf dimension Id's 
!
    integer :: gridcell_dimid         !netCDF dimension id
    integer :: landunit_dimid         !netCDF dimension id
    integer :: column_dimid           !netCDF dimension id
    integer :: pft_dimid              !netCDF dimension id
    integer :: levsoi_dimid           !netCDF dimension id
    integer :: levsno_dimid           !netCDF dimension id
    integer :: levtot_dimid           !netCDF dimension id
    integer :: levlak_dimid           !netCDF dimension id
    integer :: rtmlon_dimid           !netCDF dimension id
    integer :: rtmlat_dimid           !netCDF dimension id
!
! Netcdf variable Id's 
!
#if (defined RTM)               
    integer :: volr_id                !netCDF variable id
#endif                          
!
! Arrays needed for netCDF output
!
    integer , dimension(:)  , pointer :: ibufslc   !temporary
    real(r8), dimension(:)  , pointer :: rbufslc   !temporary 
    real(r8), dimension(:,:), pointer :: rbufmlc   !temporary  
    integer , dimension(:)  , pointer :: ibufslp   !temporary 
    real(r8), dimension(:)  , pointer :: rbufslp   !temporary 
    real(r8), dimension(:,:), pointer :: rbufmlp   !temporary 
    integer :: begc,endc,numc           !column 1d indices
    integer :: begp,endp,nump           !pft 1d indices
    character(len=16) :: namec          !column 1d name 
    character(len=16) :: namep          !pft 1d name 
!-----------------------------------------------------------------------

    ! create initial conditions file for writing

    if (masterproc) then

       loc_fn = set_init_filename()
       write(6,*)
       write(6,*)'(INICFILEMOD): Writing clm2 initial conditions dataset at ',&
            trim(loc_fn), 'at nstep = ',get_nstep()
       write(6,*)

       ! create new netCDF file (in defined mode)

       call wrap_create (trim(loc_fn), nf_clobber, ncid)

       ! set fill mode to "no fill" to optimize performance

       status = nf_set_fill (ncid, nf_nofill, omode)
       if (status /= nf_noerr) then
          write (6,*) ' netCDF error = ',nf_strerror(status)
          call endrun
       end if

       ! define dimensions 

       call wrap_def_dim (ncid, 'gridcell', grid1d%num, gridcell_dimid)
       call wrap_def_dim (ncid, 'landunit', land1d%num, landunit_dimid)
       call wrap_def_dim (ncid, 'column'  , cols1d%num, column_dimid)
       call wrap_def_dim (ncid, 'pft'     , pfts1d%num, pft_dimid)
       call wrap_def_dim (ncid, 'levsno'  , nlevsno   , levsno_dimid)
       call wrap_def_dim (ncid, 'levsoi'  , nlevsoi   , levsoi_dimid)
       call wrap_def_dim (ncid, 'levlak'  , nlevlak   , levlak_dimid)
       call wrap_def_dim (ncid, 'levtot'  , nlevsno+nlevsoi, levtot_dimid)
#if (defined RTM)
       call wrap_def_dim (ncid, 'rtmlon', rtmlon, rtmlon_dimid)
       call wrap_def_dim (ncid, 'rtmlat', rtmlat, rtmlat_dimid)
#endif

       ! define global attributes

       str = 'CF1.0'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'conventions', trim(str))
       
       call getdatetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call wrap_put_att_text (ncid, NF_GLOBAL, 'history', trim(str))
       
       call getenv ('LOGNAME', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'logname', trim(str))
       
       call getenv ('HOST', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'host', trim(str))
       
       str = 'Community Land Model: CLM2'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'source', trim(str))
       
       str = '$Name$' 
       call wrap_put_att_text (ncid, NF_GLOBAL, 'version', trim(str))
       
       str = '$Id$'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'revision_id', trim(str))
       
       str = ctitle 
       call wrap_put_att_text (ncid, NF_GLOBAL, 'case_title', trim(str))

       str = caseid
       call wrap_put_att_text (ncid, NF_GLOBAL, 'case_id', trim(str))

       str = fsurdat
       call wrap_put_att_text (ncid, NF_GLOBAL, 'surface_dataset', trim(str))

       ! Define current date

       call wrap_def_var (ncid, 'mcdate', nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            'current date as 8 digit integer (YYYYMMDD)')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       call wrap_def_var (ncid, 'mcsec',  nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            'current seconds of current date')
       call wrap_put_att_text (ncid, varid, 'units', 's')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       ! Define 1d gridcell mapping indices

       dim1_id(1) = gridcell_dimid

       call wrap_def_var(ncid, 'grid1d_lon', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            'gridcell longitude')
       call wrap_put_att_text (ncid, varid, 'units',&
            'degrees_east')

       call wrap_def_var(ncid, 'grid1d_lat', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            'gridcell latitude')
       call wrap_put_att_text (ncid, varid, 'units',&
            'degrees_north')

       call wrap_def_var(ncid, 'grid1d_ixy', nf_int, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            '2d longitude index of corresponding gridcell')

       call wrap_def_var(ncid, 'grid1d_jxy', nf_int, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', &
            '2d latitude index of corresponding gridcell')

       ! Define 1d landunit mapping indices

        dim1_id(1) = landunit_dimid

        call wrap_def_var(ncid, 'land1d_lon', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'landunit longitude')
        call wrap_put_att_text (ncid, varid, 'units', &
             'degrees_east')

        call wrap_def_var(ncid, 'land1d_lat', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'landunit latitude')

        call wrap_def_var(ncid, 'land1d_ixy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d longitude index of corresponding landunit')

        call wrap_def_var(ncid, 'land1d_jxy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d latitude index of corresponding landunit')

        call wrap_def_var(ncid, 'land1d_gi ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d grid index of corresponding landunit')

        call wrap_def_var(ncid, 'land1d_wtxy', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'landunit weight relative to corresponding gridcell')

        call wrap_def_var(ncid, 'land1d_itypwat', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'landunit water type')

       ! Define 1d column mapping indices

        dim1_id(1) = column_dimid

        call wrap_def_var(ncid, 'cols1d_lon', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'column longitude')
        call wrap_put_att_text (ncid, varid, 'units', &
             'degrees_east')

        call wrap_def_var(ncid, 'cols1d_lat', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'column latitude')
        call wrap_put_att_text (ncid, varid, 'units', &
             'degrees_north')

        call wrap_def_var(ncid, 'cols1d_ixy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d longitude index of corresponding column')

        call wrap_def_var(ncid, 'cols1d_jxy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d latitude index of corresponding column')

        call wrap_def_var(ncid, 'cols1d_gi ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d grid index of corresponding column')

        call wrap_def_var(ncid, 'cols1d_li ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d landunit index of corresponding column')

        call wrap_def_var(ncid, 'cols1d_wtxy', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'column weight relative to corresponding gridcell')

        call wrap_def_var(ncid, 'cols1d_wtlnd', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'column weight relative to corresponding landunit')

        call wrap_def_var(ncid, 'cols1d_itypwat', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'column water type')

       ! Define 1d pft mapping indices

        dim1_id(1) = pft_dimid

        call wrap_def_var(ncid, 'pfts1d_lon', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft longitude')
        call wrap_put_att_text (ncid, varid, 'units', &
             'degrees_east')

        call wrap_def_var(ncid, 'pfts1d_lat', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft latitude')
        call wrap_put_att_text (ncid, varid, 'units', &
             'degrees_north')

        call wrap_def_var(ncid, 'pfts1d_ixy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d longitude index of corresponding pft')

        call wrap_def_var(ncid, 'pfts1d_jxy', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '2d latitude index of corresponding pft')

        call wrap_def_var(ncid, 'pfts1d_gi ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d grid index of corresponding pft')

        call wrap_def_var(ncid, 'pfts1d_li ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d landunit index of corresponding pft')

        call wrap_def_var(ncid, 'pfts1d_ci ', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             '1d column index of corresponding pft')

        call wrap_def_var(ncid, 'pfts1d_wtxy', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft weight relative to corresponding gridcell')

        call wrap_def_var(ncid, 'pfts1d_wtlnd', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft weight relative to corresponding landunit')

        call wrap_def_var(ncid, 'pfts1d_wtcol', nf_double, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft weight relative to corresponding column')

        call wrap_def_var(ncid, 'pfts1d_itypwat', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid , 'long_name', &
             'pft water type')

        call wrap_def_var(ncid, 'pfts1d_itypveg', nf_int, 1, dim1_id, varid)
        call wrap_put_att_text (ncid, varid, 'long_name', &
             'pft vegetation type')

       ! define single-level fields 

       dim1_id(1) = pft_dimid
       call wrap_def_var (ncid, 'T_VEG', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'vegetation temperature (T_VEG)')
       call wrap_put_att_text (ncid, varid, 'units', 'K')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = column_dimid
       call wrap_def_var (ncid, 'T_GRND', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'ground temperature (T_GRND)')
       call wrap_put_att_text (ncid, varid, 'units', 'K')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = pft_dimid
       call wrap_def_var (ncid, 'H2OCAN', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'canopy water (H2OCAN)')
       call wrap_put_att_text (ncid, varid, 'units', 'kg/m2')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = column_dimid
       call wrap_def_var (ncid, 'H2OSNO', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'snow water (H2OSNO)')
       call wrap_put_att_text (ncid, varid, 'units', 'kg/m2')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = column_dimid
       call wrap_def_var (ncid, 'SNOWDP', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'snow depth (SNOWDP)')
       call wrap_put_att_text (ncid, varid, 'units', 'm')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = column_dimid
       call wrap_def_var (ncid, 'SNOWAGE', nf_double, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'snow age (SNOWAGE)')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim1_id(1) = column_dimid
       call wrap_def_var (ncid, 'SNLSNO', nf_int, 1, dim1_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'number of snow layers (SNLSNO)')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       ! define multi-level fields

       dim2_id(1) = levtot_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'T_SOISNO', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'soil-snow temperature')
       call wrap_put_att_text (ncid, varid, 'units', 'K')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levlak_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'T_LAKE', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'lake temperature')
       call wrap_put_att_text (ncid, varid, 'units', 'K')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levtot_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'H2OSOI_LIQ', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'liquid water')
       call wrap_put_att_text (ncid, varid, 'units', 'kg/m2')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levtot_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'H2OSOI_ICE', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name', 'ice lens')
       call wrap_put_att_text (ncid, varid, 'units', 'kg/m2')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levsno_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'ZSNO', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name','snow layer depth')
       call wrap_put_att_text (ncid, varid, 'units','m')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levsno_dimid
       dim2_id(2) = column_dimid
       call wrap_def_var (ncid, 'DZSNO', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name','snow layer thickness')
       call wrap_put_att_text (ncid, varid, 'units','m')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

       dim2_id(1) = levsno_dimid
       dim2_id(2) = column_dimid   
       call wrap_def_var (ncid, 'ZISNO', nf_double, 2, dim2_id, varid)
       call wrap_put_att_text (ncid, varid, 'long_name','snow interface depth')
       call wrap_put_att_text (ncid, varid, 'units', 'm')
       call wrap_put_att_text (ncid, varid, 'cell_method', 'time: instantaneous')

#if (defined RTM)
       dim2_id(1) = rtmlon_dimid 
       dim2_id(2) = rtmlat_dimid
       call wrap_def_var (ncid, 'RTMVOLR', nf_double, 2, dim2_id, volr_id)
       call wrap_put_att_text (ncid, volr_id, 'long_name','water volumn in cell (volr)')
       call wrap_put_att_text (ncid, volr_id, 'units', 'm3')
       call wrap_put_att_text (ncid, volr_id, 'cell_method', 'time: instantaneous')
#endif

       ! finish creating netCDF file (end define mode)

       status = nf_enddef(ncid)

       ! Write mapping incides and weights
        
       call wrap_inq_varid (ncid, 'grid1d_lon', varid)
       call wrap_put_var_realx(ncid, varid, grid1d%londeg(:))
       call wrap_inq_varid (ncid, 'grid1d_lat', varid)
       call wrap_put_var_realx(ncid, varid, grid1d%latdeg(:))
       call wrap_inq_varid (ncid, 'grid1d_ixy', varid)
       call wrap_put_var_int  (ncid, varid, grid1d%ixy(:))
       call wrap_inq_varid (ncid, 'grid1d_jxy', varid)
       call wrap_put_var_int  (ncid, varid, grid1d%jxy(:))

       call wrap_inq_varid (ncid, 'land1d_lon', varid)
       call wrap_put_var_realx(ncid, varid, land1d%londeg(:))
       call wrap_inq_varid (ncid, 'land1d_lat', varid)
       call wrap_put_var_realx(ncid, varid, land1d%latdeg(:))
       call wrap_inq_varid (ncid, 'land1d_ixy', varid)
       call wrap_put_var_int  (ncid, varid, land1d%ixy(:))
       call wrap_inq_varid (ncid, 'land1d_jxy', varid)
       call wrap_put_var_int  (ncid, varid, land1d%jxy(:))
       call wrap_inq_varid (ncid, 'land1d_gi', varid)
       call wrap_put_var_int  (ncid, varid, land1d%gindex(:))
       call wrap_inq_varid (ncid, 'land1d_wtxy', varid)
       call wrap_put_var_realx(ncid, varid, land1d%wtxy(:))
       call wrap_inq_varid (ncid, 'land1d_itypwat', varid)
       call wrap_put_var_int  (ncid, varid, land1d%itypwat(:))
        
       call wrap_inq_varid (ncid, 'cols1d_lon', varid)
       call wrap_put_var_realx(ncid, varid, cols1d%londeg(:))
       call wrap_inq_varid (ncid, 'cols1d_lat', varid)
       call wrap_put_var_realx(ncid, varid, cols1d%latdeg(:))
       call wrap_inq_varid (ncid, 'cols1d_ixy', varid)
       call wrap_put_var_int  (ncid, varid, cols1d%ixy(:))
       call wrap_inq_varid (ncid, 'cols1d_jxy', varid)
       call wrap_put_var_int  (ncid, varid, cols1d%jxy(:))
       call wrap_inq_varid (ncid, 'cols1d_gi', varid)
       call wrap_put_var_int  (ncid, varid, cols1d%gindex(:))
       call wrap_inq_varid (ncid, 'cols1d_li', varid)
       call wrap_put_var_int  (ncid, varid, cols1d%lindex(:))
       call wrap_inq_varid (ncid, 'cols1d_wtxy', varid)
       call wrap_put_var_realx(ncid, varid, cols1d%wtxy(:))
       call wrap_inq_varid (ncid, 'cols1d_wtlnd', varid)
       call wrap_put_var_realx(ncid, varid, cols1d%wtlnd(:))
       call wrap_inq_varid (ncid, 'cols1d_itypwat', varid)
       call wrap_put_var_int  (ncid, varid, cols1d%itypwat(:))
        
       call wrap_inq_varid (ncid, 'pfts1d_lon', varid)
       call wrap_put_var_realx(ncid, varid, pfts1d%londeg(:))
       call wrap_inq_varid (ncid, 'pfts1d_lat', varid)
       call wrap_put_var_realx(ncid, varid, pfts1d%latdeg(:))
       call wrap_inq_varid (ncid, 'pfts1d_ixy', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%ixy(:))
       call wrap_inq_varid (ncid, 'pfts1d_jxy', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%jxy(:))
       call wrap_inq_varid (ncid, 'pfts1d_gi', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%gindex(:))
       call wrap_inq_varid (ncid, 'pfts1d_li', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%lindex(:))
       call wrap_inq_varid (ncid, 'pfts1d_ci', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%cindex(:))
       call wrap_inq_varid (ncid, 'pfts1d_wtxy', varid)
       call wrap_put_var_realx(ncid, varid, pfts1d%wtxy(:))
       call wrap_inq_varid (ncid, 'pfts1d_wtlnd', varid)
       call wrap_put_var_realx(ncid, varid, pfts1d%wtlnd(:))
       call wrap_inq_varid (ncid, 'pfts1d_wtcol', varid)
       call wrap_put_var_realx(ncid, varid, pfts1d%wtcol(:))
       call wrap_inq_varid (ncid, 'pfts1d_itypwat', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%itypwat(:))
       call wrap_inq_varid (ncid, 'pfts1d_itypveg', varid)
       call wrap_put_var_int  (ncid, varid, pfts1d%itypveg(:))

    endif  !end of if-masterproc block

    ! Write single-level and multi-level time dependent variables

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num
    namec = cols1d%name

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num
    namep = pfts1d%name

    allocate (ibufslc(begc:endc))
    allocate (rbufslc(begc:endc))
    allocate (ibufslp(begp:endp))
    allocate (rbufslp(begp:endp))

    ! Write out zisno
    allocate (rbufmlc(-nlevsno:-1,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno:-1,ci) = c%cps%zi(-nlevsno:-1)  
    end do
    call write_ncd (ncid, 'ZISNO', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out zsno
    allocate (rbufmlc(-nlevsno+1:0,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno+1:0,ci) = c%cps%z(-nlevsno+1:0)  
    end do
    call write_ncd (ncid, 'ZSNO', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out dzsno
    allocate (rbufmlc(-nlevsno+1:0,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno+1:0,ci) = c%cps%dz(-nlevsno+1:0)  
    end do
    call write_ncd (ncid, 'DZSNO', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out h2osoi_liq
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno+1:nlevsoi,ci) = c%cws%h2osoi_liq(-nlevsno+1:nlevsoi)  
    end do
    call write_ncd (ncid, 'H2OSOI_LIQ', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out h2osoi_ice
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno+1:nlevsoi,ci) = c%cws%h2osoi_ice(-nlevsno+1:nlevsoi)  
    end do
    call write_ncd (ncid, 'H2OSOI_ICE', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out t_soisno - column
    allocate (rbufmlc(-nlevsno+1:nlevsoi,begc:endc))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(-nlevsno+1:nlevsoi,ci) = c%ces%t_soisno(-nlevsno+1:nlevsoi) 
    end do
    call write_ncd (ncid, 'T_SOISNO', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out t_lake - column 
    allocate (rbufmlc(1:nlevlak,cols1d%num))
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufmlc(1:nlevlak,ci) = c%ces%t_lake(1:nlevlak)
    end do
    call write_ncd (ncid, 'T_LAKE', rbufmlc, clmlevel=namec)
    deallocate (rbufmlc)

    ! Write out t_veg - pft
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       rbufslp(pi) = p%pes%t_veg  
    end do
    call write_ncd (ncid, 'T_VEG', rbufslp, clmlevel=namep)

    ! Write out t_grnd - column
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufslc(ci) = c%ces%t_grnd  
    end do
    call write_ncd (ncid, 'T_GRND', rbufslc, clmlevel=namec)

    ! Write out h2ocan - pft
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       rbufslp(pi) = p%pws%h2ocan  
    end do
    call write_ncd (ncid, 'H2OCAN', rbufslp, clmlevel=namep)

    ! Write out h2osno - column
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufslc(ci) = c%cws%h2osno  
    end do
    call write_ncd (ncid, 'H2OSNO', rbufslc, clmlevel=namec)

    ! Write out snowdp - column
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufslc(ci) = c%cps%snowdp
    end do
    call write_ncd (ncid, 'SNOWDP', rbufslc, clmlevel=namec)

    ! Write out snowage - column
    do ci = begc,endc
       c => cpoint%col(ci)%c
       rbufslc(ci) = c%cps%snowage  
    end do
    call write_ncd (ncid, 'SNOWAGE', rbufslc, clmlevel=namec)

    ! Write out snlsno - column
    do ci = begc,endc
       c => cpoint%col(ci)%c
       ibufslc(ci) = c%cps%snl  
    end do
    call write_ncd (ncid, 'SNLSNO', ibufslc, clmlevel=namec)

#if (defined RTM)
    ! Write out RTM volr
    if (masterproc) call wrap_put_var_realx (ncid, volr_id, volr)
#endif

    deallocate (ibufslc)
    deallocate (rbufslc)
    deallocate (ibufslp)
    deallocate (rbufslp)

! archive initial conditions dataset (Mass Store currently)

    if (masterproc) then
       call wrap_close (ncid)
       if (mss_irt > 0) then 
          rem_dir = trim(archive_dir) // '/init/'
          rem_fn = set_filename(rem_dir, loc_fn)
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .true.)
       endif
    endif

    return
  end subroutine inicwrt

!=======================================================================
! BEGIN GENERIC PROCEDURE DEFINITIONS
!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_inicwrite
!
! !INTERFACE:
  logical function do_inicwrite()
!
! !DESCRIPTION: 
! Determine if initial dataset is to be written at this time step
!
! !USES:
    use time_manager, only : get_curr_date, get_prev_date
    use clm_varctl, only : hist_crtinic
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: yr         !nstep year (0 -> ...)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: mcsec      !nstep time of day [seconds] 
    integer :: mcsecm1    !nstep-1 time of day [seconds]
!-----------------------------------------------------------------------

    ! Set calendar for current time step and previous time step

    call get_curr_date (yr, mon, day, mcsec) 
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Determine if time to write out initial dataset

    do_inicwrite = .false.
    if (hist_crtinic /= 'NONE') then
       if (hist_crtinic == 'MONTHLY') then
          if (mon /= monm1 .and. monm1 /= -1) do_inicwrite = .true.
       else if (hist_crtinic == 'YEARLY') then
          if (monm1 == 12 .and. mon == 1)  do_inicwrite = .true.
       endif
    endif

  end function do_inicwrite

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_init_filename
!
! !INTERFACE:
  character(len=256) function set_init_filename ()
!
! !DESCRIPTION: 
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine inicwrt in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec) 
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_init_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"

  end function set_init_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_real_2d
!
! !INTERFACE:
  subroutine write_real_2d (ncid, varname, rloc, clmlevel)
!
! !DESCRIPTION: 
! Writes 2d initial real field out to netCDF file
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !variable name
    real(r8), pointer :: rloc(:,:)            !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein     
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !variable id 
    integer :: lb,ub                          !bounds of first dimension 
    integer :: nsize                          !size of global array second dimension         
    real(r8), pointer :: rglob(:,:)           !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       call wrap_inq_varid (ncid, trim(varname), varid)
    endif
#if (defined SPMD)
    lb = lbound(rloc, dim=1)
    ub = ubound(rloc, dim=1)
    nsize = get_num1d (clmlevel)
    if (masterproc) then
       allocate (rglob(lb:ub,nsize))
    endif
    call gather_data_to_master (rloc, rglob, clmlevel=clmlevel)
    if (masterproc) then
       call wrap_put_var_realx (ncid, varid, rglob)
       deallocate (rglob)
    endif
#else
    call wrap_put_var_realx (ncid, varid, rloc)
#endif
    return
  end subroutine write_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_real_1d
!
! !INTERFACE:
  subroutine write_real_1d (ncid, varname, rloc, clmlevel)
!
! !DESCRIPTION: 
! Writes 1d initial real field out to netCDF file
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !variable name
    real(r8), pointer :: rloc(:)              !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !variable id 
    integer :: nsize                          !size of global array
    real(r8), pointer :: rglob(:)             !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       call wrap_inq_varid (ncid, trim(varname), varid)
    endif
#if (defined SPMD)
    nsize = get_num1d (clmlevel)
    if (masterproc) then
       allocate (rglob(nsize))
    endif
    call gather_data_to_master (rloc, rglob, clmlevel=clmlevel)
    if (masterproc) then
       call wrap_put_var_realx (ncid, varid, rglob)
       deallocate (rglob)
    endif
#else
    call wrap_put_var_realx (ncid, varid, rloc)
#endif
    return
  end subroutine write_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_int_1d
!
! !INTERFACE:
  subroutine write_int_1d (ncid, varname, iloc, clmlevel)
!
! !DESCRIPTION: 
! Writes 1d initial integer field out to netCDF file
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !variable name
    integer, pointer :: iloc(:)               !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !variable id
    integer :: nsize                          !size of global array
    integer, pointer :: iglob(:)              !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       call wrap_inq_varid (ncid, trim(varname), varid)
    endif
#if (defined SPMD)
    nsize = get_num1d (clmlevel)
    if (masterproc) then
       allocate (iglob(nsize))
    endif
    call gather_data_to_master (iloc, iglob, clmlevel=clmlevel)
    if (masterproc) then
       call wrap_put_var_int (ncid, varid, iglob)
       deallocate (iglob)
    endif
#else
    call wrap_put_var_int (ncid, varid, iloc)
#endif
    return
  end subroutine write_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_real_2d
!
! !INTERFACE:
  subroutine read_real_2d (ncid, varname, rloc, clmlevel)
!
! !DESCRIPTION: 
! Reads 2d initial real field from netCDF file.
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !input variable name
    real(r8), pointer :: rloc(:,:)            !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !netCDF variable id 
    integer :: lb,ub                          !bounds of first dimension 
    integer :: nsize                          !size of global array second dimension         
    real(r8), pointer :: rglob(:,:)           !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       lb = lbound(rloc, dim=1)
       ub = ubound(rloc, dim=1)
       nsize = get_num1d (clmlevel)
       allocate (rglob(lb:ub,nsize))
       call wrap_inq_varid (ncid, varname, varid)
       call wrap_get_var_realx (ncid, varid, rglob)
    endif
#if (defined SPMD)
    call scatter_data_from_master (rloc, rglob, clmlevel=clmlevel)
#else
    rloc(:,:) = rglob(:,:)
#endif
    if (masterproc) then
       deallocate (rglob)
    endif
    return
  end subroutine read_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_real_1d
!
! !INTERFACE:
  subroutine read_real_1d (ncid, varname, rloc, clmlevel)
!
! !DESCRIPTION: 
! Reads 1d initial real field from netCDF file.
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !input variable name
    real(r8), pointer :: rloc(:)              !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !netCDF variable id 
    integer :: nsize                          !size of global array
    real(r8), pointer :: rglob(:)             !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (rglob(nsize))
       call wrap_inq_varid (ncid, varname, varid)
       call wrap_get_var_realx (ncid, varid, rglob)
    endif
#if (defined SPMD)
    call scatter_data_from_master (rloc, rglob, clmlevel=clmlevel)
#else
    rloc(:) = rglob(:)
#endif
    if (masterproc) then
       deallocate (rglob)
    endif
    return
  end subroutine read_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_int_1d
!
! !INTERFACE:
  subroutine read_int_1d (ncid, varname, iloc, clmlevel)
!
! !DESCRIPTION: 
! Reads 1d initial integer field from netCDF file.
! Part of generic output interface
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid               !input unit
    character(len=*), intent(in) :: varname   !input variable name
    integer, pointer :: iloc(:)               !local input data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          !netCDF variable id 
    integer :: nsize                          !size of global array
    integer, pointer :: iglob(:)              !global input data
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = get_num1d (clmlevel)
       allocate (iglob(nsize))
       call wrap_inq_varid (ncid, varname, varid)
       call wrap_get_var_int (ncid, varid, iglob)
    endif
#if (defined SPMD)
    call scatter_data_from_master (iloc, iglob, clmlevel=clmlevel)
#else
    iloc(:) = iglob(:)
#endif
    if (masterproc) then
       deallocate (iglob)
    endif
    return
  end subroutine read_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_num1d
!
! !INTERFACE:
  integer function get_num1d (clmlevel)
!
! !DESCRIPTION: 
! Determine 1d size from clmlevel type
!
! !ARGUMENTS:
  implicit none
  character(len=*), intent(in) :: clmlevel      !type of clm 1d array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
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
       write(6,*) 'GET1DSIZE does not match level: ', trim(clmlevel)
       call endrun
    end select
    return
  end function get_num1d

end module inicFileMod
