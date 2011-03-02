#! /bin/csh -f

##-------------------------------------------------------------------------------
# The following script creates a CLM offline model executable
# at T42 model resolution using RTM river routing
##-------------------------------------------------------------------------------


##-------------------------------------------------------------------------------
## USER MODIFICATION:
## Define necesessary environment and build namelist
## Some of the environment variables will be used by the Makefile
##-------------------------------------------------------------------------------

# 1. Set root directory for CLM source distribution

# old value -- Paul Kishimoto
#setenv ROOTDIR /home/jscott/adam/CLM/4MIT2D
setenv ROOTDIR $PWD
setenv CLMDIR $ROOTDIR

# 2. Set directory to build the executable in

# commented -- Paul Kishimoto
#setenv MODEL_EXEDIR /home/jscott/CLM_BUILD64_TEM/exe

# 3. Set NTASKS and NTHREADS
# NTHREADS is the number of Open-MP threads.  
# ntasks is the number of MPI tasks (if set to 1, SPMD is off)   
# This will be used to also set the environment variables SMP and SPMD
# needed by the Makefile

setenv OS Linux
if ($OS == 'IRIX64') then
  # Set either NTASKS OR NTHREADS to greater than 1 (not both)     
  # Best performance is obtained in pure MPI mode (NTASKS>1, NTHREADS=1)
  setenv NTASKS  4
  setenv NTHREADS 1
else if ($OS == 'AIX') then
  setenv NTASKS  2 
  setenv NTHREADS  4 
else if ($OS == 'Linux') then
  ## On the PC assume a pure shared-memory configuration 
#CAS
#CAS -- set ntasks>1 for svante cluster parallel runs
#CAS
  setenv NTASKS 1 
  setenv NTHREADS 1 
#CAS
#CAS -- set nthreads=1 for svante cluster test
#CAS
else if ($OS == 'OSF1') then
  setenv NTASKS 2 
  setenv NTHREADS 4 
else if ($OS == 'SunOS') then
  ## Don't use OpenMP on Solaris as it currently causes problems
  setenv NTASKS 1 
  setenv NTHREADS 1 
else
   echo "ERROR: $OS is unsupported operating system"
   exit 1
endif

# 4. Set model resolution and determine 



# 5. Set required library and include paths

if ($OS == 'IRIX64') then
  # SGI-specific stuff. MIPSpro is required for f90 (NCAR specific)
  # NOTE: $MPT_SGI is obtained from invoking module mpt
  source /opt/modules/modules/init/csh
  module purge
  module load MIPSpro mpt nqe modules
  setenv INC_NETCDF  /usr/local/include
  setenv LIB_NETCDF  /usr/local/lib
  setenv INC_MPI    $MPT_SGI/usr/include
  setenv LIB_MPI    $MPT_SGI/usr/lib64
else if ($OS == 'AIX') then
  # LIB_MPI and INC_MPI do not need to be see for the IBM-SP
  setenv INC_NETCDF /usr/local/include
  setenv LIB_NETCDF /usr/local/lib32/r4i4
else if ($OS == 'Linux') then
#   setenv INC_NETCDF  /home/pkg/netcdf/netcdf-3.6.2/x86_64/include
#   setenv LIB_NETCDF  /home/pkg/netcdf/netcdf-3.6.2/x86_64/lib
#   setenv INC_MPI    /usr/home/pkg/openmpi/openmpi-1.3.3/x86_64/include
#   setenv LIB_MPI   /usr/home/pkg/openmpi/openmpi-1.3.3/x86_64/lib 
else if ($OS == 'SunOS') then
  setenv INC_MPI    /contrib/include
  setenv LIB_MPI    /contrib/lib
  setenv INC_NETCDF /usr/local/include    
  setenv LIB_NETCDF /usr/local/lib    
else if ($OS == 'OSF1') then
  # LIB_MPI and INC_MPI do not need to be see for the ORNL Compaq 
  setenv INC_NETCDF /usr/local/include
  setenv LIB_NETCDF /usr/local/lib64/r4i4
endif

# 6. Set root directory for input data

setenv CSMDATA $CLMDIR/inputdata

# 7. Create an input parameter namelist file

# commented -- Paul Kishimoto
#mkdir -p $MODEL_EXEDIR;
#cd $MODEL_EXEDIR
#cat >! lnd.stdin <<EOF
# &clmexp
# caseid         = 'clm2.1_run'
# ctitle         = 'clm2.1_run'
# offline_atmdir = '$CSMDATA/NCEPDATA'
# finidat        = ' '
# fsurdat        = '$CSMDATA/srfdata/cam/clms_64x128_T42_c020514.nc'
# fpftcon        = '$CSMDATA/pftdata/pft-physiology'
# frivinp_rtm    = '$CSMDATA/rtmdata/rdirc.05'
# nsrest         =  0
# nelapse        =  48
# dtime          =  180
# start_ymd      =  19971231
# start_tod      =  0
# irad           = -1
# wrtdia         = .true.
# mss_irt        =  0
# hist_dov2xy(1) = .true.
# hist_nhtfrq(1) =  -24
# hist_mfilt(1)  =  1
# hist_ndens(1)  =  1
# hist_crtinic   = 'MONTHLY'
# mksrf_fvegtyp  = '$CSMDATA/rawdata/mksrf_pft.nc'
# mksrf_fsoitex  = '$CSMDATA/rawdata/mksrf_soitex.10level.nc'
# mksrf_fsoicol  = '$CSMDATA/rawdata/mksrf_soicol_clm2.nc'
# mksrf_flanwat  = '$CSMDATA/rawdata/mksrf_lanwat.nc'
# mksrf_fglacier = '$CSMDATA/rawdata/mksrf_glacier.nc'
# mksrf_furban   = '$CSMDATA/rawdata/mksrf_urban.nc'
# mksrf_flai     = '$CSMDATA/rawdata/mksrf_lai.nc'
# mksrf_offline_fnavyoro = '$CSMDATA/rawdata/mksrf_navyoro_20min.nc'
# /
#EOF

##=======================================================================
## In general do not need to edit below this point
##=======================================================================

##------------------------------------------------------
## build preproc.h in ./obj directory
## define OFFLINE if offline mode
## define LSMLON  to number of longitudes on land grid
## define LSMLAT  to number of latitudes  on land grid
## define RTM     if using RTM river routing
## define STOCHASTIC if using stochastic forcing routine
## define PCP2PFT if using precip. distr. across PFTs
## define COUP_TEM if coupling with TEM
##------------------------------------------------------
#cat >! preproc.h <<EOF
#define OFFLINE
#define LSMLON $LSMLON
#define LSMLAT $LSMLAT
#define COUP_MIT2D
#define ONE_RESTART
#define PCP2PFT
#define STOCHASTIC
#define COUP_TEM
#EOF
##define COUP_TEM adam had this defined in his build??

##-----------------------------------------------------
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## OTHER COMPILE OPTIONS FOR IGSM2 USAGE
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##-----------------------------------------------------

setenv PCPTREND notdo

##------------------------------------------------------
## build misc.h in ./obj directory
## define SPMD    if running in message passing mode
##------------------------------------------------------
if ($NTASKS > 1) then
  set spmd = "#define SPMD"
else
  set spmd = "#undef  SPMD"
endif
cat >! misc.h <<EOF
$spmd
EOF

##------------------------------------------------------
## Set environment variable SMP using value of NTHREADS
## Set to TRUE to enable building in SMP mode (uses OpenMP).  
## Currently implemented for IBM, SGI, linux-pgf90. 
## (default is TRUE on IBM and linux-pgf90, and depends on SPMD setting on SGI).
##------------------------------------------------------
if ($NTHREADS > 1) then
  setenv SMP TRUE
else
  setenv SMP FALSE
endif

##------------------------------------------------------
## build Filepath 
##------------------------------------------------------
#cat >! Filepath <<EOF
#$ROOTDIR/main
#$ROOTDIR/biogeophys
#$ROOTDIR/ecosysdyn
#$ROOTDIR/riverroute
#$ROOTDIR/mksrfdata
#$ROOTDIR/atmclm_share
#$ROOTDIR/csm_share
#$ROOTDIR/utils/timing
#EOF

##------------------------------------------------------
## build executable
##------------------------------------------------------

#  Copy source files appropriate for TEM437 coupling
cp main/TEM437/histFileMod.F90 .
cp main/TEM437/program_igsm.F90 .
cp main/TEM437/program_off.F90 .

# Copy source files appropriate for pcp frequency trend uncertainty
if ( $PCPTREND == 'do' ) then
  echo ' Copying PCP trend files '
  cp PCPFREQ/atmdrvMod.F90 .
  cp PCPFREQ/clm_varsur.F90 .
  cp PCPFREQ/surfFileMod.F90 .
endif

#make

