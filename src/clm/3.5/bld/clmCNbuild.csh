#! /bin/csh -f

#-----------------------------------------------------------------------
## Default CLM Build
##------------
##
setenv MP_NODES 2
setenv MP_TASKS_PER_NODE 1
setenv MP_EUILIB us
setenv MP_RMPOOL 1

unsetenv MP_PROCS

setenv MP_STDINMODE 0

# should be set equal to (CPUs-per-node / tasks_per_node)
# Only activated if smp=on below
setenv OMP_NUM_THREADS 1

## suggestion from Jim Edwards to reintroduce XLSMPOPTS on 11/13/03
setenv XLSMPOPTS "stack=256000000"
setenv AIXTHREAD_SCOPE S
setenv MALLOCMULTIHEAP true
setenv OMP_DYNAMIC false
limit stacksize unlimited
limit datasize unlimited
## Do our best to get sufficient stack memory

## netCDF stuff
setenv INC_NETCDF /home/adam/bin/netcdf/x86-64/include
setenv LIB_NETCDF /home/adam/bin/netcdf/x86-64/lib
setenv INC_MPI /usr/pkg/x86-64/pkg/openmpi/openmpi-1.1.2/include
setenv LIB_MPI /usr/pkg/x86-64/pkg/openmpi/openmpi-1.1.2/lib

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = /home/adam/CLM/V3.5

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
## setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2       # (MAKE SURE YOU CHANGE THIS!!!)

## Configuration settings:
set spmd     = off     # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 4        # settings are 4->17               (default is 4)
set rtm      = on      # settings are [on   | off       ] (default is off)   
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
## (see below)
#--------------------------------------------------------------------------------------------

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case    = clmrun
set wrkdir  = /home/adam/CLM/V3.5
set blddir  = $wrkdir/exe
set rundir  = $wrkdir/run
set cfgdir  = $clmroot/bld
set usr_src = $clmroot/bld/usr.src

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -rtm $rtm -usr_src $usr_src -bgc cn"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"

##
## Invoke some compiler overrides
##
set flags = "$flags -fopt -O2"


echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1
## Check if config_cache.xml file exists -- if so just run make -- if NOT then run configure.
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
#--------------------------------------------------------------------------------------------
##
##  First, clear IGSM2 specific files
##
echo "Configuring default CLM compile"
rm /home/adam/CLM/V3.5/src/main/clm_stoch.F90
rm /home/adam/CLM/V3.5/src/main/pcp2pftMod.F90
rm /home/adam/CLM/V3.5/src/main/program_igsm.F90
rm /home/adam/CLM/V3.5/src/main/readIGSM2.F90
rm /home/adam/CLM/V3.5/src/main/stochrdMod.F90
cp /home/adam/CLM/V3.5/src/main/Default/* /home/adam/CLM/V3.5/src/main/.
cp /home/adam/CLM/V3.5/src/biogeophys/Default/* /home/adam/CLM/V3.5/src/biogeophys/.

if ( ! -f $blddir/config_cache.xml ) then
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 1
    echo "Building CLM in $blddir ..."
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 1
else
    echo "Re-building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! REMAKE.out      || echo "CLM build failed: see $blddir/REMAKE.out" && exit 1
endif

setenv LID "`date +%y%m%d-%H%M%S`"
echo "Completed CLM Build at $LID"

exit 0
