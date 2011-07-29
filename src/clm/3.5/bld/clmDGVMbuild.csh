#! /bin/csh -f

#-----------------------------------------------------------------------
## Default CLM Build
##------------
##
limit stacksize unlimited
limit datasize unlimited
## Do our best to get sufficient stack memory

## netCDF stuff
setenv INC_NETCDF /home/pkg/netcdf/netcdf-3.6.3/x86_64/include
setenv LIB_NETCDF /home/pkg/netcdf/netcdf-3.6.3/x86_64/lib

## MPI stuff
setenv INC_MPI /home/pkg/openmpi/openmpi-1.3/x86_64/include/
setenv LIB_MPI /home/pkg/openmpi/openmpi-1.3/x86_64/lib/

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = /home/adam/CLM/V3.5

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
## setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2       # (MAKE SURE YOU CHANGE THIS!!!)

## Configuration settings:
set spmd     = on     # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 10        # settings are 4->17               (default is 4)
set rtm      = off      # settings are [on   | off       ] (default is off)   
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
set flags = "-maxpft $maxpft -rtm $rtm -usr_src $usr_src -bgc dgvm"
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
