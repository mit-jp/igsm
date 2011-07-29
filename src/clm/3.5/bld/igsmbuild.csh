#! /bin/csh 

#-----------------------------------------------------------------------
## IGSM2 Build
##------------
##

## suggestion from Jim Edwards to reintroduce XLSMPOPTS on 11/13/03
## Do our best to get sufficient stack memory
limit stacksize unlimited
limit datasize unlimited

## module loads
#source /etc/profile.d/modules.csh
#module load pgi
#module load netcdf


## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = /home/jscott/CLM35_BUILD_mangled_11_1980

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
## setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2       # (MAKE SURE YOU CHANGE THIS!!!)

## Configuration settings:
set spmd     = off      # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 16       # settings are 4->17               (default is 4)
set rtm      = off      # settings are [on   | off       ] (default is off)   
set stoch    = on       # settings are [ on | off ]        (default is off)
set pcp2pft  = on       # settings are [ on | off ]        (default is off)
set mit2d    = on       # settings are [ on | off ]        (default is off)
set tem      = off       # settings are [ on | off ]        (default is off)

## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
## (see below)
#--------------------------------------------------------------------------------------------

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case    = clmrun
set wrkdir  = /home/jscott/CLM35_BUILD_mangled_11_1980
set blddir  = $wrkdir/exe
set rundir  = $wrkdir/run
set cfgdir  = $clmroot/bld
set usr_src = $clmroot/bld/usr.src

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -rtm $rtm -usr_src $usr_src -stoch $stoch -pcp2pft $pcp2pft -mit2d $mit2d -tem $tem"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"


echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1
## Check if config_cache.xml file exists -- if so just run make -- if NOT then run configure.
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
#--------------------------------------------------------------------------------------------
##
##  First, copy IGSM2 specific files
##
echo "Copying IGSM2 specific files to CLM source directories"
cp /home/jscott/CLM35_BUILD_mangled_11_1980/src/main/IGSM2/* /home/jscott/CLM35_BUILD_mangled_11_1980/src/main/.
cp /home/jscott/CLM35_BUILD_mangled_11_1980/src/biogeophys/IGSM2/* /home/jscott/CLM35_BUILD_mangled_11_1980/src/biogeophys/.

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
