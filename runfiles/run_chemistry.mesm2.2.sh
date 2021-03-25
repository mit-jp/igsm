#!/bin/csh
#
# A shell to generate the name list files for chemistry run
#	and run the model
#
# Chemistry runs statrs from the end (year 2005) of historical run
#================================
source /etc/profile.d/modules.csh
 module load pgi/17.3
 module load netcdf

#Kov - diffusion coefficient for ocean temperature
#sen - climate sensitivity
#FA - strengh of aerosol forcing
#rkch - of ratio coefficient for carbon to Kov
#kc - tem parameter
#kcn - tem parameter
#scenario emissions scenario
#run number is run.yy
#run0.yy0 number of histircal run from which chemistry run starts

 set Kov = 3.20
 set sen = 3.19
 set FA = -0.24
 set rkch = 0.94
 set kc =  400
 set kcn =  150
 set scenario =  e5_Ref_820
 set run0 = 524
 set run = 544
set yy = 21
set yy0 = 21


 setenv BASEDIR /home/sokolov/MESM2.2/RELEASE/igsm/

echo "run_chemistry.mesm2.2.sh "



# first year of the run
@ INYEAR=2006
# last year of the run
@ LYEAR=2100

@ YEARGT=1860
setenv CO2TEM 286.4

@ INYEARM1 = $INYEAR - 1
@ IYTEM = $INYEAR - 1
@ LYEAR0 = $INYEAR - 1
@ LYTEM= $LYEAR 
@ LYEARCLM = $LYEAR + 1
@ NYEARS = $LYEARCLM - $INYEAR

setenv MODDIR $BASEDIR

setenv OUTDIR  /scratch/sokolov
# STRDIR0 directory where historical runs is stored
setenv STRDIR0  /net/fs10/d1/sokolov/MESM2.2/RELEASE/

# STRDIR directory where chemistry runs is stored
setenv STRDIR  /net/fs10/d1/sokolov/MESM2.2/RELEASE/
 setenv OUTDIR $STRDIR

setenv DATADIR $BASEDIR/datanew/
setenv TEMDATA $DATADIR/tem/44d/
setenv CLMDATA $DATADIR/clm/3.5/
setenv EMIDIR $DATADIR/emimed/

echo $EMIDIR

 setenv EMD $EMIDIR/edaily.$scenario.dat
 setenv EMU $EMIDIR/urban.$scenario.dat
 setenv SO2ER $EMIDIR/edailyso2.$scenario

# --- indicate system name
setenv PLATFORM Linux

#
# --- you can give a specific label for the run
#	option: change it manually after creating name.dat
#
setenv RUNLABEL ${run}.$yy
setenv RUNLABEL0 ${run0}.$yy0
setenv CS $sen
setenv KV $Kov
setenv TRR  FALSE
setenv TRR  TRUE

#   RESTART is false for a new run and true if this is a restart of 
#     old run
#
setenv RESTART  TRUE
setenv RESTART  FALSE

if ($RESTART == "FALSE") then

cd $OUTDIR
if (-e $RUNLABEL) then
        echo $RUNLABEL already exists!  Exitting run_igsm!
        exit
endif
mkdir $RUNLABEL

else

setenv OUTDIR  $STRDIR

endif
setenv OUTDIR $OUTDIR/$RUNLABEL/
setenv EXEDIR $OUTDIR/

#
# --- if you have something to say (before it's too late)
#	here are your two chances
#
setenv WORDS1 "Climate chemistry model form $INYEAR to $LYEAR for $scenario "
setenv WORDS2 " S=$CS Kvh=$KV KC=$kc KCN=$kcn FA=$FA"

# === write name.dat: 
#	1. basically you can put junk inthe first line after RUNLABEL
#	2. modify parameters you'd like to change (INCLUDE/ctrparam.h too)
#	3. carefully change the directory names for initial files
#
#	================================================================
#	Description of some parameters:
#	-------------------------------
#	INYEAR:		Start year
#	LYEAR:		Last year
#	CO2:		 1: grid point gases provided by chemistry model
#			-7: off line calculation 
#	READGHG:	1 for read 0 for not read gases
#	cfdiff:		coefficient of heat diffussion into deep ocean
#	TRANSR:		if diffussion into deep ocean will be carried out
#	ISTRT1:		0 for initialize a run
#			1 for initermidiate run
#
#	================================================================


cd $EXEDIR

setenv IRST 1

	if ($RESTART == "FALSE") then

setenv IRST 0


 setenv POVDIR $STRDIR0/$RUNLABEL0

 cp -i $POVDIR/pov1.$IYTEM.$RUNLABEL0 $OUTDIR/pov1$RUNLABEL
 cp -i $POVDIR/pov2.$IYTEM.$RUNLABEL0 $OUTDIR/pov2$RUNLABEL
 cp -i $POVDIR/povocean.$IYTEM.$RUNLABEL0 $OUTDIR/povocean$RUNLABEL
 cp -i $POVDIR/tveg43d.ecd  $OUTDIR/tveg43d.ecd


 cp -i $POVDIR/T2M*  $OUTDIR/T2M.$RUNLABEL
 cp -i $POVDIR/TGO3_2D*  $OUTDIR/TGO3_2D.$RUNLABEL
 cp -i $POVDIR/QFLX_2D*  $OUTDIR/QFLX_2D.$RUNLABEL

echo $POVDIR

cat >! eppaemis.dat << EOF
 &EPPA
    FORMSO2= .FALSE.
    FORMSO2= .TRUE.
    emiss_data= '$EMD'
    SO2ERATIO= '$SO2ER'
    LYEAREM=$LYEAR
 &END

EOF


        endif
cat >! couple.nml << EOF
 &COUPLE_PARM
  dtcouple=24.0
  dtcouple=1.0
  dtatm=1.0
  dtocn=1.0
  startYear = $INYEAR
  endYear = $LYEAR
  taveDump=100
 &END

EOF


cat >! name.dat << EOF
EXP. $RUNLABEL $WORDS1
$WORDS2
 &INPUTZ
    AEXP     = $RUNLABEL,
    INYEAR   = $INYEAR,
    LYEAR    = $LYEARCLM,
    CO2      = 1.0,
    YEARGT   = $YEARGT
    AERFOR     = $FA,
    AERF4BC    = -0.35,
    CLDFEED  = .false.,
    CLDFEED  = .true.,
    fclmlice  = '$CLMDATA/clm.pctglacier.txt',
    fbaresoil = '$CLMDATA/clm.pctbaresoil.txt',
    fwmax     = '$CLMDATA/clm.wmax.txt',
    wetfrac_data     = '$CLMDATA/WETLANDS.TEM44b.dat'
    fprratio  = '$CLMDATA/pcp2lnd-data.001x046.txt',
    o3datadir = '$TEMDATA/O3/',
    clmsen   = $CS,
    TRANSR   = $TRR,
    cfdif0   =  $KV,
    rkv      = 1.00,
    cfocdif = $rkch
    CONTRR   = .true.,
    ISTRT1   = $IRST,
    ISTRTCHEM = 1,
    plotfl   = 'plot$RUNLABEL',
    nwrfl    = 'nwr$RUNLABEL',
    file1    = 'pov1$RUNLABEL',
    file2    = 'pov2$RUNLABEL',
    pov_deepo= 'pov_deepo$RUNLABEL',
    deepco2in= '$POVDIR/init_deepo.$IYTEM.$RUNLABEL0',
    fl_init_alkt= '$DATADIR/data/46lat/init_alkt_steph1'
    fl_init_salt= '$DATADIR/data/46lat/init_salt_steph1'
    caruptfile= '$OUTDIR/carupt.$RUNLABEL',
    fnememiss = '$OUTDIR/nememiss.$RUNLABEL',
    flrco2av= 'rco2av_sp$RUNLABEL',
    flin_nep = '$STRDIR0/$RUNLABEL0/init_nep.$IYTEM.$RUNLABEL0',
    last_nep = 'pov_nep$RUNLABEL',
    dirdat1  = '$DATADIR/data/',
    dirdat2  = '$DATADIR/data/46lat/',
    qffile   = '$DATADIR/data/46lat/Z1OAVE46_ZN.dat'
    zmfile   = '$DATADIR/data/46lat/Z1OMAX.dat',
    tsfile   = 'T2M.$RUNLABEL',
    bgrghg_data = '$DATADIR/forcing/ghgdata.GISS.1765-2016.dat'
    co2_data = '$DATADIR/forcing/ghgdata.GISS.1765-2016.dat'
    S0C_data = '$DATADIR/forcing/KoppLean.annual.TSI.1610-2013.dat'
    fl_volaer= '$DATADIR/forcing/STRATAER.4AR5.46.1850_2011.dat',
    o3_data  = '$DATADIR/forcing/o3.1850-2009.46x14_mass.dat'
    sulf1986 = '$DATADIR/forcing/sulfate.4x5.1986.dat',
    SO2_EM   = '$DATADIR/forcing/Total_SO2_EM.EPPA.txt',
    oco2file  = '$DATADIR/forcing/ghgdata.GISS.1765-2016.dat'
    co2rfile = '$DATADIR/data/co2_norm_46.new',
    SEN_dat =  '$DATADIR/data/SEN_TABLE.txt'
    chemdata = '$DATADIR/data/'
    chemrstfl = '$DATADIR/init/chemrenew.$IYTEM'
    chemout =  '$OUTDIR/'
 &END

EOF

cat >! name_ocean.dat << EOF

 &INPUTON
    dirdat2   = '$DATADIR/data/46lat/',
    qffile    = 'QFLX_2D.$RUNLABEL',
    kvnewfl    = '$DATADIR/data/46lat/KVNEWZON.dat',
    zmmaxfl   = 'Z1OMAX46_2D.dat',
    zmavefl   = 'Z1OAVE46_2D.dat',
    povocean  =  'povocean$RUNLABEL',
    TRRUN     = $TRR,
    CTRUN     = .true.,
    flto3     = 'TGO3_2D.$RUNLABEL',
    cdiff     = $KV
 &END


EOF
    
 @ NELAPSE = $NYEARS * 365 * (-1)
 @ STRTYMD = ( $INYEAR * 100 + 1 ) * 100 + 1
    

cat >! lnd.stdin << EOF


 &clm_inparm
 caseid         = '$RUNLABEL'
 ctitle         = '$WORDS1'
 finidat        = ''
 finidat        = '$POVDIR/$RUNLABEL0.clm2.i.$INYEAR-01-01-00000.nc'
 fsurdat        = ''
 fsurdat        = '$CLMDATA/surfdata_0046x0001_igsmveg_ar5_historical_1860.nc'
 fpftdyn        = ''
 fatmgrid        = '$CLMDATA/griddata_0046x0001.nc'
 fatmlndfrc        = '$CLMDATA/fracdata_0046x0001.nc'
 fstochdat      = '$CLMDATA/stochastic-data.001x046.GLS.nc'
 fpcp2pft       = '$CLMDATA/pcp2pft-data.001x046.AR5.nc'
 fpftcon        = '$CLMDATA/pft-physiology.ar5'
 dynamic_pft    = .true.
 dynamic_pft    = .false.
 rampYear_dynpft = $YEARGT
 rampYear_dynpft = 0
 orbitfix       = .true.
 orbityr        =   $YEARGT
 orbityr        =   2000
 offline_atmdir = '$CLMDATA/atmdat.nc'
 rpntpath       = 'pntr$RUNLABEL'
 frivinp_rtm    = '$CLMDATA/rdirc.05.061026'
 nsrest         =  0
 nrevsn         = '$STRDIR0$RUNLABEL0//$RUNLABEL0.clm2.rst.nc'
 nelapse        =  $NELAPSE
 dtime          =  3600
 start_ymd      =  $STRTYMD
 start_tod      =  0
 irad           = -1
 wrtdia         = .true.
 wrtdia         = .false.
 mss_irt        =  0
 hist_dov2xy    = .false.
 hist_dov2xy    = .true.
 hist_nhtfrq    =  0
 hist_nhtfrq    =  -8760
 hist_mfilt     =  1
 hist_ndens     =  2
 hist_crtinic   = 'YEARLY'
 /


EOF



cat >! tem.go << EOF
0
1
$POVDIR/temstate$IYTEM.kc$kc.kcn$kcn.binary.$RUNLABEL0
$LYEAR0
$INYEAR
1
46
1
92
1
1
1
1
1
1
$CO2TEM
$CO2TEM
2
1
0
1
$TEMDATA/igsmcomm44c.ecd
1
1
1
1
1
1
1
10
3000
3
0.01
1.0
0.02
$INYEAR
$LYTEM
1
0.01
20
100
$TEMDATA/tsoil43d.ecd
$TEMDATA/trotz43d.ecd
tveg43d.ecd
$TEMDATA/tleaf43d.ecd
$TEMDATA/tmcrv43d.ecd
$TEMDATA/ag43d.ecd
0
$TEMDATA/nemn2o437clm3544d.ecd
12
$TEMDATA/ice43801.dat
$TEMDATA/tpd43801.dat
$TEMDATA/tmt43801.dat
$TEMDATA/bon43801.dat
$TEMDATA/hvc43801.dat
$TEMDATA/hvd43801.dat
$TEMDATA/pwn43801.dat
$TEMDATA/cur43801.dat
$TEMDATA/duc43801.dat
$TEMDATA/gua43801.dat
$TEMDATA/tai43801.dat
$TEMDATA/gua2_43801.dat
0.001038
1
6
npp
cflux
nep
ch4flux
n2oflux
rh
NPP.GLB
CFLX.GLB
NEP.GLB
CH4.GLB
N2O.GLB
RH.GLB
1
$OUTDIR
1
temstate$LYTEM.kc$kc.kcn$kcn.binary.$RUNLABEL
$LYTEM




EOF



cp name_ocean.dat $OUTDIR/name_ocean_sp.$RUNLABEL

cp name.dat $OUTDIR/name_sp.$RUNLABEL




if ($RESTART == "FALSE") then

         cp $MODDIR/build_chemistry/src/igsm22 mod_$RUNLABEL.out

        (./mod_$RUNLABEL.out < tem.go) >>& $OUTDIR/res$RUNLABEL
#       (./mod_$RUNLABEL.out < tem.go) 

 cp -f $OUTDIR/pov_deepo$RUNLABEL $OUTDIR/init_deepo.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov_nep$RUNLABEL $OUTDIR/init_nep.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov1$RUNLABEL $OUTDIR/pov1.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov2$RUNLABEL $OUTDIR/pov2.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/povocean$RUNLABEL $OUTDIR/povocean.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/chemrenew.dat $OUTDIR/chemrenew.$LYEAR.$RUNLABEL

        cd ../
        mv -i $RUNLABEL $STRDIR
else
       ./mod_$RUNLABEL.out  >>& $OUTDIR/res$RUNLABEL

 cp -f $OUTDIR/pov_deepo$RUNLABEL $OUTDIR/init_deepo.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov_nep$RUNLABEL $OUTDIR/init_nep.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov1$RUNLABEL $OUTDIR/pov1.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov2$RUNLABEL $OUTDIR/pov2.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/povocean$RUNLABEL $OUTDIR/povocean.$LYEAR.$RUNLABEL


endif


exit 
