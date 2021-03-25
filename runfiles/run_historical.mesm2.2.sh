#!/bin/csh 
#
# A shell to generate the name list files for run with historical forcing
#	and run the model
#
#================================

source /etc/profile.d/modules.csh
 module load pgi/17.3
 module load netcdf/4 openmpi/2.1.1

#Kov - diffusion coefficient for ocean temperature
#sen - climate sensitivity
#FA - strengh of aerosol forcing
#rkch - of ratio coefficient for carbon to Kov
#kc - tem parameter
#kcn - tem parameter
#run number is run.yy


 set Kov = 3.20
 set sen = 3.19
 set FA = -0.24
 set rkch = 0.94
 set kc =  400
 set kcn =  150
 set run = 634
 set yy = 21

 setenv BASEDIR /home/sokolov/MESM2.2/RELEASE/igsm/

echo " run_pbs_mesm2.2.sp.rr.sh"


# first year of the run
@ INYEAR=1861
# last year of the run +1
@ LYEAR=2005

@ YEARGT=1860
setenv CO2EQ 286.4
setenv CO2TEM $CO2EQ

@ INYEARM1 = $INYEAR - 1
@ LYTEM= $LYEAR 
@ LYEARCLM = $LYEAR + 1
@ NYEARS = $LYEARCLM - $INYEAR

setenv MODDIR $BASEDIR

setenv OUTDIR  /scratch/sokolov/
setenv STRDIR  /net/fs10/d1/sokolov/MESM2.2/RELEASE/
 setenv OUTDIR $STRDIR
setenv DATADIR $BASEDIR/datanew/
setenv TEMDATA $DATADIR/tem/44d/
setenv CLMDATA $DATADIR/clm/3.5/


# --- indicate system name
setenv PLATFORM Linux

#
# --- you can give a specific label for the run
#	option: change it manually after creating name.dat
#
setenv RUNLABEL ${run}.$yy
setenv CS $sen
setenv KV $Kov
setenv TRR  FALSE
setenv TRR  TRUE
#TRR deffines whether heat is diffused inti deep ocean

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
#
 setenv WORDS1 "Run with histirical forcing from $INYEAR to $LYTEM "
setenv WORDS2 " S=$CS Kvh=$KV KC=$kc KCN=$kcn FA=$FA"



setenv IRST 1

	if ($RESTART == "FALSE") then

setenv IRST 0


#        initial data for 1860
 setenv POVDIR  $DATADIR/init

 cp -i $POVDIR/pov1* $OUTDIR/pov1$RUNLABEL
 cp -i $POVDIR/pov2* $OUTDIR/pov2$RUNLABEL
 cp -i $POVDIR/povocean* $OUTDIR/povocean$RUNLABEL
 cp -i $POVDIR/*clm2.i*  $OUTDIR/$RUNLABEL.clm2.i.nc
 cp -i $POVDIR/temstate.kc$kc.kcn$kcn.binary.tem44d  $OUTDIR/temstate$INYEAR.init
 cp -i $POVDIR/T2M_WL*  $OUTDIR/T2M.$RUNLABEL
 cp -i $POVDIR/TGO3_2D*  $OUTDIR/TGO3_2D.$RUNLABEL
 cp -i $POVDIR/QFLX_2D*  $OUTDIR/QFLX_2D.$RUNLABEL
 cp -i $POVDIR/ocm.climate*  $OUTDIR/ocm.climate.$RUNLABEL

cd $EXEDIR
cat >! ocup.name << EOF
 &CO2UP
    kvh0=$KV,
    kvcmin=1.0
    rkckh = $rkch
    nrun=1,
    areg='equ'
    restart=.false.
    iyeareq=$INYEAR
    co2eq=$CO2EQ
    nyeareq=10000
    co2rfile = '$DATADIR/data/co2_norm_46.new'
    dirdat2  = '$DATADIR/data/46lat/',
    eqclimate= '$OUTDIR/ocm.climate.$RUNLABEL'
    fl_init_dic= '$DATADIR/data/46lat/init_dic'
    fl_init_alkt= '$DATADIR/data/46lat/init_alkt_steph1'
    fl_init_salt= '$DATADIR/data/46lat/init_salt_steph1'
    fl_ed= '$DATADIR/data/46lat/EDZON.dat'
    fl_zoav= '$DATADIR/data/46lat/Z1OAVE46_ZN.dat'
    ocean_eq='$OUTDIR/ocean_eq_$RUNLABEL',
    deepo_eq='$OUTDIR/deepco2.eq_$RUNLABEL',
    ouptake = '$OUTDIR/ouptake_$RUNLABEL',
 &end

EOF


#spin up of ocean carbon model

         cp $MODDIR/build_historical/src/igsm22ocm $OUTDIR/ocm.$RUNLABEL.out

        ./ocm.$RUNLABEL.out  >& $OUTDIR/res$RUNLABEL
#       ./ocm.$RUNLABEL.out 



        endif

cat >! couple.nml << EOF
 &COUPLE_PARM
  dtcouple=1.0
  dtatm=1.0
  dtocn=1.0
  startYear = $INYEAR
  endYear = $LYEAR
 &END

EOF


cat >! name.dat << EOF
EXP. $RUNLABEL $WORDS1
$WORDS2
 &INPUTZ
    AEXP     = $RUNLABEL,
    INYEAR   = $INYEAR,
    LYEAR    = $LYEARCLM,
    DT       = 600.,
    CO2      = -6.0,
    S0RATE   =  -6.0,
    YEARGT   = $YEARGT
    AERFOR   = $FA,
    OBSFOR   = .true.,
    CLDFEED  = .false.,
    CLDFEED  = .true.,
    fclmlice  = '$CLMDATA/clm.pctglacier.txt',
    fbaresoil = '$CLMDATA/clm.pctbaresoil.txt',
    fwmax     = '$CLMDATA/clm.wmax.txt',
    wetfrac_data     = '$CLMDATA/WETLANDS.TEM44b.dat'
    fprratio  = '$CLMDATA/pcp2lnd-data.001x046.txt',
    o3datadir = '$TEMDATA/O3/',
    CLIMO3    = .true.
    clmsen   = $CS,
    TRANSR   = $TRR,
    cfdif0   =  $KV,
    KVCT    = .true.
    diffcar0 = 1.00,
    cfocdif = $rkch
    ISTRT1   = $IRST,
    plotfl   = '$OUTDIR/plot$RUNLABEL',
    nwrfl    = '$OUTDIR/nwr$RUNLABEL',
    file1    = '$OUTDIR/pov1$RUNLABEL',
    file2    = '$OUTDIR/pov2$RUNLABEL',
    pov_deepo= '$OUTDIR/pov_deepo$RUNLABEL',
    deepco2in= '$OUTDIR/deepco2.eq_$RUNLABEL'
    fl_init_alkt= '$DATADIR/data/46lat/init_alkt_steph1'
    fl_init_salt= '$DATADIR/data/46lat/init_salt_steph1'
    caruptfile= '$OUTDIR/carupt.$RUNLABEL',
    fnememiss = '$OUTDIR/nememiss.$RUNLABEL',
    flrco2av= '$OUTDIR/rco2av.$RUNLABEL',
    last_nep = '$OUTDIR/pov_nep$RUNLABEL',
    dirdat1  = '$DATADIR/data/',
    dirdat2  = '$DATADIR/data/46lat/',
    qffile   = '$DATADIR/data/46lat/Z1OAVE46_ZN.dat'
    zmfile   = '$DATADIR/data/46lat/Z1OMAX.dat',
    tsfile   = '$OUTDIR/T2M.$RUNLABEL',
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
 &END

EOF

cat >! name_ocean.dat << EOF

 &INPUTON
    dirdat2   = '$DATADIR/data/46lat/',
    qffile    = '$OUTDIR/QFLX_2D.$RUNLABEL',
    kvnewfl    = '$DATADIR/data/46lat/KVNEWZON.dat',
    zmmaxfl   = 'Z1OMAX46_2D.dat',
    zmavefl   = 'Z1OAVE46_2D.dat',
    povocean  =  '$OUTDIR/povocean$RUNLABEL',
    TRRUN     = $TRR,
    flto3     = '$OUTDIR/TGO3_2D.$RUNLABEL',
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
 finidat        = '$RUNLABEL.clm2.i.nc'
 fsurdat        = ''
 fsurdat        = '$CLMDATA/surfdata_0046x0001_igsmveg_ar5_historical_1860.nc'
 fpftdyn        = ' '
 fatmgrid        = '$CLMDATA/griddata_0046x0001.nc'
 fatmlndfrc        = '$CLMDATA/fracdata_0046x0001.nc'
 fstochdat      = '$CLMDATA/stochastic-data.001x046.GLS.nc'
 fpftcon        = '$CLMDATA/pft-physiology.ar5'
 offline_atmdir = '$CLMDATA/atmdat.nc'
 rpntpath       = '$OUTDIR/pntr$RUNLABEL'
 frivinp_rtm    = '$CLMDATA/rdirc.05.061026'
 fpcp2pft       = '$CLMDATA/pcp2pft-data.001x046.AR5.nc'
 dynamic_pft    = .true.
 dynamic_pft    = .false.
 rampYear_dynpft = $YEARGT
 rampYear_dynpft = 0
 orbitfix       = .false.
 orbitfix       = .true.
 orbityr        =   2000
 nsrest         =  $IRST
 nelapse        =  $NELAPSE
 dtime          =  3600
 start_ymd      =  $STRTYMD
 start_tod      =  43200
 start_tod      =  0
 irad           = -1
 wrtdia         = .true.
 wrtdia         = .false.
 mss_irt        =  0
 hist_dov2xy    = .false.
 hist_dov2xy    = .true.
 hist_nhtfrq    =  -8760
 hist_nhtfrq    =  0
 hist_mfilt     =  1
 hist_ndens     =  2
 hist_crtinic   = 'YEARLY'
 /


EOF


sed -e "s/%kcn/${kcn}/g" < $TEMDATA/tveg43dmod.ecd >! tveg43d.tmp
sed -e "s/%kc/${kc}/g" < tveg43d.tmp >! tveg43d.ecd



cat >! tem.go << EOF
0
1
$OUTDIR/temstate$INYEAR.init
$INYEARM1
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
$OUTDIR/tveg43d.ecd
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
0
1
temstate$LYTEM.kc$kc.kcn$kcn.binary.$RUNLABEL
$LYTEM




EOF







if ($RESTART == "FALSE") then

#echo "MODDIR $MODDIR " >> $outputfile
         cp $MODDIR/build_historical/src/igsm22 mod_$RUNLABEL.out

        (./mod_$RUNLABEL.out < tem.go) >>& $OUTDIR/res$RUNLABEL
#       (./mod_$RUNLABEL.out < tem.go) 

 cp -f $OUTDIR/pov_deepo$RUNLABEL $OUTDIR/init_deepo.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov_nep$RUNLABEL $OUTDIR/init_nep.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov1$RUNLABEL $OUTDIR/pov1.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov2$RUNLABEL $OUTDIR/pov2.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/povocean$RUNLABEL $OUTDIR/povocean.$LYEAR.$RUNLABEL

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
