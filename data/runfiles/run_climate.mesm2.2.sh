#!/bin/csh 
#
# A shell to generate the name list files:
# couple.nml, name.dat, name_ocean.dat, lnd.stdin
#	and run the model for PI climate run
#
#================================
source /etc/profile.d/modules.csh
 module load pgi/17.3
 module load netcdf

#Kov - diffusion coefficient for ocean temperature
#sen - climate sensitivity
#run number is run.yy

 set Kov = 3.20
 set Kov = 0.00
 set sen = 3.19
 set run = 101
 set yy = 21

 setenv BASEDIR /home/sokolov/MESM2.2/RELEASE/igsm/
 setenv BASEDIR $BASEDIR
# setenv outputfile $BASEDIR/JOBOUTS/output.$run

echo " run_pbs_mesm2.2.sp.rr.sh"


#echo "started $run $Kov $sen " >> $outputfile

# first year of the run
@ INYEAR=1861
@ INYEAR=1
# last year of the run +1
@ LYEAR=150


@ YEARGT=1860

@ INYEARM1 = $INYEAR - 1
@ LYTEM= $LYEAR 
@ LYEARCLM = $LYEAR + 1
@ NYEARS = $LYEARCLM - $INYEAR

setenv MODDIR $BASEDIR

setenv OUTDIR  /scratch/sokolov
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
 setenv WORDS1 "Climate run for $YEARGT  from $INYEAR to $LYTEM "
#setenv WORDS1 "1% per year CO2 run for $YEARGT  from $INYEAR to $LYTEM "
setenv WORDS2 " S=$CS Kvh=$KV "

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






#        for 1860
 setenv POVDIR $DATADIR/init

 cp -i $POVDIR/pov1* $OUTDIR/pov1$RUNLABEL
 cp -i $POVDIR/pov2* $OUTDIR/pov2$RUNLABEL
 cp -i $POVDIR/povocean* $OUTDIR/povocean$RUNLABEL
 cp -i $POVDIR/*clm2.i*  $OUTDIR/$RUNLABEL.clm2.i.nc
 cp -i $POVDIR/T2M*  $OUTDIR/T2M.$RUNLABEL
 cp -i $POVDIR/TGO3_2D*  $OUTDIR/TGO3_2D.$RUNLABEL
 cp -i $POVDIR/QFLX_2D*  $OUTDIR/QFLX_2D.$RUNLABEL


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
    CO2      = 2.0,
    CO2      = -5.0,
    CO2      = 1.0,
    YEARGT   = $YEARGT
    CLDFEED  = .true.,
    fclmlice  = '$CLMDATA/clm.pctglacier.txt',
    fbaresoil = '$CLMDATA/clm.pctbaresoil.txt',
    fwmax     = '$CLMDATA/clm.wmax.txt',
    fprratio  = '$CLMDATA/pcp2lnd-data.001x046.txt',
    clmsen   = $CS,
    TRANSR   = $TRR,
    cfdif0   =  $KV,
    CONTRR   = .false.,
    FORO3   = .true.
    ISTRT1   = $IRST,
    plotfl   = '$OUTDIR/plot$RUNLABEL',
    nwrfl    = '$OUTDIR/nwr$RUNLABEL',
    file1    = '$OUTDIR/pov1$RUNLABEL',
    file2    = '$OUTDIR/pov2$RUNLABEL',
    dirdat1  = '$DATADIR/data/',
    dirdat2  = '$DATADIR/data/46lat/',
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
    aveocean  =  '$OUTDIR/aveocean$RUNLABEL',
    nwravo    =  '$OUTDIR/nwravo$RUNLABEL',
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
 fpftdyn        = ''
 fatmgrid        = '$CLMDATA/griddata_0046x0001.nc'
 fatmlndfrc        = '$CLMDATA/fracdata_0046x0001.nc'
 fstochdat      = '$CLMDATA/stochastic-data.001x046.GLS.nc'
 fpcp2pft       = '$CLMDATA/pcp2pft-data.001x046.AR5.nc'
 fpftcon        = '$CLMDATA/pft-physiology.ar5'
 dynamic_pft    = .false.
 rampYear_dynpft = $YEARGT
 rampYear_dynpft = 0
 orbitfix       = .true.
 orbityr        =   2000
 offline_atmdir = '$CLMDATA/atmdat.nc'
 rpntpath       = '$OUTDIR/pntr$RUNLABEL'
 frivinp_rtm    = '$CLMDATA/rdirc.05.061026'
 nsrest         =  $IRST
 nrevsn         = '$RUNLABEL.clm2.rst'
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
 hist_nhtfrq    =  0
 hist_nhtfrq    =  -8760
 hist_mfilt     =  1
 hist_ndens     =  2
 hist_crtinic   = 'YEARLY'
 /


EOF




cp name_ocean.dat $OUTDIR/name_ocean.$RUNLABEL

cp name.dat $OUTDIR/name.$RUNLABEL




if ($RESTART == "FALSE") then


pwd

          cp $MODDIR/build_climate/src/igsm22 mod_$RUNLABEL.out
#         cp /home/sokolov/MESM2.2/RELEASE/igsm021720/build_climate/src/igsm22 mod_$RUNLABEL.out

 ./mod_$RUNLABEL.out  >>& $OUTDIR/res$RUNLABEL
#./mod_$RUNLABEL.out  

 cp -f $OUTDIR/pov1$RUNLABEL $OUTDIR/pov1.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov2$RUNLABEL $OUTDIR/pov2.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/povocean$RUNLABEL $OUTDIR/povocean.$LYEAR.$RUNLABEL

        cd ../
        mv -i $RUNLABEL $STRDIR
else
       ./mod_$RUNLABEL.out  >>& $OUTDIR/res$RUNLABEL

 cp -f $OUTDIR/pov1$RUNLABEL $OUTDIR/pov1.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/pov2$RUNLABEL $OUTDIR/pov2.$LYEAR.$RUNLABEL
 cp -f $OUTDIR/povocean$RUNLABEL $OUTDIR/povocean.$LYEAR.$RUNLABEL


endif

#echo "finished $run $Kov $sen " >> $outputfile

exit 
