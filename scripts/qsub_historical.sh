#!/bin/csh



setenv BASEDIR /home/sokolov/igsm/



set run = 621

# LEV ocean ALL median 20Y
 set Kv = 3.20
 set CS = 3.19
 set FA = -0.24
 set rkch = 0.94
 set KC =  400
 set KCN =  150



set rr = 1

 while ($rr <= 1 ) 




#New Forcing
     sbatch --export=run=$run,Kov=$Kv,sen=$CS,FA=$FA,kc=$KC,kcn=$KCN,rkch=$rkch, -J 20C.$run -n 1 -p edr --time=24:00:00 $BASEDIR/runfiles/run_historical.mesm2.2.sh >>& ../pbsouts/pbsout.$run.13

#    sbatch --export=run0=$run0,run=$run,Kov=$Kv,sen=$CS,FA=$FA,kc=$KC,kcn=$KCN,rkch=$rkch, -J 20C.$run -n 1 -p edr --time=24:00:00 $BASEDIR/RUNFILES/SVL/run_sp.mesm2.2.rstr.sh >>& ../PBSOUTS/SVL/pbsout.$run.13


   @ run = $run + 1
   @ rr = $rr + 1

end

exit
