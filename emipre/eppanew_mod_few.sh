#! /bin/csh

setenv INDIR  puts/
setenv OUTDIR  ../datanew/emimed/


 foreach  scenario ( sky2050 )

echo $INDIR/$scenario

 cp $INDIR/chm_$scenario.put eppa5chm.put



#a.out
eppaemission.out

  rm -f eppa5chm.put

mv edailytotals $OUTDIR/edaily.$scenario.fordata
mv edaily.dat $OUTDIR/edaily.$scenario.dat
mv edaily.so2 $OUTDIR/edailyso2.$scenario

 end

exit
