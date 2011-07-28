#!/bin/sh
# Remove some files left over after a run.

FILES="budget globalmean globalmean2 meta monthly monthly2 photodiag \
	pov_deepo* pov_nep* temstate* \
	cflx.glb ch4.glb n2o.glb npp.glb \
	carupt_chem*.out nememiss_chem*.out nwr*.out plot*.out rco2_chem*.out \
	chemrenew.dat name.tmp mittem4.log"

rm -f $FILES
