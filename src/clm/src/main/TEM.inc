!
!       CLM/IGSM to TEM/NEM variables
!
!	NOTES: 
!             - Filter out Antarctica for now...
!             - Western boundary -180, Eastern boundary 180
!	      - TEM vector of land points starts at SP, then goes north to NP
!  			and then back down to SP (stepping east)
!	
    integer, parameter :: mxmsaics = 35
    integer, parameter :: mxmthdys = 31
    integer, parameter :: mxdayhrs = 24
    integer, parameter :: max3hrs = 8
    integer, parameter :: mxnlayers = 6
    integer, parameter :: totlayers = 10
    integer, parameter :: igsmlat = 46
    real*8 :: mitco2(igsmlat)
    real*8 :: mito3(max3hrs,igsmlat)
    real*8 :: mittemp(igsmlat)
    real*8 :: mitdaytemp(mxmthdys,igsmlat)
    real*8 :: mitswrs(igsmlat)
    real*8 :: mitpre(igsmlat)
    real*8 :: mitstrmdur(mxmthdys,igsmlat)
    real*8 :: mitqstrm(mxmthdys,igsmlat)
    real*8 :: mitpet(mxmsaics,igsmlat)
    real*8 :: mitaet(mxmsaics,igsmlat)
    real*8 :: mitsh2o1m(mxmsaics,igsmlat)
    real*8 :: mitsh2o2m(mxmsaics,igsmlat)
    real*8 :: mitswe(mxmsaics,igsmlat)
    real*8 :: mitsfr(mxmsaics,igsmlat)
    real*8 :: mitdrn(mxmsaics,igsmlat)
    real*8 :: mitdaytsoil(mxmthdys,mxmsaics,igsmlat,totlayers)
    real*8 :: mitdaysh2o(mxmthdys,mxmsaics,igsmlat,totlayers)
    real*8 :: mithrsh2o(mxdayhrs,mxmthdys,mxmsaics,igsmlat,mxnlayers)
    common/climate4tem/mitco2,mito3,mittemp,mitdaytemp,mitswrs,mitpre,mitstrmdur,&
                       mitqstrm,mitpet,mitaet,mitsh2o1m,mitsh2o2m,mitswe,mitsfr,mitdrn,&
                       mitdaytsoil,mitdaysh2o,&
                       mithrsh2o

    common/vegcoverclm/ vegfrac
    real*8 :: vegfrac(mxmsaics,igsmlat)
    common/landcover4tem/ cellarea4tem, cohortarea4tem, landfrac4tem, Incohortarea, &
                          pctclay4tem, pctsand4tem, pctsilt4tem, &
                          soilph4tem, ksat4tem, por4tem, zthick4tem  
    real*8 :: cellarea4tem(igsmlat)
    integer*8 :: cohortarea4tem(mxmsaics,igsmlat)
    integer*8 :: landfrac4tem(igsmlat)
    integer*8 :: Incohortarea(mxmsaics,igsmlat)
    real*8 :: pctclay4tem(igsmlat)
    real*8 :: pctsand4tem(igsmlat)
    real*8 :: pctsilt4tem(igsmlat)
    real*8 :: soilph4tem(igsmlat)
    real*8 :: ksat4tem(totlayers,igsmlat)
    real*8 :: por4tem(totlayers,igsmlat)
    real*8 :: zthick4tem(totlayers,igsmlat)
    real*8 :: croptotal(igsmlat)

