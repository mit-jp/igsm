!
!       TEM variables
!
!	NOTES: 
!             - Filter out Antarctica for now...
!             - Western boundary -180, Eastern boundary 180
!	      - TEM vector of land points starts at SP, then goes north to NP
!  			and then back down to SP (stepping east)
!	      - Pass water table depth
!	      - Fraction of grid area saturated (water table at sfc.)
!
!	      - Change lsmlat to numgrids (total land-points only)
!	
    integer, parameter :: mxmsaics = 35
    integer, parameter :: mxmthdys = 31
    integer, parameter :: mxdayhrs = 24
    integer, parameter :: max3hrs = 8
    integer, parameter :: mxnlayers = 6
    integer, parameter :: totlayers = 10
    real(r8) :: mitco2(lsmlat)
    real(r8) :: mito3(max3hrs,lsmlat)
    real(r8) :: mittemp(lsmlat)
    real(r8) :: mitdaytemp(mxmthdys,lsmlat)
    real(r8) :: mitswrs(lsmlat)
    real(r8) :: mitpre(lsmlat)
    real(r8) :: mitstrmdur(mxmthdys,lsmlat)
    real(r8) :: mitqstrm(mxmthdys,lsmlat)
    real(r8) :: mitaet(mxmsaics,lsmlat)
    real(r8) :: mitsh2o1m(mxmsaics,lsmlat)
    real(r8) :: mitsh2o2m(mxmsaics,lsmlat)
    real(r8) :: mitswe(mxmsaics,lsmlat)
    real(r8) :: mitsfr(mxmsaics,lsmlat)
    real(r8) :: mitdrn(mxmsaics,lsmlat)
    real(r8) :: mitdaytsoil(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mitdaysh2o(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mithrsh2o(mxdayhrs,mxmthdys,mxmsaics,lsmlat,mxnlayers)
    real(r8) :: wtable(mxmthdys,numgrids)
    real(r8) :: fsat(mxmthdys, numgrids)
    common/climate4tem/mitco2,mito3,mittemp,mitdaytemp,mitswrs,mitpre,mitstrmdur,&
                       mitqstrm,mitaet,mitsh2o1m,mitsh2o2m,mitswe,mitsfr,mitdrn,&
                       mitdaytsoil,mitdaysh2o,wtable,fsat,&
                       mithrsh2o
