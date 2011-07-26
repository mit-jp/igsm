C $Header: /home/igsm-repository/INC/CLM.h,v 1.1 2010/08/31 19:15:27 jscott Exp $
C $Name:  $

      INTEGER lsmlon, lsmlat, nlevsoi
      PARAMETER(lsmlon=1,lsmlat=N_LAT,nlevsoi=10)
      COMMON/mit2din/ ps4clm, pcpl4clm, pcpc4clm, tpr4clm,
     &  tsl4clm, qs4clm, ws4clm , us4clm, vs4clm, dsw4clm,
     &  dlw4clm, pco24clm, swinr4clm, swvis4clm, dt2mland
      REAL*8 ps4clm (lsmlon,lsmlat)
      REAL*8 pcpl4clm (lsmlon,lsmlat)
      REAL*8 pcpc4clm (lsmlon,lsmlat)
      REAL*8 tpr4clm (lsmlon,lsmlat)
      REAL*8 tsl4clm (lsmlon,lsmlat)
      REAL*8 qs4clm (lsmlon,lsmlat)
      REAL*8 ws4clm (lsmlon,lsmlat)
      REAL*8 us4clm (lsmlon,lsmlat)
      REAL*8 vs4clm (lsmlon,lsmlat)
      REAL*8 dsw4clm (lsmlon,lsmlat)
      REAL*8 dlw4clm (lsmlon,lsmlat)
      REAL*8 pco24clm (lsmlon,lsmlat)
      REAL*8 swinr4clm (lsmlon,lsmlat)
      REAL*8 swvis4clm (lsmlon,lsmlat)
      REAL*8 dt2mland
!
! MIT 2D COMMON BLOCK
!
      COMMON/clm4mit/ lhfclm, shfclm, tauxclm, tauyclm,
     &    asdirclm, aldirclm, asdifclm, aldifclm, sroclm,
     &    ssrclm, glrclm, vetclm, sevclm, cevclm, lwuclm,
     &    tref2mclm, tflxclm, tgndclm, snwdclm, snwcclm,
     &    tsoiclm, h2oiclm, h2olclm
      REAL*8 lhfclm(lsmlon,lsmlat)
      REAL*8 shfclm(lsmlon,lsmlat)
      REAL*8 tauxclm(lsmlon,lsmlat)
      REAL*8 tauyclm(lsmlon,lsmlat)
      REAL*8 asdirclm(lsmlon,lsmlat)
      REAL*8 aldirclm(lsmlon,lsmlat)
      REAL*8 asdifclm(lsmlon,lsmlat)
      REAL*8 aldifclm(lsmlon,lsmlat)
      REAL*8 sroclm(lsmlon,lsmlat)
      REAL*8 ssrclm(lsmlon,lsmlat)
      REAL*8 glrclm(lsmlon,lsmlat)
      REAL*8 vetclm(lsmlon,lsmlat)
      REAL*8 sevclm(lsmlon,lsmlat)
      REAL*8 cevclm(lsmlon,lsmlat)
      REAL*8 lwuclm(lsmlon,lsmlat)
      REAL*8 tref2mclm(lsmlon,lsmlat)
      REAL*8 tflxclm(lsmlon,lsmlat)
      REAL*8 tgndclm(lsmlon,lsmlat)
      REAL*8 snwdclm(lsmlon,lsmlat)
      REAL*8 snwcclm(lsmlon,lsmlat)
      REAL*8 tsoiclm(lsmlon,lsmlat,nlevsoi)
      REAL*8 h2oiclm(lsmlon,lsmlat,nlevsoi)
      REAL*8 h2olclm(lsmlon,lsmlat,nlevsoi)

      COMMON/clm4r95/ ncallclm
      INTEGER  ncallclm

      COMMON/clmprtr/ PRTREND
      LOGICAL  PRTREND
 
!     common/prratio/prlnd2total(lsmlat,12)
