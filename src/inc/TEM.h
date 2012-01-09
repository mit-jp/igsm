c   These two common block are use to pass data from climate model to
C        TEM and back
      INTEGER no3, ndperm, ncoh, nhor, mnlay, mnlay1
      PARAMETER (no3=8,ndperm=31,ncoh=35,nhor=24,mnlay=6,mnlay1=10)
#if ( defined CLM35 )
      common/wetfrac/wetfrac_data
      character *256 wetfrac_data
#endif

      COMMON/climate4tem/ co24tem, o34tem, temp4tem, dtem4tem, sws4tem,
     &   pre4tem, strmdur, qstrm, 
#if ( defined CLM35 )
     &   pet, 
#endif

     &   aet, sh2o1m, sh2o2m, swe, 
     &   sfr, drn,
     &   daytsoil, daysho, rsh2o
!    &           ,istorms(ndperm,jm0),
!    &            istrmdry(ndperm,jm0)
      REAL*8 co24tem(jm0)
      REAL*8 o34tem(no3,jm0)
      REAL*8 temp4tem(jm0)
      REAL*8 dtem4tem(ndperm,jm0)
      REAL*8 sws4tem(jm0)
      REAL*8 pre4tem(jm0)
      REAL*8 strmdur(ndperm,jm0)
      REAL*8 qstrm(ndperm,jm0)
#if ( defined CLM35 )
      REAL*8 pet(ncoh,jm0)
#endif
      REAL*8 aet(ncoh,jm0)
      REAL*8 sh2o1m(ncoh,jm0)
      REAL*8 sh2o2m(ncoh,jm0)
      REAL*8 swe(ncoh,jm0)
      REAL*8 sfr(ncoh,jm0)
      REAL*8 drn(ncoh,jm0)
      REAL*8 daytsoil(ndperm,ncoh,jm0,mnlay1)
      REAL*8 daysho(ndperm,ncoh,jm0,mnlay1)
      REAL*8 rsh2o(nhor,ndperm,ncoh,jm0,mnlay)

#if ( defined CLM35 )
      COMMON/landcover4tem/ cellarea, vegfrac, landfrac
     &  ,cohortarea
      REAL*8 cellarea(jm0)
      REAL*8 landfrac(jm0)
      REAL*8 vegfrac(ncoh,jm0)
      REAL*8 croptotal(jm0)
      Integer*8 cohortarea(ncoh,jm0),landarea(jm0)
#endif


      COMMON/upt4chem/ temco2, temch4, temn2o
      REAL*8 temco2(jm0)
      REAL*8 temch4(jm0)
      REAL*8 temn2o(jm0)

      COMMON/temtake/ adupt, antemnep, temnep
      REAL*8 adupt
      REAL*8 antemnep(jm0)
      REAL*8 temnep(12,jm0)

C      real  nepan,nepav,n2oann

      COMMON/ozon4tem/ obso3, o3datadir, CLIMO3
      REAL*8 obso3(no3,jm0,12)
      CHARACTER*120 o3datadir
      LOGICAL CLIMO3
C
