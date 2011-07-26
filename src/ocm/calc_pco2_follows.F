c       ==========================================================
c
c       CALC_PCO2_FOLLOWS:  Subroutine for calculating pH
c                following Follows, Ito and Dutkiewicz (2006)
c                Ocean Modelling, 12, 290-301.
c            -------------------------------------------------
c

c
      subroutine calc_pco2_follows(
     I                       diclocal,pt,sit,ta,
     I                       k1local,k2local,
     I                       k1plocal,k2plocal,k3plocal,
     I                       kblocal,kwlocal,
     I                       ksilocal,fflocal,btlocal,
     U                       pHlocal,co2slocal,co3local)
c     ==========================
       implicit none
C     == Routine arguments ==
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
        real  pt, sit, ta
        real  diclocal 
        real  fflocal, btlocal
        real  k1local, k2local
        real  k1plocal, k2plocal, k3plocal
        real  kblocal, kwlocal, ksilocal
        real  pHlocal,co2slocal,co3local

C     == Local variables ==
        real  phguess
        real  cag
        real  bohg
        real  hguess
        real  stuff
        real  gamm
        real  hnew
        real  co2s
        real  h3po4g, h2po4g, hpo4g, po4g   
        real  siooh3g
        real  permil


c ---------------------------------------------------------------------
C Change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
c where the ocean's mean surface density is 1024.5 kg/m^3
c Note: mol/kg are actually what the body of this routine uses
c for calculations.  Units are reconverted back to mol/m^3 at the
c end of this routine.
c To convert input in mol/m^3 -> mol/kg
        permil=1.0 / 1024.5
        pt=pt*permil
        sit=sit*permil
        ta=ta*permil
        diclocal=diclocal*permil
c ---------------------------------------------------------------------
c set first guess for [H+] solvers
c first guess (for newton-raphson)
        phguess = phlocal
cmick - new approx method
cmick - make estimate of htotal (hydrogen ion conc) using
cmick   appromate estimate of CA, carbonate alkalinity
        hguess = 10.0**(-phguess)
cmick - first estimate borate contribution using guess for [H+]
        bohg = btlocal*kblocal/(hguess+kblocal)
cmick - first estimate of contribution from phosphate
cmick based on Dickson and Goyet
        stuff = hguess*hguess*hguess
     &           + (k1plocal*hguess*hguess)
     &           + (k1plocal*k2plocal*hguess)
     &           + (k1plocal*k2plocal*k3plocal)
        h3po4g = (pt*hguess*hguess*hguess) / stuff
        h2po4g = (pt*k1plocal*hguess*hguess) / stuff
        hpo4g  = (pt*k1plocal*k2plocal*hguess) / stuff
        po4g   = (pt*k1plocal*k2plocal*k3plocal) / stuff

cmick - estimate contribution from silicate
cmick based on Dickson and Goyet
        siooh3g = sit*ksilocal / (ksilocal + hguess)

cmick - now estimate carbonate alkalinity
        cag = ta - bohg - (kwlocal/hguess) + hguess
     &           - hpo4g - 2.0*po4g + h3po4g
     &           - siooh3g

cmick - now evaluate better guess of hydrogen ion conc
cmick   htotal = [H+], hydrogen ion conc
        gamm  = diclocal/cag
        stuff = (1.0-gamm)*(1.0-gamm)*k1local*k1local
     &          - 4.0*k1local*k2local*(1.0-2.0*gamm)
        hnew  = 0.5*( (gamm-1.0)*k1local + sqrt(stuff) )
cmick - now determine [CO2*]
        co2slocal  = diclocal/
     &   (1.0 + (k1local/hnew) + (k1local*k2local/(hnew*hnew)))
cmick - return update pH to main routine
        phlocal = -log10(hnew)

c NOW EVALUATE CO32-, carbonate ion concentration
c used in determination of calcite compensation depth
c Karsten Friis & Mick - Sep 2004
        co3local = k1local*k2local*diclocal / 
     &         (hnew*hnew + k1local*hnew + k1local*k2local)

C ----------------------------------------------------------------
c Reconvert from mol/kg -> mol/m^3
        pt=pt/permil
        sit=sit/permil
        ta=ta/permil
        diclocal=diclocal/permil
        co2slocal=co2slocal/permil
        co3local=co3local/permil
c
        return
        end
