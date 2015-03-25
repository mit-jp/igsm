
      SUBROUTINE ORBIT (OBLIQ,ECCN,OMEGT,DAY,SDIST,SIND,COSD,LAMBDA)    8201.
C****                                                                   8202.
C**** ORBIT receives the orbital parameters and time of year, and       8203.
C**** returns the distance from the sun and its declination angle.      8204.
C**** The reference for the following caculations is: V.M.Blanco        8205.
C**** and S.W.McCuskey, 1961, "Basic Physics of the Solar System",      8206.
C**** pages 135 - 151.                                                  8207.
C****                                                                   8208.
C**** Program authors: Gary L. Russell and Robert J. Suozzo, 12/13/85   8209.
C****                                                                   8210.
C****        All computations are in double-precision;                  8211.
C****        but the arguments are single-precision.                    8212.
C**** Input: OBLIQ = latitude of tropics in degrees                     8213.
C****        ECCEN = eccentricity of the orbital ellipse                8214.
C****        OMEGT = angle from vernal equinox to perihelion in degrees 8215.
C****        DAY   = day of the year in days; 0 = Jan 1, hour 0         8216.
C****                                                                   8217.
C**** Constants: EDAYPY = Earth days per year = 365                     8218.
C****            VERQNX = occurence of vernal equinox = day 79 = Mar 21 8219.
C****                                                                   8220.
C**** Intermediate quantities:                                          8221.
C****    PERIHE = perihelion during the year in temporal radians        8222.
C****    MA     = mean anomaly in temporal radians = 2J DAY/365 - PERIHE8223.
C****    EA     = eccentric anomaly in radians                          8224.
C****    TA     = true anomaly in radians                               8225.
C****    BSEMI  = semi minor axis in units of the semi major axis       8226.
C****    GREENW = longitude of Greenwich in the Earth's reference frame 8227.
C****                                                                   8228.
C**** Output: DIST = distance to the sun in units of the semi major axis8229.
C****        SDIST = square of DIST                                     8229.5
C****         SIND = sine of the declination angle                      8230.
C****         COSD = cosine of the declination angle                    8231.
C****       LAMBDA = sun longitude in Earth's rotating reference frame  8232.
C****                                                                   8233.
      IMPLICIT REAL*8 (A-H,O-Z)                                         8234.
      REAL*8 MA                                                         8235.
C     REAL*4 SIND,COSD,SDIST,LAMBDA,OBLIQ,ECCN,OMEGT,DAY                8236.
C****                                                                   8237.
      PI = 3.14159265358979D0                                           8238.
      EDAYPY = 365.                                                     8239.
      VERQNX = 79.                                                      8240.
      OMEGA=OMEGT*(PI/180.D0)                                           8241.
      DOBLIQ=OBLIQ*(PI/180.D0)                                          8242.
      ECCEN=ECCN                                                        8243.
C****                                                                   8244.
C**** Determine time of perihelion using Kepler's equation:             8245.
C**** PERIHE-VERQNX = OMEGA - ECCEN sin(OMEGA)                          8246.
C****                                                                   8247.
      PERIHE = OMEGA-ECCEN*SIN(OMEGA)+VERQNX*2.*PI/365.                 8248.
C     PERIHE = DMOD(PERIHE,2.*PI)                                       8249.
      MA = 2.*PI*DAY/365.-PERIHE                                        8250.
      MA = DMOD(MA,2.*PI)                                               8251.
C****                                                                   8252.
C**** Numerically solve Kepler's equation: MA = EA - ECCEN sin(EA)      8253.
C****                                                                   8254.
      EA = MA+ECCEN*(SIN(MA)+ECCEN*SIN(2.*MA)/2.)                       8255.
  110 DEA = (MA-EA+ECCEN*SIN(MA))/(1.-ECCEN*COS(EA))                    8256.
      EA = EA+DEA                                                       8257.
      IF (DABS(DEA).GT.1.D-8)  GO TO 110                                8258.
C****                                                                   8259.
C**** Calculate the distance to the sun and the true anomaly            8260.
C****                                                                   8261.
      BSEMI = DSQRT(1.-ECCEN*ECCEN)                                     8262.
      COSEA = COS(EA)                                                   8263.
      SINEA = SIN(EA)                                                   8264.
      SDIST  = (1.-ECCEN*COSEA)*(1.-ECCEN*COSEA)                        8265.
      TA = DATAN2(SINEA*BSEMI,COSEA-ECCEN)                              8266.
C****                                                                   8267.
C**** Change the reference frame to be the Earth's equatorial plane     8268.
C**** with the Earth at the center and the positive x axis parallel to  8269.
C**** the ray from the sun to the Earth were it at vernal equinox.      8270.
C**** The distance from the current Earth to that ray (or x axis) is:   8271.
C**** DIST sin(TA+OMEGA).  The sun is located at:                       8272.
C****                                                                   8273.
C**** SUN    = (-DIST cos(TA+OMEGA),                                    8274.
C****           -DIST sin(TA+OMEGA) cos(OBLIQ),                         8275.
C****            DIST sin(TA+OMEGA) sin(OBLIQ))                         8276.
C**** SIND   = sin(TA+OMEGA) sin(OBLIQ)                                 8277.
C**** COSD   = sqrt(1-SIND**2)                                          8278.
C**** LAMBDA = atan[tan(TA+OMEGA) cos(OBLIQ)] - GREENW                  8279.
C**** GREENW = 2*3.14159 DAY (EDAYPY-1)/EDAYPY                          8280.
C****                                                                   8281.
      SINDD = SIN(TA+OMEGA)*SIN(DOBLIQ)                                 8282.
      COSD = DSQRT(1.-SINDD*SINDD)                                      8283.
      SIND = SINDD                                                      8284.
C     GREENW = 2.*PI*(DAY-VERQNX)*(EDAYPY+1.)/EDAYPY                    8285.
C     SUNX = -COS(TA+OMEGA)                                             8286.
C     SUNY = -SIN(TA+OMEGA)*COS(DOBLIQ)                                 8287.
C     LAMBDA = DATAN2(SUNY,SUNX)-GREENW                                 8288.
C     LAMBDA = DMOD(LAMBDA,2.*PI)                                       8289.
C****                                                                   8290.
      RETURN                                                            8291.
      END                                                               8292.
