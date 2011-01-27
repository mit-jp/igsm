
C  common data for ML Ocean model

      DATA JDPERY/365/,JMPERY/12/,EDPERY/365./
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,    
     *  ALAMI/2.1762/,TFO/-1.56/
      DATA Z1I/.1/,Z2OIM/0.4/,Z2OIX/6.9/,Z2OIMAX /5.0/
      DATA TTRUNC/0./
