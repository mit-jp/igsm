C $Header: /home/igsm-repository/INC/TSRF.COM,v 1.1 2010/08/31 19:27:01 jscott Exp $
C $Name:  $

      common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DT2MGL(JM0)
     *,cfcld(JM0,3)
      common/TSLD/TLANDC(JM0,0:13),TLANDT(JM0),TLANDD(JM0),DT2MLD(JM0)
