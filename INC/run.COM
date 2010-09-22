C $Header: /home/igsm-repository/INC/run.COM,v 1.1 2010/08/31 19:27:01 jscott Exp $
C $Name:  $

      common/rlength/AEXPA,INDAYA,IMONTHA,INYEARA,
     &  LDAYA,LMONTHA,LYEARA,IYEARA,TAUEM, IRESTART
      common/ATMTIME/TAUATM,MONTHATM,JDATEATM,JYEARATM
      character *4 MONTHATM
