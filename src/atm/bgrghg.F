
#include "ctrparam.h"

      SUBROUTINE bgrghg(YEART)
      common /bghgdata/bgrghg_data
      common/ BACKGRGHG/GHGBGR(5)
      character * 256 bgrghg_data
      iyeart=yeart
      nyd=300
      open (unit=561,file=bgrghg_data,
     &  status='OLD',form='formatted')
   50 continue
      read (561,*,end=500),ighgyr,GHGBGR
      if(ighgyr.eq.iyeart) then 
        close (561)
        return
      endif
      go to 50
  500 continue
      print *,'From bgrghg wrong YEART ',YEART
      stop
      return
      end
