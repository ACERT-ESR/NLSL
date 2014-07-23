c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                        function VCHANGE
c                    =========================
c
c check if the difference between two FP variables is significant
c .true. means they are different.  Must be double precision v1,v2.
c----------------------------------------------------------------------
      function vchange(v1,v2)
c
      implicit none
      logical vchange
      double precision v1,v2
      include 'rndoff.inc'
c
      vchange=.false.
      if (abs((v1+v2)/2.0D0) .gt. rndoff) then
       vchange=(abs(v1-v2)/(abs(v1+v2)/2.0D0) .gt. 0.000001*abs(v1+v2))
      end if
      return
      end
