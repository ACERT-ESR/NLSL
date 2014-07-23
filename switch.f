c NLSPMC Version 1.0 2/5/99
c*********************************************************************
c
c                         SUBROUTINE SWITCH
c                         =================
c
c        This subroutine shifts the second half of the FT,
c        making it the first half.
c
c*********************************************************************
c
      subroutine switch(temp,npt)
c
      integer i,j,k
c
      complex*16 temp(npt),tmp
c
c#####################################################################
c
      k=npt/2
      do 10 i=1,k
        tmp=temp(i+k)
        temp(i+k)=temp(i)
 10     temp(i)=tmp
c
      return
      end
