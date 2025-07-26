!        NLS VERSION 20 May 1992
!**********************************************************************
!
!                        SUBROUTINE : CGLTRI
!                        -------------------
!
!       This subroutine will generate the "Lanczos" tridiagonal 
!       from some of the quantities which arise naturally in the
!       course of the conjugate gradients algorithm (in the absence
!       of preconditioning).  This is accomplished by multiplying 
!       out the factors of the Cholesky decomposition of the Lanczos 
!       tridiagonal matrix which are implicitly constructed during 
!       the conjugate gradients procedure.  See pp. 370-1 in "Matrix
!       Computations", G. Golub and C. Van Loan, Johns Hopkins Univ.
!       Press, 1983 and pp.  in " ", J. Cullum and R. Willoughby,
!       Springer-Verlag, 1985 for more information.
!
!       written by DJS 20-SEP-87
!
!       Includes:
!               nlsdim.inc
!               rndoff.inc
!
!       Uses:
!
!**********************************************************************
!
      subroutine cgltri(ndone,a,b,shiftr,shifti)
!
      use nlsdim
      use rnddbl
!
      integer ndone
      double precision shiftr,shifti
      double precision a,b
      dimension a(2,MXSTEP),b(2,MXSTEP)
!
      integer i,j
      double precision amp,phase,tr,ti,trhalf,tihalf
      double precision c
      dimension c(2,MXSTEP)
!
!######################################################################
!
      do 10 i=1,ndone
        c(1,i)=a(1,i)
        c(2,i)=a(2,i)
 10   continue
!
      do 20 i=2,ndone
        j=i-1
        tr=b(1,i)
        ti=b(2,i)
!
        amp=sqrt(sqrt(tr*tr+ti*ti))
        if (amp.gt.rndoff) then
          phase=0.5D0*datan2(ti,tr)
          trhalf=amp*cos(phase)
          tihalf=amp*sin(phase)
          if (abs(trhalf).lt.rndoff) trhalf=0.0D0
          if (abs(tihalf).lt.rndoff) tihalf=0.0D0
        else
          trhalf=0.0D0
          tihalf=0.0D0
        end if
!
        b(1,j)=c(1,j)*trhalf-c(2,j)*tihalf
        b(2,j)=c(1,j)*tihalf+c(2,j)*trhalf
        if (b(1,j).lt.0.0D0) then
          b(1,j)=-b(1,j)
          b(2,j)=-b(2,j)
        end if
        a(1,i)=c(1,i)+c(1,j)*tr-c(2,j)*ti
        a(2,i)=c(2,i)+c(1,j)*ti+c(2,j)*tr
 20   continue
!
      do 30 i=1,ndone
        a(1,i)=a(1,i)-shiftr
        a(2,i)=a(2,i)-shifti
 30   continue
!
      b(1,ndone)=0.0D0
      b(2,ndone)=0.0D0
!
!----------------------------------------------------------------------
!     return to caller
!----------------------------------------------------------------------
!
      return 
      end
