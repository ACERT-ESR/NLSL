c     Version 1.5  5/2/94 
c**********************************************************************
c
c                        SUBROUTINE : CGLTRI
c                        -------------------
c
c       This subroutine will generate the "Lanczos" tridiagonal 
c       from some of the quantities which arise naturally in the
c       course of the conjugate gradients algorithm (in the absence
c       of preconditioning).  This is accomplished by multiplying 
c       out the factors of the Cholesky decomposition of the Lanczos 
c       tridiagonal matrix which are implicitly constructed during 
c       the conjugate gradients procedure.  See pp. 370-1 in "Matrix
c       Computations", G. Golub and C. Van Loan, Johns Hopkins Univ.
c       Press, 1983 and pp.  in " ", J. Cullum and R. Willoughby,
c       Springer-Verlag, 1985 for more information.
c
c       written by DJS 20-SEP-87
c
c       Includes:
c               stddim.inc
c               rndoff.inc
c
c       Uses:
c
c**********************************************************************
c
      subroutine cgltri(ndone,a,b,shiftr,shifti)
c
      include 'stddim.inc'
      include 'rndoff.inc'
c
      integer ndone
      double precision shiftr,shifti
      double precision a,b
      dimension a(2,MXSTEP),b(2,MXSTEP)
c
      integer i,j
      double precision amp,phase,tr,ti,trhalf,tihalf
      double precision c
      dimension c(2,MXSTEP)
c
c######################################################################
c
      do i=1,ndone
        c(1,i)=a(1,i)
        c(2,i)=a(2,i)
      end do
c
      do i=2,ndone
        j=i-1
        tr=b(1,i)
        ti=b(2,i)
c
        amp=sqrt(sqrt(tr*tr+ti*ti))
        if (amp.gt.RNDOFF) then
          phase=0.5D0*datan2(ti,tr)
          trhalf=amp*cos(phase)
          tihalf=amp*sin(phase)
          if (abs(trhalf).lt.RNDOFF) trhalf=0.0D0
          if (abs(tihalf).lt.RNDOFF) tihalf=0.0D0
        else
          trhalf=0.0D0
          tihalf=0.0D0
        end if
c
        b(1,j)=c(1,j)*trhalf-c(2,j)*tihalf
        b(2,j)=c(1,j)*tihalf+c(2,j)*trhalf
        if (b(1,j).lt.0.0D0) then
          b(1,j)=-b(1,j)
          b(2,j)=-b(2,j)
        end if
        a(1,i)=c(1,i)+c(1,j)*tr-c(2,j)*ti
        a(2,i)=c(2,i)+c(1,j)*ti+c(2,j)*tr
      end do
c
      do i=1,ndone
        a(1,i)=a(1,i)-shiftr
        a(2,i)=a(2,i)-shifti
      end do
c
      b(1,ndone)=0.0D0
      b(2,ndone)=0.0D0
c
      return 
      end
