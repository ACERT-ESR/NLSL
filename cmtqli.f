c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c
c       This subroutine will diagonalize a complex symmetric 
c       tridiagonal matrix using the QL algorithm with implicit
c       shifts.  It also accumulates the first row of the eigenvector
c       matrix which is required for calculating ESR spectra from
c       the eigenvalues of the stochastic Liouville matrix.  The
c       algorithm used here is a slight modification of the algorithm
c       discussed CMTQL1 in "Lanczos Algorithms for Large Symmetric
c       Eigenvalue Computations", vol. 2, J. Cullum and R. Willoughby,
c       Birkhauser, 1985.
c
c       written by DJS 26-NOV-87
c
c       Includes:
c               nlsdim.inc
c               rndoff.inc
c
c       Uses:
c
c**********************************************************************
c
      subroutine cmtqli(n,d,e,z,isort,ierr)
c
      include 'limits.inc'
      include 'rndoff.inc'
c
      integer n,ierr
      complex*16 d,e,z
      dimension d(mxstep),e(mxstep),z(mxstep)
c
      integer i,j,l,m,mml,isort
      double precision temp,t0,t1,eps
      complex*16 b,c,f,g,p,r,s,w
c
      complex*16 czero,cone
      parameter (czero=(0.0D0,0.0D0),cone=(1.0D0,0.0D0))
c
c######################################################################
c
      eps=100.0D0*rndoff
c
c----------------------------------------------------------------------
c       intialize first row of eigenvector matrix
c----------------------------------------------------------------------
c
      z(1)=cone
      do 10 i=2,n
        z(i)=czero
 10   continue
c
      ierr=0
      if (n.eq.1) go to 180
c
c======================================================================
c               loop over eigenvalues
c======================================================================
c
      e(n)=czero
c
      do 140 l=1,n
c
        j=0
c
c----------------------------------------------------------------------
c               find a small subdiagonal matrix element
c----------------------------------------------------------------------
c
 20     continue
c
        do 30 m=l,n-1
          temp=cdabs(d(m))+cdabs(d(m+1))
          if (cdabs(e(m)).le.temp*rndoff) go to 40
 30     continue
c
        m=n
c
 40     p=d(l)
        if (m.ne.l) then
          if (j.eq.100) go to 170
          j=j+1
c
c----------------------------------------------------------------------
c               form shift 
c----------------------------------------------------------------------
c
          g=(d(l+1)-p)*0.5D0
          t0=cdabs(g)
          t1=cdabs(e(l))
c
          if (t0.le.t1) then
            w=g/e(l)
            r=cdsqrt(cone+w*w)
            t0=cdabs(w+r)
            t1=cdabs(w-r)
            if (t1.le.t0) then
              g=d(m)-p+e(l)/(w+r)
            else
              g=d(m)-p+e(l)/(w-r)
            end if
          else
            w=e(l)/g
            r=cdsqrt(cone+w*w)
            t0=cdabs(cone+r)
            t1=cdabs(cone-r)
            if (t1.le.t0) then
              g=d(m)-p+w*e(l)/(cone+r)

            else
              g=d(m)-p+w*e(l)/(cone-r)
            end if
          end if        
c
          s=cone
          c=-cone
          p=czero
          mml=m-l
c
          do 90 i=m-1,l,-1
            f=s*e(i)
            b=-c*e(i)
            t0=cdabs(g)
            t1=cdabs(f)
            if (t1.le.t0) then
              w=f/g
              r=cdsqrt(cone+w*w)
              e(i+1)=g*r
              c=cone/r
              s=w*c
            else
              w=g/f
              r=cdsqrt(cone+w*w)
              e(i+1)=f*r
              s=cone/r
              c=w*s
            end if
            temp=1.0D0+cdabs(w)**2
            t0=dsqrt(temp)
            t1=cdabs(r)
c
            if (t1.le.eps*t0) then
              ierr=-l
              go to 180
            else
              ierr=0
            end if
c
c----------------------------------------------------------------------
c               finish up loop over i
c----------------------------------------------------------------------
c
            g=d(i+1)-p
            r=(d(i)-g)*s+2.0D0*c*b
            p=s*r
            d(i+1)=g+p
            g=b-c*r
c
c----------------------------------------------------------------------
c               keep track of effect of rotations on first row
c----------------------------------------------------------------------
c
            w=z(i+1)
            z(i+1)=s*z(i)+c*w
            if (cdabs(z(i+1)).lt.rndoff) z(i+1)=czero
            z(i)=c*z(i)-s*w
            if (cdabs(z(i)).lt.rndoff) z(i)=czero
c
c----------------------------------------------------------------------
c               end of loop to restore to tridiagonal form
c----------------------------------------------------------------------
c
 90       continue
c
c----------------------------------------------------------------------
c               go back to find another small off-diagonal element
c----------------------------------------------------------------------
c
          d(l)=d(l)-p
          e(l)=g
          e(m)=czero
c
          go to 20
        end if
c
c----------------------------------------------------------------------
c               end of loop over eigenvalues
c----------------------------------------------------------------------
c
 140  continue
c
      do 145 l=1,n
 145	z(l)=z(l)*z(l)
c
c----------------------------------------------------------------------
c      sort in the order of increasing magnitude of the eigenvalues
c      or decreasing order of weighting factor
c----------------------------------------------------------------------
c
c                           *** sort by eigenvalue ***
      if (isort.eq.0) then
        do 150 l=1,n-1
          t0=cdabs(d(l))
          do 150 m=l+1,n
            t1=cdabs(d(m))
            if (t1.lt.t0) then
              t0=t1
              w=d(m)
              d(m)=d(l)
              d(l)=w
              w=z(m)
              z(m)=z(l)
              z(l)=w
            end if
  150   continue
      else
c                           *** sort by weighting factor ***
        do 160 l=1,n-1
          t0=cdabs(z(l))
          do 160 m=l+1,n
            t1=cdabs(z(m))
            if (t1.gt.t0) then
              t0=t1
              w=d(m)
              d(m)=d(l)
              d(l)=w
              w=z(m)
              z(m)=z(l)
              z(l)=w
            end if
  160   continue
      end if     
c
      go to 180
c
 170  ierr=l
c
c----------------------------------------------------------------------
c       return to caller
c----------------------------------------------------------------------
c
 180  return
      end
