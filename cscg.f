c Version 1.6 10/20/94
c**********************************************************************
c
c                         SUBROUTINE : CSCG
c                         -----------------
c
c       (C)omplex (S)ymmetric (C)onjugate (G)radient Algorithm 
c
c       This subroutine attempts a solution of an Ax=b type problem 
c       where A is a known sparse complex symmetric matrix, b is a 
c       known vector and the approximate solution vector is obtained 
c       via a simple complex symmetric conjugate gradient procedure.  
c       The algorithm used here is based on algorithm 10.3-1 in 
c       "Matrix Computations", first edition, by G. Golub and C. 
c       Van Loan, Johns Hopkins Univ. Press, 1983.  The "Lanczos" 
c       tridiagonal matrix in constructed using formula 10.2-14 in 
c       the reference above.
c
c       It is assumed in this routine that the right hand vector is 
c       normalized and the solution vector is initialized at the start
c       of the routine.  The matrix and index arrays are passed to 
c       this routine through the inclusion of the file eprmat.inc.
c
c       arguments 
c       ---------
c
c             b     : right hand vector
c             ndim  : number of rows and columns in the matrix 
c             mxcgs : maximum number of conjugate gradient steps 
c                     allowed
c             cgtol : maximum residual allowed for termination
c             shift : a complex shift value used to avoid extraneous
c                     divisions by zero 
c             x     : approximate solution vector 
c             al,bl : quantities that can be used to construct the
c                     "Lanczos" tridiagonal matrix
c             ndone : number of cg steps executed (positive if 
c                     converged, negative otherwise)
c             error : pseudonorm of the residual vector at return
c
c             a     : diagonal elements of tridiagonal matrix
c             b     : off-diagonal elements of tridiagonal matrix
c             z     : elements of input matrix
c             izmat : column number and indicator for ends of rows and
c                     real and imaginary parts of matrix elements of
c                     stored in z
c                         if izmat(i)<0 then z(i) is pure real
c                         if izmat(i)=0 then end of row
c                         if izmat(i)>0 then z(i) is pure imaginary
c             zdiag : the diagonal elements of the matrix
c
c
c       written by DJS 20-SEP-87
c
c       Includes:
c               stddim.inc
c               rndoff.inc
c               eprmat.inc
c
c       Uses:
c               zdotu.f
c               zaypx.f
c               scmvm.f
c               zaxpy.f
c                               
c**********************************************************************
c
      subroutine cscg(b,ndim,mxcgs,cgtol,shift,x,al,bl,ndone,error)
c
      include 'stddim.inc'
      include 'rndoff.inc'
      include 'eprmat.inc'
      include 'cgdata.inc'
c
      integer ndim,mxcgs,ndone
      double precision cgtol,error
      complex*16 shift
      complex*16 x,b,al,bl
      dimension x(MXDIM),b(MXDIM),al(MXSTEP),bl(MXSTEP)
c
      integer i,nstep
      double precision ZERO
      complex*16 alpha,beta,rho1,rho2,CZERO,CONE
      parameter (ZERO=0.0D0,CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
c
      complex*16 zdotu
      external zdotu
c
c######################################################################
c
c
c----------------------------------------------------------------------
c    Set up r,x and p
c----------------------------------------------------------------------
c
      do i=1,ndim
        r(i)=b(i)
      end do
c
      do i=1,ndim
        p(i)=b(i)
      end do
c
      do i=1,ndim
        x(i)=czero
      end do
c
      rho1=zdotu(r,r,ndim)
      rho2=czero
c
c======================================================================
c     begin loop over CG steps
c======================================================================
c
      nstep=0
 100  nstep=nstep+1
c
c----------------------------------------------------------------------
c     calculate beta 
c----------------------------------------------------------------------
c
      if (nstep.eq.1) then
        beta=czero
      else
        rho2=rho1
        rho1=zdotu(r,r,ndim)
        beta=rho1/rho2
      end if
c
      bl(nstep)=beta
c
c----------------------------------------------------------------------
c     check for convergence
c----------------------------------------------------------------------
c
      error=sqrt(abs(rho1))
c
      if (error.lt.cgtol) go to 200
c
c----------------------------------------------------------------------
c     update p  ( p <- r+beta*p )
c----------------------------------------------------------------------
c
      if (nstep.ne.1) call zaypx(r,p,beta,ndim)
c
c----------------------------------------------------------------------
c     calculate w ( w <- A*p )                
c----------------------------------------------------------------------
c
      call scmvm(p,w,ndim)
      call zaxpy(p,w,shift,ndim)
c
c----------------------------------------------------------------------
c     calculate alpha 
c----------------------------------------------------------------------
c
      alpha=rho1/zdotu(p,w,ndim)
c
      al(nstep)=1.0D0/alpha
c
c----------------------------------------------------------------------
c     update x ( x <- x+alpha*p )
c----------------------------------------------------------------------
c
      call zaxpy(p,x,alpha,ndim)
c
c----------------------------------------------------------------------
c     update r ( r <- r-alpha*w )
c----------------------------------------------------------------------
c
      call zaxpy(w,r,-alpha,ndim)
c
c======================================================================
c     end of loop over CG steps
c======================================================================
c
      if (nstep.lt.mxcgs) go to 100
c
c----------------------------------------------------------------------
c     construct Lanczos tridiagonal matrix from CG quantities
c----------------------------------------------------------------------
c
 200  continue
c
      call cgltri(nstep,al,bl,dreal(shift),dimag(shift))
c
c----------------------------------------------------------------------
c     return error code if not converged
c----------------------------------------------------------------------
c
      if (error.le.cgtol) then
        ndone=nstep
      else
        ndone=-nstep
      end if
c
      return
      end
