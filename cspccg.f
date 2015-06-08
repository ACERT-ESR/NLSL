c     Version 1.5  5/2/94
c**********************************************************************
c
c                         SUBROUTINE : CSPCCG
c                         -------------------
c
c               (C)omplex (S)ymmetric (P)re(C)onditioned 
c                   (C)onjugate (G)radient Algorithm 
c
c       This subroutine attempts a solution of an Ax=b type problem 
c       where A is a known sparse complex symmetric matrix, b is a 
c       known vector and the approximate solution vector is obtained 
c       via a simple diagonal preconditioned complex symmetric 
c       conjugate gradient procedure.  The algorithm used here is 
c       based in algorithm 10.3-3 in "Matrix Computations", first 
c       edition, by G. Golub and C. Van Loan, Johns Hopkins Univ.
c       Press, 1983.  The real, diagonal part of the matrix is used
c       as a preconditioner here.  The real part of the shift 
c       argument should be chosen to ensure that the preconditioner 
c       is positive definite.  The "Lanczos" tridiagonal matrix is 
c       not formed since the matrix is preconditioned.
c
c       arguments :
c       ---------
c               b      : right hand vector
c               ndim   : number of rows and columns in the matrix 
c               mxcgs  : maximum number of conjugate gradient steps 
c                        allowed
c               cgtol  : maximum residual allowed for termination
c               shift  : a complex shift value used to avoid extraneous
c               invrs  : Inverse of pre-conditioning matrix, taken
c                        to be diagonal
c                        divisions by zero 
c               x      : approximate solution vector 
c               ndone  : number of cg steps executed (positive if 
c                        converged, negative otherwise)
c               error  : pseudonorm of the residual vector at return
c               pcflag : input as .true. if preconditioning is to be performed
c
c       The matrix and index arrays are passed to this routine through 
c       the inclusion of the file eprmat.inc.
c
c       written by DJS 20-SEP-87
c       Modified by SHL and DEB 5-MAY-92 
c         1) Inverse of diagonal part of Z matrix now passed
c            to routine to avoid re-calculating it each iteration
c         2) Preconditioning using this inverse is now an option 
c            specified by setting the pcflag argument to .true.
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
      subroutine cspccg(b,ndim,mxcgs,cgtol,shift,invrs,x,ndone,error,
     #                  pcflag)
c
      implicit none
c
      include 'stddim.inc'
      include 'rndoff.inc'
      include 'eprmat.inc'
      include 'cgdata.inc'
c
      integer ndim,mxcgs,ndone
      double precision cgtol,error,invrs(MXDIM)
      logical pcflag
c
      complex*16 shift
      complex*16 x,b
      dimension x(MXDIM),b(MXDIM)
c
      integer i,nstep
      complex*16 rho1,rho2,alpha,beta,CZERO
      parameter (CZERO=(0.0D0,0.0D0))
c
      complex*16 zdotu
      external zdotu
c
c######################################################################
c
c----------------------------------------------------------------------
c     setup r and rho
c----------------------------------------------------------------------
c
      do i=1,ndim
        r(i)=b(i)
      end do
c
      alpha=CZERO
      beta=CZERO
c
      rho1=zdotu(r,r,ndim)
      rho2=rho1
c
c======================================================================
c     begin loop over CG steps
c======================================================================
c
      nstep=0
 100  nstep=nstep+1
c
c----------------------------------------------------------------------
c     check for convergence
c----------------------------------------------------------------------
c
      error=sqrt(abs(zdotu(r,r,ndim)))
c
      if ((error.lt.cgtol).and.(nstep.gt.1)) then
        nstep=nstep-1
        go to 200
      end if
c
c----------------------------------------------------------------------
c     solve diagonal system of equations
c----------------------------------------------------------------------
c
      if (pcflag) then
         do i=1,ndim
           z(i)=r(i)*invrs(i)
         end do
      else
        do i=1,ndim
          z(i)=r(i)
        end do
      end if 
c
c----------------------------------------------------------------------
c     calculate beta 
c----------------------------------------------------------------------
c
      if (nstep.eq.1) then
        beta=czero
      else
        rho2=rho1
        rho1=zdotu(r,z,ndim)
        beta=rho1/rho2
      end if
c
c----------------------------------------------------------------------
c     update p  ( p <- r+beta*p )
c----------------------------------------------------------------------
c
      if (nstep.eq.1) then
        do i=1,ndim
          p(i)=z(i)
        end do
      else
        call zaypx(z,p,beta,ndim)
      end if
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
c     return error code if not converged
c----------------------------------------------------------------------
c
 200  continue
c
      if (error.le.cgtol) then
        ndone=nstep
      else
        ndone=-nstep
      end if
c
      return
      end
