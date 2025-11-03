! NLSL Version 1.5.1 beta 11/25/95
!**********************************************************************
!
!                         SUBROUTINE : CSCG
!                         -----------------
!
!       (C)omplex (S)ymmetric (C)onjugate (G)radient Algorithm 
!
!       This subroutine attempts a solution of an Ax=b type problem 
!       where A is a known sparse complex symmetric matrix, b is a 
!       known vector and the approximate solution vector is obtained 
!       via a simple complex symmetric conjugate gradient procedure.  
!       The algorithm used here is based on algorithm 10.3-1 in 
!       "Matrix Computations", first edition, by G. Golub and C. 
!       Van Loan, Johns Hopkins Univ. Press, 1983.  The "Lanczos" 
!       tridiagonal matrix in constructed using formula 10.2-14 in 
!       the reference above.
!
!       It is assumed in this routine that the right hand vector is 
!       normalized and the solution vector is initialized at the start
!       of the routine.  The matrix and index arrays are passed to 
!       this routine through the inclusion of the file eprmat.inc.
!
!       arguments 
!       ---------
!
!             b     : right hand vector
!             ndim  : number of rows and columns in the matrix 
!             mxcgs : maximum number of conjugate gradient steps 
!                     allowed
!             cgtol : maximum residual allowed for termination
!             shift : a complex shift value used to avoid extraneous
!                     divisions by zero 
!             x     : approximate solution vector 
!             al,bl : quantities that can be used to construct the
!                     "Lanczos" tridiagonal matrix
!             ndone : number of cg steps executed (positive if 
!                     converged, negative otherwise)
!             error : pseudonorm of the residual vector at return
!
!             a     : diagonal elements of tridiagonal matrix
!             b     : off-diagonal elements of tridiagonal matrix
!             z     : elements of input matrix
!             izmat : column number and indicator for ends of rows and
!                     real and imaginary parts of matrix elements of
!                     stored in z
!                         if izmat(i)<0 then z(i) is pure real
!                         if izmat(i)=0 then end of row
!                         if izmat(i)>0 then z(i) is pure imaginary
!             zdiag : the diagonal elements of the matrix
!
!
!       Includes:
!               nlsdim.inc
!               rndoff.inc
!               eprmat.inc
!
!       Uses:
!               zdotu.f
!               zaypx.f
!               scmvm.f
!               zaxpy.f
!                               
!**********************************************************************
!
      subroutine cscg(b,ndim,mxcgs,cgtol,shift,
     #                x,al,bl,ndone,error)
!
      use nlsdim
      use rnddbl
      use eprmat
      use stdio
!
      integer ndim,mxcgs,ndone
      double precision cgtol,error
      complex*16 shift
      complex*16 x,b,al,bl
      dimension x(MXDIM),b(MXDIM),al(MXSTEP),bl(MXSTEP)
      integer i,nstep
!
      complex*16 rho1,rho2,CZERO
      parameter (CZERO=(0.0D0,0.0D0))
!
      complex*16 p,r,w,z,alpha,beta
      common /cgdata/ p(MXDIM),r(MXDIM),w(MXDIM),z(MXDIM)
!
      complex*16 zdotu
      external zdotu
!
!######################################################################
!
!.......... Debugging purposes only ..........
!      integer j
!      open (unit=20,file='nlcg.tst',status='unknown',
!     #     access='sequential',form='formatted')
!      write (20,7002)
!      do i=1,ndim
!         if (abs(b(i)).gt. RNDOFF) write (20,7003) i,real(b(i))
!      end do
!
!      write (20,7001) 'Imaginary'
!      do i=1,ndim
!         write (20,7000) i,i,zdiag(2,i)
!         write (20,7000) (i,izmat(j),zmat(j),j=jzmat(i),jzmat(i+1)-1)
!      end do
!
!      write (20,7001) 'Real'
!      do i=1,ndim
!         write (20,7000) i,i,zdiag(1,i)
!         do j=kzmat(i),kzmat(i+1)-1
!            k=MXEL-j+1
!            write (20,7000) i,izmat(k),zmat(k)
!         end do
!      end do
!      close(20)
! 7000 format(i5,',',i5,2x,g14.7)
! 7001 format('c --- ',a,' matrix elements')
! 7002 format('c --- vector elements')
! 7003 format(i5,2x,g14.7)
!........................................
!
!----------------------------------------------------------------------
!               setup r,x and p
!----------------------------------------------------------------------
!
      do i=1,ndim
        r(i)=b(i)
      end do
!
      do i=1,ndim
        p(i)=b(i)
      end do
!
      do i=1,ndim
        x(i)=CZERO
        end do
!
      rho1=zdotu(r,r,ndim)
      rho2=CZERO
!
!======================================================================
!     begin loop over CG steps
!======================================================================
!
      nstep=0
 100  nstep=nstep+1
!
!----------------------------------------------------------------------
!     calculate beta 
!----------------------------------------------------------------------
!
      if (nstep.eq.1) then
        beta=CZERO
      else
        rho2=rho1
        rho1=zdotu(r,r,ndim)
        beta=rho1/rho2
      end if
!
      bl(nstep)=beta
!
!----------------------------------------------------------------------
!     check for convergence
!----------------------------------------------------------------------
!
      error=sqrt(abs(rho1))
!
      if (error.lt.cgtol) go to 200
!
!----------------------------------------------------------------------
!     update p  ( p <- r+beta*p )
!----------------------------------------------------------------------
!
      if (nstep.ne.1) call zaypx(r,p,beta,ndim)
!
!----------------------------------------------------------------------
!     calculate w ( w <- A*p )                
!----------------------------------------------------------------------
!
!
!     Check for window input and/or user halt before proceeding
!
      call wpoll
      if (hltcmd.ne.0 .or. hltfit.ne.0) return

      call scmvm(p,w,ndim)
      call zaxpy(p,w,shift,ndim)
!
!----------------------------------------------------------------------
!     calculate alpha 
!----------------------------------------------------------------------
!
      alpha=rho1/zdotu(p,w,ndim)
!
      al(nstep)=1.0D0/alpha
!
!----------------------------------------------------------------------
!     update x ( x <- x+alpha*p )
!----------------------------------------------------------------------
!
      call zaxpy(p,x,alpha,ndim)
!
!----------------------------------------------------------------------
!     update r ( r <- r-alpha*w )
!----------------------------------------------------------------------
!
      call zaxpy(w,r,-alpha,ndim)
!
!======================================================================
!     end of loop over CG steps
!======================================================================
!
      if (nstep.lt.mxcgs) go to 100
!
!----------------------------------------------------------------------
!     construct Lanczos tridiagonal matrix from CG quantities
!----------------------------------------------------------------------
!
 200  continue
!
      call cgltri(nstep,al,bl,dreal(shift),dimag(shift))
!
!----------------------------------------------------------------------
!     return error code if not converged
!----------------------------------------------------------------------
!
      if (error.le.cgtol) then
        ndone=nstep
      else
        ndone=-nstep
      end if
!
!----------------------------------------------------------------------
!     return to calling routine
!----------------------------------------------------------------------
!
      return
      end
