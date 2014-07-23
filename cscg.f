c  VERSION 1.0  (NLSPMC version)   2/5/99 
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
c       Modified by SHL :  Lanczos vector is stored in array 'lnzvec'.
c       They are used in the conversion of the tridiagonal eigenvector
c	into the eigenvector in the original ritz basis.  The role of
c	the old cgltri.f routine (calculation of the Lanczos tridiagonal
c	from the quantities which arise naturally in the course of the
c	CG algorithm) is also incorporated in this routine to resolve
c	the sign ambiguity that comes in when one takes sqrt(complex #).
c	This happens twice here.
c	  i)  converting beta(CG) into beta(LA)
c	  ii) normalizing the residual vector
c	As long as the sign is consistent between these, the resulting
c	eigenvectors and eigenvalues are the same, even though the sign
c	of eigenvectors and beta's of the tridiagonal matrix may be
c 	different from that of LA.
c
c       arguments 
c       ---------
c
c             b     : right hand vector
c             ndim  : number of rows and columns in the matrix 
c             mxcgs : maximum number of conjugate gradient steps allowed
c             cgtol : maximum residual allowed for termination
c             shift : a complex shift value used to avoid extraneous
c                     divisions by zero 
c             ter,r,p,w : work spaces
c      returns
c      -------
c             x     : approximate solution vector 
c             lnzvec: Lanczos vectors that defines the transformation
c                     from the original to the tridiagonal matrix
c             al    : diagonal elements of tridiagonal matrix
c             bl    : off-diagonal elements of tridiagonal matrix
c             ter   : b-A*x
c             ndone : number of cg steps executed (positive if 
c                     converged, negative otherwise)
c             terror: true error of the solution vector estimated
c                     Termination of the CG steps is based on this
c                     true error.
c
c       written by DJS 20-SEP-87
c	modified by SHL 22-JUL-92
c
c       Includes:
c               nlsdim.inc
c               parcom.inc
c               stdio.inc
c               rndoff.inc
c       Uses:
c               zdotu2.f
c               zaypx.f
c               scmvm.f
c               zaxpy2.f
c                               
c**********************************************************************
c
      subroutine cscg(b,ndim,mxcgs,cgtol,shift,
     #     x,lnzvec,al,bl,ter,r,p,w,ndone,terror,trmin,ntrmin)
c
      include 'limits.inc'
      include 'parms.inc'
      include 'stdio.inc'
      include 'rndoff.inc'
c
      integer ndim,mxcgs,ndone,ntrmin
      double precision cgtol,error,terror,trmin
      complex*16 alpha,alph1,alph2,beta,x,b,lnzvec,al,bl
      dimension x(mxdim),b(mxdim),lnzvec(mxdim,mxstep),
     #          al(mxstep),bl(mxstep)
c
      integer i,nstpcg
      double precision zero,one,amp,phase,tr,ti,trhalf,tihalf
      complex*16 rho1,rho2,rho1h,rho2h,czero,ci,shift
      parameter (zero=0.0D0,one=1.0D0,czero=(0.0D0,0.0D0),
     #           ci=(0.0D0,1.0D0))
c
      complex*16 ter,p,r,w
      dimension ter(mxdim),p(mxdim),r(mxdim),w(mxdim)
c
      complex*16 zdotu2
      external zdotu2
c
c######################################################################
c
      if (idebug.ne.0) write (ludeb,1010)
c
c----------------------------------------------------------------------
c    Initialize r,x and p
c----------------------------------------------------------------------
c
      do 10 i=1,ndim
        r(i)=b(i)
        p(i)=b(i)
 10   continue
c
      rho1=zdotu2(r,r,ndim)
c
      alph1=one
      rho2=one
      rho2h=one
c
c======================================================================
c     begin loop over CG steps
c======================================================================
c
      nstpcg=0
 100  nstpcg=nstpcg+1
c
c----------------------------------------------------------------------
c     calculate beta 
c----------------------------------------------------------------------
c
      if (nstpcg.eq.1) then
        beta=czero
      else
        rho2=rho1
        rho2h=rho1h
        rho1=zdotu2(r,r,ndim)
        beta=rho1/rho2
c
      end if
c
c----------------------------------------------------------------------
c   calculate the pseudo norm of r (sqrt(rho1))
c----------------------------------------------------------------------
c
      tr=dreal(rho1)
      ti=dimag(rho1)
      amp=sqrt(sqrt(tr*tr+ti*ti))
      if (amp.gt.rndoff) then
        phase=0.5D0*datan2(ti,tr)
        trhalf=amp*cos(phase)
        tihalf=amp*sin(phase)
        if (abs(trhalf).lt.rndoff) trhalf=zero
        if (abs(tihalf).lt.rndoff) tihalf=zero
      else
        trhalf=zero
        tihalf=zero
      end if
      rho1h=trhalf+ci*tihalf
c
c----------------------------------------------------------------------
c     store the residual vectors to get the Lanczos vector
c----------------------------------------------------------------------
c
      do 20 i=1,ndim
 20   lnzvec(i,nstpcg)=r(i)/rho1h
c
c----------------------------------------------------------------------
c     update p  ( p <- r+beta*p )
c----------------------------------------------------------------------
c
      if (nstpcg.ne.1) call zaypx(r,p,beta,ndim)
c
c----------------------------------------------------------------------
c     calculate w ( w <- A*p )                
c----------------------------------------------------------------------
c
      call scmvm(p,w,ndim)
c
c----------------------------------------------------------------------
c     add diagonal shift term to matrix-vector product
c----------------------------------------------------------------------
c
      call zaxpy2(p,w,shift,ndim)
c
c----------------------------------------------------------------------
c     calculate alpha and convert into LA tridiagonal form
c----------------------------------------------------------------------
c
      alph2=alph1
      alph1=rho1/zdotu2(p,w,ndim)
c
      al(nstpcg)=one/alph1+beta/alph2-shift
c
      bl(nstpcg)=-beta*rho2h/(alph2*rho1h)
c
c----------------------------------------------------------------------
c     update x ( x <- x+alpha*p )
c----------------------------------------------------------------------
c
      call zaxpy2(p,x,alph1,ndim)
c
c----------------------------------------------------------------------
c     update r ( r <- r-alpha*w )
c----------------------------------------------------------------------
c
      call zaxpy2(w,r,-alph1,ndim)
c
c----------------------------------------------------------------------
c     calculate error and terror & check for convergence
c----------------------------------------------------------------------
c
      call scmvm(x,ter,ndim)
      call zaxpy2(x,ter,shift,ndim)
c
      terror=zero
      do 30 i=1,ndim
        ter(i)=b(i)-ter(i)
        tr=dreal(ter(i))
        ti=dimag(ter(i))
        terror=terror+tr*tr+ti*ti
 30   continue
c
      terror=dsqrt(terror)
c
      error=dsqrt(abs(rho1))
c
      if (error.lt.1.D-17) go to 200
      if (terror.lt.cgtol) go to 200
c record minimum so far found
      if (terror.lt.trmin)then 
	trmin=terror
	ntrmin=nstpcg
      end if
c
      if (idebug.ne.0) then
         if ( (nstpcg/10)*10 .eq. nstpcg ) 
     #         write (ludeb,1020) nstpcg,error,terror
      end if
c
      if (nstpcg.lt.mxcgs) go to 100
c
c======================================================================
c     end of loop over CG steps
c======================================================================
c
 200  continue	! exit loop
c
      if (idebug.ne.0) write (ludeb,1020) nstpcg,error,terror
c
c   shift the index of beta by 1 to match convention
c
      do 40 i=1,nstpcg-1
 40     bl(i)=bl(i+1)
      bl(nstpcg)=czero
c
c----------------------------------------------------------------------
c     return error code if not converged
c----------------------------------------------------------------------
c
      if (terror.le.cgtol) then
c        ndone=nstpcg-1
        ndone=nstpcg   ! test this version, re neg. real parts.
c this choice provides more stable operaition and maybe no neg reals!
      else
c        ndone=-nstpcg+1
        ndone=-nstpcg
      end if
c
c----------------------------------------------------------------------
c     return to calling routine
c----------------------------------------------------------------------
c
      return
c
c=====================================================================
c     format statements
c=====================================================================
c
 1010 format(' ** CG calculation **'/,3x,'step',5x,'error',3x,
     #       'true error')
 1020 format(2x,i4,3x,g12.5,2x,g12.5)
c
      end
