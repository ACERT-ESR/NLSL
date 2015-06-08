c      Version 1.5  5/2/94
c**********************************************************************
c
c			===================
c			SUBROUTINE : CSLNZS
c			===================
c
c	This subroutine will perform the desired number of steps of 
c	the complex symmetric Lanczos algorithm on the matrix passed
c	throught the common block defined in the include file 
c	eprmat.inc.  The dimension of the arrays are defined in 
c	the file stddim.inc.  The diagonal and off-diagonal matrix 
c	elements of the Lanczos tridiagoanl matrix are returned in 
c	the arrays a and b, respectively.  The algorithm used here is 
c	based on algorithm 9.1-1 in "Matrix Computations", first 
c	edition, by G. Golub and C. Van Loan, Johns Hopkins Univ. 
c	Press, 1983.
c
c
c	arguments
c	---------
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
c	      zdiag : the diagonal elements of the matrix
c
c             nstep : the number ot Lanczos steps to execute
c	      ndim  : the dimension of the matrix
c
c
c       Includes:
c               stddim.inc
c               eprdat.inc
c
c       Uses:
c               scmvma.f
c               zaxpy.f
c               znormu.f
c               zscsw.f
c
c       Written by DJS 26-AUG-87
c
c**********************************************************************
c
        subroutine cslnzs(x,y,nstep,ndim,a,b)
c
        include 'stddim.inc'
        include 'eprmat.inc'
c
        integer nstep,ndim
        complex*16 x,y,a,b         
        dimension x(mxdim),y(mxdim),a(mxstep),b(mxstep)
c
        integer nl
        complex*16 u
        dimension u(mxdim)
c
        complex*16 zdotu,znormu
        external zdotu,znormu
c
c######################################################################
c
c=======================================================================
c               loop over Lanczos steps
c=======================================================================
c
        do 100 nl=1,nstep
c
c----------------------------------------------------------------------
c             update y by doing matrix-vector multiply then adding
c             old y vector to the result, i.e.,
c                       y <- A*x+y = u+y
c----------------------------------------------------------------------
c
         call scmvm(x,u,ndim)
         call zaxpy(u,y,dcmplx(1.0D0,0.0D0),ndim)
c
c----------------------------------------------------------------------
c             calculate alpha
c----------------------------------------------------------------------
c
         a(nl)=zdotu(x,y,ndim)
c
c----------------------------------------------------------------------
c             update y
c----------------------------------------------------------------------
c
         call zaxpy(x,y,-a(nl),ndim)
c
c----------------------------------------------------------------------
c             calculate beta
c----------------------------------------------------------------------
c
         b(nl)=znormu(y,ndim)
c
c----------------------------------------------------------------------
c             rearrange vectors
c----------------------------------------------------------------------
c
         call zscsw(y,x,b(nl),ndim)
c
c=======================================================================
c             end of loop over Lanczos steps
c=======================================================================
c
100     continue
c
c----------------------------------------------------------------------
c             return to calling program
c----------------------------------------------------------------------
c
        return
        end
