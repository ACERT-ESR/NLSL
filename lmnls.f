c NLSPMC Version 1.0 2/5/99
c----------------------------------------------------------------------
c                    =========================
c                      subroutine LMNLS
c                    =========================
c     (L)evenberg-(M)arquardt (N)onlinear (L)east (S)quares
c
c This is a modification of the original lmder subroutine from the 
c MINPACK subroutine library. It uses a Levenberg-Marquardt nonlinear
c least-squares algorithm modified to carry out a local optimization
c constrained to lie within a "trust region" defined by a step bound
c delta using scaling of the variables.
c
c For a description of the trust region approach for least squares 
c problems, see J.E. Dennis and R.B. Schnabel, Numerical Methods for
c Unconstrained Optimization and Nonlinear Equations, Prentice-Hall,
c Englewood Cliffs, NJ (1983), sections 6.4, 7.1, and 10.2.
c
c Modified to use limits on the variables as set in the calling 
c program.  x when called is the physical real parameter, convert it
c to an infinite range parameter for real range defined by prmin,
c prmax.  Convert back to physical parameters on call to simulation
c program.  
c 
c----------------------------------------------------------------------
      subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
     *                 maxfev,maxitr,diag,scale,factor,nprint,itrace,
     *                 jacobi,info,nfev,njev,ipvt,qtf,gnvec,gradf,
     *                 wa1,wa2,wa3,wa4)
      integer m,n,ldfjac,maxfev,maxitr,nprint,info,istep,nfev,njev,
     *        itrace,jacobi
      integer ipvt(n)
      double precision ftol,xtol,gtol,factor
      double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),scale(n),
     *                 qtf(n),gnvec(n),gradf(n),wa1(n),wa2(n),wa3(n),
     *                 wa4(m)
      external fcn
c
c----------------------------------------------------------------------
c
c     The purpose of LMDER is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the Levenberg-Marquardt algorithm. The user must provide a
c     subroutine which calculates the functions and the Jacobian.
c
c     The subroutine statement is
c
c       subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
c                        maxfev,maxitr,diag,scale,factor,nprint,info,
c                        nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       FCN is the name of the user-supplied subroutine which
c         calculates the functions and the Jacobian. FCN must
c         be declared in an external statement in the user
c         calling program, and should be written as follows:
c
c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
c         integer m,n,ldfjac,iflag
c         double precision x(n),fvec(m),fjac(ldfjac,n)
c         ----------
c         If iflag=1 calculate the functions at x and
c         return this vector in fvec. Do not alter fjac.
c         If iflag=2 calculate the Jacobian at x and
c         return this matrix in fjac. Do not alter fvec.
c         ----------
c         return
c         end
c
c         The value of IFLAG should not be changed by FCN unless
c         the user wants to terminate execution of LMDER.
c         In this case set iflag to a negative integer.
c
c       M is a positive integer input variable set to the number
c         of functions.
c
c       N  is a positive integer input variable set to the number
c          of variables. N must not exceed M.
c
c       X is an array of length N. On input X must contain
c         an initial estimate of the solution vector. On output X
c         contains the final estimate of the solution vector.
c
c       FVEC is an output array of length M which contains
c         the functions evaluated at the output X.
c
c       FJAC is an output M by N array. the upper N by N submatrix
c         of FJAC contains an upper triangular matrix R with
c         diagonal elements of nonincreasing magnitude such that
c
c                T     T           T
c               P *(JAC *JAC)*P = R *R,
c
c         where P is a permutation matrix and JAC is the final
c         calculated Jacobian. column j of P is column IPVT(j)
c         (see below) of the identity matrix. The lower trapezoidal
c         part of FJAC contains information generated during
c         the computation of R.
c
c       LDFJAC is a positive integer input variable not less than M
c         which specifies the leading dimension of the array FJAC.
c
c       FTOL is a nonnegative input variable. Termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most FTOL.
c         Therefore, FTOL measures the relative error desired
c         in the sum of squares.
c
c       XTOL is a nonnegative input variable. Termination
c         occurs when the relative error between two consecutive
c         iterates is at most XTOL. Therefore, XTOL measures the
c         relative error desired in the approximate solution.
c
c       GTOL is a nonnegative input variable. Termination
c         occurs when the cosine of the angle between FVEC and
c         any column of the Jacobian is at most GTOL in absolute
c         value. therefore, GTOL measures the orthogonality
c         desired between the function vector and the columns
c         of the Jacobian.
c
c       MAXFEV is a positive integer input variable. Termination
c         occurs when the number of calls to FCN with IFLAG=1
c         has reached MAXFEV.
c
c       SCALE is an array of length N containing multiplicative scale
c         factors for each of the variables in X. If an element of SCALE
c         is non-positive, it will be reset internally to unity. 
c         Positive entries in the SCALE array will be retained as 
c         user-specified scaling factors for the trust-region search 
c         of the algorithm. The step for the Ith parameter will be scaled 
c         by SCALE(I) times the norm of the Ith column of the Jacobian. 
c         The default value for all parameters is unity (i.e., the
c         column norms of the Jacobian will be used).
c          
c         NB: This convention differs from the original
c         specifications of LMDER in MINPACK.
c
c       FACTOR is a positive input variable used in determining the
c         initial trust region bound. This bound is set to the product of
c         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else
c         to FACTOR itself. In most cases FACTOR should lie in the
c         interval (.1,100.). 100 is a generally recommended value.
c
c       NPRINT is an integer input variable that enables controlled
c         printing of iterates if it is positive. In this case,
c         fcn is called with IFLAG=0 at the beginning of the first
c         iteration and every NPRINT iterations thereafter and
c         immediately prior to return, with X, FVEC, and FJAC
c         available for printing. FVEC and FJAC should not be
c         altered. If NPRINT is not positive, no special calls
c         of FCN with IFLAG=0 are made.
c
c       INFO is an integer output variable. If the user has
c         terminated execution, INFFO is set to the (negative)
c         value of IFLAG. See description of FCN. Otherwise,
c         INFO is set as follows:
c
c         INFO=0  Improper input parameters.
c
c         INFO=1  Both actual and predicted relative reductions
c                   in the sum of squares are at most FTOL.
c
c         INFO=2  Relative error between two consecutive iterates
c                   is at most XTOL.
c
c         INFO=3  conditions for INFO=1 and INFO=2 both hold.
c
c         INFO=4  The cosine of the angle between FVEC and any
c                   column of the Jacobian is at most GTOL in
c                   absolute value.
c
c         INFO=5  number of calls to FCN with IFLAG=1 has
c                   reached MAXFEV.
c
c         INFO=6  FTOL is too small. No further reduction in
c                   the sum of squares is possible.
c
c         INFO=7  XTOL is too small. No further improvement in
c                   the approximate solution X is possible.
c
c         INFO=8  GTOL is too small. FVEC is orthogonal to the
c                   columns of the Jacobian to machine precision.
c
c       NFEV is an integer output variable set to the number of
c         calls to FCN with IFLAG=1.
c
c       NJEV is an integer output variable set to the number of
c         calls to NJEV with IFLAG=2.
c
c       IPVT is an integer output array of length N. IPVT
c         defines a permutation matrix p such that JAC*P=Q*R,
c         where JAC is the final calculated Jacobian, Q is
c         orthogonal (not stored), and R is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of P is column IPVT(j) of the identity matrix.
c
c       QTF is an output array of length N which contains
c         the first N elements of the vector (Q transpose)*FVEC.
c
c       WA1, WA2, and WA3 are work arrays of length N.
c
c       WA4 is a work array of length M.
c
c     Subprograms called
c
c       User-supplied ...... FCN
c
c       MINPACK-supplied ... DPMPAR,ENORM,LMPAR,QRFAC
c
c       FORTRAN-supplied ... DABS,DMAX1,DMIN1,DSQRT,MOD
c
c     Argonne National Laboratory. MINPACK Project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c----------------------------------------------------------------------
      integer i,iflag,j,l,ld
c
      double precision actred,delta,dirder,epsmch,fnorm1,gnorm,
     *                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
     *                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c
c ----added for EPR NLS
      include 'limits.inc'
      include 'parms.inc'
c      include 'iterat.inc'
c
c     epsmch is the machine precision.
c
      epsmch=dpmpar(1)
c
c on entry, convert physical xxx into infinite range lmnls variables.
c
      call mapxxx(x,n,-1)
      info=0
      iflag=0
      nfev=0
      njev=0
c
c----------------------------------------------------------------------
c     Check the input parameters for errors.
c----------------------------------------------------------------------
      if (n.le.0 .or. m.lt.n .or. ldfjac.lt.m
     *    .or. ftol.lt.zero .or. xtol.lt.zero .or. gtol.lt.zero
     *    .or. maxfev.le.0 .or.maxitr.le.0 .or. factor.le.zero) 
     *    go to 300
c
c**********************************************************************
c
c----------------------------------------------------------------------
c     Evaluate the function at the starting point
c     and calculate its norm.
c----------------------------------------------------------------------
      iflag=1
      call mapxxx(x,n,1)	! map to physical variables
      call fcn(m,n,x,fvec,fjac,ldfjac,iflag)	! pfun
      if (ihltcmd.ne.0) return
      call mapxxx(x,n,-1)	! map to lmnls variables
      nfev=1
c      seteval=.true.	! got eigenvalues. - not used anymore.
      if (iflag.lt.0) go to 300
      fnorm=enorm(m,fvec)
c
c----------------------------------------------------------------------
c     Initialize Levenberg-Marquardt parameter and iteration counter
c----------------------------------------------------------------------
      par=zero
      iter=1
c
c----------------------------------------------------------------------
c ********* Beginning of the outer loop *******************************
c----------------------------------------------------------------------
c
   30 continue
c----------------------------------------------------------------------
c        Calculate the Jacobian matrix (iflag=2)
c        and signal fcn to output it if necessary (iflag=3) 
c----------------------------------------------------------------------
         iflag=2
         if (jacobi.ne.0) iflag=3
         call mapxxx(x,n,1)	! map to physical variables
         call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         if (ihltcmd.ne.0) return
         call mapxxx(x,n,-1)	! map to lmnls variables
c
c       call ftest1(fjac(1,1),16384,'in lmnls, fjac(1)',17)
c       call ftest1(fjac(1,2),16384,'in lmnls, fjac(2)',17)
         njev=njev+1
         nfev=nfev+n
         if (iflag.lt.0) go to 300
c----------------------------------------------------------------------
c        If requested, call fcn to enable printing of iterates
c----------------------------------------------------------------------
         if (nprint.gt.0) then
           iflag=0
           if (mod(iter-1,nprint).eq.0)then
              call mapxxx(x,n,1)	! map to physical variables
              call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
              if (ihltcmd.ne.0) return
              call mapxxx(x,n,-1)	! map to lmnls variables
           end if
           if (iflag.lt.0) go to 300
         end if
c----------------------------------------------------------------------
c        Compute the QR factorization of the Jacobian
c----------------------------------------------------------------------
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c ***** wa2 here is large?
c----------------------------------------------------------------------
c        On the first iteration, set each non-positive element of the
c        SCALE scaling array according to the norms of the columns of 
c        the initial Jacobian
c----------------------------------------------------------------------
         if (iter.eq.1) then
             do 50 j=1, n
               if (scale(j).le.zero) then
                  diag(j)=wa2(j)
               else
                  diag(j)=wa2(j)/scale(j)
               end if
               if (diag(j).eq.zero) diag(j)=one
   50          continue
c
            if (itrace.ne.0) then
               write(itrace,1000) (tag(j),j=1,n)
               write(itrace,1001) (wa2(j),j=1,n)
               write(itrace,1002) (diag(j),j=1,n)
            end if
c
c----------------------------------------------------------------------
c        On the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta
c----------------------------------------------------------------------
           do 70 j=1, n
             wa3(j)=diag(j)*x(j)
   70        continue
           xnorm=enorm(n,wa3)
           delta=factor*xnorm
           if (delta.eq.zero) delta=factor
           if (itrace.ne.0) write(itrace,1003) xnorm,delta,factor
         end if
c
c----------------------------------------------------------------------
c        Form (Q transpose)*fvec and store the first n components in
c        QtF.
c----------------------------------------------------------------------
         do 90 i=1, m
           wa4(i)=fvec(i)
   90      continue
         do 130 j=1, n
c
           if (fjac(j,j).ne.zero) then
              sum=zero
              do 100 i=j, m
                sum=sum + fjac(i,j)*wa4(i)
  100           continue
              temp=-sum/fjac(j,j)
              do 110 i=j, m
                wa4(i)=wa4(i) + fjac(i,j)*temp
  110           continue
            end if
c
            fjac(j,j)=wa1(j)
            qtf(j)=wa4(j)
  130       continue
c
c----------------------------------------------------------------------
c        Compute the norm of the scaled gradient.
c----------------------------------------------------------------------
         gnorm=zero
         if (fnorm.ne.zero) then
           do 160 j=1, n
             l=ipvt(j)
             if (wa2(l).ne.zero) then
               sum=zero
               do 140 i=1, j
                 sum=sum + fjac(i,j)*(qtf(i)/fnorm)
  140            continue
              gnorm=dmax1(gnorm,dabs(sum/wa2(l)))
            end if
  160       continue
          end if
c
c----------------------------------------------------------------------
c        Test for convergence of the gradient norm
c----------------------------------------------------------------------
         if (gnorm.le.gtol) info=4
         if (info.ne.0) go to 300
c
c----------------------------------------------------------------------
c        Rescale diag array
c----------------------------------------------------------------------
           do 180 j=1,n
               if (scale(j).le.zero) then
                  temp=wa2(j)
               else
                  temp=wa2(j)/scale(j)
               end if
                 diag(j)=dmax1(diag(j),temp)
  180        continue
c
c
c----------------------------------------------------------------------
c  ******** Beginning of the inner loop ******************************
c----------------------------------------------------------------------
	istep=1
  200   continue
c----------------------------------------------------------------------
c           Determine the Levenberg-Marquardt parameter.
c----------------------------------------------------------------------
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4,gnvec,gradf)
c----------------------------------------------------------------------
c           Store the direction p and X + p. Calculate the norm of p.
c----------------------------------------------------------------------
c
            do 210 j=1, n
               wa1(j)=-wa1(j)
               wa2(j)=x(j) + wa1(j)
               wa3(j)=diag(j)*wa1(j)
 210        continue
            pnorm=enorm(n,wa3)
c
            if (itrace.ne.0 .and. istep.eq.1) then
               write(itrace,1004) iter,fnorm,(tag(j),j=1,n)
               write(itrace,1005) (x(j),j=1,n)
               write(itrace,1006) (diag(j),j=1,n)
               if (par.ne.zero) write(itrace,1007) (gradf(j),j=1,n)
               write(itrace,1008) (-gnvec(j),j=1,n)
            end if
c
c----------------------------------------------------------------------
c        On the first iteration, adjust the initial trust region bound
c        to the size of the initial step.
c----------------------------------------------------------------------
            if (iter.eq.1) then
              if (itrace.ne.0) then
                 write(itrace,1010) pnorm
                 if (delta.gt.pnorm) write (itrace,1011) 
              end if
              delta=dmin1(delta,pnorm)
            end if
c
c----------------------------------------------------------------------
c           Evaluate the function at x + p and calculate its norm.
c----------------------------------------------------------------------
            iflag=1
            call mapxxx(wa2,n,1)	! map to physical variables
            call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
            if (ihltcmd.ne.0) return
            call mapxxx(wa2,n,-1)	! map to lmnls variables
            nfev=nfev+1
            if (iflag.lt.0) go to 300
            fnorm1=enorm(m,wa4)
c
            if (itrace.ne.0) write (itrace,1012) istep,par,
     #                                           delta,fnorm1*fnorm1
            istep=istep+1
c
c----------------------------------------------------------------------
c           Compute the scaled actual reduction.
c----------------------------------------------------------------------
            actred=-one
            if (p1*fnorm1.lt.fnorm) actred=one - (fnorm1/fnorm)**2
c
c----------------------------------------------------------------------
c           Compute the scaled predicted reduction and
c           the scaled directional derivative.
c----------------------------------------------------------------------
            do 230 j=1,n
               wa3(j)=zero
               l=ipvt(j)
               temp=wa1(l)
               do 220 i=1, j
                  wa3(i)=wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1=enorm(n,wa3)/fnorm
            temp2=(dsqrt(par)*pnorm)/fnorm
            prered=temp1**2 + temp2**2/p5
            dirder=-(temp1**2 + temp2**2)
c
c----------------------------------------------------------------------
c           Compute the ratio of the actual to the predicted
c           reduction.
c----------------------------------------------------------------------
            ratio=zero
            if (prered.ne.zero) ratio=actred/prered
c
c----------------------------------------------------------------------
c           Update the step bound.
c----------------------------------------------------------------------
            if (ratio.le.p25) then
c
c----------------------------------------------------------------------
c             If actual reduction is too much smaller than the predicted
c             reduction (i.e. actred/prered ratio is too small)
c             the function is not well-approximated by a quadratic
c             equation. Reduce the size of the trust region by a
c             factor of 0.1 to 0.5 and increase the L-M parameter.
c----------------------------------------------------------------------
c
              if (actred.ge.zero) temp=p5
              if (actred.lt.zero) temp=p5*dirder/(dirder + p5*actred)
              if (p1*fnorm1.ge.fnorm .or. temp.lt.p1) temp=p1
	      if (p1*delta.gt.pnorm) then
                 if (itrace.ne.0) write(itrace,1013) one/p1,pnorm/p1
                 delta=pnorm/p1
              endif
              delta=delta*temp
              par=par/temp
c
            else
c
c----------------------------------------------------------------------
c             If ratio of actual to predicted reduction is close to 1,
c             the quadratic model is a good approximation to the function,
c             and we can try increasing the trust region by a factor of
c             two to see if a better solution is available. Otherwise, 
c             the size of the step bound is left unchanged.
c----------------------------------------------------------------------
c
              if (par.eq.zero .or. ratio.ge.p75) then
                delta=pnorm/p5
                par=p5*par
                temp=one/p5
              else
                 temp=one
              end if
            end if
c
            if (itrace.ne.0) write (itrace,1014) ratio,temp
            if (itrace.ne.0) write (itrace,1015) (wa2(j)-x(j),j=1,n)
c
c----------------------------------------------------------------------
c           Test for successful iteration.
c----------------------------------------------------------------------
            if (ratio.ge.p0001) then
c----------------------------------------------------------------------
c           Successful iteration. Update X, FVEC, and their norms.
c----------------------------------------------------------------------
              do 270 j=1, n
                x(j)=wa2(j)
                wa2(j)=diag(j)*x(j)
  270           continue
              do 280 i=1, m
                fvec(i)=wa4(i)
  280           continue
              xnorm=enorm(n,wa2)
              fnorm=fnorm1
              iter=iter+1
            end if
c----------------------------------------------------------------------
c           Tests for convergence.
c----------------------------------------------------------------------
            info=0
            if (dabs(actred).le.ftol .and. prered.le.ftol
     *          .and. p5*ratio.le.one) info=1
            if (delta.le.xtol*xnorm) info=info + 2
            if (info.ne.0) go to 300
c----------------------------------------------------------------------
c           Tests for termination and stringent tolerances.
c----------------------------------------------------------------------
            if (nfev.ge.maxfev) info=5
            if (iter.ge.maxitr) info=6
            if (dabs(actred).le.epsmch .and. prered.le.epsmch
     *          .and. p5*ratio.le.one) info=7
            if (delta.le.epsmch*xnorm) info=8
            if (gnorm.le.epsmch) info=9
            if (info.ne.0) go to 300
c----------------------------------------------------------------------
c           End of the inner loop. Repeat if iteration unsuccessful.
c----------------------------------------------------------------------
            if (ratio.lt.p0001) go to 200
c----------------------------------------------------------------------
c        End of the outer loop.
c----------------------------------------------------------------------
         go to 30
  300 continue
c----------------------------------------------------------------------
c     Termination, either normal or user-imposed.
c----------------------------------------------------------------------
      if (iflag.lt.0) then
         info=10
         nprint=0
      end if
      iflag=0
      call mapxxx(x,n,1)	! map to physical variables
      if (nprint.gt.0) call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      return
c
c ##### format statements for trace printout ########################
c
 
 1000 format(65('#')/10x,'TRACE OF LEVENBERG-MARQUARDT MINIMIZATION',
     #     /65('#')//
     #       'INITIAL VALUES:'/16x,10(3x,a6,3x))
 1001 format('Col norms of J:',10(2x,g10.4)/)
 1002 format(9x,'Scale:',10(2x,g10.4))
 1003 format(/10x,'Scaled X norm: ',g11.5/5x,'Trust region bound:',
     #     g11.5,'  =(',g9.3,'*Xnorm)')
 1004 format(/65('#')/'Iteration',i3,': Chi-Sq=',g11.5//12x,
     #     10(3x,a6,3x))
 1005 format(7x,'Xvec:',10(2x,g10.4))
 1006 format(6x,'Scale:',10(2x,g10.4))
 1007 format(3x,'Gradient:',10(2x,g10.4))
 1008 format(1x,'G-N vector:',10(2x,g10.4))
 1010 format(/3x,'Initial step size=',g10.4)
 1011 format(3x,'(TR bound has been reduced to the initial ',
     #     'step size)'/)
 1012 format(5x,65('-')/5x,'Step',i3,'; LMpar=',g10.4,' TR bound=',
     #     g10.4,'; Chi-Sq=',g12.5)
 1013 format(5x,'(TR bound reduced to ',f4.1,'*step length =',g11.5,')')
 1014 format(5x,'Actual/predicted Chi-Sq reduction=',g10.4,
     #     '; TR scaled by ',f4.1/5x,65('-'))
 1015 format(12x,10(2x,g10.4))
 2000 format(/6x,'*** Jacobian calculation : ',i3,' ***'/6x,
     #     'trans(J)*fvec',15x,'trans(J)*J'/)
 2001 format(6x,g12.5,5x,6(g12.5,1x))
      end
