c  NLSL Version 1.5 beta 11/25/95
c----------------------------------------------------------------------
c                    =========================
c                      subroutine LMNLS
c                    =========================
c> @brief (L)evenberg-(M)arquardt (N)onlinear (L)east (S)quares
c> ============================================================
c> This is a modification of the original lmder subroutine from the 
c> MINPACK subroutine library. It uses a Levenberg-Marquardt nonlinear
c> least-squares algorithm modified to carry out a local optimization
c> constrained to lie within a "trust region" defined by a step bound
c> delta using scaling of the variables.
c> @details
c> For a description of the trust region approach for least squares 
c> problems, see J.E. Dennis and R.B. Schnabel, Numerical Methods for
c> Unconstrained Optimization and Nonlinear Equations, Prentice-Hall,
c> Englewood Cliffs, NJ (1983), sections 6.4, 7.1, and 10.2.
c----------------------------------------------------------------------
      subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
     *                 maxfev,maxitr,diag,scale,factor,nprint,
     *                 info,nfev,njev,ipvt,qtf,gnvec,gradf,
     *                 wa1,wa2,wa3,wa4)
c
c ----added for EPR NLS -----
c   (also note that fnorm and iter have been moved to common /iterat/)
      use rnddbl
      use nlsdim
      use parcom
      use iterat
      use stdio
c
      integer m,n,ldfjac,maxfev,maxitr,nprint,info,istep,nfev,njev
      integer ipvt(n)
      double precision ftol,xtol,gtol,factor
      double precision x(n),fvec(m),fjac(ldfjac,njcol),diag(njcol),
     *       scale(njcol),qtf(njcol),gnvec(njcol),gradf(njcol),
     *       wa1(njcol),wa2(njcol),wa3(njcol),wa4(m)
      external fcn
c
c----------------------------------------------------------------------
c
c> @details
c>    The purpose of LMDER is to minimize the sum of the squares of
c>    m nonlinear functions in n variables by a modification of
c>    the Levenberg-Marquardt algorithm. The user must provide a
c>    subroutine which calculates the functions and the Jacobian.
c>
c>    The subroutine statement is
c>
c>      subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
c>                       maxfev,maxitr,diag,scale,factor,nprint,info,
c>                       nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c>@param FCN is the name of the user-supplied subroutine which
c>       calculates the functions and the Jacobian. FCN must
c>       be declared in an external statement in the user
c>       calling program, and should be written as follows:
c>>
c>>      subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
c>>      integer m,n,ldfjac,iflag
c>>      double precision x(n),fvec(m),fjac(ldfjac,n)
c>>      ----------
c>>      If iflag=1 calculate the functions at x and
c>>      return this vector in fvec. Do not alter fjac.
c>>      If iflag=2 calculate the Jacobian at x and
c>>      return this matrix in fjac. Do not alter fvec.
c>>      ----------
c>>      return
c>>      end
c>>
c>       The value of IFLAG should not be changed by FCN unless
c>       the user wants to terminate execution of LMDER.
c>       In this case set iflag to a negative integer.
c>
c>@param M is a positive integer input variable set to the number
c>       of functions.
c>
c>@param N  is a positive integer input variable set to the number
c>        of variables. N must not exceed M.
c>
c>@param X is an array of length N. On input X must contain
c>       an initial estimate of the solution vector. On output X
c>       contains the final estimate of the solution vector.
c>
c>@param FVEC is an output array of length M which contains
c>       the functions evaluated at the output X.
c>
c>@param FJAC is an output M by N array. the upper N by N submatrix
c>       of FJAC contains an upper triangular matrix R with
c>       diagonal elements of nonincreasing magnitude such that
c>>
c>>             T     T           T
c>>            P *(JAC *JAC)*P = R *R,
c>>
c>       where P is a permutation matrix and JAC is the final
c>       calculated Jacobian. column j of P is column IPVT(j)
c>       (see below) of the identity matrix. The lower trapezoidal
c>       part of FJAC contains information generated during
c>       the computation of R.
c>
c>@param LDFJAC is a positive integer input variable not less than M
c>       which specifies the leading dimension of the array FJAC.
c>
c>@param FTOL is a nonnegative input variable. Termination
c>       occurs when both the actual and predicted relative
c>       reductions in the sum of squares are at most FTOL.
c>       Therefore, FTOL measures the relative error desired
c>       in the sum of squares.
c>
c>@param XTOL is a nonnegative input variable. Termination
c>       occurs when the relative error between two consecutive
c>       iterates is at most XTOL. Therefore, XTOL measures the
c>       relative error desired in the approximate solution.
c>
c>@param GTOL is a nonnegative input variable. Termination
c>       occurs when the cosine of the angle between FVEC and
c>       any column of the Jacobian is at most GTOL in absolute
c>       value. therefore, GTOL measures the orthogonality
c>       desired between the function vector and the columns
c>       of the Jacobian.
c>
c>@param MAXFEV is a positive integer input variable. Termination
c>       occurs when the number of calls to FCN with IFLAG=1
c>       has reached MAXFEV.
c>
c>@param SCALE is an array of length N containing multiplicative scale
c>       factors for each of the variables in X. If an element of SCALE
c>       is non-positive, it will be reset internally to unity. 
c>       Positive entries in the SCALE array will be retained as 
c>       user-specified scaling factors for the trust-region search 
c>       of the algorithm. The step for the Ith parameter will be scaled 
c>       using the norm of the Ith column of the Jacobian *divided* by
c>       SCALE(I). This produces larger steps along the ith dimension,
c>       at least initialy. The default value for all parameters is unity 
c>       (i.e., the column norms of the Jacobian will be used).
c>>      NB: This convention differs from the original
c>>      specifications of LMDER in MINPACK.
c>
c>@param FACTOR is a positive input variable used in determining the
c>       initial trust region bound. This bound is set to the product of
c>       FACTOR and the Euclidean norm of DIAG*X if nonzero, or else
c>       to FACTOR itself. In most cases FACTOR should lie in the
c>       interval (.1,100.). 100 is a generally recommended value.
c>
c>@param NPRINT is an integer input variable that enables controlled
c>       printing of iterates if it is positive. In this case,
c>       fcn is called with IFLAG=0 at the beginning of the first
c>       iteration and every NPRINT iterations thereafter and
c>       immediately prior to return, with X, FVEC, and FJAC
c>       available for printing. FVEC and FJAC should not be
c>       altered. If NPRINT is not positive, no special calls
c>       of FCN with IFLAG=0 are made.
c>
c>@param INFO @parblock is an integer output variable. If the user has
c>       terminated execution, INFFO is set to the (negative)
c>       value of IFLAG. See description of FCN. Otherwise,
c>       INFO is set as follows:
c>
c>       INFO=0  Improper input parameters.
c>
c>       INFO=1  Both actual and predicted relative reductions
c>                 in the sum of squares are at most FTOL.
c>
c>       INFO=2  Relative error between two consecutive iterates
c>                 is at most XTOL.
c>
c>       INFO=3  conditions for INFO=1 and INFO=2 both hold.
c>
c>       INFO=4  The cosine of the angle between FVEC and any
c>                 column of the Jacobian is at most GTOL in
c>                 absolute value.
c>
c>       INFO=5  number of calls to FCN with IFLAG=1 has
c>                 reached MAXFEV.
c>
c>       INFO=6  FTOL is too small. No further reduction in
c>                 the sum of squares is possible.
c>
c>       INFO=7  XTOL is too small. No further improvement in
c>                 the approximate solution X is possible.
c>
c>       INFO=8  GTOL is too small. FVEC is orthogonal to the
c>                 columns of the Jacobian to machine precision.
c>       @endparblock
c>@param NFEV is an integer output variable set to the number of
c>       calls to FCN with IFLAG=1.
c>
c>@param NJEV is an integer output variable set to the number of
c>       calls to NJEV with IFLAG=2.
c>
c>@param IPVT is an integer output array of length N. IPVT
c>       defines a permutation matrix p such that JAC*P=Q*R,
c>       where JAC is the final calculated Jacobian, Q is
c>       orthogonal (not stored), and R is upper triangular
c>       with diagonal elements of nonincreasing magnitude.
c>       column j of P is column IPVT(j) of the identity matrix.
c>
c>@param QTF is an output array of length N which contains
c>       the first N elements of the vector (Q transpose)*FVEC.
c>
c>@param WA1, WA2, and WA3 are work arrays of length N.
c>
c>@param WA4 is a work array of length M.
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
cc
      integer i,iflag,j,l
      double precision actred,delta,dirder,epsmch,fnorm1,gnorm,par,
     *                 pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm
      PARAMETER (EPSMCH=RNDOFF)
      double precision enorm
c
      double precision ONE,P1,P5,P25,P75,P0001,ZERO
c
      character*1 trstr
      logical grdclc
c
      integer itrim
      external itrim
c ---------------------------
      data ONE,P1,P5,P25,P75,P0001,ZERO
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c
c######################################################################
c
c     epsmch is the machine precision.
c
c
      info=0
      iflag=0
      nfev=0
      njev=0
c
c----------------------------------------------------------------------
c     Check the input parameters for errors.
c----------------------------------------------------------------------
      if (n.le.0 .or. m.lt.n .or. ldfjac.lt.m
     *    .or. ftol.lt.ZERO .or. xtol.lt.ZERO .or. gtol.lt.ZERO
     *    .or. maxfev.le.0 .or.maxitr.le.0 .or. factor.le.ZERO) 
     *    go to 300
c
c**********************************************************************
c
c----------------------------------------------------------------------
c     Initialize Levenberg-Marquardt parameter and iteration counter
c----------------------------------------------------------------------
      par=ZERO
      iter=1
c
c----------------------------------------------------------------------
c     Evaluate the function at the starting point
c     and calculate its norm.
c----------------------------------------------------------------------
      iflag=1
      call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      nfev=1
      if (iflag.lt.0) go to 300
      fnorm=enorm(m,fvec)
c
c----------------------------------------------------------------------
c ********* Beginning of the outer loop *******************************
c----------------------------------------------------------------------
c
   30 continue
c
c----------------------------------------------------------------------
c        Calculate the Jacobian matrix (iflag=2)
c----------------------------------------------------------------------
         iflag=2
         call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         njev=njev+1
         nfev=nfev+n*njev
         if (iflag.lt.0) go to 300
c
c----------------------------------------------------------------------
c        If requested, call fcn to enable printing of iterates
c----------------------------------------------------------------------
         if (nprint.gt.0) then
           iflag=0
           if (mod(iter-1,nprint).eq.0)
     *        call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
           if (iflag.lt.0) go to 300
         end if
c
c----------------------------------------------------------------------
c        Compute the QR factorization of the Jacobian
c----------------------------------------------------------------------

         call qrfac(m,njcol,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c
c----------------------------------------------------------------------
c        On the first iteration, set each non-positive element of the
c        SCALE scaling array according to the norms of the columns of 
c        the initial Jacobian
c----------------------------------------------------------------------
         if (iter.eq.1) then
             do j=1, n
               if (scale(j).le.ZERO) then
                  diag(j)=wa2(j)
               else
                  diag(j)=wa2(j)/scale(j)
               end if
               if (diag(j).eq.ZERO) diag(j)=ONE
            end do
c
            if (itrace.ne.0) then
               write(itrace,1000) (tag(j)(:itrim(tag(j))),j=1,n)
               write(itrace,1001) (wa2(j),j=1,n)
               write(itrace,1002) (diag(j),j=1,n)
            end if
c
c----------------------------------------------------------------------
c        On the first iteration, calculate the norm of the scaled x
c        and initialize the trust region bound delta
c----------------------------------------------------------------------
           do j=1, n
             wa3(j)=diag(j)*x(j)
          end do
          xnorm=enorm(n,wa3)
          delta=factor*xnorm
          if (delta.eq.ZERO) delta=factor
          if (itrace.ne.0) write(itrace,1003) xnorm,delta,factor
       end if
c
c----------------------------------------------------------------------
c        Form (Q transpose)*fvec and store the first njcol components in
c        QtF.
c----------------------------------------------------------------------
         do i=1, m
            wa4(i)=fvec(i)
        end do
c
        do j=1, njcol
           if (fjac(j,j).ne.ZERO) then
              sum=ZERO
              do i=j, m
                 sum=sum + fjac(i,j)*wa4(i)
              end do
              temp=-sum/fjac(j,j)
              do i=j, m
                wa4(i)=wa4(i) + fjac(i,j)*temp
             end do
          end if
c
          fjac(j,j)=wa1(j)
          qtf(j)=wa4(j)
       end do
c
c----------------------------------------------------------------------
c        Compute the norm of the scaled gradient.
c----------------------------------------------------------------------
         gnorm=ZERO
         if (fnorm.ne.ZERO) then
            do j=1, n
               l=ipvt(j)
               if (wa2(l).ne.ZERO) then
                  sum=ZERO
                  do i=1, j
                     sum=sum + fjac(i,j)*(qtf(i)/fnorm)
                  end do
                  gnorm=dmax1(gnorm,dabs(sum/wa2(l)))
               end if
            end do
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
           do j=1,n
               if (scale(j).gt.ZERO) then
                  temp=wa2(j)/scale(j)
                  diag(j)=dmax1(diag(j),temp)
               end if
            end do
c
c----------------------------------------------------------------------
c  ******** Beginning of the inner loop ******************************
c----------------------------------------------------------------------
        istep=0
        grdclc=.false.
  200   continue
c
c----------------------------------------------------------------------
c           Determine the Levenberg-Marquardt parameter.
c----------------------------------------------------------------------
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4,gnvec,gradf)
c
            grdclc=grdclc.or.(par.ne.ZERO)
c----------------------------------------------------------------------
c           Store the direction p and X + p. Calculate the norm of p.
c----------------------------------------------------------------------
            do j=1, n
               wa1(j)=-wa1(j)
               wa2(j)=x(j) + wa1(j)
               wa3(j)=diag(j)*wa1(j)
            end do
            pnorm=enorm(n,wa3)
c
c----------------------------------------------------------------------
c        On the first iteration, adjust the initial trust region bound
c        to the size of the initial step.
c----------------------------------------------------------------------
            trstr=' '
            if (iter.eq.1) then
               if (delta.gt.pnorm) then
                  trstr='*'
               else
                  trstr='s'
               end if
               delta=dmin1(delta,pnorm)
           end if
c
           if (istep.eq.0 .and. itrace.ne.0) then
              write (itrace,1012) iter,(tag(j),j=1,n)
              write (itrace,1013) fnorm,(x(j),j=1,n)
           end if
c
c----------------------------------------------------------------------
c           Evaluate the function at x + p and calculate its norm.
c----------------------------------------------------------------------
            iflag=1
            call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
            nfev=nfev+1
            if (iflag.lt.0) go to 300
            fnorm1=enorm(m,wa4)
            istep=istep+1
c
c----------------------------------------------------------------------
c           Compute the scaled actual reduction.
c----------------------------------------------------------------------
            actred=-ONE
            if (P1*fnorm1.lt.fnorm) actred=ONE-(fnorm1/fnorm)**2
c
c----------------------------------------------------------------------
c           Compute the scaled predicted reduction and
c           the scaled directional derivative.
c----------------------------------------------------------------------
            do j=1,n
               wa3(j)=ZERO
               l=ipvt(j)
               temp=wa1(l)
               do i=1, j
                  wa3(i)=wa3(i) + fjac(i,j)*temp
               end do
            end do
            temp1=enorm(n,wa3)/fnorm
            temp2=(dsqrt(par)*pnorm)/fnorm
            prered=temp1**2 + temp2**2/P5
            dirder=-(temp1**2 + temp2**2)
c
c----------------------------------------------------------------------
c           Compute the ratio of the actual to the predicted
c           reduction.
c----------------------------------------------------------------------
            ratio=ZERO
            if (prered.ne.ZERO) ratio=actred/prered
c
c----------------------------------------------------------------------
c           Update the step bound.
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c             If actual reduction is too much smaller than the predicted
c             reduction (i.e. actred/prered ratio is too small)
c             the function is not well-approximated by a quadratic
c             equation. Reduce the size of the trust region by a
c             factor of 0.1 to 0.5 and increase the L-M parameter.
c----------------------------------------------------------------------
            if (ratio.le.P25) then
              if (actred.ge.ZERO) temp=P5
              if (actred.lt.ZERO) temp=P5*dirder/(dirder + P5*actred)
              if (P1*fnorm1.ge.fnorm .or. temp.lt.P1) temp=P1
              if (delta.gt.pnorm/P1) then
                 trstr='x'
                 delta=pnorm/P1
              endif
              delta=delta*temp
              par=par/temp
c
            else
c
c----------------------------------------------------------------------
c             If ratio of actual to predicted reduction is close to 1,
c             the quadratic model is a good approximation to the function,
c             and we can try increasing the trust region to twice the
c             last step size in order to check whether a better solution
c             is available. Otherwise, the size of the trust region
c             is left unchanged.
c----------------------------------------------------------------------
              if (par.eq.ZERO .or. ratio.ge.P75) then
                 delta=pnorm/P5
                 par=P5*par
                 temp=ONE/P5
              else
                 temp=ONE
              end if
            end if
c
            if (itrace.ne.0) write (itrace,1014) istep,par,ratio,
     #                       ONE/temp,trstr,fnorm1,(wa2(j)-x(j),j=1,n)
c
c----------------------------------------------------------------------
c           Test for successful iteration.
c----------------------------------------------------------------------
            if (ratio.ge.P0001) then
c
c----------------------------------------------------------------------
c           Successful iteration. Update X, FVEC, and their norms.
c----------------------------------------------------------------------
               do j=1, n
                  x(j)=wa2(j)
                  wa2(j)=diag(j)*x(j)
               end do
c
               do i=1, m
                  fvec(i)=wa4(i)
               end do
               xnorm=enorm(n,wa2)
               fnorm=fnorm1
               iter=iter+1
               if (itrace.ne.0) then
                  write(itrace,1006) (diag(j),j=1,n)
                  if (grdclc) write(itrace,1007) (gradf(j),j=1,n)
                  write(itrace,1008) (gnvec(j),j=1,n)
                  write (itrace,1015) delta
               end if
c
            end if
c----------------------------------------------------------------------
c           Tests for convergence.
c----------------------------------------------------------------------
            info=0
            if (dabs(actred).le.ftol .and. prered.le.ftol
     *          .and. P5*ratio.le.ONE) info=1
            if (delta.le.xtol*xnorm) info=info + 2
            if (info.ne.0) go to 300
c
c----------------------------------------------------------------------
c           Tests for termination and stringent tolerances.
c----------------------------------------------------------------------
            if (nfev.ge.maxfev) info=5
            if (iter.ge.maxitr) info=6
            if (dabs(actred).le.epsmch .and. prered.le.epsmch
     *          .and. P5*ratio.le.ONE) info=7
            if (delta.le.epsmch*xnorm) info=8
            if (gnorm.le.epsmch) info=9
            if (info.ne.0) go to 300
c----------------------------------------------------------------------
c           End of the inner loop. Repeat if iteration unsuccessful.
c----------------------------------------------------------------------
            if (ratio.lt.P0001) go to 200
c----------------------------------------------------------------------
c        End of the outer loop.
c----------------------------------------------------------------------
         go to 30
  300 continue
c----------------------------------------------------------------------
c     Termination, either normal or user-imposed.
c----------------------------------------------------------------------
      if (iflag.lt.0) info=10
      iflag=3
      if (info.gt.4) iflag=-iflag
      if (info.ne.0.and.info.ne.10.and.nprint.gt.0) 
     #     call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      return
c
c ##### format statements for trace printout ########################
c
 1000 format(/10x,41('=')/10x,
     #      'TRACE OF LEVENBERG-MARQUARDT MINIMIZATION'/
     #     10x,41('=')//'INITIAL SCALING:'/13x,10(1x,a9,2x))
 1001 format('Col norms of J:',10(2x,g10.4)/)
 1002 format(9x,'Scale:',10(2x,g10.4))
 1003 format(/10x,'Scaled X norm: ',g11.5/6x,'Trust region (TR):',
     #     g11.5,'  =(Xnorm*',g9.3,')')
 1006 format(79('-')/t26,'Scale:',10(1x,g10.4))
 1007 format(t23,'Gradient:',10(1x,g10.4))
 1008 format(t21,'G-N vector:',10(1x,g10.4))
 1012 format(/'##### Iteration',i3,1x,59('#')/
     #'Stp LMpar Ratio Trscl  Fnorm    ',12(2x,a9))
 1013 format(79('-')/'  0',t22,g11.5,12g11.4)
 1014 format(i3,2f6.2,f5.2,a1,g11.5,sp,12g11.4)
 1015 format(t23,'Final TR:    ',g10.4)
      end
