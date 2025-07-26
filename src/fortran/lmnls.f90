!  NLSL Version 1.5 beta 11/25/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine LMNLS
!                    =========================
!> @brief (L)evenberg-(M)arquardt (N)onlinear (L)east (S)quares
!> ============================================================
!> This is a modification of the original lmder subroutine from the 
!> MINPACK subroutine library. It uses a Levenberg-Marquardt nonlinear
!> least-squares algorithm modified to carry out a local optimization
!> constrained to lie within a "trust region" defined by a step bound
!> delta using scaling of the variables.
!> @details
!> For a description of the trust region approach for least squares 
!> problems, see J.E. Dennis and R.B. Schnabel, Numerical Methods for
!> Unconstrained Optimization and Nonlinear Equations, Prentice-Hall,
!> Englewood Cliffs, NJ (1983), sections 6.4, 7.1, and 10.2.
!----------------------------------------------------------------------
      subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
     *                 maxfev,maxitr,diag,scale,factor,nprint,
     *                 info,nfev,njev,ipvt,qtf,gnvec,gradf,
     *                 wa1,wa2,wa3,wa4)
!
! ----added for EPR NLS -----
!   (also note that fnorm and iter have been moved to common /iterat/)
      use rnddbl
      use nlsdim
      use parcom
      use iterat
      use stdio
!
      integer m,n,ldfjac,maxfev,maxitr,nprint,info,istep,nfev,njev
      integer ipvt(n)
      double precision ftol,xtol,gtol,factor
      double precision x(n),fvec(m),fjac(ldfjac,njcol),diag(njcol),
     *       scale(njcol),qtf(njcol),gnvec(njcol),gradf(njcol),
     *       wa1(njcol),wa2(njcol),wa3(njcol),wa4(m)
      external fcn
!
!----------------------------------------------------------------------
!
!> @details
!>    The purpose of LMDER is to minimize the sum of the squares of
!>    m nonlinear functions in n variables by a modification of
!>    the Levenberg-Marquardt algorithm. The user must provide a
!>    subroutine which calculates the functions and the Jacobian.
!>
!>    The subroutine statement is
!>
!>      subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!>                       maxfev,maxitr,diag,scale,factor,nprint,info,
!>                       nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!>@param FCN is the name of the user-supplied subroutine which
!>       calculates the functions and the Jacobian. FCN must
!>       be declared in an external statement in the user
!>       calling program, and should be written as follows:
!>>
!>>      subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!>>      integer m,n,ldfjac,iflag
!>>      double precision x(n),fvec(m),fjac(ldfjac,n)
!>>      ----------
!>>      If iflag=1 calculate the functions at x and
!>>      return this vector in fvec. Do not alter fjac.
!>>      If iflag=2 calculate the Jacobian at x and
!>>      return this matrix in fjac. Do not alter fvec.
!>>      ----------
!>>      return
!>>      end
!>>
!>       The value of IFLAG should not be changed by FCN unless
!>       the user wants to terminate execution of LMDER.
!>       In this case set iflag to a negative integer.
!>
!>@param M is a positive integer input variable set to the number
!>       of functions.
!>
!>@param N  is a positive integer input variable set to the number
!>        of variables. N must not exceed M.
!>
!>@param X is an array of length N. On input X must contain
!>       an initial estimate of the solution vector. On output X
!>       contains the final estimate of the solution vector.
!>
!>@param FVEC is an output array of length M which contains
!>       the functions evaluated at the output X.
!>
!>@param FJAC is an output M by N array. the upper N by N submatrix
!>       of FJAC contains an upper triangular matrix R with
!>       diagonal elements of nonincreasing magnitude such that
!>>
!>>             T     T           T
!>>            P *(JAC *JAC)*P = R *R,
!>>
!>       where P is a permutation matrix and JAC is the final
!>       calculated Jacobian. column j of P is column IPVT(j)
!>       (see below) of the identity matrix. The lower trapezoidal
!>       part of FJAC contains information generated during
!>       the computation of R.
!>
!>@param LDFJAC is a positive integer input variable not less than M
!>       which specifies the leading dimension of the array FJAC.
!>
!>@param FTOL is a nonnegative input variable. Termination
!>       occurs when both the actual and predicted relative
!>       reductions in the sum of squares are at most FTOL.
!>       Therefore, FTOL measures the relative error desired
!>       in the sum of squares.
!>
!>@param XTOL is a nonnegative input variable. Termination
!>       occurs when the relative error between two consecutive
!>       iterates is at most XTOL. Therefore, XTOL measures the
!>       relative error desired in the approximate solution.
!>
!>@param GTOL is a nonnegative input variable. Termination
!>       occurs when the cosine of the angle between FVEC and
!>       any column of the Jacobian is at most GTOL in absolute
!>       value. therefore, GTOL measures the orthogonality
!>       desired between the function vector and the columns
!>       of the Jacobian.
!>
!>@param MAXFEV is a positive integer input variable. Termination
!>       occurs when the number of calls to FCN with IFLAG=1
!>       has reached MAXFEV.
!>
!>@param SCALE is an array of length N containing multiplicative scale
!>       factors for each of the variables in X. If an element of SCALE
!>       is non-positive, it will be reset internally to unity. 
!>       Positive entries in the SCALE array will be retained as 
!>       user-specified scaling factors for the trust-region search 
!>       of the algorithm. The step for the Ith parameter will be scaled 
!>       using the norm of the Ith column of the Jacobian *divided* by
!>       SCALE(I). This produces larger steps along the ith dimension,
!>       at least initialy. The default value for all parameters is unity 
!>       (i.e., the column norms of the Jacobian will be used).
!>>      NB: This convention differs from the original
!>>      specifications of LMDER in MINPACK.
!>
!>@param FACTOR is a positive input variable used in determining the
!>       initial trust region bound. This bound is set to the product of
!>       FACTOR and the Euclidean norm of DIAG*X if nonzero, or else
!>       to FACTOR itself. In most cases FACTOR should lie in the
!>       interval (.1,100.). 100 is a generally recommended value.
!>
!>@param NPRINT is an integer input variable that enables controlled
!>       printing of iterates if it is positive. In this case,
!>       fcn is called with IFLAG=0 at the beginning of the first
!>       iteration and every NPRINT iterations thereafter and
!>       immediately prior to return, with X, FVEC, and FJAC
!>       available for printing. FVEC and FJAC should not be
!>       altered. If NPRINT is not positive, no special calls
!>       of FCN with IFLAG=0 are made.
!>
!>@param INFO @parblock is an integer output variable. If the user has
!>       terminated execution, INFFO is set to the (negative)
!>       value of IFLAG. See description of FCN. Otherwise,
!>       INFO is set as follows:
!>
!>       INFO=0  Improper input parameters.
!>
!>       INFO=1  Both actual and predicted relative reductions
!>                 in the sum of squares are at most FTOL.
!>
!>       INFO=2  Relative error between two consecutive iterates
!>                 is at most XTOL.
!>
!>       INFO=3  conditions for INFO=1 and INFO=2 both hold.
!>
!>       INFO=4  The cosine of the angle between FVEC and any
!>                 column of the Jacobian is at most GTOL in
!>                 absolute value.
!>
!>       INFO=5  number of calls to FCN with IFLAG=1 has
!>                 reached MAXFEV.
!>
!>       INFO=6  FTOL is too small. No further reduction in
!>                 the sum of squares is possible.
!>
!>       INFO=7  XTOL is too small. No further improvement in
!>                 the approximate solution X is possible.
!>
!>       INFO=8  GTOL is too small. FVEC is orthogonal to the
!>                 columns of the Jacobian to machine precision.
!>       @endparblock
!>@param NFEV is an integer output variable set to the number of
!>       calls to FCN with IFLAG=1.
!>
!>@param NJEV is an integer output variable set to the number of
!>       calls to NJEV with IFLAG=2.
!>
!>@param IPVT is an integer output array of length N. IPVT
!>       defines a permutation matrix p such that JAC*P=Q*R,
!>       where JAC is the final calculated Jacobian, Q is
!>       orthogonal (not stored), and R is upper triangular
!>       with diagonal elements of nonincreasing magnitude.
!>       column j of P is column IPVT(j) of the identity matrix.
!>
!>@param QTF is an output array of length N which contains
!>       the first N elements of the vector (Q transpose)*FVEC.
!>
!>@param WA1, WA2, and WA3 are work arrays of length N.
!>
!>@param WA4 is a work array of length M.
!
!     Subprograms called
!
!       User-supplied ...... FCN
!
!       MINPACK-supplied ... DPMPAR,ENORM,LMPAR,QRFAC
!
!       FORTRAN-supplied ... DABS,DMAX1,DMIN1,DSQRT,MOD
!
!     Argonne National Laboratory. MINPACK Project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!----------------------------------------------------------------------
!c
      integer i,iflag,j,l
      double precision actred,delta,dirder,epsmch,fnorm1,gnorm,par,
     *                 pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm
      PARAMETER (EPSMCH=RNDOFF)
      double precision enorm
!
      double precision ONE,P1,P5,P25,P75,P0001,ZERO
!
      character*1 trstr
      logical grdclc
!
      integer itrim
      external itrim
! ---------------------------
      data ONE,P1,P5,P25,P75,P0001,ZERO
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!######################################################################
!
!     epsmch is the machine precision.
!
!
      info=0
      iflag=0
      nfev=0
      njev=0
!
!----------------------------------------------------------------------
!     Check the input parameters for errors.
!----------------------------------------------------------------------
      if (n.le.0 .or. m.lt.n .or. ldfjac.lt.m
     *    .or. ftol.lt.ZERO .or. xtol.lt.ZERO .or. gtol.lt.ZERO
     *    .or. maxfev.le.0 .or.maxitr.le.0 .or. factor.le.ZERO) 
     *    go to 300
!
!**********************************************************************
!
!----------------------------------------------------------------------
!     Initialize Levenberg-Marquardt parameter and iteration counter
!----------------------------------------------------------------------
      par=ZERO
      iter=1
!
!----------------------------------------------------------------------
!     Evaluate the function at the starting point
!     and calculate its norm.
!----------------------------------------------------------------------
      iflag=1
      call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      nfev=1
      if (iflag.lt.0) go to 300
      fnorm=enorm(m,fvec)
!
!----------------------------------------------------------------------
! ********* Beginning of the outer loop *******************************
!----------------------------------------------------------------------
!
   30 continue
!
!----------------------------------------------------------------------
!        Calculate the Jacobian matrix (iflag=2)
!----------------------------------------------------------------------
         iflag=2
         call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         njev=njev+1
         nfev=nfev+n*njev
         if (iflag.lt.0) go to 300
!
!----------------------------------------------------------------------
!        If requested, call fcn to enable printing of iterates
!----------------------------------------------------------------------
         if (nprint.gt.0) then
           iflag=0
           if (mod(iter-1,nprint).eq.0)
     *        call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
           if (iflag.lt.0) go to 300
         end if
!
!----------------------------------------------------------------------
!        Compute the QR factorization of the Jacobian
!----------------------------------------------------------------------

         call qrfac(m,njcol,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!----------------------------------------------------------------------
!        On the first iteration, set each non-positive element of the
!        SCALE scaling array according to the norms of the columns of 
!        the initial Jacobian
!----------------------------------------------------------------------
         if (iter.eq.1) then
             do j=1, n
               if (scale(j).le.ZERO) then
                  diag(j)=wa2(j)
               else
                  diag(j)=wa2(j)/scale(j)
               end if
               if (diag(j).eq.ZERO) diag(j)=ONE
            end do
!
            if (itrace.ne.0) then
               write(itrace,1000) (tag(j)(:itrim(tag(j))),j=1,n)
               write(itrace,1001) (wa2(j),j=1,n)
               write(itrace,1002) (diag(j),j=1,n)
            end if
!
!----------------------------------------------------------------------
!        On the first iteration, calculate the norm of the scaled x
!        and initialize the trust region bound delta
!----------------------------------------------------------------------
           do j=1, n
             wa3(j)=diag(j)*x(j)
          end do
          xnorm=enorm(n,wa3)
          delta=factor*xnorm
          if (delta.eq.ZERO) delta=factor
          if (itrace.ne.0) write(itrace,1003) xnorm,delta,factor
       end if
!
!----------------------------------------------------------------------
!        Form (Q transpose)*fvec and store the first njcol components in
!        QtF.
!----------------------------------------------------------------------
         do i=1, m
            wa4(i)=fvec(i)
        end do
!
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
!
          fjac(j,j)=wa1(j)
          qtf(j)=wa4(j)
       end do
!
!----------------------------------------------------------------------
!        Compute the norm of the scaled gradient.
!----------------------------------------------------------------------
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
!
!----------------------------------------------------------------------
!        Test for convergence of the gradient norm
!----------------------------------------------------------------------
         if (gnorm.le.gtol) info=4
         if (info.ne.0) go to 300
!
!----------------------------------------------------------------------
!        Rescale diag array
!----------------------------------------------------------------------
           do j=1,n
               if (scale(j).gt.ZERO) then
                  temp=wa2(j)/scale(j)
                  diag(j)=dmax1(diag(j),temp)
               end if
            end do
!
!----------------------------------------------------------------------
!  ******** Beginning of the inner loop ******************************
!----------------------------------------------------------------------
        istep=0
        grdclc=.false.
  200   continue
!
!----------------------------------------------------------------------
!           Determine the Levenberg-Marquardt parameter.
!----------------------------------------------------------------------
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4,gnvec,gradf)
!
            grdclc=grdclc.or.(par.ne.ZERO)
!----------------------------------------------------------------------
!           Store the direction p and X + p. Calculate the norm of p.
!----------------------------------------------------------------------
            do j=1, n
               wa1(j)=-wa1(j)
               wa2(j)=x(j) + wa1(j)
               wa3(j)=diag(j)*wa1(j)
            end do
            pnorm=enorm(n,wa3)
!
!----------------------------------------------------------------------
!        On the first iteration, adjust the initial trust region bound
!        to the size of the initial step.
!----------------------------------------------------------------------
            trstr=' '
            if (iter.eq.1) then
               if (delta.gt.pnorm) then
                  trstr='*'
               else
                  trstr='s'
               end if
               delta=dmin1(delta,pnorm)
           end if
!
           if (istep.eq.0 .and. itrace.ne.0) then
              write (itrace,1012) iter,(tag(j),j=1,n)
              write (itrace,1013) fnorm,(x(j),j=1,n)
           end if
!
!----------------------------------------------------------------------
!           Evaluate the function at x + p and calculate its norm.
!----------------------------------------------------------------------
            iflag=1
            call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
            nfev=nfev+1
            if (iflag.lt.0) go to 300
            fnorm1=enorm(m,wa4)
            istep=istep+1
!
!----------------------------------------------------------------------
!           Compute the scaled actual reduction.
!----------------------------------------------------------------------
            actred=-ONE
            if (P1*fnorm1.lt.fnorm) actred=ONE-(fnorm1/fnorm)**2
!
!----------------------------------------------------------------------
!           Compute the scaled predicted reduction and
!           the scaled directional derivative.
!----------------------------------------------------------------------
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
!
!----------------------------------------------------------------------
!           Compute the ratio of the actual to the predicted
!           reduction.
!----------------------------------------------------------------------
            ratio=ZERO
            if (prered.ne.ZERO) ratio=actred/prered
!
!----------------------------------------------------------------------
!           Update the step bound.
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!             If actual reduction is too much smaller than the predicted
!             reduction (i.e. actred/prered ratio is too small)
!             the function is not well-approximated by a quadratic
!             equation. Reduce the size of the trust region by a
!             factor of 0.1 to 0.5 and increase the L-M parameter.
!----------------------------------------------------------------------
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
!
            else
!
!----------------------------------------------------------------------
!             If ratio of actual to predicted reduction is close to 1,
!             the quadratic model is a good approximation to the function,
!             and we can try increasing the trust region to twice the
!             last step size in order to check whether a better solution
!             is available. Otherwise, the size of the trust region
!             is left unchanged.
!----------------------------------------------------------------------
              if (par.eq.ZERO .or. ratio.ge.P75) then
                 delta=pnorm/P5
                 par=P5*par
                 temp=ONE/P5
              else
                 temp=ONE
              end if
            end if
!
            if (itrace.ne.0) write (itrace,1014) istep,par,ratio,
     #                       ONE/temp,trstr,fnorm1,(wa2(j)-x(j),j=1,n)
!
!----------------------------------------------------------------------
!           Test for successful iteration.
!----------------------------------------------------------------------
            if (ratio.ge.P0001) then
!
!----------------------------------------------------------------------
!           Successful iteration. Update X, FVEC, and their norms.
!----------------------------------------------------------------------
               do j=1, n
                  x(j)=wa2(j)
                  wa2(j)=diag(j)*x(j)
               end do
!
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
!
            end if
!----------------------------------------------------------------------
!           Tests for convergence.
!----------------------------------------------------------------------
            info=0
            if (dabs(actred).le.ftol .and. prered.le.ftol
     *          .and. P5*ratio.le.ONE) info=1
            if (delta.le.xtol*xnorm) info=info + 2
            if (info.ne.0) go to 300
!
!----------------------------------------------------------------------
!           Tests for termination and stringent tolerances.
!----------------------------------------------------------------------
            if (nfev.ge.maxfev) info=5
            if (iter.ge.maxitr) info=6
            if (dabs(actred).le.epsmch .and. prered.le.epsmch
     *          .and. P5*ratio.le.ONE) info=7
            if (delta.le.epsmch*xnorm) info=8
            if (gnorm.le.epsmch) info=9
            if (info.ne.0) go to 300
!----------------------------------------------------------------------
!           End of the inner loop. Repeat if iteration unsuccessful.
!----------------------------------------------------------------------
            if (ratio.lt.P0001) go to 200
!----------------------------------------------------------------------
!        End of the outer loop.
!----------------------------------------------------------------------
         go to 30
  300 continue
!----------------------------------------------------------------------
!     Termination, either normal or user-imposed.
!----------------------------------------------------------------------
      if (iflag.lt.0) info=10
      iflag=3
      if (info.gt.4) iflag=-iflag
      if (info.ne.0.and.info.ne.10.and.nprint.gt.0) 
     #     call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      return
!
! ##### format statements for trace printout ########################
!
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
