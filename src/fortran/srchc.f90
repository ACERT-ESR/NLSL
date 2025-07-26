! Version 1.5.1b 11/6/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine SRCHC
!                    =========================
!
!   search  <parameter> { xtol <xtol> step <step> bound <bound> }
!
!       xtol:   Convergence tolerance for fitting parameter
!       step:   Size of first step to take away from initial parameter value  
!       bound:  Boundary for search (may not exceed initial value
!               +/- bound
!       maxfun: Maximum number of function evaluations
!
!   Uses
!      mnbrak
!      brent 
!      l1pfun
!      setprm
!      gettkn
!      ftoken
!      ipfind
!      spcpar
!      tcheck
!      indtkn
!      itrim
!----------------------------------------------------------------------
      subroutine srchc( line )
!
      use nlsdim
      use eprprm
      use expdat
      use lmcom
      use parcom
      use iterat
      use errmsg
      use lpnam
      use stdio
!
      implicit none
      character line*80
!
      integer i,itmp,bflag,jx,jx1,jx2,lth,lu
      double precision ax,bx,cx,fa,fb,fc,fmin,fval,pmin
      character token*30,prmID*9,tagtmp*9
!
      integer NKEYWD
      parameter(NKEYWD=6)
!
      integer ipfind,indtkn,itrim
      double precision l1pfun,brent,getprm,gammq,residx,corrl,wtdres
      logical ftoken,tcheck,spcpar
      external l1pfun,brent,ftoken,spcpar,tcheck,
     #     ipfind,indtkn,itrim,getprm,residx,corrl,wtdres
!
      character*8 keywrd(NKEYWD)
      data keywrd / 'XTOL','FTOL','STEP','BOUND','MAXFUN','SRANGE'/
!
!######################################################################
!
!----------------------------------------------------------------------
! Get the name of the parameter
!----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
!
      ixp1p=ipfind(token,lth)
!
!----------------------------------------------------------------------
!     Check whether parameter may be varied 
!----------------------------------------------------------------------
      if (ixp1p.eq.0 .or. ixp1p.gt.NFPRM) then
         write(luttyo,1002) token(:lth)
         return
      end if
!
      if (ixp1p.eq.iser) then
         write(luttyo,1011) token(:lth)
         return
      end if
!
      if (ixp1p.lt.-100) then 
        prmID=alias2( -99-(IWXX+ixp1p) )
      else if (ixp1p.lt.0) then
        prmID=alias1( 1-(IWXX+ixp1p) )
      else
        prmID=parnam(ixp1p)
      end if
!
!     --- Get secondary index
!
      ixs1p=indtkn( line )
!
!   --- Set bounds for site index of search parameter
!
      if (ixs1p.le.0) then
!
!        --- if parameter is to be searched globally for all sites, 
!            check that the initial values of the parameter are the
!            same for all sites
!
         do i=2,nsite
            if (getprm(ixp1p,i).ne.getprm(ixp1p,1)) then
               write (luout,1004) prmid(:itrim(prmID))
               if (luout.ne.luttyo) 
     #            write (luttyo,1004) prmid(:itrim(prmID))
               return
            end if
         end do
!
         jx1=1
         jx2=MXSITE
         if (spcpar(ixp1p)) jx2=MXSPC
      else
         jx1=ixs1p
         jx2=ixs1p
         write (prmID(itrim(prmid)+1:),1005) ixs1p
      end if
!
!  --- Check whether proper symmetry has been specified for
!      tensor components 
!
      lu=luttyo
      do jx=jx1,jx2
         if (.not.tcheck(ixp1p,ixs1p,prmID,lu)) return
         lu=0
      end do
!
      ixp1p=iabs( mod(ixp1p,100) )
!         
!----------------------------------------------------------------------
!  Look for a keyword
!----------------------------------------------------------------------
!
 13   call gettkn(line,token,lth)
      lth=min(lth,8)
!
!------------------------------------------------
!************************************************
! No more keywords: call the line search routine
!************************************************
!------------------------------------------------
      if (lth.eq.0) then
!
!        ---------------------------------------------------
!        Find a starting place for the bracketing procedure
!        ---------------------------------------------------
!
         if (ixs1p.le.0) then
            ax=fparm(ixp1p,1)
         else
            ax=fparm(ixp1p,ixs1p)
         end if
         bx=ax+pstep
!
!        -------------------------------------------------------------------
!        Swap indices for the search parameter into the first elements of
!        the parameter and site index arrays for the NLS x-vector
!        (This is so x can be used by the lfun routine as in a NLS procedure)
!        -------------------------------------------------------------------
!        NOTE: Nonsense is saved if no paramaters have been varied!
!        The swap-out should occur only if nprm .gt. 1
!
         itmp=ixpr(1)
         ixpr(1)=ixp1p
         ixp1p=itmp
         itmp=ixst(1)
         ixst(1)=ixs1p
         ixs1p=itmp
         tagtmp=tag(1)
         tag(1)=prmID
!
!        --------------------
!        Bracket the minimum
!        --------------------
! 
         call catchc( hltfit )
         warn=.true.
!
         write (luout,1009) prmID(:itrim(prmID))
         if (luout.ne.luttyo) write (luttyo,1009) prmID(:itrim(prmID))
!
         bflag=1
         iter=0
         fnmin=0.0d0
         call mnbrak(ax,bx,cx,fa,fb,fc,l1pfun,bflag,pbound)
!
!        ---------------------------------
!        User-terminated during bracketing
!        ---------------------------------
         if (bflag.lt.0) then
            write (luttyo,1007)
            go to 15
         else if (bflag.eq.NOMIN) then
            write (luttyo,1012) prmID(:itrim(prmID)),cx
            call setprm( ixpr(1), ixst(1), cx )
            fnorm=fc
            go to 14
         end if
!
!        -----------------------------------
!        Use Brent's method for minimization 
!        -----------------------------------
!            
         write (luout,1010) ax,cx
         if (luout.ne.luttyo) write (luttyo,1010) ax,cx
         bflag=1
         iter=0
         fnmin=fb
         fmin=brent(ax,bx,cx,fb,l1pfun,ptol,pftol,pmin,bflag)
!
         if (bflag.lt.0) then
!
!           ---------------------------------
!           User-terminated:report best value
!           ---------------------------------
            write (luttyo,1008) prmID,pmin
            if (luout.ne.luttyo) write (luout,1008) prmID,pmin
         else
!
!           ------------------------
!           Minimum found: report it   
!           ------------------------
            write (luttyo,1006) prmID,pmin
            if (luout.ne.luttyo) write (luout,1006) prmID,pmin
         end if
!
         fnorm=fmin
         call setprm( ixpr(1), ixst(1), pmin )
!
 14      if (iwflag.ne.0) then
            chisqr=fnorm*fnorm
            rdchsq=chisqr/float(ndatot-nprm)
            qfit=gammq( 0.5d0*(ndatot-nprm), 0.5*chisqr )
         else
            chisqr=wtdres( fvec,ndatot,nspc,ixsp,npts,rmsn )**2
            rdchsq=chisqr/float(ndatot)
         end if
!
         write(luout,2036) fnorm,chisqr,rdchsq,corrl(),residx()
         if (iwflag.ne.0) write (luout,2038) qfit
! 
         if (luout.ne.luttyo)
     #         write(luttyo,2036) fnorm,chisqr,rdchsq,corrl(),
     #                            residx()
!
!        --------------------------------------------------------
!        Restore the parameter/site index arrays and I.D. string
!        --------------------------------------------------------
 15      itmp=ixpr(1)
         ixpr(1)=ixp1p
         ixp1p=itmp
!
         itmp=ixst(1)
         ixst(1)=ixs1p
         ixs1p=itmp
         tag(1)=tagtmp
!
         call uncatchc( hltfit )
!
         return
      end if
!
!     -------------------
!      Check keyword list
!     -------------------
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) goto 16
      end do
!                                        *** Unrecognized keyword
      write (luttyo,1000) token(:lth)
      return
!
!      ----------------------------------------------------------
!      Keyword found: for keywords requiring an argument, convert 
!      next token and assign appropriate value
!      ------------------------------------------------------------
 16   call gettkn(line,token,lth)
!                                        *** No value given
      if (lth.eq.0) then
         write(luttyo,1003) keywrd(i)
         return
      end if
!
      if (ftoken(token,lth,fval)) then
!                                          *** XTOL keyword
         if (i.eq.1) then 
            ptol=fval
!                                          *** FTOL keyword
         else if (i.eq.2) then 
            pftol=fval
!                                          *** STEP keyword
         else if (i.eq.3) then 
            pstep=fval
!                                          *** BOUND keyword
         else if (i.eq.4) then
            pbound=fval
!                                          *** MAXFUN keyword
         else if (i.eq.5) then
            mxpitr=int(fval)
!                                          *** MAXFUN keyword
         else if (i.eq.6) then
            srange=fval/100.0
         end if
!                                      *** Illegal numeric value
      else
         write(luttyo,1001) token(:lth)
      end if
      go to 13
!
! ##### Format statements ########################################
!
 1000 format('*** Unrecognized SEARCH keyword: ''',a,''' ***')
 1001 format('*** Numeric value expected: ''',a,''' ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** No value given for ''',a,''' ***')
 1004 format('*** ',a,' is not the same for all currently defined',
     #       ' sites ***')
 1005 format('(',i1,')')
 1006 format(/2x,'Minimum found at ',a,'= ',g12.6/)
 1007 format('*** Terminated by user during bracketing of minimum ***')
 1008 format('*** Terminated by user ***'/'Best point is ',a,'= ',f9.5/)
 1009 format(/2x,'Bracketing the minimum in ',a)
 1010 format(/2x,'Minimum is between ',f9.5,' and ',f9.5/)
 1011 format('*** ',a,' is series parameter: cannot be searched ***')
 1012 format('*** Minimum is at step bound ***'/2x,a,'=',f9.5/)
 2036 format(10x,'Residual norm= ',g13.6/
     #       10x,'Chi-squared=',g13.6,5x,'Reduced Chi-sq=',g13.6/
     #       10x,'Correlation = ',f8.5,5x,'Residual index =',f8.5/)
 2038 format(12x,'Goodness of fit from chi-squared (Q)=',g13.6)
      end
