c Version 1.5.1b 11/6/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SRCHC
c                    =========================
c
c   search  <parameter> { xtol <xtol> step <step> bound <bound> }
c
c       xtol:   Convergence tolerance for fitting parameter
c       step:   Size of first step to take away from initial parameter value  
c       bound:  Boundary for search (may not exceed initial value
c               +/- bound
c       maxfun: Maximum number of function evaluations
c
c   Uses
c      mnbrak
c      brent 
c      l1pfun
c      setprm
c      gettkn
c      ftoken
c      ipfind
c      spcpar
c      tcheck
c      indtkn
c      itrim
c----------------------------------------------------------------------
      subroutine srchc( line )
      implicit none
      character line*80
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'lmcom.inc'
      include 'parcom.inc'
      include 'iterat.inc'
      include 'errmsg.inc'
      include 'lpnam.inc'
      include 'stdio.inc'
c
      integer i,itmp,bflag,jx,jx1,jx2,lth,lu
      double precision ax,bx,cx,fa,fb,fc,fmin,fval,pmin
      character token*30,prmID*9,tagtmp*9
c
      integer NKEYWD
      parameter(NKEYWD=6)
c
      integer ipfind,indtkn,itrim
      double precision l1pfun,brent,getprm,gammq,residx,corrl,wtdres
      logical ftoken,tcheck,spcpar
      external l1pfun,brent,ftoken,spcpar,tcheck,
     #     ipfind,indtkn,itrim,getprm,residx,corrl,wtdres
c
      character*8 keywrd(NKEYWD)
      data keywrd / 'XTOL','FTOL','STEP','BOUND','MAXFUN','SRANGE'/
c
c######################################################################
c
c----------------------------------------------------------------------
c Get the name of the parameter
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
c
      ixp1p=ipfind(token,lth)
c
c----------------------------------------------------------------------
c     Check whether parameter may be varied 
c----------------------------------------------------------------------
      if (ixp1p.eq.0 .or. ixp1p.gt.NFPRM) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
      if (ixp1p.eq.iser) then
         write(luttyo,1011) token(:lth)
         return
      end if
c
      if (ixp1p.lt.-100) then 
        prmID=alias2( -99-(IWXX+ixp1p) )
      else if (ixp1p.lt.0) then
        prmID=alias1( 1-(IWXX+ixp1p) )
      else
        prmID=parnam(ixp1p)
      end if
c
c     --- Get secondary index
c
      ixs1p=indtkn( line )
c
c   --- Set bounds for site index of search parameter
c
      if (ixs1p.le.0) then
c
c        --- if parameter is to be searched globally for all sites, 
c            check that the initial values of the parameter are the
c            same for all sites
c
         do i=2,nsite
            if (getprm(ixp1p,i).ne.getprm(ixp1p,1)) then
               write (luout,1004) prmid(:itrim(prmID))
               if (luout.ne.luttyo) 
     #            write (luttyo,1004) prmid(:itrim(prmID))
               return
            end if
         end do
c
         jx1=1
         jx2=MXSITE
         if (spcpar(ixp1p)) jx2=MXSPC
      else
         jx1=ixs1p
         jx2=ixs1p
         write (prmID(itrim(prmid)+1:),1005) ixs1p
      end if
c
c  --- Check whether proper symmetry has been specified for
c      tensor components 
c
      lu=luttyo
      do jx=jx1,jx2
         if (.not.tcheck(ixp1p,ixs1p,prmID,lu)) return
         lu=0
      end do
c
      ixp1p=iabs( mod(ixp1p,100) )
c         
c----------------------------------------------------------------------
c  Look for a keyword
c----------------------------------------------------------------------
c
 13   call gettkn(line,token,lth)
      lth=min(lth,8)
c
c------------------------------------------------
c************************************************
c No more keywords: call the line search routine
c************************************************
c------------------------------------------------
      if (lth.eq.0) then
c
c        ---------------------------------------------------
c        Find a starting place for the bracketing procedure
c        ---------------------------------------------------
c
         if (ixs1p.le.0) then
            ax=fparm(ixp1p,1)
         else
            ax=fparm(ixp1p,ixs1p)
         end if
         bx=ax+pstep
c
c        -------------------------------------------------------------------
c        Swap indices for the search parameter into the first elements of
c        the parameter and site index arrays for the NLS x-vector
c        (This is so x can be used by the lfun routine as in a NLS procedure)
c        -------------------------------------------------------------------
c        NOTE: Nonsense is saved if no paramaters have been varied!
c        The swap-out should occur only if nprm .gt. 1
c
         itmp=ixpr(1)
         ixpr(1)=ixp1p
         ixp1p=itmp
         itmp=ixst(1)
         ixst(1)=ixs1p
         ixs1p=itmp
         tagtmp=tag(1)
         tag(1)=prmID
c
c        --------------------
c        Bracket the minimum
c        --------------------
c 
         call catchc( hltfit )
         warn=.true.
c
         write (luout,1009) prmID(:itrim(prmID))
         if (luout.ne.luttyo) write (luttyo,1009) prmID(:itrim(prmID))
c
         bflag=1
         iter=0
         fnmin=0.0d0
         call mnbrak(ax,bx,cx,fa,fb,fc,l1pfun,bflag,pbound)
c
c        ---------------------------------
c        User-terminated during bracketing
c        ---------------------------------
         if (bflag.lt.0) then
            write (luttyo,1007)
            go to 15
         else if (bflag.eq.NOMIN) then
            write (luttyo,1012) prmID(:itrim(prmID)),cx
            call setprm( ixpr(1), ixst(1), cx )
            fnorm=fc
            go to 14
         end if
c
c        -----------------------------------
c        Use Brent's method for minimization 
c        -----------------------------------
c            
         write (luout,1010) ax,cx
         if (luout.ne.luttyo) write (luttyo,1010) ax,cx
         bflag=1
         iter=0
         fnmin=fb
         fmin=brent(ax,bx,cx,fb,l1pfun,ptol,pftol,pmin,bflag)
c
         if (bflag.lt.0) then
c
c           ---------------------------------
c           User-terminated:report best value
c           ---------------------------------
            write (luttyo,1008) prmID,pmin
            if (luout.ne.luttyo) write (luout,1008) prmID,pmin
         else
c
c           ------------------------
c           Minimum found: report it   
c           ------------------------
            write (luttyo,1006) prmID,pmin
            if (luout.ne.luttyo) write (luout,1006) prmID,pmin
         end if
c
         fnorm=fmin
         call setprm( ixpr(1), ixst(1), pmin )
c
 14      if (iwflag.ne.0) then
            chisqr=fnorm*fnorm
            rdchsq=chisqr/float(ndatot-nprm)
            qfit=gammq( 0.5d0*(ndatot-nprm), 0.5*chisqr )
         else
            chisqr=wtdres( fvec,ndatot,nspc,ixsp,npts,rmsn )**2
            rdchsq=chisqr/float(ndatot)
         end if
c
         write(luout,2036) fnorm,chisqr,rdchsq,corrl(),residx()
         if (iwflag.ne.0) write (luout,2038) qfit
c 
         if (luout.ne.luttyo)
     #         write(luttyo,2036) fnorm,chisqr,rdchsq,corrl(),
     #                            residx()
c
c        --------------------------------------------------------
c        Restore the parameter/site index arrays and I.D. string
c        --------------------------------------------------------
 15      itmp=ixpr(1)
         ixpr(1)=ixp1p
         ixp1p=itmp
c
         itmp=ixst(1)
         ixst(1)=ixs1p
         ixs1p=itmp
         tag(1)=tagtmp
c
         call uncatchc( hltfit )
c
         return
      end if
c
c     -------------------
c      Check keyword list
c     -------------------
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) goto 16
      end do
c                                        *** Unrecognized keyword
      write (luttyo,1000) token(:lth)
      return
c
c      ----------------------------------------------------------
c      Keyword found: for keywords requiring an argument, convert 
c      next token and assign appropriate value
c      ------------------------------------------------------------
 16   call gettkn(line,token,lth)
c                                        *** No value given
      if (lth.eq.0) then
         write(luttyo,1003) keywrd(i)
         return
      end if
c
      if (ftoken(token,lth,fval)) then
c                                          *** XTOL keyword
         if (i.eq.1) then 
            ptol=fval
c                                          *** FTOL keyword
         else if (i.eq.2) then 
            pftol=fval
c                                          *** STEP keyword
         else if (i.eq.3) then 
            pstep=fval
c                                          *** BOUND keyword
         else if (i.eq.4) then
            pbound=fval
c                                          *** MAXFUN keyword
         else if (i.eq.5) then
            mxpitr=int(fval)
c                                          *** MAXFUN keyword
         else if (i.eq.6) then
            srange=fval/100.0
         end if
c                                      *** Illegal numeric value
      else
         write(luttyo,1001) token(:lth)
      end if
      go to 13
c
c ##### Format statements ########################################
c
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
