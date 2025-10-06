! NLSL Version 1.5 beta 11/24/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine VARYC
!                    =========================
!
! vary <name> {  minimum <minval> maximum <maxval> scale <accuracy>
!                fdstep <step> }
!
!      name    : name of the parameter to be varied
!                Vary parameter individually for each spectrum in 
!                a series
!      minval  : Minimum permissible value for parameter
!      maxval  : Maximum permissible value for parameter
!      scale   : Factor by which to scale search vector for this parameter
!      step    : Relative step size for forward-differences approximation
!
!  NOTE: Minimum and maximum ARE NOT YET IMPLEMENTED
!
!      
! Special rules:
!
!  (1) Spectrum parameters may only be varied for all sites associated with
!      a given spectrum. These include:
!        PSI   -- tilt angle
!        PHASE -- spectral phase (absorption vs. dispersion)
!        LB    -- gaussian inhomogeneous broadening (orientation-independent)
!
!  (2) Note that B0, also a spectrum parameter may not be varied at all,
!      but it may be a series variable.
!
!  (3) GIB2 may only be varied for a PSI series or for a MOMD calculation

!  Other rules that will be implemented:
!
!    Shifting is disabled if the average g-value is allowed to float
!
!----------------------------------------------------------------------
      subroutine varyc( line )
!
      use nlsdim
      use eprprm
!      use prmeqv
      use parcom
      use lpnam
      use stdio
!
      implicit none
      character*80 line
!
      integer i,ibd,ix,ix2,lth
      double precision fval,prmn,prmx,prsc,step
      character token*30,prmID*9
      logical gib2OK
!
      integer NKEYWD
      parameter(NKEYWD=4)
!
      double precision ONE,ZERO
      parameter (ONE=1.0d0,ZERO=0.0d0)
!
      integer ipfind,indtkn,itrim
      logical ftoken
      external ftoken,ipfind,indtkn,itrim
      character*8 keywrd(NKEYWD)
      data keywrd /'MINIMUM','MAXIMUM','SCALE','FDSTEP'/
!
!    -------------------------------
!     Get the name of the parameter
!    -------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
!
      if (lth.le.0) then
         write(luttyo,1004)
         return
      end if
!
      call touppr(token,lth)
!
 1    ix=ipfind(token,lth)
!
!     ---------------------------------------
!      Check whether parameter may be varied 
!     ---------------------------------------
      if (ix.eq.0 .or. ix.gt.NVPRM) then
         write(luttyo,1002) token(:lth)
         return
      end if
!
      if (ix.lt.-100) then 
        prmID=alias2( -99-(IWXX+ix) )
      else if (ix.lt.0) then
        prmID=alias1( 1-(IWXX+ix) )
      else
        prmID=parnam(ix)
      end if
!
!     --- Get secondary index
!
      ix2=indtkn( line )
!
!     --------------------------------------------------
!      GIB2 may only be varied for MOMD calculations or 
!      if PSI is the series variable.
!     --------------------------------------------------
      if (ix.eq.IGIB2) then
         gib2OK=iser.eq.IPSI
         if (ix2.le.0) then
            do i=1,nsite
               gib2OK=gib2OK .or. (iparm(INORT,i).gt.1)
            end do
         else
            gib2OK=gib2OK .or. (iparm(INORT,ix2).gt.1)
         end if
!
         if (.not. gib2OK) then
            write(luttyo,1006)
            return
         end if
      end if
!
!    -------------
!    set defaults
!    -------------
      prmn=ZERO
      prmx=ZERO
      prsc=ONE
      step=1.0D-6
      ibd=0
!
!     --------------------
!      Look for a keyword
!     --------------------
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
!
!     ------------------------------------------------
!      No more tokens: vary the last parameter and exit
!     ------------------------------------------------
      if (lth.eq.0) then
         call addprm(ix,ix2,ibd,prmn,prmx,prsc,step,prmID)
         return
      end if
!
!     --------------------
!      Check keyword list
!     --------------------
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) goto 16
      end do
!
!     ----------------------------------------------------------------
!      Word was not recognized: add last parameter specified
!      and treat present token as a possible new parameter name 
!     ----------------------------------------------------------------
      call addprm(ix,ix2,ibd,prmn,prmx,prsc,step,prmID)
      go to 1
!
!     ----------------------------------------------------------------
!      Keyword found: convert next token and assign appropriate value
!     ----------------------------------------------------------------
   16 call gettkn(line,token,lth)
!                                        *** No value given
      if (lth.eq.0) then
        write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
        return
      end if
!
      if (ftoken(token,lth,fval)) then
!                                        *** MINIMUM keyword
        if (i.eq.1) then
          prmn=fval
          if (mod(ibd,2).eq.0) ibd=ibd+1
!                                        *** MAXIMUM keyword
        else if (i.eq.2) then
          prmx=fval
          ibd=ibd+2
!                                        *** SCALE keyword
        else if (i.eq.3) then
          prsc=fval
!                                        *** FDSTEP keyword
        else if (i.eq.4) then
          step=fval
        end if
!                               *** Illegal real value
      else
        write(luttyo,1001) token(:lth)
      end if
!
      go to 14
!
! ###### format statements ########################################
!
 1001 format('*** Real value expected: ''',a,''' ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** No value given for ''',a,''' ***')
 1004 format('*** Parameter name expected ***')
 1006 format('*** GIB2 may only be varied for a series of PSI angles',
     #       ' ***')

      end
