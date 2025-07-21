c NLSL Version 1.5 beta 11/24/95
c----------------------------------------------------------------------
c                    =========================
c                      subroutine VARYC
c                    =========================
c
c vary <name> {  minimum <minval> maximum <maxval> scale <accuracy>
c                fdstep <step> }
c
c      name    : name of the parameter to be varied
c                Vary parameter individually for each spectrum in 
c                a series
c      minval  : Minimum permissible value for parameter
c      maxval  : Maximum permissible value for parameter
c      scale   : Factor by which to scale search vector for this parameter
c      step    : Relative step size for forward-differences approximation
c
c  NOTE: Minimum and maximum ARE NOT YET IMPLEMENTED
c
c      
c Special rules:
c
c  (1) Spectrum parameters may only be varied for all sites associated with
c      a given spectrum. These include:
c        PSI   -- tilt angle
c        PHASE -- spectral phase (absorption vs. dispersion)
c        LB    -- gaussian inhomogeneous broadening (orientation-independent)
c
c  (2) Note that B0, also a spectrum parameter may not be varied at all,
c      but it may be a series variable.
c
c  (3) GIB2 may only be varied for a PSI series or for a MOMD calculation

c  Other rules that will be implemented:
c
c    Shifting is disabled if the average g-value is allowed to float
c
c----------------------------------------------------------------------
      subroutine varyc( line )
c
      use nlsdim
      use eprprm
c      use prmeqv
      use parcom
      use lpnam
      use stdio
c
      implicit none
      character*80 line
c
      integer i,ibd,ix,ix2,lth
      double precision fval,prmn,prmx,prsc,step
      character token*30,prmID*9
      logical gib2OK
c
      integer NKEYWD
      parameter(NKEYWD=4)
c
      double precision ONE,ZERO
      parameter (ONE=1.0d0,ZERO=0.0d0)
c
      integer ipfind,indtkn,itrim
      logical ftoken
      external ftoken,ipfind,indtkn,itrim
      character*8 keywrd(NKEYWD)
      data keywrd /'MINIMUM','MAXIMUM','SCALE','FDSTEP'/
c
c    -------------------------------
c     Get the name of the parameter
c    -------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
c
      if (lth.le.0) then
         write(luttyo,1004)
         return
      end if
c
      call touppr(token,lth)
c
 1    ix=ipfind(token,lth)
c
c     ---------------------------------------
c      Check whether parameter may be varied 
c     ---------------------------------------
      if (ix.eq.0 .or. ix.gt.NVPRM) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
      if (ix.lt.-100) then 
        prmID=alias2( -99-(IWXX+ix) )
      else if (ix.lt.0) then
        prmID=alias1( 1-(IWXX+ix) )
      else
        prmID=parnam(ix)
      end if
c
c     --- Get secondary index
c
      ix2=indtkn( line )
c
c     --------------------------------------------------
c      GIB2 may only be varied for MOMD calculations or 
c      if PSI is the series variable.
c     --------------------------------------------------
      if (ix.eq.IGIB2) then
         gib2OK=iser.eq.IPSI
         if (ix2.le.0) then
            do i=1,nsite
               gib2OK=gib2OK .or. (iparm(INORT,i).gt.1)
            end do
         else
            gib2OK=gib2OK .or. (iparm(INORT,ix2).gt.1)
         end if
c
         if (.not. gib2OK) then
            write(luttyo,1006)
            return
         end if
      end if
c
c    -------------
c    set defaults
c    -------------
      prmn=ZERO
      prmx=ZERO
      prsc=ONE
      step=1.0D-6
      ibd=0
c
c     --------------------
c      Look for a keyword
c     --------------------
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
c
c     ------------------------------------------------
c      No more tokens: vary the last parameter and exit
c     ------------------------------------------------
      if (lth.eq.0) then
         call addprm(ix,ix2,ibd,prmn,prmx,prsc,step,prmID)
         return
      end if
c
c     --------------------
c      Check keyword list
c     --------------------
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) goto 16
      end do
c
c     ----------------------------------------------------------------
c      Word was not recognized: add last parameter specified
c      and treat present token as a possible new parameter name 
c     ----------------------------------------------------------------
      call addprm(ix,ix2,ibd,prmn,prmx,prsc,step,prmID)
      go to 1
c
c     ----------------------------------------------------------------
c      Keyword found: convert next token and assign appropriate value
c     ----------------------------------------------------------------
   16 call gettkn(line,token,lth)
c                                        *** No value given
      if (lth.eq.0) then
        write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
        return
      end if
c
      if (ftoken(token,lth,fval)) then
c                                        *** MINIMUM keyword
        if (i.eq.1) then
          prmn=fval
          if (mod(ibd,2).eq.0) ibd=ibd+1
c                                        *** MAXIMUM keyword
        else if (i.eq.2) then
          prmx=fval
          ibd=ibd+2
c                                        *** SCALE keyword
        else if (i.eq.3) then
          prsc=fval
c                                        *** FDSTEP keyword
        else if (i.eq.4) then
          step=fval
        end if
c                               *** Illegal real value
      else
        write(luttyo,1001) token(:lth)
      end if
c
      go to 14
c
c ###### format statements ########################################
c
 1001 format('*** Real value expected: ''',a,''' ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** No value given for ''',a,''' ***')
 1004 format('*** Parameter name expected ***')
 1006 format('*** GIB2 may only be varied for a series of PSI angles',
     #       ' ***')

      end
