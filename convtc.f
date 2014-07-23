c NLSPMC VERSION (VERSION 1.0) 2/5/99
c----------------------------------------------------------------------
c                         ====================
c                          subroutine CONVTC
c                         ====================
c
c     Subroutine to convert the tensor named on the passed command line
c     according to the symmetry specified by iflg. The following are
c     valid symmetries: 
c         1:  Cartesian
c         2:  Spherical
c         3:  Axial
c     See source file tensym.f for an explanation of these symmetries.
c
c     Uses:
c         gettkn, touppr (in strutl.f)
c         ipfind.f
c         rmvprm.f
c         tocart, tosphr, toaxil (in tensym.f)
c
c     Includes 
c         limits.inc
c         simparm.inc
c         parms.inc
c         stdio.inc
c         miscel.inc
c         names.inc
c         lpnam.inc
c
c 7/31/98 - modified to not use ipfind but rather to just for A, G or R.
c 	    Previous version didn't handle G correctly, (saw GIB rather 
c	    than GXX)
c 	Symmetries apply to all spectra, all components.
c
c---------------------------------------------------------------------- 

      subroutine convtc( line,iflg )
      implicit none
      include 'limits.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'stdio.inc'
      include 'miscel.inc'
      include 'names.inc'
      include 'lpnam.inc'
c
      character line*(LINELG), token*(WORDLG)
      integer i,iflg,ix,ixa,ixf,ixten,j,k,lth,lu,ispec,icomp,ident
      integer CARTESIAN,SPHERICAL,AXIAL
      parameter (CARTESIAN=1,SPHERICAL=2,AXIAL=3)
c
      integer ipfind
      external ipfind
c
c######################################################################
c
      if (iflg.lt.CARTESIAN .or. iflg.gt.AXIAL) iflg=CARTESIAN
c
      call gettkn(line,token,lth)
      call touppr(token,lth)
c
c                                          *** No tensor specified
      if (lth.eq.0) then
         write (luttyo,1000)
         return
      end if
c
c Check for specific tensors:
c
      if (token(:1).eq.'A') then
        ixa=IAXX
      elseif (token(:1).eq.'G') then
        ixa=IGXX
      elseif (token(:1).eq.'R') then
        ixa=IDX
      else
        write (luout,1001) token(1:1)
        if (luout.ne.luttyo) write (luttyo,1001) token(1:1)
        return
      end if
c
c----------------------------------------------------------------------
c     Tensor found: Check existing symmetry: a) all equal? b) set?
c----------------------------------------------------------------------
c      else
         ixf=IIGFLG+(ixa-IGXX)/3	! 0,1 or 2 added to IIGFLG
         ixten=IGXX+3*(ixf-IIGFLG)	! 0,3 or 6 added to IGXX
c----------------------------------------------------------------------
c check existing symmetry, insist equal for all sites, spectra:
c----------------------------------------------------------------------
         i=iparm(ixf,1,1)
         do 11 icomp=1,ncomps
           do 11 ispec=1,nspectra
             if (iparm(ixf,ispec,icomp) .ne. i) then
               write(luttyo,*)'in convtc, error, multiple symmetries '
               return
             end if
 11      continue
c----------------------------------------------------------------------
c        Symmetry not set yet: set it and exit
c----------------------------------------------------------------------
         if (iparm(ixf,1,1) .eq. 0) then
	   do 20 icomp=1,ncomps
	     do 20 ispec=1,nspectra
 20            iparm(ixf,ispec,icomp)=iflg
           write(luttyo,1002) token(1:1),symstr(iflg)
           return
c
c----------------------------------------------------------------------
c        Symmetry same as the one specified: exit
c----------------------------------------------------------------------
         else if (iparm(ixf,1,1).eq.iflg) then
	   write(luttyo,*)'in convtc, symmetry unchanged, 
     #                     nothing changed'
           return
         end if
c
c----------------------------------------------------------------------
c        Here when want new symmetry:
c        Remove any tensor components of any symmetry from
c        the list of variable parameters
c----------------------------------------------------------------------
         do 10 i=0,2	!all components of tensor to be changed
           j=ixten+i
	   do 30 icomp=1,ncomps
	     do 30 ispec=1,nspectra
               if (ixx(j,ispec,icomp).ne.0) then	! if we vary this one
                 ident=ixx(j,ispec,icomp)		! get xxx index
	         call rmvprm(ident)       ! remove it
c rmvprm when asked to remove a variable, removes all specs/comps
c where that variable applies, so this call should only happen for the 
c first occasion. 
	       end if
 30	   continue
 10      continue
c     
c----------------------------------------------------------------------
c        Now...convert the tensor symmetry!
c----------------------------------------------------------------------
         do 40 icomp=1,ncomps
           do 40 ispec=1,nspectra
            if (iflg.eq.2) then
              call tosphr( fparm(ixten,ispec,icomp),
     #		iparm(ixf,ispec,icomp) )
            else if (iflg.eq.3) then
              call toaxil( fparm(ixten,ispec,icomp),
     #		iparm(ixf,ispec,icomp) )
            else
              call tocart( fparm(ixten,ispec,icomp),
     #		iparm(ixf,ispec,icomp) )
            end if
              iparm(ixf,ispec,icomp)=iflg
c           end if
 40	 continue
      write (luout,1003) token(1:1),symstr(iflg)
      if (luout.ne.luttyo) write (luttyo,1003) token(1:1),symstr(iflg)
c     
c----------------------------------------------------------------------
c        Report new values for site 1 spectrum 1
c----------------------------------------------------------------------
      lu=luout
 300  continue
      write(lu,*)'in convtc ',ncomps,iflg
      do 400 icomp=1,ncomps
        write (lu,*)'For component # ',icomp
        if (ixf.eq.IIGFLG) then
c                                                  * g tensor
          if (iflg.eq.1) then
            write (lu,1011) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else if (iflg.eq.2) then
            write (lu,1012) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else
            write (lu,1013) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          end if
        else if (ixf.eq.IIGFLG+1) then
c                                                  * hf tensor
          if (iflg.eq.1) then
            write (lu,1021) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else if (iflg.eq.2) then
            write (lu,1022) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else
            write (lu,1023) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          end if
        else
c                                                  * r tensor
          if (iflg.eq.1) then
            write (lu,1031) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else if (iflg.eq.2) then
            write (lu,1032) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          else
            write (lu,1033) icomp,(fparm(ixten+i,1,icomp),i=0,2)
          end if
        end if
 400  continue
      if (lu.ne.luttyo) then
         lu=luttyo
         go to 300
      end if
c
      return
c     
 1000 format('*** No tensor specified ***')
 1001 format('*** Unknown tensor: ''',a,''' ***')
 1002 format('*** ',a1,' tensor set for all spec,comps to ',a,' ***')
 1003 format('*** ',a1,' tensor converted to ',a,' ***')
 1011 format('spec 1, site ',i2,': gxx,gyy,gzz = ',2(f9.6,','),f9.6)
 1012 format('spec 1, site ',i2,': g1,g2,g3 = ',2(f9.6,','),f9.6)
 1013 format('spec 1, site ',i2,': gprp,grhm,gpll = ',2(f9.6,','),f9.6)
 1021 format('spec 1, site ',i2,': Axx,Ayy,Azz = ',2(f9.6,','),f9.6)
 1022 format('spec 1, site ',i2,': A1,A2,A3 = ',2(f9.6,','),f9.6)
 1023 format('spec 1, site ',i2,': Aprp,Arhm,Apll = ',2(f9.6,','),f9.6)
 1031 format('spec 1, site ',i2,': log(Rx,Ry,Rz) = ',2(f9.6,','),f9.6)
 1032 format('spec 1, site ',i2,': log(Rbar,N,Nxy) = ',2(f9.6,','),
     #		f9.6)
 1033 format('spec 1, site ',i2,': log(Rprp,Rrhm,Rpll) = ',
     #		2(f9.6,','),f9.6)
      end
