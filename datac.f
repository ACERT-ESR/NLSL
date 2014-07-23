c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine DATAC
c                    =========================
c
c Interprets a "data" command from the given line. The format of
c the command is as follows:
c
c data <dataid> { ascii|binary scale|noscale shift|noshift store}
c
c ** Special Caution necessary for NLSPMC **
c   The order of <dataid> in the fit of series spectra SHOULD
c   match the order of parameters specified in series command.
c
c   dataid : base name for datafile and associated output files:
c            <dataid>.SPC -- datafile
c            <dataid>.FIT -- fit spectrum
c   ascii|binary : specifies the format of the datafile
c            default (ascii)
c   scale|noscale : specifies that NLSPMC program should scale the
c            corresponding spectrum to fit.  The flag idepnd are
c            set according to the choice.  Refer comments in sscale.f.
c            default (scale)
c   shift|noshift : specifies that NLSPMC program should shift the
c            corresponding spectrum to fit.  The flag sishft is
c            set according to the choice.  Refer comments in pfunnew.f
c            before calling xshft.  default (shift)
c   store|nostore :  specifies that the actual data (may have been 
c            spline-smoothened) be stored in the file <datid>.nsp
c            in 2D plotting format.
c
c----------------------------------------------------------------------
      subroutine datac( line )
      implicit none
      character line*80,token*30
c
      include 'limits.inc'
      include 'datas.inc'
      include 'lmcomm.inc'
      include 'parms.inc'
      include 'names.inc'
      include 'wkspcm.inc'
      include 'stdio.inc'
c
      integer nkeywd
      parameter (nkeywd=8)
c
      integer i,iret,ispc,ival,ix,ixcw,j,lth,iform,istore,lu,
     #        isp,ispi,mpts,mptscw
      double precision  inif2,res2,inif1,res1,amin,amax,scfac,tog
c
      character*8 keywrd(nkeywd)
      data keywrd /'ASCII','BINARY','SCALE','NOSCALE',
     #             'SHIFT','NOSHIFT','STORE','NOSTORE'/
c
      character*20 wndoid(MXSPEC)
      logical itoken
      integer getd2d
      external itoken,getd2d
c
c----------------------------------------------------------------------
c  Check if the series parameters are properly set
c----------------------------------------------------------------------
      
      if (nsparm.lt.9) then
         write(luttyo,1100)nsparm
         return
      end if
c
c----------------------------------------------------------------------
c  Get the name of the datafile
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c  If (1) no data name was specified, or (2) the number of datafiles
c  already read in equals the number of spectra to be calculated (nser)
c  then reset the data buffer
c----------------------------------------------------------------------
      if (lth.eq.0.or.nspc.ge.nser) then
         nspc=0
         nptot=0
         write(luttyo,1020)
         if (lth.eq.0) return
      end if
c
      nspc=nspc+1
      dataid(nspc)=token
c
c ....initialization of default values
c
      iform=0
      ndata(nspc)=snpt1(nspc)*snpt2(nspc)
      idepnd(nspc)=0
      idepnd(nspc+1)=0
      sishft(nspc)=1
      istore=0
c
c----------------------------------------------------------------------
c     Look for a keyword
c----------------------------------------------------------------------
 5    call gettkn(line,token,lth)
c
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do 6 i=1,nkeywd
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
 6       continue 
c
         write (luttyo,1000) token(:lth)
         go to 5
c
c----------------------------------------------------------------------
c     Keyword found: assign appropriate value using next token
c----------------------------------------------------------------------
c
c --- Non-argument keywords: 
c     ASCII, BINARY, SCALE, NOSCALE, SHIFT, NOSHIFT
c
c                                         *** ASCII keyword
 7       if (i.eq.1) then
            iform=0
c                                         *** BINARY keyword
         else if (i.eq.2) then
            iform=1
            write (luttyo,1030)
            iform=0
c                                         *** SCALE keyword
         else if (i.eq.3) then
            idepnd(nspc)=0
c                                         *** NOSCALE keyword
         else if (i.eq.4) then
            idepnd(nspc)=1
c                                         *** SHIFT keyword
         else if (i.eq.5) then
            sishft(nspc)=1
c                                         *** NOSHIFT keyword
         else if (i.eq.6) then
            sishft(nspc)=0
c                                         *** STORE keyword
         else if (i.eq.7) then
            istore=1
c                                         *** NOSTORE keyword
         else if (i.eq.8) then
            istore=0
         end if
c
         go to 5
c
c----------------------------------------------------------------------
c     Read in datafile (end of data command)
c----------------------------------------------------------------------
      else
c
         call setdat( dataid(nspc) )
c
c    --- Check whether there is enough storage for the new data
c
         ix=nptot+1
         ixcw=nptotcw+1
         if ( (nptot+ndata(nspc)).gt.MXPT ) then
            write(luttyo,1050) MXPT	    
            nspc=nspc-1
            return
         end if
c
         lu=0
c
         iret=getd2d (dtname,iform,snpt1(nspc),snpt2(nspc),lu,
     #           siexp(nspc),sicomb(nspc),sstept1(nspc),sstept2(nspc),
     #           sratio(nspc),data(ix),cwdata(ixcw),rwsp,cwsp1,cwsp2,
     #           cwsp3,c2wsp1,c2wsp2)
c
         if (iret.ne.0) then
c                               *** Error opening/reading datafile
            if (iret.eq.-1) write(luttyo,1060) dtname(:lthdnm)
            if (iret.eq.-2) write(luttyo,1070) dtname(:lthdnm)
            if (iret.eq.-3) write(luttyo,1080) dtname(:lthdnm)
            if (iret.eq.-4) write(luttyo,1090) dtname(:lthdnm)
            nspc=nspc-1
            return
         end if
c
         if (istore.eq.1) then
c                               *** store the splined datafile
c	    tog=dble(snpt1(nspc)-1)/dble(snpt1(nspc))
c            inif1=-5.0d2*tog/sstept1(nspc)
c            res1=-2.0d0*inif1/dble(snpt1(nspc)-1)
c	    inif1=-5.0d2/sstept1(nspc)
c	    tog=dble(snpt2(nspc)-1)/dble(snpt2(nspc))
c            inif2=-5.0d2*tog/sstept2(nspc)
c            res2=-2.0d0*inif2/dble(snpt2(nspc)-1)
c	    inif2=-5.0d2/sstept2(nspc)
            amin=0.0d0
            amax=0.0d0
            do 45 j=ix,ix+ndata(nspc)-1
               if ( data(j).gt.amax ) amax=data(j)
c               if ( data(j).lt.amin ) amin=data(j)
 45         continue
c            call wrfit( data(ix),snpt1(nspc),snpt2(nspc),
c     #            inif2,res2,inif1,res1,amin,amax,nsname,lthdnm )
            call wrfit( data(ix),nspc,amin,amax,nsname )
c
         end if
c
         ixsp(nspc)=ix
         ixspcw(nspc)=ixcw
         sshft(nspc)=0.0d0
         wndoid(nspc)=dataid(nspc)
         wndoid(nspc)(10:10)=char(0)
         nptot=nptot+ndata(nspc)
         nptotcw=nptotcw+snpt2(nspc)
c
c----------------------------------------------------------------------
c        Scale the series spectra within reasonable range (smax=1000)
c----------------------------------------------------------------------
c if we have all the data sets:
         if (nspc.eq.nser) then
c
            isp=1
 100        if (uniflg) then

               ispi=isp
               ix=ixsp(isp)
               ixcw=ixspcw(isp)
               mpts=ndata(isp)
               mptscw=snpt1(isp)
c                              * search next independent spectra
 50            isp=isp+1
               if (idepnd(isp).eq.1) then
                  mpts=mpts+ndata(isp)
                  mptscw=mptscw+snpt1(isp)
                  go to 50
               end if
            else
               ispi=1
               isp=nspc+1
               ix=1
               mpts=nptot
               ixcw=1
               mptscw=nptotcw
            end if
c                              * obtain scale factor
            amax=0.0d0
            do 52 j=ispi,isp-1
 52            if (sratio(j).gt.amax) amax=sratio(j)
            scfac=1.0d3/amax
            do 54 j=ispi,isp-1
 54            sratio(j)=scfac
c                              * scale the spectra
            do 56 j=ix,ix+mpts-1
 56            data(j)=data(j)*scfac
            do 58 j=ixcw,ixcw+mptscw-1
 58            cwdata(j)=cwdata(j)*scfac
c
            if (isp.gt.nspc) then
                dataOK=.true.
                return
             end if
c
            go to 100
c
         end if
c
c######################################################################
c
      end if
c
      return
c
c #### format statements ########################################
c
 1000 format('*** Unrecognized DATA keyword: ''',a,''' ***')
 1020 format('*** Data buffer has been reset *** ')
 1030 format('*** Binary format for input data is not supported ***')
 1050 format('*** Maximum number of data points (',i4,') exceeded ***')
 1060 format(13x,'*** Error opening or reading datafile ''',a,
     #       ''' ***')
 1070 format(13x,'*** Too many data points in datafile ''',a,''' ***')
 1080 format(13x,'*** Data pts outside input range in ''',a,''' ***')
 1090 format(13x,'*** Can not extract cw-equivalent spectrum ***')
 1100 format('*** Series parameters not properly set ***',i3)
      end
