c----------------------------------------------------------------------
c                         ====================
c                          common block DATAS
c                         ====================
c     Experimental data used for least-squares fitting and associated
c     parameters.  It also keeps the variable in series option.
c
c     NOTE: File 'limits.inc' should be included before this file
c           in the program source code.
c
c  ** User-specified parameters **
c
c    *Set by data command
c     data           : Experimental data array
c     idepnd         : Dependency flag for scale factor (sfac)
c     sishft         : shift flag
c
c    *Set by series command
c     siexp          : Experimental type (0=FID, 1=COSY, 2=SECSY, 3=ELDOR
c                                         4=Echo-ELDOR, 5=Stim. SECSY)
c     sicomb         : Spectral combination
c     snpt1,snpt2    : Number of data points in t1 & t2
c     sinit1,sstept1 : Initial and step values in t1 (ns)
c     sinit2,sstept2 : Initial and step values in t2 (ns)
c     stfix          : Fixed time (mixing time in ELDOR)
c     datwgt         : Weight for this data set (1=standard weight)
c
c  ** parameters calculated during fit **
c
c     cwdata         : cw-equivalent spectra extracted from 2D data array
c     cwsim          : cw-equivalent spectra simulated
c     cwsimtot       : sum simulated cw-equivalent spectra over sites
c     sfac           : Scale factor for each spectrum in series option
c     sshft          : Field offset (shift spectrum)
c     ndata          : Number of data points in each spectrum
c     ixsp           : Starting index of each spectrum in data array
c     ixspcw         : Starting index for cw-eq. spectrum in cwdata array
c     nspc           : Number of spectra
c     nptot          : Total number of points (all spectra)
c     nptotcw        : Total number of cw data points (all spectra)
c     dataid         : Identification string for each data file
c     sratio         : scaling ratio for each independent experimental
c                      spectrum to have same maximum intensity
c     stfft1	     : specify if fft is to be zero filled in t1
c     stfft2	     : same in t2.
c
c----------------------------------------------------------------------
      character dataid*30
      integer siexp,sicomb,snpt1,snpt2,ndata,ixsp,ixspcw,idepnd,
     #        sishft,nspc,nptot,nptotcw
c    #        ,ssfft1,ssfft2

      double precision data,cwdata,cwsim,cwsimtot,sinit1,sstept1,
     #             sinit2,sstept2,stfix,datwgt,sfac,sshft,sratio
      complex*16 ccwsim,ccwsimtot
c 
      common /expdat/ data(mxpt),cwdata(mxnpt*MXSPEC),
     #	    cwsim(mxnpt*MXSPEC,MXSITE),cwsimtot(mxnpt*MXSPEC),
     #	    ccwsim(mxnpt*MXSPEC,MXSITE),ccwsimtot(mxnpt*MXSPEC),
     #      sinit1(MXSPEC),sstept1(MXSPEC),sinit2(MXSPEC),
     #      sstept2(MXSPEC),stfix(MXSPEC),datwgt(MXSPEC),sfac(MXSPEC),
     #	    sshft(MXSPEC),sratio(MXSPEC),
     #      siexp(MXSPEC),sicomb(MXSPEC),snpt1(MXSPEC),snpt2(MXSPEC),
     #      ndata(MXSPEC),ixsp(MXSPEC),ixspcw(MXSPEC),idepnd(MXSPEC+1),
     #      sishft(MXSPEC),nspc,nptot,nptotcw,dataid(MXSPEC)

c     #	    ssfft1(MXSPEC),ssfft2(MXSPEC),

c
