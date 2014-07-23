c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       =================
c                       SUBROUTINE CONVFT
c                       =================
c
c     This routine performs convolution with Gaussian inhomogeneous
c     broadening and subsequent FT on comlex 2D time-domain data to
c     generate frequency domain spectrum.  It is intended for use
c     with nonlinear least-squares applications.
c
c     On Entry :
c        xspec :  complex 2D spectrum in time domain
c
c     On Exit  :
c        cspec  :  complex spectrum in freq. domain
c
c     Includes :
c               nlsdim.inc
c               eprprm.inc
c               wkspcm.inc
c               physcn.inc
c
c     Uses :
c               fft.f
c               switch.f
c
c*********************************************************************
c
      subroutine convft(xspec,cspec,npt1,npt2)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'wkspcm.inc'
      include 'physcn.inc'
c     
      double precision zero,gtf,cfact,twopi
      parameter (zero=0.0d0)
c
      integer i,j,npt1,npt2,k
      double precision gib0,lib0,tdabs,expnt,gconv,t1,t2,
     #        linit1,lstept1,linit2,lstept2,frst,tdt
c
      complex*16 xspec(npt1,npt2),cspec(npt1,npt2),czero
      parameter (czero=(0.0d0,0.0d0))
c      real*8 spec(npt1,npt2)
c
c#####################################################################
c
      twopi=8.0d0*datan(1.0d0)
      gtf=g0*betae/(hbar*twopi)*1.0d-6
      cfact=twopi*gtf
c
      linit1=init1*1.0d-3
      lstept1=stept1*1.0d-3
      if (iexp.ge.1) then
         linit2=init2*1.0d-3
         lstept2=stept2*1.0d-3
      end if
c
      gib0=dabs(gib*cfact/twopi)
      lib0=dabs(lib*cfact)
c
c---------------------------------------------------------------------
c     Convolution with Gaussian and/or Lorentzian inhomogeneous broad.
c---------------------------------------------------------------------
c
 400  if ((gib0.gt.1.0d-13).or.(lib0.gt.1.0d-13)) then
c
         gib0=twopi*twopi*gib0*gib0/2.0d0
         do 410 i=1,npt1
            t1=linit1+(i-1)*lstept1
c                                     ** FID **             
            if (iexp.eq.0) then
               c2wsp1(i,1)=xspec(i,1)*dexp(-gib0*t1*t1)*dexp(-lib0*t1)
               go to 410
            end if
c
            if ((iexp.eq.2).or.(iexp.eq.4).or.(iexp.eq.5)) t1=zero
            do 420 j=1,npt2             
               t2=linit2+(j-1)*lstept2
               tdabs=dabs(t2+icomb*t1)
               expnt=gib0*tdabs*tdabs
               gconv=zero
               if (expnt.le.5.0d1)  gconv=dexp(-expnt)
               c2wsp1(i,j)=xspec(i,j)*gconv*dexp(-lib0*tdabs)
 420        continue
 410     continue
c
      else
         do 430 j=1,npt2
         do 430 i=1,npt1
            c2wsp1(i,j)=xspec(i,j)
 430     continue
      end if
      if ((iexp.eq.2).or.(iexp.eq.4)) then
        ipt=0
        do 470 i=1,npt1
           tdt=init1+dble(i-1)*stept1
c               * skip the region with asymmetric dead time
           if (ipt.eq.0) then
              if ((tdt.ge.init2).or.((init2-tdt).lt.1.0d-15))
     #        then
                 ipt=i
              else
                 go to 470
              end if
            end if
c         * get the threshold index when t2 > t1 starts
            do 480 j=1,npt2
               frst=init2+(j-1)*stept2
               if ((frst.ge.tdt).or.((tdt-frst).lt.1.d-4))
     #            go to 440
 480        continue
c                    * copy FID and zero-fill
 440        do 450 k=j,npt2
 450           c2wsp1(i-ipt+1,k-j+1)=c2wsp1(i,k)
            do 460 k=npt2-j+2,npt2
 460           c2wsp1(i-ipt+1,k)=czero
c                
 470    continue
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c---------------------------------------------------------------------
c     Frequency spectrum by Fourier transform
c---------------------------------------------------------------------
c
      if (iexp.ge.1) then
c                                 **  along t2
         do 500 i=1,npt1
            do 510 j=1,npt2
 510           cwsp3(j)=c2wsp1(i,j)
c divide first pt by 2 for baseline
	    cwsp3(1)=cwsp3(1)/2.0d0
            call fft(cwsp3,npt2)
            call switch(cwsp3,npt2)
            do 520 j=1,npt2
 520           c2wsp1(i,j)=cwsp3(j)
 500     continue
c                                 **  along t1
         do 550 j=1,npt2
            do 560 i=1,npt1
 560           cwsp3(i)=c2wsp1(i,j)
c divide first pt by 2 for baseline
	    cwsp3(1)=cwsp3(1)/2.0d0
            call fft(cwsp3,npt1)
            call switch(cwsp3,npt1)
            do 570 i=1,npt1
 570           cspec(i,j)=cwsp3(i)
 550     continue
c
      else
c
         do 580 i=1,npt1
 580        cwsp3(i)=c2wsp1(i,1)
c divide first pt by 2 for baseline
         cwsp3(1)=cwsp3(1)/2.0d0
         call fft(cwsp3,npt1)
         call switch(cwsp3,npt1)
         do 590 i=1,npt1
c
c not a real CW spectrum, but something to compare with the diagonal 
c slice from the 2D magnitude plot.  (iexp=0 case)
c
           cspec(i,1)=cwsp3(i)
 590     continue
c
      end if
c
      return
      end
