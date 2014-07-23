c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       ================= 
c                       SUBROUTINE SPCALC
c                       =================
c
c     Subroutine version of EPRPS spectral calculation routine.
c
c     This routine is intended for use with nonlinear least-squares
c     applications.  The routine calculates 2D-spectrum for the 
c     parameters given in /eprprm/ using the eigenvalues and eigenvectors
c     already obtained by evesa routine.
c
c     !!!  No nuclear modulation allowed here.  !!!
c
c     On Entry :
c
c        scal  :  scale factor for MOMD
c                 (1/2 for end points in Simpson's integration)
c        evalx,evalz,evecx,evecz   :  eigenvalues and eigenvectors
c                 for off-diagonal & diagonal subspace
c
c     On Exit  :
c        xspec :  complex 2D spectrum in time domain
c
c     Includes :
c               nlsdim.inc
c               eprprm.inc
c               stvcom.inc
c               wkspcm.inc
c               physcn.inc
c
c     Uses :
c               zgemul, zgemv   (ESSL library routines)
c               ZGEMV and ZGEMM from IMSL replace the ESSL libs.
c
c*********************************************************************

c
      subroutine spcalc(scal,evalx,evalz,evecx,evecz,xspec,
     #                  npt1,npt2)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'stvcom.inc'
      include 'wkspcm.inc'
      include 'physcn.inc'
c     
      double precision zero,gtf,cfact,scal,twopi,omshft
      parameter (zero=0.0d0)
      complex*16 czero,ci,cone
      parameter (czero=(0.0D0,0.0D0),ci=(0.0D0,1.0D0),
     #           cone=(1.0d0,0.0d0))
c
      integer i,j,k,it1,it2,npt1,npt2,ntmp
      double precision wid,t1,t2,f1,f2,ltfix,
     #       linit1,lstept1,linit2,lstept2
c
      complex*16 sum1,sum2
      complex*16 evalx(mxegv),evalz(mxegv),
     #           levalx(mxegv),levaly(mxegv),levalz(mxegv),
     #           evecx(mxdim,mxegv),evecy(mxdim,mxegv),
     #           evecz(mxdim,mxegv),xspec(npt1,npt2),
     #           plsmix(mxegv,mxegv)
c
c#####################################################################
c
      twopi=8.0d0*datan(1.0d0)
      gtf=g0*betae/(hbar*twopi)*1.0d-6
      cfact=twopi*gtf
      omshft=twopi*shft
c
      ltfix=tfix*1.0d-3
      linit1=init1*1.0d-3
      lstept1=stept1*1.0d-3
      if (iexp.ge.1) then
         linit2=init2*1.0d-3
         lstept2=stept2*1.0d-3
      end if
c
c RC modification 11/18/96 to allow T2 to vary with mixing time
c with exponential dependence.  To do this, the old variables
c hwid and mwid have been redefined.   The old equation was:
c      wid=hwid+(mwid*twopi/gtf)*ltfix
c In the revised equation below,  hwid becomes the increase in line
c broadening (Gauss) between Tmix=0 and Tmix=infinity in usec and
c mwid describes how fast that change takes place with Tmix.
c mwid is in micro seconds.  To provide for a constant
c wid, use t2edi which has a similar effect.
c       given in run file. [Gauss]
      if (abs(mwid).gt.1.0D-10) then
        wid=abs(hwid)*(1.0d0-dexp(-ltfix/abs(mwid)))
      else
        wid=0.0d0
      endif
c
      do 10 i=1,nevo
         levalx(i)=(evalx(i)+wid)*cfact-ci*omshft
 10      levaly(i)=dconjg(levalx(i))
c
      do 20 j=1,nevo
      do 20 i=1,ndimo
 20      evecy(i,j)=dconjg(evecx(i,j))
c
      if (iexp.ge.3) then         
         do 30 i=1,nevd
 30         levalz(i)=evalz(i)*cfact
      end if
c
c---------------------------------------------------------------------
c     calculate the projection of each eigenvector on starting vector
c---------------------------------------------------------------------
c
      if (icomb.eq.1) then
         do 50 i=1,nevo
            levaly(i)=levalx(i)
            do 55 j=1,ndimo
 55            evecy(j,i)=evecx(j,i)
 50      continue
      end if
c
      do 60 i=1,nevo
         sum1=czero
         sum2=czero
         do 65 j=1,ndimo
            sum1=sum1+evecx(j,i)*stvo(j)
 65         sum2=sum2+evecy(j,i)*stvo(j)
         cwsp1(i)=sum1
 60      cwsp2(i)=sum2
c
c---------------------------------------------------------------------
c     calculate constant part in 2-D calculation for the given
c     experimental type
c---------------------------------------------------------------------
c
c                 *** 2 Pulse exp. : COSY or SECSY ***
c
      if ((iexp.eq.1).or.(iexp.eq.2)) then
c
c         call zgemul(evecy,mxdim,'t',evecx,mxdim,'n',
c     #        plsmix,mxegv,nevo,ndimo,nevo)
c
         call zgemm('t','n',nevo,nevo,ndimo,cone,evecy,mxdim,evecx,
     #	mxdim,czero,plsmix,mxegv)
c
      else if (iexp.ge.3) then
c
         if ((iexp.eq.3).or.(iexp.eq.4)) then
c
c                 *** 3 Pulse exp. : ELDOR, Echo-ELDOR ***
c
            do 90 i=1,nevd
 90            levalz(i)=cdexp(-levalz(i)*ltfix)
c
c            call zgemul(evecx,mxdim,'t',evecz,mxdim,'n',
c     #                  c2wsp1,mxdim,nevo,ndimo,nevd)
c
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecx,mxdim,evecz,
     #			mxdim,czero,c2wsp1,mxdim)

c
c            call zgemul(evecy,mxdim,'t',evecz,mxdim,'n',
c     #                  c2wsp2,mxegv,nevo,ndimo,nevd)
c
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecy,mxdim,evecz,
     #			mxdim,czero,c2wsp2,mxegv)
c

            do 230 i=1,nevo
            do 230 j=1,nevo
               sum1=czero
               do 235 k=1,nevd
 235              sum1=sum1+c2wsp2(i,k)*levalz(k)*c2wsp1(j,k)
 230           plsmix(i,j)=sum1*icomb
c
         else
c
c                 *** 3 Pulse exp. : Stim.-SECSY ***
c
            do 240 i=1,nevo
 240           levaly(i)=cdexp(-levaly(i)*tfix)
c
            do 245 i=1,nevo
 245           cwsp3(i)=cwsp2(i)*levaly(i)
c
c            call zgemul(evecy,mxdim,'t',evecz,mxdim,'n',
c     #                  c2wsp1,mxdim,nevo,ndimo,nevd)
c
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecy,mxdim,evecz,
     #			mxdim,czero,c2wsp1,mxdim)
c
c
            do 260 i=1,nevd
               sum1=czero
               do 265 j=1,nevo
 265              sum1=sum1+c2wsp1(j,i)*cwsp3(j)
 260           cwsp2(i)=sum1*icomb
c
c            call zgemul(evecz,mxdim,'t',evecx,mxdim,'n',
c     #                  plsmix,mxegv,nevd,ndimo,nevo)
c
            call zgemm('t','n',nevd,nevo,ndimo,cone,evecz,mxdim,evecx,
     #			mxdim,czero,plsmix,mxegv)

c
         end if
c
      end if
c
c=====================================================================
c     calculate 2-D spectrum in time-domain
c=====================================================================
c
c---------------------------------------------------------------------
c     step out t1 values
c---------------------------------------------------------------------
c
      do 300 it1=1,npt1
         t1=linit1+(it1-1)*lstept1
c
         if (iexp.ne.5) then
            do 305 i=1,nevo
 305           c2wsp1(i,it1)=cdexp(-levaly(i)*t1)*cwsp2(i)
         else
            do 307 i=1,nevd
 307           c2wsp1(i,it1)=cdexp(-levalz(i)*t1)*cwsp2(i)
         end if
 300  continue
c
      ntmp=nevo
      if (iexp.eq.5) ntmp=nevd
c
      if (iexp.ge.1) then
c         call zgemul(plsmix,mxegv,'t',c2wsp1,mxdim,'n',
c     #               c2wsp2,mxegv,nevo,ntmp,npt1)
         call zgemm('t','n',nevo,npt1,ntmp,cone,plsmix,mxegv,c2wsp1,
     #			mxdim,czero,c2wsp2,mxegv)
      else
c                                      *** iexp = 0: CW 1D-FID ***
c         call zgemv('t',nevo,npt1,cone,c2wsp1,mxdim,
c     #              cwsp2,1,czero,plsmix,1)	! returns plsmix
c						! = cone*c2wsp1'*cwsp2
c
         call zgemv('t',nevo,npt1,cone,c2wsp1,mxdim,
     #              cwsp2,1,czero,plsmix,1)	! returns plsmix
c						! = cone*c2wsp1'*cwsp2
c
         npt2=1
         go to 400
      end if
c
c---------------------------------------------------------------------
c     step out t2 values
c---------------------------------------------------------------------
c
c               *** t2 independent of t1 : no SECSY or echo-ELDOR ***
c
      if ((iexp.eq.1).or.(iexp.eq.3).or.(iexp.eq.5)) then
c
         do 330 it2=1,npt2
            t2=linit2+(it2-1)*lstept2
            if (iexp.eq.5) t2=ltfix+t2
c
            do 335 i=1,nevo
 335           c2wsp1(i,it2)=cwsp1(i)*cdexp(-levalx(i)*t2)
 330     continue
c
c         call zgemul(c2wsp2,mxegv,'t',c2wsp1,mxdim,'n',
c     #               plsmix,mxegv,npt1,nevo,npt2)
c
         call zgemm('t','n',npt1,npt2,nevo,cone,c2wsp2,mxegv,c2wsp1,
     #			mxdim,czero,plsmix,mxegv)
c
c
c               *** t2 dependent of t1 : SECSY or echo-ELDOR ***
c
      else if ((iexp.eq.2).or.(iexp.eq.4)) then
c
         do 350 it1=1,npt1
            t1=0.0d0                     
            do 360 it2=1,npt2
               t2=linit2+(it2-1)*lstept2
               t2=t2+t1
               do 365 i=1,nevo
 365              c2wsp1(i,it2)=cwsp1(i)*cdexp(-levalx(i)*t2)
 360        continue
c
c            call zgemv('t',nevo,npt2,cone,c2wsp1,mxdim,
c     #                 c2wsp2(1,it1),1,czero,plsmix(it1,1),mxegv)
            call zgemv('t',nevo,npt2,cone,c2wsp1,mxdim,
     #                 c2wsp2(1,it1),1,czero,plsmix(it1,1),mxegv)
 350     continue
c
      end if
c
c---------------------------------------------------------------------
c     Add spectrum
c---------------------------------------------------------------------
c
 400  continue
c
      do 410 j=1,npt2
      do 410 i=1,npt1
 410     xspec(i,j)=xspec(i,j)+scal*plsmix(i,j)
c
      return
      end
