c  VERSION 1.1  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       ==================
c                       SUBROUTINE : EVCGF
c                       ==================
c
c     Subroutine version of EPRCGF program. (intended for NLSPMC)
c
c     It calculates the eigenvalues and eigenvectors for the general 
c     complex matrix using Conjugate gradient algorithm.
c
c     On Entry :
c
c        Parameters are passed through the common block /simparm/.
c        Pulse propagator is passed through the common block /basis/.
c        Off-diagonal basis vectors are passed through the common
c               block /basis/.
c
c     On Exit  :
c
c        nev   : number of important eigenvalues
c        eval  : important eigenvalues
c        evec  : important eigenvectors
c        ierr  : error flag
c
c     Includes :
c               nlsdim.inc
c               stdio.inc
c               eprprm.inc
c               parcom.inc
c               eprmat.inc
c               indexf.inc
c               wkspcm.inc
c
c     Uses :
c               matrxo.f
c               matrxd.f
c               cscg.f
c               csval.f
c               csvec.f
c
c*********************************************************************
c
      subroutine evcgf(psval,stvx,nev,eval,evec,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'eprmat.inc'
      include 'basis.inc'
      include 'wkspcm.inc'
c
      integer i,j,k,psval,ndim,nev,ierr,mxit,ival,ndone,nmin
      double precision zr,zi,zero,bshft,wtolcg,error,errorv,ermin
      complex*16 czero,cone,eval,evec,evtmp,stvx,y,w
      parameter (zero=0.D0,czero=(0.D0,0.D0),cone=(1.D0,0.D0))
      dimension eval(mxstep),evec(mxdim,mxegv),stvx(mxdim),y(mxdim),
     #          w(mxdim)
c
c#####################################################################
c
c---------------------------------------------------------------------
c     calculate the matrix elements
c---------------------------------------------------------------------
c

      ierr=0
      if (psval.ne.0) then
         ndim=ndimo
         call matrxo(ierr)
         if (idebug.ne.0) then
c                                ** store starting vector **
            open (unit=ludiskb,file='stvec.stvx',status='unknown',
     #           access='sequential',form='unformatted')
            write (ludiskb) (stvx(i),i=1,ndim)
            close (unit=ludiskb)
c                                ** store matrix elements **
            open (unit=ludiskb,file='matrx.mtxx',status='unknown',
     #           access='sequential',form='unformatted')
            write (ludiskb) (zdiag(1,i),zdiag(2,i),i=1,ndim)
            write (ludiskb) (jzmat(i),i=1,ndim+1) 
            write (ludiskb) (izmat(i),zmat(i),i=1,nelim)
            write (ludiskb) (kzmat(i),i=1,ndim+1)
            write (ludiskb) (izmat(mxel-i+1),zmat(mxel-i+1),i=1,nelre)
            close (unit=ludiskb)
         end if
      else
         ndim=ndimd
         call matrxd(ierr)
         if (idebug.ne.0) then
c                                ** store starting vector **
            open (unit=ludiskb,file='stvec.stvz',status='unknown',
     #           access='sequential',form='unformatted')
            write (ludiskb) (stvx(i),i=1,ndim)
            close (unit=ludiskb)
c                                ** store matrix elements **
            open (unit=ludiskb,file='matrx.mtxz',status='unknown',
     #           access='sequential',form='unformatted')
            write (ludiskb) (zdiag(1,i),zdiag(2,i),i=1,ndim)
            write (ludiskb) (jzmat(i),i=1,ndim+1) 
            write (ludiskb) (izmat(i),zmat(i),i=1,nelim)
            write (ludiskb) (kzmat(i),i=1,ndim+1)
            write (ludiskb) (izmat(mxel-i+1),zmat(mxel-i+1),i=1,nelre)
            close (unit=ludiskb)
         end if
      end if
c
      if (ierr.ne.0) then
         write (lulog,1000)
         write (luttyo,1000)
c                  ** make ierr negative for critical error **
         ierr=-abs(ierr)
         return
      end if
c
c---------------------------------------------------------------------
c     perform CG calc to obtain eigenvalues and weighting factors
c---------------------------------------------------------------------
c
      do 20 i=1,mxdim
 20      y(i)=czero
c                          * precondition the matrix for CG
      wtolcg=wtol
      if(psval.eq.1) then
         bshft=shifti+b0
      else if(psval.eq.-1) then
         bshft=shifti-b0
         write(*,*)'should never get here in evcgf'
         stop
      else
         bshft=shifti
         wtolcg=wtol*0.2d0
      end if
c
c call changed from nstep+1 to nstep.  RC 2/3/99.
      ermin=1.0d20
      call cscg(stvx,ndim,nstep,cgtol,dcmplx(shiftr,bshft),y,c2wsp1,
     #    cwsp1,cwsp2,cwsp3,cwsp4,cwsp5,w,ndone,error,ermin,nmin)
c
      if (ndone.le.0) then	! not converged
         ndone=abs(ndone)
         write (lulog,1010) abs(ndone),error
         write (lulog,*) 'CG min was, at step ',ermin,nmin
         write (*,*) 'CG min was at step ',ermin,nmin
         if(lulog.ne.luttyo)write(luttyo,1010)abs(ndone),error
      end if
c
      call csval(cwsp1,cwsp2,cwsp3,cwsp4,cwsp5,y,eval,w,ndone,ierr)
c
c---------------------------------------------------------------------
c     keep the eigenvalues with
c     (i)  weighting factor greater than wtolcg,
c     (ii) real part is smaller than 20.0 gauss
c         (fastly decaying eigenmodes do not contribute to the signal)
c---------------------------------------------------------------------
c
      nev=0
      do 30 i=1,ngood
        if ((cdabs(w(i)).gt.wtolcg).and.(dreal(eval(i)).lt.20.0d0))
     #                         then
          nev=nev+1
          eval(nev)=eval(i)
          w(nev)=w(i)
        end if
 30   continue
c
c---------------------------------------------------------------------
c     calculate the eigenvectors for nev important eigenvalues
c---------------------------------------------------------------------
c
      if ((idebug.ne.0).and.(ievec.eq.0)) write (ludeb,1020)
      do 100 ival=1,nev
c
         evtmp=eval(ival)
c
         mxit=10
c         if(ival.eq.3)write(*,*)'E2.05 ',cwsp1(1),cwsp2(1),
c     #		cwsp3(1),cwsp4(1),cwsp5(1),ndiag
c         if(ival.eq.3)write(*,*)'y2.5.1 ',cwsp1(1),cwsp1(2)
c         if(ival.eq.3)write(*,*)'y2.5.2 ',cwsp2(1),cwsp2(2)
c         if(ival.eq.3)write(*,*)'y2.5.3 ',cwsp3(1),cwsp3(2)
c         if(ival.eq.3)write(*,*)'y2.5.4 ',cwsp4(1),cwsp4(2)
c         if(ival.eq.3)write(*,*)'y2.5.5 ',cwsp5(5),cwsp5(2)
c ok at this point ...
         call csvec(cwsp1,cwsp2,cwsp3,y,evtmp,error,errorv,
     #        cwsp4,cwsp5,iwsp,ndiag,mxit)
c
         do 60 i=1,ndim
            evec(i,ival)=czero
            do 60 j=1,ndiag
               evec(i,ival)=evec(i,ival)+c2wsp1(i,j)*y(j)
 60      continue
         if (psval.eq.0) then
            k=1
            do 70 j=1,ndimo
               if (pid(j).eq.1) then
                  evec(j,ival)=-pp(k)*evec(k,ival)
                  k=k+1
               else if (pid(j).eq.2) then
                  evec(j,ival)=-pp(k)*evec(k,ival)
     #                         -pp(k+1)*evec(k+1,ival)
                  k=k+2
               else 
                  write(*,*)'in evcgf, pid error stopping. ',j,pid(j)
                  stop
               end if
 70         continue
         end if
c
         if (idebug.ne.0) then
            if (ievec.eq.0) then
               write (ludeb,1025) ival,eval(ival),w(ival)
            else
               write (ludeb,1030) ival,eval(ival),w(ival)
               write (ludeb,1040) (j,evec(j,ival),j=1,ndimo)
            end if
         end if
c
 100  continue
c
      return
c
c=====================================================================
c     format statements
c=====================================================================
c
 1000 format(2x,'*** ERROR in generating matrix elements ***')
 1010 format(2x,' CG did not converge in ',i4,' steps : error = ',
     #       g14.7)
 1020 format(/' ** Eigenvalues **'/)
 1025 format(2x,i3,2x,g14.7,1x,g14.7,2x,g14.7,1x,g14.7)
 1030 format(/,' ** Eigenvalue ',i3,' : ',g12.5,1x,g12.5,2x,g12.5,
     #       1x,g12.5,/)
 1040 format(5x,i3,2x,g14.7,1x,g14.7)
c
      end
