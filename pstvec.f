c NLSL Version 1.4 10/10/94
c**********************************************************************
c                    =========================
c                       subroutine PSTVEC
c                    =========================
c                      *** Pruning version ***
c
c       Calculates the starting vector for EPRLL
c
c       This differs from the standalone version in the addition of
c       the parameter ierr. It is used to report non-convergence of
c       Bessel function evaluations instead of halting when such
c       errors occur.   
c       
c       Notes:
c               1) The orientational integral is evaluated numerically 
c                  using the Curtis-Clenshaw-Romberg algorithm.
c                  For details, see Bruno's thesis, p. 554.
c
c       written by DJS 10-SEP-87
c
c NB:   Modified so that nelv is kept in common eprdat instead of
c       passed as a subroutine argument. DEB 22-OCT-91
c
c       Uses:
c               ipar.f
c               ccrint.f
c
c**********************************************************************
c
      subroutine pstvec(bss,v,ierr)
      implicit none
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'errmsg.inc'
      include 'rndoff.inc'
c
      integer lr,kr,iberr
      common /ifzdat/lr,kr,iberr
c
      integer bss(5,MXDIM),ierr
      double precision v(2,MXDIM)
c
      integer i,id,ioldlr,ioldkr,ioldmr,ioldpnr,iparlr,ipnr,ipnrsg,
     #        iqnr,jkr,jmr,krsgn,mr,mrsgn,nrow,nup
      double precision cnl,stv,stvlk,stvm,factor,dsq2,vnorm
      logical evenl,evenk,newlr,newkr,newmr,newpnr
c
      double precision ACCRCY,SMALL
      parameter (ACCRCY=1.0D-8,SMALL=1.0D-10)
c
      double precision ONE,ZERO
      parameter(ONE=1.0D0,ZERO=0.0D0)
c
      integer ipar
      double precision fz
      external fz,ipar
c
c######################################################################
c
c.... Debugging purposes only!
c      open (unit=20,file='nlpvec.tst',status='unknown',
c     #     access='sequential',form='formatted')
c............................
c
      dsq2=sqrt(2.0D0)
c
      nrow=0
      nelv=0
      vnorm=ZERO
c
c----------------------------------------------------------------------
c     *** loop over rows ***
c----------------------------------------------------------------------
c
      ioldlr=-1
      ioldkr=kmn-1
      ioldmr=mmn-1
      ioldpnr=-in2-1
c
      do 100 nrow=1,ndim
         lr=bss(1,nrow)
         krsgn=bss(2,nrow)
         mrsgn=bss(3,nrow)
         ipnrsg=bss(4,nrow)
         iqnr=bss(5,nrow)
c
         newlr=lr.ne.ioldlr
         newkr=newlr.or.(krsgn.ne.ioldkr)
         newmr=newkr.or.(mrsgn.ne.ioldmr)
         newpnr=newmr.or.(ipnrsg.ne.ioldpnr)
c
         if(newlr) then
            ioldlr=lr
            iparlr=ipar(lr)
            if (iparlr.eq.1) then
               evenl=.true.
               cnl=dsqrt(dble(2*lr+1))
            else
               evenl=.false.
               cnl=ZERO
            end if
         end if
c
        if (newkr) then
           ioldkr=krsgn
           kr=abs(krsgn)
           jkr=isign(1,krsgn)
           if (krsgn.eq.0) jkr=iparlr
           evenk=ipar(kr).eq.1
c
           if (evenl.and.evenk) then
c
c             --- No potential:  Only elements are for L,K=0
c
              if(lptmx.eq.0) then
                 if((lr.eq.0).and.(kr.eq.0)) then
                    stvlk=ONE
                 else
                    stvlk=ZERO
                 end if
c
c             --- Axial potential: Only elements are for K=0
c
              else if((kptmx.eq.0).and.(kr.ne.0)) then
                 stvlk=ZERO
c
c             --- Integration to find vector element
c
              else
                 call ccrint(ZERO,ONE,ACCRCY,SMALL,stvlk,nup,fz,id)
                 if (iberr.ne.0) ierr=BADBESS
                 if(kr.ne.0) then
                    factor=ONE
                    do 500 i=lr-kr+1,lr+kr
                       factor=factor*dble(i)
 500                continue
                    factor=ONE/dsqrt(factor)
                 else
                    factor=ONE/dsq2
                 end if
                 stvlk=factor*stvlk*cnl
              end if
           else
              stvlk=ZERO
           end if
        end if
c
        if (newmr) then
           ioldmr=mrsgn
           mr=abs(mrsgn)
           jmr=isign(1,mrsgn)
           if (mrsgn.eq.0) jmr=0
c
           if (mr.ne.0) then
              stvm=ZERO
           else
              stvm=stvlk
           end if
        end if
c
        if (newpnr) then
           ioldpnr=ipnrsg
           if (mr.eq.0) then
              ipnr=abs(ipnrsg)
              jmr=isign(1,ipnrsg)
              if (ipnr.eq.0) jmr=iparlr
           else
              ipnr=ipnrsg
           end if
c
           if (ipnr.ne.0.or.jmr.ne.1 .or.jkr.ne.1) then
              stv=ZERO
           else
              stv=stvm
           end if
        end if
c
        v(1,nrow)=stv
        vnorm=vnorm+stv*stv
c
c.................... Debugging purposes only!
c        if (abs(v(1,nrow)).gt.RNDOFF) then
c           write(20,7000) lr,krsgn,mrsgn,ipnrsg,iqnr,v(1,nrow)
c 7000      format('Re <v|',4(i3,','),i3,'> = ',2g14.7)
c        end if
c.............................................
c
 100  continue
c
c----------------------------------------------------------------------
c     normalize starting vector and zero out imaginary part
c----------------------------------------------------------------------
c
      nelv=0
      vnorm=ONE/dsqrt(vnorm)
      do i=1,ndim
         v(1,i)=v(1,i)*vnorm
         if(abs(v(1,i)).gt.RNDOFF) then
            nelv=nelv+1
         else
            v(1,i)=ZERO
         end if
         v(2,i)=ZERO
      end do
c
c----------------------------------------------------------------------
c     zero out remainder of vector
c----------------------------------------------------------------------
c
      do i=ndim+1,MXDIM
         v(1,i)=ZERO
         v(2,i)=ZERO
      end do
c
c..........Debugging purposes only!
c      close(20)
c...................................
c
      return
      end
