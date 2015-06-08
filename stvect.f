c Version 1.6  10/25/94
c**********************************************************************
c                    =========================
c                       subroutine STVECT
c                    =========================
c                      *** Pruning version ***
c
c       Calculates the starting vector for EPRLL
c       NEW BASIS: Uses antisymmetric K- and M- combinations!
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
c       Includes:
c               stddim.inc
c               eprdat.inc
c               indexl.inc
c               rndoff.inc
c
c       Uses:
c               ipar.f
c               ccrint.f
c
c**********************************************************************
c
      subroutine stvect(v)
      implicit none
c
      include 'stddim.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
      include 'rndoff.inc'
c
      double precision v(2,MXDIM)
c
      logical newlr,newkr,newmr,newpnr
      integer i,id,iparkr,iparlr,ioldlr,ioldkr,ioldmr,ioldpnr,
     #        ipnr,ipnrsg,iqnr,jkr,jmr,krsgn,mr,mrsgn,nrow,nup
      double precision cnl,stv,stvlk,stvm,factor,dsq2,vnorm
c
      double precision ACCRCY,SMALL
      parameter (ACCRCY=1.0D-8,SMALL=1.0D-10)
c
      double precision ONE,ZERO
      parameter(ONE=1.0D0,ZERO=0.0D0)
c
      integer lr,kr
      common /ifzdat/lr,kr
c
      integer ipar
      double precision fz
      external fz,ipar
c
c######################################################################
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
         lr=l1(nrow)
         krsgn=k1(nrow)
         mrsgn=m1(nrow)
         ipnrsg=pi1(nrow)
         iqnr=qi1(nrow)
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
               cnl=dsqrt(dble(2*lr+1))
            else
               cnl=ZERO
            end if
         end if
c
        if (newkr) then
           ioldkr=krsgn
           kr=abs(krsgn)
           iparkr=ipar(kr)
           jkr=isign(1,krsgn)
           if (krsgn.eq.0) jkr=iparlr
c
c         --- Only vector elements are for even L, K, and jK=+1
c
           if (iparlr.eq.1 .and. iparkr.eq.1 .and.jkr.eq.1) then
c
c             --- No potential:  Only vector elements are for L,K=0
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
                 if(kr.ne.0) then
                    factor=ONE
                    do i=lr-kr+1,lr+kr
                       factor=factor*dble(i)
                    end do
                    factor=ONE/dsqrt(factor)
                 else
                    factor=ONE/dsq2
                 end if
                 stvlk=factor*stvlk*cnl
              end if
c
c          --- L or K is odd: no starting vector element
c
           else
              stvlk=ZERO
           end if
        end if
c
        if (newmr) then
           ioldmr=mrsgn
           mr=abs(mrsgn)
           jmr=isign(1,mrsgn)
           if (mr.ne.0) then
              stvm=ZERO
           else
              stvm=stvlk
           end if
        end if
c
        if (newpnr) then
           ioldpnr=ipnrsg
c
           if (mr.eq.0) then
              ipnr=abs(ipnrsg)
              jmr=isign(1,ipnrsg)
              if (ipnr.eq.0) jmr=iparlr
           else
              ipnr=ipnrsg
           end if
           if (ipnr.ne.0.or.jmr.ne.1) then
              stv=ZERO
           else
              stv=stvm
           end if
        end if
c
        v(1,nrow)=stv
        vnorm=vnorm+stv*stv
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
      return
      end
