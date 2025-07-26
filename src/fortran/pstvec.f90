! NLSL Version 1.4 10/10/94
!**********************************************************************
!                    =========================
!                       subroutine PSTVEC
!                    =========================
!                      *** Pruning version ***
!
!       Calculates the starting vector for EPRLL
!
!       This differs from the standalone version in the addition of
!       the parameter ierr. It is used to report non-convergence of
!       Bessel function evaluations instead of halting when such
!       errors occur.   
!       
!       Notes:
!               1) The orientational integral is evaluated numerically 
!                  using the Curtis-Clenshaw-Romberg algorithm.
!                  For details, see Bruno's thesis, p. 554.
!
!       written by DJS 10-SEP-87
!
! NB:   Modified so that nelv is kept in common eprdat instead of
!       passed as a subroutine argument. DEB 22-OCT-91
!
!       Uses:
!               ipar.f
!               ccrint.f
!
!**********************************************************************
!
      subroutine pstvec(bss,v,ierr)
!
      use nlsdim
      use eprprm
      use errmsg
      use rnddbl
!
      implicit none
      integer lr,kr,iberr
      common /ifzdat/lr,kr,iberr
!
      integer bss(5,MXDIM),ierr
      double precision v(2,MXDIM)
!
      integer i,id,ioldlr,ioldkr,ioldmr,ioldpnr,iparlr,ipnr,ipnrsg,
     #        iqnr,jkr,jmr,krsgn,mr,mrsgn,nrow,nup
      double precision cnl,stv,stvlk,stvm,factor,dsq2,vnorm
      logical evenl,evenk,newlr,newkr,newmr,newpnr
!
      double precision ACCRCY,SMALL
      parameter (ACCRCY=1.0D-8,SMALL=1.0D-10)
!
      double precision ONE,ZERO
      parameter(ONE=1.0D0,ZERO=0.0D0)
!
      integer ipar
      double precision fz
      external fz,ipar

      cnl=ZERO
      stv=ZERO
      stvlk=ZERO
      stvm=ZERO
      evenl=.false.
      evenk=.false.
      mr=0
      jmr=0
      jkr=0
      iparlr=0
!
!######################################################################
!
!.... Debugging purposes only!
!      open (unit=20,file='nlpvec.tst',status='unknown',
!     #     access='sequential',form='formatted')
!............................
!
      dsq2=sqrt(2.0D0)
!
      nrow=0
      nelv=0
      vnorm=ZERO
!
!----------------------------------------------------------------------
!     *** loop over rows ***
!----------------------------------------------------------------------
!
      ioldlr=-1
      ioldkr=kmn-1
      ioldmr=mmn-1
      ioldpnr=-in2-1
!
      do 100 nrow=1,ndim
         lr=bss(1,nrow)
         krsgn=bss(2,nrow)
         mrsgn=bss(3,nrow)
         ipnrsg=bss(4,nrow)
         iqnr=bss(5,nrow)
!
         newlr=lr.ne.ioldlr
         newkr=newlr.or.(krsgn.ne.ioldkr)
         newmr=newkr.or.(mrsgn.ne.ioldmr)
         newpnr=newmr.or.(ipnrsg.ne.ioldpnr)
!
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
!
        if (newkr) then
           ioldkr=krsgn
           kr=abs(krsgn)
           jkr=isign(1,krsgn)
           if (krsgn.eq.0) jkr=iparlr
           evenk=ipar(kr).eq.1
!
           if (evenl.and.evenk) then
!
!             --- No potential:  Only elements are for L,K=0
!
              if(lptmx.eq.0) then
                 if((lr.eq.0).and.(kr.eq.0)) then
                    stvlk=ONE
                 else
                    stvlk=ZERO
                 end if
!
!             --- Axial potential: Only elements are for K=0
!
              else if((kptmx.eq.0).and.(kr.ne.0)) then
                 stvlk=ZERO
!
!             --- Integration to find vector element
!
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
!
        if (newmr) then
           ioldmr=mrsgn
           mr=abs(mrsgn)
           jmr=isign(1,mrsgn)
           if (mrsgn.eq.0) jmr=0
!
           if (mr.ne.0) then
              stvm=ZERO
           else
              stvm=stvlk
           end if
        end if
!
        if (newpnr) then
           ioldpnr=ipnrsg
           if (mr.eq.0) then
              ipnr=abs(ipnrsg)
              jmr=isign(1,ipnrsg)
              if (ipnr.eq.0) jmr=iparlr
           else
              ipnr=ipnrsg
           end if
!
           if (ipnr.ne.0.or.jmr.ne.1 .or.jkr.ne.1) then
              stv=ZERO
           else
              stv=stvm
           end if
        end if
!
        v(1,nrow)=stv
        vnorm=vnorm+stv*stv
!
!.................... Debugging purposes only!
!        if (abs(v(1,nrow)).gt.RNDOFF) then
!           write(20,7000) lr,krsgn,mrsgn,ipnrsg,iqnr,v(1,nrow)
! 7000      format('Re <v|',4(i3,','),i3,'> = ',2g14.7)
!        end if
!.............................................
!
 100  continue
!
!----------------------------------------------------------------------
!     normalize starting vector and zero out imaginary part
!----------------------------------------------------------------------
!
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
!
!----------------------------------------------------------------------
!     zero out remainder of vector
!----------------------------------------------------------------------
!
      do i=ndim+1,MXDIM
         v(1,i)=ZERO
         v(2,i)=ZERO
      end do
!
!..........Debugging purposes only!
!      close(20)
!...................................
!
      return
      end
