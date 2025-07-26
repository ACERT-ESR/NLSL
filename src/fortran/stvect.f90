! NLSL Version 1.4 10/10/94
!**********************************************************************
!                       ===================
!                       SUBROUTINE : STVECT
!                       ===================
!
!       subroutine for calculating the starting vector for eprll
!       
!       This differs from the standalone version in the addition of
!       the parameter ierr. It is used to report non-convergence of
!       Bessel function evaluations instead of halting when such
!       errors occur.   
!       
!       Notes:
!               1) The orientational integral is evaluated numerically 
!                  using the Curtis-Clenshaw-Romberg extrapolation algorithm.
!                  For details, see Bruno's thesis.
!
!       written by DJS 10-SEP-87
!
!       Includes:
!               nlsdim.inc
!               eprprm.inc
!               rndoff.inc
!
!       Uses:
!               ipar.f
!               ccrint.f
!
!**********************************************************************
!
      subroutine stvect(v,ierr)
!
      use nlsdim
      use eprprm
      use errmsg
      use rnddbl
!
      implicit none
      double precision v(2,MXDIM)
      integer ierr
!
      logical evenk,fmpi0
      integer i,iparkr,iparlr,jkr,jmr,nup,id,nrow,krmn,krmx,krsgn,
     #     ipnr,ipnrmx,ipnrmn,ipnrsg,iqnr,iqnrmx,iqnrmn,mr,mrmn,mrmx,
     #     mrsgn
!
      double precision cnl,stv,stvlk,stvm,factor,dsq2,vnorm
!
      double precision ACCRCY,SMALL
      parameter (ACCRCY=1.0D-8,SMALL=1.0D-10)
!
      double precision ONE,ZERO
      parameter(ONE=1.0D0,ZERO=0.0D0)
!
      integer lr,kr,iberr
      common /ifzdat/lr,kr,iberr
!
      integer ipar
      double precision fz
      external fz,ipar
!
!######################################################################
!
!.... Debugging purposes only!
!      open (unit=20,file='nlvec.tst',status='unknown',
!     #     access='sequential',form='formatted')
!............................
!
      do i=1,ndim
         v(1,i)=ZERO
         v(2,i)=ZERO
      end do
!
      dsq2=sqrt(2.0D0)
!
      nrow=0
      nelv=0
      vnorm=ZERO
!
!     --------------------
!     *** loop over lr ***
!     --------------------
      do 100 lr=0,lemx,ldelta
         iparlr=ipar(lr)
         if((lr.gt.lomx).and.(iparlr.ne.1)) go to 100
!
         if (iparlr.eq.1) then
            cnl=dsqrt(dble(2*lr+1))
         else
            cnl=ZERO
         end if
!
!        --------------------
!        *** loop over kr ***
!        --------------------
         krmx=min(kmx,lr)
         krmn=max(kmn,-lr)
         do 110 krsgn=krmn,krmx,kdelta
            kr=abs(krsgn)
            iparkr=ipar(kr)
            jkr=isign(1,krsgn)
            if(kr.eq.0) jkr=iparlr
            if(jkr.lt.jkmn) goto 110
!
!        --- Only vector elements are for even L, K and jK=1
!
            if(iparlr.eq.1 .and. iparkr.eq.1 .and. jkr.eq.1) then
!
!            --- No potential: only elements are for L,K=0
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
                     do i=lr-kr+1,lr+kr
                        factor=factor*dble(i)
                     end do
                     factor=ONE/dsqrt(factor)
                  else
                     factor=ONE/dsq2
                  end if
                  stvlk=factor*stvlk*cnl
               end if
            else
               stvlk=ZERO
            end if
!
!           ----------------------
!           *** loop over mr ***
!           ----------------------
            mrmx=min0(mmx,lr)
            mrmn=max0(mmn,-lr)
            do 120 mrsgn=mrmn,mrmx
               mr=abs(mrsgn)
               jmr=isign(1,mrsgn)
!
               if (mr.ne.0) then
                  stvm=ZERO
               else
                  stvm=stvlk
               end if
!
               if(ipsi0.eq.0) then
                  ipnrmn=mr
                  ipnrmx=mr
               else
                  ipnrmx=min0(in2,ipnmx)
                  if(mr.eq.0 .and. jmmn.eq.1) then
                     ipnrmn=0
                  else
                     ipnrmn=-ipnrmx
                  end if
               end if
!
!              ----------------------
!              *** loop over ipnr ***
!              ----------------------
               do 130 ipnrsg=ipnrmn,ipnrmx
                  ipnr=abs(ipnrsg)
                  if (mr.eq.0) then
                     jmr=isign(1,ipnrsg)
                     if (ipnrsg.eq.0) jmr=iparlr
                  end if
                  if(jmr.lt.jmmn) goto 130
!
                  if (ipnr.ne.0 .or. jmr.ne.1) then
                     stv=ZERO
                  else
                     stv=stvm
                  end if
!
                  iqnrmx=in2-iabs(ipnr)
                  iqnrmn=-iqnrmx
!
!                 ----------------------
!                 *** loop over iqnr ***
!                 ----------------------
                  do 140 iqnr=iqnrmn,iqnrmx,2
!
!======================================================================
!                  center of loops
!======================================================================
!
                     nrow=nrow+1
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
!======================================================================
!
!
 140              continue
 130           continue
 120        continue
 110     continue
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
         if(abs(v(1,i)).gt.rndoff) then
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
      do i=ndim+1,mxdim
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
