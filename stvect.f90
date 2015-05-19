c NLSL Version 1.4 10/10/94
c**********************************************************************
c                       ===================
c                       SUBROUTINE : STVECT
c                       ===================
c
c       subroutine for calculating the starting vector for eprll
c       
c       This differs from the standalone version in the addition of
c       the parameter ierr. It is used to report non-convergence of
c       Bessel function evaluations instead of halting when such
c       errors occur.   
c       
c       Notes:
c               1) The orientational integral is evaluated numerically 
c                  using the Curtis-Clenshaw-Romberg extrapolation algorithm.
c                  For details, see Bruno's thesis.
c
c       written by DJS 10-SEP-87
c
c       Includes:
c               nlsdim.inc
c               eprprm.inc
c               rndoff.inc
c
c       Uses:
c               ipar.f
c               ccrint.f
c
c**********************************************************************
c
      subroutine stvect(v,ierr)
c
      use nlsdim
      use eprprm
      use errmsg
      use rnddbl
c
      implicit none
      double precision v(2,MXDIM)
      integer ierr
c
      logical evenk,fmpi0
      integer i,iparkr,iparlr,jkr,jmr,nup,id,nrow,krmn,krmx,krsgn,
     #     ipnr,ipnrmx,ipnrmn,ipnrsg,iqnr,iqnrmx,iqnrmn,mr,mrmn,mrmx,
     #     mrsgn
c
      double precision cnl,stv,stvlk,stvm,factor,dsq2,vnorm
c
      double precision ACCRCY,SMALL
      parameter (ACCRCY=1.0D-8,SMALL=1.0D-10)
c
      double precision ONE,ZERO
      parameter(ONE=1.0D0,ZERO=0.0D0)
c
      integer lr,kr,iberr
      common /ifzdat/lr,kr,iberr
c
      integer ipar
      double precision fz
      external fz,ipar
c
c######################################################################
c
c.... Debugging purposes only!
c      open (unit=20,file='nlvec.tst',status='unknown',
c     #     access='sequential',form='formatted')
c............................
c
      do i=1,ndim
         v(1,i)=ZERO
         v(2,i)=ZERO
      end do
c
      dsq2=sqrt(2.0D0)
c
      nrow=0
      nelv=0
      vnorm=ZERO
c
c     --------------------
c     *** loop over lr ***
c     --------------------
      do 100 lr=0,lemx,ldelta
         iparlr=ipar(lr)
         if((lr.gt.lomx).and.(iparlr.ne.1)) go to 100
c
         if (iparlr.eq.1) then
            cnl=dsqrt(dble(2*lr+1))
         else
            cnl=ZERO
         end if
c
c        --------------------
c        *** loop over kr ***
c        --------------------
         krmx=min(kmx,lr)
         krmn=max(kmn,-lr)
         do 110 krsgn=krmn,krmx,kdelta
            kr=abs(krsgn)
            iparkr=ipar(kr)
            jkr=isign(1,krsgn)
            if(kr.eq.0) jkr=iparlr
            if(jkr.lt.jkmn) goto 110
c
c        --- Only vector elements are for even L, K and jK=1
c
            if(iparlr.eq.1 .and. iparkr.eq.1 .and. jkr.eq.1) then
c
c            --- No potential: only elements are for L,K=0
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
c
c           ----------------------
c           *** loop over mr ***
c           ----------------------
            mrmx=min0(mmx,lr)
            mrmn=max0(mmn,-lr)
            do 120 mrsgn=mrmn,mrmx
               mr=abs(mrsgn)
               jmr=isign(1,mrsgn)
c
               if (mr.ne.0) then
                  stvm=ZERO
               else
                  stvm=stvlk
               end if
c
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
c
c              ----------------------
c              *** loop over ipnr ***
c              ----------------------
               do 130 ipnrsg=ipnrmn,ipnrmx
                  ipnr=abs(ipnrsg)
                  if (mr.eq.0) then
                     jmr=isign(1,ipnrsg)
                     if (ipnrsg.eq.0) jmr=iparlr
                  end if
                  if(jmr.lt.jmmn) goto 130
c
                  if (ipnr.ne.0 .or. jmr.ne.1) then
                     stv=ZERO
                  else
                     stv=stvm
                  end if
c
                  iqnrmx=in2-iabs(ipnr)
                  iqnrmn=-iqnrmx
c
c                 ----------------------
c                 *** loop over iqnr ***
c                 ----------------------
                  do 140 iqnr=iqnrmn,iqnrmx,2
c
c======================================================================
c                  center of loops
c======================================================================
c
                     nrow=nrow+1
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
c======================================================================
c
c
 140              continue
 130           continue
 120        continue
 110     continue
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
         if(abs(v(1,i)).gt.rndoff) then
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
      do i=ndim+1,mxdim
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
