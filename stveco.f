c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c                       ===================
c                       SUBROUTINE : STVECO
c                       ===================
c
c     This subroutine calculates the starting vector for the off-diagonal
c     space.  It is intended for use with nonlinear lease-squares
c     applications.
c
c     On Entry :
c
c        Parameters are passed through the common block /eprprm/.
c        Basis set is passed through the common block /indexf/.
c
c     On Exit  :
c        The common block /stvcom/ contains the starting vector for
c        off-diagonal space.
c
c     Notes:
c        1) The orientational integral is evaluated numerically 
c           using the Curtis-Clenshaw-Romberg extrapolation algorithm.
c           For details, see Bruno's thesis.
c
c     10-APR-93 Sanghyuk Lee
c
c     Includes :
c        nlsdim.inc
c        eprprm.inc
c        indexf.inc
c        stvcom.inc
c        rndoff.inc
c
c     Uses:
c        ipar.f
c        ccrint.f
c
c**********************************************************************
c
      subroutine stveco
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'basis.inc'
      include 'stvcom.inc'
      include 'rndoff.inc'
c
      logical flr,fkr
c
      integer i,iparlr,nup,id,nrow,ipar,
     #        iper,jqer,lr,jkr,kr,jmr,mr,ipnr,iqnr
c
      double precision cnl,stvec,factor,dsq2,vnorm
c
      double precision acc,sml,fz,one,zero,two
      parameter (acc=1.0D-8,sml=1.0D-10)
      parameter (one=1.0D0,zero=0.0D0,two=2.0D0)
      complex*16 czero
      parameter (czero=(0.0d0,0.0d0))
c
      common /ifzdat/lr,kr
c
      external fz,ipar
c
c######################################################################
c
      dsq2=sqrt(2.0D0)
c
c----------------------------------------------------------------------
c             initialize counters
c----------------------------------------------------------------------
c
      vnorm=zero
c
c----------------------------------------------------------------------
c             *** loop over rows ***
c----------------------------------------------------------------------
c
      do 100 nrow=1,ndimo
c
        stvo(nrow)=czero
c
        jqer=jqe1(nrow)
        lr=l1(nrow)
        jkr=jk1(nrow)
        kr=k1(nrow)
        jmr=jm1(nrow)
        mr=m1(nrow)
        ipnr=pi1(nrow)
        iqnr=qi1(nrow)
c
        if ( (jqer.ne.0).or.(jkr.ne.1).or.(jmr.ne.1).or.
     #      (mr.ne.0).or.(ipnr.ne.0) ) go to 100
c
        iparlr=ipar(lr)
        flr=iparlr.eq.1
        cnl=zero
        if(flr) cnl=dsqrt(dble(2*lr+1))
c
        fkr=ipar(kr).eq.1
        stvec=zero
        if(flr.and.fkr) then
          if(lptmx.eq.0) then
            if((lr.eq.0).and.(kr.eq.0)) stvec=one
          else if((kptmx.eq.0).and.(kr.ne.0)) then
            stvec=zero
          else
            call ccrint(zero,one,acc,sml,stvec,nup,fz,id)
            if(kr.ne.0) then
              factor=one
              do 500 i=lr-kr+1,lr+kr
                factor=factor*dble(i)
 500          continue
              factor=one/dsqrt(factor)
            else
              factor=one/dsq2
            end if
            stvec=factor*stvec*cnl
          end if
        end if
c
c      **  center of loops  **
c
        stvo(nrow)=dcmplx(stvec,zero)
        vnorm=vnorm+stvec*stvec
c
 100  continue
c
c----------------------------------------------------------------------
c     normalize starting vector and zero out imaginary part
c----------------------------------------------------------------------
c
      nelv=0
      vnorm=one/dsqrt(vnorm)
      do 200 i=1,ndimo
        stvo(i)=stvo(i)*vnorm
        if(cdabs(stvo(i)).gt.rndoff) then
          nelv=nelv+1
        else
          stvo(i)=dcmplx(zero,zero)
        end if
 200  continue
c
c----------------------------------------------------------------------
c     zero out remainder of vector
c----------------------------------------------------------------------
c
c      do 210 i=ndimo+1,mxdim
c 210     stvo(i)=dcmplx(zero,zero)
c
c----------------------------------------------------------------------
c     return to caller
c----------------------------------------------------------------------
c

      return
      end
