c Version 1.6  8/12/94
c----------------------------------------------------------------------
c                    =========================
c                         subroutine MTSL  
c                    =========================
c
c   Given an array representing the maximum projection of each basis
c   vector on the solution vector (array basis), the number of elements
c   (ndim) and a minimum projection defining the level at which a basis
c   vector is to be considered significant (bcut), this routine returns 
c   parameters defining the MTS (Minimum Truncation Set) of basis vector
c   indices; i.e. the minimum and maximum values of each index within
c   the set of "significant" basis vectors. The routine also returns the 
c   number of vectors qualifying as "significant" in nbas.
c
c   If the value of nmts is negative upon entry, the routine will also
c   count the number of vectors defined by the MTS parameters and return
c   it in nmts. In general, nmts need not equal nbas; i.e., the MTS still 
c   contains some insignificant vectors.
c
c   This routine replaces the function of program TNLL in EPRLL versions
c   through 1.3. It is intended to be called primarily by routine EPRBL.
c
c   Written by David Budil, 18 May 1993, Cornell Dept. of Chemistry
c
c   Uses:
c     ipar.f
c
c   Includes:
c     maxl.inc 
c     stddim.inc
c     eprdat.inc
c     indexl.inc
c----------------------------------------------------------------------
      subroutine mtsl( bsswt,nbasis,bcut,maxle,minle,maxlo,
     #     minlo,maxk,mink,maxm,minm,maxpi,minpi,nbas,nmts )
c
      implicit none
      integer nbasis,maxle,minle,maxlo,minlo,maxk,mink,maxm,
     #     minm,maxpi,minpi,nbas,nmts
c
      include 'maxl.inc'
      include 'stddim.inc'
      include 'indexl.inc'
      include 'eprdat.inc'
c
      integer i,ipnr,ipnrmn,ipnrmx,iqnr,iqnrmn,iqnrmx,kr,krmn,krmx,
     #        lr,mr,mrmn,mrmx,n,iparlr
      double precision bcut,bsswt(ndim)
c
      integer ipar
      external ipar
c
c----------------------------------------------------------------------
c
      minle=MXLVAL+1
      maxle=-MXLVAL-1
      minlo=MXLVAL
      maxlo=-MXLVAL
      mink=minle
      maxk=maxle
      minm=minle
      maxm=maxle
      minpi=in2+1
      maxpi=-in2-1
c
      nbas=0
      do i=1,nbasis 
         if (bsswt(i).ge.bcut) then
            nbas=nbas+1
            if (ipar(l1(i)).eq.1) then
               if (l1(i).ge.maxle) maxle=l1(i)
               if (l1(i).le.minle) minle=l1(i)
            else
               if (l1(i).ge.maxlo) maxlo=l1(i)
               if (l1(i).le.minlo) minlo=l1(i)
            end if
c
            if (k1(i).ge.maxk) maxk=k1(i)
            if (k1(i).le.mink) mink=k1(i)
c
            if (m1(i).ge.maxm) maxm=m1(i)
            if (m1(i).le.minm) minm=m1(i)
c
            if (pi1(i).ge.maxpi) maxpi=pi1(i)
            if (pi1(i).le.minpi) minpi=pi1(i)
         end if
      end do
c
c----------------------------------------------------------------------
c     Now count the number of elements in the MTS if needed 
c----------------------------------------------------------------------
      if (nmts.lt.0) then
         nmts=0
c
c     *** loop over lr ***
c
         do 100 lr=0,maxle,ldelta
            iparlr=ipar(lr)
            if((iparlr.ne.1).and.((lr.gt.maxlo).or.(lr.lt.minlo)) )
     #         go to 100
c
c
c     *** loop over kr ***
c
            krmx=min(maxk,lr)
            krmn=max(mink,-lr)
            do 110 kr=krmn,krmx,kdelta
               if((kr.eq.0).and.(iparlr.lt.jmmn)) go to 110
c     
c     *** loop over mr ***
c
               mrmx=min(maxm,lr)
               mrmn=max(minm,-lr)
               do 120 mr=mrmn,mrmx
                  if(ipsi0.eq.0) then
                     ipnrmn=abs(mr)
                     ipnrmx=abs(mr)
                  else
                     ipnrmx=min(in2,ipnmx)
                     if(mr.eq.0.and.jmmn.eq.1) then
                        ipnrmn=0
                     else
                        ipnrmn=-ipnrmx
                     end if
                  end if
c
c     *** loop over ipnr ***
c
                  do 130 ipnr=ipnrmn,ipnrmx
                     if((mr.eq.0).and.(ipnr.eq.0).and.(iparlr.lt.jmmn)) 
     #                    go to 130
c
c     *** loop over iqnr ***
c
                     iqnrmx=in2-iabs(ipnr)
                     iqnrmn=-iqnrmx
                     do 140 iqnr=iqnrmn,iqnrmx,2
c
                        nmts=nmts+1
c
c      *** end loop over rows ***
c
 140                 continue
 130              continue
 120           continue
 110        continue
 100     continue
      end if
c     
      return
      end
