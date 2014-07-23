c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c
c                       ===================
c                       SUBROUTINE : FBASIS
c                       ===================
c     *** COMBINATION OF FBASO & FBASD ***
c
c          This routine builds a list of the basis set indices for
c     both diagonal and off-diagonal spaces in common block /indexf/
c     for the truncation parameters set in /eprprm/.
c          If the file, <fileid>.ind exists, it simply reads the basis
c     set.  This basis would be the pruned basis obtained by running
c     program eprbf.  If the file does not exist, it generates full
c     basis set within the given MTS parameters.
c
c     Includes:
c        nlsdim.inc
c        eprprm.inc
c        basis.inc
c        stdio.inc
c
c     Uses:
c        ipar.f
c
c**********************************************************************
c
      subroutine fbasis(ixname,igen,ierr)
c
c igen=0 to read from file, else generate full.
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'basis.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c
      logical fexist
      character*(NAMELG) ixname
      integer lu,igen,ierrtmp
c
      integer lr,jkr,kr,krmx,jmr,mr,mrmx,iper,iqer,iqermn,iqermx,
     #     ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,iqnrmn,nrow,i,j,ierr,ilm
c
      double precision one,sqrt2
c
      integer ipar
      external ipar
c
c######################################################################
c
      ierr=0
      ierrtmp=0
c
      if (igen.ne.0) go to 100	! generate full
c read from file:
      inquire(file=ixname,exist=fexist)
      if (fexist) then
c
c----------------------------------------------------------------------
c     Read basis set information from index file
c----------------------------------------------------------------------
c
	write(*,*)'opening ixname ',ludiskb, ixname
         open (unit=ludiskb,file=ixname,status='old',
     #         access='sequential',form='unformatted')
	write(*,*)'opened it'
         read (ludiskb) ndimoss(nbasis)
c fill next area:
         read (ludiskb) (mpi1(i),mqi1(i),ml1(i),mjk1(i),mk1(i),mjm1(i),
     #      mm1(i),i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1)
         if (nbasis+1 .le. nspectra*ncomps) then   ! ptr to next basis
           pidptr(nbasis+1)=pidptr(nbasis)+ndimoss(nbasis)
         end if
         close (unit=ludiskb)
         go to 200
      end if
c
      ierrtmp=1		! error, gen full basis (assumes have lr, etc.)
c
c----------------------------------------------------------------------
c     Check the magnetic parameters set in cmds.f
c----------------------------------------------------------------------
c
c  setbas should be false here since are asking for a basis.
      setbas=.false.
      ndimo=ndimoss(nbasis)
 100  call pcheck(luttyo,ierr)
c
      if (ierr.lt.0) then
         ierr=3
         return
      end if
c
      nrow=pidptr(nbasis)-1
c
c----------------------------------------------------------------------
c     *** loop over lr ***
c----------------------------------------------------------------------
c
      do 110 lr=0,lemx,ldelta
        if((ipar(lr).ne.1).and.(lr.gt.lomx)) go to 110
c
c----------------------------------------------------------------------
c     *** loop over jkr ***
c----------------------------------------------------------------------
c
        do 120 jkr=jkmn,1,2
c
c----------------------------------------------------------------------
c     *** loop over kr ***
c----------------------------------------------------------------------
c
        krmx=min(kmx,lr)
        do 130 kr=0,krmx,kdelta
          if((kr.eq.0).and.(ipar(lr).ne.jkr)) go to 130
c
c----------------------------------------------------------------------
c     *** loop over jmr ***
c----------------------------------------------------------------------
c
          do 140 jmr=jmmn,1,2
c
c----------------------------------------------------------------------
c     *** loop over mr ***
c----------------------------------------------------------------------
c
            mrmx=min(mmx,lr)
            do 150 mr=0,mrmx
c
c----------------------------------------------------------------------
c     *** loop over ipnr ***
c----------------------------------------------------------------------
c
              ipnrmx=min(in2,ipnmx)
              if (mr.eq.0) then
                ipnrmn=0
              else
                ipnrmn=-ipnrmx
              end if                
c
              do 160 ipnr=ipnrmn,ipnrmx
               if((mr.eq.0).and.(ipnr.eq.0).and.(ipar(lr).ne.jmr))
     #                                 go to 160
                if((ipsi0.eq.0).and.(ipnr.ne.mr)) go to 160
c
c----------------------------------------------------------------------
c     *** loop over iqnr ***
c----------------------------------------------------------------------
c
                iqnrmx=in2-iabs(ipnr)
                iqnrmn=-iqnrmx
                do 170 iqnr=iqnrmn,iqnrmx,2
c
                  nrow=nrow+1
                  mjqe1(nrow)=0
                  ml1(nrow)=lr
                  mjk1(nrow)=jkr
                  mk1(nrow)=kr
                  mjm1(nrow)=jmr
                  mm1(nrow)=mr
                  mpi1(nrow)=ipnr
                  mqi1(nrow)=iqnr
c
c----------------------------------------------------------------------
c     end loop over rows
c----------------------------------------------------------------------
c
 170              continue
 160            continue
 150          continue
 140        continue
 130      continue
 120    continue
 110  continue
c
      if (nbasis+1 .le. nspectra*ncomps) then   ! ptr to next basis
        pidptr(nbasis+1)=nrow+1		! location of next basis
      end if
        ndimoss(nbasis)=nrow-pidptr(nbasis)+1	! size of this one

c
 200  continue
c
c      ndimo=nrow
c
c**********************************************************************
c     Basis index for diagonal space (FBASD)
c**********************************************************************
c
      one=1.0D0
      sqrt2=dsqrt(2.0D0)
c
c----------------------------------------------------------------------
c       convert pruned off-diagonal basis set into diagonal one.
c----------------------------------------------------------------------
c
      j=dpidptr(nbasis)-1
      do 210 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
         j=j+1
         mdl1(j)=ml1(i)
         mdjk1(j)=mjk1(i)
         mdk1(j)=mk1(i)
         mdjm1(j)=0
         mdm1(j)=mm1(i)
         mdjqe1(j)=mjm1(i)
         mdpi1(j)=mpi1(i)
         mdqi1(j)=mqi1(i)
         if ((mpi1(i).eq.0).and.(mm1(i).eq.0)) then
            mpid(i)=1
            mpp(j)=sqrt2
         else
            mpp(j)=one
c
            j=j+1
            mdl1(j)=ml1(i)
            mdjk1(j)=mjk1(i)
            mdk1(j)=mk1(i)
            mdjm1(j)=0
            mdm1(j)=-mm1(i)
            mdjqe1(j)=mjm1(i)
            mdpi1(j)=-mpi1(i)
            mdqi1(j)=mqi1(i)
            mpid(i)=2
            ilm=ipar(ml1(i)+mm1(i))
            if (ilm.eq.mjm1(i)) then
               mpp(j)=one
            else
               mpp(j)=-one
            end if
         end if
 210  continue
         do 211 i=pidptr(nbasis)+ndimoss(nbasis),j
           mpid(i)=-1	! test for problems
 211     continue
      if (nbasis+1 .le. nspectra*ncomps) then   ! ptr to next basis
        dpidptr(nbasis+1)=j+1		! location of next basis
      end if
      ndimdss(nbasis)=j-dpidptr(nbasis)+1	! length of dp
c
      if (ierrtmp.eq.1) ierr=1
      if ((pidptr(nbasis).gt.mmxdim).or.(dpidptr(nbasis)
     #		.gt. 2*mmxdim)) ierr=2
c
      return
      end
