c Version 1.6  8/12/94
c**********************************************************************
c                       ===================
c                       SUBROUTINE : LBASIX
c                       ===================
c
c       This routine builds a list of the basis set indices in common 
c       block /indexl/ for the truncation parameters set in /eprdat/ 
c       or optionally reads them from the index file "<name>.ind"
c       which contains a pruned basis set produced by the EPRBL program. 
c       If one of the following conditions is true:
c         (1) The "new" flag is set nonzero on input
c         (2) <name>.ind does not exist 
c         (3) an error occurs reading <name.ind>
c       the routine builds the full basis set indicated
c       by the {lemx,lomx,kmn,kmx,mmn,mmx,ipnmx} MTS parameters.
c
c       Includes:
c               stddim.inc
c               stdio.inc
c               eprdat.inc
c               indexl.inc
c
c       Uses:
c               ipar.f
c
c**********************************************************************
c
      subroutine lbasix( ixname,lth,new )
      implicit none
      character*30 ixname
      integer lth,new,i
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
c
      integer lr,kr,krmn,krmx,mr,mrmn,mrmx,
     #        ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,
     #        iqnrmn,nrow,iparlr
c
      logical fexist
      integer ipar
      external ipar
c
c######################################################################
c
      if (new.eq.0) then
c
c----------------------------------------------------------------------
c     Check whether the basis set indices may be found in a file
c----------------------------------------------------------------------
c
         inquire(file=ixname(:lth),exist=fexist)
         if (fexist) then
c
c           -------------------------------------------
c           Read basis set information from index file
c           -------------------------------------------
            open (unit=ludisk,file=ixname(:lth),status='old',
     #           access='sequential',form='unformatted',err=9)
            read (ludisk,err=9) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
            read (ludisk,err=9) (l1(i),k1(i),m1(i),pi1(i),qi1(i),
     #           i=1,ndim)
            close (unit=ludisk)
            return
c
c           ----------------------------------------------------------
c           Error return: close file and continue to generate indices
c           ----------------------------------------------------------
 9          close (unit=ludisk)
         end if
c
      end if
c
c----------------------------------------------------------------------
c     Loop through all index limits to generate basis set
c----------------------------------------------------------------------
c
      nrow=0
c     ---------------------
c     *** loop over lr ***
c     ---------------------
      do 100 lr=0,lemx,ldelta
         iparlr=ipar(lr)
         if((iparlr.ne.1).and.(lr.gt.lomx)) go to 100
c
c     ---------------------
c     *** loop over kr ***
c     ---------------------
         krmx=min(kmx,lr)
         krmn=max(kmn,-lr)
         do 110 kr=krmn,krmx,kdelta
            if((kr.eq.0).and.(iparlr.lt.jkmn)) go to 110
c
c          ----------------------
c           *** loop over mr ***
c          ----------------------
            mrmx=min(mmx,lr)
            mrmn=max(mmn,-lr)
            do 120 mr=mrmn,mrmx
c
               if(ipsi0.eq.0) then
                  ipnrmn=abs(mr)
                  ipnrmx=abs(mr)
               else
                  ipnrmx=min(in2,ipnmx)
                  if(mr.eq.0 .and. jmmn.eq.1) then
                     ipnrmn=0
                  else
                     ipnrmn=-ipnrmx
                  end if
               end if
c
c           -------------------------    
c            *** loop over ipnr ***
c           -------------------------    
            do 130 ipnr=ipnrmn,ipnrmx
               if (mr.eq.0.and.ipnr.eq.0.and.iparlr.lt.jmmn) goto 130
c
c                -------------------------
c                 *** loop over iqnr ***
c                -------------------------
                  iqnrmx=in2-iabs(ipnr)
                  iqnrmn=-iqnrmx
                  do 140 iqnr=iqnrmn,iqnrmx,2
                     nrow=nrow+1
                     if (nrow.le.MXDIM) then
                        l1(nrow)=lr
                        k1(nrow)=kr
                        m1(nrow)=mr
                        pi1(nrow)=ipnr
                        qi1(nrow)=iqnr
                     end if
c
 140              continue
 130           continue
 120        continue
 110     continue
 100  continue
c
      ndim=nrow
c
      return
      end
