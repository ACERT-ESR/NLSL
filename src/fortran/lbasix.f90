c  Version 1.4 (NLS)
c**********************************************************************
c                       ===================
c                       SUBROUTINE : LBASIX
c                       ===================
c
c       This routine builds a list of the basis set indices in common 
c       block /indexl/ for the truncation parameters set in /eprdat/ 
c       or reads them from the index file "<name>.ind" which contains 
c       a pruned basis set produced by the EPRBL program. 
c
c       If one of the following conditions is true:
c         (1) The "new" flag is set nonzero on input
c         (2) the input filename is blank
c
c       this routine routine builds the full basis set indicated
c       by the truncation parameters specified in the MTS array.
c       These are as follows:
c          mts(1)  lemx    Maximum even L
c          mts(2)  lomx    Maximum odd L
c          mts(3)  kmn     Minimum K (K<0 => jK = -1)
c          mts(4)  kmx     Maximum K
c          mts(5)  mmn     Minimum M (M<0 => jM = -1)
c          mts(6)  mmx     Maximum M
c          mts(7)  ipnmx   Maximum pI
c          mts(8)  ldelta  L increment
c          mts(9)  kdelta  K increment
c          mts(10) ipsi0   Tilt angle flag
c          mts(11) in2     Nuclear spin
c          mts(12) jkmn    Minimum jK
c          mts(13) jmmn    Minimum jM
c
c       Uses:
c               ipar.f
c
c**********************************************************************
c
      subroutine lbasix( ixname,bss,mts,lthb,maxb,new,ierr )
c
      use nlsdim
      use eprprm
      use mtsdef
      use stdio
c
      implicit none
      character*30 ixname,fname
      integer lthb,maxb,new,ierr,bss(5,lthb),mts(MXMTS)
c
      integer i,ioerr,ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,iqnrmn,j,kr,
     #        krmx,lth,lr,mr,mrmx,nrow,iparlr,krmn,mrmn
c
      logical fexist
      integer ipar,iroot
      external ipar,iroot
c
c######################################################################
c
      ierr=0
      if (new.eq.0) then
c
c----------------------------------------------------------------------
c     Look for the basis set indices in the specified file
c----------------------------------------------------------------------
c
         fname=ixname
         lth=iroot(fname)
         fname=fname(:lth)//'.ind'
         lth=lth+4
         inquire(file=fname(:lth),exist=fexist)
         if (.not.fexist) then
c                                           *** file not found
            write(luout,1003) fname(:lth)
            if (luout.ne.luttyo) write (luttyo,1003) fname(:lth)
            ierr=-1
            return
         end if
c
c----------------------------------------------------------------------
c     Read basis set information from index file
c----------------------------------------------------------------------
c
         open (unit=ludisk,file=fname(:lth),status='old',
     #         access='sequential',form='unformatted',iostat=ioerr)
         read (ludisk,iostat=ioerr) lthb,(mts(i),i=1,NTRC)
         if (ioerr.ne.0 .or. lthb.gt.maxb) then
            if (ioerr.ne.0) then
c                                          *** error reading file
               write(luout,1002) fname(:lth)
               if (luout.ne.luttyo) write(luttyo,1002) fname(:lth)
            else
c                                           *** insufficient room 
               write(luout,1001) fname(:lth)
               if (luout.ne.luttyo) write(luttyo,1001) fname(:lth)
            end if
            ierr=-1
            return
         else
            read (ludisk,iostat=ioerr) ((bss(i,j),i=1,5),j=1,lthb)
         end if
c                                  *** error reading file
         if (ioerr.ne.0) then
            write(luout,1002) fname(:lth)
            if (luout.ne.luttyo) write(luttyo,1002) fname(:lth)
            ierr=-1
            return
         end if
c
         close(ludisk)
c
         return
      end if
c
c----------------------------------------------------------------------
c     Routine was called with blank filename or "new" flag nonzero: 
c     Construct a basis using specified MTS. First check the MTS
c----------------------------------------------------------------------
c
      if (mts(NLEMX).le.0 .or. mts(NLOMX).lt.0 .or.
     #    mts(NKMX).lt.0 .or.  mts(NMMX).lt.0 .or.
     #    mts(NIPNMX).lt.0 .or. mts(NKMN).gt.mts(NKMX) .or.
     #    mts(NMMN).gt.mts(NMMX)) then
          write (luout,1004) (mts(i),i=1,NTRC)
         ierr=-1
         return
      end if
c
c *** loop over lr ***
c
      nrow=0
      do 100 lr=0,mts(NLEMX),mts(NLDEL)
        iparlr=ipar(lr)
        if((iparlr.ne.1).and.(lr.gt.mts(NLOMX))) go to 100
c
c   *** loop over kr *** 
c
        krmx=min(mts(NKMX),lr)
        krmn=max(mts(NKMN),-lr)
        do 110 kr=krmn,krmx,mts(NKDEL)
          if((kr.eq.0).and.(iparlr.lt.mts(NJKMN))) go to 110
c
c        *** loop over mr ***
c
          mrmx=min(mts(NMMX),lr)
          mrmn=max(mts(NMMN),-lr)
          do 120 mr=mrmn,mrmx
            if(mts(NIPSI0).eq.0) then
              ipnrmn=mr
              ipnrmx=mr
            else
              ipnrmx=min(mts(NIN2),mts(NIPNMX))
              if(mr.eq.0.and.mts(NJMMN).eq.1) then
                ipnrmn=0
              else
                ipnrmn=-ipnrmx
              end if
            end if
c
c       *** loop over ipnr ***
c
            do 130 ipnr=ipnrmn,ipnrmx
              if((mr.eq.0).and.(ipnr.eq.0).and.(iparlr.lt.mts(NJMMN))) 
     #             go to 130
c
c         *** loop over iqnr ***
c
              iqnrmx=mts(NIN2)-iabs(ipnr)
              iqnrmn=-iqnrmx
              do 140 iqnr=iqnrmn,iqnrmx,2
c
                 nrow=nrow+1
                 if (nrow.le.maxb) then
                    bss(1,nrow)=lr
                    bss(2,nrow)=kr
                    bss(3,nrow)=mr
                    bss(4,nrow)=ipnr
                    bss(5,nrow)=iqnr
                 end if
c
 140          continue
 130        continue
 120      continue
 110    continue
 100  continue
c
      lthb=nrow
      return
c
 1001 format('*** Insufficient room for basis set ''',a,''' ***')
 1002 format('*** Error reading basis set ''',a,''' ***')
 1003 format('*** file ''',a,''' not found ***')
 1004 format('*** lbasix called with illegal MTS: (',6(i3,','),
     #       i3,') ***')
      end
