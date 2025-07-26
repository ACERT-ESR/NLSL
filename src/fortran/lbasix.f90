!  Version 1.4 (NLS)
!**********************************************************************
!                       ===================
!                       SUBROUTINE : LBASIX
!                       ===================
!
!       This routine builds a list of the basis set indices in common 
!       block /indexl/ for the truncation parameters set in /eprdat/ 
!       or reads them from the index file "<name>.ind" which contains 
!       a pruned basis set produced by the EPRBL program. 
!
!       If one of the following conditions is true:
!         (1) The "new" flag is set nonzero on input
!         (2) the input filename is blank
!
!       this routine routine builds the full basis set indicated
!       by the truncation parameters specified in the MTS array.
!       These are as follows:
!          mts(1)  lemx    Maximum even L
!          mts(2)  lomx    Maximum odd L
!          mts(3)  kmn     Minimum K (K<0 => jK = -1)
!          mts(4)  kmx     Maximum K
!          mts(5)  mmn     Minimum M (M<0 => jM = -1)
!          mts(6)  mmx     Maximum M
!          mts(7)  ipnmx   Maximum pI
!          mts(8)  ldelta  L increment
!          mts(9)  kdelta  K increment
!          mts(10) ipsi0   Tilt angle flag
!          mts(11) in2     Nuclear spin
!          mts(12) jkmn    Minimum jK
!          mts(13) jmmn    Minimum jM
!
!       Uses:
!               ipar.f
!
!**********************************************************************
!
      subroutine lbasix( ixname,bss,mts,lthb,maxb,new,ierr )
!
      use nlsdim
      use eprprm
      use mtsdef
      use stdio
!
      implicit none
      character*30 ixname,fname
      integer lthb,maxb,new,ierr,bss(5,lthb),mts(MXMTS)
!
      integer i,ioerr,ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,iqnrmn,j,kr,
     #        krmx,lth,lr,mr,mrmx,nrow,iparlr,krmn,mrmn
!
      logical fexist
      integer ipar,iroot
      external ipar,iroot
!
!######################################################################
!
      ierr=0
      if (new.eq.0) then
!
!----------------------------------------------------------------------
!     Look for the basis set indices in the specified file
!----------------------------------------------------------------------
!
         fname=ixname
         lth=iroot(fname)
         fname=fname(:lth)//'.ind'
         lth=lth+4
         inquire(file=fname(:lth),exist=fexist)
         if (.not.fexist) then
!                                           *** file not found
            write(luout,1003) fname(:lth)
            if (luout.ne.luttyo) write (luttyo,1003) fname(:lth)
            ierr=-1
            return
         end if
!
!----------------------------------------------------------------------
!     Read basis set information from index file
!----------------------------------------------------------------------
!
         open (unit=ludisk,file=fname(:lth),status='old',
     #         access='sequential',form='unformatted',iostat=ioerr)
         read (ludisk,iostat=ioerr) lthb,(mts(i),i=1,NTRC)
         if (ioerr.ne.0 .or. lthb.gt.maxb) then
            if (ioerr.ne.0) then
!                                          *** error reading file
               write(luout,1002) fname(:lth)
               if (luout.ne.luttyo) write(luttyo,1002) fname(:lth)
            else
!                                           *** insufficient room 
               write(luout,1001) fname(:lth)
               if (luout.ne.luttyo) write(luttyo,1001) fname(:lth)
            end if
            ierr=-1
            return
         else
            read (ludisk,iostat=ioerr) ((bss(i,j),i=1,5),j=1,lthb)
         end if
!                                  *** error reading file
         if (ioerr.ne.0) then
            write(luout,1002) fname(:lth)
            if (luout.ne.luttyo) write(luttyo,1002) fname(:lth)
            ierr=-1
            return
         end if
!
         close(ludisk)
!
         return
      end if
!
!----------------------------------------------------------------------
!     Routine was called with blank filename or "new" flag nonzero: 
!     Construct a basis using specified MTS. First check the MTS
!----------------------------------------------------------------------
!
      if (mts(NLEMX).le.0 .or. mts(NLOMX).lt.0 .or.
     #    mts(NKMX).lt.0 .or.  mts(NMMX).lt.0 .or.
     #    mts(NIPNMX).lt.0 .or. mts(NKMN).gt.mts(NKMX) .or.
     #    mts(NMMN).gt.mts(NMMX)) then
          write (luout,1004) (mts(i),i=1,NTRC)
         ierr=-1
         return
      end if
!
! *** loop over lr ***
!
      nrow=0
      do 100 lr=0,mts(NLEMX),mts(NLDEL)
        iparlr=ipar(lr)
        if((iparlr.ne.1).and.(lr.gt.mts(NLOMX))) go to 100
!
!   *** loop over kr *** 
!
        krmx=min(mts(NKMX),lr)
        krmn=max(mts(NKMN),-lr)
        do 110 kr=krmn,krmx,mts(NKDEL)
          if((kr.eq.0).and.(iparlr.lt.mts(NJKMN))) go to 110
!
!        *** loop over mr ***
!
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
!
!       *** loop over ipnr ***
!
            do 130 ipnr=ipnrmn,ipnrmx
              if((mr.eq.0).and.(ipnr.eq.0).and.(iparlr.lt.mts(NJMMN))) 
     #             go to 130
!
!         *** loop over iqnr ***
!
              iqnrmx=mts(NIN2)-iabs(ipnr)
              iqnrmn=-iqnrmx
              do 140 iqnr=iqnrmn,iqnrmx,2
!
                 nrow=nrow+1
                 if (nrow.le.maxb) then
                    bss(1,nrow)=lr
                    bss(2,nrow)=kr
                    bss(3,nrow)=mr
                    bss(4,nrow)=ipnr
                    bss(5,nrow)=iqnr
                 end if
!
 140          continue
 130        continue
 120      continue
 110    continue
 100  continue
!
      lthb=nrow
      return
!
 1001 format('*** Insufficient room for basis set ''',a,''' ***')
 1002 format('*** Error reading basis set ''',a,''' ***')
 1003 format('*** file ''',a,''' not found ***')
 1004 format('*** lbasix called with illegal MTS: (',6(i3,','),
     #       i3,') ***')
      end
