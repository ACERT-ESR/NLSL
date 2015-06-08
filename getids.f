c  Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                       subroutine GETIDS
c                    =========================
c
c  Returns a list of filenames in array <files> in common /filcom/
c  Optionally uses UNIX command line parameter convention
c
c  Includes:
c     stddim.inc
c     stdio.inc
c     fnames.inc
c
c  Written by DEB and SHL.
c----------------------------------------------------------------------
      subroutine getids(n,limit)
      integer n,limit
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'fnames.inc'
c
      integer nargs
      character deflgs*10,arg*30
c
c...............................................
c Code between dotted lines is intended for 
c running these programs in a UNIX environment
c
      integer i,iargc
c      external iargc,getarg
c
c
c  set the initial defaults for each flag
c
      rdflg=.false.
      stmflg=.false.
      stvflg=.false.
      tdlflg=.false.
      savflg=.false.
c
      deflgs=' '
      nargs=iargc()
      n=0
c
      if (nargs.gt.0) then
         do 5 i=1,nargs
            call getarg(i,arg)
            if (arg(1:1) .eq. '-' .or. arg(1:1).eq.'+') then
               deflgs=arg
            else
               n=n+1
               if (n.gt.limit) return
               files(n)=arg
               flags(n)=deflgs
            endif
  5       continue
      else
c..................................................
c
 15     write (luttyo,1020)
        read (luttyi,1030,err=15) arg
        if (arg.eq.' ') return
        if (arg(1:1) .eq. '-' .or. arg(1:1) .eq.'+') then
           deflgs=arg
        else
           n=n+1
           if (n.gt.limit) return
           files(n)=arg
           flags(n)=deflgs
        endif
        goto 15
c..................................................
      end if
c..................................................
c
      return
c
c######################################################################
c  Format statements
c######################################################################
c
 1020 format(2x,'Please enter file identifier (<CR> to end): ',$)
 1030 format(a)
      end
