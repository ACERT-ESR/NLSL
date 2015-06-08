c   Version 1.5 5/2/94
c----------------------------------------------------------------------
c                    =====================
c                     subroutine SETFLG
c                    =====================
c
c      This subroutine interprets <flag> string and sets various flags
c      The flags active are :
c
c		-d : reset to the default flags, which calculates the
c		     SLE matrix elements and starting vector elements
c		     and does not store these into a file.
c		-r : read the matrix and vector elements from old files
c		     rather than calculate again.  If the file doesn't
c		     exist, it calculates those elements.
c		-m : store the matrix elements into file with extension
c		     <.mtx>.
c		-v : store the vector elements into file with extension
c		     <.stv>.
c		-t : store the tridiagonal matrix into file with extension
c		     <.tdl>.
c               -s : save the spectra calculated for individual orientations
c                    in a MOMD calculation, in files <name>.001, <name.002> ...
c               -a : distribute the MOMD orientations using constant angle step
c                    (default is constant cos(angle) )
c
c  Includes
c     stddim.inc   Definition of mxcalc
c     fnames.inc   Definition of names in common /fnames/
c
c----------------------------------------------------------------------
      subroutine setflg( flgstr )
      implicit none
c
      include 'stddim.inc'
      include 'fnames.inc'
c
      integer i
      character*5 flgstr*5,chr*1
      logical value
c----------------------------------------------------------------------
c
      value = (flgstr(1:1).eq.'+')
c
      do 10 i=2,mxflgs
         chr=flgstr(i:i)
         if (chr.eq.' ') return
         if ((chr.eq.'d').or.(chr.eq.'D')) then
            rdflg=.false.
            stmflg=.false.
            stvflg=.false.
            tdlflg=.false.
            savflg=.false.
            angflg=.false.
         end if
         if ((chr.eq.'r').or.(chr.eq.'R')) rdflg=value
         if ((chr.eq.'m').or.(chr.eq.'M')) stmflg=value
         if ((chr.eq.'v').or.(chr.eq.'V')) stvflg=value
         if ((chr.eq.'t').or.(chr.eq.'T')) tdlflg=value
         if ((chr.eq.'s').or.(chr.eq.'S')) savflg=value
         if ((chr.eq.'a').or.(chr.eq.'A')) angflg=value
 10   continue
c
      return
      end
