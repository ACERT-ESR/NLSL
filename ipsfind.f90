c NLSL Version 1.9.0 beta 2/7/15
c----------------------------------------------------------------------
c                    =========================
c                         function IPFIND
c                    =========================
c
c  Search for the given token in one of the lists of the names of EPRLL
c  spectral calculation parameters, which are defined in module lpnam. 
c  Function return values are as follows:
c
c -200 < ipfind < -100   Specified name was the (axial tensor) alias 
c                        of a floating point parameter. 
c                        Returned value is the negative of the index
c                        of the parameter in the fepr array, minus 100. 
c
c -100 < ipfind < 0      Specified name was the (spherical tensor) alias 
c                        of a floating point parameter. 
c                        Returned value is the negative of the index of
c                        the parameter in the fepr array. 
c
c    0 < ipfind < 100    Specified name was a floating point parameter.
c                        Returned value is the index of the parameter
c                        in the fepr array
c
c        ipfind > 100    Specified name was an integer parameter.
c                        Returned value minus 100 is the index of the
c                        parameter in the iepr array
c                       
c
c  NOTE: The order of names in this list MUST correspond to the order
c   of parameters in modules eprprm and lpnam for this routine to work 
c   properly.
c 
c----------------------------------------------------------------------
      function ipfind(token,lth)

      use nlsdim
      use parcom
      use eprprm
      use lpnam

      implicit none
      integer ipfind,lth
      character token*30
c
      integer i
c
c----------------------------------------------------------------------
c     Search the list of floating point parameter names
c----------------------------------------------------------------------
      do i=1,nfprm
        if (token(:lth).eq.parnam(i)(:lth)) then
          ipfind=i
          return
        end if
      end do
c
c----------------------------------------------------------------------
c     Search both lists of floating point parameter aliases 
c     (these are all names of spherical tensor components)
c     Return negative index if found in alias1
c     Return negative index-100 if found in alias2
c----------------------------------------------------------------------
      do i=1,nalias
        if (token(:lth).eq.alias1(i)(:lth)) then
           ipfind=1-(IWXX+i)
           return
        else if (token(:lth).eq.alias2(i)(:lth)) then
           ipfind=-(99+IWXX+i)
           return
        end if 
      end do
c
c----------------------------------------------------------------------
c     Search the list of integer parameter names
c----------------------------------------------------------------------
      do i=1,niprm
        if (token(:lth).eq.iprnam(i)(:lth)) then
          ipfind=100+i
          return
        end if
      end do
c
c----------------------------------------------------------------------
c     Token was not found
c----------------------------------------------------------------------
      ipfind=0
      return
      end function ipfind

c----------------------------------------------------------------------
c                    =========================
c                        function ISFIND
c                    =========================
c
c  Search for the given token in one of several lists of strings.
c  Potential matches are stored in the arrays that hold the datafile,
c  basisid, and symbolic names. Datafile and basisid names are defined
c  via I/O within the datac and basisc routines. Symbolic names are
c  defined in the initialization segment of module lpnam--formerly
c  these were encoded in the block data section of nlstxt.f.
c
c  If a match is found, ISFIND returns an index into the array which 
c  contains the match. The index is positive if the match is in the
c  symbol array. The index is multipled by -1 if the match is in the
c  datafile or basisid array. Otherwise, 0 is returned.
c----------------------------------------------------------------------
      function isfind(token,lth)
c
      use nlsdim
      use expdat
      use basis
      use lpnam
      use stdio
      implicit none
c
      integer isfind,lth
      character token*30,tmpstr*30
c
      integer i,nlth
c
c----------------------------------------------------------------------
c     Search the list of datafile names
c----------------------------------------------------------------------
      do i=1,nspc
         tmpstr=dataid(i)
         nlth=itrim(tmpstr)
         if (nlth.ge.lth) then
            if (token(:lth).eq.tmpstr(:lth)) then
               isfind=-i
               return
            end if
         end if
      end do
c
c----------------------------------------------------------------------
c     Search the list of basis set names
c----------------------------------------------------------------------
      do i=1,nbas
         tmpstr=basisid(i)
         nlth=itrim(tmpstr)
         if (nlth.ge.lth) then
            if (token(:lth).eq.tmpstr(:lth)) then
               isfind=-i
               return
            end if
         end if
      end do
c
c
c----------------------------------------------------------------------
c     Search the list of symbolic names
c----------------------------------------------------------------------
      tmpstr=token
      call touppr(tmpstr,lth)
      do i=1,nsymbl
         if (tmpstr(:lth).eq.symbol(i)(:lth)) then
            isfind=i
            return
         end if
      end do
c
      isfind=0
      return
      end
