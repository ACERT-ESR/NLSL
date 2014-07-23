c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine IPFIND
c                    =========================
c
c  Search for the given token in a list of the names of EPRESA spectral
c  calculation parameters, which is defined in the BLOCK DATA segment 
c  coded below. iprtn returns a value as follows:
c
c -200 < iprtn < -100   Specified name was the (axial tensor) alias 
c                        of a floating point parameter. 
c                        Returned value is the negative of the index
c                        of the parameter in the fepr array, minus 100. 
c
c -100 < iprtn < 0      Specified name was the (spherical tensor) alias 
c                        of a floating point parameter. 
c                        Returned value is the negative of the index of the
c                        parameter in the fepr array. 
c
c    0 < iprtn < 100    Specified name was a floating point parameter.
c                        Returned value is the index of the parameter
c                        in the fepr array
c
c        iprtn > 100    Specified name was an integer parameter.
c                        Returned value minus 100 is the index of the parameter
c                        in the iepr array
c
c
c   NOTE: The order of names in this list MUST correspond to the
c         order of parameters in EPRPRM.INC for this routine to work 
c         properly.
c 
c varspec is set to 1 if the selected parameter requires a new matrix,
c  set to 2 if new Eiginvalues and Eigenvectors are required for a change in
c  the specified variable and set to 3 for variables requiring new
c  xspec, set to 4 in the case of gib or lib requiring 
c  only a new time domain spec calculated from xspec, and 6 if nothing 
c  is required.  
c
c   Includes
c      nlsdim.inc
c      eprprm.inc
c      parcom.inc
c      lpnam.inc
c
c----------------------------------------------------------------------
      subroutine ipfind(iprtn,token,lth,varspec)
      implicit none
      integer iprtn,lth,varspec
c
      include 'limits.inc'
      include 'parms.inc'
      include 'simparm.inc'
      include 'lpnam.inc'
      character token*(WORDLG)
c
      integer i
c
c----------------------------------------------------------------------
c     Search the list of floating point parameter names
c----------------------------------------------------------------------
      do 10 i=1,NFPRM
        if (token(:lth).eq.parnam(i)(:lth)) then
          iprtn=i
          varspec=specid(i)
          return
        end if
   10   continue
c
c----------------------------------------------------------------------
c     Search the list of floating point parameter aliases 
c     (these are all names of spherical tensor components)
c     Return negative index if found in alias1
c     Return negative index-100 if found in alias2
c----------------------------------------------------------------------
      do 11 i=1,nalias
        if (token(:lth).eq.alias1(i)(:lth)) then
           iprtn=1-(IGXX+i)
           varspec=1
           return
        else if (token(:lth).eq.alias2(i)(:lth)) then
           iprtn=-(99+IGXX+i)
           varspec=1
           return
        end if 
  11    continue
c
c----------------------------------------------------------------------
c     Search the list of integer parameter names
c----------------------------------------------------------------------
      do 12 i=1,NIPRM
        if (token(:lth).eq.iprnam(i)(:lth)) then
          iprtn=100+i
          return
        end if
   12   continue
c
c----------------------------------------------------------------------
c     Token was not found
c----------------------------------------------------------------------
      iprtn=0
      varspec=-1
      return
      end


c----------------------------------------------------------------------
c                    =========================
c                       function ISFIND
c                    =========================
c
c  Search for the given token in a list of the symbolic names defined
c  in the BLOCK DATA segment coded below. If a match is found, ISFIND
c  returns an index into the symval array which contains the values 
c  symbolized by each element of the symbol array. Otherwise, -1 is
c  returned.
c----------------------------------------------------------------------
      function isfind(token,lth)
      implicit none
      integer isfind,lth
      character token*30
c
      include 'limits.inc'
      include 'lpnam.inc'
c
      integer i
c
c----------------------------------------------------------------------
c     Search the list of symbolic names
c----------------------------------------------------------------------
      do 10 i=1,NSYMBL
        if (token(:lth).eq.symbol(i)(:lth)) then
          isfind=i
          return
        end if
   10   continue
        isfind=-1
        return
        end


c----------------------------------------------------------------------
c                    =========================
c                       block data LPINIT
c                    =========================
c
c     Initialization of the strings used for pre-defined parameter
c     names for the EPRP family of programs.
c
c     NOTE: The order of names in this list MUST correspond to the
c      order of parameters in EPRPRM.INC for this routine to work 
c      properly.
c----------------------------------------------------------------------
      block data
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'lpnam.inc'
c
c----------------------------------------------------------------------
c
c ########## Definition of names of EPRESA parameters ###############
c see lpnam.inc
c
      data parnam/'GIB',   'LIB',   'HWID',  'MWID',  'GXX',
     #            'GYY',   'GZZ',   'AXX',   'AYY',   'AZZ',
     #            'RX',    'RY',    'RZ',    'PML',   'PMXY',
     #            'PMZZ',  'DJF',   'DJFPRP','OSS',   'PSI',
     #            'ALD',   'BED',   'GAD',   'ALM',   'BEM',
     #            'GAM',   'C20',   'C22',   'C40',   'C42',
     #            'C44',   'T2EDI', 'T2NDI', 'T2EFI', 'T1EDI',
     #            'T1NDI', 'SITEWG',
     #            'B0',    'GAMMAN','SHFT',  'FIELDI',
     #            'FIELDF','CGTOL', 'SHIFTR','SHIFTI','PTOL',
     #           'INIT1', 'STEPT1','INIT2', 'STEPT2','TFIX','WEIGHT'/
c
      data iprnam/'IN2',   'IPDF',  'IST',   'ITYPE', 'NORT',
     #            'IGSPH', 'IASPH', 'IRSPH', 'ML',    'MXY',
     #            'MZZ',   'LEMX',  'LOMX',  'KMX',   'MMX', 
     #            'IPNMX', 'JKMN',  'JMMN',  'NSTEP', 'NFIELD',
     #            'IEXP',  'ICOMB', 'SPTST', 'DUM1',  'DUM2' / 
c
      data alias1/'G1',    'G2',    'G3',  
     #            'A1',    'A2',    'A3', 
     #            'RBAR',  'N',     'NXY' /
c
      data alias2/'GPRP',  'GRHM',  'GPLL',
     #            'APRP',  'ARHM',  'APLL',
     #            'RPRP',  'RRHM',  'RPLL'/
c
c specid: 0=new matrix required, 1=use evals,evecs, 2=use xspec, 3=have spec,
c 52 values:
c
      data specid/2,2,1,1,0,	!gib,lib,hwid,mwid,gxx
     #	0,0,0,0,0,		!gyy,gzz,axx,ayy,azz
     #  0,0,0,0,0,		!rx,ry,rz,pml,pmxy
     #  0,0,0,0,0,		!pmzz,djf,djfprp,oss,psi
     #  0,0,0,0,0,		!ald,bed,gad,alm,bem
     #  0,0,0,0,0,		!gam,c20,c22,c40,c42
     #	0,0,0,0,0,		!c44,t2edi,t2ndi,t2efi,t1edi
     #  0,4,			!t1ndi,sw1 - site weight
     #  0,3,1,3,		!b0,gamman,shft,fieldi
     #  3,0,0,0,3,		!fieldf,cgtol,shiftr,shifti,ptol
     #  1,1,1,1,1,1/		!init1,stept1,init2,stept2,tfix,weight
c
      data symstr/'CARTESIAN','SPHERICAL','AXIAL'/
c
      data symbol/'BROWNIAN','NONBROWNIA','ANISOVISCO','FREE','JUMP'/
      data symval/ 0,         1,           2,           1,     0     /
c
      end
