c  NLSPMC VERSION 1.0 2/5/99
c*********************************************************************
c
c	declarations of parameters for new series and component fits.
c
c	Notes:
c	  1) The parameters define here are used as
c
c	    nspectra - specify number of spectra expected (SPECTRA)
c	    ncomps - specify number of components in each spectra (COMP)
c           basinfo - 1=basis set index, 2=have matrix and tridiag,
c			3=have spectra, 4=undefined yet
c           nbasis - current number of basis sets (max=MXSPEC*MXSITE)
c           ndimoss - off diagonal dimension
c           ndimdss - diagonal dimension
c           cursite - current site under consideration in simulation
c                     get from basinfo(1)
c c           setbass - got this basis set
c
c To use the basis set location information, use basinfo(index,spec,site)
c  to identify which basis (index=1) the given spec and site require.
c With this information, ndimoss and ndimdss can be used in conjunction with
c pidptr and dpidptr, in basis.inc, to locate the specific basis information.
c
c      logical setbass
c
      
      integer nspectra,ncomps,basinfo,nbasis,
     #		ndimoss,ndimdss,cursite
c
c initflg is set true after program starts and user enters type of fit
c   information
c
      common /mulcomp/ nspectra,ncomps,basinfo(5,MXSPEC,MXSITE),
     #	nbasis,ndimoss(MXSPEC*MXSITE),ndimdss(MXSPEC*MXSITE)

c     #		setbass(MXSPEC*MXSITE)

      
