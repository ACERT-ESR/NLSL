c     Version 1.5  5/2/94
c----------------------------------------------------------------------
c                   =========================
c                      common block FNAMES
c                   =========================
c
c Common block for all the filenames that may possible be associated
c with a given file i.d.  The filenames are defined as follows:
c
c prname : <fileid>.PAR   Parameter file (unformatted)
c fmname : <fileid>.FMT   Parameter listing (formatted)
c bsname : <fileid>.BSS   Basis set datafile
c spname : <fileid>.SPC   Calculated spectrum
c rlname : <fileid>.RLF   Data from swept-field CG calculation
c tdname : <fileid>.TDL   Tridiagonal SLE matrix (unformatted)
c tfname : <fileid>.TDF   Tridiagonal SLE matrix (formatted)
c tsname : <fileid>.MTS   MTS database file
c mtname : <fileid>.MTX   SLE matrix (unformatted)
c mfname : <fileid>.MTF   SLE matrix (formatted)
c egname : <fileid>.EGV   SLE eigenvalues/weights
c ixname : <fileid>.IND   basis set index list (unformatted)
c vtname : <fileid>.STV   Starting vector
c
c NOTE: file 'stddim.inc' must be included before this file.
c----------------------------------------------------------------------
c
      integer MXLTH,MXFLGS
      parameter(MXLTH=30,MXFLGS=10)
c
      character*30 fileid,prname,fmname,bsname,spname,rlname,tdname,
     #             tfname,mtname,mfname,vtname,tsname,egname,ixname
      logical nameOK
      integer namlth
c
      common /fnames/ nameOK,namlth,fileid,prname,fmname,bsname,
     #                spname,rlname,tdname,tfname,mtname,mfname,
     #                vtname,tsname,egname,ixname
c
      character files*30,flags*10
      common /filcom/ files(MXCALC),flags(MXCALC)
c
      logical rdflg,stmflg,stvflg,tdlflg,savflg,angflg
      common /flgcom/ rdflg,stmflg,stvflg,tdlflg,savflg,angflg
