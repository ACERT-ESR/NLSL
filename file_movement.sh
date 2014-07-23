### marks things that I should LOOK AT LATER
# "very similar" MEANS THAT THERE A MOSTLY STYLISTIC / FORTRAN-VERSION TYPE CHANGES
# "modified extensively" are files that are called the same thing, but likely do different things now
# ccrints.f is just a duplicate file
# qrsolv.f is just a duplicate file
# covar.f is just a duplicate file
# plgndr.f is just a duplicate file
# scmvm.f is just a duplicate file
# anxlk.f is just a duplicate file
# w3j.f is just a duplicate file
# ipar.f is just a duplicate file
# lmpar.f is just a duplicate file
# qrfac.f is just a duplicate file
# helpc.f is just a duplicate file
# bessel.f is just a duplicate file
# setnm.f is just a duplicate file
# zdotu2 is very similar to zdotu, so:
git mv zdotu.f zdotu2.f
# ordrpr.f is just a duplicate file
# strutl.f is just a duplicate file
# fz.f is just a duplicate file
# lmnls.f is just a duplicate file
# cd2km.f is just a duplicate file
# enorm.f is just a duplicate file
# zaxpy2 is very similar to zaxpy so:
git mv zaxpy.f zaxpy2.f
# tensym.f is just a duplicate file
# zaypx.f is just a duplicate file
# nlspmc is different, but clearly derived from nlsl:
git mv nlsl.f nlspmc.f
# xshft.f is basically the same as ftfuns.f, but with more documentation:
git mv ftfuns.f xshft.f
# ipfind.f has apparently been modified
### has some block data that might be better off in an include file
# brent.f has been modified
### some commentary got stripped out
# matrxo.f and matrxd.f are apparently similar to each other
# matrxo.f and matrxd.f are apparently similar to each other
# getdat.f has been modified extensively
# therefore, afterwards, I will need to run:
# mv getdat.f getd2d.f # before checkin of movements
# grep *.f makefile -Rile "\<getdat\>" | xargs -i@ sed -e "s/\<getdat\>/getd2d/g" -i @
# after movement checkin
# cscg.f has been modified
# comps.f and spectra.f are apparently similar to each other
# comps.f and spectra.f are apparently similar to each other
# stveco.f and stvect.f are not similar enough to warrant a rename
# pcheck.f and lcheck.f are similar, and the description is similar enough:
git mv lcheck.f pcheck.f
# datac.f has been modified extensively
### this likely will require some manual merging down the line
# addprm.f has been modified 
### this likely will require some manual merging down the line
# srchc.f has been modified  -- hard to tell what exactly is up here
# rmvprm.f has been modified a bit but should be doing roughly the same
# letcmc.f appears to be derived from letc.f
git mv letc.f letcmc.f
# dpmpar.f and enorm.f are not related
# mnbrak.f has been modified but seems to be doing the same thing
# p1pfun.f and mnbrak.f are not related
# convtc.f has been modified -- clearly, it's derived from the original file, but I'm not sure whether or not it's doing the same thing

#now the .inc files
# rndoff.inc is just a duplicate file
# maxl.inc is just a duplicate file
# parmequ.inc seems to supersede prmeqv.inc
git mv prmeqv.inc parmequ.inc
#eprmat.inc is just a duplicate file
#stdio.inc is just a duplicate file
#physcn.inc is just a duplicate file
# names.inc seems to be an upgrade of nlsnam.inc
git mv nlsnam.inc names.inc
# parms seems to be an upgrade of parcom
git mv parcom.inc parms.inc
# pidef.inc is just a duplicate file
# lmcomm.inc is an upgrade of lmcom.inc
git mv lmcom.inc lmcomm.inc
# lpnam has been significantly modified, but is the same
# simparm.inc is a modified version of eprprm.inc
git mv eprprm.inc simparm.inc
# a bunch of the header files within nlspmcL appear to have duplicate info -- ignore
# expdat.inc (old) and datas.inc (new) both define the expdat common block
git mv expdat.inc datas.inc
