########################################################################
#
#	     ====================================================
#	     Make descriptor file for NLSPC and related programs
#	     ====================================================
#
#  VERSION FOR RUNNING UNDER windows 2000
#
########################################################################
#
.SUFFIXES: .f .inc .o .c

# See: https://computing.llnl.gov/tutorials/openMP/
#F77=gfortran -std=legacy
F77=gfortran
#FFLAGS = -c -O2 -g
UNAME := $(shell uname)
$(info $$UNAME is [${UNAME}])
ifneq (,$(findstring MINGW32_NT,$(UNAME)))
OS := "windows"
$(info in windows)
FFLAGS = -c -O2 -g -fopenmp -std=gnu
else
OS := "linux"
$(info not in windows)
FFLAGS = -c -O2 -g -fopenmp -std=gnu -m64 -mcmodel=medium
endif
LIB = -L /usr/X11/lib64 -lX11 -L/usr/X11R6/lib64 

LIBS2  = C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsl.lib C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsls_err.lib

#	
#-----------------------------------------------------------------------
#			Definitions
#-----------------------------------------------------------------------
#
NLSPI = limits.inc names.inc parms.inc stdio.inc miscel.inc 
NLSPI2 = parmequ.inc simparm.inc
CMDI = $(NLSPI) $(NLSPI2) lpnam.inc names.inc basis.inc datas.inc 
CMDI2 = $(CMDI) lmcomm.inc parms.inc rndoff.inc basis.inc
VCH = vchange.o

LETI = limits.inc simparm.inc datas.inc parms.inc lpnam.inc
NLSI = limits.inc simparm.inc datas.inc lmcomm.inc stdio.inc parms.inc
FUNI = limits.inc simparm.inc datas.inc parms.inc 
MATI = limits.inc simparm.inc basis.inc eprmat.inc stvcom.inc
# CMDI = $(NLSI) lmcomm.inc parms.inc lpnam.inc
DATI = limits.inc datas.inc names.inc stdio.inc parms.inc 
CGO  = cscg.o zaypx.o zdotu2.o zaxpy2.o scmvm.o znormu.o
EVCG = csval.o compev.o lump.o comgap.o isoev.o inverr.o cmtqli.o csvec.o
EPO1 = sim2d.o evcgf.o cd2km.o anxlk.o w3j.o ipar.o fbasis.o
EPO2 = matrxo.o matrxd.o stveco.o stvecd.o ccrints.o fz.o bessel.o plgndr.o
EPO3 = spcalc.o convft.o fft.o switch.o zgemm.o zgemv.o lsame.o xerbla.o
NEWSTF = spectra.o comps.o vchange.o
NLSC = cmds.o datac.o letcmc.o addprm.o rmvprm.o srchc.o tensym.o convtc.o
NLSS = strutl.o lprmpt.o helpc.o ipfind.o
NLSD = setnm.o getd2d.o wrfit.o
NLSH = xshft.o sscale.o mnbrak.o brent.o p1pfun.o
NLSF = fitp.o pfunnew.o pcheck.o $(NLSH) $(EPO1) $(EPO2) $(EPO3) $(EVCG) $(CGO)
NLSB = lmnls.o enorm.o dpmpar.o qrfac.o lmpar.o qrsolv.o covar.o mapxxx.o
NLSO = nlspmc.o $(NLSC) $(NLSS) $(NLSD) $(NLSF) $(NLSB) nlsinit.o ordrpr.o

#-----------------------------------------------------------------------
#		Object files
#-----------------------------------------------------------------------

ifeq ($(OS),"windows")
all : nlspmc.exe
clean:
	echo "running clean for windows"
	$(RM) *.o *.mod nlspmc.exe
else
all             : nlspmc
clean           :
	echo "running clean for linux"
	$(RM) *.o *.mod nlspmc_cmd
endif

addprm.o	: addprm.f $(CMDI)
anxlk.o		: anxlk.f limits.inc simparm.inc rndoff.inc
bessel.o	: bessel.f rndoff.inc pidef.inc
brent.o		: brent.f limits.inc stdio.inc parms.inc
cmds.o		: cmds.f $(CMDI2) 
cd2km.o		: cd2km.f rndoff.inc pidef.inc
cmtqli.o	: cmtqli.f limits.inc rndoff.inc
comgap.o	: comgap.f limits.inc
compev.o	: compev.f limits.inc simparm.inc
comps.o		: comps.f $(CMDI2)
convft.o	: convft.f limits.inc simparm.inc wkspcm.inc physcn.inc
convtc.o	: convtc.f limits.inc simparm.inc parms.inc stdio.inc lpnam.inc
cscg.o 		: cscg.f limits.inc parms.inc stdio.inc rndoff.inc
csval.o		: csval.f limits.inc stdio.inc simparm.inc rndoff.inc
csvec.o		: csvec.f limits.inc rndoff.inc
datac.o		: datac.f $(DATI) lmcomm.inc wkspcm.inc
evcgf.o		: evcgf.f limits.inc stdio.inc simparm.inc parms.inc basis.inc wkspcm.inc
fbasis.o	: fbasis.f limits.inc simparm.inc basis.inc stdio.inc
fft.o		: fft.f
fitp.o  	: fitp.f $(NLSI) names.inc tdspec.inc lmtxt.inc parms.inc
ftest.o		: ftest.f
fz.o		: fz.f limits.inc simparm.inc
getd2d.o	: getd2d.f limits.inc stdio.inc
helpc.o		: helpc.f stdio.inc
inverr.o	: inverr.f limits.inc simparm.inc
ipar.o		: ipar.f
ipfind.o	: ipfind.f limits.inc simparm.inc parms.inc lpnam.inc
isoev.o		: isoev.f limits.inc simparm.inc
letcmc.o	: letcmc.f $(VCH) $(LETI) stdio.inc
lmnls.o		: lmnls.f mapxxx.f limits.inc parms.inc parms.inc
lump.o		: lump.f limits.inc simparm.inc
mapxxx.o	: mapxxx.f parms.inc
matrxo.o	: matrxo.f $(MATI) maxl.inc rndoff.inc physcn.inc
matrxd.o	: matrxd.f $(MATI) maxl.inc rndoff.inc physcn.inc
mnbrak.o	: mnbrak.f limits.inc parms.inc
nlsinit.o 	: nlsinit.f $(DATI) parms.inc simparm.inc parmequ.inc lmcomm.inc
nlspmc.o	: nlspmc.f $(NLSPI) $(NLSPI2)
ordrpr.o	: ordrpr.f rndoff.inc
p1pfun.o	: p1pfun.f limits.inc parms.inc datas.inc lmcomm.inc
pcheck.o	: pcheck.f $(NLSI) names.inc maxl.inc rndoff.inc
pfunnew.o	: pfunnew.f $(VCH) $(FUNI) lpnam.inc tdspec.inc wkspcm.inc stdio.inc parms.inc
rmvprm.o	: rmvprm.f limits.inc parms.inc simparm.inc stdio.inc
scmvm.o 	: scmvm.f limits.inc rndoff.inc eprmat.inc
scspec.o	: scspec.f stdio.inc
setnm.o		: setnm.f limits.inc names.inc
sim2d.o		: sim2d.f $(MATI) parmequ.inc parms.inc egvcom.inc stdio.inc
spcalc.o	: spcalc.f limits.inc simparm.inc stvcom.inc wkspcm.inc physcn.inc
spectra.o	: spectra.f $(NLSPI) $(NLSPI2) lpnam.inc
srchc.o		: srchc.f $(VCH) limits.inc simparm.inc parms.inc lpnam.inc stdio.inc
sscale.o	: sscale.f limits.inc datas.inc
strutl.o	: strutl.f stdio.inc
stveco.o 	: stveco.f limits.inc simparm.inc basis.inc stvcom.inc rndoff.inc
stvecd.o 	: stvecd.f limits.inc simparm.inc basis.inc stvcom.inc wkspcm.inc
tensym.o	: tensym.f limits.inc simparm.inc parms.inc lpnam.inc stdio.inc
tstjac.o	: tstjac.f stdio.inc
vchange.o	: vchange.f rndoff.inc
w3j.o		: w3j.f maxl.inc
whris.o		: whris.f
wrfit.o		: wrfit.f limits.inc datas.inc stdio.inc
xshft.o		: xshft.f
zaypx.o 	: zaypx.f limits.inc rndoff.inc
zaxpy2.o 	: zaxpy2.f limits.inc rndoff.inc
znormu.o	: znormu.f limits.inc rndoff.inc
zscsw.o		: zscsw.f limits.inc rndoff.inc
zdotu2.o		: zdotu2.f limits.inc rndoff.inc

lsame.o         : lsame.f
zgemm.o         : zgemm.f
zgemv.o         : zgemv.f
xerbla.o        : xerbla.f

#-----------------------------------------------------------------------
#		Executable files
#-----------------------------------------------------------------------

ifeq ($(OS),"windows")
nlspmc.exe	: $(NLSO) $(NEWSTF) ftest.o
	$(F77) -fopenmp $(LOADFLG) $(NEWSTF) ftest.o -o $@ $(NLSO) $(LIB) $(LIB2) -lc
else
nlspmc	: $(NLSO) $(NEWSTF) ftest.o
	$(F77) -fopenmp $(LOADFLG) $(NEWSTF) ftest.o -o $@ $(NLSO) $(LIB) $(LIB2) -lc
endif


#-----------------------------------------------------------------------
#			Default actions
#-----------------------------------------------------------------------

.c.o   :
	cc $(CFLAGS) $*.c

.f.o   :
	$(F77) $(FFLAGS) $*.f

.inc.o :
	$(F77) $(FFLAGS) $*.f
