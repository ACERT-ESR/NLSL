########################################################################
#
#	     ====================================================
#	     Make descriptor file for NLSL and related programs
#	     ====================================================
#
########################################################################

#Suffixes are now replaced by pattern rules
.PRECIOUS: .f90 .mod .f .inc .h .c
.PHONY: all clean

########################################################################

RM = rm -f

F77 = gfortran
F90 = $(F77)
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
#Note, compiler bug in gfortran<4.8: -mcmodel=medium causes spurious warnings
#Message says, "Warning: ignoring incorrect section type for .lbss"
endif

CC = gcc
CFLAGS = -c -O2 -DADD_UNDERSCORE

FLINK = gfortran
LIBS = -L/usr/X11R6/lib -lX11 -lc

#LIBS = dfor.lib
#LIBS  = C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsl.lib C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsls_err.lib 

#F77 = g77
#FFLAGS = -c -O9 -fomit-frame-pointer -ffast-math -malign-double\
	-funroll-loops -march=i586
#FFLAGS = -c -g

########################################################################

#	
#-----------------------------------------------------------------------
#			Definitions
#-----------------------------------------------------------------------
#

MOD1 = eprprm.o parcom.o nlsdim.o errmsg.o lpnam.o expdat.o basis.o
MOD2 = symdef.o maxl.o bincom.o dfunc.o eprmat.o ftwork.o iterat.o lmcom.o mspctr.o mtsdef.o 
MOD3 = nlsnam.o parsav.o physcn.o pidef.o rnddbl.o stdio.o timer.o tridag.o 
SPCI = nlsdim.mod eprprm.mod expdat.mod parcom.mod tridag.mod basis.mod pidef.mod
BASI = nlsdim.mod eprprm.mod parcom.mod basis.mod stdio.mod
LETI = nlsdim.mod eprprm.mod expdat.mod parcom.mod lpnam.mod
NLSI = nlsdim.mod nlsnam.mod eprprm.mod expdat.mod lmcom.mod parcom.mod symdef.mod stdio.mod parcom.mod
PRMI = nlsdim.mod eprprm.mod lmcom.mod parcom.mod lpnam.mod
FUNI = nlsdim.mod eprprm.mod expdat.mod tridag.mod parcom.mod iterat.mod
FNI2 = basis.mod pidef.mod stdio.mod mspctr.mod ftwork.mod 
FITI = nlsnam.mod eprmat.mod errmsg.mod mspctr.mod iterat.mod
MATI = nlsdim.mod eprprm.mod eprmat.mod rnddbl.mod
CHKI = maxl.mod nlsdim.mod eprprm.mod rnddbl.mod errmsg.mod
DATI = nlsdim.mod expdat.mod mspctr.mod nlsnam.mod stdio.mod parcom.mod 
SCLI = nlsdim.mod expdat.mod lmcom.mod parcom.mod mspctr.mod iterat.mod 
CGO  = cscg.o zaypx.o cgltri.o zdotu.o zaxpy.o scmvm.o tdsqz.o tdchek.o
EPO1 = momdls.o eprls.o matrll.o cd2km.o anxlk.o w3j.o stvect.o cfs.o $(CGO)
EPO2 = pstvec.o pmatrl.o ccrints.o fz.o bessel.o plgndr.o ipar.o
NLSB = lbasix.o setmts.o 
NLSC1 = datac.o fitc.o letc.o addprm.o rmvprm.o statc.o tensym.o writec.o
NLSC2 = convtc.o helpc.o srchc.o setprm.o assgnc.o basisc.o shiftc.o scalec.o
NLSC3 = covrpt.o stats.o confc.o series.o sitec.o varyc.o fixc.o parc.o
NLSC = $(NLSC1) $(NLSC2) $(NLSC3)
NLSD = setnm.o getdat.o genio.o fmomnt.o correl.o ftfuns.o gconvl.o
NLSH = sshift.o sscale.o lgrint.o mnbrak.o brent.o l1pfun.o qrutil.o writr.o
NLSF = fitl.o lfun.o lcheck.o setspc.o $(NLSH) $(EPO1) $(EPO2) ordrpr.o
NLSL = lmnls.o enorm.o qrfac.o lmpar.o qrsolv.o covar.o
NLSN = dchex.o daxpy.o dcopy.o ddot.o drotg.o
NLSP = pltx.o 
NLSS = strutl1.o strutl2.o lprmpt.o catch.o ipsfind.o 
ifeq ($(OS),"windows")
NLSO = nlsl.o pltx_dummy.o $(MOD1) $(MOD2) $(MOD3) $(NLSC) $(NLSS) $(NLSD) $(NLSF) $(NLSL) $(NLSN) $(NLSB) $(NLSW)
else
NLSO = nlsl.o $(MOD1) $(MOD2) $(MOD3) $(NLSC) $(NLSS) $(NLSD) $(NLSF) $(NLSL) $(NLSP) $(NLSN) $(NLSB) $(NLSW)
endif

#-----------------------------------------------------------------------
#		Object files
#-----------------------------------------------------------------------

ifeq ($(OS),"windows")
all : nlsl.exe
clean :
#	del *.o nlsl.exe
	rm *.o nlsl.exe
else
all             : nlsl
clean           :
			rm -f *.o nlsl
endif
addprm.o	: addprm.f90 nlsdim.mod eprprm.mod expdat.mod parcom.mod lpnam.mod\
                  lmcom.mod stdio.mod rnddbl.mod
assgnc.o	: assgnc.f90 $(BASI)
bessel.o	: bessel.f90 rnddbl.mod
basisc.o	: basisc.f90 $(BASI)
catch.o         : catch.c fortrancall.h
cd2km.o         : cd2km.f90
cgltri.o	: cgltri.f90 nlsdim.mod rnddbl.mod
confc.o		: confc.f90 iterat.mod stdio.mod
convlv.o	: convlv.f90 nlsdim.mod
convtc.o	: convtc.f90 nlsdim.mod eprprm.mod stdio.mod lpnam.mod parcom.mod
covrpt.o	: covrpt.f90 nlsdim.mod eprprm.mod expdat.mod parcom.mod\
		  lmcom.mod lpnam.mod mspctr.mod iterat.mod rnddbl.mod\
		  stdio.mod
cscg.o		: cscg.f90 nlsdim.mod eprmat.mod rnddbl.mod 
cslnzs.o 	: cslnzs.f90 nlsdim.mod eprmat.mod
datac.o		: datac.f90 $(DATI) mspctr.mod lmcom.mod rnddbl.mod
eprls.o		: eprls.f90 $(MATI) ftwork.mod tridag.mod pidef.mod
fitc.o		: fitc.f90 nlsdim.mod lmcom.mod parcom.mod stdio.mod iterat.mod
fitl.o  	: fitl.f90 $(NLSI) $(FITI) rnddbl.mod timer.mod
fixc.o		: fixc.f90 nlsdim.mod eprprm.mod parcom.mod lpnam.mod stdio.mod
ftfuns.o	: ftfuns.f90
fz.o		: fz.f90 eprprm.mod
gconvl.o	: gconvl.f90 nlsdim.mod eprprm.mod ftwork.mod pidef.mod
genio.o         : genio.c fortrancall.h
getids.o	: getids.f90 nlsdim.mod stdio.mod nlsnam.mod
getdat.o	: getdat.f90 nlsdim.mod stdio.mod
helpc.o		: helpc.f90 stdio.mod
ipsfind.o 	: ipsfind.f90 $(BASI) lpnam.mod expdat.mod
lbasix.o        : lbasix.f90 nlsdim.mod eprprm.mod mtsdef.mod stdio.mod
lcheck.o	: lcheck.f90 $(CHKI)
letc.o		: $(LETI) stdio.mod
lfun.o		: lfun.f90 $(FUNI) $(FNI2) 
lgrint.o	: lgrint.f90
lmnls.o		: lmnls.f90 iterat.mod
lprmpt.o        : lprmpt.c fortrancall.h
l1pfun.o	: l1pfun.f90 nlsdim.mod expdat.mod parcom.mod
matrll.o	: matrll.f90 $(MATI) maxl.mod physcn.mod
mnbrak.o	: mnbrak.f90 errmsg.mod
momdls.o	: momdls.f90 $(MATI) pidef.mod stdio.mod dfunc.mod
nlsl.o		: nlsl.f90 $(NLSI) basis.mod iterat.mod mspctr.mod tridag.mod
#nlstxt.o	: nlstxt.f90 nlsdim.mod eprprm.mod lpnam.mod errmsg.mod
ordrpr.o	: ordrpr.f90 rnddbl.mod dfunc.mod pidef.mod
parc.o		: parc.f90 nlsdim.mod eprprm.mod expdat.mod parcom.mod lpnam.mod stdio.mod\
                  rnddbl.mod symdef.mod mtsdef.mod
pltx.o		: pltx.c fortrancall.h
pmatrl.o	: pmatrl.f90 $(MATI) maxl.mod physcn.mod
pstvec.o 	: pstvec.f90 nlsdim.mod eprprm.mod errmsg.mod rnddbl.mod
rdpar.o		: rdpar.f90 eprprm.mod stdio.mod
rms.o		: rms.f90 nlsdim.mod mspctr.mod parcom.mod rnddbl.mod
scalec.o        : scalec.f90 $(SCLI) stdio.mod
scmvm.o 	: scmvm.f90 nlsdim.mod rnddbl.mod eprmat.mod eprprm.mod stdio.mod
series.o	: series.f90 nlsdim.mod eprprm.mod expdat.mod mspctr.mod parcom.mod stdio.mod
setmts.o	: setmts.f90 nlsdim.mod eprprm.mod maxl.mod rnddbl.mod mtsdef.mod
setnm.o		: setnm.f90 nlsdim.mod nlsnam.mod
setprm.o	: setprm.f90 $(BASI) $(PRMI) iterat.mod stdio.mod symdef.mod
setspc.o	: setspc.f90 $(SPCI) iterat.mod
shiftc.o	: shiftc.f90 nlsdim.mod expdat.mod parcom.mod stdio.mod
sitec.o		: sitec.f90 nlsdim.mod expdat.mod parcom.mod tridag.mod basis.mod stdio.mod
srchc.o 	: srchc.f90 nlsdim.mod eprprm.mod lmcom.mod parcom.mod\
                  iterat.mod errmsg.mod lpnam.mod stdio.mod expdat.mod
sscale.o	: sscale.f90 nlsdim.mod
statc.o		: statc.f90 $(BASI) nlsdim.mod eprprm.mod expdat.mod parcom.mod\
                  lmcom.mod mtsdef.mod errmsg.mod
stats.o		: stats.f90
strutl1.o	: strutl1.f90 stdio.mod
strutl2.o	: strutl2.f90
stvect.o 	: stvect.f90 nlsdim.mod eprprm.mod errmsg.mod rnddbl.mod
tdchek.o        : tdchek.f90 $(BASI) tridag.mod errmsg.mod mtsdef.mod
tdqz.o		: tdsqz.f90 nlsdim.mod expdat.mod tridag.mod stdio.mod
varyc.o		: varyc.f90 nlsdim.mod eprprm.mod parcom.mod lpnam.mod stdio.mod
writr.o		: writr.f90 nlsdim.mod expdat.mod lmcom.mod parcom.mod stdio.mod\
                  tridag.mod mspctr.mod iterat.mod
writec.o	: writec.f90 nlsdim.mod expdat.mod parcom.mod lmcom.mod mspctr.mod stdio.mod
w3j.o		: w3j.f90 maxl.mod stdio.mod bincom.mod
zaxpy.o 	: zaxpy.f90 nlsdim.mod rnddbl.mod
qzaypx.o 	: zaypx.f90 nlsdim.mod rnddbl.mod
znormu.o	: znormu.f90 nlsdim.mod rnddbl.mod
zscsw.o		: zscsw.f90 nlsdim.mod rnddbl.mod
zdotu.o		: zdotu.f90 nlsdim.mod rnddbl.mod
nlsdim.o nlsdim.mod: nlsdim.f90
parcom.o parcom.mod: parcom.f90 nlsdim.mod
eprprm.o eprprm.mod: eprprm.f90 parcom.mod nlsdim.mod
errmsg.o errmsg.mod: errmsg.f90
lpnam.o lpnam.mod: lpnam.f90 nlsdim.mod
expdat.o expdat.mod: expdat.f90 nlsdim.mod
basis.o basis.mod: basis.f90 nlsdim.mod
stdio.o stdio.mod: stdio.f90
maxl.o maxl.mod: maxl.f90
bincom.o bincom.mod: bincom.f90 maxl.mod
eprmat.o eprmat.mod: eprmat.f90 nlsdim.mod
ftwork.o ftwork.mod: ftwork.f90 nlsdim.mod
iterat.o iterat.mod: iterat.f90
lmcom.o lmcom.mod: lmcom.f90 nlsdim.mod
mspctr.o mspctr.mod: mspctr.f90 nlsdim.mod
mtsdef.o mtsdef.mod: mtsdef.f90
nlsnam.o nlsnam.mod: nlsnam.f90 nlsdim.mod stdio.mod
parsav.o parsav.mod: parcom.f90 nlsdim.mod
physcn.o physcn.mod: physcn.f90
pidef.o pidef.mod: pidef.f90
rnddbl.o rnddbl.mod: rnddbl.f90
timer.o timer.mod: timer.f90
tridag.o tridag.mod: tridag.f90 nlsdim.mod
testmods.o: testmods.f90 nlsdim.mod parcom.mod eprprm.mod errmsg.mod lpnam.mod

#-----------------------------------------------------------------------
#		Executable files
#-----------------------------------------------------------------------

ifeq ($(OS),"windows")
nlsl.exe: $(NLSO) 
	$(FLINK) -o nlsl $(NLSO) 
else
nlsl: $(NLSO) 
	$(FLINK) -o $@ $(NLSO) $(LIBS)
endif

clean:
	$(RM) *.o *.mod nlsl testmods

TESTS = testmods.o eprprm.o parcom.o nlsdim.o errmsg.o lpnam.o ipsfind.o expdat.o basis.o strutl2.o
EXTRAS = strutl1.o symdef.o maxl.o bincom.o dfunc.o eprmat.o ftwork.o iterat.o lmcom.o mspctr.o mtsdef.o nlsnam.o parsav.o physcn.o pidef.o rnddbl.o timer.o tridag.o
# no tests yet for $(EXTRAS)... compile them anyway

testmods: $(TESTS) $(EXTRAS)
	$(FLINK) -o $@ $(TESTS)

testclean:
	$(RM) $(TESTS) $(EXTRAS) *.mod testmods

#-----------------------------------------------------------------------
#			Default actions
#-----------------------------------------------------------------------

%.o : %.c
	$(CC) $(CFLAGS) $*.c

%.o %.mod : %.f90
	if [ -a $*.mod ] ; then rm -f $*.mod ; fi
	$(F90) $(FFLAGS) -ffixed-form $*.f90
