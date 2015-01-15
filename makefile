########################################################################
#
#	     ====================================================
#	     Make descriptor file for NLSL and related programs
#	     ====================================================
#
########################################################################

.SUFFIXES:
.SUFFIXES: .f90 .f .inc .o .c 
.PRECIOUS: .f90 .f .inc .h .c

########################################################################

F77=gfortran
FFLAGS = -c -O2 -g -fopenmp -std=gnu -m64 -mcmodel=medium
LIB = -L /usr/X11/lib64 -lX11 -L/usr/X11R6/lib64

#CC = gcc 
#LIBS = dfor.lib
#LIBS2  = C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsl.lib C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsls_err.lib 

# gcc on i386 Linux
 CC = gcc
#F77 = g77
#FFLAGS = -c -O9 -fomit-frame-pointer -ffast-math -malign-double\
	-funroll-loops -march=i586
#FFLAGS = -c -g
 CFLAGS = -c -O2 -DADD_UNDERSCORE
 FLINK = gfortran
 LIB = -L/usr/X11R6/lib
########################################################################

#	
#-----------------------------------------------------------------------
#			Definitions
#-----------------------------------------------------------------------
#


SPCI = nlsdim.inc eprprm.inc expdat.inc parcom.inc tridag.inc basis.inc pidef.inc
BASI = nlsdim.inc eprprm.inc parcom.inc basis.inc stdio.inc
LETI = nlsdim.inc eprprm.inc expdat.inc prmeqv.inc parcom.inc lpnam.inc
NLSI = nlsdim.inc eprprm.inc expdat.inc lmcom.inc stdio.inc parcom.inc
PRMI = nlsdim.inc eprprm.inc lmcom.inc parcom.inc lpnam.inc
FUNI = nlsdim.inc eprprm.inc expdat.inc tridag.inc parcom.inc iterat.inc
FNI2 = basis.inc pidef.inc stdio.inc prmeqv.inc mspctr.inc ftwork.inc 
FITI = nlsnam.inc eprmat.inc errmsg.inc mspctr.inc iterat.inc
MATI = nlsdim.inc eprprm.inc eprmat.inc rndoff.inc
CHKI = maxl.inc nlsdim.inc eprprm.inc rndoff.inc prmeqv.inc errmsg.inc
DATI = nlsdim.inc expdat.inc mspctr.inc nlsnam.inc stdio.inc parcom.inc 
SCLI = nlsdim.inc expdat.inc lmcom.inc parcom.inc mspctr.inc iterat.inc 
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
NLSS = strutl.o lprmpt.o catch.o ipfind.o nlstxt.o  
NLSO = nlsl.o $(NLSC) $(NLSS) $(NLSD) $(NLSF) $(NLSL) $(NLSP) $(NLSN) $(NLSB)\
       $(NLSW)

#-----------------------------------------------------------------------
#		Object files
#-----------------------------------------------------------------------

all             : nlsl
clean           :
			rm -f *.o nlsl
addprm.o	: addprm.f nlsdim.inc eprprm.inc expdat.inc parcom.inc lpnam.inc\
                  lmcom.inc stdio.inc rndoff.inc prmeqv.inc
assgnc.o	: assgnc.f $(BASI)
bessel.o	: bessel.f rndoff.inc
basisc.o	: basisc.f $(BASI)
catch.o         : catch.c fortrancall.h
cd2km.o         : cd2km.f
cgltri.o	: cgltri.f nlsdim.inc rndoff.inc
confc.o		: confc.f iterat.inc stdio.inc
convlv.o	: convlv.f nlsdim.inc
convtc.o	: convtc.f nlsdim.inc eprprm.inc stdio.inc lpnam.inc parcom.inc
covrpt.o	: covrpt.f nlsdim.inc eprprm.inc expdat.inc parcom.inc\
		  lmcom.inc lpnam.inc mspctr.inc iterat.inc rndoff.inc\
		  stdio.inc
cscg.o		: cscg.f nlsdim.inc eprmat.inc rndoff.inc 
cslnzs.o 	: cslnzs.f nlsdim.inc eprmat.inc
datac.o		: datac.f $(DATI) mspctr.inc lmcom.inc rndoff.inc
eprls.o		: eprls.f $(MATI) ftwork.inc tridag.inc pidef.inc
fitc.o		: fitc.f nlsdim.inc lmcom.inc parcom.inc stdio.inc iterat.inc
fitl.o  	: fitl.f $(NLSI) $(FITI) rndoff.inc timer.inc
fixc.o		: fixc.f nlsdim.inc eprprm.inc prmeqv.inc parcom.inc lpnam.inc stdio.inc
ftfuns.o	: ftfuns.f
fz.o		: fz.f eprprm.inc
gconvl.o	: gconvl.f nlsdim.inc eprprm.inc ftwork.inc pidef.inc
genio.o         : genio.c fortrancall.h
getids.o	: getids.f nlsdim.inc stdio.inc nlsnam.inc
getdat.o	: getdat.f nlsdim.inc stdio.inc
helpc.o		: helpc.f stdio.inc
ipfind.o	: ipfind.f $(BASI) eprprm.inc parcom.inc lpnam.inc
lbasix.o        : lbasix.f nlsdim.inc eprprm.inc mtsdef.inc stdio.inc
lcheck.o	: lcheck.f $(CHKI) prmeqv.inc
letc.o		: $(LETI) stdio.inc
lfun.o		: lfun.f $(FUNI) $(FNI2) 
lgrint.o	: lgrint.f
lmnls.o		: lmnls.f iterat.inc
lprmpt.o        : lprmpt.c fortrancall.h
l1pfun.o	: l1pfun.f nlsdim.inc expdat.inc parcom.inc
matrll.o	: matrll.f $(MATI) maxl.inc physcn.inc
mnbrak.o	: mnbrak.f errmsg.inc
momdls.o	: momdls.f $(MATI) pidef.inc stdio.inc dfunc.inc
nlsl.o		: nlsl.f $(NLSI) parcom.inc tridag.inc basis.inc iterat.inc
nlstxt.o	: nlstxt.f nlsdim.inc eprprm.inc lpnam.inc errmsg.inc
ordrpr.o	: ordrpr.f rndoff.inc dfunc.inc pidef.inc
parc.o		: parc.f nlsdim.inc eprprm.inc expdat.inc parcom.inc lpnam.inc stdio.inc\
                  rndoff.inc symdef.inc mtsdef.inc
pltx.o		: pltx.c fortrancall.h
pmatrl.o	: pmatrl.f $(MATI) maxl.inc physcn.inc
pstvec.o 	: pstvec.f nlsdim.inc eprprm.inc errmsg.inc rndoff.inc
rdpar.o		: rdpar.f eprprm.inc stdio.inc
rms.o		: rms.f nlsdim.inc mspctr.inc parcom.inc rndoff.inc
scalec.o        : scalec.f $(SCLI) stdio.inc
scmvm.o 	: scmvm.f nlsdim.inc rndoff.inc eprmat.inc eprprm.inc stdio.inc
series.o	: series.f nlsdim.inc eprprm.inc expdat.inc mspctr.inc parcom.inc stdio.inc
setmts.o	: setmts.f nlsdim.inc eprprm.inc maxl.inc rndoff.inc mtsdef.inc
setnm.o		: setnm.f nlsdim.inc nlsnam.inc
setprm.o	: setprm.f $(BASI) $(PRMI) iterat.inc stdio.inc symdef.inc
setspc.o	: setspc.f $(SPCI) iterat.inc
shiftc.o	: shiftc.f nlsdim.inc expdat.inc parcom.inc stdio.inc
sitec.o		: sitec.f nlsdim.inc expdat.inc parcom.inc tridag.inc basis.inc stdio.inc
srchc.o 	: srchc.f nlsdim.inc eprprm.inc lmcom.inc parcom.inc\
                  iterat.inc errmsg.inc lpnam.inc stdio.inc expdat.inc
sscale.o	: sscale.f nlsdim.inc
statc.o		: statc.f $(BASI) nlsdim.inc eprprm.inc expdat.inc parcom.inc\
                  lmcom.inc mtsdef.inc
stats.o		: stats.f
strutl.o	: strutl.f stdio.inc
stvect.o 	: stvect.f nlsdim.inc eprprm.inc errmsg.inc rndoff.inc
tdchek.o        : tdchek.f $(BASI) tridag.inc errmsg.inc mtsdef.inc
tdqz.o		: tdsqz.f nlsdim.inc expdat.inc tridag.inc
varyc.o		: varyc.f nlsdim.inc eprprm.inc prmeqv.inc parcom.inc lpnam.inc stdio.inc
writr.o		: writr.f nlsdim.inc expdat.inc lmcom.inc parcom.inc stdio.inc\
                  tridag.inc mspctr.inc iterat.inc
writec.o	: writec.f nlsdim.inc expdat.inc parcom.inc lmcom.inc mspctr.inc stdio.inc
w3j.o		: w3j.f maxl.inc bincom.inc stdio.inc
zaxpy.o 	: zaxpy.f nlsdim.inc rndoff.inc
zaypx.o 	: zaypx.f nlsdim.inc rndoff.inc
znormu.o	: znormu.f nlsdim.inc rndoff.inc
zscsw.o		: zscsw.f nlsdim.inc rndoff.inc
zdotu.o		: zdotu.f nlsdim.inc rndoff.inc
nlsdim.o	: nlsdim.f90
parcom.o	: parcom.f90
eprprm.o	: eprprm.f90

#-----------------------------------------------------------------------
#		Executable files
#-----------------------------------------------------------------------

nlsl	: $(NLSO) 
	$(FLINK) -o $@ $(NLSO) $(LIB) -lX11 -lc

#-----------------------------------------------------------------------
#			Default actions
#-----------------------------------------------------------------------

testmods: nlsdim.o parcom.o eprprm.o testmods.f90
	$(F77) -g testmods.f90 nlsdim.o parcom.o eprprm.o -o testmods

.c.o   :
	$(CC) $(CFLAGS) $*.c

.f.o   :
	$(F77) $(FFLAGS) $*.f

.f90.o   :
	$(F77) $(FFLAGS) -ffixed-form $*.f90
