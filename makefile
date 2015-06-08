########################################################################
#
#	 ====================================================
#	  MAKE descriptor file for EPRLL and related programs
#	 ====================================================
#
# VERSION 1.6  8/12/94
# Version for xlf Fortran under IBM AIX system for RISC/6000
# Adapted from NMAKE makefile for Microsoft Fortran version 5.1 under MS-DOS.
#
########################################################################
#

.SUFFIXES:
.SUFFIXES: .f .inc .o .c
.PRECIOUS: .f .inc .h .c

LIB = -L /usr/X11/lib64 -lX11 -L/usr/X11R6/lib64 -llapack -lblas 
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
F77 = gfortran
g77 = gfortran
CFLAGS = -c -O2 -DADD_UNDERSCORE
FLINK = gfortran
CC = gcc

#-----------------------------------------------------------------------
#			Definitions
#-----------------------------------------------------------------------
#

BSSLOBJS   = bssl.o getids.o ipar.o lbasix.o setnam.o setflg.o

EPRBLOBJS  = eprbl.o anxlk.o bessi.o bessi0.o bessi1.o cd2km.o\
        ccrint.o cgltri.o cscg.o cspccg.o dtime.o eprfsl.o fz.o\
        getids.o ipar.o lbasix.o matrll.o mtsl.o plgndr.o\
        prmpt.o rddat.o scmvm.o setflg.o setnam.o stvect.o\
        w3j.o wrdat.o wrparm.o zaxpy.o zaypx.o zdotu.o znormu.o

EPRCGLOBJS = eprcgl.o anxlk.o bessi.o bessi0.o bessi1.o cd2km.o\
        ccrint.o cgltri.o cscg.o dtime.o fz.o getids.o ipar.o\
        lbasix.o matrll.o plgndr.o prmpt.o rddat.o scmvm.o\
        setflg.o setnam.o stvect.o w3j.o wrdat.o wrparm.o\
        zaxpy.o zaypx.o znormu.o zdotu.o zscsw.o

EPRLLOBJS  = eprll.o anxlk.o bessi.o bessi0.o bessi1.o cd2km.o\
        ccrint.o cslnzs.o dtime.o fz.o getids.o ipar.o lbasix.o\
        matrll.o plgndr.o prmpt.o rddat.o scmvm.o setflg.o\
        setnam.o stvect.o w3j.o wrdat.o wrparm.o zaxpy.o\
        znormu.o zscsw.o zdotu.o

LBLLOBJS  = lbll.o cd2km.o getids.o ipar.o lbasix.o prmpt.o\
        rddat.o setflg.o setnam.o wrdat.o wrparm.o wrlbas.o

LMATRXOBJS = lmatrx.o getids.o ipar.o lbasix.o rddat.o setnam.o

LVECTROBJS  = lvectr.o getids.o ipar.o lbasix.o rddat.o setnam.o\
	stvect.o bessi.o bessi0.o bessi1.o fz.o ccrint.o\
	plgndr.o

MOMDLOBJS  = momdl.o anxlk.o bessi.o bessi0.o bessi1.o cd2km.o\
        ccrint.o cfvd.o cgltri.o cscg.o dtime.o fft.o gconvl.o\
        fz.o getids.o ipar.o lbasix.o matrll.o plgndr.o prmpt.o\
        rddat.o scmvm.o setflg.o setnam.o stvect.o w3j.o\
        wrdat.o wrparm.o zaxpy.o zaypx.o znormu.o zdotu.o\
        zscsw.o

ORDEROBJS  = order.o ccrint.o ccrin1.o

PRUNELOBJS  = prunel.o getids.o ipar.o lbasix.o prmpt.o rddat.o\
        setflg.o setnam.o

TDLLOBJS  = tdll.o cfvd.o cmtqli.o fft.o gconvl.o getids.o\
        setflg.o setnam.o prmpt.o rddat.o

#
#  Temporary or test files
#

LMATOBJS = lmat.o getids.o ipar.o lbasix.o rddat.o setnam.o


.suffixes:
.suffixes: .o .f

all: 
ifeq ($(OS),"windows")
all : bssl.exe eprbl.exe eprcgl.exe eprll.exe lbll.exe lmatrx.exe lvectr.exe momdl.exe order.exe prunel.exe tdll.exe
clean:
	echo "running clean for windows"
	$(RM) *.o *.mod bssl.exe eprbl.exe eprcgl.exe eprll.exe lbll.exe lmatrx.exe lvectr.exe momdl.exe order.exe prunel.exe tdll.exe
else
all             : bssl eprbl eprcgl eprll lbll lmatrx lvectr momdl order prunel tdll
clean           :
	echo "running clean for linux"
	$(RM) *.o *.mod 
endif

#
#-----------------------------------------------------------------------
#		object files
#-----------------------------------------------------------------------
#

anxlk.o  : anxlk.f rndoff.inc eprdat.inc
bessi.o  : bessi.f rndoff.inc
bessi0.o : bessi0.f rndoff.inc pidef.inc
bessi1.o : bessi1.f rndoff.inc pidef.inc
bssl.o   : bssl.f stdio.inc stddim.inc eprdat.inc fnames.inc indexl.inc\
           baswt.inc rndoff.inc version.inc
cd2km.o  : cd2km.f pidef.inc rndoff.inc
ccrint.o : ccrint.f
ccrin1.o : ccrin1.f
cfvd.o   : cfvd.f
cgltri.o : cgltri.f stddim.inc rndoff.inc
cmtqli.o : cmtqli.f rndoff.inc
cscg.o   : cscg.f stddim.inc rndoff.inc eprmat.inc cgdata.inc
cslnzs.o : cslnzs.f stddim.inc eprmat.inc
cspccg.o : cspccg.f stddim.inc rndoff.inc eprmat.inc cgdata.inc
#dtime.o  : dtime.c
eprbl.o  : eprbl.f stddim.inc stdio.inc eprdat.inc eprmat.inc spectr.inc\
             baswt.inc pidef.inc fnames.inc indexl.inc version.inc timer.inc
eprcgl.o : eprcgl.f stddim.inc stdio.inc eprdat.inc eprmat.inc tridag.inc\
             vectrs.inc fnames.inc version.inc timer.inc
eprfsl.o : eprfsl.f stddim.inc stdio.inc eprdat.inc vectrs.inc pidef.inc\
             fnames.inc indexl.inc timer.inc
eprll.o  : eprll.f stddim.inc stdio.inc eprdat.inc eprmat.inc tridag.inc\
             vectrs.inc fnames.inc version.inc timer.inc
fft.o    : fft.f
fz.o     : fz.f eprdat.inc
gconvl.o : gconvl.f stddim.inc eprdat.inc pidef.inc
getids.o : getids.f stddim.inc stdio.inc fnames.inc
ipar.o   : ipar.f
lbasix.o : lbasix.f stddim.inc stdio.inc eprdat.inc indexl.inc
lbll.o   : lbll.f stddim.inc maxl.inc stdio.inc eprdat.inc fnames.inc\
             indexl.inc rndoff.inc version.inc
lmatrx.o : lmatrx.f stddim.inc stdio.inc eprmat.inc eprdat.inc indexl.inc\
             fnames.inc rndoff.inc version.inc
lvectr.o : lvectr.f stddim.inc stdio.inc eprdat.inc indexl.inc fnames.inc\
	     rndoff.inc version.inc
matrll.o : matrll.f stddim.inc rndoff.inc eprmat.inc eprdat.inc\
             indexl.inc physcn.inc maxl.inc
momdl.o  : momdl.f stddim.inc stdio.inc eprdat.inc eprmat.inc tridag.inc\
	     vectrs.inc spectr.inc fnames.inc pidef.inc version.inc timer.inc\
             version.inc
mtsl.o   : mtsl.f maxl.inc stddim.inc indexl.inc eprdat.inc
order.o  : order.f
plgndr.o : plgndr.f
prmpt.o  : prmpt.f stdio.inc
prunel.o : prunel.f stdio.inc eprdat.inc stddim.inc fnames.inc indexl.inc\
	     rndoff.inc version.inc
rddat.o  : rddat.f stdio.inc eprdat.inc version.inc
scmvm.o  : scmvm.f stddim.inc rndoff.inc eprmat.inc
setflg.o : setflg.f stddim.inc fnames.inc
setnam.o : setnam.f stddim.inc fnames.inc
stvect.o : stvect.f stddim.inc eprdat.inc indexl.inc rndoff.inc
tdll.o   : tdll.f stddim.inc stdio.inc eprdat.inc fnames.inc tridag.inc\
             pidef.inc rndoff.inc version.inc
w3j.o    : w3j.f
wrdat.o  : wrdat.f stdio.inc eprdat.inc version.inc
wrlbas.o : wrlbas.f stddim.inc eprdat.inc indexl.inc
wrparm.o : wrparm.f rndoff.inc eprdat.inc
zaxpy.o  : zaxpy.f stddim.inc rndoff.inc
zaypx.o  : zaypx.f stddim.inc rndoff.inc
znormu.o : znormu.f stddim.inc rndoff.inc
zdotu.o  : zdotu.f stddim.inc rndoff.inc
zscsw.o  : zscsw.f stddim.inc

lmat : lmat.f stddim.inc stdio.inc eprmat.inc eprdat.inc indexl.inc\
             fnames.inc rndoff.inc version.inc

#
#-----------------------------------------------------------------------
#			executable files
#-----------------------------------------------------------------------
#

bssl :  $(BSSLOBJS)
	$(F77) $(LOADFLG) -o $@ $(BSSLOBJS) $(LIB)

eprbl : $(EPRBLOBJS)
	$(F77) $(LOADFLG) -o $@ $(EPRBLOBJS) $(LIB)


eprll : $(EPRLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(EPRLLOBJS) $(LIB)


eprcgl : $(EPRCGLOBJS)
	$(F77) $(LOADFLG) -o $@ $(EPRCGLOBJS) $(LIB)


lbll : $(LBLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(LBLLOBJS) $(LIB)


lmatrx : $(LMATRXOBJS)
	$(F77) $(LOADFLG) -o $@ $(LMATRXOBJS) $(LIB)


lvectr : $(LVECTROBJS)
	$(F77) $(LOADFLG) -o $@ $(LVECTROBJS) $(LIB)

momdl : $(MOMDLOBJS)
	$(F77) $(LOADFLG) -o $@ $(MOMDLOBJS) $(LIB)


order : $(ORDEROBJS)
	$(F77) $(LOADFLG) -o $@ $(ORDEROBJS) $(LIB)


prunel : $(PRUNELOBJS)
	$(F77) $(LOADFLG) -o $@ $(PRUNELOBJS) $(LIB)


tdll : $(TDLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(TDLLOBJS) $(LIB)

eprll.exe : $(EPRLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(EPRLLOBJS) $(LIB)


eprcgl.exe : $(EPRCGLOBJS)
	$(F77) $(LOADFLG) -o $@ $(EPRCGLOBJS) $(LIB)


lbll.exe : $(LBLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(LBLLOBJS) $(LIB)


lmatrx.exe : $(LMATRXOBJS)
	$(F77) $(LOADFLG) -o $@ $(LMATRXOBJS) $(LIB)


lvectr.exe : $(LVECTROBJS)
	$(F77) $(LOADFLG) -o $@ $(LVECTROBJS) $(LIB)

momdl.exe : $(MOMDLOBJS)
	$(F77) $(LOADFLG) -o $@ $(MOMDLOBJS) $(LIB)


order.exe : $(ORDEROBJS)
	$(F77) $(LOADFLG) -o $@ $(ORDEROBJS) $(LIB)


prunel.exe : $(PRUNELOBJS)
	$(F77) $(LOADFLG) -o $@ $(PRUNELOBJS) $(LIB)


tdll.exe : $(TDLLOBJS)
	$(F77) $(LOADFLG) -o $@ $(TDLLOBJS) $(LIB)

#
#  Temporary and test files
#

lmat: $(LMATOBJS)
	$(F77) $(LOADFLG) -o $@ $(LMATOBJS) $(LIB)


#
#-----------------------------------------------------------------------
#			default actions
#-----------------------------------------------------------------------
#
.f.o :
	$(F77) $(FFLAGS) $*.f

.inc.o :
	$(F77) $(FFLAGS) $*.f

.c.o   :
	cc -c $*.c
