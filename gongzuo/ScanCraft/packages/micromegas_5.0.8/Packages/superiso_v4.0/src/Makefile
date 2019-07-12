
.KEEP_STATE:

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libsuperiso.a

-include FlagsForMake

.c.o:
	$(CC) -c $(CFLAGS) $(CFLAGS_MP) $<
.c.a:
	$(CC) -c $(CFLAGS) $(CFLAGS_MP) $<
	$(AR) $(ARFLAGS) $@ $*.o; rm $*.o

all: libsuperiso.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) $(CFLAGS_MP) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac
	   
clean:
	rm -f *.a isajet.h softsusy.h spheno.h suspect.h 2hdmc.h nmssmtools.h;

resetpaths:
	rm -f isajet.h softsusy.h spheno.h suspect.h 2hdmc.h;
	
distclean: 
	rm -f *.a *.o *.x isajet.h softsusy.h spheno.h suspect.h 2hdmc.h;
	
libsuperiso.a: libsuperiso.a(softsusy.o) libsuperiso.a(leshouches.o)\
	  libsuperiso.a(isospin.o) libsuperiso.a(isajet.o) libsuperiso.a(wilson.o)\
	  libsuperiso.a(bsgamma.o) libsuperiso.a(alphas.o) libsuperiso.a(general.o)\
	  libsuperiso.a(excluded_masses.o) libsuperiso.a(quark_masses.o) libsuperiso.a(gmuon.o)\
	  libsuperiso.a(btaunu.o) libsuperiso.a(bdtaunu.o) libsuperiso.a(kmunu.o)\
	  libsuperiso.a(bsmumu.o) libsuperiso.a(dslnu.o) libsuperiso.a(dmunu.o)\
	  libsuperiso.a(spheno.o) libsuperiso.a(suspect.o) libsuperiso.a(2hdmc.o)\
	  libsuperiso.a(fleshouches.o) libsuperiso.a(bsll.o) libsuperiso.a(bsll_extra.o)\
	  libsuperiso.a(bkll.o) libsuperiso.a(bkstarll.o) libsuperiso.a(bkstargamma.o) libsuperiso.a(bsphill.o)\
	  libsuperiso.a(chi2.o) libsuperiso.a(zwidths.o) libsuperiso.a(nmssmtools.o)\
	  $(RANL)

libsuperiso.a(softsusy.o): softsusy.c include.h softsusy.h
libsuperiso.a(leshouches.o): leshouches.c include.h
libsuperiso.a(isospin.o): isospin.c include.h
libsuperiso.a(isajet.o): isajet.c include.h isajet.h
libsuperiso.a(wilson.o): wilson.c include.h
libsuperiso.a(bsgamma.o): bsgamma.c include.h
libsuperiso.a(alphas.o): alphas.c include.h
libsuperiso.a(quark_masses.o): quark_masses.c include.h
libsuperiso.a(general.o): general.c include.h
libsuperiso.a(excluded_masses.o): excluded_masses.c include.h
libsuperiso.a(gmuon.o): gmuon.c include.h
libsuperiso.a(btaunu.o): btaunu.c include.h
libsuperiso.a(bdtaunu.o): bdtaunu.c include.h
libsuperiso.a(kmunu.o): kmunu.c include.h
libsuperiso.a(bsmumu.o): bsmumu.c include.h
libsuperiso.a(dslnu.o): dslnu.c include.h
libsuperiso.a(dmunu.o): dmunu.c include.h
libsuperiso.a(spheno.o): spheno.c include.h spheno.h
libsuperiso.a(suspect.o): suspect.c include.h suspect.h
libsuperiso.a(2hdmc.o): 2hdmc.c include.h 2hdmc.h
libsuperiso.a(nmssmtools.o): nmssmtools.c nmssmtools.h include.h
libsuperiso.a(fleshouches.o): fleshouches.c include.h
libsuperiso.a(bsll.o): bsll.c include.h
libsuperiso.a(bsll_extra.o): bsll_extra.c include.h
libsuperiso.a(bkll.o): bkll.c include.h
libsuperiso.a(bkstarll.o): bkstarll.c include.h
libsuperiso.a(bkstargamma.o): bkstargamma.c include.h
libsuperiso.a(bsphill.o): bsphill.c include.h
libsuperiso.a(chi2.o): chi2.c include.h experimental_input.h correlations_exp_input.h
libsuperiso.a(zwidths.o): zwidths.c include.h

softsusy.h:
	@echo \#define SOFTSUSY \"$(SOFTSUSY)\" > softsusy.h;

isajet.h:
	@echo \#define ISAJET \"$(ISAJET)\" > isajet.h;

spheno.h:
	@echo \#define SPHENO \"$(SPHENO)\" > spheno.h;

suspect.h:
	@echo \#define SUSPECT \"$(SUSPECT)\" > suspect.h;
	
2hdmc.h:
	@echo \#define THDMC \"$(THDMC)\" > 2hdmc.h;

nmssmtools.h:
	@echo \#define NMSSMTools \"$(NMSSMTOOLS)\" > nmssmtools.h;
