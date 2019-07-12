.KEEP_STATE:

VERSION = v4.0
YEAR = 2018

# Choose your compilers here (in general gcc/gfortran on Linux systems):
CC = gcc
CFLAGS = -O3 -pipe -fomit-frame-pointer
CFLAGS_MP = 
#CFLAGS_MP = -fopenmp

#CC = icc
#CFLAGS = -O3
#CFLAGS_MP = -qopenmp

MAKE = make
MAKE_MP=
#MAKE_MP= -j$(nproc)

AR = ar
ARFLAGS = rcsU

# Add the link to Softsusy, Isajet, SuSpect, SPheno and NMSSMTools main programs, if available.
# Otherwise, comment them in the main programs */
# If you change the path, you must "make distclean" first.
SOFTSUSY = ~/softsusy/softpoint.x
ISAJET = ~/isajet/isasugra.x
SUSPECT = ~/suspect/suspect2
SPHENO = ~/spheno/bin/SPheno
NMSSMTOOLS = ~/nmssmtools/main
# Add the link to 2HDMC directory, if available.
THDMC = ~/2HDMC

############## DO NOT MODIFY THE FOLLOWING ##############

CINCLUDE= -I./src -L./src

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libsuperiso.a

all: libsuperiso.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.c" into "name.x"';\
	echo ' '

%:: %.c libsuperiso.a
	$(CC) -c $(CFLAGS) $(CFLAGS_MP) $*.c;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lsuperiso -lm;
	@rm -f $*.o;
	@touch $*.x

%.c:: %.c libsuperiso.a
	$(CC) -c $(CFLAGS_MP) $(CFLAGS) $@;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lsuperiso -lm;
	@rm -f $*.o;
	@touch $*.x

save:
	rm -f superiso_$(VERSION).tgz;\
	mkdir superiso_$(VERSION);\
	cp -p README superiso_$(VERSION)/;\
	cp -p example.lha superiso_$(VERSION)/;\
	cp -p amsb.c superiso_$(VERSION)/;\
	cp -p cnmssm.c superiso_$(VERSION)/;\
	cp -p hcamsb.c superiso_$(VERSION)/;\
	cp -p mmamsb.c superiso_$(VERSION)/;\
	cp -p gmsb.c superiso_$(VERSION)/;\
	cp -p cmssm.c superiso_$(VERSION)/;\
	cp -p ngmsb.c superiso_$(VERSION)/;\
	cp -p nnuhm.c superiso_$(VERSION)/;\
	cp -p nuhm.c superiso_$(VERSION)/;\
	cp -p slha.c superiso_$(VERSION)/;\
	cp -p slha_chi2.c superiso_$(VERSION)/;\
	cp -p sm.c superiso_$(VERSION)/;\
	cp -p sm_chi2.c superiso_$(VERSION)/;\
	cp -p thdm.c superiso_$(VERSION)/;\
	cp -p flha.c superiso_$(VERSION)/;\
	cp -p Makefile superiso_$(VERSION)/;\
	mkdir superiso_$(VERSION)/src;\
	cp -p src/*.h superiso_$(VERSION)/src/;\
	cp -p src/*.c superiso_$(VERSION)/src/;\
	rm -f superiso_$(VERSION)/src/softsusy.h;\
	rm -f superiso_$(VERSION)/src/nmssmtools.h;\
	rm -f superiso_$(VERSION)/src/isajet.h;\
	rm -f superiso_$(VERSION)/src/spheno.h;\
	rm -f superiso_$(VERSION)/src/suspect.h;\
	rm -f superiso_$(VERSION)/src/2hdmc.h;\
	cp -p src/Makefile superiso_$(VERSION)/src/;\
	tar czvf superiso_$(VERSION).tgz superiso_$(VERSION);\
	rm -rf superiso_$(VERSION)
	
libsuperiso.a:
	@echo;
	@echo SuperIso $(VERSION) - F. Mahmoudi $(YEAR);
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo CFLAGS_MP = $(CFLAGS_MP) >> src/FlagsForMake;\
	echo MAKE = $(MAKE) >> src/FlagsForMake;\
	echo AR = $(AR) >> src/FlagsForMake;\
	echo ARFLAGS = $(ARFLAGS) >> src/FlagsForMake;
	echo SOFTSUSY = $(SOFTSUSY) >> src/FlagsForMake;\
	echo ISAJET = $(ISAJET) >> src/FlagsForMake;\
	echo SPHENO = $(SPHENO) >> src/FlagsForMake;\
	echo SUSPECT = $(SUSPECT) >> src/FlagsForMake;\
	echo THDMC = $(THDMC) >> src/FlagsForMake;\
	echo NMSSMTOOLS = $(NMSSMTOOLS) >> src/FlagsForMake;\
	$(MAKE) $(MAKE_MP) -C src/ libsuperiso.a

clean:
	rm -rf tmp *.x *.bin .*.bin *.tmplha *.is? *.ss *.sptmp *.sutmp libsuperiso.a;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ clean;
	
distclean:
	rm -rf tmp *.x *.bin .*.bin *.tmplha *.fh *.is? *.ss *.sptmp *.sutmp libsuperiso.a;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ distclean;
