#
# Unix/Linux Makefile Makefile for MATLAB interface to sparseLM
#


# HAVE_LARGE_ARRAYS should be set to 'yes' if your matlab version is 7.3 (R2006b) or later.
# Starting with this version, sparse matrix indices are 64 (rather than 32) bits long.
# Your matlab's version can be determined by typing 'version' at the matlab prompt.
#
# More details are at http://www.mathworks.com/support/solutions/en/data/1-5C27B9/?solution=1-5C27B9
#
#HAVE_LARGE_ARRAYS=no
HAVE_LARGE_ARRAYS=yes


MEX=mex
ifeq ($(strip $(HAVE_LARGE_ARRAYS)), yes)
  MEX:=$(MEX) -largeArrayDims -DUSE_LARGEARRAYDIMS
endif

MEXCFLAGS=-I.. -O #-g
# WHEN USING LAPACK, CHANGE THE NEXT TWO LINES TO WHERE YOUR COMPILED LAPACK/BLAS & F2C LIBS ARE!
LAPACKBLASLIBS_PATH=/home/sysbio/msrfby/GotoBLAS2
F2CLIBS_PATH=/home/sysbio/msrfby/libf2c

INTFACESRC=sparselm.c
INTFACESOBJ=sparselm.mexglx
UTILSRC=sparselm_spsolvr.c
UTILOBJ=sparselm_spsolvr.mexglx
# CURRENTLY, ONLY CHOLMOD IS SUPPORTED
#SPARSELIBS=/home/lourakis/SuiteSparse-3.4.0/CHOLMOD/Lib/libcholmod.a /home/lourakis/SuiteSparse-3.4.0/CAMD/Lib/libcamd.a\
/home/lourakis/SuiteSparse-3.4.0/AMD/Lib/libamd.a /home/lourakis/SuiteSparse-3.4.0/COLAMD/Lib/libcolamd.a\
/home/lourakis/SuiteSparse-3.4.0/CCOLAMD/Lib/libccolamd.a /usr/local/lib/libmetis.a
SPARSELIBS=/home/sysbio/msrfby/SuiteSparse/CHOLMOD/Lib/libcholmod.a /home/sysbio/msrfby/SuiteSparse/CAMD/Lib/libcamd.a\
/home/sysbio/msrfby/SuiteSparse/AMD/Lib/libamd.a /home/sysbio/msrfby/SuiteSparse/COLAMD/Lib/libcolamd.a\
/home/sysbio/msrfby/SuiteSparse/CCOLAMD/Lib/libccolamd.a /home/sysbio/msrfby/SuiteSparse/metis-4.0/libmetis.a


SPLMLIB=../libsplm.a
LAPACKLIBS=$(LAPACKBLASLIBS_PATH)/libgoto2.a $(F2CLIBS_PATH)/libf2c.a

LIBS=$(SPLMLIB) $(SPARSELIBS) $(LAPACKLIBS)

dummy: $(INTFACESOBJ) $(UTILOBJ)

$(INTFACESOBJ): $(INTFACESRC) $(SPLMLIB)
	$(MEX) $(MEXCFLAGS) $(INTFACESRC) $(LIBS)

$(UTILOBJ): $(UTILSRC)
	$(MEX) $(MEXCFLAGS) $(UTILSRC)

clean:
	@rm -f $(INTFACESOBJ) $(UTILOBJ)

realclean cleanall: clean

depend:
	makedepend -f Makefile $(INTFACESRC) $(UTILOBJ)

# DO NOT DELETE THIS LINE -- make depend depends on it.

