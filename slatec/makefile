#
#
# SLATEC libraries for non-negative least squares
#
#FFLAGS = -g -nofor_main 
#CFLAGS = -g
#
# $Id: makefile,v 1.1 2004/01/25 03:07:20 becker Exp becker $
#
#
# to compile, edit i1mach.f, r1mach.f and d1mach.f for individual
# machine constants and change the compiler options below
#
# SGI
#CFLAGS = -Ofast=$(ARCH)
#LOPT = -IPA
# LINUX
#CFLAGS = -O2
LOPT = $(MATHLIBS)
#LOPT = -lblas
#LOPT = -L/usr/local/src/atlas_Linux_P4SSE2/lib/  -lf77blas -latlas
#LOPT = -L/opt/lib/ -lblas 
#LOPT = -L/usr/local/lib -lblas
#
CFLAGS = -O2
FFLAGS = -O2 -fdollar-ok  -x f77-cpp-input -x f95-cpp-input -ffixed-line-length-132
F90FLAGS = -O2 -fdollar-ok  -x f77-cpp-input -x f95-cpp-input -ffixed-line-length-132
#FFLAGS = $(CFLAGS)
#> nnls 12 6 6
#i686/nnls: m(rows): 12 n(columns): 6, no NN rows
#i686/nnls: completely unconstrained least squares
#i686/nnls: dwnnls: rnorm:    0.857243 mode 0
# 4.3197307974e-01
# 3.5234583071e-01
# 2.9140619210e-01
# 1.8767460515e-01
# 6.2197116693e-02
#-1.5879688554e-01
#
# and 
#> nnls 12 6 0
#i686/nnls: m(rows): 12 n(columns): 6, NN: rows 1 to 6
#i686/nnls: completely non-negative  least squares
#i686/nnls: dwnnls: rnorm:    0.864629 mode 0
# 5.4529868644e-01
# 2.7165906671e-01
# 2.5196205254e-01
# 9.3772309486e-02
# 2.8522624145e-02
# 0.0000000000e+00
#
# same for the single precision version.
#

BIN = ./

L_SGL_FLAGS = -L.  -I. -lslatec_sgl $(LOPT) -DSINGLE_PREC

L_DBL_FLAGS = -L. -I.  -lslatec_dbl $(LOPT) 


SLATEC_DBL_OBJS = d1mach.o dh12.o dwnlit.o dwnlsm.o \
	dwnlt1.o \
	dwnlt2.o dwnlt3.o dwnnls.o fdump.o \
	i1mach.o j4save.o \
	xercnt.o xerhlt.o xermsg.o xerprn.o \
	xersve.o xgetua.o nnls_driver.o

SLATEC_SGL_OBJS = r1mach.o h12.o wnlit.o wnlsm.o \
	wnlt1.o \
	wnlt2.o wnlt3.o wnnls.o fdump.o \
	i1mach.o j4save.o \
	xercnt.o xerhlt.o xermsg.o xerprn.o \
	xersve.o xgetua.o nnls_driver_sgl.o


all:	obj_dirs \
	libslatec_sgl.a libslatec_dbl.a \
	 nnls   nnls_sgl

obj_dirs:
	if [ ! -s ./ ]; then\
		mkdir $(ARCH);\
	fi;

nnls: nnls.o libslatec_dbl.a
	$(CC) $(CFLAGS) nnls.o  -o $(BIN)/nnls $(L_DBL_FLAGS)  

nnls_sgl: nnls_sgl.o libslatec_sgl.a
	$(CC) $(CFLAGS) nnls_sgl.o  -o $(BIN)/nnls_sgl $(L_SGL_FLAGS)  


nnls_sgl.o: nnls.c
	$(CC) -c $(CFLAGS) $< -o nnls_sgl.o -DSINGLE_PREC

nnls.o: nnls.c
	$(CC) -c $(CFLAGS) $< -o nnls.o

nnls_driver_sgl.o: nnls_driver.c
	$(CC) -c $(CFLAGS) $< -o nnls_driver_sgl.o -DSINGLE_PREC

nnls_driver.o: nnls_driver.c
	$(CC) -c $(CFLAGS) $< -o nnls_driver.o


%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $*.o


%.o: %.f
	$(F77) -c $(FFLAGS) -fno-range-check $< -o $*.o

libslatec_sgl.a: $(SLATEC_SGL_OBJS)
	$(AR) -rv libslatec_sgl.a $(SLATEC_SGL_OBJS)

libslatec_dbl.a: $(SLATEC_DBL_OBJS)
	$(AR) -rv libslatec_dbl.a $(SLATEC_DBL_OBJS)


clean:
	rm *.o *.a
