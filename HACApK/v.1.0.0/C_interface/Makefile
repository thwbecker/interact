#gnu compiler
#F90= mpifort -fbacktrace -g -fallow-argument-mismatch -ffree-form -ffree-line-length-none -O3 -march=x86-64-v3 -fopenmp

#Sherlock (intel MPI)
#F90 = mpiifort -O3 -qopenmp -heap-arrays -ip -g -traceback


CC=mpicc
F90=mpifort
F77=mpif77
LD=$(F90)

OPT = -O3 -march=x86-64-v3 -fno-strict-aliasing
#OPT = -g

CFLAGS = $(OPT)
F90FLAGS=-cpp -fbacktrace -fallow-argument-mismatch -ffree-form -ffree-line-length-none $(OPT) -fopenmp
F77FLAGS = $(OPT) -fopenmp
LDFLAGS= -fopenmp -llapack

#F90FLAGS = $(OPTFLAG)
#LDFLAGS = -mkl=parallel

ODIR=objects/

HACAPACK_OBJS= $(ODIR)/m_HACApK_calc_entry_ij.o $(ODIR)/HACApK_lib.o $(ODIR)/m_HACApK_base.o \
	$(ODIR)/m_HACApK_solve.o $(ODIR)/m_HACApK_use.o \
	$(ODIR)/HACApK_c_interface.o $(ODIR)/HACApK_testing_routines.o


all: libhacapk.a test_hacpak_c

libhacapk.a: $(HACAPACK_OBJS)
	$(AR) rv libhacapk.a $(HACAPACK_OBJS)

# C interface
test_hacpak_c: $(OBJS) libhacapk.a $(ODIR)/test_hacpak_c.o 
	$(LD) -o $@ $(OBJS) $(ODIR)/test_hacpak_c.o -L. -lhacapk $(LDFLAGS)

$(ODIR)/%.o: %.f90 HACApK_c_interface.h
	$(F90) -c $< $(F90FLAGS)  -o $(ODIR)/$*.o

$(ODIR)/%.o: %.c HACApK_c_interface.h
	$(CC) -c $< $(CFLAGS)  -o $(ODIR)/$*.o


clean:
	rm -f $(ODIR)/*.o *.mod 

