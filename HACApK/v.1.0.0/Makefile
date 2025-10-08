#gnu compiler
#F90= mpifort -fbacktrace -g -fallow-argument-mismatch -ffree-form -ffree-line-length-none -O3 -march=native -fopenmp

#Sherlock (intel MPI)
#F90 = mpiifort -O3 -qopenmp -heap-arrays -ip -g -traceback

CC=mpicc
F90=mpifort

#OPT = -O3 -march=native -fno-strict-aliasing
OPT = -g

CFLAGSE = $(OPT)
F90FLAGS=-cpp -fbacktrace -fallow-argument-mismatch -ffree-form -ffree-line-length-none $(OPT) -fopenmp

F77=mpif77
F77FLAGS = $(OPT) -fopenmp

LDFLAGS= -fopenmp -llapack
#wisteria (Fujitsu compiler)
#F90=mpifrtpx -Kfast -Kopenmp

#F90FLAGS = $(OPTFLAG)
#LDFLAGS = -mkl=parallel

LD=$(F90)


OBJS= HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base.o m_HACApK_solve.o m_HACApK_use.o HACApK_c_interface.o 


all: test_hacpak_f90 test_hacpak_c

.SUFFIXES: .o .f90

test_hacpak_f90: $(OBJS) test_hacpak_f90.o
	$(LD) -o $@ $(OBJS) test_hacpak_f90.o $(LDFLAGS)

test_hacpak_c: $(OBJS) test_hacpak_c.o
	$(LD) -o $@ $(OBJS) test_hacpak_c.o $(LDFLAGS)

.f90.o: *.f90
	$(F90) -c $< $(F90FLAGS)
.c.o: *.c
	$(CC) -c $< $(CFLAGS)


clean:
	rm -f *.o *.mod $(TARGET)

rmod:
	rm -f m_*.o *.mod
