#gnu compiler
#F90= mpifort -fbacktrace -g -fallow-argument-mismatch -ffree-form -ffree-line-length-none -O3 -march=native -fopenmp

#Sherlock (intel MPI)
#F90 = mpiifort -O3 -qopenmp -heap-arrays -ip -g -traceback

F90=mpifort

#OPT = -O3 -march=native
OPT = -g

F90FLAGS=-cpp -fbacktrace -fallow-argument-mismatch -ffree-form -ffree-line-length-none $(OPT) -fopenmp
F77=mpif77
F77FLAGS = $(OPT) -fopenmp
LDFLAGS= -fopenmp -llapack
#wisteria (Fujitsu compiler)
#F90=mpifrtpx -Kfast -Kopenmp

#F90FLAGS = $(OPTFLAG)
#LDFLAGS = -mkl=parallel

LD=$(F90)


OBJS= HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base.o m_HACApK_solve.o m_HACApK_use.o test_hacpak.o


TARGET=test_hacpak
.SUFFIXES: .o .f90

$(TARGET): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

.f90.o: *.f90
	$(F90) -c $< $(F90FLAGS)
clean:
	rm -f *.o *.mod $(TARGET)

rmod:
	rm -f m_*.o *.mod
