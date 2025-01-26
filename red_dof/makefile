# makefile for interact and related programs
# Thorsten Becker, February 2000
# $Id: makefile,v 1.11 2000/08/04 20:57:19 becker Exp becker $
CC = cc
F77 = f77
# machine dependent optimization flags are in
# makefile.machine
# here, set other flags:
# set -DUSE_DOUBLE_PRECISION    for double precision, if left away, uses single prec
# set -DBINARY_PATCH_EVENT_FILE for binary output of single event files in interact
# set -DNO_OPENING_MODES        for restricted operation in strike and dip direction only
#ARGS = 
ARGS =  -DBINARY_PATCH_EVENT_FILE -DNO_OPENING_MODES  
# choice of Okada routine,
# comment out if modified routine is used
OKROUTINE = dc3d.o
# use this for modified Okada routine
#ARGS = -O -DUSEMYDC3D
#OKROUTINE = mydc3d.o
# debugging flag, comment out the above and use the following
# instead for debugging
#ARGS = -g -DDEBUG
# machine dependent flags
include makefile.sgi
# add this for pgplot support, otherwise comment it out
include makefile.pgplot
# assemble all arguments, ie. add compiler flags,
# system dependent flags, and pgplot flags if set
FLAGS = $(ARGS) $(PARGS)
# C and FORTRAN compiler specific flags, add them to the 
# other flags
CFLAGS = $(FLAGS) $(SCARGS)
FFLAGS = $(FLAGS) $(SFARGS)

MATRIX_SOLVER_OBJS = numrec_svd.o nnls_lawson.o
COMPUTATIONAL_LIBS = 

ALL_INTERACT_OBJS = $(OKROUTINE) eval_green.o geometry.o input.o output.o \
	init.o interact.o calc_stress.o rupture.o svd.o \
	quake.o fracture_criterion.o sincos.o solve.o \
	fltcopy.o optimize.o myrand.o nnls.o myopen.o \
	$(MATRIX_SOLVER_OBJS)
RANDOMFLT_OBJS = randomflt.o divide_fault_in_patches.o \
	sincos.o geometry.o myrand.o tritri.o far_enough.o

RANDOMIZE_STRIKE_OBJS = randomize_strike.o \
	divide_fault_in_patches.o tritri.o \
	sincos.o geometry.o myrand.o far_enough.o

OBJ = main.o $(ALL_INTERACT_OBJS)
TOBJ = test_stuff.o $(ALL_INTERACT_OBJS)

INC = interact.h precision_single.h precision_double.h
LIBS = -lm $(COMPUTATIONAL_LIBS)




interact: $(OBJ) $(INC) $(POBJS)
	$(F77) $(FFLAGS) $(POBJS) $(OBJ) -o interact \
	$(PGLIBS) $(LIBS) 

test_stuff: $(TOBJ) $(INC) $(POBJS)
	$(F77) $(FFLAGS) $(POBJS) $(TOBJ) -o test_stuff \
	$(PGLIBS) $(LIBS)

randomflt: $(RANDOMFLT_OBJS)  $(INC)
	$(CC) $(CFLAGS) $(RANDOMFLT_OBJS) \
	-o randomflt -lm

patch2xyz: divide_fault_in_patches.o sincos.o \
	geometry.o patch2xyz.o  $(INC)
	$(CC) $(CFLAGS) divide_fault_in_patches.o \
		sincos.o geometry.o patch2xyz.o \
		-o patch2xyz -lm
patch2geom: divide_fault_in_patches.o sincos.o \
	geometry.o patch2geom.o  $(INC)
	$(CC) $(CFLAGS) divide_fault_in_patches.o \
		sincos.o geometry.o patch2geom.o \
		-o patch2geom -lm
patch2bc: divide_fault_in_patches.o sincos.o \
	geometry.o patch2bc.o  $(INC)
	$(CC) $(CFLAGS) divide_fault_in_patches.o \
		sincos.o geometry.o patch2bc.o \
		-o patch2bc -lm

randomize_strike: $(RANDOMIZE_STRIKE_OBJS)  $(INC)
	$(CC) $(CFLAGS) $(RANDOMIZE_STRIKE_OBJS) \
	-o randomize_strike -lm

makefault: makefault.o divide_fault_in_patches.o sincos.o\
	geometry.o	$(INC)
	$(CC) $(CFLAGS) makefault.o \
	divide_fault_in_patches.o sincos.o geometry.o \
	-o makefault -lm

project_stress: project_stress.o geometry.o  $(INC)
	$(CC) $(CFLAGS) project_stress.o  geometry.o \
	-o project_stress -lm

all: interact test_stuff randomflt  makefault \
	project_stress randomize_strike \
	patch2xyz patch2geom patch2bc

clean: 
	rm -rf *.o *~ *.ps *.dat rii_files/

makefault.o: makefault.c $(INC) properties.h
	$(CC) $(CFLAGS) -c makefault.c 
project_stress.o: project_stress.c $(INC)
	$(CC) $(CFLAGS) -c project_stress.c 
fltcopy.o: fltcopy.c $(INC)
	$(CC) $(CFLAGS) -c fltcopy.c 
myopen.o: myopen.c $(INC)
	$(CC) $(CFLAGS) -c myopen.c 

patch2xyz.o: patch2xyz.c $(INC)
	$(CC) $(CFLAGS) -c patch2xyz.c 

patch2geom.o: patch2geom.c $(INC)
	$(CC) $(CFLAGS) -c patch2geom.c 

patch2bc.o: patch2bc.c $(INC)
	$(CC) $(CFLAGS) -c patch2bc.c 

far_enough.o: far_enough.c $(INC)
	$(CC) $(CFLAGS) -c far_enough.c 

randomflt.o: randomflt.c $(INC)
	$(CC) $(CFLAGS) -c randomflt.c 
randomize_strike.o: randomize_strike.c $(INC)
	$(CC) $(CFLAGS) -c randomize_strike.c 
main.o: main.c $(INC) properties.h
	$(CC) -c $(CFLAGS) main.c 
rupture.o: rupture.c $(INC) 
	$(CC) -c $(CFLAGS) rupture.c 
fracture_criterion.o: fracture_criterion.c $(INC) \
	properties.h
	$(CC) -c $(CFLAGS) fracture_criterion.c 
eval_green.o: eval_green.c $(INC) properties.h
	$(CC) -c $(CFLAGS) eval_green.c
geometry.o: geometry.c $(INC)
	$(CC) -c $(CFLAGS) geometry.c
input.o: input.c $(INC) properties.h
	$(CC) -c $(CFLAGS) input.c
output.o: output.c $(INC) properties.h
	$(CC) -c $(CFLAGS) output.c
init.o: init.c $(INC) properties.h
	$(CC) -c $(CFLAGS) init.c
svd.o: svd.c $(INC) svd.h
	$(CC) -c $(CFLAGS) svd.c
nnls.o: nnls.c $(INC) nnls.h
	$(CC) -c $(CFLAGS) nnls.c
optimize.o: optimize.c $(INC) svd.h
	$(CC) -c $(CFLAGS) optimize.c
solve.o: solve.c $(INC) 
	$(CC) -c $(CFLAGS) solve.c
numrec_svd.o: numrec_svd.F $(INC)
	$(F77) -c $(FFLAGS) numrec_svd.F 
nnls_lawson.o: nnls_lawson.F $(INC) 
	$(F77) -c $(FFLAGS) nnls_lawson.F
sincos.o: sincos.c $(INC) 
	$(CC) -c $(CFLAGS) sincos.c
quake.o: quake.c $(INC) 
	$(CC) -c $(CFLAGS) quake.c
divide_fault_in_patches.o: divide_fault_in_patches.c $(INC) \
	properties.h
	$(CC) -c $(CFLAGS) divide_fault_in_patches.c
calc_stress.o: calc_stress.c $(INC) properties.h
	$(CC) -c $(CFLAGS) calc_stress.c
myrand.o: myrand.c $(INC) properties.h
	$(CC) -c $(CFLAGS) myrand.c
mydc3d.o: mydc3d.F $(INC)
	$(F77) -c $(FFLAGS) mydc3d.F
dc3d.o: dc3d.F $(INC)
	$(F77) -c $(FFLAGS) dc3d.F
test_stuff.o: test_stuff.c $(INC)
	$(CC) -c $(CFLAGS) test_stuff.c

segseg.o: segseg.c $(INC)
	$(CC) -c $(CFLAGS) segseg.c


# these are only for PGPLOT support, we can leave them in
# anyhows
plotting.o: plotting.c $(INC)
	$(CC) $(CFLAGS) -c plotting.c 
plotting_palette.o: plotting_palette.F $(INC)
	$(F77) $(FFLAGS) -c plotting_palette.F
plotting_palette_.o: plotting_palette_.F $(INC)
	$(F77) $(FFLAGS) -c plotting_palette_.F

# some generic rules

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(F77) $(FFLAGS) -c $<
