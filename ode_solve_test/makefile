MY_PRECISION = -DUSE_DOUBLE_PRECISION


# d
include ../config/makefile.gcc
include ../config/makefile.petsc
ifndef MPILD
MPILD = $(LD)
endif

LOCAL_INCLUDES = -I../src/includes/ -I../src/util/ -I../src/block/

INCLUDES = $(LOCAL_INCLUDES) $(PETSC_INCLUDES)
BDIR = bin/
ODIR = objects/

# petsc h2opus and htool examples
all: $(BDIR)/ode_solve_test


$(BDIR)/ode_solve_test: $(ODIR)/ode_solve_test.o $(LIBLIST) 
	$(MPILD) -o $(BDIR)/ode_solve_test $(ODIR)/ode_solve_test.o \
	$(PETSC_LIBS) $(LIBS) $(SUPERLU_LIBS)  $(SLATEC_LIBS)   $(LDFLAGS)




clean: 
	rm -rf $(ODIR)/*.o $(ODIR)/*.a  
dist_clean:
	rm -rf $(BDIR)/* auto_proto.h auto_proto.sgl.h

obj_directories:
	if [ ! -s $(ODIR) ];then \
		mkdir -p $(ODIR);\
	fi;\
	if [ ! -s $(BDIR) ];then\
		mkdir -p $(BDIR);\
	fi;\

#
$(ODIR)/%.o:	%.c 
	$(CC) $(INCLUDES) $(CFLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f 
	$(F77) $(INCLUDES) $(FFLAGS) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.F 
	$(F77) $(INCLUDES) $(FFLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f90 
	$(F90)   $(INCLUDES) $(F90FLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

