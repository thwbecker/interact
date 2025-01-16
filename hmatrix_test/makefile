MY_PRECISION = -DUSE_DOUBLE_PRECISION


# noise level for random interact version. this version 
# d
include ../makefile.gcc
include ../makefile.petsc
ifndef MPILD
MPILD = $(LD)
endif

BDIR = bin/
ODIR = objects/


all: $(BDIR)/ex21

$(BDIR)/ex21: $(ODIR)/ex21.o $(GEN_P_INC) $(LIBLIST) 
	$(MPILD) -o $(BDIR)/ex21 $(ODIR)/ex21.o $(PETSC_LIBS) $(LIBS) $(SUPERLU_LIBS)  $(SLATEC_LIBS)   $(LDFLAGS)



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
$(ODIR)/%.o:	%.c $(GEN_P_INC)
	$(CC) $(PETSC_INCLUDES) $(CFLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f $(GEN_P_INC)
	$(F77) $(PETSC_INCLUDES) $(FFLAGS) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.F $(GEN_P_INC)
	$(F77) $(PETSC_INCLUDES) $(FFLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f90 $(GEN_P_INC)
	$(F90)   $(PETSC_INCLUDES) $(F90FLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

