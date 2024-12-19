MY_PRECISION = -DUSE_DOUBLE_PRECISION



# noise level for random interact version. this version 
# d
include ../makefile.gcc
include ../makefile.petsc
ifndef MPILD
MPILD = $(LD)
endif

GEN_P_INC = $(PETSC_INCLUDES)
ODIR = objects/
BDIR = bin/

all: $(BDIR)/ex21

clean: 
	rm -rf $(ODIR)/*.o $(ODIR)/*.a  


$(BDIR)/ex21: $(ODIR)/ex21.o 
	$(MPILD) $(ODIR)/ex21.o -o $(BDIR)/ex21 \
		$(PETSC_LIBS) $(LIBS)  $(LDFLAGS)

$(ODIR)/ex21.o:	ex21.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c ex21.c -o $(ODIR)/ex21.o

$(ODIR)/%.o:	%.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o
