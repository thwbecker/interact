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

all: $(BDIR)/ex21 $(BDIR)/ex21-kdtree $(BDIR)/ex21-htool

clean: 
	rm -rf $(ODIR)/*.o $(ODIR)/*.a  


$(BDIR)/ex21: $(ODIR)/ex21.o 
	$(MPILD) $(ODIR)/ex21.o -o $(BDIR)/ex21 \
		$(PETSC_LIBS) $(LIBS)  $(LDFLAGS)

$(BDIR)/ex21-htool: $(ODIR)/ex21-htool.o 
	$(MPILD) $(ODIR)/ex21-htool.o -o $(BDIR)/ex21-htool \
		$(PETSC_LIBS) $(LIBS)  $(LDFLAGS)

$(BDIR)/ex21-kdtree: $(ODIR)/ex21-kdtree.o 
	$(MPILD) $(ODIR)/ex21-kdtree.o -o $(BDIR)/ex21-kdtree \
		$(PETSC_LIBS) $(LIBS)  $(LDFLAGS)

$(ODIR)/ex21.o:	ex21.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c ex21.c -o $(ODIR)/ex21.o

$(ODIR)/ex21-kdtree.o:	ex21-kdtree.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c ex21-kdtree.c -o $(ODIR)/ex21-kdtree.o

$(ODIR)/ex21-htool.o:	ex21-htool.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c ex21-htool.c -o $(ODIR)/ex21-htool.o

$(ODIR)/%.o:	%.c 
	$(CC) $(CFLAGS)  $(GEN_P_INC) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o
