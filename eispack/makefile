#
#
# makefile for EISPACK rs symmetric matrix eigenvalues and eigenvectors routine
# in single and double precision. to avoid confusion, all single precision 
# subroutines within the .sgl.F files have a s_ before the normal name. e.g.
# rs(...) turns into s_rs
#
# set architecture flag
#ARCH=$(shell uname -m | awk '{print(tolower($$1))}')


#ODIR = $(ARCH)/
ODIR = .

all: directories libraries

libraries: $(ODIR)/libmyeis_dbl.a $(ODIR)/libmyeis_sgl.a \
	$(ODIR)/test_dbl $(ODIR)/test_sgl

clean: 
	rm $(ODIR)/*.o $(ODIR)/*.a 

SGL_OBJS =  $(ODIR)/rs.sgl.o $(ODIR)/rg.sgl.o $(ODIR)/pythag.sgl.o \
	$(ODIR)/tql1.sgl.o $(ODIR)/tql2.sgl.o $(ODIR)/tred1.sgl.o \
	$(ODIR)/tred2.sgl.o $(ODIR)/hqr2.sgl.o $(ODIR)/hqr.sgl.o \
	$(ODIR)/balanc.sgl.o $(ODIR)/balbak.sgl.o $(ODIR)/eltran.sgl.o \
	$(ODIR)/cdiv.sgl.o $(ODIR)/elmhes.sgl.o 

DBL_OBJS =  $(ODIR)/rs.dbl.o $(ODIR)/rg.dbl.o $(ODIR)/pythag.dbl.o \
	$(ODIR)/tql1.dbl.o $(ODIR)/tql2.dbl.o $(ODIR)/tred1.dbl.o \
	$(ODIR)/tred2.dbl.o $(ODIR)/hqr2.dbl.o $(ODIR)/hqr.dbl.o \
	$(ODIR)/balanc.dbl.o $(ODIR)/balbak.dbl.o $(ODIR)/eltran.dbl.o \
	$(ODIR)/cdiv.dbl.o $(ODIR)/elmhes.dbl.o 


directories:
	if [ ! -s $(ODIR)/ ];then \
		mkdir $(ODIR);\
	fi



#
# LIBRARIES
#
# double precision objects
$(ODIR)/libmyeis_dbl.a: $(DBL_OBJS) 
	$(AR) -rv  $(ODIR)/libmyeis_dbl.a  $(DBL_OBJS) 
# single precision objects
$(ODIR)/libmyeis_sgl.a: $(SGL_OBJS) 
	$(AR) -rv  $(ODIR)/libmyeis_sgl.a  $(SGL_OBJS)
# both
$(ODIR)/libmyeis.a: $(SGL_OBJS)  $(DBL_OBJS) 
	$(AR) -rv  $(ODIR)/libmyeis.a  $(DBL_OBJS)  $(SGL_OBJS) 


#
# rules
#
$(ODIR)/%.sgl.o: %.sgl.F
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.dbl.o: %.dbl.F
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.dbl.o

$(ODIR)/test_dbl: $(ODIR)/test.dbl.o $(ODIR)/libmyeis.a
	$(F77) $(ODIR)/test.dbl.o $(ODIR)/libmyeis.a $(FFLAGS)  \
	 -o $(ODIR)/test_dbl \
	$(LDFLAGS)

$(ODIR)/test_sgl: $(ODIR)/test.sgl.o $(ODIR)/libmyeis.a
	$(F77) $(ODIR)/test.sgl.o $(ODIR)/libmyeis.a $(FFLAGS) \
	-o $(ODIR)/test_sgl \
	$(LDFLAGS)








