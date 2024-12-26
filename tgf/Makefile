ODIR = objects

#
# TGF accuracy settings for quadrature (see readme)
TGF_ACC_SETTINGS = -DTGF_ELHTRIAADAP_EPS=1e-12 -DTGF_ELHTRIAADAP_NQ=20 

FLAGS =

F90FLAGS = -Wall -O3 -cpp $(TGF_ACC_SETTINGS)

LIB_TGF_OBJS=$(ODIR)/mod_tgf.o $(ODIR)/el3dtria.o $(ODIR)/el3green.o $(ODIR)/el3greenh.o \
	$(ODIR)/elh3dtria.o $(ODIR)/legewhts.o $(ODIR)/mindlinpieces.o  \
	$(ODIR)/tria3adamar.o $(ODIR)/triagauc.o $(ODIR)/triahquad.o $(ODIR)/triasymq.o $(ODIR)/trilib.o \
	  $(ODIR)/el3dtria_omp.o $(ODIR)/elh3dtria_omp.o 
TGF_TEST_OBJS=test_driver/objects/dc3.o test_driver/objects/prini.o test_driver/objects/test_dc3_elh.o \
	test_driver/objects/test_dc3_elh_grid.o test_driver/objects/test_dc3_elt_grid.o test_driver/objects/test_tgf_main.o

all: tgf_library

tests:  tgf_main_test 

clean:
	rm -rf $(ODIR)/*.o $(ODIR)/*.a  


tgf_library: $(ODIR)/libtgf.a


$(ODIR)/libtgf.a: $(LIB_TGF_OBJS)
	$(AR) rv $(ODIR)/libtgf.a $(LIB_TGF_OBJS)

tgf_main_test:  $(ODIR)/libtgf.a  $(TGF_TEST_OBJS)
	$(F90) $(F90FLAGS) $(TGF_TEST_OBJS) -o tgf_main_test  $(ODIR)/libtgf.a

$(ODIR)/%.o:	source/%.f 
	$(F90) $(F90FLAGS) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

test_driver/objects/%.o:	test_driver/%.f 
	$(F90) $(F90FLAGS) -c $< -o test_driver/objects/$*.o
