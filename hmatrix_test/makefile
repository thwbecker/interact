MY_PRECISION = -DUSE_DOUBLE_PRECISION



# noise level for random interact version. this version 
# d
include ../makefile.gcc
include ../makefile.petsc
ifndef MPILD
MPILD = $(LD)
endif


ODIR = objects/

all: ex21



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

saved:
	cp $(BDIR)/$(INTERACT_BINARY_NAME) $(BDIR)/interact_saved;\
	cp $(BDIR)/$(INTERACT_BINARY_NAME).sgl $(BDIR)/interact_saved.sgl;\
	cp $(BDIR)/randomize_strike $(BDIR)/randomize_strike_saved
	if [ -s $(BDIR)/blockinvert_sph ];then \
		cp $(BDIR)/blockinvert_sph $(BDIR)/blockinvert_sph_saved;\
	fi;


noisefile:
	echo $(NOISELEVEL) > noise.dat

# individual programs

$(BDIR)/$(INTERACT_BINARY_NAME): $(OBJ) $(GEN_P_INC) $(LIBLIST) 
	$(MPILD) $(OBJ) -o $(BDIR)/$(INTERACT_BINARY_NAME) \
		$(PETSC_LIBS) $(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS)   $(LDFLAGS)

$(BDIR)/$(INTERACT_BINARY_NAME).sgl: $(OBJ_SGL) $(GEN_P_INC) $(LIBLIST_SGL) 
	$(MPILD) $(OBJ_SGL) -o $(BDIR)/$(INTERACT_BINARY_NAME).sgl \
		$(PETSC_LIBS) $(LIBS_SGL) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS)   $(LDFLAGS) 

$(BDIR)/interact_noise.$(NOISELEVEL): $(NOBJ) $(GEN_P_INC) $(LIBLIST) 
	$(MPILD)  $(NOBJ) -o $(BDIR)/interact_noise.$(NOISELEVEL) \
		$(PETSC_LIBS) $(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS)  $(LDFLAGS)


$(ODIR)/test_stuff: $(TOBJ) $(GEN_P_INC)  $(LIBLIST)  
	$(LD) $(LDFLAGS) $(TOBJ) -o $(BDIR)/test_stuff \
	$(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 

$(ODIR)/test_triangle_stress: $(T2OBJ) $(GEN_P_INC)  $(LIBLIST)  $(TRI_GREEN_OBJS)
	$(LD) $(LDFLAGS) $(T2OBJ) $(TRI_GREEN_OBJS) -o $(BDIR)/test_triangle_stress \
	$(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/randomflt: $(RANDOMFLT_OBJS)  $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(RANDOMFLT_OBJS) \
	-o $(BDIR)/randomflt $(LIBS)  $(LDFLAGS)

$(BDIR)/generate_random_2d: $(GENERATE_RANDOM_2D_OBJS)  $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(GENERATE_RANDOM_2D_OBJS) \
	-o $(BDIR)/generate_random_2d $(LIBS) $(LDFLAGS)


$(BDIR)/patchquad2patchtri: $(ODIR)/patchquad2patchtri.o $(GEN_P_INC) $(ODIR)/read_geometry.o \
	 $(LIBLIST) 
	$(MPILD) $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patchquad2patchtri.o \
		-o $(BDIR)/patchquad2patchtri  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)  $(LDFLAGS)  

$(BDIR)/patch2xyz: $(ODIR)/patch2xyz.o $(GEN_P_INC) $(ODIR)/read_geometry.o \
	 $(LIBLIST) 
	$(MPILD)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2xyz.o \
		-o $(BDIR)/patch2xyz  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)  $(LDFLAGS) 

$(BDIR)/patch2xyzvec: $(ODIR)/patch2xyzvec.o $(GEN_P_INC) \
	$(ODIR)/read_geometry.o $(LIBLIST) 
	$(MPILD)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2xyzvec.o \
		-o $(BDIR)/patch2xyzvec  $(LIBS) \
	$(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)  $(LDFLAGS)

$(BDIR)/patch2poly3d: $(ODIR)/patch2poly3d.o $(GEN_P_INC) $(ODIR)/read_geometry.o $(LIBLIST)
	$(MPILD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2poly3d.o \
		-o $(BDIR)/patch2poly3d \
		$(PETSC_LIBS) $(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 

$(BDIR)/patch2dis3d: $(ODIR)/patch2dis3d.o $(GEN_P_INC) \
	$(ODIR)/read_geometry.o $(LIBLIST) 
	$(MPILD) $(LDFLAGS)  $(ODIR)/read_geometry.o $(ODIR)/patch2dis3d.o \
		-o $(BDIR)/patch2dis3d  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/points2patch: $(ODIR)/points2patch.o  $(ODIR)/fit_plane.o $(ODIR)/libpatchio.a $(GEN_P_INC)  \
	 $(LIBLIST) 
	$(MPILD) $(ODIR)/points2patch.o $(ODIR)/libpatchio.a \
		$(ODIR)/fit_plane.o \
		-o $(BDIR)/points2patch  $(LIBS)  $(BLASLIB) $(LDFLAGS)   

$(BDIR)/create_random_stress_file: $(ODIR)/create_random_stress_file.o \
	$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(ODIR)/create_random_stress_file.o \
		-o $(BDIR)/create_random_stress_file  \
	$(LIBS)  $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)

$(BDIR)/calc_stress_stat: $(ODIR)/calc_stress_stat.o $(ODIR)/read_geometry.o \
	$(GEN_P_INC)  $(LIBLIST) $(ODIR)/calc_spatial_correlation.o
	$(MPILD)   $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/calc_spatial_correlation.o \
		$(ODIR)/calc_stress_stat.o -o $(BDIR)/calc_stress_stat  \
	$(LIBS)  $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)  

$(BDIR)/create_random_mu_file: $(ODIR)/create_random_mu_file.o \
	 $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(ODIR)/create_random_mu_file.o \
		 -o $(BDIR)/create_random_mu_file  \
	$(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)  

$(BDIR)/tri2patch: $(ODIR)/tri2patch.o $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)   $(ODIR)/tri2patch.o \
		-o $(BDIR)/tri2patch  $(LIBS)  $(LDFLAGS)


$(BDIR)/patch2geom: $(ODIR)/read_geometry.o $(ODIR)/patch2geom.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(MPILD)   $(ODIR)/read_geometry.o  $(ODIR)/patch2geom.o \
		-o $(BDIR)/patch2geom  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)

$(BDIR)/patch2vtk: $(ODIR)/read_geometry.o $(ODIR)/patch2vtk.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(MPILD) $(ODIR)/read_geometry.o  $(ODIR)/patch2vtk.o \
		-o $(BDIR)/patch2vtk  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)  $(LDFLAGS) 

$(BDIR)/patch2bc: $(ODIR)/patch2bc.o $(ODIR)/read_geometry.o $(GEN_P_INC)  \
	$(LIBLIST) 
	$(MPILD)  $(ODIR)/read_geometry.o \
		$(ODIR)/patch2bc.o -o $(BDIR)/patch2bc  $(LIBS) \
		 $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)	  $(LDFLAGS)

$(BDIR)/patch2corners:  $(ODIR)/read_geometry.o $(ODIR)/patch2corners.o \
	$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)    $(ODIR)/read_geometry.o $(ODIR)/patch2corners.o \
		-o $(BDIR)/patch2corners  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)

$(BDIR)/patch2group:   $(ODIR)/patch2group.o $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)   $(ODIR)/patch2group.o \
		-o $(BDIR)/patch2group  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)  $(LDFLAGS)

$(BDIR)/randomize_strike: $(RANDOMIZE_STRIKE_OBJS)  $(GEN_P_INC)
	$(MPILD)  $(RANDOMIZE_STRIKE_OBJS) \
	-o $(BDIR)/randomize_strike $(LIBS) $(LDFLAGS)


$(BDIR)/makefault: $(ODIR)/makefault.o   $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(ODIR)/makefault.o  \
	-o $(BDIR)/makefault $(LIBS)  $(LDFLAGS)

$(BDIR)/sort_events: $(ODIR)/sort_events.o $(ODIR)/myopen.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/sort_events.o  $(ODIR)/myopen.o \
	-o $(BDIR)/sort_events $(LIBS)

$(BDIR)/check_feedback: $(ODIR)/coulomb_stress.o $(ODIR)/check_feedback.o  $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(ODIR)/check_feedback.o \
	$(ODIR)/coulomb_stress.o \
		-o $(BDIR)/check_feedback  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS) 

$(BDIR)/mspectral: $(ODIR)/mspectral.o  interact.h $(ODIR)/myopen.o $(ODIR)/period.o
	$(LD) $(LDFLAGS)  $(ODIR)/mspectral.o $(ODIR)/myopen.o $(ODIR)/period.o \
		-o $(BDIR)/mspectral  $(LIBS)

$(BDIR)/calc_interaction_matrix: $(ODIR)/coulomb_stress.o \
	$(ODIR)/calc_interaction_matrix.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(MPILD)  $(ODIR)/coulomb_stress.o \
	$(ODIR)/calc_interaction_matrix.o  \
	-o $(BDIR)/calc_interaction_matrix $(LIBS) $(SUPERLU_LIBS) \
		$(PGLIBS) $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/calc_design_matrix: $(ODIR)/calc_design_matrix.o  \
	$(ODIR)/coulomb_stress.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(MPILD)  $(ODIR)/calc_design_matrix.o   \
	$(ODIR)/coulomb_stress.o \
	-o $(BDIR)/calc_design_matrix $(LIBS) $(SUPERLU_LIBS) \
		$(PGLIBS) $(SLATEC_LIBS) $(LDFLAGS)


$(BDIR)/test_sparse: $(ODIR)/test_sparse.o $(ODIR)/coulomb_stress.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/test_sparse.o $(ODIR)/coulomb_stress.o \
	-o $(BDIR)/test_sparse $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LIBS)

$(BDIR)/test_optimize: $(ODIR)/test_optimize.o $(ODIR)/optimize.o  \
	$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/test_optimize.o \
	 $(ODIR)/optimize.o  \
	-o $(BDIR)/test_optimize $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/test_solvers: $(ODIR)/test_solvers.o $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(LDFLAGS)  $(ODIR)/test_solvers.o \
	-o $(BDIR)/test_solvers $(PETSC_LIBS) $(LIBS) $(SUPERLU_LIBS) \
	$(SLATEC_LIBS) $(PGLIBS) 

$(BDIR)/ex_dense: $(ODIR)/ex_dense_v2.o $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(LDFLAGS)  $(ODIR)/ex_dense_v2.o \
	-o $(BDIR)/ex_dense $(PETSC_LIBS) $(LIBS) $(SUPERLU_LIBS) \
	$(SLATEC_LIBS) $(PGLIBS) 

$(BDIR)/project_stress: $(ODIR)/project_stress.o $(ODIR)/mysincos.o \
	$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(ODIR)/project_stress.o  $(ODIR)/llgeo.o \
	$(ODIR)/mysincos.o $(ODIR)/geometry.o \
	-o $(BDIR)/project_stress $(LIBS) $(LDFLAGS)

$(BDIR)/generate_slipdia: $(ODIR)/generate_slipdia.o $(ODIR)/mysincos.o \
	$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(ODIR)/generate_slipdia.o  $(ODIR)/llgeo.o \
	$(ODIR)/mysincos.o $(ODIR)/geometry.o \
	-o $(BDIR)/generate_slipdia $(LIBS) $(LDFLAGS)

$(BDIR)/calc_eigen_from_cart_stress: $(ODIR)/calc_eigen_from_cart_stress.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(ODIR)/calc_eigen_from_cart_stress.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/calc_eigen_from_cart_stress $(LIBS) $(EISPACK_LIB) \
		 $(LDFLAGS)

$(BDIR)/calc_cart_from_eigen_stress: $(ODIR)/calc_cart_from_eigen_stress.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(ODIR)/calc_cart_from_eigen_stress.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/calc_cart_from_eigen_stress $(LIBS) $(EISPACK_LIB) \
		$(LDFLAGS)

$(BDIR)/fit_simple_stress_from_cart: $(ODIR)/fit_simple_stress_from_cart.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)  $(ODIR)/fit_simple_stress_from_cart.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/fit_simple_stress_from_cart $(LIBS) $(EISPACK_LIB) \
		 $(LDFLAGS)

$(BDIR)/plotevents: $(ODIR)/plotevents.o $(ODIR)/plotting.o \
	$(ODIR)/plotting_palette.o $(ODIR)/mysincos.o $(ODIR)/llgeo.o \
	$(ODIR)/geometry.o $(GEN_P_INC)  $(LIBLIST) 
	$(MPILD) $(ODIR)/plotevents.o $(ODIR)/llgeo.o \
		$(ODIR)/geometry.o $(ODIR)/plotting_palette.o \
		$(ODIR)/myopen.o $(ODIR)/plotting.o \
		-o  $(BDIR)/plotevents  $(PGLIBS) $(LIBS) $(LDFLAGS)

$(BDIR)/read_bin_events: $(ODIR)/read_bin_events.o $(GEN_P_INC)  \
	 $(LIBLIST) 
	$(MPILD) $(ODIR)/read_bin_events.o \
		-o  $(BDIR)/read_bin_events $(LIBS) $(PGLIBS) \
		$(SUPERLU_LIBS) $(SLATEC_LIBS)  $(LDFLAGS) 


$(BDIR)/blockinvert_sph: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST) \
		$(BLOCKINVERT_SPH_OBJS)  $(ODIR)/blockinvert.sph.o  
	$(MPILD)  $(BLOCKINVERT_SPH_OBJS) $(ODIR)/blockinvert.sph.o  \
		$(MY_LIBDIR_SPEC)$(ODIR)/ $(MATHLIB)   $(GEOPROJECT_OBJS) \
		-o  $(BDIR)/blockinvert_sph  $(LIBS)		\
		$(PETSC_LIBS)	$(GEOPROJECT_LIBS) \
		$(DEBUG_LIBS) $(EISPACK_LIB) $(PGLIBS) 	\
		$(COMPUTATIONAL_LIBS)  $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/fstress2hor: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o
	$(MPILD)  $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2hor   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/fit_mean_stress: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o
	$(MPILD) $(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fit_mean_stress   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)   $(DEBUG_LIBS) $(EISPACK_LIB)  $(LDFLAGS) 

$(BDIR)/fstress2eig: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o
	$(MPILD)  $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2eig   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/block_evaluate_solution: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST)	\
		$(BLOCK_EVALUATE_SOLUTION_OBJS) 		\
		 $(ODIR)/block_evaluate_solution.o
	$(MPILD)  $(BLOCK_EVALUATE_SOLUTION_OBJS)		\
		$(ODIR)/block_evaluate_solution.o 	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/block_evaluate_solution   -linput -lpatchio	\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/block_checkflt: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST) \
		$(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o
	$(MPILD) $(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)			\
		-o  $(BDIR)/block_checkflt   -linput -lpatchio		\
		$(GEOPROJECT_LIBS) $(PETSC_LIBS)				\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB)  $(PGLIBS) $(LDFLAGS) 

$(BDIR)/geo_okada: $(ODIR)/geo_okada.o $(ODIR)/coulomb_stress.o $(GEN_P_INC)  \
	$(GEOPROJECT_OBJS) $(LIBLIST) 
	$(MPILD)  $(ODIR)/geo_okada.o $(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS) \
		-o  $(BDIR)/geo_okada  $(ODIR)/coulomb_stress.o -lpatchio \
		$(GEOPROJECT_LIBS)	$(COMPUTATIONAL_LIBS) $(DEBUG_LIBS) $(MATHLIB)  $(LDFLAGS)

#
# C function prototyper
#

proto: 	auto_proto.h auto_proto.sgl.h

auto_proto.h: 
	rm -f auto_proto.h 2> /dev/null;\
	touch auto_proto.h;\
	cproto  $(DEFINE_FLAGS)  $(GEOPROJECT_INCLUDES) \
		$(PGPLOT_DEFINES) $(PGPLOT_INCLUDES) $(MY_PRECISION) \
		$(SLATEC_INCLUDES)  $(SUPERLU_INCLUDES)  -f2 -q *.c 2> /dev/null | grep -v "void main("  | grep -v "int main(" > tmp.h; \
		mv tmp.h auto_proto.h



auto_proto.sgl.h: 
	rm -f auto_proto.sgl.h 2> /dev/null;\
	touch auto_proto.sgl.h;\
	cproto  $(DEFINE_FLAGS)  $(GEOPROJECT_INCLUDES) \
		$(PGPLOT_DEFINES) $(PGPLOT_INCLUDES)  \
		$(SLATEC_INCLUDES)  $(SUPERLU_INCLUDES)  -f2 -q *.c 2> /dev/null  | \
		grep -v "void main("  | grep -v "int main(" > tmp.h;\
	mv tmp.h auto_proto.sgl.h


#
# libraries
#
libraries: $(LIBLIST)

$(EISPACK_DIR)/libmyeis.a:
	cd $(EISPACK_DIR); make ; cd -

$(SLATEC_DIR)/libslatec_dbl.a:
	cd $(SLATEC_DIR); make ; cd -

$(SLATEC_DIR)/libslatec_sgl.a:
	cd $(SLATEC_DIR); make ; cd -

$(ODIR)/libpatchio.a: $(PATCH_IO_OBJS)
	$(AR) rv $(ODIR)/libpatchio.a $(PATCH_IO_OBJS)

$(ODIR)/libinput.a: $(INPUT_OBJS)
	$(AR) rv $(ODIR)/libinput.a $(INPUT_OBJS)

$(ODIR)/libpatchio.sgl.a: $(PATCH_IO_OBJS_SGL)
	$(AR) rv $(ODIR)/libpatchio.sgl.a $(PATCH_IO_OBJS_SGL)

$(ODIR)/libinput.sgl.a: $(INPUT_OBJS_SGL)
	$(AR) rv $(ODIR)/libinput.sgl.a $(INPUT_OBJS_SGL)


#
# source code
#
#
#
# files with special compiler options
#
#
$(ODIR)/numrec_svd_routines.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F $(MY_PRECISION) \
	-o $(ODIR)/numrec_svd_routines.o
$(ODIR)/numrec_svd_routines.sgl.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F -o $(ODIR)/numrec_svd_routines.sgl.o
#
# other specialities
#
$(ODIR)/solve_mode_dependend_1.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1 $(MY_PRECISION) -o  $(ODIR)/solve_mode_dependend_1.o
$(ODIR)/solve_mode_dependend_2.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2 $(MY_PRECISION) -o  $(ODIR)/solve_mode_dependend_2.o
$(ODIR)/solve_mode_dependend_3.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_3 $(MY_PRECISION) -o  $(ODIR)/solve_mode_dependend_3.o
$(ODIR)/solve_mode_dependend_4.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4 $(MY_PRECISION) -o  $(ODIR)/solve_mode_dependend_4.o
#
# single prec (leave out my_precision)
$(ODIR)/solve_mode_dependend_1.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1  -o  $(ODIR)/solve_mode_dependend_1.sgl.o
$(ODIR)/solve_mode_dependend_2.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2  -o  $(ODIR)/solve_mode_dependend_2.sgl.o
$(ODIR)/solve_mode_dependend_3.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_3  -o  $(ODIR)/solve_mode_dependend_3.sgl.o
$(ODIR)/solve_mode_dependend_4.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4  -o  $(ODIR)/solve_mode_dependend_4.sgl.o
#
# spherical versions of the blockinvert code
$(ODIR)/blockinvert.sph.o:	blockinvert.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL -c $< $(MY_PRECISION) \
	-o  $(ODIR)/blockinvert.sph.o
$(ODIR)/block_matrix.sph.o:	block_matrix.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_matrix.sph.o
$(ODIR)/read_stress_observations.sph.o:	read_stress_observations.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/read_stress_observations.sph.o
$(ODIR)/block_read_gps.sph.o:	block_read_gps.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_gps.sph.o
$(ODIR)/block_output.sph.o:	block_output.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_output.sph.o
$(ODIR)/block_solve.sph.o:	block_solve.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_solve.sph.o
$(ODIR)/block_read_bflt.sph.o:	block_read_bflt.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_bflt.sph.o
$(ODIR)/block_stress.sph.o:	block_stress.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_stress.sph.o
$(ODIR)/block_read_euler.sph.o:	block_read_euler.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_euler.sph.o
#
# single
$(ODIR)/blockinvert.sph.sgl.o:	blockinvert.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL -c $<  -o  $(ODIR)/blockinvert.sph.sgl.o
$(ODIR)/block_matrix.sph.sgl.o:	block_matrix.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_matrix.sph.sgl.o
$(ODIR)/read_stress_observations.sph.sgl.o:	read_stress_observations.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/read_stress_observations.sph.sgl.o
$(ODIR)/block_read_gps.sph.sgl.o:	block_read_gps.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_gps.sph.sgl.o
$(ODIR)/block_output.sph.sgl.o:	block_output.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_output.sph.sgl.o
$(ODIR)/block_solve.sph.sgl.o:	block_solve.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_solve.sph.sgl.o
$(ODIR)/block_read_bflt.sph.sgl.o:	block_read_bflt.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_bflt.sph.sgl.o
$(ODIR)/block_stress.sph.sgl.o:	block_stress.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_stress.sph.sgl.o
$(ODIR)/block_read_euler.sph.sgl.o:	block_read_euler.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_euler.sph.sgl.o

#
# the random noise version of coulomb stress
$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o: coulomb_stress.c $(GEN_P_INC) noise.dat
	$(CC) $(CFLAGS) -c coulomb_stress.c \
	-DADD_COULOMB_STRESS_NOISE=$(NOISELEVEL) \
	$(MY_PRECISION) -o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o

$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o: coulomb_stress.c $(GEN_P_INC) noise.dat
	$(CC) $(CFLAGS) -c coulomb_stress.c \
	-DADD_COULOMB_STRESS_NOISE=$(NOISELEVEL) \
	-o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o
# 

$(ODIR)/terminate.o:	terminate.c $(GEN_P_INC)
	$(MPICC) $(CFLAGS)  $(MY_PRECISION) -c terminate.c -o $(ODIR)/terminate.o

$(ODIR)/test_solvers.o:	test_solvers.c $(GEN_P_INC)
	$(MPICC) $(CFLAGS)  $(MY_PRECISION) -c test_solvers.c -o $(ODIR)/test_solvers.o

$(ODIR)/ex_dense.o:	ex_dense.c $(GEN_P_INC)
	$(MPICC) $(CFLAGS)  $(MY_PRECISION) -c ex_dense.c -o $(ODIR)/ex_dense.o

$(ODIR)/ex_dense_v2.o:	ex_dense_v2.c $(GEN_P_INC)
	$(MPICC) $(CFLAGS)  $(MY_PRECISION) -c ex_dense_v2.c -o $(ODIR)/ex_dense_v2.o



#
# some generic rules with normal dependencies
#
$(ODIR)/%.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f90 $(GEN_P_INC)
	$(F90)   $(F90FLAGS) $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

# single prec versions
$(ODIR)/%.sgl.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) -c $<  -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.f90 $(GEN_P_INC)
	$(F90)  $(F90FLAGS) -c $<  -o $(ODIR)/$*.sgl.o
