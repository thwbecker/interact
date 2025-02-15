################################################################################
#
#
#  (C) Thorsten Becker, thwbecker@post.harvard.edu, 1999-2025
#
#  makefile for interact and related programs
#
#
#################################################################################
#
# read through this makefile to make changes adjusting to your environment
# (especially the PGPLOT support has to be adapted or commented out. so far
# the sun makefiles have no PGPLOT, the SGI and LINUX files have PGPLOT support). 
#
# Then type "make" and "interact -h" after succesful compilation for the help page.
#
# 
# for calc_eigen_from_cart_stress you need EISPACK in the modified
# form as supplied in the myeispack.tar.gz file.  Extract this files
# to the location as specified at the end of the makefile.gcc (by
# default eispack/) and run "make" in that directory to create the
# modified EISPACK libraries first
#
#
# WARNING:
#
# If you change the flags, type "make clean" and "make all". This will create a new
# C function prototype file and recompile everything, a necessary step.
#
# A list of flags follows:
#
#  -DUSE_DOUBLE_PRECISION    for double precision, if left 
#                            away, uses single prec.
#                            
#  -DBINARY_PATCH_EVENT_FILE for binary output of single 
#                            event files in interact (smaller file)
#
#  -DTWO_DIM	             limit the analysis to a 2-D plane strain case
#
#  -DNO_OPENING_MODES        for restricted operation in 
#                            strike and dip direction only (saves memory)
#
#  -DLATENCY=0.5             define this to allow for a 
#                            certain period of time before patches 
#                            can be reactivated (didn't do much good)
#
#  -DNO_COHESION             remove all references to a non-zero 
#                            cohesion (improved speed)
#
#  -DCHECK_CI_ONE_WAY        normally, randomize_strike and other programs will
#                            check if the self-interaction is smaller than the Coulomb
#			     stress change on other patches for positive 
#    			     feedback loops for pairs of
#                            patches. if this switch is set, will also remove the patch
#                            when condition is only met one way
#
#  -DALLOW_NON_3DQUAD_GEOM   allows for different fault element types, so far rectangular and 
#                            point source are implemented. if all faults are rectangular, leave
#		  	     away for improved speed 
#
#
#  -DUSE_NUMREC_SVD          use SVD decomposition and backsubstitution from Numerical
#                            recipes. Default is from LAPACK since numerical reps
#                            didn't work for -O:fast=ip27 compilation on IRIX systems
#                            If this is not set, you will have to specify the flags
#                            for the inclusion of LAPACK in the makefile.gcc files.
#
#  -DUSE_GEOPROJECT          use geoproject and produce code that allows for I/O - Flag should be set in makefile.geoproject
#                            using geographic coordinates and projections. for this to
#                            work, you will have to have GMT (www.gmt.soest.hawaii.edu)
#                            installed - THIS USES GMT < VERSION 5 - this is set in makefile.geoproject
# 			     which is commented out by default
# 		             
#
#
#  -DUSE_PETSC               use the Petsc libraries and run the code in parallel, where only matrix inversions are parallelized right now
#                            comments out makefile.petsc if no petsc support required. Flag should be set in makefile.petsc
# 			     possible LU solvers are
# 				-pc_factor_mat_solver_type scalapack -mat_type scalapack
#  				-pc_factor_mat_solver_type elemental -mat_type elemental
# 		             iterative solver
# 		 	         -ksp_type fgmres -pc_type none -ksp_max_it 10000 -ksp_rtol 1.0e-8 
# 		 	         -ksp_type fgmres -pc_type jacobi -ksp_max_it 10000 -ksp_rtol 1.0e-8
#                            for this to work, you will have to have  $(PETSC_DIR) and $(PETSC_ARCH) defined
#
# to run in parallel for example
#
#  mpirun -np 8 interact -pc_factor_mat_solver_type scalapack -mat_type scalapack
# to debug:
# mpirun -np 1 valgrind --tool=memcheck -q --num-callers=20 --log-file=pmia.%p.log interact -fpetsc -malloc=off
#
# ____________________________________________________________________________________________________
#
# 		              NOTE: CERTAIN PROGRAMS USE petsc_settings.yaml FOR DEFAULTS 
# ____________________________________________________________________________________________________

#  Example settings:
#
MY_PRECISION = -DUSE_DOUBLE_PRECISION
#
# if we want to test the HBI version of the triangular stress routine
# -DUSE_HBI_TDDEF
#
# if we want to use the TGF triangle (NOT TESTED yet)
# -DUSE_TGF_TRIANGLE
#TGF_LIB = tgf/objects/libtgf.a
TGF_LIB = 
#
COMMON_DEFINES =  -DBINARY_PATCH_EVENT_FILE -DCHECK_CI_ONE_WAY  \
	-DALLOW_NON_3DQUAD_GEOM 

#
# noise level for random interact version. this version 
# doesn't get compiled automatically.
NOISELEVEL=1e-08
#
MAIN_DEFINES = $(COMMON_DEFINES) 

#
#
# directory for object files
ODIR = objects
# directory for binaries
BDIR = bin
#
# choice of Okada routine
OKROUTINE = $(ODIR)/dc3d.o	# my modified version
OKROUTINE_DEBUG = $(ODIR)/dc3d.dbg.o	# my modified version
#OKROUTINE = $(ODIR)/dc3d_original_double.o	# original, changed to double

# 
# include the machine dependent flags
# 
include makefile.gcc
#include makefile.icc
#include makefile.mixed_mkl
#include makefile.mixed
#include makefile.icc_frontera
#
# add this for pgplot support, otherwise comment it out
# you will use runtime plotting capabilities
#include makefile.pgplot
#
# add this for slatec NNLS routine support, otherwise comment it out
# you will use NNLS solving capabilities
#include makefile.slatec

#
# petsc, will override some of the flags
# comment out if not needed
include makefile.petsc
ifndef MPILD
MPILD = $(LD)
endif

#
# add this for superlu support, otherwise comment it out
# you will loose sparse matrix SuperLU LU solver capabilities
#include makefile.superlu
#
# add this for geoprojection support
#
include makefile.geoproject

# add up all define flags
DEFINE_FLAGS = $(MAIN_DEFINES) $(SLATEC_DEFINES) $(PETSC_DEFINES) \
	$(PGPLOT_DEFINES) $(SUPERLU_DEFINES)  \
	$(GEOPROJECT_DEFINES)
# defines and pgplot flags
FLAGS = $(DEFINE_FLAGS) $(PGPLOT_INCLUDES) $(SLATEC_INCLUDES) \
	$(SUPERLU_INCLUDES) $(GEOPROJECT_INCLUDES) $(PETSC_INCLUDES) 
# C and FORTRAN compiler specific flags, add them to the 


# other flags
CFLAGS =   $(FLAGS) $(SCARGS)   $(MACHINE_DEFINES)
FFLAGS =   $(FLAGS) $(SFARGS)   $(MACHINE_DEFINES)
F90FLAGS = $(FLAGS) $(SF90ARGS) $(MACHINE_DEFINES)


#
# list of object files
#
# solve mode dependend, ie. one for each mode of accesing the
# interaction matrix (improved speed, one would hope)
SMD_OBJS = $(ODIR)/solve_mode_dependend_1.o	\
	$(ODIR)/solve_mode_dependend_2.o	\
	$(ODIR)/solve_mode_dependend_3.o	\
	$(ODIR)/solve_mode_dependend_4.o 

SMD_OBJS_DEBUG = $(SMD_OBJS:.o=.dbg.o)
#
# matrix solver 
#
MATRIX_SOLVER_OBJS = $(ODIR)/numrec_svd_routines.o $(ODIR)/nnls_lawson.o	\
	$(ODIR)/nnls.o $(ODIR)/svd.o $(ODIR)/solve.o		\
	$(ODIR)/lusolve.o $(ODIR)/sparse_solve.o \
	$(ODIR)/ilaenv_wrapper.o

MATRIX_SOLVER_OBJS_DEBUG = $(MATRIX_SOLVER_OBJS:.o=.dbg.o)

#
# triangular dislocation routines
#

TRI_GREEN_OBJS = $(ODIR)/tdstresshs.o  $(ODIR)/tddisphs.o   $(ODIR)/tdstress_hbi.o 

TRI_GREEN_OBJS_DEBUG = $(TRI_GREEN_OBJS:.o=.dbg.o)

# list of objects for patch i/o. these also include
# all object files that deal with interaction coefficients
# since we might want to check for fatal interaction 
# setups. 
#
PATCH_IO_OBJS = $(ODIR)/divide_fault_in_patches.o $(ODIR)/sparse.o	\
	$(ODIR)/sparse_nr.o $(ODIR)/randomize_list.o			\
	$(ODIR)/print_patch_geometry.o $(SMD_OBJS) $(ODIR)/llgeo.o	\
	$(ODIR)/geometry.o $(ODIR)/fltcopy.o		\
	$(ODIR)/check_interaction.o $(ODIR)/eval_2dsegment.o		\
	$(ODIR)/eval_okada.o $(ODIR)/tdd_coeff.o $(ODIR)/rhs.o		\
	$(ODIR)/eval_green.o $(ODIR)/eval_triangle_nw.o			\
	$(ODIR)/eval_iquad.o $(ODIR)/eval_triangle_tgf.o			\
	$(ODIR)/interact.o   $(ODIR)/mysincos.o 	\
	$(OKROUTINE) $(ODIR)/fracture_criterion.o			\
	$(ODIR)/myopen.o  $(ODIR)/randgen.o \
	$(ODIR)/string_compare.o  $(TRI_GREEN_OBJS)

PATCH_IO_OBJS_DEBUG = $(PATCH_IO_OBJS:.o=.dbg.o)
PATCH_IO_OBJS_SGL = $(PATCH_IO_OBJS:.o=.sgl.o)
#
#
#
# objects for the main input routines of interact
# by some chance, the solver routines are referenced, too
#
# will go to libinput.a
#
INPUT_OBJS = $(ODIR)/read_boundary_conditions.o $(ODIR)/read_fltdat.o \
	$(ODIR)/read_geometry.o $(ODIR)/init.o $(ODIR)/help_and_comments.o	\
	$(ODIR)/output.o $(ODIR)/input.o $(ODIR)/matrixio.o			\
	$(ODIR)/calc_spatial_correlation.o	$(ODIR)/stress_aux.o		\
	$(ODIR)/quake.o  $(ODIR)/restart.o			\
	$(ODIR)/coulomb_stress.o $(POBJS) 
INPUT_OBJS_SGL = $(INPUT_OBJS:.o=.sgl.o)
INPUT_OBJS_DEBUG = $(INPUT_OBJS:.o=.dbg.o)

#
# list of objects for the main program, interact
INTERACT_OBJS = $(ODIR)/rupture.o $(ODIR)/adjust_time_step.o  $(ODIR)/terminate.o \
	$(ODIR)/calc_stress.o $(MATRIX_SOLVER_OBJS)

INTERACT_OBJS_DEBUG = $(ODIR)/rupture.dbg.o $(ODIR)/adjust_time_step.dbg.o  $(ODIR)/terminate.dbg.o \
	$(ODIR)/calc_stress.dbg.o $(MATRIX_SOLVER_OBJS_DEBUG)


# this is a random noise added version
INTERACT_NOISE_OBJS = $(ODIR)/rupture.o	$(ODIR)/calc_stress.o \
	$(ODIR)/adjust_time_step.o $(MATRIX_SOLVER_OBJS)	\
	$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o $(ODIR)/terminate.o
#
# objects for randomflt
RANDOMFLT_OBJS = $(ODIR)/randomflt.o $(ODIR)/compare_fault.o	\
	$(ODIR)/tritri.o $(ODIR)/far_enough.o	\
	$(ODIR)/coulomb_stress.o
#
# objects for generate_random_2d
GENERATE_RANDOM_2D_OBJS = $(ODIR)/generate_random_2d.o $(ODIR)/compare_fault.o	\
	$(ODIR)/tritri.o $(ODIR)/far_enough.o	\
	$(ODIR)/coulomb_stress.o

#
# objects for randomize_strike
RANDOMIZE_STRIKE_OBJS = $(ODIR)/randomize_strike.o		\
	$(ODIR)/tritri.o $(ODIR)/far_enough.o	\
	$(ODIR)/compare_fault.o $(ODIR)/coulomb_stress.o
#
# objects for the GPS velocity/crustal block inversion program
# blockinvert
BLOCKINVERT_OBJS = $(ODIR)/block_read_bflt.o $(ODIR)/coulomb_stress.o	\
	$(ODIR)/block_matrix.o $(ODIR)/block_read_euler.o	\
	$(ODIR)/read_stress_observations.o $(ODIR)/block_read_gps.o $(ODIR)/svd.o \
	$(ODIR)/eigensystem.o $(ODIR)/block_solve.o $(ODIR)/solve.o $(ODIR)/nnls.o		\
	$(ODIR)/levmarq_numrec.o $(ODIR)/block_output.o  $(ODIR)/numrec_svd_routines.o \
	$(ODIR)/block_stress.o $(ODIR)/block_levmarq.o 

BLOCKINVERT_SPH_OBJS = $(ODIR)/block_read_bflt.sph.o \
	$(ODIR)/coulomb_stress.o $(ODIR)/solve.o  $(ODIR)/nnls.o  $(ODIR)/svd.o	 \
	$(ODIR)/block_read_euler.sph.o  $(ODIR)/lusolve.o \
	$(ODIR)/block_matrix.sph.o $(ODIR)/numrec_svd_routines.o	\
	$(ODIR)/block_read_gps.sph.o $(ODIR)/nnls_lawson.o \
	$(ODIR)/read_stress_observations.sph.o	\
	$(ODIR)/block_solve.sph.o 	$(ODIR)/ilaenv_wrapper.o	\
	$(ODIR)/eigensystem.o 	 \
	$(ODIR)/levmarq_numrec.o $(ODIR)/block_output.sph.o  \
	$(ODIR)/block_stress.sph.o $(ODIR)/block_levmarq.o 

FSTRESS2HOR_OBJS = $(ODIR)/block_read_gps.sph.o $(ODIR)/svd.o $(ODIR)/numrec_svd_routines.o \
	$(ODIR)/block_output.o $(ODIR)/block_matrix.o	$(ODIR)/ilaenv_wrapper.o \
	$(ODIR)/eigensystem.o $(ODIR)/block_read_bflt.o  \
	$(ODIR)/block_stress.o $(ODIR)/block_solve.o	

#
# real interact and test program
OBJ = $(ODIR)/regular_interact_main.o $(INTERACT_OBJS) 
OBJ_SGL = $(OBJ:.o=.sgl.o)
OBJ_DEBUG = $(OBJ:.o=.dbg.o)
NOBJ = $(ODIR)/regular_interact_main.o $(INTERACT_NOISE_OBJS)
# test programs
TOBJ = $(ODIR)/test_stuff.o $(INTERACT_OBJS)
T2OBJ = $(ODIR)/test_triangle_stress.o $(INTERACT_OBJS) 
#

#
# 
# include dependencies for all source codes
#
GEN_P_INC = interact.h precision_single.h precision_double.h \
	structures.h macros.h auto_proto.h auto_proto.sgl.h fortran_proto.h \
	filenames.h properties.h blockinvert.h constants.h

LIBLIST = $(ODIR)/libpatchio.a $(ODIR)/libinput.a $(EISPACK_DIR)/libmyeis.a $(TGF_LIB)
LIBLIST_DEBUG = $(ODIR)/libpatchio.dbg.a $(ODIR)/libinput.dbg.a $(EISPACK_DIR)/libmyeis.a $(TGF_LIB)
LIBLIST_SGL = 	$(ODIR)/libpatchio.sgl.a $(ODIR)/libinput.sgl.a

#
# libraries, also linker flags
#
LIBS = $(MY_LIBDIR_SPEC)$(ODIR)/    -linput -lpatchio  $(TGF_LIB) \
	$(COMPUTATIONAL_LIBS) $(MATHLIB) 

LIBS_SGL = $(MY_LIBDIR_SPEC)$(ODIR)/     -linput.sgl -lpatchio.sgl $(TGF_LIB) \
	$(COMPUTATIONAL_LIBS) $(MATHLIB) 

LIBS_DEBUG = $(MY_LIBDIR_SPEC)$(ODIR)/     -linput.dbg -lpatchio.dbg $(TGF_LIB) \
	$(COMPUTATIONAL_LIBS) $(MATHLIB) 

#
# list of all programs in groups
#

all: obj_directories libraries main_prog \
	tools converters geom_converters

really_all: obj_directories debug_libraries libraries main_prog \
	tools converters geom_converters $(BDIR)/compress_interaction_matrix.dbg \
	inoise analysis geographic_tools $(BDIR)/$(INTERACT_BINARY_NAME).dbg
#	pgplot_progs 

main_prog: $(BDIR)/$(INTERACT_BINARY_NAME) $(BDIR)/$(INTERACT_BINARY_NAME).sgl $(BDIR)/rsf_solve

inoise: noisefile $(BDIR)/interact_noise.$(NOISELEVEL)

tools: misc_tools random_geom_tools random_prop_tools

misc_tools: $(BDIR)/makefault $(BDIR)/calc_interaction_matrix $(BDIR)/calc_design_matrix \
		$(BDIR)/project_stress $(BDIR)/calc_eigen_from_cart_stress   \
	$(BDIR)/check_feedback $(BDIR)/fit_simple_stress_from_cart  $(BDIR)/petsc_simple_solve \
	$(BDIR)/calc_cart_from_eigen_stress $(BDIR)/compress_interaction_matrix \
	$(BDIR)/sort_events $(BDIR)/generate_slipdia

random_geom_tools:  $(BDIR)/randomflt  $(BDIR)/generate_random_2d \
	$(BDIR)/randomize_strike 

random_prop_tools: $(BDIR)/create_random_stress_file \
	$(BDIR)/create_random_mu_file $(BDIR)/calc_stress_stat


geographic_tools: $(BDIR)/blockinvert_sph $(BDIR)/geo_okada  \
	$(BDIR)/block_checkflt $(BDIR)/block_evaluate_solution \
	$(BDIR)/fstress2hor	$(BDIR)/fstress2eig \
	$(BDIR)/fit_mean_stress

analysis: $(BDIR)/mspectral
#$(BDIR)/patch2geom : geomview, outdated

converters: $(BDIR)/patch2xyz   $(BDIR)/patch2vtk $(BDIR)/patch2bc \
	$(BDIR)/patch2vertices $(BDIR)/patch2group  \
	$(BDIR)/patch2xyzvec # $(BDIR)/patch2poly3d $(BDIR)/patch2dis3d

geom_converters: $(BDIR)/points2patch $(BDIR)/tri2patch  $(BDIR)/patchquad2patchtri 

test:    $(ODIR)/test_stuff  $(ODIR)/test_triangle_stress

matrix_test_progs: $(BDIR)/test_sparse $(BDIR)/test_optimize $(BDIR)/test_solvers \
	$(BDIR)/ex_dense

pgplot_progs:  $(BDIR)/plotevents $(BDIR)/read_bin_events 

#
# some special targets
#

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

$(BDIR)/$(INTERACT_BINARY_NAME).dbg: $(OBJ_DEBUG) $(GEN_P_INC) $(LIBLIST_DEBUG) 
	$(MPILD) $(OBJ_DEBUG) -o $(BDIR)/$(INTERACT_BINARY_NAME).dbg \
		$(PETSC_LIBS) $(LIBS_DEBUG) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS)   $(LDFLAGS) 

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

$(BDIR)/patch2vertices:  $(ODIR)/read_geometry.o $(ODIR)/patch2vertices.o \
	$(GEN_P_INC)  $(LIBLIST) 
	$(MPILD)    $(ODIR)/read_geometry.o $(ODIR)/patch2vertices.o \
		-o $(BDIR)/patch2vertices  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LDFLAGS)

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

$(BDIR)/compress_interaction_matrix: $(ODIR)/compress_interaction_matrix.o $(ODIR)/coulomb_stress.o \
	$(ODIR)/interact.o $(GEN_P_INC) $(LIBLIST) 
	$(MPILD)   $(ODIR)/compress_interaction_matrix.o    $(ODIR)/petsc_interact.o $(ODIR)/coulomb_stress.o $(ODIR)/interact.o  \
	-o $(BDIR)/compress_interaction_matrix $(LIBS) $(SUPERLU_LIBS) \
			$(PETSC_LIBS) $(PGLIBS) $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/petsc_simple_solve: $(ODIR)/petsc_simple_solve.o $(ODIR)/coulomb_stress.o $(ODIR)/interact.o  $(ODIR)/petsc_interact.o \
	$(GEN_P_INC) $(LIBLIST) 
	$(MPILD)   $(ODIR)/petsc_simple_solve.o $(ODIR)/petsc_interact.o  $(ODIR)/coulomb_stress.o $(ODIR)/interact.o  \
	-o $(BDIR)/petsc_simple_solve $(LIBS) $(SUPERLU_LIBS) \
			$(PETSC_LIBS) $(PGLIBS) $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/rsf_solve: $(ODIR)/rsf_solve.o $(ODIR)/coulomb_stress.o $(ODIR)/interact.o  $(ODIR)/petsc_interact.o \
	$(GEN_P_INC) $(LIBLIST) 
	$(MPILD)   $(ODIR)/rsf_solve.o $(ODIR)/petsc_interact.o  $(ODIR)/coulomb_stress.o $(ODIR)/interact.o  \
	-o $(BDIR)/rsf_solve $(LIBS) $(SUPERLU_LIBS) \
			$(PETSC_LIBS) $(PGLIBS) $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/compress_interaction_matrix.dbg: $(ODIR)/compress_interaction_matrix.dbg.o $(ODIR)/coulomb_stress.dbg.o  \
	$(ODIR)/interact.dbg.o $(ODIR)/petsc_interact.dbg.o $(GEN_P_INC) $(LIBLIST) 
	$(MPILD)   $(ODIR)/compress_interaction_matrix.dbg.o  $(ODIR)/petsc_interact.dbg.o \
	$(ODIR)/coulomb_stress.dbg.o $(ODIR)/interact.dbg.o \
	-o $(BDIR)/compress_interaction_matrix.dbg $(LIBS) $(SUPERLU_LIBS) \
			$(PETSC_LIBS) $(PGLIBS) $(SLATEC_LIBS)  $(LDFLAGS)

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
		 $(EISPACK_LIB) $(PGLIBS) 	\
		$(COMPUTATIONAL_LIBS)  $(SLATEC_LIBS)  $(LDFLAGS)

$(BDIR)/fstress2hor: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o
	$(MPILD)  $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2hor   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)   $(TGF_LIB) $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/fit_mean_stress: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o
	$(MPILD) $(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fit_mean_stress   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)	$(TGF_LIB)				\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)    $(EISPACK_LIB)  $(LDFLAGS) 

$(BDIR)/fstress2eig: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o
	$(MPILD)  $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2eig   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)   $(TGF_LIB) $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/block_evaluate_solution: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST)	\
		$(BLOCK_EVALUATE_SOLUTION_OBJS) 		\
		 $(ODIR)/block_evaluate_solution.o
	$(MPILD)  $(BLOCK_EVALUATE_SOLUTION_OBJS)		\
		$(ODIR)/block_evaluate_solution.o 	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/block_evaluate_solution   -linput -lpatchio	\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)  $(TGF_LIB)  $(EISPACK_LIB) $(LDFLAGS)

$(BDIR)/block_checkflt: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST) \
		$(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o
	$(MPILD) $(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)			\
		-o  $(BDIR)/block_checkflt   -linput -lpatchio		\
		$(GEOPROJECT_LIBS) $(PETSC_LIBS)				\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)    $(TGF_LIB) $(EISPACK_LIB)  $(PGLIBS) $(LDFLAGS) 

$(BDIR)/geo_okada: $(ODIR)/geo_okada.o $(ODIR)/coulomb_stress.o $(GEN_P_INC)  \
	$(GEOPROJECT_OBJS) $(LIBLIST) 
	$(MPILD)  $(ODIR)/geo_okada.o $(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS) \
		-o  $(BDIR)/geo_okada  $(ODIR)/coulomb_stress.o -lpatchio \
		$(GEOPROJECT_LIBS)	$(COMPUTATIONAL_LIBS)  $(MATHLIB)  $(LDFLAGS)

#
# C function prototyper
#

proto: 	auto_proto.h auto_proto.sgl.h

auto_proto.h: 
	rm -f auto_proto.h 2> /dev/null;\
	touch auto_proto.h;\
	cproto  $(DEFINE_FLAGS)  $(GEOPROJECT_INCLUDES) \
		$(PGPLOT_DEFINES) $(PGPLOT_INCLUDES) $(MY_PRECISION) \
		$(SLATEC_INCLUDES)  $(SUPERLU_INCLUDES)  -f2 -q *.c 2> /dev/null | \
		grep -v "void main("  | grep -v "int main(" > tmp.h; \
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
debug_libraries: $(LIBLIST_DEBUG)

tgf/objects/libtgf.a:
	cd tgf; make ; cd -

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

$(ODIR)/libpatchio.dbg.a: $(PATCH_IO_OBJS_DEBUG)
	$(AR) rv $(ODIR)/libpatchio.dbg.a $(PATCH_IO_OBJS_DEBUG)

$(ODIR)/libinput.dbg.a: $(INPUT_OBJS_DEBUG)
	$(AR) rv $(ODIR)/libinput.dbg.a $(INPUT_OBJS_DEBUG)

#
# source code
#
#
#
# files with special compiler options
#
#
$(ODIR)/numrec_svd_routines.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F $(MY_PRECISION) $(OPTIM_FLAGS) \
	-o $(ODIR)/numrec_svd_routines.o

$(ODIR)/numrec_svd_routines.dbg.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F $(MY_PRECISION) $(DEBUG_FLAGS) \
	-o $(ODIR)/numrec_svd_routines.dbg.o

$(ODIR)/numrec_svd_routines.sgl.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) $(OPTIM_FLAGS) numrec_svd_routines.F -o $(ODIR)/numrec_svd_routines.sgl.o

#
# other specialities
#
$(ODIR)/solve_mode_dependend_1.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1 $(MY_PRECISION) $(OPTIM_FLAGS)  -o  $(ODIR)/solve_mode_dependend_1.o
$(ODIR)/solve_mode_dependend_2.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2 $(MY_PRECISION) $(OPTIM_FLAGS)   -o  $(ODIR)/solve_mode_dependend_2.o
$(ODIR)/solve_mode_dependend_3.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_3 $(MY_PRECISION) $(OPTIM_FLAGS)  -o  $(ODIR)/solve_mode_dependend_3.o
$(ODIR)/solve_mode_dependend_4.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4 $(MY_PRECISION) $(OPTIM_FLAGS)  -o  $(ODIR)/solve_mode_dependend_4.o
# debug
$(ODIR)/solve_mode_dependend_1.dbg.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1 $(MY_PRECISION) $(DEBUG_FLAGS)  -o  $(ODIR)/solve_mode_dependend_1.dbg.o
$(ODIR)/solve_mode_dependend_2.dbg.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2 $(MY_PRECISION) $(DEBUG_FLAGS)   -o  $(ODIR)/solve_mode_dependend_2.dbg.o
$(ODIR)/solve_mode_dependend_3.dbg.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_3 $(MY_PRECISION) $(DEBUG_FLAGS)  -o  $(ODIR)/solve_mode_dependend_3.dbg.o
$(ODIR)/solve_mode_dependend_4.dbg.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4 $(MY_PRECISION) $(DEBUG_FLAGS)  -o  $(ODIR)/solve_mode_dependend_4.dbg.o
#
# single prec (leave out my_precision)
$(ODIR)/solve_mode_dependend_1.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1  -o  $(ODIR)/solve_mode_dependend_1.sgl.o
$(ODIR)/solve_mode_dependend_2.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2  -o  $(ODIR)/solve_mode_dependend_2.sgl.o
$(ODIR)/solve_mode_dependend_3.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS)  -c solve_mode_dependend.c \
	-DCOMP_MODE_3  -o  $(ODIR)/solve_mode_dependend_3.sgl.o
$(ODIR)/solve_mode_dependend_4.sgl.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4  -o  $(ODIR)/solve_mode_dependend_4.sgl.o
#
# spherical versions of the blockinvert code
$(ODIR)/blockinvert.sph.o:	blockinvert.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL -c $< $(MY_PRECISION) \
	-o  $(ODIR)/blockinvert.sph.o
$(ODIR)/block_matrix.sph.o:	block_matrix.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_matrix.sph.o
$(ODIR)/read_stress_observations.sph.o:	read_stress_observations.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/read_stress_observations.sph.o
$(ODIR)/block_read_gps.sph.o:	block_read_gps.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_gps.sph.o
$(ODIR)/block_output.sph.o:	block_output.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_output.sph.o
$(ODIR)/block_solve.sph.o:	block_solve.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_solve.sph.o
$(ODIR)/block_read_bflt.sph.o:	block_read_bflt.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_bflt.sph.o
$(ODIR)/block_stress.sph.o:	block_stress.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_stress.sph.o
$(ODIR)/block_read_euler.sph.o:	block_read_euler.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS) -DBLOCK_SPHERICAL  -c $< $(MY_PRECISION) -o  \
	$(ODIR)/block_read_euler.sph.o
#
# single
$(ODIR)/blockinvert.sph.sgl.o:	blockinvert.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL -c $<  -o  $(ODIR)/blockinvert.sph.sgl.o
$(ODIR)/block_matrix.sph.sgl.o:	block_matrix.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_matrix.sph.sgl.o
$(ODIR)/read_stress_observations.sph.sgl.o:	read_stress_observations.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/read_stress_observations.sph.sgl.o
$(ODIR)/block_read_gps.sph.sgl.o:	block_read_gps.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_gps.sph.sgl.o
$(ODIR)/block_output.sph.sgl.o:	block_output.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_output.sph.sgl.o
$(ODIR)/block_solve.sph.sgl.o:	block_solve.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_solve.sph.sgl.o
$(ODIR)/block_read_bflt.sph.sgl.o:	block_read_bflt.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_bflt.sph.sgl.o
$(ODIR)/block_stress.sph.sgl.o:	block_stress.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_stress.sph.sgl.o
$(ODIR)/block_read_euler.sph.sgl.o:	block_read_euler.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -DBLOCK_SPHERICAL  -c $< -o  $(ODIR)/block_read_euler.sph.sgl.o

#
# the random noise version of coulomb stress
$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o: coulomb_stress.c $(GEN_P_INC) noise.dat
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -c coulomb_stress.c \
	-DADD_COULOMB_STRESS_NOISE=$(NOISELEVEL) \
	$(MY_PRECISION) -o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o

$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o: coulomb_stress.c $(GEN_P_INC) noise.dat
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -c coulomb_stress.c \
	-DADD_COULOMB_STRESS_NOISE=$(NOISELEVEL) \
	-o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o
# 

$(ODIR)/terminate.o:	terminate.c $(GEN_P_INC)
	$(MPICC) $(OPTIM_FLAGS)  $(CFLAGS)  $(MY_PRECISION) -c terminate.c -o $(ODIR)/terminate.o

$(ODIR)/test_solvers.o:	test_solvers.c $(GEN_P_INC)
	$(MPICC)  $(OPTIM_FLAGS) $(CFLAGS)  $(MY_PRECISION) -c test_solvers.c -o $(ODIR)/test_solvers.o

$(ODIR)/ex_dense.o:	ex_dense.c $(GEN_P_INC)
	$(MPICC)  $(OPTIM_FLAGS) $(CFLAGS)  $(MY_PRECISION) -c ex_dense.c -o $(ODIR)/ex_dense.o

$(ODIR)/ex_dense_v2.o:	ex_dense_v2.c $(GEN_P_INC)
	$(MPICC)  $(OPTIM_FLAGS) $(CFLAGS)  $(MY_PRECISION) -c ex_dense_v2.c -o $(ODIR)/ex_dense_v2.o



#
# some generic rules with normal dependencies
#
$(ODIR)/%.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(OPTIM_FLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) $(OPTIM_FLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) $(OPTIM_FLAGS)  $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f90 $(GEN_P_INC)
	$(F90)   $(F90FLAGS) $(OPTIM_FLAGS)  $(MY_PRECISION) -c $<  -o $(ODIR)/$*.o

# debug versions
$(ODIR)/%.dbg.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  $(DEBUG_FLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) $(DEBUG_FLAGS)  $(MY_PRECISION) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) $(DEBUG_FLAGS)  $(MY_PRECISION) -c $<  -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o:	%.f90 $(GEN_P_INC)
	$(F90)   $(F90FLAGS) $(DEBUG_FLAGS)  $(MY_PRECISION) -c $<  -o $(ODIR)/$*.dbg.o

# single prec versions
$(ODIR)/%.sgl.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS) $(OPTIM_FLAGS)  -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS)  $(OPTIM_FLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS)  $(OPTIM_FLAGS) -c $<  -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.f90 $(GEN_P_INC)
	$(F90)  $(F90FLAGS)  $(OPTIM_FLAGS) -c $<  -o $(ODIR)/$*.sgl.o
