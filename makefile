################################################################################
#
#
#  (C) Thorsten Becker, thwbecker@post.harvard.edu
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
# YOU WILL HAVE TO HAVE THE FOLLOWING ENVIRONMENT VARIABLES DEFINED:
#
# ARCH = architecture of the machine, e.g. ip27 or i686, that's what you get with
# uname -m | awk '{print(tolower($1))}'
#
# all machine dependent compiler flags are in makefile.$(ARCH)
#  where ARCH is 
#  ip27 or ip32 for SGI
#  i686 or something for LINUX 
#  sun4d or sun4u for SUN and so on
#
# you should edit the makefile.$(ARCH) files for machine specific settings such
# as compiler flags, memory limits and the like
#
# other flags (switches) should be set here, ie. added to the 
# MAIN_DEFINES line as in the examples below. 
#
# for calc_eigen_from_cart_stress you need EISPACK in the modified form as 
# supplied in the myeispack.tar.gz file. 
# Extract this files to the location as specified at the end
# of the makefile.$(ARCH) (by default $(HOME)/progs/src/eispack) 
# and run "make" in that directory to create the modified
# EISPACK libraries first
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
#  -DALLOW_NON_3DQUAD_GEOM      allows for different fault element types, so far rectangular and 
#                            point source are implemented. if all faults are rectangular, leave
#		  	     away for improved speed 
#
#
#  -DUSE_NUMREC_SVD          use SVD decomposition and backsubstitution from Numerical
#                            recipes. Default is from LAPACK since numerical reps
#                            didn't work for -O:fast=ip27 compilation on IRIX systems
#                            If this is not set, you will have to specify the flags
#                            for the inclusion of LAPACK in the makefile.$ARCH files.
#
#  -DUSE_GEOPROJECT          use geoproject and produce code that allows for I/O
#                            using geographic coordinates and projections. for this to
#                            work, you will have to have GMT (www.gmt.soest.hawaii.edu)
#                            installed - THIS USES GMT < VERSION 5
#  Example settings:
#
#  example with latency:
#              MAIN_DEFINES =   \
#                      -DBINARY_PATCH_EVENT_FILE -DUSE_DOUBLE_PRECISION  -DLATENCY=0.5  
#  example without:
#              DEFINE_FLAGS =  \
#                      -DBINARY_PATCH_EVENT_FILE -DUSE_DOUBLE_PRECISION 
#
COMMON_DEFINES =  -DBINARY_PATCH_EVENT_FILE -DCHECK_CI_ONE_WAY 
#
# noise level for random interact version. this version 
# doesn't get compiled automatically.
NOISELEVEL=1e-08
#
MAIN_DEFINES = $(COMMON_DEFINES)

#
#
# directory for object files
ODIR = objects/$(ARCH)/
# directory for binaries
BDIR = bin/$(ARCH)/
#
# choice of Okada routine, comment out if modified routine is used
OKROUTINE = $(ODIR)/dc3d.o
# 
# include the machine dependent flags
# 
include makefile.$(ARCH)
#
# add this for pgplot support, otherwise comment it out
# you will use runtime plotting capabilities
include makefile.pgplot
#
# add this for slatec NNLS routine support, otherwise comment it out
# you will use NNLS solving capabilities
include makefile.slatec
#
# add this for superlu support, otherwise comment it out
# you will loose sparse matrix SuperLU LU solver capabilities
#include makefile.superlu
#
# add this for geoprojection support
#
include makefile.geoproject

# add up all define flags
DEFINE_FLAGS = $(MAIN_DEFINES) $(SLATEC_DEFINES) \
	$(PGPLOT_DEFINES) $(SUPERLU_DEFINES)  \
	$(GEOPROJECT_DEFINES)
# defines and pgplot flags
FLAGS = $(DEFINE_FLAGS) $(PGPLOT_INCLUDES) $(SLATEC_INCLUDES) \
	$(SUPERLU_INCLUDES) $(GEOPROJECT_INCLUDES)
# C and FORTRAN compiler specific flags, add them to the 
# other flags
CFLAGS = $(FLAGS) $(SCARGS)
FFLAGS = $(FLAGS) $(SFARGS)




#
# list of object files
#
# solve mode dependend, ie. one for each mode of accesing the
# interaction matrix (improved speed, one would hope)
SMD_OBJS = $(ODIR)/solve_mode_dependend_1.o	\
	$(ODIR)/solve_mode_dependend_2.o	\
	$(ODIR)/solve_mode_dependend_3.o	\
	$(ODIR)/solve_mode_dependend_4.o 
#
# matrix solver 
#
MATRIX_SOLVER_OBJS = $(ODIR)/numrec_svd_routines.o $(ODIR)/nnls_lawson.o	\
	$(ODIR)/nnls.o $(ODIR)/svd.o $(ODIR)/solve.o		\
	$(ODIR)/lusolve.o $(ODIR)/sparse_solve.o \
	$(ODIR)/ilaenv_wrapper.o
#
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
	$(ODIR)/eval_okada.o $(ODIR)/tdd_coeff.o						\
	$(ODIR)/eval_green.o $(ODIR)/eval_triangle.o			\
	$(ODIR)/interact.o	$(ODIR)/mysincos.o 			\
	$(OKROUTINE) $(ODIR)/fracture_criterion.o			\
	$(ODIR)/myopen.o  $(ODIR)/randgen.o \
	$(ODIR)/string_compare.o 

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
	$(ODIR)/calc_stress.o		\
	$(ODIR)/read_geometry.o $(ODIR)/init.o $(ODIR)/help_and_comments.o	\
	$(ODIR)/output.o $(ODIR)/input.o $(ODIR)/matrixio.o			\
	$(ODIR)/calc_spatial_correlation.o	$(ODIR)/stress_aux.o		\
	$(ODIR)/quake.o $(MATRIX_SOLVER_OBJS) $(ODIR)/restart.o			\
	$(ODIR)/coulomb_stress.o $(POBJS) 
INPUT_OBJS_SGL = $(INPUT_OBJS:.o=.sgl.o)
#
# list of objects for the main program, interact
INTERACT_OBJS = $(ODIR)/rupture.o $(ODIR)/adjust_time_step.o 
# this is a random noise added version
INTERACT_NOISE_OBJS = $(ODIR)/rupture.o	$(ODIR)/adjust_time_step.o	\
	$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o
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
	$(ODIR)/read_stress_observations.o $(ODIR)/block_read_gps.o \
	$(ODIR)/eigensystem.o $(ODIR)/block_solve.o			\
	$(ODIR)/levmarq_numrec.o $(ODIR)/block_output.o  \
	$(ODIR)/block_stress.o $(ODIR)/block_levmarq.o 

BLOCKINVERT_SPH_OBJS = $(ODIR)/block_read_bflt.sph.o \
	$(ODIR)/coulomb_stress.o	 \
	$(ODIR)/block_read_euler.sph.o\
	$(ODIR)/block_matrix.sph.o	\
	$(ODIR)/block_read_gps.sph.o \
	$(ODIR)/read_stress_observations.sph.o	\
	$(ODIR)/block_solve.sph.o		\
	$(ODIR)/eigensystem.o 			\
	$(ODIR)/levmarq_numrec.o $(ODIR)/block_output.sph.o  \
	$(ODIR)/block_stress.sph.o $(ODIR)/block_levmarq.o 

FSTRESS2HOR_OBJS = $(ODIR)/block_read_gps.sph.o \
	$(ODIR)/block_output.o $(ODIR)/block_matrix.o\
	$(ODIR)/eigensystem.o $(ODIR)/block_read_bflt.o  \
	$(ODIR)/block_stress.o $(ODIR)/block_solve.o	

#
# real interact and test program
OBJ = $(ODIR)/main.o $(INTERACT_OBJS)
OBJ_SGL = $(OBJ:.o=.sgl.o)
NOBJ = $(ODIR)/main.o $(INTERACT_NOISE_OBJS)
TOBJ = $(ODIR)/test_stuff.o $(INTERACT_OBJS) 
#

#
# 
# include dependencies for all source codes
#
GEN_P_INC = interact.h precision_single.h precision_double.h \
	structures.h macros.h auto_proto.h auto_proto.sgl.h fortran_proto.h \
	filenames.h properties.h blockinvert.h

LIBLIST = $(ODIR)/libpatchio.a $(ODIR)/libinput.a 
LIBLIST_SGL = 	$(ODIR)/libpatchio.sgl.a $(ODIR)/libinput.sgl.a

#
# libraries, also linker flags
#
LIBS = $(MY_LIBDIR_SPEC)$(ODIR)/   -linput -lpatchio \
	$(COMPUTATIONAL_LIBS) $(DEBUG_LIBS) $(MATHLIB) 

LIBS_SGL = $(MY_LIBDIR_SPEC)$(ODIR)/   -linput.sgl -lpatchio.sgl \
	$(COMPUTATIONAL_LIBS) $(DEBUG_LIBS) $(MATHLIB) 


#
# list of all programs in groups
#

all: obj_directories libraries main_prog \
	tools converters geom_converters

really_all: obj_directories libraries main_prog \
	tools converters geom_converters \
	inoise analysis geographic_tools \
	 pgplot_progs 

main_prog: $(BDIR)/$(INTERACT_BINARY_NAME) $(BDIR)/$(INTERACT_BINARY_NAME).sgl

inoise: noisefile $(BDIR)/interact_noise.$(NOISELEVEL)

tools: misc_tools random_geom_tools random_prop_tools

misc_tools: $(BDIR)/makefault $(BDIR)/calc_interaction_matrix \
		$(BDIR)/project_stress $(BDIR)/calc_eigen_from_cart_stress \
	$(BDIR)/check_feedback $(BDIR)/fit_simple_stress_from_cart \
	$(BDIR)/calc_cart_from_eigen_stress \
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

converters: $(BDIR)/patch2xyz $(BDIR)/patch2geom  $(BDIR)/patch2vtk $(BDIR)/patch2bc \
	$(BDIR)/patch2corners $(BDIR)/patch2group \
	$(BDIR)/patch2xyzvec $(BDIR)/patch2poly3d $(BDIR)/patch2dis3d

geom_converters: $(BDIR)/points2patch $(BDIR)/tri2patch

test:    $(ODIR)/test_stuff 

matrix_test_progs: $(BDIR)/test_sparse $(BDIR)/test_optimize

pgplot_progs:  $(BDIR)/plotevents $(BDIR)/read_bin_events 

#
# some special targets
#

clean: 
	rm -rf $(ODIR)/*.o $(ODIR)/*.a  $(ODIR)/rii_files/ *.ps *.dat rii_files/ \
	auto_proto.h auto_proto.sgl.h

dist_clean:
	rm -rf $(BDIR)/*

obj_directories:
	if [ ! -s ./objects/ ]; then\
		mkdir objects;\
	fi;
	if [ ! -s ./objects/ ]; then\
		mkdir objects;\
	fi;
	if [ ! -s $(ODIR) ];then \
		mkdir $(ODIR);\
	fi;\
	if [ ! -s ./bin/ ];then\
		mkdir bin;\
	fi;\
	if [ ! -s ./bin ];then\
		mkdir bin;\
	fi;\
	if [ ! -s bin/$(ARCH)/ ];then \
		mkdir bin/$(ARCH);\
	fi;

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
	$(LD) $(LDFLAGS) $(OBJ) -o $(BDIR)/$(INTERACT_BINARY_NAME) \
		$(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 

$(BDIR)/$(INTERACT_BINARY_NAME).sgl: $(OBJ_SGL) $(GEN_P_INC) $(LIBLIST_SGL) 
	$(LD) $(LDFLAGS) $(OBJ_SGL) -o $(BDIR)/$(INTERACT_BINARY_NAME).sgl \
		$(LIBS_SGL) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 

$(BDIR)/interact_noise.$(NOISELEVEL): $(NOBJ) $(GEN_P_INC) $(LIBLIST) 
	$(LD) $(LDFLAGS) $(NOBJ) -o $(BDIR)/interact_noise.$(NOISELEVEL) \
		$(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 


$(ODIR)/test_stuff: $(TOBJ) $(GEN_P_INC)  $(LIBLIST)  
	$(LD) $(LDFLAGS) $(TOBJ) -o $(ODIR)/test_stuff \
	$(LIBS) $(PGLIBS)  $(SUPERLU_LIBS)  $(SLATEC_LIBS) 

$(BDIR)/randomflt: $(RANDOMFLT_OBJS)  $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(RANDOMFLT_OBJS) \
	-o $(BDIR)/randomflt $(LIBS)

$(BDIR)/generate_random_2d: $(GENERATE_RANDOM_2D_OBJS)  $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(GENERATE_RANDOM_2D_OBJS) \
	-o $(BDIR)/generate_random_2d $(LIBS)


$(BDIR)/patch2xyz: $(ODIR)/patch2xyz.o $(GEN_P_INC) $(ODIR)/read_geometry.o \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2xyz.o \
		-o $(BDIR)/patch2xyz  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2xyzvec: $(ODIR)/patch2xyzvec.o $(GEN_P_INC) \
	$(ODIR)/read_geometry.o $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2xyzvec.o \
		-o $(BDIR)/patch2xyzvec  $(LIBS) \
	$(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2poly3d: $(ODIR)/patch2poly3d.o $(GEN_P_INC) \
	$(ODIR)/read_geometry.o $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2poly3d.o \
		-o $(BDIR)/patch2poly3d  $(LIBS) \
	$(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2dis3d: $(ODIR)/patch2dis3d.o $(GEN_P_INC) \
	$(ODIR)/read_geometry.o $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/patch2dis3d.o \
		-o $(BDIR)/patch2dis3d  $(LIBS) \
	$(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/points2patch: $(ODIR)/points2patch.o  $(ODIR)/fit_plane.o $(ODIR)/libpatchio.a $(GEN_P_INC)  \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/points2patch.o $(ODIR)/libpatchio.a \
		$(ODIR)/fit_plane.o \
		-o $(BDIR)/points2patch  $(LIBS)  $(BLASLIB)

$(BDIR)/create_random_stress_file: $(ODIR)/create_random_stress_file.o \
	$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/create_random_stress_file.o \
		-o $(BDIR)/create_random_stress_file  \
	$(LIBS)  $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/calc_stress_stat: $(ODIR)/calc_stress_stat.o $(ODIR)/read_geometry.o \
	$(GEN_P_INC)  $(LIBLIST) $(ODIR)/calc_spatial_correlation.o
	$(LD) $(LDFLAGS)    $(ODIR)/read_geometry.o \
		$(ODIR)/libpatchio.a $(ODIR)/calc_spatial_correlation.o \
		$(ODIR)/calc_stress_stat.o -o $(BDIR)/calc_stress_stat  \
	$(LIBS)  $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/create_random_mu_file: $(ODIR)/create_random_mu_file.o \
	 $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/create_random_mu_file.o \
		 -o $(BDIR)/create_random_mu_file  \
	$(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/tri2patch: $(ODIR)/tri2patch.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/tri2patch.o \
		-o $(BDIR)/tri2patch  $(LIBS) 


$(BDIR)/patch2geom: $(ODIR)/read_geometry.o $(ODIR)/patch2geom.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o  $(ODIR)/patch2geom.o \
		-o $(BDIR)/patch2geom  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2vtk: $(ODIR)/read_geometry.o $(ODIR)/patch2vtk.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o  $(ODIR)/patch2vtk.o \
		-o $(BDIR)/patch2vtk  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2bc: $(ODIR)/patch2bc.o $(ODIR)/read_geometry.o $(GEN_P_INC)  \
	$(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/read_geometry.o \
		$(ODIR)/patch2bc.o -o $(BDIR)/patch2bc  $(LIBS) \
		 $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)	

$(BDIR)/patch2corners:  $(ODIR)/read_geometry.o $(ODIR)/patch2corners.o \
	$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/read_geometry.o $(ODIR)/patch2corners.o \
		-o $(BDIR)/patch2corners  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/patch2group:   $(ODIR)/patch2group.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)   $(ODIR)/patch2group.o \
		-o $(BDIR)/patch2group  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) 

$(BDIR)/randomize_strike: $(RANDOMIZE_STRIKE_OBJS)  $(GEN_P_INC)
	$(LD) $(LDFLAGS)  $(RANDOMIZE_STRIKE_OBJS) \
	-o $(BDIR)/randomize_strike $(LIBS)


$(BDIR)/makefault: $(ODIR)/makefault.o   $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/makefault.o  \
	-o $(BDIR)/makefault $(LIBS)

$(BDIR)/sort_events: $(ODIR)/sort_events.o $(ODIR)/myopen.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/sort_events.o  $(ODIR)/myopen.o \
	-o $(BDIR)/sort_events $(LIBS)

$(BDIR)/check_feedback: $(ODIR)/coulomb_stress.o $(ODIR)/check_feedback.o  $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/check_feedback.o \
	$(ODIR)/coulomb_stress.o \
		-o $(BDIR)/check_feedback  $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/mspectral: $(ODIR)/mspectral.o  interact.h $(ODIR)/myopen.o $(ODIR)/period.o
	$(LD) $(LDFLAGS)  $(ODIR)/mspectral.o $(ODIR)/myopen.o $(ODIR)/period.o \
		-o $(BDIR)/mspectral  $(LIBS)

$(BDIR)/calc_interaction_matrix: $(ODIR)/coulomb_stress.o \
	$(ODIR)/calc_interaction_matrix.o $(GEN_P_INC) \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/coulomb_stress.o \
	$(ODIR)/calc_interaction_matrix.o  \
	-o $(BDIR)/calc_interaction_matrix $(LIBS) $(SUPERLU_LIBS) \
		$(PGLIBS) $(SLATEC_LIBS) 


$(BDIR)/test_sparse: $(ODIR)/test_sparse.o $(ODIR)/coulomb_stress.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/test_sparse.o $(ODIR)/coulomb_stress.o \
	-o $(BDIR)/test_sparse $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS) $(LIBS)

$(BDIR)/test_optimize: $(ODIR)/test_optimize.o $(ODIR)/optimize.o  \
	$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS)  $(ODIR)/test_optimize.o \
	 $(ODIR)/optimize.o  \
	-o $(BDIR)/test_optimize $(LIBS) $(SUPERLU_LIBS) $(SLATEC_LIBS) $(PGLIBS)

$(BDIR)/project_stress: $(ODIR)/project_stress.o $(ODIR)/mysincos.o \
	$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/project_stress.o  $(ODIR)/llgeo.o \
	$(ODIR)/mysincos.o $(ODIR)/geometry.o \
	-o $(BDIR)/project_stress $(LIBS)

$(BDIR)/generate_slipdia: $(ODIR)/generate_slipdia.o $(ODIR)/mysincos.o \
	$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/generate_slipdia.o  $(ODIR)/llgeo.o \
	$(ODIR)/mysincos.o $(ODIR)/geometry.o \
	-o $(BDIR)/generate_slipdia $(LIBS)

$(BDIR)/calc_eigen_from_cart_stress: $(ODIR)/calc_eigen_from_cart_stress.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/calc_eigen_from_cart_stress.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/calc_eigen_from_cart_stress $(LIBS) $(EISPACK_LIB)

$(BDIR)/calc_cart_from_eigen_stress: $(ODIR)/calc_cart_from_eigen_stress.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/calc_cart_from_eigen_stress.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/calc_cart_from_eigen_stress $(LIBS) $(EISPACK_LIB)

$(BDIR)/fit_simple_stress_from_cart: $(ODIR)/fit_simple_stress_from_cart.o $(ODIR)/mysincos.o \
		$(ODIR)/llgeo.o $(ODIR)/geometry.o  $(ODIR)/eigensystem.o \
		$(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/fit_simple_stress_from_cart.o  $(ODIR)/llgeo.o \
		$(ODIR)/mysincos.o $(ODIR)/geometry.o $(ODIR)/eigensystem.o \
		-o $(BDIR)/fit_simple_stress_from_cart $(LIBS) $(EISPACK_LIB)

$(BDIR)/plotevents: $(ODIR)/plotevents.o $(ODIR)/plotting.o \
	$(ODIR)/plotting_palette.o $(ODIR)/mysincos.o $(ODIR)/llgeo.o \
	$(ODIR)/geometry.o $(GEN_P_INC)  $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/plotevents.o $(ODIR)/llgeo.o \
		$(ODIR)/geometry.o $(ODIR)/plotting_palette.o \
		$(ODIR)/myopen.o $(ODIR)/plotting.o \
		-o  $(BDIR)/plotevents  $(PGLIBS) $(LIBS) 

$(BDIR)/read_bin_events: $(ODIR)/read_bin_events.o $(INTERACT_OBJS) $(GEN_P_INC)  \
	 $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/read_bin_events.o  $(INTERACT_OBJS) \
		-o  $(BDIR)/read_bin_events $(LIBS) $(PGLIBS) \
		$(SUPERLU_LIBS) $(SLATEC_LIBS) 


$(BDIR)/blockinvert_sph: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST) \
		$(BLOCKINVERT_SPH_OBJS)  $(ODIR)/blockinvert.sph.o
	$(LD) $(LDFLAGS) $(BLOCKINVERT_SPH_OBJS) $(ODIR)/blockinvert.sph.o \
		$(MY_LIBDIR_SPEC)$(ODIR)/ $(MATHLIB)   $(GEOPROJECT_OBJS)			\
		-o  $(BDIR)/blockinvert_sph   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(DEBUG_LIBS) $(EISPACK_LIB) $(PGLIBS) 		$(COMPUTATIONAL_LIBS) 

$(BDIR)/fstress2hor: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o
	$(LD) $(LDFLAGS) $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2hor.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2hor   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB)

$(BDIR)/fit_mean_stress: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o
	$(LD) $(LDFLAGS) $(FSTRESS2HOR_OBJS) $(ODIR)/fit_mean_stress.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fit_mean_stress   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB)   $(DEBUG_LIBS) $(EISPACK_LIB)

$(BDIR)/fstress2eig: $(GEN_P_INC)  $(LIBLIST) $(GEOPROJECT_OBJS)\
		$(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o
	$(LD) $(LDFLAGS) $(FSTRESS2HOR_OBJS) $(ODIR)/fstress2eig.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/fstress2eig   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB)

$(BDIR)/block_evaluate_solution: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST)	\
		$(BLOCK_EVALUATE_SOLUTION_OBJS) 		\
		 $(ODIR)/block_evaluate_solution.o
	$(LD) $(LDFLAGS) $(BLOCK_EVALUATE_SOLUTION_OBJS)		\
		$(ODIR)/block_evaluate_solution.o 	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)	\
		-o  $(BDIR)/block_evaluate_solution   -linput -lpatchio	\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB)

$(BDIR)/block_checkflt: $(GEN_P_INC)  $(GEOPROJECT_OBJS) $(LIBLIST) \
		$(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o
	$(LD) $(LDFLAGS) $(BLOCKINVERT_SPH_OBJS) $(ODIR)/block_checkflt.o	\
			$(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS)			\
		-o  $(BDIR)/block_checkflt   -linput -lpatchio		\
		$(GEOPROJECT_LIBS)					\
		$(COMPUTATIONAL_LIBS) $(MATHLIB) $(DEBUG_LIBS)   $(EISPACK_LIB)

$(BDIR)/geo_okada: $(ODIR)/geo_okada.o $(ODIR)/coulomb_stress.o $(GEN_P_INC)  \
	$(GEOPROJECT_OBJS) $(LIBLIST) 
	$(LD) $(LDFLAGS) $(ODIR)/geo_okada.o $(MY_LIBDIR_SPEC)$(ODIR)/ $(GEOPROJECT_OBJS) \
		-o  $(BDIR)/geo_okada  $(ODIR)/coulomb_stress.o -lpatchio \
		$(GEOPROJECT_LIBS)	$(COMPUTATIONAL_LIBS) $(DEBUG_LIBS) $(MATHLIB)  

#
# C function prototyper
#

proto: 	auto_proto.h auto_proto.sgl.h

auto_proto.h: 
	rm -f auto_proto.h 2> /dev/null;\
	touch auto_proto.h;\
	cproto  $(DEFINE_FLAGS)  $(GEOPROJECT_INCLUDES) \
		$(PGPLOT_DEFINES) $(PGPLOT_INCLUDES) -DUSE_DOUBLE_PRECISION \
		$(SLATEC_INCLUDES)  $(SUPERLU_INCLUDES)  -f2 -q \
		`ls *.c | grep -v solve_mode_dependent.c | grep -v geoproject.c `  | \
		grep -v "void main("  | grep -v "int main(" > tmp.h;\
	mv tmp.h auto_proto.h
auto_proto.sgl.h: 
	rm -f auto_proto.sgl.h 2> /dev/null;\
	touch auto_proto.sgl.h;\
	cproto  $(DEFINE_FLAGS)  $(GEOPROJECT_INCLUDES) \
		$(PGPLOT_DEFINES) $(PGPLOT_INCLUDES)  \
		$(SLATEC_INCLUDES)  $(SUPERLU_INCLUDES)  -f2 -q \
		`ls *.c | grep -v solve_mode_dependent.c | grep -v geoproject.c `  | \
		grep -v "void main("  | grep -v "int main(" > tmp.h;\
	mv tmp.h auto_proto.sgl.h


#
# libraries
#
libraries: $(LIBLIST)

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
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F -DUSE_DOUBLE_PRECISION -o $(ODIR)/numrec_svd_routines.o
$(ODIR)/numrec_svd_routines.sgl.o: numrec_svd_routines.F $(GEN_P_INC)
	$(F77) -c  $(FFLAGS) numrec_svd_routines.F -o $(ODIR)/numrec_svd_routines.sgl.o
#
# other specialities
#
$(ODIR)/solve_mode_dependend_1.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_1 -DUSE_DOUBLE_PRECISION -o  $(ODIR)/solve_mode_dependend_1.o
$(ODIR)/solve_mode_dependend_2.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_2 -DUSE_DOUBLE_PRECISION -o  $(ODIR)/solve_mode_dependend_2.o
$(ODIR)/solve_mode_dependend_3.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_3 -DUSE_DOUBLE_PRECISION -o  $(ODIR)/solve_mode_dependend_3.o
$(ODIR)/solve_mode_dependend_4.o: solve_mode_dependend.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c solve_mode_dependend.c \
	-DCOMP_MODE_4 -DUSE_DOUBLE_PRECISION -o  $(ODIR)/solve_mode_dependend_4.o
#
# single prec
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
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/blockinvert.sph.o
$(ODIR)/block_matrix.sph.o:	block_matrix.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_matrix.sph.o
$(ODIR)/read_stress_observations.sph.o:	read_stress_observations.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/read_stress_observations.sph.o
$(ODIR)/block_read_gps.sph.o:	block_read_gps.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_read_gps.sph.o
$(ODIR)/block_output.sph.o:	block_output.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_output.sph.o
$(ODIR)/block_solve.sph.o:	block_solve.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_solve.sph.o
$(ODIR)/block_read_bflt.sph.o:	block_read_bflt.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_read_bflt.sph.o
$(ODIR)/block_stress.sph.o:	block_stress.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_stress.sph.o
$(ODIR)/block_read_euler.sph.o:	block_read_euler.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -DBLOCK_SPHERICAL  -c $< -DUSE_DOUBLE_PRECISION -o  $(ODIR)/block_read_euler.sph.o
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
	-DUSE_DOUBLE_PRECISION -o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).o

$(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o: coulomb_stress.c $(GEN_P_INC) noise.dat
	$(CC) $(CFLAGS) -c coulomb_stress.c \
	-DADD_COULOMB_STRESS_NOISE=$(NOISELEVEL) \
	-o  $(ODIR)/coulomb_noise_stress.$(NOISELEVEL).sgl.o

#
# some generic rules with normal dependencies
#
$(ODIR)/%.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS)  -DUSE_DOUBLE_PRECISION -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) -DUSE_DOUBLE_PRECISION -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) -DUSE_DOUBLE_PRECISION -c $<  -o $(ODIR)/$*.o
# single prec versions
$(ODIR)/%.sgl.o:	%.c $(GEN_P_INC)
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.f $(GEN_P_INC)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.sgl.o

$(ODIR)/%.sgl.o:	%.F $(GEN_P_INC)
	$(F77) $(FFLAGS) -c $<  -o $(ODIR)/$*.sgl.o
