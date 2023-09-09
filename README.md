#  interact: model fault interactions using dislocations in an elastic  medium
Copyright (C) Thorsten W. Becker 2000 - 2023                                

Interact uses dc3d.f as provided by Y. Okada as in Okada (BSSA, 1992)
and might come with linear algebra routines from SLATEC, LAPACK, and
EISPACK included. 
                                                                              
Petsc parallel solver interface built with much help from Dave May
(UCSD).
                                                                              
Might contain material copyrighted by others (e.g. Numerical Recipes
  etc).
                                                                              
                                                                              
See this readme, help.txt, and the makefile for documentation and
COPYRIGHT and COPYING for the license and warranty disclaimers.

If you use interact in real work, a citation to use would be:               

  @Article{becker02c,                                                         
    author = 	 {Becker, T. W. and Schott, B.},                               
    title = 	 {On boundary-element models of elastic fault interaction      
                  (abstract)},                                                 
    journal = 	 {EOS Trans. AGU},                                             
    year = 	 2002,                                                         
    volume =	 83,                                                           
    number =	 47,                                                           
    pages =	 {NG62A-0925}                                                  
  }                                                                           
                                                                              
thwbecker@post.harvard.edu                                                  
                                                                              

This directory contains source codes needed to compile the elastic 2-D
and half-space dislocation program "interact" as well as several tools
that go with it.

To compile, edit the makefile as described in the INSTALLATION
document and type 'make'. By default, we use the GNU compilers, with
options set in makefile.gcc

There are a number of optional additions
such as petsc, pgplot, GMT projection interfaces, and such. Those are
not important for the general operation, so if there are issues with
any of these, edit the makefile to comment out
include makefile.pgplot
include makefile.petsc
include makefile.geoproject

General settings as to the operation of interact are set in
"makefile", material properties such as the shear modulus are set in
"properties.h", and all system-dependent settings should be changed in
"makefile.gcc" or similar.  The reasoning behind this setup is that
you are not likely to change the, say, shear modulus very often
because it simply affects the scaling. Likewise, you will also not
want to change the way interact works that frequently, so we hardwire
several options to avoid calls to IF statements for efficiency. After
compilation, the operation of interact is affected by command line
arguments and input files.

For documentation, use 'program -h' after compilation, and refer to
the numerous comments in the source code as well as in the
makefiles. The output of 'interact -h' is supplied as the file
'help.txt'.

The example.?.txt files describe some simple example calculations to,
hopefully, elucidate the usage of interact and affiliated programs.
(See our Fall 2002 AGU poster at
http:///www-udc.ig.utexas.edu/external/becker/agu02_disp.pdf).

Typing 'make' will compile most programs, see the makefile. If you
want all tools, type 'make really_all'. 


A list of the programs that will be compiled and several batch scripts
that are provided in addition follows with short descriptions what
they do. For each program, more information can be obtained by typing
'program -h'.

IN PARTICULAR, THE INPUT FORMATS ETC. ARE DESCRIBED WHEN TYPING
interact -h

1) main program:

interact: 
	reads in fault geometry and boundary conditions. Can do a
	one-step calculation where a system of equations is solved
	based on boundary conditions in stress and/or displacement on
	the rectangular fault patches. BCs can have non-negativity
	constraints, e.g. motion only in one direction. Program can
	also simulate a loading experiment where patches behave
	according to a friction criterion (so far, only Coulomb with
	static and kinetic friction is implemented).

2) generation of input/output files for interact:

makefault:
	produces regular fault segments based on input
	parameters such as strike that can be selected on the command
	line with arguments like '-strike 32'. output is in the
	geom.in format.

points2patch:
	convert four points defining a rectangle in 3-D to a best
	fitting, Okada type fault patch

patch2geom:
	converts a patch file (with fault segments) into a file that
	holds the patch geometry as a geomview GEOM file

patch2xyz:
	converts a patch file (with fault segments) into a file that
	holds the geometry in a xyz file to be read with GMT (<5)

create_random_stress_file:
	create fsi.in input file for interact which holds randomly
	distributed initial values for the strike, dip and normal
	stress components for each patch of the fault system.

create_random_mu_file:
	create fp.in input file for interact which holds randomly
	distributed initial values for static and dynamic friction
	coefficients each patch of the fault system.

randomflt:
	produces randomly oriented and sized fault patches in patch format

randomize_strike:
	reads in locations for future faults and randomizes the
	strike, optionally the dip and depth alignment to. output in
	patch format.

3) main plotting

pgeom:  (GMT batch script)
	plots the GMT xyz file produced by patch2xyz in 3D using a GMT batch
	command file

pgeom2d: (GMT batch script)
	 plots 2D projection of the GMT xyz file produced by patch2xyz

4) main other tools

project_stress:
	projects a stress tensor (six components) into a fault local
	system where the traction vector is resolved onto the local
	strike, dip, and normal vector

calc_eigen_from_cart_stress:
	calculates the eigensystem of a given cartesian stress matrix

calc_cart_from_eigen_stress:
	convert a given stress state in principal axis to a cartesian
	stress matrix


5) other programs

i) plotting or data analysis

pdisp:  (batch script)
	plots displacement fields as produced by one-step calculations
	with interact

peventfile: (batch/gnuplot script)
	plots the cevents.dat file that holds individual rupture
	events produced by a loading simulation with interact

pfstress:
	plots the stress and slip on a fault as a function of time,
	read flkt.*.dat files

plotevents: (only if PGPLOT support included)
	reads a binary event file produced by interact and plots the
	activations using PGPLOT

plotgr: (batch script)
	reads the cevents.dat file that holds individual rupture
	events produced by a loading simulation with interact and
	plots a Gutenberg Richter statistic

pslip:  (GMT script)
	plots slip on the patches that belong to a certain fault
	group. useful for one-step calculation when solving for slip
	subject to stress boundary conditions.

ii) input/output conversion

patch2bc:
	converts a patch file (with fault segments) into a file that
	hold boundary condition input based on the orientation of
	faults.


patch2corners:
	converts a patch file (with fault segments) into a file that
	holds the four corners of the fault patch as (x,y,z) triples


patch2group:
	converts a patch file (with fault segments) into a file that
	holds the group geometry (can be made up of several patches)
	in patch format


patch2xyzvec:
	converts a patch file (with fault segments) into a file that
	holds the corner coordinates and basis (strike, normal, dip) vectors

patch2poly3d:
	converts a patch and a boundary condition file to Poly3D
	format

read_bin_events:
	reads in binary events from an interact loading simulation and
	writes the ASCII equivalent to stdout.

sort_events:
	reads in cevents.dat file and sorts for fore- and 
	aftershocks

tri2patch: (not really implemented yet)
	convert three points defining a triangle in 3-D to internal format.

iii) interaction matrix related

calc_interaction_matrix:
	calculates the interaction matrix which holds the Green's
	functions coefficients for stress change of type l on fault i
	due to slip of fault j of type k. input is a file in the
	'patch' format. writes the result to file, and prints to the
        screen if matrix is small (nrflts<100)


check_feedback:
	calculates the interaction matrix for a fault system as
	described by a file in the 'patch' format and checks for
	interactions between patches that would lead to a positive
	feedback loop. this can arise when the coulomb stress changes
	on another patch are larger than the reduction of coulomb
	stress on the patch itself. this is an effect of the
	discretization, on average stresses of groups of patches will
	always be OK.


________________________________________________________________________________

Output from interact -h follows (run code for updated versions!)

________________________________________________________________________________
main: binary: interact
main: compiled for double precision on Apr  9 2023 at 19:29:33, initializing
init: version: $Id: init.c,v 2.51 2011/01/09 02:02:43 becker Exp $
init: compiled on Apr  9 2023 19:29:28

Interact: calculate fault stresses and displacements in a half space or in 2-D
          using a boundary element approach.

Program reads in fault patch geometries (usually rectangular) and either solves
for stresses and/or displacements in a one-step problem for various boundary
conditions (sign-constrained or unconstrained), or for a simulated loading
cycle where each patch follows a frictional constitutive law (e.g. Navier-Coulomb)
and repetitive rupture is allowed during continuous "plate-tectonic" loading.

WHEN STRESSES (STRESS TENSOR COMPONENTS) OR PRESSURE ARE REFERRED TO,
POSITIVE VALUES MEAN EXTENSIONAL AND COMPRESSIONAL REGIMES,
RESPECTIVELY (PHYSICS CONVENTION).

Output is written to files and in real-time to X11 if PGPLOT support was compiled in.
(PGPlot support was compiled in.)

(1) The fault geometry is input via the file "geom.in",
    a list of fault patches.
    This file has the following ("patch") format for rectangular faults or point sources
    (free format ASCII list of parameters):

    x_0 y_0 z_0 strike_0 dip_0 hlength_0 hwidth_0 group_0 ...
    x_1 y_1 z_1 strike_1 dip_1 hlength_1 hwidth_1 group_1 ...
    ...
    x_N-1 y_N-1 z_N-1 strike_N-1 dip_N-1 hlength_N-1 hwidth_N-1 group_N-1 ...

    Here, N is the total number of patches and

    - x,y,z are the coordinates of the center of the rectangular fault patch or
      in the case of a point source, the point source location.

    - strike and dip are fault plane orientation angles in degree
      strike is defined in degrees clockwise from North
      dip    is defined in degrees downward from the horizontal, 90 degrees means vertical fault
      note that right now, the geometrical rake has to be 90 or 0 degrees, i.e. no tilted rectangles
      slip can, however, have a rake as given by different along-strike and along-dip values, see below

    - hlength (L) and hwidth (W) are the patch's -->HALF<-- length and half width in 
      strike and dip direction, respectively. The patch area is thus 4LW.

    - group is the patch's fault group, a number that is used to define faults that consist of
      several patches (discretizing fault planes by rectangular elements) (0 <= group <= N-1))
    - ... at the end mean possible additional input, see below.
    If you would like to use the 2-D mode, point sources or triangles, recompile with ALLOW_NON_3DQUAD_GEOM flag set.


    Examples:

    0 0 -1 30 90 4 2 0

    for a rectangular, vertical fault patch centered at depth -1, striking 30 degrees from North
    with length 8 and width 4, patch group number is 0.

    Please note that there are various programs to convert to and from the above patch format,
    e.g. patch2xyz to go to GMT style xyz coordinates for the edges of each patch,
    patch2geom to go from patch format to the Geometry center's geomview OFF format, or
    points2patch to convert a set of four points in FE ordering to the above patch format.


(2) Boundary conditions and the operational mode are input via the "bc.in" file
    which has the following format (free format ASCII list of parameters):

    operational_mode:
      1: one step calculation
      2: loading simulation, and
      3: loading simulation and x window plotting.


 IF ONE STEP CALCULATION IS CHOSEN, the following input fields should be

     print_bulk_fields

      If print_bulk_fields is set to 0, no stress or displacement fields will be output.
      If print_bulk_fields is set to 1, the program will read in the boundaries from the next line

       xmin xmax n ymin ymax m zmin zmax o

       where n, m, and o are the number of samples between the given limits.

         If o is set to -1, the output will be in a plane with the
          average strike and dip direction of fault group 0. x/y min/max are then the
          limits in strike and dip direction.
         If o is set to -2, the output will be in a plane with the
          average strike and normal direction of fault group 0. x/y min/max are then the
          limits in strike and normal direction.
         In these cases, the coordinates in the output files will correspond to the local,
         rotated systems, e.g. for -2, x will go along strike, and y along the normal direction
         with respect to the average fault patch geometry, with origin in the center of the fault
         group. HOWEVER, the displacements and stresses will still be given in the old system!
         If the -pc switch is set, the global coordinates will be printed to "plane.xyz"

      If print_bulk_fields is set to 2, the program will read the output locations from file "oloc.dat"
      This file has the x, y, and z coordinates in ASCII as an unformatted list.


     patch_number boundary_code boundary_value

      patch_number runs from 0 to N-1, where N is determined from the number of patches
      (fully read lines, unformatted) in the geometry file.

      If patch_number < 0, then the line should instead read

       -increment start_fault stop_fault boundary_code boundary_value ... ...

      and boundary code boundary_code will be assigned with value
      boundary_value from start_fault to stop_fault with increments increment.

      If start_fault < 0 and stop_fault < 0, will select all patches.

     The patch_number ... boundary_value line can be repeated as often as necessary (say, twice for each fault
     if the strike and dip movement modes are to be activated the same time).


   Possible boundary conditions (as indicated by their integer codes) are:

   (A) 0, 1, 2: slip on fault specified (if called several times, will add up)
    0 means strike, 1 dip, and 2 normal direction

    For strike: positive values of slip mean left-lateral fault, negative right-lateral;
    	this corresponds to a resulting drop resp. increase in the strike component of stress.
    	(We use stress values with signs, that is a negative stress drop corresponds to a
    	reduction of positive shear stress, a positive stress drop to reduction of shear
    	stress with negative sign. Both times, the absolute value of the shear stress is
    	reduced).
    For dip:    positive values of slip mean thrust fault, negative normal fault;
    	this corresponds to a resulting drop resp. increase in the dip component of stress.
    For normal: positive values of slip mean explosive source, negative implosive;
    	this corresponds to a resulting drop resp. increase in the normal component of stress.


   (B) 10, 11, 12 or 20, 21, 22 or 30, 31, 32: stresses are specified

    where the following codes indicate the degree of freedom for the faults to
     achieve that stress starting from 0
    10	strike slip, either way
    11	strike slip, left lateral   (stress drop bc usually < 0 for activation)
    12	strike slip, right lateral  (stress drop bc usually > 0 for activation)
    20	dip direction, either way
    21	dip direction, thrust       (stress drop bc usually < 0 for activation)
    22	dip direction, normal       (stress drop bc usually > 0 for activation)
    30	normal direction, either way
    31	normal direction, explosive (stress drop bc usually < 0 for activation)
    32	normal direction, implosive (stress drop bc usually > 0 for activation)

    If any of the strike or dip modes (10, 11, 12 or 20, 21, or 22) are given as 
    (110, 111, 112 or 120, 121, or 122) instead (add 100) this means
    that the target stress should be corrected by the change in normal
    stress times the dynamic coefficient of friction (makes sense for Coulomb type rupture).

   (C) 50, 51, 52: read in background loading (see -l, NOT -s), reduce resolved stresses on faults

    50	such that faults are shear-stress free after slip in strike and dip directions
    51	or such that the shear stress is given by the static friction coefficient
      	times the normal stress, after slip in strike and dip direction
    52	or such that faults are shear stress free after slip in strike direction only
    For the above two cases, the BC value (the last number) is meaningless.


    NOTE: The final stresses, e.g. in flt.dat, after slip will correspond to the stress drops
    needed to reduce the background load as specified above. The stresses will NOT include
    the resolved background stresses since those are only used to initialize the problem
    as boundary conditions. I.e. add the resolved background shear stress to the shear stress
    specified in flt.dat stress to get shear stress free conditions for the 50 case, for example.
    For debugging, the resolved shear, dip, and normal stresses for each patch will be written to "rstress.dat".

    Mixed boundary conditions: If you want to mix specified slip boundary conditions with stress
    boundary conditions, slip has to be completely specified first such that subsequent stress conditions
    can take the stress changes induced by the slip into account.


    Example:
    1 1
    -5 5 51 -5 5 51 -1 -1 1
    -1 -1 -1 10 1


 IF LOADING SIMULATION IS CHOSEN (with or without X output)

   the first three numbers are:

    time_step print_interval stop_time

   where time runs from 0 to stop_time in base time steps of time_step
   and the program will print to fault files every print_interval intervals.

   The time_step will normally be adjusted (but see -ct) based on the
   program's prediction of the next rupture time or change in shear stress sign.
   The max step is given by print_interval to ensure frequent
   enough output to the flt.*.dat files.
   The minimum time step was set to 1e-15.
   Time has to monotonously increase from zero.

   If X windows output is chosen, the next line (two numbers) should read

    x_plot_interval x_scroll_interval

   where x_plot_interval and x_scroll_intervals give the timing for the map view
   and time vs. x plots, respectively.

   All following lines assign active modes to faults and are of format:

    patch_number activation_mode_code

   where fault numbers should be given from 0 to N-1, N being the number of patches.
   If patch_number < 0 then give

    -increment start_fault stop_fault activation_mode_code

   instead. If start_fault<0 and stop_fault<0 will select all patches, as above.

   Activation_mode_codes are:
      0		inactive(default)
     10		strike slip,   either way
     110	strike slip,   Coulomb correction, either way
     11		strike slip,   left lateral
     111	strike slip,   Coulomb correction, left lateral
     12		strike slip,   right lateral
     112	strike slip,   Coulomb correction, right lateral
     20		dip direction, either way
     120	dip direction, Coulomb correction, either way
     21		dip direction, thrust
     121	dip direction, Coulomb correction, thrust
     22		dip direction, normal
     122	dip direction, Coulomb correction, normal
     30		normal direction, either way (not implemented yet)
     31		normal direction, explosive  (not implemented yet)
     32		normal direction, implosive  (not implemented yet)
     40		maximum stress direction in strike/dip plane
     140	maximum stress direction in strike/dip plane, Coulomb correction

   If Coulomb correction is set, the modes that involve slip in strike and dip direction will
   result in a correction for normal stress changes, assuming that the target
   stress drop is modified by the product of normal stress change during slip and
   the dynamic friction coefficient.


 NOTICE:
 Program maintains normal (opening) modes. If this is not needed, recompile
 with the NO_OPENING_MODES switch set to save memory.


PARAMETERS:

 Elastic and frictional law properties are hard-coded and should be set in the "properties.h" file.
 Whenever you want to change those, you will have to recompile the program. This was done to improve
 the execution speed, assuming that material properties will not change between different model runs.

 In the following, we will refer to two different matrices, A and I. The interaction, or I, matrix is
 assembled once for the loading simulations and holds all possible stress influence coefficients for all
 combinations of patch slips for a given geometry. I is typically held in memory and can be quite big, depending
 on the problem. Hence, there is the option to store a reduced, sparse version of I. If patches are activated
 (slipping) during a loading simulation (according to their frictional law), a system of equations
 A x = b is set up to solve for the slip vector x that corresponds to the target stress drop. This A matrix
 is also the only matrix that is assembled for a one-step calculation as we specifically know the patch/boundary
 conditions we have to consider in this case. There are several possible modes for solving completely
 unconstrained A x = b systems; constrained systems (where some entries in x have to be >=0 for unidirectionally
 constrained slip) are always solved with a direct implementation of Lawson and Hanson's solver.

 The following parameters can be set to non-default values on the command line:

 -c  value      sets the cohesion term of the friction law to value, by default: 5
                If the cohesion=0 at all times, recompile with NO_COHESION set
                to improve speed. (Stress drop is set to 1, which is 0.0001 times mu).

 -pr value      set the background pressure to value "value".
                The background pressure amplitude is 10 by default, both for one-step and loading
                simulations. The pressure is constant with depth (see properties.h).

 -ms value      set minimum stress drop to value "value". (Characteristic stress drops are
                set to 1.) The minimum stress drop is zero by default, ie.
                allowing stress drops of all sizes (loading simulation).

 -ei value      set matrix cutoff value for sparse storage of I or A in absolute values
                             (5.000000e-05 by default.)
                (only used if I matrix is stored as sparse, see -nd option, or for sparse Ax=b solver.)

 -wc value      will use value as the maximum fraction of singular value for an SVD solution (default: 7.00e-15)



OPTIONS:

 The following operational modes (switches) can be switched on or off by
 using the toggle options listed below (type, eg, interact -s)
 Please read carefully which switches are ON by default. If those
 active switches are selected, they will be switched OFF instead.

 -s  read in initial stress values for each fault from file "fsi.in", format:

     s_strike^0 s_dip^0 s_normal^0
     s_strike^1 s_dip^1 s_normal^1
     ...

     (repeated for each patch). Here, the s_i are the initial values of the stress 
     tensor resolved on the fault slip directions. All values are ADDED to the
     homogeneous background stress as set in the program, ie. use all zeros for normal
     initial stress state. (loading simulation)
 --> WARNING: This option will only work for the loading simulation mode, for one-step calculations
     specify fault stresses as individual boundary conditions, or use -l to get a resolved background
     stress.

     ON by default, if switch is set will be OFF.

 -f  read in fault properties from file "fp.in", so far only parameters for the friction law
     are supported. Format:

     mu_static^0 mu_dynamic^0
     mu_static^1 mu_dynamic^1
     ...

     (repeated for each fault). Static and dynamic coefficients of friction for
     Coulomb.
     ON by default, if switch is set will be OFF.

 -i  suppress all interactions between fault patches except the self-effect of slip on any
     patches within the group that the patch belongs to and the effect of slip
     on the number 0 master group (this should be the main fault).
     (Meant as a experimental feature for LOADING SIMULATIONS ONLY!)
     OFF by default, if switch is set will be ON. See -ni for suppressed interactions in one step.

 -ni suppress all interactions between faults except the interactions of fault patches

 -ct use constant time steps for the loading simulation
     OFF by default, if switch is set will be ON.

 -b  suppress output of all fault stress and slip data. Normally,
     individual output files for each
     group of faults are generated when the number of groups is smaller than 50.
     If -b is set, this maximum number will be set to zero.

 -p  prints information about the cumulative slip of every patch in a group to the
     slipline.*.dat files if the fault files are used (see -b). (loading simulation)
     OFF by default, if switch is set will be ON.

 -w  whole fault activation mode (loading simulation)
     If set, one patch will activate all patches in the group for sliding.
     Also, if one patch is switched off because of low normal stress, all
     patches in the group will be switched off as well (*unless* -v is set).
     OFF by default, if switch is set will be ON.

 -k  keep slipping mode (loading simulation)
     If set, the first critical patches can trigger slip on other patches.
     If so, the initial patches will keep slipping, and the solution is obtained
     at each iteration until no more triggering takes place. If not set, slip
     on triggering patches is exhaustive, and iteration checks for triggering step
     by step (faster, since old solution doesn't have to be removed from stress
     field). The exhaustive method is affected by tiny differences in critical stress!
     If keep slipping, critical stress is -7.000000e-15, else -7.000000e-15.
     ON by default, if switch is set will be OFF.

 -v  whole fault deactivation mode (loading simulation)
     If set (and -w is *not* set), the deactivation (low compressive normal stress) of a patch
     will lead to a deactivation of the whole fault group.
     If -w is on (whole fault activation mode), this is the default. If you use -v in this
     case, the whole fault de-activations will be suppressed.
     OFF by default, if switch is set will be ON.

 -l  reads a/b factors for the linear background loading stress relationship.
     If set, the file "smlin.in" is read, it should hold 6 pairs of a_i b_i factors
     where the *non-hydrostatic* background stress at time t is given by s_i=a_i+b_i*t
     s_i = {s_xx, s_xy, s_xz, s_yy, s_yz, s_zz}, so that the format is

     a_xx b_xx
     a_xy b_xy
     ...
     a_zz b_zz

     To these factors, a depth independent hydrostatic pressure of 10
     will be added. You can change the amplitude of the pressure (e.g. to zero) with -pr.

     If no factors are read in (if "smlin.in" is not found), the code will use simple shear loading 
     with b_xy=1 plus the depth independent hydrostatic pressure of 10.

     (Modify pressure amp. with -pr). If you want to specify resolved background
      stresses in a one-step calculation, simply use the a_i factors to specify the stress 
      at time t=0 but mind the hydrostatic part.

     ON by default, if switch is set will be OFF.

 -nd switches to sparse matrix I storage, no matter what the I matrix size is.
     By default, the program goes from dense to sparse only if the
     matrix is bigger than 8000 MB.
     The sparse matrix storage implies that matrix elements whose absolute
     value is less than 5.000e-05
     are discarded. This behavior can be changed by specifying
     a different cutoff value using the -ei parameter option.
     (Make sure to use -si if you want to keep the matrix.)
     OFF by default, if switch is set will be ON.
     note that sparse here refers to only the storage mode of the I matrix,
     not the solution method of A.x=b. (see -ss).

 -ci checks the interaction matrix for positive Coulomb stress feedback loops, once calculated
     This may be useful for loading simulations.
     OFF by default, if switch is set will be ON.

 -si saves the interaction (I) matrix during a loading simulation,
     if I is written to a file during runtime and not held in memory.
     Else, "/tmp/i" and "/tmp/i.hdr" will be deleted.
     OFF by default, if switch is set will be ON.

 -oi reads old interaction (I) matrix from files "/tmp/i.dat" and "/tmp/i.hdr".
     Saves time when the geometry
     does not change between subsequent model runs (this is not automatically checked!)
     (Make sure to use -si if you want to keep the matrix.)
     OFF by default, if switch is set will be ON.

 -sa saves the A x = b A matrix to "a.dat" and "a.dat.hdr" during a one-step calculation.
     Binary format, with dimensions (n x m) and precision (p, in bytes) in the header file as n p m
     OFF by default, if switch is set will be ON.

 -oa reads A matrix from files "a.dat" and "a.dat.hdr" during a one-step calculation.
     Saves time when A
     does not change between subsequent model runs (this is not automatically checked!)
     (Make sure to use -sa if you want to keep the matrix.)
     NOTE that A is only constant for a constant geometry and the boundary conditions are
     either pure slip or pure stress components.
     If you are specifying a Coulomb stress change, elements of A will be modified
     by the boundary conditions, and not only b as one might think!
     OFF by default, if switch is set will be ON.

 -r  attempts to restart a loading simulation model run based on previous results
     WARNING: all result files will be overwritten NOT appended, so save before.
     For this to work, "geom.in" and "bc.in" (and possibly "fsi.in" and "fp.in"
     should not be substantially changed.

     Slip events from the previous run in the format as in the output file "events.bin" will be read in
     from the restart event file "restart.events.bin"

 -gv will output files with slip on faults in the Geometry center's geomview COFF format,
     with a CQUAD file for each fault group. Slip values are scaled to the global maxima.
     Note that the fltdat2cquad script can produce CQAUD files from geom.in and flt.dat files.
     OFF by default, if switch is set will be ON.

 -sn will print NaN values in displacement and stress outputs of the one-step calculation.
     Normally, those locations (end patch/segments) yielding one or more inf values in the
     stress matrix or the displacement vector will be omitted during output.
     OFF by default, if switch is set will be ON.

 -pc will print the global coordinates if fault-local planes are chosen for bulk stress
     and displacement output (-1 or -2 for o in zmin zmax o), writes xyz to "plane.xyz"
     ON by default, if switch is set will be OFF.

 -sv uses SVD solver for the unconstrained A x = b system.
     This approach should reliably find the least-squares minimum norm solution
     if the A matrix is singular.
     Note that LU can be orders of magnitude faster than SVD but use caution.
     If non-negative solutions are sought (as for specified slip directions),
     the program will automatically default back to the NNLS solver.

 -sl uses LU solver for the unconstrained A x = b system (default).
     LU can be orders of magnitude faster than SVD, but use caution.
     If non-negative solutions are sought (as for specified slip directions),
     the program will automatically default back to the NNLS solver.

 -d  run in debug mode
     OFF by default, if switch is set will be ON.


 -h  prints out this help message and exits to the operating system

(C) Thorsten Becker, thwbecker@post.harvard.edu, 1999 - 2021
    interact - boundary element code for elastic half-spaces
    Main 3-D dislocation code based on dc3d.f by Y. Okada, as of Okada, BSSA, 1992
    2-D segment slip solution from Crouch and Starfield (1973)
    May include routines based on copyrighted software of others.
    Distributed under the GNU public license, see "COPYING".

