#include "interact.h"
#include "properties.h"
/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.havard.edu

*/
/*
  print some documentation
*/

#define PE(x) {fprintf(stderr,"%s\n",x);}
void phelp(void)
{
  PE("");
  PE("Interact: calculate fault stresses and displacements in a half space or in 2-D");
  PE("          using a boundary element approach.");
  PE("");
  PE("Program reads in fault patch geometries (usually rectangular) and either solves");
  PE("for stresses and/or displacements in a one-step problem for various boundary");
  PE("conditions (sign-constrained or unconstrained), or for a simulated loading");
  PE("cycle where each patch follows a frictional constitutive law (e.g. Navier-Coulomb)");
  PE("and repetitive rupture is allowed during continuous \"plate-tectonic\" loading.");
  PE("");
  PE("When stresses (stress tensor components) or pressure are referred to, positive values mean");
  PE("extensional and compressional regimes, respectively (physics convention).");
  PE("");
  PE("Output is written to files and in real-time to X11 if PGPLOT support was compiled in.");
#ifdef USE_PGPLOT
  PE("PGPlot support was compiled in. ");
#else
  PE("(No PGPlot support, X11 output unavailable. Compile with USE_PGPLOT flag set, if wanted.)");
#endif
#ifdef USE_PETSC
  PE("");
  PE("PetSc support was compiled in, providing limited access to parallel solves for one-step problems.");
  PE("      For LU solve, use    \"-pc_factor_mat_solver_type scalapack -mat_type scalapack\" or");
  PE("                           \"-pc_factor_mat_solver_type elemental -mat_type elemental\".");
  PE("      For iterative solve \"-ksp_type fgmres -pc_type none   -ksp_max_it 10000 -ksp_rtol 1.0e-8\" or");
  PE("                          \"-ksp_type fgmres -pc_type jacobi -ksp_max_it 10000 -ksp_rtol 1.0e-8\".");
  PE("      Check the makefile for other solver options and MPI settings.");
  PE("      When running a one-step computation, will also compute stress fields in parallel.");

#else
  PE("(Petsc support not compiled in, if parallel support is desired, check makefile.)");
#endif
  PE("");
  fprintf(stderr,"(1) The fault geometry is input via the file \"%s\",\n    a list of fault patches.\n",
	  GEOMETRY_FILE);
  PE("    This file has the following (\"patch\") format for rectangular faults or point sources");
  PE("    (free format ASCII list of parameters):\n");
  PE("    x_0 y_0 z_0 strike_0 dip_0 hlength_0 hwidth_0 group_0 ...");
  PE("    x_1 y_1 z_1 strike_1 dip_1 hlength_1 hwidth_1 group_1 ...");
  PE("    ...");
  PE("    x_N-1 y_N-1 z_N-1 strike_N-1 dip_N-1 hlength_N-1 hwidth_N-1 group_N-1 ...\n");
  PE("    Here, N is the total number of patches and\n");
  PE("    - x,y,z are the coordinates of the center of the rectangular fault patch or");
  PE("      in the case of a point source, the point source location.");
  PE("");
  PE("    - strike and dip are fault plane orientation angles in degree");
  PE("      strike is defined in degrees clockwise from North");
  PE("      dip    is defined in degrees downward from the horizontal, 90 degrees means vertical fault");
  PE("      note that right now, the geometrical rake has to be 90 or 0 degrees, i.e. no tilted rectangles");
  PE("      slip can, however, have a rake as given by different along-strike and along-dip values, see below");
  PE("")
  PE("    - hlength (L) and hwidth (W) are the patch's -->HALF<-- length and half width in ");
  PE("      strike and dip direction, respectively. The patch area is thus 4LW.");
  PE("");
  PE("    - group is the patch's fault group, a number that is used to define faults that consist of");
  PE("      several patches (discretizing fault planes by rectangular elements) (0 <= group <= N-1))");
  PE("    - ... at the end mean possible additional input, see below.");
#ifdef ALLOW_NON_3DQUAD_GEOM
  PE("");
  PE("    Since the \"ALLOW_NON_3DQUAD_GEOM\" flag was used during runtime, the program can deal with");
  PE("    some additional patch geometries besides rectangular. Those are selected as follows:");
  PE("");
  PE("    - point source in half-space");
  PE("      If ONLY fault half-length is set to a negative value, a point source will be assumed instead");
  PE("      of a rectangular fault patch. In this case, width should be set to the `fault' area");
  PE("      and half-length is equivalent to -aspect_ratio, where aspect_ratio is some equivalent L/W.");
  PE("");
  PE("      WARNING: Not properly tested yet.");
  PE("");
  PE("    - triangular in half-space:");
  PE("      If BOTH fault half-width and length are negative, then the patch is a triangular element.");
  PE("      In this case, x, y, z, dip, strike, length, and width have no meaning but will be reassigned.");
  PE("      The next nine numbers in the input geometry line will be the coordinates of the three nodes of the");
  PE("      triangle, however. FE-style counterclockwise numbering, input is then in the format:\n");
  PE("      999 999 999 999 999 -1 -1 group_0 x_x^1 x_y^1 x_z^1 x_x^2 x_y^2 x_z^2 x_x^3 x_y^3 x_z^3\n");
  PE("      where exponents indicate the local number of the node, and 999 are place holder values, not used.");
  PE("");
  PE("      WARNING: This is not properly implemented yet.");
  PE("");
  PE("    - segments in two dimensions in the x-y plane");
  PE("      If fault half-width is zero, then dip should be 90, z=0, and program will run in 2-D mode.");
  PE("      In this case, all z coordinates should still be specified (e.g. in the grid output requests)");
  fprintf(stderr,"      but z should be set to zero. Computation is plane-%s unless changed by -ps.\n",
	  (TWO_DIM_APPROX_IS_PLANE_STRESS_DEF)?("stress"):("strain"));
  PE("      For plane stress, all output values of u[Z] (displacements into stress free direction) are");
  PE("      given as strains e_zz.");
  fprintf(stderr,"      Calculations are performed in a %s plane by default, can be changed with the -hp option.\n",
	  (HALF_PLANE_DEF)?("half"):("full"));
  PE("");
  PE("    If all your patches are rectangular in 3-D, you could recompile\n    without the ALLOW_NON_3DQUAD_GEOM set to save space.");
#else
  PE("    If you would like to use the 2-D mode, point sources or triangles, recompile with ALLOW_NON_3DQUAD_GEOM flag set.");
#endif  
  PE("");
  PE("");
  PE("    Examples:");
  PE("");
  PE("    0 0 -1 30 90 4 2 0");
  PE("");
  PE("    for a rectangular, vertical fault patch centered at depth -1, striking 30 degrees from North");
  PE("    with length 8 and width 4, patch group number is 0.");
#ifdef ALLOW_NON_3DQUAD_GEOM 
  PE("");
  PE("    0 0 -1 30 90 -1 32 0");
  PE("");
  PE("    is the same kind of fault, only as a point source. Note that 32 is the fault area from above.");
  PE("");
  PE("    999 999 999 999 999 -1 -1 0   2.0 3.4641 0  -2.0 -3.4641 0  0 0 -4");
  PE("");
  PE("    would be for a triangular element with roughly the same geometry.");
#endif
  PE("");
  PE("    Please note that there are various programs to convert to and from the above patch format,");
  PE("    e.g. patch2xyz to go to GMT style xyz coordinates for the edges of each patch,");
  PE("    patch2vtk to go from patch format to (paraview) VTK format, or");
  PE("    points2patch to convert a set of four points in FE ordering to the above patch format.");
  PE("");
  PE("");
  
  fprintf(stderr,"(2) Boundary conditions and the operational mode are input via the \"%s\" file\n",
	  BC_FILE);
  PE("    which has the following format (free format ASCII list of parameters):\n");
  PE("    operational_mode:");
  fprintf(stderr,"      %i: one step calculation\n",
	  ONE_STEP_CALCULATION);
  fprintf(stderr,"      %i: loading simulation, and\n",
	  SIMULATE_LOADING);
  fprintf(stderr,"      %i: loading simulation and x window plotting.\n",
	  SIMULATE_LOADING_AND_PLOT);
  PE("");
  PE("");
  PE(" IF ONE STEP CALCULATION IS CHOSEN, the following input fields should be");
  PE("");
  PE("     print_bulk_fields");
  PE("");
  PE("      If print_bulk_fields is set to 0, no stress or displacement fields will be output.");
  PE("      If print_bulk_fields is set to 1, the program will read in the boundaries from the next line");
  PE("");
  PE("       xmin xmax n ymin ymax m zmin zmax o");
  PE("");
  PE("       where n, m, and o are the number of samples between the given limits.");
  PE("");
  PE("         If o is set to -1, the output will be in a plane with the");
  PE("          average strike and dip direction of fault group 0. x/y min/max are then the");
  PE("          limits in strike and dip direction.");
  PE("         If o is set to -2, the output will be in a plane with the");
  PE("          average strike and normal direction of fault group 0. x/y min/max are then the");
  PE("          limits in strike and normal direction.");
  PE("         In these cases, the coordinates in the output files will correspond to the local,");
  PE("         rotated systems, e.g. for -2, x will go along strike, and y along the normal direction");
  PE("         with respect to the average fault patch geometry, with origin in the center of the fault");
  PE("         group. HOWEVER, the displacements and stresses will still be given in the old system!");
  fprintf(stderr,"         If the -pc switch is set, the global coordinates will be printed to \"%s\"\n",
	  PLANE_COORD_FILE);
  PE("");
  fprintf(stderr,
	  "      If print_bulk_fields is set to 2, the program will read the output locations from file \"%s\"\n",
	  OLOC_FILE);
  PE("      This file has the x, y, and z coordinates in ASCII as an unformatted list.");
  PE("");
  PE("");
  PE("     patch_number boundary_code boundary_value");
  PE("");
  PE("      patch_number runs from 0 to N-1, where N is determined from the number of patches");
  PE("      (fully read lines, unformatted) in the geometry file.");
  PE("");
  PE("      If patch_number < 0, then the line should instead read");
  PE("");
  PE("       -increment start_fault stop_fault boundary_code boundary_value ... ...");
  PE("");
  PE("      and boundary code boundary_code will be assigned with value");
  PE("      boundary_value from start_fault to stop_fault with increments increment.");
  PE("");
  PE("      If start_fault < 0 and stop_fault < 0, will select all patches.");
  PE("");
  PE("     The patch_number ... boundary_value line can be repeated as often as necessary (say, twice for each fault");
  PE("     if the strike and dip movement modes are to be activated the same time).");
#ifdef ALLOW_NON_3DQUAD_GEOM
  PE("");
  fprintf(stderr,"     If the slip was specified (codes %i, %i, or %i) and one or more fault patches are\n",
	  STRIKE,DIP,NORMAL);
  PE("     triangular, then two more fields have to be given in the end (indicated by dots),");
  PE("     the strike and dip dip angle of the global reference frame to which the prescribed slip values");
  PE("     refer. The slip boundary conditions will then be converted to the local reference");
  PE("     frame of the triangular element.");
#endif
  PE("");
  PE("");
  PE("   Possible boundary conditions (as indicated by their integer codes) are:");
  PE("");
  fprintf(stderr,"   (A) %i, %i, %i: slip on fault specified (if called several times, will add up)\n",
	  STRIKE,DIP,NORMAL);
  fprintf(stderr,"    %i means strike, %i dip, and %i normal direction\n",
	  STRIKE,DIP,NORMAL);
  PE("");
  PE("    For strike: positive values of slip mean left-lateral fault, negative right-lateral;");
  PE("    \tthis corresponds to a resulting drop resp. increase in the strike component of stress.");
  PE("    \t(We use stress values with signs, that is a negative stress drop corresponds to a");
  PE("    \treduction of positive shear stress, a positive stress drop to reduction of shear");
  PE("    \tstress with negative sign. Both times, the absolute value of the shear stress is");
  PE("    \treduced).");
  PE("    For dip:    positive values of slip mean thrust fault, negative normal fault;");
  PE("    \tthis corresponds to a resulting drop resp. increase in the dip component of stress.");
  PE("    For normal: positive values of slip mean explosive source, negative implosive;");
  PE("    \tthis corresponds to a resulting drop resp. increase in the normal component of stress.\n");
  PE("");
  fprintf(stderr,"   (B) %i, %i, %i or %i, %i, %i or %i, %i, %i: stresses are specified\n",
	  STRIKE_SLIP,STRIKE_SLIP_LEFTLATERAL,STRIKE_SLIP_RIGHTLATERAL,
	  DIP_SLIP,DIP_SLIP_UPWARD,DIP_SLIP_DOWNWARD, 
	  NORMAL_SLIP,NORMAL_SLIP_OUTWARD,NORMAL_SLIP_INWARD);
  PE("");
  PE("    where the following codes indicate the degree of freedom for the faults to\n     achieve that stress starting from 0");
  fprintf(stderr,"    %2i\tstrike slip, either way\n",
	  STRIKE_SLIP);
  fprintf(stderr,"    %2i\tstrike slip, left lateral   (stress drop bc usually < 0 for activation)\n",
	  STRIKE_SLIP_LEFTLATERAL);
  fprintf(stderr,"    %2i\tstrike slip, right lateral  (stress drop bc usually > 0 for activation)\n",
	  STRIKE_SLIP_RIGHTLATERAL);
  fprintf(stderr,"    %2i\tdip direction, either way\n",
	  DIP_SLIP);
  fprintf(stderr,"    %2i\tdip direction, thrust       (stress drop bc usually < 0 for activation)\n",
	  DIP_SLIP_UPWARD);
  fprintf(stderr,"    %2i\tdip direction, normal       (stress drop bc usually > 0 for activation)\n",
	  DIP_SLIP_DOWNWARD);
#ifndef NO_OPENING_MODES
  fprintf(stderr,"    %2i\tnormal direction, either way\n",
	  NORMAL_SLIP);
  fprintf(stderr,"    %2i\tnormal direction, explosive (stress drop bc usually < 0 for activation)\n",
	  NORMAL_SLIP_OUTWARD);
  fprintf(stderr,"    %2i\tnormal direction, implosive (stress drop bc usually > 0 for activation)\n",
	  NORMAL_SLIP_INWARD);
#endif 
  PE("");
  fprintf(stderr,
	  "    If any of the strike or dip modes (%i, %i, %i or %i, %i, or %i) are given as \n",
	  STRIKE_SLIP,STRIKE_SLIP_LEFTLATERAL,STRIKE_SLIP_RIGHTLATERAL,
	  DIP_SLIP,DIP_SLIP_UPWARD,DIP_SLIP_DOWNWARD);
  fprintf(stderr,"    (%i, %i, %i or %i, %i, or %i) instead (add %i) this means\n",
	  COULOMB_STRIKE_SLIP,COULOMB_STRIKE_SLIP_LEFTLATERAL,COULOMB_STRIKE_SLIP_RIGHTLATERAL,
	  COULOMB_DIP_SLIP,COULOMB_DIP_SLIP_UPWARD,COULOMB_DIP_SLIP_DOWNWARD,OS_C_OFFSET);
	  
  PE("    that the target stress should be corrected by the change in normal");
  PE("    stress times the dynamic coefficient of friction (makes sense for Coulomb type rupture).");
  PE("");
  fprintf(stderr,"   (C) %i, %i, %i: read in background loading (see -l, NOT -s), reduce resolved stresses on faults\n",
	  SHEAR_STRESS_FREE, SHEAR_STRESS_FRICTION,SHEAR_STRESS_FREE_STRIKE_ONLY);
  PE("");
  fprintf(stderr,"    %2i\tsuch that faults are shear-stress free after slip in strike and dip directions\n",
	  SHEAR_STRESS_FREE);
  fprintf(stderr,"    %2i\tor such that the shear stress is given by the static friction coefficient\n",
	  SHEAR_STRESS_FRICTION);
  fprintf(stderr,"      \ttimes the normal stress, after slip in strike and dip direction\n");
  fprintf(stderr,"    %2i\tor such that faults are shear stress free after slip in strike direction only\n",
	  SHEAR_STRESS_FREE_STRIKE_ONLY);
  PE("    For the above two cases, the BC value (the last number) is meaningless.\n");
  PE("");
  PE("    NOTE: The final stresses, e.g. in flt.dat, after slip will correspond to the stress drops");
  PE("    needed to reduce the background load as specified above. The stresses will NOT include");
  PE("    the resolved background stresses since those are only used to initialize the problem");
  PE("    as boundary conditions. I.e. add the resolved background shear stress to the shear stress");
  fprintf(stderr,"    specified in flt.dat stress to get shear stress free conditions for the %2i case, for example.",
     SHEAR_STRESS_FREE);
  PE("");
  fprintf(stderr,"    For debugging, the resolved shear, dip, and normal stresses for each patch will be written to \"%s\".\n",
	  RES_STRESS_FILE);
  PE("");
  PE("    Mixed boundary conditions: If you want to mix specified slip boundary conditions with stress");
  PE("    boundary conditions, slip has to be completely specified first such that subsequent stress conditions");
  PE("    can take the stress changes induced by the slip into account.");
  PE("");
#ifdef ALLOW_NON_3DQUAD_GEOM
  PE("    WARNING: Stress boundary conditions will NOT work for point sources, and are");
  PE("    unreliable for triangular sources right now.");
#endif
  PE("");
  PE("    Example:");
  PE("    1 1");
  PE("    -5 5 51 -5 5 51 -1 -1 1");
  PE("    -1 -1 -1 10 1");
  PE("");
  PE("");
  PE(" IF LOADING SIMULATION IS CHOSEN (with or without X output)");
  PE("");
  PE("   the first three numbers are:");
  PE("");
  PE("    time_step print_interval stop_time");
  PE("");
  PE("   where time runs from 0 to stop_time in base time steps of time_step");
  PE("   and the program will print to fault files every print_interval intervals.");
  PE("");
  PE("   The time_step will normally be adjusted (but see -ct) based on the");
  PE("   program's prediction of the next rupture time or change in shear stress sign.");
  PE("   The max step is given by print_interval to ensure frequent\n   enough output to the flt.*.dat files.");
  fprintf(stderr,"   The minimum time step was set to %g.\n",MIN_TIME_STEP);
  PE("   Time has to monotonously increase from zero.");
  PE("");
  PE("   If X windows output is chosen, the next line (two numbers) should read\n");
  PE("    x_plot_interval x_scroll_interval\n");
  PE("   where x_plot_interval and x_scroll_intervals give the timing for the map view");
  PE("   and time vs. x plots, respectively.");
  PE("");
  PE("   All following lines assign active modes to faults and are of format:");
  PE("");
  PE("    patch_number activation_mode_code");
  PE("");
  PE("   where fault numbers should be given from 0 to N-1, N being the number of patches.");
  PE("   If patch_number < 0 then give");
  PE("");
  PE("    -increment start_fault stop_fault activation_mode_code");
  PE("");
  PE("   instead. If start_fault<0 and stop_fault<0 will select all patches, as above.");
  PE("");
  PE("   Activation_mode_codes are:");
  fprintf(stderr,"     %2i\t\tinactive(default)\n",
	  INACTIVE);
  fprintf(stderr,"     %2i\t\tstrike slip,   either way\n",
	  STRIKE_SLIP);
  fprintf(stderr,"     %2i\tstrike slip,   Coulomb correction, either way\n",
	  COULOMB_STRIKE_SLIP);
  fprintf(stderr,"     %2i\t\tstrike slip,   left lateral\n",
	  STRIKE_SLIP_LEFTLATERAL);  
  fprintf(stderr,"     %2i\tstrike slip,   Coulomb correction, left lateral\n",
	  COULOMB_STRIKE_SLIP_LEFTLATERAL);
  fprintf(stderr,"     %2i\t\tstrike slip,   right lateral\n",
	  STRIKE_SLIP_RIGHTLATERAL);
  fprintf(stderr,"     %2i\tstrike slip,   Coulomb correction, right lateral\n",
	  COULOMB_STRIKE_SLIP_RIGHTLATERAL);
  fprintf(stderr,"     %2i\t\tdip direction, either way\n",
	  DIP_SLIP);
  fprintf(stderr,"     %2i\tdip direction, Coulomb correction, either way\n",
	  COULOMB_DIP_SLIP);
  fprintf(stderr,"     %2i\t\tdip direction, thrust\n",
	  DIP_SLIP_UPWARD);
  fprintf(stderr,"     %2i\tdip direction, Coulomb correction, thrust\n",
	  COULOMB_DIP_SLIP_UPWARD);
  fprintf(stderr,"     %2i\t\tdip direction, normal\n",
	  DIP_SLIP_DOWNWARD);
  fprintf(stderr,"     %2i\tdip direction, Coulomb correction, normal\n",
	  COULOMB_DIP_SLIP_DOWNWARD);
#ifndef NO_OPENING_MODES
  fprintf(stderr,"     %2i\t\tnormal direction, either way (not implemented yet)\n",
	  NORMAL_SLIP);
  fprintf(stderr,"     %2i\t\tnormal direction, explosive  (not implemented yet)\n",
	  NORMAL_SLIP_OUTWARD);
  fprintf(stderr,"     %2i\t\tnormal direction, implosive  (not implemented yet)\n",
	  NORMAL_SLIP_INWARD);
#endif
  fprintf(stderr,"     %2i\t\tmaximum stress direction in strike/dip plane\n",
	  MAXSDIR_SLIP);
  fprintf(stderr,"     %2i\tmaximum stress direction in strike/dip plane, Coulomb correction\n",
	  COULOMB_MAXSDIR_SLIP);
  PE("");
  PE("   If Coulomb correction is set, the modes that involve slip in strike and dip direction will");
  PE("   result in a correction for normal stress changes, assuming that the target");
  PE("   stress drop is modified by the product of normal stress change during slip and");
  PE("   the dynamic friction coefficient.");
  PE("");
  PE("");
  PE(" NOTICE:");
#ifdef NO_OPENING_MODES
  PE(" Program was compiled with the NO_OPENING_MODES switch set, therefore all");
  PE(" normal (opening) modes of slip are suppressed to save memory");
  PE("");
#else
  PE(" Program maintains normal (opening) modes. If this is not needed, recompile");
  PE(" with the NO_OPENING_MODES switch set to save memory.");
#endif
  PE("");
  PE("");
  PE("PARAMETERS:");
  PE("");
  PE(" Elastic and frictional law properties are hard-coded and should be set in the \"properties.h\" file.");
  PE(" Whenever you want to change those, you will have to recompile the program. This was done to improve");
  PE(" the execution speed, assuming that material properties will not change between different model runs.");
  PE("");
  PE(" In the following, we will refer to two different matrices, A and I. The interaction, or I, matrix is");
  PE(" assembled once for the loading simulations and holds all possible stress influence coefficients for all");
  PE(" combinations of patch slips for a given geometry. I is typically held in memory and can be quite big, depending");
  PE(" on the problem. Hence, there is the option to store a reduced, sparse version of I. If patches are activated");
  PE(" (slipping) during a loading simulation (according to their frictional law), a system of equations"); 
  PE(" A x = b is set up to solve for the slip vector x that corresponds to the target stress drop. This A matrix");
  PE(" is also the only matrix that is assembled for a one-step calculation as we specifically know the patch/boundary");
  PE(" conditions we have to consider in this case. There are several possible modes for solving completely");
  PE(" unconstrained A x = b systems; constrained systems (where some entries in x have to be >=0 for unidirectionally");
  PE(" constrained slip) are always solved with a direct implementation of Lawson and Hanson's solver.");
  PE("");
  PE(" The following parameters can be set to non-default values on the command line:");
  PE("");
#ifdef NO_COHESION 
  PE(" -c option for cohesion was switched off, always assuming that cohesion C=0");
#else
  fprintf(stderr," -c  value      sets the cohesion term of the friction law to value, by default: %g\n",
	  COHESION_DEF);
  fprintf(stderr,"                If the cohesion=0 at all times, recompile with NO_COHESION set\n");
  fprintf(stderr,"                to improve speed. (Stress drop is set to %g, which is %g times mu).\n",
	  STRESS_DROP,STRESS_DROP/SHEAR_MODULUS);
#endif
  PE("");
  PE(" -pr value      set the background pressure to value \"value\".");
  fprintf(stderr,"                The background pressure amplitude is %g by default, both for one-step and loading\n",
	  PRESSURE_DEF);
#ifdef HYDROSTATIC_PRESSURE
  fprintf(stderr,"                simulations. The pressure varies with depth on length scales of %g (see properties.h).\n",
	  HYDROSTATIC_PRESSURE);
#else
  fprintf(stderr,"                simulations. The pressure is constant with depth (see properties.h).\n");
#endif
  PE("");
  PE(" -ms value      set minimum stress drop to value \"value\". (Characteristic stress drops are");
  fprintf(stderr,"                set to %g.) The minimum stress drop is zero by default, ie.\n",
	  STRESS_DROP);
  PE("                allowing stress drops of all sizes (loading simulation).");
  PE("");
  
  PE(" -ei value      set matrix cutoff value for sparse storage of I or A in absolute values");
  fprintf(stderr,"                             (%e by default.)\n",I_MAT_CUTOFF_DEF);
  PE("                (only used if I matrix is stored as sparse, see -nd option, or for sparse Ax=b solver.)");
  PE("");
  fprintf(stderr," -wc value      will use value as the maximum fraction of singular value for an SVD solution (default: %.2e)\n",SVD_THRESHOLD);
  PE("");
  PE("");
  PE("");
  PE("OPTIONS:");
  PE("");
  PE(" The following operational modes (switches) can be switched on or off by");
  PE(" using the toggle options listed below (type, eg, interact -s)");
  PE(" Please read carefully which switches are ON by default. If those");
  PE(" active switches are selected, they will be switched OFF instead.");
  PE("");
  
  fprintf(stderr," -s  read in initial stress values for each fault from file \"%s\", format:\n",
	  FAULT_STRESS_INIT_FILE);
  PE("");
  PE("     s_strike^0 s_dip^0 s_normal^0\n     s_strike^1 s_dip^1 s_normal^1\n     ...\n");
  PE("     (repeated for each patch). Here, the s_i are the initial values of the stress ");
  PE("     tensor resolved on the fault slip directions. All values are ADDED to the");
  PE("     homogeneous background stress as set in the program, ie. use all zeros for normal");
  PE("     initial stress state. (loading simulation)");
  PE(" --> WARNING: This option will only work for the loading simulation mode, for one-step calculations");
  PE("     specify fault stresses as individual boundary conditions, or use -l to get a resolved background");
  PE("     stress.");
  PE("")
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(READ_INITIAL_FAULT_STRESS_DEF),
	  name_boolean(TOGV(READ_INITIAL_FAULT_STRESS_DEF)));
  PE("");

  fprintf(stderr," -f  read in fault properties from file \"%s\", so far only parameters for the friction law\n",
	  FAULT_PROP_FILE);
  PE("     are supported. Format:\n");
  PE("     mu_static^0 mu_dynamic^0\n     mu_static^1 mu_dynamic^1\n     ...\n");
  PE("     (repeated for each fault). Static and dynamic coefficients of friction for");
  PE("     Coulomb.");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(READ_FAULT_PROPERTIES_DEF),name_boolean(TOGV(READ_FAULT_PROPERTIES_DEF)));
  PE("");

  PE(" -i  suppress all interactions between fault patches except the self-effect of slip on any");
  fprintf(stderr,"     patches within the group that the patch belongs to and the effect of slip\n");
  fprintf(stderr,"     on the number %i master group (this should be the main fault).\n",
	  MASTER_FAULT_GROUP);
  fprintf(stderr,"     (Meant as a experimental feature for LOADING SIMULATIONS ONLY!)\n");
  fprintf(stderr,"     %s by default, if switch is set will be %s. See -ni for suppressed interactions in one step.\n",
	  name_boolean(SUPPRESS_INTERACTIONS_DEF),name_boolean(TOGV(SUPPRESS_INTERACTIONS_DEF)));
  
  PE("");
  PE(" -ni suppress all interactions between faults except the interactions of fault patches");
  PE("");
  PE(" -ct use constant time steps for the loading simulation");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(CONSTANT_TIME_STEP_DEF),
	  name_boolean(TOGV(CONSTANT_TIME_STEP_DEF)));
  PE("");
  fprintf(stderr," -b  suppress output of all fault stress and slip data. Normally,\n     individual output files for each\n");
  fprintf(stderr,"     group of faults are generated when the number of groups is smaller than %i.\n",
	  MAX_NR_FLT_FILES_DEF);
  fprintf(stderr,"     If -b is set, this maximum number will be set to zero.\n");
  PE("");

  PE(" -p  prints information about the cumulative slip of every patch in a group to the");
  PE("     slipline.*.dat files if the fault files are used (see -b). (loading simulation)");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(PRINT_SLIPLINE_DEF),name_boolean(TOGV(PRINT_SLIPLINE_DEF)));
  PE("");

  fprintf(stderr," -w  whole fault activation mode (loading simulation)\n");
  fprintf(stderr,"     If set, one patch will activate all patches in the group for sliding.\n");
  fprintf(stderr,"     Also, if one patch is switched off because of low normal stress, all\n");
  fprintf(stderr,"     patches in the group will be switched off as well (*unless* -v is set).\n");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(WHOLE_FAULT_MODE_DEF),name_boolean(TOGV(WHOLE_FAULT_MODE_DEF)));
  PE("");

  fprintf(stderr," -k  keep slipping mode (loading simulation)\n");
  fprintf(stderr,"     If set, the first critical patches can trigger slip on other patches.\n");
  fprintf(stderr,"     If so, the initial patches will keep slipping, and the solution is obtained\n");
  fprintf(stderr,"     at each iteration until no more triggering takes place. If not set, slip\n");
  fprintf(stderr,"     on triggering patches is exhaustive, and iteration checks for triggering step\n");
  fprintf(stderr,"     by step (faster, since old solution doesn't have to be removed from stress\n");
  fprintf(stderr,"     field). The exhaustive method is affected by tiny differences in critical stress!\n");
  fprintf(stderr,"     If keep slipping, critical stress is %e, else %e.\n",
	  CRITICAL_STRESS_EPS,EXHAUSTIVE_CRITICAL_STRESS_EPS);
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(KEEP_SLIPPING_DEF),name_boolean(TOGV(KEEP_SLIPPING_DEF)));
  PE("");
  
  PE(" -v  whole fault deactivation mode (loading simulation)");
  PE("     If set (and -w is *not* set), the deactivation (low compressive normal stress) of a patch");
  PE("     will lead to a deactivation of the whole fault group.");
  PE("     If -w is on (whole fault activation mode), this is the default. If you use -v in this");
  PE("     case, the whole fault de-activations will be suppressed.");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(WHOLE_FAULT_DEACTIVATIONS_DEF),
	  name_boolean(TOGV(WHOLE_FAULT_DEACTIVATIONS_DEF)));
  PE("");

  fprintf(stderr," -l  reads a/b factors for the linear background loading stress relationship.\n");
  fprintf(stderr,"     If set, the file \"%s\" is read, it should hold 6 pairs of a_i b_i factors\n",
	  STRESS_RELATION_FILE);
  fprintf(stderr,"     where the *non-hydrostatic* background stress at time t is given by s_i=a_i+b_i*t\n");
  fprintf(stderr,"     s_i = {s_xx, s_xy, s_xz, s_yy, s_yz, s_zz}, so that the format is\n\n");
  PE("     a_xx b_xx");
  PE("     a_xy b_xy");
  PE("     ...");
  PE("     a_zz b_zz");
  PE("");
  if(PRESSURE_DEF != 0.0){
#ifdef HYDROSTATIC_PRESSURE
    fprintf(stderr,"     To these factors, a depth dependent hydrostatic pressure with length-scale %g and amplitude %g\n",
	    HYDROSTATIC_PRESSURE, PRESSURE_DEF);
#else
    fprintf(stderr,"     To these factors, a depth independent hydrostatic pressure of %g\n",PRESSURE_DEF);
#endif
    fprintf(stderr,"     will be added. You can change the amplitude of the pressure (e.g. to zero) with -pr.\n");
  }
  fprintf(stderr,"\n     If no factors are read in (if \"%s\" is not found), the code will use simple shear loading \n",
	  STRESS_RELATION_FILE);
#ifdef HYDROSTATIC_PRESSURE
  fprintf(stderr,"     with b_xy=%g plus the depth dependent hydrostatic pressure with length-scale %g and amplitude %g\n",
	  STRESSING_RATE,HYDROSTATIC_PRESSURE, PRESSURE_DEF);
#else
  fprintf(stderr,"     with b_xy=%g plus the depth independent hydrostatic pressure of %g.\n",STRESSING_RATE,PRESSURE_DEF);
#endif
  PE("\n     (Modify pressure amp. with -pr). If you want to specify resolved background");
  PE("      stresses in a one-step calculation, simply use the a_i factors to specify the stress ");
  PE("      at time t=0 but mind the hydrostatic part.");
  fprintf(stderr,"\n     %s by default, if switch is set will be %s.\n",
	  name_boolean(READ_STRESS_RELATION_FACTORS_DEF),
	  name_boolean(TOGV(READ_STRESS_RELATION_FACTORS_DEF)));
  PE("");
  
  fprintf(stderr," -nd switches to sparse matrix I storage, no matter what the I matrix size is.\n");
  fprintf(stderr,"     By default, the program goes from dense to sparse only if the\n     matrix is bigger than %g MB.\n",
	  IMAT_SPARSE_LIM);
  fprintf(stderr,"     The sparse matrix storage implies that matrix elements whose absolute\n     value is less than %7.3e\n",
	  I_MAT_CUTOFF_DEF);
  PE("     are discarded. This behavior can be changed by specifying");
  PE("     a different cutoff value using the -ei parameter option.\n     (Make sure to use -si if you want to keep the matrix.)");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(USE_SPARSE_STORAGE_DEF),
	  name_boolean(TOGV(USE_SPARSE_STORAGE_DEF)));
  PE("     note that sparse here refers to only the storage mode of the I matrix,");
  PE("     not the solution method of A.x=b. (see -ss).");
  PE("");
  PE(" -ci checks the interaction matrix for positive Coulomb stress feedback loops, once calculated");
  PE("     This may be useful for loading simulations.");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(CHECK_FOR_INTERACTION_FEEDBACK_DEF),
	  name_boolean(TOGV(CHECK_FOR_INTERACTION_FEEDBACK_DEF)));
  PE("");
  fprintf(stderr," -si saves the interaction (I) matrix during a loading simulation,\n");
  fprintf(stderr,"     if I is written to a file during runtime and not held in memory.\n     Else, \"%s\" and \"%s.hdr\" will be deleted.\n",
	  INTERACTION_MATRIX_FILE,INTERACTION_MATRIX_FILE);
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(SAVE_IMAT_DEF),
	  name_boolean(TOGV(SAVE_IMAT_DEF)));
  PE("");
  fprintf(stderr," -oi reads old interaction (I) matrix from files \"%s.dat\" and \"%s.hdr\".\n     Saves time when the geometry\n",
	  INTERACTION_MATRIX_FILE,INTERACTION_MATRIX_FILE);
  PE("     does not change between subsequent model runs (this is not automatically checked!)\n     (Make sure to use -si if you want to keep the matrix.)");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(USE_OLD_IMAT_DEF),
	  name_boolean(TOGV(USE_OLD_IMAT_DEF)));
  PE("");
  fprintf(stderr,
	  " -sa saves the A x = b A matrix to \"%s\" and \"%s.hdr\" during a one-step calculation.\n",
	  A_MATRIX_FILE,A_MATRIX_FILE);
  PE("     Binary format, with dimensions (n x m) and precision (p, in bytes) in the header file as n p m");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(SAVE_AMAT_DEF),
	  name_boolean(TOGV(SAVE_AMAT_DEF)));
  PE("");
  fprintf(stderr," -oa reads A matrix from files \"%s\" and \"%s.hdr\" during a one-step calculation.\n     Saves time when A\n",
	  A_MATRIX_FILE,A_MATRIX_FILE);
  PE("     does not change between subsequent model runs (this is not automatically checked!)\n     (Make sure to use -sa if you want to keep the matrix.)");
  PE("     NOTE that A is only constant for a constant geometry and the boundary conditions are");
  PE("     either pure slip or pure stress components.");
  PE("     If you are specifying a Coulomb stress change, elements of A will be modified");
  PE("     by the boundary conditions, and not only b as one might think!");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(USE_OLD_AMAT_DEF),
	  name_boolean(TOGV(USE_OLD_AMAT_DEF)));
  PE("");



  PE(" -r  attempts to restart a loading simulation model run based on previous results");
  PE("     WARNING: all result files will be overwritten NOT appended, so save before.");
  fprintf(stderr,"     For this to work, \"%s\" and \"%s\" (and possibly \"%s\" and \"%s\"\n     should not be substantially changed.\n\n",
	  GEOMETRY_FILE, BC_FILE,FAULT_STRESS_INIT_FILE,FAULT_PROP_FILE);
  fprintf(stderr,"     Slip events from the previous run in the format as in the output file \"%s\" will be read in\n",
#ifdef BINARY_PATCH_EVENT_FILE
	  EVENT_FILE_BINARY);
  fprintf(stderr,"     from the restart event file \"%s\"\n",
	  RESTART_EVENT_FILE_BINARY);
#else
          EVENT_FILE_ASCII);
  fprintf(stderr,"     from the restart event file \"%s\"\n",
	  RESTART_EVENT_FILE_ASCII);
#endif
  PE(""); 

  PE(" -sn will print NaN values in displacement and stress outputs of the one-step calculation.");
  PE("     Normally, those locations (end patch/segments) yielding one or more inf values in the");
  PE("     stress matrix or the displacement vector will be omitted during output.");
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(TOGV(SUPPRESS_NAN_OUTPUT_DEF)),
	  name_boolean(SUPPRESS_NAN_OUTPUT_DEF));
  PE("");
  PE(" -pc will print the global coordinates if fault-local planes are chosen for bulk stress");
  fprintf(stderr,"     and displacement output (-1 or -2 for o in zmin zmax o), writes xyz to \"%s\"\n",
	  PLANE_COORD_FILE);
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(TOGV(PRINT_PLANE_COORD_DEF)),
	  name_boolean(PRINT_PLANE_COORD_DEF));
  PE("");
  PE(" -sv uses SVD solver for the unconstrained A x = b system.");
  PE("     This approach should reliably find the least-squares minimum norm solution");
  PE("     if the A matrix is singular.");
  PE("     Note that LU can be orders of magnitude faster than SVD but use caution.");
  PE("     If non-negative solutions are sought (as for specified slip directions),");
  PE("     the program will automatically default back to the NNLS solver.");
  PE("");
  PE(" -sl uses LU solver for the unconstrained A x = b system (default).");
  PE("     LU can be orders of magnitude faster than SVD, but use caution.");
  PE("     If non-negative solutions are sought (as for specified slip directions),");
  PE("     the program will automatically default back to the NNLS solver.");
  PE("");
  PE(" -d  debug output, if you want extra checks, compile with -DDEBUG")
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(DEBUG_DEF),
	  name_boolean(TOGV(DEBUG_DEF)));
  PE("");
#ifdef USE_SUPERLU
  PE(" -ss uses LU solver and sparse storage/solution method for unconstrained A x = b system");
  PE("     WARNING:");
  PE("     This is meant for the one-step calculation only, and involves cutting off");
  PE("     small interaction matrix values as specified by -ei (bad idea?!)");
  PE("");
#endif
#ifdef USE_PETSC
  PE(" -fpetsc    force Petsc solvers even for serial runs (else LAPACK for LU)");
  PE("")
#endif
#ifdef ALLOW_NON_3DQUAD_GEOM
  fprintf(stderr," -ps Change the default 2-D elastic approximation for segments from plane-%s to plane-%s.\n",
	  (TWO_DIM_APPROX_IS_PLANE_STRESS_DEF)?("stress"):("strain"),
	  (!TWO_DIM_APPROX_IS_PLANE_STRESS_DEF)?("stress"):("strain"));
  fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	  name_boolean(TWO_DIM_APPROX_IS_PLANE_STRESS_DEF),
	  name_boolean(TOGV(TWO_DIM_APPROX_IS_PLANE_STRESS_DEF)));
 PE("");
 fprintf(stderr," -hp run a 2-D plane strain calculation in a half-plane (y<=0) instead of plane\n");
 fprintf(stderr,"     %s by default, if switch is set will be %s.\n",
	name_boolean(HALF_PLANE_DEF),
	name_boolean(TOGV(HALF_PLANE_DEF)));
#endif
  PE("");
  PE(" -h  prints out this help message and exits to the operating system");
  PE("");
  PE("(C) Thorsten Becker, thwbecker@post.harvard.edu, 1999 - 2023)");
  PE("    interact - boundary element code for elastic half-spaces");
  PE("    Main 3-D dislocation code based on dc3d.f by Y. Okada, as of Okada, BSSA, 1992");
  PE("    2-D segment slip solution from Crouch and Starfield (1973)");
  PE("    May include routines based on copyrighted software of others.");
  PE("    Distributed under the GNU public license, see \"COPYING\".");
  PE("");
}
#undef PE
/*


  called with code as integer, will return a description in words
  used by read_boundary_conditions

  these are the slip modes for loading simulation runs

*/
char *comment_on_code(short int code)
{
switch(code){
case INACTIVE:{return("inactive");}

case STRIKE_SLIP:{return("pure strike slip, unconstrained");}
case STRIKE_SLIP_LEFTLATERAL:{return("pure strike slip, leftlateral");}
case STRIKE_SLIP_RIGHTLATERAL:{return("pure strike slip, rightlateral");}

case DIP_SLIP:{return("pure dip slip, unconstrained");}
case DIP_SLIP_UPWARD:{return("pure dip slip, upward");}
case DIP_SLIP_DOWNWARD:{return("pure dip slip, downard");}

case NORMAL_SLIP:{return("pure normal slip, unconstrained");}
case NORMAL_SLIP_OUTWARD:{return("pure normal slip, explosion");}
case NORMAL_SLIP_INWARD:{return("pure normal slip, implosion");}

case MAXSDIR_SLIP:{return("mix of strike and dip slip");}

case COULOMB_STRIKE_SLIP:{return("pure strike slip, Coulomb corr.");}
case COULOMB_STRIKE_SLIP_LEFTLATERAL:{return("pure strike slip, leftlateral, Coulomb corr.");}
case COULOMB_STRIKE_SLIP_RIGHTLATERAL:{return("pure strike slip, rightlateral, Coulomb corr.");}

case COULOMB_DIP_SLIP:{return("pure dip slip, Coulomb corr.");}
case COULOMB_DIP_SLIP_UPWARD:{return("pure dip slip, upward, Coulomb corr.");}
case COULOMB_DIP_SLIP_DOWNWARD:{return("pure dip slip, downard, Coulomb corr.");}
case COULOMB_MAXSDIR_SLIP:{return("mix of strike and dip slip, Coulomb corr.");}

  default:{return("undefined");}
  }
}
/*
  
  same as above but for one-step boundary conditions

 */
char *comment_on_code_bc(short int code, COMP_PRECISION bc)
{
  switch(code){
    // slip codes
  case STRIKE:{return((bc>=0)?("left-lateral transform"):("right-lateral transform"));}
  case DIP:{return((bc>0)?("thrust fault"):("normal fault"));}
#ifndef NO_OPENING_MODES
  case NORMAL:{return((bc>=0)?("explosive source"):("implosive source"));}
#endif
  case STRIKE_SLIP_LEFTLATERAL:{return("left-lateral strike-slip stress");}
  case STRIKE_SLIP_RIGHTLATERAL:{return("right-lateral strike-slip stress");}
  case STRIKE_SLIP:{return("general strike-slip stress");}
  case COULOMB_STRIKE_SLIP_LEFTLATERAL:{return("left-lateral strike-slip stress, Coulomb corr.");}
  case COULOMB_STRIKE_SLIP_RIGHTLATERAL:{return("right-lateral strike-slip stress, Coulomb corr.");}
  case COULOMB_STRIKE_SLIP:{return("general strike-slip stress, Coulomb corr.");}
  case DIP_SLIP_UPWARD:{return("dip-slip thrust stress");}
  case DIP_SLIP_DOWNWARD:{return("dip-slip normal stress");}
  case DIP_SLIP:{return("dip-slip general stress");}
  case COULOMB_DIP_SLIP_UPWARD:{return("dip-slip thrust stress, Coulomb corr.");}
  case COULOMB_DIP_SLIP_DOWNWARD:{return("dip-slip normal stress, Coulomb corr.");}
  case COULOMB_DIP_SLIP:{return("dip-slip general stress, Coulomb corr.");}
#ifndef NO_OPENING_MODES
  case NORMAL_SLIP_OUTWARD:{return("normal explosive stress");}
  case NORMAL_SLIP_INWARD:{return("normal implosive stress");}
  case NORMAL_SLIP:{return("normal stress");}
#endif
  case SHEAR_STRESS_FREE:{return("shear stress  from ext. load, strike and dip");}
  case SHEAR_STRESS_FREE_STRIKE_ONLY:{return("shear stress  from ext. load, strike only");}
  case SHEAR_STRESS_FRICTION:{return("friction modified shear stress from ext. load");}
  default:{return("undefined");}
  }
}



