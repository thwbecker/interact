################################################################################
#                                                                              #
#  interact: model fault interactions using dislocations in an elastic         #
#            medium                                                            #
#	                                                                       #
#  Copyright (C) Thorsten W. Becker 2000 - 2021                                #
#                                                                              #
#  Interact uses dc3d.f as provided by Y. Okada as in Okada (BSSA, 1992) and   #
#  linear algebra routines from SLATEC, LAPACK, and EISPACK.                   #
#  Might contain material copyrighted by others (e.g. Numerical Recipes etc)   #
#                                                                              #
#  See INSTALLATION, README.TXT, and help.txt for documentation and            #
#  COPYRIGHT and COPYING for the license and warranty disclaimers.             #
#                                                                              #
#  If you use interact in real work, a citation to use would be:               #
#                                                                              #
#  @Article{becker02c,                                                         #
#    author = 	 {Becker, T. W. and Schott, B.},                               #
#    title = 	 {On boundary-element models of elastic fault interaction      #
#                  (abstract)},                                                # 
#    journal = 	 {EOS Trans. AGU},                                             #
#    year = 	 2002,                                                         #
#    volume =	 83,                                                           #
#    number =	 47,                                                           #
#    pages =	 {NG62A-0925}                                                  #
#  }                                                                           #
#                                                                              #
#  thwbecker@post.harvard.edu                                                  #
#                                                                              #
################################################################################


This directory contains source codes needed to compile the elastic 2-D
and half-space dislocation program "interact" as well as several tools
that go with it. 

To compile, edit the makefile as described in the INSTALLATION
document and type 'make all'. (The normal make might give you
problems, the GNU gmake is therefore to be preferred.) 

General settings as to the operation of interact are set in
"makefile", material properties such as the shear modulus are set in
"properties.h", and all system-dependent settings should be changed in
"makefile.$ARCH".  The reasoning behind this setup is that you are not
likely to change the, say, shear modulus very often because it simply
affects the scaling. Likewise, you will also not want to change the
way interact works that frequently, so we hardwire several options to
avoid calls to IF statements for efficiency. After compilation, the
operation of interact is affected by command line arguments and input
files.

For documentation, use 'program -h' after compilation, and refer to
the numerous comments in the source code as well as in the
makefiles. The output of 'interact -h' is supplied as the file
'help.txt'.

The example.?.txt files describe some simple example calculations to,
hopefully, elucidate the usage of interact and affiliated programs.
(See also our Fall 2002 AGU poster at
http://www.seismology.harvard.edu/~becker/agu02_disp.pdf).

Typing 'make' will compile most programs, see the makefile. If you
want all tools, type 'make really_all'. 


A list of the programs that will be compiled and several batch scripts
that are provided in addition follows with short descriptions what
they do. For each program, more information can be obtained by typing
'program -h'.


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













