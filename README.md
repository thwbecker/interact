# interact

`interact` models fault interactions using dislocations in an elastic medium.

`interact` uses `dc3d.f` as provided by Y. Okada as in Okada (BSSA,
1992) and  linear algebra routines from SLATEC, LAPACK, and
EISPACK. Might contain material copyrighted by others (e.g. Numerical
Recipes). Petsc implementation based on Dave May's examples and
assistance. 

See files `INSTALLATION`, `README.md`, and `help.txt` for documentation and `COPYRIGHT` and `COPYING` for the license and warranty disclaimers.

If you use `interact` please cite the following

>@Article{becker02c,                         
>  author = {Becker, T. W. and Schott, B.},
>  title = {On boundary-element models of elastic fault interaction (abstract)},
>    journal = {EOS Trans. AGU},
>    year = 2002,
>    volume = 83,
>    number = 47,
>    pages = {NG62A-0925}
>  }

Copyright (C) Thorsten W. Becker 2000 - 2023

thwbecker@post.harvard.edu

## Overview

This directory contains source codes needed to compile the elastic 2-D and half-space dislocation program `interact` as well as several tools that go with it. 

To compile, edit the `makefile` as described in the INSTALLATION document and type `make all`. The normal make might give you problems, the GNU gmake is therefore to be preferred.

General settings as to the operation of interact are set in `makefile`, material properties such as the shear modulus are set in `properties.h`, and all system-dependent settings should be changed in `makefile.gcc` or similar.  The reasoning behind this setup is that you are not likely to change the, say, shear modulus very often because it simply affects the scaling. Likewise, you will also not want to change the way interact works that frequently, so we hardwire several options to avoid calls to IF statements for efficiency. After compilation, the operation of interact is affected by command line arguments and input files.

For documentation, use `program -h` after compilation, and refer to the numerous comments in the source code as well as in the makefiles. The output of `interact -h` is supplied as the file `help.txt`.

The `example.?.txt` files describe some simple example calculations to, hopefully, elucidate the usage of interact and affiliated programs.
See our [Fall 2002 AGU poster](http:///www-udc.ig.utexas.edu/external/becker/agu02_disp.pdf).

Typing `make` will compile most programs, see `makefile`. If you want all tools, type `make really_all`. 

A list of the programs that will be compiled and several batch scripts that are provided in addition follows with short descriptions what they do. For each program, more information can be obtained by typing `program -h`. The input formats etc. are described when executing the command

```
interact -h
```

## Parameters

Elastic and frictional law properties are hard-coded and should be set in the `properties.h` file. Whenever you want to change those, you will have to recompile the program. This was done to improve the execution speed, assuming that material properties will not change between different model runs.

 In the following, we will refer to two different matrices, $\mathbf A$ and $\mathbf I$. The interaction, or $\mathbf I$, matrix is assembled once for the loading simulations and holds all possible stress influence coefficients for all combinations of patch slips for a given geometry. $\mathbf I$ is typically held in memory and can be quite big, depending on the problem. Hence, there is the option to store a reduced, sparse version of $\mathbf I$. If patches are activated (slipping) during a loading simulation (according to their frictional law), a system of equations $\mathbf A \mathbf x = \mathbf b$ is set up to solve for the slip vector x that corresponds to the target stress drop. This $\mathbf A$ matrix is also the only matrix that is assembled for a one-step calculation as we specifically know the patch/boundary conditions we have to consider in this case. There are several possible modes for solving completely unconstrained $\mathbf A \mathbf x = \mathbf b$ systems; constrained systems (where some entries in $\mathbf x$ have to be >=0 for unidirectionally constrained slip) are always solved with a direct implementation of Lawson and Hanson's solver.

 The following parameters can be set to non-default values on the command line:

* `-c value` Sets the cohesion term of the friction law to `value`. Default value is 5. If the cohesion=0 at all times, recompile with NO_COHESION set to improve speed. Stress drop is set to 1, which is $0.0001 \times \mu$.

* `-pr value` Set the background pressure to value given by `value`. The background pressure amplitude is 10 by default, both for one-step and loading simulations. The pressure is constant with depth (see `properties.h`).

* `-ms value` Sets the minimum stress drop to value `value`. Characteristic stress drops are set to 1. The minimum stress drop is zero by default, ie. allowing stress drops of all sizes (loading simulation).

* `-ei value` Set matrix cutoff value for sparse storage of $\mathbf I$ or $\mathbf A$ in absolute values Default is 5.000000e-05. Only used if I matrix is stored as sparse, see `-nd` option, or for sparse $\mathbf A \mathbf x = \mathbf b$ solver.

* ` -wc value` Will use `value` as the maximum fraction of singular value for an SVD solution. Default value is 7.00e-15.


## Options

 The following operational modes (switches) can be switched on or off by using the toggle options listed below (type, eg, `interact -s`) Please read carefully which switches are ON by default. If those active switches are selected, they will be switched OFF instead.

*  `-s`  Read in initial stress values for each fault from file `fsi.in`, format:
  ```
  s_strike^0 s_dip^0 s_normal^0
   s_strike^1 s_dip^1 s_normal^1
   ...  
  ```
        
   (repeated for each patch). Here, the `s_i` are the initial values of the stress  tensor resolved on the fault slip directions. All values are ADDED to the homogeneous background stress as set in the program, ie. use all zeros for normal initial stress state. (loading simulation)
   
   **WARNING** This option will only work for the loading simulation mode, for one-step calculations specify fault stresses as individual boundary conditions, or use `-l` to get a resolved background stress.

     ON by default, if switch is set will be OFF.

*  `-f`  Read in fault properties from file `fp.in`, so far only parameters for the friction law are supported. Format:
  ```
     mu_static^0 mu_dynamic^0
     mu_static^1 mu_dynamic^1
     ...
  ```
     (repeated for each fault). Static and dynamic coefficients of friction for Coulomb.
     
     ON by default, if switch is set will be OFF.

* `-i`  Suppress all interactions between fault patches except the self-effect of slip on any patches within the group that the patch belongs to and the effect of slip on the number 0 master group (this should be the main fault).
     
     **Meant as a experimental feature for LOADING SIMULATIONS ONLY!**
     
     OFF by default, if switch is set will be ON. See `-ni` for suppressed interactions in one step.

* `-ni` Suppress all interactions between faults except the interactions of fault patches.

* `-ct` use constant time steps for the loading simulation
     
     OFF by default, if switch is set will be ON.

* `-b`  Suppress output of all fault stress and slip data. Normally,
     individual output files for each group of faults are generated when the number of groups is smaller than 50. If `-b` is set, this maximum number will be set to zero.
     
*  `-p`  Print information about the cumulative slip of every patch in a group to the `slipline.*.dat` files if the fault files are used (see `-b`). (loading simulation)

     OFF by default, if switch is set will be ON.

* `-w`  Whole fault activation mode (loading simulation).
If set, one patch will activate all patches in the group for sliding. Also, if one patch is switched off because of low normal stress, all patches in the group will be switched off as well (*unless* `-v` is set).

     OFF by default, if switch is set will be ON.
     
* `-k`  Keep slipping mode (loading simulation).
If set, the first critical patches can trigger slip on other patches. If so, the initial patches will keep slipping, and the solution is obtained at each iteration until no more triggering takes place. If not set, slip on triggering patches is exhaustive, and iteration checks for triggering step by step (faster, since old solution doesn't have to be removed from stress field). The exhaustive method is affected by tiny differences in critical stress! If keep slipping, critical stress is -7.000000e-15, else -7.000000e-15.
     
     ON by default, if switch is set will be OFF.

* `-v`  Whole fault deactivation mode (loading simulation).
     If set (and `-w` is **not** set), the deactivation (low compressive normal stress) of a patch will lead to a deactivation of the whole fault group. If `-w` is on (whole fault activation mode), this is the default. If you use `-v` in this case, the whole fault de-activations will be suppressed.
     
     OFF by default, if switch is set will be ON.

* `-l`  Reads a/b factors for the linear background loading stress relationship.
     If set, the file `smlin.in` is read, it should hold 6 pairs of a_i b_i factors where the *non-hydrostatic* background stress at time t is given by s_i=a_i+b_i*t,  `s_i = {s_xx, s_xy, s_xz, s_yy, s_yz, s_zz}`, so that the format is
  ```
     a_xx b_xx
     a_xy b_xy
     ...
     a_zz b_zz
  ```
     To these factors, a depth independent hydrostatic pressure of 10
     will be added. You can change the amplitude of the pressure (e.g. to zero) with `-pr`.

     If no factors are read in (if `smlin.in` is not found), the code will use simple shear loading with b_xy=1 plus the depth independent hydrostatic pressure of 10.

   Modify pressure amp. with `-pr`. If you want to specify resolved background stresses in a one-step calculation, simply use the a_i factors to specify the stress at time t=0 but mind the hydrostatic part.

     ON by default, if switch is set will be OFF.

* `-nd` Switches to sparse matrix $\mathbf I$ storage, no matter what the $\mathbf I$ matrix size is.
     By default, the program goes from dense to sparse only if the matrix is bigger than 8000 MB. The sparse matrix storage implies that matrix elements whose absolute value is less than 5.000e-05 are discarded. This behavior can be changed by specifying a different cutoff value using the `-ei` parameter option. (Make sure to use `-si` if you want to keep the matrix.)
     
     OFF by default, if switch is set will be ON.
     
     Note that sparse here refers to only the storage mode of the $\mathbf I$ matrix, not the solution method of $\mathbf A \mathbf x = \mathbf b$. (see `-ss`).

* `-ci` Checks the interaction matrix for positive Coulomb stress feedback loops, once calculated. This may be useful for loading simulations.

     OFF by default, if switch is set will be ON.

* `-si` Saves the interaction ($\mathbf I$) matrix during a loading simulation, if $\mathbf I$ is written to a file during runtime and not held in memory. Else, `/tmp/i` and `/tmp/i.hdr` will be deleted.

     OFF by default, if switch is set will be ON.

* `-oi` Reads old interaction ($\mathbf I$) matrix from files `/tmp/i.dat` and `/tmp/i.hdr`.
     Saves time when the geometry does not change between subsequent model runs (this is not automatically checked!). Make sure to use `-si` if you want to keep the matrix.)
     
     OFF by default, if switch is set will be ON.

* `-sa` Saves the $\mathbf A \mathbf x = \mathbf b$ $\mathbf A$ matrix to `a.dat` and `a.dat.hdr` during a one-step calculation.
     Binary format, with dimensions (n x m) and precision (p, in bytes) in the header file as n p m
     
     OFF by default, if switch is set will be ON.

* `-oa` Reads $\mathbf A$ matrix from files `a.dat` and `a.dat.hdr` during a one-step calculation.
     Saves time when $\mathbf A$ does not change between subsequent model runs (this is not automatically checked!). Make sure to use `-sa` if you want to keep the matrix.
     NOTE that $\mathbf A$ is only constant for a constant geometry and the boundary conditions are either pure slip or pure stress components. If you are specifying a Coulomb stress change, elements of $\mathbf A$ will be modified by the boundary conditions, and not only $\mathbf b$ as one might think!
     
     OFF by default, if switch is set will be ON.

* `-r`  Attempts to restart a loading simulation model run based on previous results.
     **WARNING** All result files will be overwritten NOT appended, so save before. For this to work, `geom.in` and `bc.in` (and possibly `fsi.in` and `fp.in`should not be substantially changed.

     Slip events from the previous run in the format as in the output file `events.bin` will be read in from the restart event file `restart.events.bin`.

* `-gv` Output files with slip on faults in the Geometry center's geomview COFF format, with a CQUAD file for each fault group. Slip values are scaled to the global maxima. Note that the `fltdat2cquad` script can produce CQAUD files from `geom.in` and flt.dat files.

     OFF by default, if switch is set will be ON.

* `-sn` Print NaN values in displacement and stress outputs of the one-step calculation.
     Normally, those locations (end patch/segments) yielding one or more inf values in the stress matrix or the displacement vector will be omitted during output.
     
     OFF by default, if switch is set will be ON.

* `-pc` Print the global coordinates if fault-local planes are chosen for bulk stress and displacement output (-1 or -2 for o in zmin zmax o), writes xyz to `plane.xyz`.

     ON by default, if switch is set will be OFF.

* `-sv` Use SVD solver for the unconstrained $\mathbf A \mathbf x = \mathbf b$ system.
     This approach should reliably find the least-squares minimum norm solution if the $\mathbf A$ matrix is singular. Note that LU can be orders of magnitude faster than SVD but use caution. If non-negative solutions are sought (as for specified slip directions), the program will automatically default back to the NNLS solver.

* `-sl` Use LU solver for the unconstrained $\mathbf A \mathbf x = \mathbf b$ system (default).
     LU can be orders of magnitude faster than SVD, but use caution. If non-negative solutions are sought (as for specified slip directions), the program will automatically default back to the NNLS solver.

* `-d`  Run in debug mode.

     OFF by default, if switch is set will be ON.

* ` -h`  Print this help message and exits.

