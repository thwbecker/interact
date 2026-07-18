/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

This package provides:
- C++ routines and a program to create an H-matrix approximation to a matrix.
- C++ routines to compute matrix-vector products (MVP) and related operations
  using this H-matrix approximation.
- An optional Matlab interface to create the input files to build the H-matrix.
- An optional Matlab interface to compute MVP.
- A limited C interface to the MVP routines.
- A limited Fortran interface to the MVP routines.
- A limited C++ MPI interface to the MVP routines.


Installation
------------
- If you are working in a Unix-like developer environment, follow these
  instructions. I don't have anything for Windows or Macs that lack Unix-like
  developer tools. But you can still follow the Makefile as a guide to how to
  write the appropriate build script for your environment.
- Edit the top of Makefile with compiler and LAPACK/BLAS information. Set the
  optimization level (opt =) and serial/parallel mode (mode =, with
  (s)erial, OpenMP (omp), and MPI (mpi)).
- On the command line, type 'make'.
- (Optional.) Start Matlab and cd to 'matlab/'. Type 'make'. It seems Macs have
  a problem with recent versions of Matlab; follow these instructions:
      http://www.mathworks.com/matlabcentral/answers/94092


Usage
-----
- See examples/ex.m for Matlab usage. ('addpath examples/' to access it in
  Matlab.) Even if you don't intend to use hmmvp in Matlab, reading through ex.m
  can be useful.
- In Matlab, type
      addpath matlab; help hmmvp;
- In a terminal, type
      ./bin/hmmvpbuild_mode help
  and
      ./bin/hmmvpbuild_mode help compress
  where mode is s, omp, or mpi from Makefile.
- As explained in greater detail in ex.m and the hmmvpbuild help, you need to
  create a key-value file to describe your problem. This can be done in Matlab
  using matlab/kvf.m (no mex file needed) or in C++ using
  util/include/KeyValueFile.hpp
- See examples/mvp_*.c* for example matrix-vector product usage.
- To add a new Green's function, follow the instructions in src/hmmvpbuild.cpp
  tagged with 'DEV'. Use src/GreensFnInverseR.cpp as an example.
- To apply an H-matrix, link against lib/libhmmvp_*.a, lapack, and blas. Follow
  the example in the Makefile target 'mvp'.


You may find my report "H-Matrix and Block Error Tolerances"
    http://arxiv.org/abs/1110.2807
helpful when setting the error tolerance.


If you use this software in research that you publish, I would appreciate your
citing this software package as follows:

    A. M. Bradley, Software for efficient static dislocation-traction
    calculations in fault simulators, Seis. Res. Lett. 85(6), 2014.

Feel free to contact me at (my permanent email address) ambrad@cs.stanford.edu
to check for a new version of this software. In addition, current versions are
maintained on these sites:
    pangea.stanford.edu/research/CDFM/software
    cs.stanford.edu/~ambrad


Limitations
-----------
- Not all transpose MVP functionality is implemented.
- The MPI code currently assumes all machines have the same endianness; I should
  use Boost.MPI. This affects only the compressor code.


Acknowledgements
----------------
Paul Segall, my postdoc advisor, and I gratefully acknowledge support from the
National Science Foundation (EAR-0838267), the U.S. Geological Survey
(G12AP20030), and SCEC (13096). I would also like to thank Paul for his support
of this work. (Any errors are mine.)

Thanks also to early users for their feedback, bug reports, and help with
platform issues: Michael Barall, Sylvain Barbot, Junichi Fukuda, Kaj Johnson,
Mark McClure.
