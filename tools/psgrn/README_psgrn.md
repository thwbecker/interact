# PSGRN/PSCMP for interact: reference solutions and 3-D kernel generation

`install_psgrn [dir] [fc]` clones and compiles Wang's PSGRN/PSCMP
2020 (github.com/RongjiangWang/PSGRN-PSCMP_2020) with gfortran and
smoke-tests the binaries.  Two roles in this repository:

1. REFERENCE SOLUTIONS.  PSGRN/PSCMP is the community standard for
   co-/postseismic deformation of layered viscoelastic-gravitational
   half-spaces.  test_relax/inplane_ve_proto/psgrn_compare.py
   cross-validates the in-plane plate-over-Maxwell(+gravity) kernel
   machinery against it (postseismic uplift-change profiles agree to
   1-3 percent of peak, gravity on and off; input templates in
   test_relax/inplane_ve_proto/psgrn_inputs/).  A long PSCMP fault
   evaluated at mid-strike gives the 2-D plane-strain limit.

2. 3-D KERNEL GENERATION for rsf_solve -ve_mode 3 (-ve_prony_file).
   The file contract is generator-agnostic: shared relaxation times
   tau_p plus per-pair amplitude matrices C_p with K(t) = C_inf +
   sum_p C_p exp(-t/tau_p) and C_inf implicit (the assembled elastic
   operator stays exact), plus a high-k elastic K0 block for the
   runtime consistency gate against Is (off-diagonal entries;
   interact assembles Is in the stress-DROP convention, so kernels
   from stress-from-slip codes need a global sign flip).  The 2-D
   reference implementation of the whole pipeline is
   test_relax/inplane_ve_proto/bp3_ve_kernels.py: sample the
   RELAXATION part dK_ij(t) = K_ij(t) - K_ij(0) per source patch
   (differences cancel the singular direct terms), fit fixed-ladder
   amplitudes with the relaxed limit as a constraint and a held-out
   sample as the gate, write the file.  For 3-D: replace the 2-D
   sampler with PSGRN/PSCMP runs per source patch (pscmp gives
   stress tensors at receiver patches for a slip step; resolve onto
   receiver slip directions), keep everything else identical.
   Caveats for that step: PSGRN's sources are POINT dislocations
   interpolated from a depth/distance grid, so sample at receiver
   separations of a few patch sizes and let the assembled Is carry
   the near field (exactly what the implicit C_inf and the
   K0-gate's near-diagonal skip are designed for); and check the
   Burgers/rake conventions against the K0 gate before trusting any
   physics.
