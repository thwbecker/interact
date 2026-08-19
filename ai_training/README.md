# hbi_pair_tests: HBI M_vs2square / pair000 on a triangular mesh

Runs the HBI "pair000" scenario (20 x 10 km vertical strike slip fault,
2500 triangular cells, smoothed a field with two VS squares inside the VW
zone) with rsf_solve, following the layout and conventions of
seas_bp_tests.  The purpose is a code comparison against the HBI reference
run of the same case, so the setup mirrors the HBI configuration for this
particular scenario; other cases or configurations may need changes.

## Files

    pair000/geometry.stl    HBI fault mesh (from the HBI case directory)
    pair000/param.dat       HBI per-cell parameters (same source)
    pair000/stl2geom.awk    ASCII STL -> interact triangular patch format,
                            km -> m scaling, group 1 for a < b (VW) when
                            given the parameter file
    pair000/param2rsf.awk   HBI parameter file -> rsf_ab (a b, -rsf_file)
                            and ic (tau[Pa] v[m/s], -rsf_ic_file); FAILS if
                            any column assumed constant actually varies
    pair000/prepare_files   runs both converters plus consistency checks
    create_run              as in seas_bp_tests: model name picks the
                            evolution law and imex, e.g.
                            ./create_run pair000_prz_imex 0 8
    run_many                loops create_run over a set of laws

## Usage

    cd pair000; ./prepare_files; cd ..
    ./create_run pair000 0 8        # aging law, dense matrix, 8 cores

## Notes and caveats, as of 2026-08-18

- Mesh: the STL is a 50 x 25 grid of 0.4 x 0.4 km quads, each split into
  two congruent right triangles (legs 400 m, sqrt(area) 283 m), not the
  "50x50 cells of 0.4 x 0.2 km" of the HBI case README.  The per triangle
  area matches a 0.4 x 0.2 cell, which may be where that shorthand comes
  from.
- Constants: dc 0.018, f0 0.6, sigma 58 MPa, vpl 1e-9 m/s, v0 1e-6
  (HBI vref default), G 32.04 GPa and vs 3464 m/s (HBI m_const.f90
  defaults, not set in perf.in).  param2rsf.awk verifies the constancy of
  everything it does not convert and prints the flag values it finds.
- ICs: tau is per patch steady state at vel = vpl (param2rsf.awk checks
  this against the plain log law to ~3e-7 MPa; HBI writes the regularized
  arcsinh form plus radiation damping, both negligible at 1e-9 m/s).
  There is no forced nucleation patch, so in both codes first events grow
  from roundoff scale perturbations and exact first event times should not
  be expected to match; event statistics are the comparison target.
- Resolution: L_b = G dc / (b sigma) = 497 m against 283-400 m element
  dimensions, i.e. roughly one element per cohesive zone, and
  h* (RA, largest b - a = 0.01) = 1266 m, roughly 3-4 elements.  This is
  the discretization the HBI reference run used, so the comparison is like
  for like, but results at this resolution should not be treated as
  converged for either code without a refinement check.
- tmax defaults to 1000 yr to match the HBI perf.in for this case.
- Groups: 0 = VS (1896), 1 = VW (604, includes the smoothed VW rim), for
  -rsf_monitor_by_group.  The two VS converted squares inside the VW zone
  are not separable from the outer VS region by parameter values alone and
  are not tagged separately.
- Signs: all facet normals point to -y, giving local strike +x uniformly;
  HBI rake 0 slip sense conventions may differ from interact by a global
  sign, so compare |v| and slip magnitudes.
