# Replication tests: Kato (EPS 2002) and Miyake & Noda (EPS 2019)

Driver: `run_noda_test` (about 6-8 minutes; gates below, then renders
`noda_kato02_test.png` and `noda_mn19_test.png` in `tmp_noda/`).
Generators: `gen_kato02.py`, `gen_mn19.py`.  Independent arbiters:
`kato02_exact_chain.py`, `mn19_exact_chain.py` (see below).  The broader
context of all antiplane VE validations is in `README_ve_antiplane.md`.

## The two models

Kato (2002): vertical strike-slip fault cutting an elastic layer
(h = 20/15/12 km) over a Maxwell half-space, backslip loading at
Vpl = 35 mm/yr, composite Kato-Tullis friction (state_law 5, Vc =
0.01 um/s), sigma_eff = (rho - rho_w) g z, G = 30 GPa, c_s = 3.27 km/s,
L = 5 cm, a(z) and a-b(z) DIGITIZED from his Fig. 2 (see gen_kato02.py
header; the digitization was pixel-measured and tick-calibrated from
the article PDF, with an overlay QC figure delivered outside the repo).
His kernel (his Eq. 2) is exactly the Bonafede et al. (1984) two-family
image construction with Erlang time weights that this repo's layered
machinery implements; his tr = 2 eta/G = 2x our -ve_tmaxwell_yr.
NOTE his Eq. (3) as printed starts the sum at l = 1 (the Erlang CDF
requires l = 0); integrating that misprint changes cycles by ~6 percent
and does not otherwise matter (mode `kato3` of the chain).

Miyake & Noda (2019): rate-weakening patch (2R = 1 km, a = 0.02,
a/b = 0.8 smoothed boxcar, b = 0 outside per their Fig. 2, aging law,
sigma = 100 MPa, f0 = 0.6 at V0 = Vpl = 1e-9 m/s) on an antiplane fault
in a uniform Maxwell medium, constant far-field stress (equivalent to
backslip), t_c scanned against t_load = 2 a sigma R/(mu Vpl) = 2.11 yr.
Replica approximations: quasi-dynamic (they are fully dynamic), buried
finite strip (16R, welded exterior, free surface 52R away) instead of
their 4.096R-periodic infinite fault; `gen_mn19.py [.. npatch]` places
npatch patches every 4.096R to approximate the periodicity.

## Independent arbiters (no Prony fit, no shared code)

`kato02_exact_chain.py M tr tend [n] [dc] [mode]` integrates Kato's
hereditary formulation exactly by converting each image order m into an
Erlang memory chain (S_m' = (S_{m-1}-S_m)/tr, exact since
A_m' = (A_{m-1}-A_m)/tr), with the same friction and regularization,
its own adaptive RK3(2).  Kernel modes: `exact`, `kato3` (his Eq. 3
misprint), `static` (images permanently on), `current` (history-free
shortcut), `single`/`singleB` (one image family only).  An optional
probe file records v(t) at 15 and 18 km (his Fig. 4 depths).

`mn19_exact_chain.py tc tend [ds] [RcR] [lenR]` integrates the uniform
Maxwell kernel exactly via one memory state per cell
(K(t) = K0 e^(-t/tc); tc <= 0 gives the elastic reference).

## Results and what the gates check

1. Kato ELASTIC: Tcy = 137.87 yr and us = 3.88 m with ZERO free
   parameters against his 138 yr and 3.8 m; gate on Tcy.  This
   validates the digitized profile, the friction, and the loading.

2. Kato VISCOELASTIC: his Table 1 (Tcy = 138 for tr = 1..300, i.e. no
   viscoelastic effect at all) does NOT replicate.  The full-resolution
   sweep (virgin start, cycles 2-6) gives

       tr [yr]     1      3     10     30    100    300   h15,tr10 h12,tr10
       Tcy [yr]  105.5  106.9  110.3  116.9  125.6  131.9   94.6     93.6
       change    -23%   -22%   -20%   -15%    -9%    -4%    -31%     -32%
       us [m]     3.89   3.89   3.89   3.88   3.88   3.88  (elastic 3.88)

   The gate cross-checks rsf_solve against the exact chain on a reduced
   problem (tr = 1: both about 0.79 x elastic, within 8 percent of each
   other; currently 1.4 percent).

3. The explanation for his null: integrating his own equations with
   only image family A (depths 2mh + z' and mirror) while dropping
   family B (depths 2mh - z', the near-field substrate reflections
   that hug the plate bottom below the afterslip zone) gives cycles
   IDENTICAL to elastic at every tr (gate: within 3 percent), with us
   unchanged, and postseismic fault slip rates that reproduce his
   Fig. 4 exactly in character (tr = 3 and tr = 30 curves nearly
   coincident, the shorter tr peak "slightly larger"; the exact kernel
   instead separates them clearly).  Regenerate that comparison with
   the chain's probe output, e.g.
   `kato02_exact_chain.py 10 3 165 144 0.07 exact probe.dat` for
   elastic/exact/single at tr = 3 and 30, plotting v(15 km) and
   v(18 km) for 20 yr after the second event.  This
   is the SAME single-family construction this project found and fixed
   in the repo's own layered antiplane kernels: family A alone is
   dynamically inert on the fault while still producing order-unity
   viscoelastic SURFACE deformation, which is why his Figs. 5-7 (slip
   histories post-processed through medium Green functions) look
   sensible.  Every other candidate explanation was tested and fails:
   his m <= 10 truncation (1-3 percent), his Eq. 3 misprint (-8 to
   -14 percent, tr-dependent), a permanently relaxed kernel
   (-14 percent, tr-independent but not elastic), the history-free
   shortcut (-5 to -19 percent), initial conditions and protocol
   (shortening appears in the FIRST cycle from his virgin start).

4. M&N ELASTIC patch (Rc_EL/R = 0.35): t_rec = 12.4 yr, rsf_solve and
   chain agreeing to 3 digits (their FD periodic value is 7.15 yr at
   Rc/R = 0.4; QD plus geometry offset expected).

5. M&N VISCOELASTIC (t_c/t_load = 1.0): t_rec = 28.5 yr (2.2x
   elastic), rsf_solve -ve_mode 1 vs chain within 3 percent.  The full
   scan reproduces their Fig. 5 phenomenology: t_rec rises and
   DIVERGES as t_c decreases, fit by alpha/(x - t_cr) + beta with
   t_cr(QD, strip) = 0.86, vs their 1.518 (FD, periodic).  Approximate
   periodicity (4 patches every 4.096R) raises t_cr to about 1.1; the
   remaining gap is inertia (fully dynamic events deepen the
   post-event stress hole and starve the patch earlier), exactly the
   QD caveat they state.  QD also predicts a mild t_rec MINIMUM
   (about 0.75x elastic near t_c/t_load = 2.4) before the divergence,
   which their FD scan does not show; whether that is QD-only or
   outside their scan is an open question.

6. M&N STUCK class: t_c/t_load = 0.8 produces exactly one event and
   locks permanently (their ST), confirmed by the chain.

## The two machinery findings these tests forced (fixed in this patch set)

- The VE sink's within-step forcing was lagged to the last accepted
  step (frozen sink), an error OUTSIDE the TS controller; harmless for
  the validated layered configurations (1-2 percent) but biasing
  near-critical uniform-Maxwell recurrence by tens of percent (the
  first M&N scan produced a spurious "runaway repeater" where the
  converged answer is stable lengthened cycles).  Fixed by
  stage-consistent sink forcing (-ve_h_stage 1, the default; 0
  restores the cheap lagged mode, np fewer matvecs per RHS).
- Restart chaining re-wrote a checkpoint AT the restart step itself
  with a not-yet-valid time step (poisoning the file with dt = NaN);
  fixed by suppressing the write at the restored step.

## Caveats to keep in mind

- Through-plate faults (Kato's geometry) have a stress-free uniform
  slip mode in the relaxed limit: no steady cycle attractor exists.
  Protocol-matched EARLY cycles are meaningful (and chain-validated);
  multi-kyr trajectories integrate any small kernel-representation
  bias along the neutral mode (the exact chain drifts back toward the
  elastic recurrence over ~100 cycles, the 6-pole Prony ladder, which
  cannot carry multi-century Erlang tails, drifts differently).  For
  physics keep the fault bottom above the plate bottom; those
  configurations converge and are initial-condition independent.
- rsf_solve is quasi-dynamic; near-critical M&N-type problems are
  sensitive to inertia (t_cr shifts).
