# `rsf_solve -fd_mode`: fully dynamic earthquake sequences

`rsf_solve` integrates rate-and-state earthquake sequences quasi-dynamically
(QD): the elastostatic interaction operators plus radiation damping. With
`-fd_mode 1` the coseismic phases are instead computed fully dynamically
(FD), with inertia and wave propagation, on the same fault geometry, friction
and loading, and the quasi-dynamic integrator resumes afterwards. Runs
without the option are unchanged, bit for bit.

This document gives the formulation, what it was checked against, the
options, and the current limits. Development notes and the Python prototype
that the C code was verified against are in `fully_dynamic/` (status summary
in `fully_dynamic/fd_plan.md`, section 0; how to resume the work in
`fully_dynamic/HANDOFF.md`).

Capabilities at a glance (2026-10-25):

- 2-D (in-plane, antiplane): one fault or parallel segments with
  normal-stress coupling; rate-and-state and slip-weakening friction;
  second-order time marching (default) with variable substeps and an exact
  coarse history; mirror free surface for antiplane. Benchmarked against
  exact solutions, TPV205-2D, Harris-Day and SEAS BP1.
- 3-D, two parallel planes (stepovers, strike slip): one lattice for both, whole
  space or half space by the mirror image, exact statics (Is, and In with
  -calc_sigma_dot) for both planes; check with
  src/testing/rsf_solve/test_fd3d_pair.sh. The planes must be farther apart than
  the P-wave reach of one step (about 0.9 element at the default CFL).
- 3-D (one planar fault of squares, strike or dip slip): full space; half
  space with a mirror image for vertical strike-slip faults (BP5-FD first
  event within 0.1 s and 2 percent of BICyclE), an optional opening layer
  that makes the mirror exact, and a full virtual free surface for faults of
  any dip with normal-stress coupling (validated for buried faults; a
  resolution-dependent error at the trace for faults that reach the
  surface). Fixed fine steps; quasi-steady hand-back for speed.

## 1. Formulation

Boundary integral form of the elastodynamic traction on the fault, with the
static operators of `rsf_solve` kept and a hereditary dynamic correction
added. For element i at time t, with slip rate history V,

    tau_i(t) = tau0_i + loading_i t + sum_j Ks[i,j] delta_j(t)
               + sum_{l=1}^{L} sum_j W[l][i,j] V_j(t - l dt) - eta V_i(t),

eta = G/(2 c_s) the radiation damping that QD uses alone, Ks the static
operator (`rsf_solve`'s `Is`), and W[l] = K[l] - Ks dt the difference between
the elastodynamic response at lag l dt to unit slip rate during one time step
and its static limit. W decays to zero behind the waves, which is what makes
the window L finite; before the first arrival it equals -Ks dt exactly.
Radiation damping is therefore not an approximation in this scheme but the
instantaneous part of the kernel. The normal traction has the same structure
(Kn, Wn) and feeds the effective normal stress when `-calc_sigma_dot` is on.

Kernels. 2-D plane strain (in-plane, mode II, element type 30) and antiplane
(mode III, type 33), homogeneous medium, from the closed-form time-domain
primitives of Tada and Madariaga (2001) and Ando et al. (2007), with a
large-time series beyond beta t / r = 10. The antiplane free surface (half
plane, y = 0) is added by image sources. The tables are stored per pair of
parallel segments on the element-offset (Toeplitz) axis and applied by FFT.

Time stepping during a phase. Fixed fine step dt_min = CFL dx / c_s (CFL 0.5)
with end-of-step collocation; the current step's slip rate enters through the
first-lag operator K[1], which couples the neighbours reached by the P wave
within one step, and that coupled system is solved exactly (banded) together
with the friction law: a Newton iteration in log V for rate and state. The
step is lengthened to k dt_min (k a power of two up to `-fd_kmax`) whenever
xi D_c / V_max allows, using the cumulative kernels sum_{l<=k} K[l]; the
history stays on the fine grid, so no re-interpolation is needed. The state
variable is advanced with the exact solution of the aging law over the step
(Lapusta et al. 2000); an explicit step was found to produce early healing and
40 percent too little shallow slip in BP1.

Switching. A phase starts when the maximum slip rate crosses `-fd_vsw`
(1e-3 m/s), detected as a TS event; the slip-rate history before the switch is
rebuilt from the last accepted QD steps. The phase ends when the slip rate has
stayed below `-fd_vback` (1e-4 m/s) for one full window, i.e. when the
hereditary term has flushed; psi, tau, sigma and slip are then written back
into the QD solution and the TS continues from the end of the phase. The
result is insensitive to `-fd_vsw` between 1e-3 and 1e-2 m/s (1 percent in
slip).

Window. Lags are kept while |W| exceeds `-fd_tol` (1e-2) times the static
interaction envelope; tightening to 1e-3 changed BP1-FD results by less than
the printed precision. `-fd_lmax` caps the window (memory), which is
appropriate when the natural window is set by distant pairs.

## 2. Verification

The C implementation reproduces the Python prototype to roundoff (TPV205-2D
station series to 1e-10 relative; BP1-FD and the single-fault sequence to the
printed precision of steps, durations, peak slip rates and slip per event),
and the prototype was checked against:

- exact solutions: the self-similar (Kostrov/Broberg) crack in modes II and
  III to 1e-3; static limits and the plane-wave identity;
- SCEC TPV205-2D against E. Dunham's MDSBI reference at 6.25 m: main-front
  arrivals within about 1.5 time steps (first order in time) at 100 and
  50 m, final slip within 0.4 percent;
- Harris and Day (1993) parallel stepovers, case D reproduced (trigger
  delays and locations); supershear transition lengths of Andrews (1976b)
  and Dunham (2007) when the cohesive zone is resolved by about 40 elements;
- SEAS BP1 in QD and FD at 50 and 25 m: event times, recurrence (105.7 and
  118.7 yr in FD), slip per event and station histories against the
  community solutions (sem2dpack, Jiang's SBIM); peak slip rates within a
  few percent at the shallow stations at 25 m. `rsf_solve -fd_mode 1` at
  25 m gave FD events at 196.36, 302.10 and 420.80 yr, the prototype's
  196.36, 302.11 and 420.81 yr, with the same step counts and peak slip
  rates (both before the history alignment fix of 2026-10-20, which moves
  the 50 m event times by 0.1 yr; fully_dynamic/fd_c_design.md section 13).

`rsf_solve`'s own QD on the BP1 half-plane geometry agrees with the prototype
and with group 4's reference to 0.01 yr, and the FD tables' static limit
agrees with the assembled `Is`/`In` to 1e-16 (checked at every start,
`static gate` lines).

Added since (details and numbers in `fully_dynamic/fd_plan.md`):

- second-order time marching (Noda, Sato and Kurihara 2020): equal to the
  prototype's version to 1e-10; BP1-FD converges 3 to 5 times faster than
  with the first-order scheme; BP5-FD at 1 km matches Lambert's BICyclE
  arrivals to 0.1 s and peak slip rates to 2 percent;
- 3-D kernels (Tada 2006) against an independent point-source quadrature
  to 1e-8 or better, C against Python at 600 points; TPV3 and TPV101 runs
  (TPV101 converging at first order, no station data available yet);
- free surface: mirror plus opening layer against Okada's half space to
  4.5e-4; the full virtual surface for buried dipping faults against Okada
  (1e-3 to 1e-4) and against the independent 2-D plane-strain prototype
  (slip to 0.3 to 0.6 percent, identical rupture times, with normal-stress
  coupling; `fd3d_test q2d`); BP4-QD and BP5 first events against the
  published runs.

## 3. Usage

    rsf_solve <usual options> -fd_mode 1 [-fd_lmax 6000] [-fd_stations 49,149,249]

Options (defaults in parentheses):

- `-fd_mode <0|1|-1>` off (0); dynamic phases (1); -1 only reports the switch
  threshold crossings.
- `-fd_vsw <m/s>` (1e-3) switch to FD; `-fd_vback <m/s>` (1e-4) return to QD
  after one window below it.
- `-fd_tol <f>` (1e-2) window tolerance; `-fd_lmax <n>` (0 = natural bound)
  cap on the window in fine-grid lags; `-fd_cfl <f>` (0.5) c_s dt_min / dx.
- `-fd_xi <f>` (0.25) and `-fd_kmax <n>` (64): variable step dt <= xi D_c /
  V_max, in powers of two of dt_min.
- `-fd_et <f>` (1.0) collocation offset (1 end of step; about 2 is the
  Peirce-Siebrits scheme, halves grid-scale noise at a 0.5 percent cost in
  accuracy).
- `-fd_gate_tol <f>` (0.05) tolerance of the static gate (abort at 5x).
- `-fd_length_scale <m>` (1) geometry unit in metres.
- `-fd_max_steps <n>` (2e6) abort a phase beyond this; `-fd_dt_restart <s>`
  (1e-2) first QD step after a phase.
- `-fd_conv <0|1|2>` (2) history convolution: direct sum; partitioned FFT
  (about twice the memory); hybrid, partitioned for short substeps and direct
  for long ones (about three times the memory, fastest). Results agree to
  roundoff. `-fd_conv_block <n>` (0 = automatic) partition length in lags.
- `-fd_coarse <0|1>` (1) steps of `-fd_kmax` substeps start at multiples of
  `-fd_kmax`; once a window holds only such steps, the history term is summed
  over the aligned blocks with block-aggregated kernels, which is exact and
  `-fd_kmax` times cheaper. The alignment rule changes BP1 event times by
  about 0.003 yr and peak slip rates by about 1 percent at 50 m.
- `-fd_back_eps <f>` (0) return to QD once every slip rate has changed by
  less than this fraction over one window and the maximum is below
  `-fd_back_vmax` (half of `-fd_vsw`), instead of waiting until all slip
  rates are below `-fd_vback`. Justified because the hereditary term of a
  slip rate that is constant over a window equals the QD radiation damping;
  0.05 shortened BP1 phases by about 40 percent and moved event times by
  about 1e-5 relative.
- `-fd_stations i,j,...` global patch indices written during every phase to
  `fd_phase_<n>_st<i>.dat` (time since switch, slip, slip rate, shear
  stress, fault-wide max slip rate, substep multiple); the final slip of each
  phase goes to `fd_phase_<n>_slip.dat`.

Output: one line per phase (`rsf_fd: phase n: t, steps, duration, peak slip
rate, mean slip, band, wall time`); the usual `rsf_events.dat`, monitor,
stations and catalog continue across phases (the tracker's onset line is the
pre-phase state, its arrest line the post-phase state).

Memory and cost. Tables take 16 L M bytes per table (two for in-plane), L
the window in lags and M the padded FFT length (next power of two above 2 N -
1 for one segment): BP1 at 50 m needs 0.5 GB with `-fd_lmax 6000`, 25 m
about 2.3 GB with `-fd_lmax 12000` with `-fd_conv 0`, about three times
that with the default hybrid convolution (use `-fd_conv 0` if memory is
short). A BP1 phase takes under a minute at 50 m on one core with the default
settings (it took 5 min before the speed-ups of 2026-10-20 and -21), and
about 45 min at 25 m before the coarse history; the history convolution dominates and scales as L M per
step, i.e. fourfold per halving of the element size; resolution of 5 or more elements per cohesive zone is needed for
slip and timing, about 40 for supershear transitions or near-threshold
triggering.

3-D (2026-10-23): with all elements 3-D rectangles forming a complete grid
of equal squares on one plane, `-fd_mode 1` runs the phases with the 3-D
marcher (Tada 2006 quadrant kernels). Whole space (`-full_space 1`, slip
along strike or dip) or half space for a vertical strike-slip fault whose
top row reaches the surface, where the free surface enters the dynamic part
through a mirror image and the static part uses interact's exact half-space
operator. 3-D phases use fixed fine steps by default; `-fd_back_eps 0.05` (quasi-steady
hand-back, as in 2-D) shortens the afterslip tail, about 4 times faster on
the BP5-FD first event at 2 km with the later slip unchanged to 1e-4 m.
`-fd_layer n` (3-D half space, vertical strike-slip fault reaching the
surface) completes the mirror image with a strip of opening elements on the
surface, +-n fault elements across the trace (`-fd_layer_m`, default 2:
elements of dx/2; `-fd_layer_margin`, default n, beyond the fault ends).
It removes the normal-traction error of the mirror; for BP5 the effect on
the first event is below 0.5 percent (fd_plan.md). Off by default.
In a half space, a planar fault that is not a vertical strike-slip fault (or any
fault with `-fd_surface 1`) uses a full virtual free surface: surface elements of
dx/4 (`-fd_surface_m`) over the fault's surface projection plus
`-fd_surface_ext` fault elements (default the fault width), three unknowns each,
zero traction. With `-calc_sigma_dot` the fault's normal stress follows the exact
static In plus the surface's dynamic part. This is new and tested on small
problems only (see fd_plan.md, steps B2 and B3). The surface's shear modes are
unstable for cS dt/hs between about 0.85 and 1.05; by default the surface is refined
to keep cS dt/hs >= 2 (dx/4 at the default CFL). Rows are graded by default
(`-fd_surface_grade 1`): widths of at most half the distance to the fault, from
dx/m near the trace up to dx far away. A fault whose top edge lies in the surface
continues one row above it (`-fd_surface_extend`, tied to the top row). Check with
src/testing/rsf_solve/test_fd3d_dip.sh. For faults that reach the surface the full
surface converges slowly at the trace (fd_plan.md, step B5); for vertical strike-slip
faults the mirror image is exact there and should be used. `-fd_table_statics 1`
replaces interact's exact Is by the tables' own statics (experiments only).
For 3-D sequences, `-fd_back_dyn 1e4` ends a dynamic phase once the fully dynamic
traction differs from the quasi-dynamic one by less than 10 kPa everywhere (for
`-fd_back_nwin` steps, default 10) and the peak rate is below the switching
threshold; on BP5-FD this halves the phase with afterslip unchanged to 1e-4 m
(fd_plan.md, L2). In 3-D the defaults are now `-fd_back_eps 0.05 -fd_back_dyn 1e4
-fd_reentry 10` (after a phase, the next one opens only above 10 vsw until the peak
rate has fallen below vback; fd_plan.md, L3); 2-D defaults are unchanged. The
traction criterion applies up to the re-entry threshold, so slow nucleation above vsw
is carried by QD until the rate exceeds 10 vsw (fd_plan.md, P4). For 3-D models at
500 m or finer, `-ts_rk_type 3bs` was needed for the QD part (the default scheme can
stop in a rejection storm during QD coseismic phases).
`-fd_kmax3 1` enables steps of several substeps in 3-D (experimental, see
fd_plan.md). The 3-D phase station files have a fifth column, the number
of substeps of the step.

`-fd_nl 1` (default since 2026-10-25, 2-D and 3-D) is second-order time
marching after Noda, Sato and Kurihara (2020): midpoint slip rates in the
history and the radiation damping applied to the slip rate at the
collocation time. It is more accurate than the earlier first-order scheme
(`-fd_nl 0`) for event times, peak slip rates and rupture arrivals in every
benchmark tried (BP1-FD, BP5-FD, TPV205, a rate-and-state rupture); near
the jump threshold of a stepover, triggering delays at 100 m may be too
short and should be checked at finer resolution. See fd_plan.md.

Restrictions in this version, refused with a message at start: one MPI
rank; `-ve_mode 0`; aging law (`-state_law 1`); 2-D elements only, forming
straight runs of equal, collinear elements with one common element size,
all parallel; full plane for in-plane, full or half plane (vertical
segments) for antiplane. H-matrix backends for `Is` are fine (the static
part goes through `MatMult`); the dynamic tables are dense per segment pair.

## 4. Example

Ready-made drivers in `fully_dynamic/c_port/` (run from the repository
root): `run_bp1_fd.sh [dz] [tend] [lmax]` writes the BP1 files, runs QD and
FD, prints the event times and phase lines; `run_step_fd.sh W overlap margin
[tend] [qd|fd|both]` runs a two-segment stepover sequence and
`summarize_rsf_step.py` tabulates the jumps (delay of the second segment's
rupture after the first, binned as immediate, minutes, hours, days).

### Stepover sequence

    fully_dynamic/c_port/run_step_fd.sh 1.0 2.0 0.5 1400 both
    python3 fully_dynamic/c_port/summarize_rsf_step.py step_*_W1.0_ov2.0_vs0.5

Two parallel 24 km in-plane segments (100 m elements, velocity-weakening
except for the margins), step width W in km (positive W puts the second
segment at +y, which is the dilational side for interact's left-lateral
positive slip), overlap and margin in km, end time in years, mode. With 1
km width and 0.5 km margins the quasi-dynamic run gives the second segment
3 h after the first in every cycle; the fully dynamic run gives joint
ruptures and 25 to 50 min delays in alternating cycles (fd_plan.md has the
full table). Two further arguments, the length of the second segment in
km and b - a on its weakening part, make the segments unequal, so that the
receiver's position in its own cycle varies from event to event; the
summary then reports the fraction of source events followed by a receiver
rupture within one day and within 30 days:

    fully_dynamic/c_port/run_step_fd.sh 1.0 2.0 0.5 2500 both 16 0.004

### BP1 by hand

SEAS BP1 (antiplane half plane, 40 km, 50 m):

    python3 - << 'EOF'
    # geom_bp1.in / rsf_bp1.dat / ic_bp1.in as in fully_dynamic/bp1.py
    EOF
    rsf_solve -geom_file geom_bp1.in -rsf_file rsf_bp1.dat -rsf_ic_file ic_bp1.in \
       -dc 0.008 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.204e10 \
       -s_wave_speed 3464 -rtol 1e-6 -stop_time_yr 520 -track_events 1 -vel_event 1e-3 \
       -fd_mode 1 -fd_lmax 6000 -fd_stations 49,149,249

gives events at 196.36, 301.69 and 420.48 yr (QD: 196.36, 271.96, 349.77,
...), with per-event mean slips 1.561, 1.303 and 1.461 m.

## 5. References

Ando, R., Kame, N., Yamashita, T. (2007), Earth Planets Space 59, 173.
Andrews, D. J. (1976b), J. Geophys. Res. 81, 5679.
Dunham, E. M. (2007), J. Geophys. Res. 112, B07302.
Erickson, B. A. et al. (2020, 2023), Seismol. Res. Lett. 91, 874; Bull.
Seismol. Soc. Am. 113, 499.
Harris, R. A., Day, S. M. (1993), J. Geophys. Res. 98, 4461.
Lapusta, N., Rice, J. R., Ben-Zion, Y., Zheng, G. (2000), J. Geophys. Res.
105, 23765.
Tada, T., Madariaga, R. (2001), Int. J. Numer. Meth. Engng. 50, 227.
