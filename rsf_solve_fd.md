# `rsf_solve -fd_mode`: fully dynamic earthquake sequences

`rsf_solve` integrates rate-and-state earthquake sequences quasi-dynamically
(QD): the elastostatic interaction operators plus radiation damping. With
`-fd_mode 1` the coseismic phases are instead computed fully dynamically
(FD), with inertia and wave propagation, on the same fault geometry, friction
and loading, and the quasi-dynamic integrator resumes afterwards. Runs
without the option are unchanged, bit for bit.

This document gives the formulation, what it was checked against, the
options, and the current limits. Development notes and the Python prototype
that the C code was verified against are in `fully_dynamic/`.

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
short). A BP1 phase takes about 3 min at 50 m on one core with the default
settings (2 min with `-fd_back_eps 0.05`), and about an hour at 25 m before
the M6 speed-ups; the history convolution dominates and scales as L M per
step, i.e. fourfold per halving of the element size; resolution of 5 or more elements per cohesive zone is needed for
slip and timing, about 40 for supershear transitions or near-threshold
triggering.

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
