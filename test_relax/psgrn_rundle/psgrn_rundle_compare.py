#!/usr/bin/env python3
"""psgrn_rundle_compare.py: read the PSCMP output of run_psgrn_rundle.sh
and compare with Rundle (1982) Figs. 2 and 3 and with the 2-D codes.

usage: psgrn_rundle_compare.py [dir] [lengths] [gravs] [ve_thrust_relax.py]
       dir       where out_g<g>_L<L>/ sit                      (.)
       lengths   "200 600 2000"
       gravs     "0 1"
       script    path to ve_thrust_relax.py; if found it is run with the
                 matching parameters (H 30, mu 30 GPa, rho 3300/3800,
                 g 9.8, times 10 90) to provide the 2-D curves

Produces psgrn_rundle.png (uplift change at 10 and 90 tM, gravity off and
on, every fault length, 2-D curve, Rundle's digitised points) and
psgrn_rundle_horizontal.png (same for the horizontal component towards
the hanging wall), and prints: basin depth and zero crossings per
length, finite-length factors relative to the longest fault, rms against
the 2-D curve for the longest fault, rms against Rundle's points for
L = 200 km without gravity, and the coseismic field of the longest
fault against the 2-D code (a geometry and sign check of the PSCMP
setup).  Conventions: x positive towards the hanging wall (south in the
PSCMP frame) in units of H = 30 km; displacements x 100 per metre of
slip; uplift positive up; horizontal positive towards the hanging wall.
"""
import sys, os, subprocess
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else '.'
lengths = [int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else '200 600 2000').split()]
gravs = [int(v) for v in (sys.argv[3] if len(sys.argv) > 3 else '0 1').split()]
script = sys.argv[4] if len(sys.argv) > 4 else 've_thrust_relax.py'
H = 30.0
KM = 111.195
lat = np.linspace(90.0 / KM, -160.0 / KM, 101)
xH = -lat * KM / H                       # towards the hanging wall, in H
times = (10.0, 90.0)

dig = {10.0: np.array([[-3.0, 5], [-2.0, 4], [-1.5, 2], [-0.85, 0], [-0.5, -6], [0.0, -15],
                       [0.5, -24], [0.65, -26], [1.0, -23], [1.5, -12], [2.0, -4], [2.2, 0],
                       [3.0, 5], [4.0, 6], [5.0, 5]], float),
       90.0: np.array([[-3.0, 5], [-2.5, 4], [-2.0, 0], [-1.5, -13], [-1.0, -30], [-0.5, -57],
                       [0.0, -70], [0.5, -75], [1.0, -70], [1.5, -52], [2.0, -30], [2.5, -13],
                       [3.0, -2], [3.5, 5], [4.0, 9], [4.5, 11], [5.0, 13]], float)}
rundle_fig3_basin = -45.0                 # with gravity, 45 tau_a


def read(g, L):
    out = os.path.join(d, f'out_g{g}_L{L}')
    U = {c: np.loadtxt(os.path.join(out, f'U_{c}.dat'), skiprows=1) for c in ('down', 'north')}
    t = U['down'][:, 0] / 365.25
    res = {'t': t}
    for c in ('down', 'north'):
        u = U[c][:, 1:102]
        res[c + '0'] = u[0]
        for tt in times:
            res[(c, tt)] = np.array([np.interp(tt, t, u[:, i]) for i in range(101)]) - u[0]
    # uplift (up) and horizontal towards hanging wall (south = -north)
    res['uz0'] = -res['down0'] * 100; res['ux0'] = -res['north0'] * 100
    for tt in times:
        res[('uz', tt)] = -res[('down', tt)] * 100
        res[('ux', tt)] = -res[('north', tt)] * 100
    return res


def two_d(g):
    """2-D curves from ve_thrust_relax.py (uplift and horizontal change)."""
    if not os.path.isfile(script):
        return None
    out = f'ref2d_g{g}'
    if not os.path.isfile(out + '.npz'):
        subprocess.run([sys.executable, script, '--H', '30', '--mu1', '30', '--nu', '0.25',
                        '--rho1', '3300', '--rho2', '3800', '--g', '9.8', '--gravity', str(g),
                        '--times', '10', '90', '--xmax', '100', '--nx', '1601', '--out', out,
                        '--no-plot'], check=True, stdout=subprocess.DEVNULL)
    r = np.load(out + '.npz')
    x = r['x']
    o = {'x': x, 'uz0': np.interp(xH, x, r['uz'][:, 0] * 100), 'ux0': np.interp(xH, x, r['ux'][:, 0] * 100)}
    for it, tt in enumerate(times):
        o[('uz', tt)] = np.interp(xH, x, (r['uz'][:, it + 1] - r['uz'][:, 0]) * 100)
        o[('ux', tt)] = np.interp(xH, x, (r['ux'][:, it + 1] - r['ux'][:, 0]) * 100)
    return o


def zeros(y):
    j = np.where(np.diff(np.sign(y)) != 0)[0]
    return [round(float(xH[i] - y[i] * (xH[i + 1] - xH[i]) / (y[i + 1] - y[i])), 2) for i in j]


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Lmax = max(lengths)
for comp, fname in (('uz', 'psgrn_rundle.png'), ('ux', 'psgrn_rundle_horizontal.png')):
    fig, ax = plt.subplots(len(gravs), 2, figsize=(13, 4.8 * len(gravs)), squeeze=False)
    for ig, g in enumerate(gravs):
        ref = two_d(g)
        R = {L: read(g, L) for L in lengths}
        if comp == 'uz':
            print(f'\n=== gravity {"on" if g else "off"}')
            if ref is not None:
                # the receiver on the trace (x = 0) sits on the coseismic
                # step (45 vertical, 85 horizontal x 1e-2 slip) and the two
                # codes resolve it differently; exclude |x| < 0.1 H
                off = np.abs(xH) > 0.1
                e0 = (R[Lmax]['uz0'] - ref['uz0'])[off]; ex = (R[Lmax]['ux0'] - ref['ux0'])[off]
                print(f'coseismic, L = {Lmax} km vs 2-D (|x| > 0.1 H): uz max |d| {np.abs(e0).max():.2f} '
                      f'rms {np.sqrt(np.mean(e0**2)):.2f}, ux max |d| {np.abs(ex).max():.2f} '
                      f'rms {np.sqrt(np.mean(ex**2)):.2f} (x100 slip; 2-D peaks uz {np.abs(ref["uz0"][off]).max():.1f}, '
                      f'ux {np.abs(ref["ux0"][off]).max():.1f})')
        for it, tt in enumerate(times):
            a = ax[ig, it]
            for L in lengths:
                y = R[L][(comp, tt)]
                a.plot(xH, y, '-', lw=1.3 if L == Lmax else 1.0, label=f'PSCMP L = {L} km')
                if comp == 'uz':
                    s = f'  t = {tt:3.0f} tM  L = {L:5d}: basin {y.min():7.1f} at x = {xH[y.argmin()]:5.2f} H, zeros {zeros(y)}'
                    if L != Lmax:
                        s += f', depth ratio to L = {Lmax}: {y.min() / R[Lmax][(comp, tt)].min():.3f}'
                    if ref is not None and L == Lmax:
                        s += f', rms vs 2-D {np.sqrt(np.mean((y - ref[(comp, tt)])**2)):.2f}'
                    if g == 0 and L == 200:
                        p = np.interp(dig[tt][:, 0], xH, y)
                        s += f', rms vs Rundle {np.sqrt(np.mean((p - dig[tt][:, 1])**2)):.1f} (peak {abs(dig[tt][:, 1]).max():.0f})'
                    print(s)
            if ref is not None:
                a.plot(xH, ref[(comp, tt)], 'k--', lw=1.2, label='2-D (ve_thrust_relax)')
            if comp == 'uz' and g == 0:
                a.plot(dig[tt][:, 0], dig[tt][:, 1], 'o', mfc='none', color='r', label='Rundle 1982 Fig. 2')
            if comp == 'uz' and g == 1 and tt == 90.0:
                a.axhline(rundle_fig3_basin, color='r', ls=':', label='Rundle Fig. 3 basin')
            a.axhline(0, color='0.7', lw=0.6); a.grid(alpha=0.3)
            a.set_title(f'gravity {"on" if g else "off"}, t = {tt:.0f} tM ({tt / 2:.0f} tau_a)', fontsize=10)
            a.set_xlabel('x / H (towards hanging wall)')
            a.set_ylabel(('uplift' if comp == 'uz' else 'horizontal') + ' change x 100 / slip')
            a.legend(fontsize=8)
    fig.suptitle('PSGRN/PSCMP, 30 deg surface-breaking thrust to 0.5 H, H = 30 km, mu = lam = 30 GPa', fontsize=11)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    print('wrote', fname)
