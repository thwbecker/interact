#!/usr/bin/env python3
"""psgrn_convergence_compare.py: table and figure for run_psgrn_convergence.sh.

usage: psgrn_convergence_compare.py [dir] [variants] [ve_thrust_relax.py]

For every variant (and its _fine PSCMP patch variant where present) it
reports the basin depth of the uplift change at 10 and 90 tM, its
difference to the 2-D code in percent, the rms difference over the
profile, and the zero crossings; then plots PSCMP minus 2-D for all
variants so that the parameter responsible for the offset stands out.
Variants ending in _g1 are compared with the gravity-on 2-D run, all
others with gravity off.  Writes psgrn_convergence.png.
"""
import sys, os, subprocess
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else '.'
variants = (sys.argv[2] if len(sys.argv) > 2 else 'base acc01 nr252 nz29 r1200 all acc005 nr504 all_g1').split()
script = sys.argv[3] if len(sys.argv) > 3 else 've_thrust_relax.py'
H = 30.0; KM = 111.195
lat = np.linspace(90.0 / KM, -160.0 / KM, 101)
xH = -lat * KM / H
times = (10.0, 90.0)


def two_d(g):
    out = f'ref2d_g{g}'
    if not os.path.isfile(out + '.npz'):
        subprocess.run([sys.executable, script, '--H', '30', '--mu1', '30', '--nu', '0.25',
                        '--rho1', '3300', '--rho2', '3800', '--g', '9.8', '--gravity', str(g),
                        '--times', '10', '90', '--xmax', '100', '--nx', '1601', '--out', out,
                        '--no-plot'], check=True, stdout=subprocess.DEVNULL)
    r = np.load(out + '.npz'); x = r['x']
    return {tt: np.interp(xH, x, (r['uz'][:, it + 1] - r['uz'][:, 0]) * 100) for it, tt in enumerate(times)}


def read(out):
    U = np.loadtxt(os.path.join(d, out, 'U_down.dat'), skiprows=1)
    t = U[:, 0] / 365.25; u = U[:, 1:102]
    return {tt: -(np.array([np.interp(tt, t, u[:, i]) for i in range(101)]) - u[0]) * 100 for tt in times}


def zeros(y):
    j = np.where(np.diff(np.sign(y)) != 0)[0]
    return [round(float(xH[i] - y[i] * (xH[i + 1] - xH[i]) / (y[i + 1] - y[i])), 2) for i in j]


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
print(f'{"variant":12s} {"t":>3s} {"basin":>7s} {"2-D":>7s} {"diff%":>6s} {"rms":>6s}  zero crossings (PSCMP | 2-D)')
rows = []
for v in variants:
    g = 1 if v.endswith('_g1') else 0
    ref = two_d(g)
    for out in (f'out_{v}', f'out_{v}_fine'):
        if not os.path.isfile(os.path.join(d, out, 'U_down.dat')):
            continue
        P = read(out)
        for it, tt in enumerate(times):
            y, r = P[tt], ref[tt]
            diff = 100 * (y.min() / r.min() - 1)
            rms = np.sqrt(np.mean((y - r)**2))
            print(f'{out[4:]:12s} {tt:3.0f} {y.min():7.1f} {r.min():7.1f} {diff:+6.1f} {rms:6.2f}  {zeros(y)} | {zeros(r)}')
            ax[it].plot(xH, y - r, label=out[4:])
for it, tt in enumerate(times):
    ax[it].axhline(0, color='0.6', lw=0.6); ax[it].grid(alpha=0.3)
    ax[it].set_title(f'PSCMP minus 2-D, uplift change, t = {tt:.0f} tM', fontsize=10)
    ax[it].set_xlabel('x / H (towards hanging wall)'); ax[it].set_ylabel('difference x 100 / slip')
    ax[it].legend(fontsize=8)
fig.suptitle('PSGRN/PSCMP numerical parameters, L = 2000 km thrust, H = 30 km', fontsize=11)
fig.tight_layout(); fig.savefig('psgrn_convergence.png', dpi=130)
print('wrote psgrn_convergence.png')
