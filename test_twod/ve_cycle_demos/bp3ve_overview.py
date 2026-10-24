#!/usr/bin/env python3
"""bp3ve_overview.py: event size vs time, elastic and end-member Maxwell
times, thrust and normal (BP3 dip 60, 25 m, sigma coupled, gravity on).
usage: bp3ve_overview.py [root] [out]"""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
root = sys.argv[1] if len(sys.argv) > 1 else '.'
out = sys.argv[2] if len(sys.argv) > 2 else 'bp3ve_overview'
def cat(c, tag):
    d = np.loadtxt(f'{root}/{c}/cat_{tag}.dat', comments='#', ndmin=2); return d[:, 1], d[:, 6]
rows = [('ve_g1', 'el', 'thrust, elastic', '4-cycle {63, 87, 91, 98}'),
        ('ve_g1', 'tm15', 'thrust, t_M = 15 yr (0.2 T)', '1-cycle 88'),
        ('ve_g1', 'tm100', 'thrust, t_M = 100 yr (1.2 T)', '3-cycle {65, 86, 90} from 2600 yr'),
        ('ve_g1', 'tm1000', 'thrust, t_M = 1000 yr (12 T)', 'same 3-cycle, from 7700 yr'),
        ('ve_g1_normal', 'el', 'normal, elastic', '1-cycle 95'),
        ('ve_g1_normal', 'tm50', 'normal, t_M = 50 yr (0.5 T)', '2-cycle {75, 99} from 950 yr'),
        ('ve_g1_normal', 'tm500', 'normal, t_M = 500 yr (5 T)', 'same 2-cycle, from 10200 yr')]
fig, ax = plt.subplots(len(rows), 2, figsize=(14, 1.55 * len(rows)), sharey='row',
                       gridspec_kw=dict(width_ratios=[3, 1], wspace=0.04, hspace=0.35))
for i, (c, tag, lab, note) in enumerate(rows):
    on, sl = cat(c, tag)
    col = 'k' if tag == 'el' else ('C0' if 'normal' not in c else 'C3')
    for j, (x0, x1) in enumerate(((0, 15000), (13500, 15000))):
        a = ax[i, j]; m = (on >= x0) & (on <= x1)
        a.vlines(on[m], 0, sl[m], color=col, lw=0.6 if j == 0 else 1.2)
        a.plot(on[m], sl[m], '.', color=col, ms=2 if j == 0 else 4)
        a.set_xlim(x0, x1); a.set_ylim(1.5, 3.3); a.grid(alpha=0.3)
        if j == 0:
            a.set_ylabel('slip [m]'); a.text(0.01, 0.95, f'{lab}: {note}', transform=a.transAxes, va='top', fontsize=9,
                                             bbox=dict(fc='white', ec='none', alpha=0.8))
        if i < len(rows) - 1: a.set_xticklabels([])
ax[-1, 0].set_xlabel('time [yr]'); ax[-1, 1].set_xlabel('time [yr], last 1500 yr')
fig.suptitle('BP3 dip-60 fault, 25 m, 40 km plate over Maxwell half-space, sigma coupling and gravity on: '
             'mean slip per event vs onset time (T = elastic mean recurrence 85 / 95 yr)', fontsize=10)
fig.savefig(out + '.png', dpi=130, bbox_inches='tight'); print('wrote', out + '.png')
