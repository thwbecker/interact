#!/usr/bin/env python3
"""bp3ve_analyze.py: recurrence and cycle properties of the BP3 VE sweep
from the catalogs alone (cat_<tag>.dat, one row per event), so that the
monitor and station files need not be moved around.

usage: bp3ve_analyze.py [root] [twin] [out] [min_ncell] [tol] [maxper]
       root      directory holding the ve_* case directories    (default .)
       twin      fixed-window start, yr, for the comparison column (2500)
       out       prefix for outputs                                (bp3ve)
       min_ncell drop events rupturing fewer cells than this       (0)
       tol       interval tolerance for cycle detection, yr        (0.5)
       maxper    longest period searched                           (40)

Catalog columns (rsf_solve SEAS-style): 0 index, 1 onset, 2 arrest,
3 duration, 4 n cells, 5 area, 6 mean slip, 7 max slip, 8 mean drop,
9 max drop, 10 peak slip rate, 11 M0, 12 Mw.

Settled state: for each period p = 1..maxper the longest tail of the
interval series dt with |dt[i] - dt[i+p]| < tol is found; the p with
the longest tail wins (ties to the smaller p).  The tail is accepted as
settled if it holds at least two full periods; its first onset is
t_settle.  Means over the settled state use whole periods only.  A
second pass at 4 tol flags runs whose cycle is recognisable but whose
members still drift ("creep": the laminar phase of an intermittent run,
or a slow approach to the attractor); the drift of the fastest-moving
member is reported in yr per cycle.

Outputs (prefix out):
  <out>_summary.txt / .csv   one row per case and tag
  <out>_recurrence.png       mean recurrence change vs tM, settled and
                             fixed window, all cases
  <out>_members.png          attractor members vs tM per case
  <out>_intervals_<case>.png interval series dt(t), one panel per tM
  <out>_props.png            cycle properties vs tM relative to elastic
  <out>_settle.png           t_settle and number of settled cycles vs tM
  <out>_cats.tgz             the catalogs, branch.txt, stop_*.yr and the
                             summary, small enough to send
"""
import sys, os, re, glob, tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = sys.argv[1] if len(sys.argv) > 1 else '.'
twin = float(sys.argv[2]) if len(sys.argv) > 2 else 2500.0
out = sys.argv[3] if len(sys.argv) > 3 else 'bp3ve'
min_ncell = int(sys.argv[4]) if len(sys.argv) > 4 else 0
tol = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
maxper = int(sys.argv[6]) if len(sys.argv) > 6 else 40

COLS = dict(onset=1, arrest=2, dur=3, ncell=4, area=5, slip=6, maxsl=7,
            drop=8, maxdr=9, vmax=10, M0=11, Mw=12)
PROPS = ['slip', 'maxsl', 'drop', 'maxdr', 'ncell', 'Mw', 'dur', 'vmax']


def tmval(tag):
    return 0.0 if tag == 'el' else float(tag[2:])


def load_cat(fn):
    d = np.loadtxt(fn, comments='#', ndmin=2)
    if d.shape[1] < 13:
        raise ValueError(f'{fn}: {d.shape[1]} columns, expected 13')
    d = d[np.argsort(d[:, COLS['onset']])]
    if min_ncell > 0:
        d = d[d[:, COLS['ncell']] >= min_ncell]
    return d


def settled_tail(dt, tol, maxper):
    """Return (p, nint) for the longest periodic tail of dt; p = 0 if none."""
    n = len(dt)
    best_p, best_n = 0, 0
    for p in range(1, maxper + 1):
        if n < 2 * p:
            break
        k = 0                       # number of satisfied comparisons
        while k < n - p and abs(dt[n - 1 - k] - dt[n - 1 - k - p]) < tol:
            k += 1
        nint = k + p                # intervals covered by the tail
        if nint < 2 * p:
            continue
        if nint > best_n + p:       # clearly longer than what we have
            best_p, best_n = p, nint
    return best_p, best_n


def analyse(d):
    on = d[:, COLS['onset']]
    dt = np.diff(on)
    r = dict(nev=len(on), t_first=on[0], t_last=on[-1])
    # fixed window
    w = on[1:] > twin
    r['n_win'] = int(w.sum())
    r['T_win'] = dt[w].mean() if w.any() else np.nan
    r['T_win_std'] = dt[w].std() if w.any() else np.nan
    # settled tail; a second, loose pass (4 tol) catches a cycle whose
    # members are still creeping (laminar phase of an intermittent run)
    p, nint = settled_tail(dt, tol, maxper)
    pl, nl = settled_tail(dt, 4 * tol, maxper)
    r['per'] = p
    r['creep'] = ''
    if pl and (nl > 2 * max(nint, 2 * pl)):
        tl = dt[-nl:]
        # drift of each phase over the loose tail, yr per cycle
        drift = [(tl[k::pl][-1] - tl[k::pl][0]) / max(1, len(tl[k::pl]) - 1) for k in range(pl)]
        j = int(np.argmax(np.abs(drift)))
        r['creep'] = f'creep p={pl} {nl // pl}cyc member {tl[j::pl][-1]:.0f} {drift[j]:+.2f}/cyc'
    if p:
        nfull = (nint // p) * p
        tail = dt[-nfull:]
        r['n_settled'] = nint
        r['n_cycles'] = nint // p
        r['t_settle'] = on[len(on) - 1 - nint]
        r['T_set'] = tail.mean()
        r['members'] = np.sort(dt[-p:])
        # intervals' end events: events on[-nfull:]
        ev = d[len(on) - nfull:]
    else:
        r['n_settled'] = 0
        r['n_cycles'] = 0
        r['t_settle'] = np.nan
        r['T_set'] = np.nan
        r['members'] = None
        ev = d[1:][w] if w.any() else d[1:]
    for k in PROPS:
        r[k] = ev[:, COLS[k]].mean() if len(ev) else np.nan
    r['on'] = on
    r['dt'] = dt
    return r


# ---------------------------------------------------------------- read
cases = sorted(c for c in os.listdir(root)
               if os.path.isdir(os.path.join(root, c)) and c.startswith('ve_'))
if not cases:
    sys.exit(f'{root}: no ve_* directories')
res = {}
for c in cases:
    fns = glob.glob(os.path.join(root, c, 'cat_*.dat'))
    tags = sorted((tmval(re.sub(r'.*cat_(.*)\.dat$', r'\1', f)), f) for f in fns)
    res[c] = {}
    for tm, f in tags:
        try:
            d = load_cat(f)
        except Exception as e:
            print(f'skip {f}: {e}', file=sys.stderr)
            continue
        if len(d) < 3:
            print(f'skip {f}: {len(d)} events', file=sys.stderr)
            continue
        r = analyse(d)
        sf = os.path.join(root, c, 'stop_' + os.path.basename(f)[4:-4] + '.yr')
        r['stop'] = float(open(sf).read().split()[0]) if os.path.isfile(sf) else np.nan
        res[c][tm] = r
    if not res[c]:
        del res[c]

# ---------------------------------------------------------------- tables
hdr = (f'{"case":<16}{"tag":>7}{"nev":>5}{"t_last":>9}{"stop":>7}'
       f'{"per":>4}{"ncyc":>5}{"t_settle":>9}{"T_set":>8}{"dT_set":>8}'
       f'{"T_win":>8}{"dT_win":>8}'
       + ''.join(f'{k:>8}' for k in PROPS) + '  members')
lines = [f'# bp3ve_analyze: root={root} twin={twin:g} min_ncell={min_ncell} '
         f'tol={tol:g} maxper={maxper}',
         '# dT_* in percent relative to the el run of the same case; '
         'properties are means over the settled state (or the window)',
         hdr]
csv = ['case,tm,nev,t_last,stop,per,ncyc,t_settle,T_set,dT_set,T_win,dT_win,'
       + ','.join(PROPS) + ',members']
for c in res:
    ref = res[c].get(0.0)
    for tm in sorted(res[c]):
        r = res[c][tm]
        tag = 'el' if tm == 0 else f'tm{tm:g}'
        r['dT_set'] = 100 * (r['T_set'] / ref['T_set'] - 1) if ref and ref['per'] else np.nan
        r['dT_win'] = 100 * (r['T_win'] / ref['T_win'] - 1) if ref else np.nan
        mem = ' '.join(f'{x:.1f}' for x in r['members']) if r['members'] is not None else 'aperiodic'
        if r['creep']:
            mem += '  [' + r['creep'] + ']'
        f = lambda x, w=8, p=2: f'{x:{w}.{p}f}' if np.isfinite(x) else ' ' * (w - 1) + '-'
        lines.append(f'{c:<16}{tag:>7}{r["nev"]:5d}{r["t_last"]:9.1f}{f(r["stop"], 7, 0)}'
                     f'{r["per"]:4d}{r["n_cycles"]:5d}{f(r["t_settle"], 9, 0)}{f(r["T_set"])}'
                     f'{f(r["dT_set"])}{f(r["T_win"])}{f(r["dT_win"])}'
                     + ''.join(f(r[k], 8, 3 if k in ('slip', 'maxsl', 'drop', 'Mw') else 1)
                               for k in PROPS)
                     + f'  {mem}')
        csv.append(','.join([c, f'{tm:g}', str(r['nev']), f'{r["t_last"]:.2f}', f'{r["stop"]:g}',
                             str(r['per']), str(r['n_cycles']), f'{r["t_settle"]:.1f}',
                             f'{r["T_set"]:.3f}', f'{r["dT_set"]:.2f}', f'{r["T_win"]:.3f}',
                             f'{r["dT_win"]:.2f}']
                            + [f'{r[k]:.4g}' for k in PROPS] + ['"' + mem + '"']))
txt = '\n'.join(lines) + '\n'
open(out + '_summary.txt', 'w').write(txt)
open(out + '_summary.csv', 'w').write('\n'.join(csv) + '\n')
print(txt)

# ---------------------------------------------------------------- plots
style = {'ve_g1': ('C0', 'o', '-'), 've_g0': ('C1', 's', '-'), 've_shear': ('C2', '^', '-'),
         've_g1_normal': ('C0', 'o', '--'), 've_g0_normal': ('C1', 's', '--'),
         've_shear_normal': ('C2', '^', '--')}


def sty(c):
    return style.get(c, ('k', 'x', ':'))


def kw(c):
    """plot kwargs: normal-branch cases get open markers and dashed lines"""
    col, mk, ls = sty(c)
    d = dict(color=col, marker=mk, ls=ls)
    if c.endswith('_normal'):
        d['mfc'] = 'none'
    return d


def tms_of(c):
    return sorted(t for t in res[c] if t > 0)


# recurrence change vs tM
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for c in res:
    col, mk, ls = sty(c)
    tms = tms_of(c)
    if not tms:
        continue
    y1 = [res[c][t]['dT_set'] for t in tms]
    y2 = [res[c][t]['dT_win'] for t in tms]
    ax[0].plot(tms, y1, label=c, **kw(c))
    ax[1].plot(tms, y2, label=c, **kw(c))
    # mark unsettled runs
    bad = [t for t in tms if not res[c][t]['per']]
    if bad:
        ax[0].plot(bad, [0] * len(bad), color=col, marker='x', ls='none', ms=9)
for a, t in zip(ax, ['settled state (x: aperiodic, not plotted)',
                     f'fixed window t > {twin:g} yr']):
    a.axhline(0, color='k', lw=0.5)
    a.set_xscale('log')
    a.set_xlabel('t_M [yr]')
    a.set_title(t)
    a.grid(alpha=0.3)
ax[0].set_ylabel('mean recurrence change [%]')
ax[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out + '_recurrence.png', dpi=130)
plt.close(fig)

# attractor members vs tM
nc = len(res)
fig, ax = plt.subplots(1, nc, figsize=(3.6 * nc, 4), sharey=True, squeeze=False)
for a, c in zip(ax[0], res):
    col, mk, ls = sty(c)
    ref = res[c].get(0.0)
    if ref and ref['members'] is not None:
        for m in ref['members']:
            a.axhline(m, color='gray', lw=0.8, ls=':')
    for t in tms_of(c):
        r = res[c][t]
        if r['members'] is not None:
            a.plot([t] * len(r['members']), r['members'], ls='none', **{k: v for k, v in kw(c).items() if k != 'ls'})
            a.plot(t, r['T_set'], color=col, marker='_', ms=14, mew=2, ls='none')
        else:
            dtw = r['dt'][r['on'][1:] > twin]
            a.plot([t] * len(dtw), dtw, color='r', marker='.', ls='none', alpha=0.4)
    a.set_xscale('log')
    a.set_title(c, fontsize=10)
    a.set_xlabel('t_M [yr]')
    a.grid(alpha=0.3)
ax[0][0].set_ylabel('interval [yr]')
fig.suptitle('attractor members (dotted: elastic members; bar: settled mean; '
             'red dots: aperiodic, all intervals in window)', fontsize=10)
fig.tight_layout()
fig.savefig(out + '_members.png', dpi=130)
plt.close(fig)

# interval series dt(t): one figure per case, one panel per tM, el behind
for c in res:
    tms = tms_of(c)
    if not tms:
        continue
    ref = res[c].get(0.0)
    ncol = 3
    nrow = (len(tms) + ncol - 1) // ncol
    fig, ax = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 2.3 * nrow),
                           sharex=True, sharey=True, squeeze=False)
    for a, t in zip(ax.ravel(), tms):
        r = res[c][t]
        if ref:
            a.plot(ref['on'][1:], ref['dt'], color='0.75', lw=0.7)
        a.plot(r['on'][1:], r['dt'], color=sty(c)[0], lw=0.8)
        if r['per']:
            a.axvline(r['t_settle'], color='r', lw=0.8, ls=':')
            lab = f'p={r["per"]}, {r["n_cycles"]} cycles, T={r["T_set"]:.1f} ({r["dT_set"]:+.1f}%)'
        else:
            lab = 'aperiodic'
        a.set_title(f't_M = {t:g} yr:  {lab}', fontsize=9)
        a.grid(alpha=0.3)
    for a in ax.ravel()[len(tms):]:
        a.set_visible(False)
    for a in ax[-1]:
        a.set_xlabel('onset time [yr]')
    for a in ax[:, 0]:
        a.set_ylabel('interval [yr]')
    fig.suptitle(f'{c}: interevent intervals (gray: elastic; red dotted: t_settle)', fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{out}_intervals_{c}.png', dpi=120)
    plt.close(fig)

# cycle properties relative to elastic
fig, ax = plt.subplots(2, 4, figsize=(14, 6.5), sharex=True)
for a, k in zip(ax.ravel(), PROPS):
    for c in res:
        col, mk, ls = sty(c)
        ref = res[c].get(0.0)
        tms = tms_of(c)
        if not ref or not tms:
            continue
        y = [100 * (res[c][t][k] / ref[k] - 1) for t in tms]
        a.plot(tms, y, label=c, ms=4, **kw(c))
    a.axhline(0, color='k', lw=0.5)
    a.set_xscale('log')
    a.set_title(f'{k} change [%]', fontsize=10)
    a.grid(alpha=0.3)
for a in ax[1]:
    a.set_xlabel('t_M [yr]')
ax[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(out + '_props.png', dpi=130)
plt.close(fig)

# settling
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for c in res:
    col, mk, ls = sty(c)
    tms = tms_of(c)
    if not tms:
        continue
    ax[0].plot(tms, [res[c][t]['t_settle'] for t in tms], label=c, **kw(c))
    ax[1].plot(tms, [res[c][t]['n_cycles'] for t in tms], label=c, **kw(c))
for c in res:
    for t in tms_of(c):
        if np.isfinite(res[c][t]['stop']):
            ax[0].axhline(res[c][t]['stop'], color='gray', lw=0.5)
            break
ax[0].set_ylabel('t_settle [yr]  (gray: run length)')
ax[1].set_ylabel('settled cycles (whole periods)')
for a in ax:
    a.set_xscale('log')
    a.set_xlabel('t_M [yr]')
    a.grid(alpha=0.3)
ax[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out + '_settle.png', dpi=130)
plt.close(fig)

# ---------------------------------------------------------------- bundle
with tarfile.open(out + '_cats.tgz', 'w:gz') as tf:
    for c in res:
        for pat in ('cat_*.dat', 'branch.txt', 'stop_*.yr', 'bp3_rsf.dat'):
            for f in glob.glob(os.path.join(root, c, pat)):
                tf.add(f, arcname=os.path.join(c, os.path.basename(f)))
    for ext in ('_summary.txt', '_summary.csv'):
        tf.add(out + ext, arcname=os.path.basename(out + ext))
print(f'wrote {out}_summary.txt/.csv, {out}_{{recurrence,members,props,settle}}.png, '
      f'{out}_intervals_<case>.png, '
      f'{out}_cats.tgz ({os.path.getsize(out + "_cats.tgz") / 1e6:.2f} MB)')
