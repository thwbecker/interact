#!/usr/bin/env python3
"""
rsf_series.py: read the all-cell binary time series written by rsf_solve
(-rsf_series_file), the event catalogs, and the geometry of one run.

File layout (little endian, written by rsf_output.c):
  int32 magic 0x52534631, int32 ncell, int32 nfield (5; 4 in older files,
  6 with two state variables, the sixth field log10 theta2 [s])
  records: float64 t[s], float32[ncell] log10|v| [m/s], float32[ncell]
           tau [MPa], float32[ncell] log10 theta [s], float32[ncell] slip [m],
           float32[ncell] sigma [MPa] (nfield 5; older files have 4)

Also provides resampling to a uniform grid and an event-aligned
extraction, the two forms intended for training.
"""
import json
import os
import warnings

import numpy as np

SEC_PER_YEAR = 365.25 * 86400.0
MAGIC = 0x52534631
FIELDS = ("logv", "tau", "logtheta", "slip", "sigma")


class Run:
    def __init__(self, rdir):
        self.dir = rdir
        with open(os.path.join(rdir, "params.json")) as f:
            self.par = json.load(f)
        g = np.loadtxt(os.path.join(rdir, "geom.in"), ndmin=2)
        self.x, self.y, self.z = g[:, 0], g[:, 1], g[:, 2]
        self.group = g[:, 7].astype(int)
        self.ngroup = self.group.max() + 1
        ab = np.loadtxt(os.path.join(rdir, "rsf.dat"), ndmin=2)
        self.a, self.b = ab[:, 0], ab[:, 1]
        self.read_series()
        self.read_catalogs()

    def read_series(self, fname="series.bin"):
        with open(os.path.join(self.dir, fname), "rb") as f:
            hdr = np.fromfile(f, dtype="<i4", count=3)
            if hdr[0] != MAGIC:
                raise IOError("bad magic in %s" % fname)
            n, nf = int(hdr[1]), int(hdr[2])
            if n != len(self.x):
                raise IOError("series has %i cells, geometry %i" % (n, len(self.x)))
            rec = np.dtype([("t", "<f8"), ("f", "<f4", (nf, n))])
            d = np.fromfile(f, dtype=rec)
        self.t = d["t"] / SEC_PER_YEAR          # [yr]
        self.logv = d["f"][:, 0, :]
        self.tau = d["f"][:, 1, :]
        self.logtheta = d["f"][:, 2, :]
        self.slip = d["f"][:, 3, :]
        # normal stress [MPa] was added as a fifth field (compression
        # positive); older files have four fields and a constant sigma
        if nf >= 5:
            self.sigma = d["f"][:, 4, :]
        else:
            self.sigma = np.full_like(self.tau, self.par.get("sigma", 0.0) / 1e6)
        # with two state variables (rsf_solve -rsf_two_state) a sixth
        # field holds log10 theta2 [s]; logtheta is then theta1
        self.logtheta2 = d["f"][:, 5, :] if nf >= 6 else None
        self.nfield = nf
        self.nrec = len(self.t)

    def read_catalogs(self):
        """global catalog and per-group catalogs (if -rsf_monitor_by_group)"""
        self.cat = _load_cat(os.path.join(self.dir, "rsf_catalog.dat"))
        self.gcat = []
        for k in range(self.ngroup):
            fn = os.path.join(self.dir, "rsf_catalog.g%03i.dat" % k)
            self.gcat.append(_load_cat(fn) if os.path.exists(fn) else None)

    def group_stats(self, field):
        """mean and std over the cells of each group, shape (ngroup, nrec)"""
        f = getattr(self, field)
        mean = np.empty((self.ngroup, self.nrec))
        std = np.empty((self.ngroup, self.nrec))
        for k in range(self.ngroup):
            sel = self.group == k
            mean[k] = f[:, sel].mean(axis=1)
            std[k] = f[:, sel].std(axis=1)
        return mean, std

    def resample(self, dt_yr, fields=FIELDS, t0=0.0, t1=None):
        """linear interpolation of all cells onto a uniform grid; returns
        (t, dict of arrays (nt, ncell)).  log quantities are interpolated
        as logs.  Event-scale detail is lost unless dt_yr is tiny; use
        event_windows for that."""
        if t1 is None:
            t1 = self.t[-1]
        tu = np.arange(t0, t1, dt_yr)
        out = {}
        for fld in fields:
            f = getattr(self, fld)
            o = np.empty((len(tu), f.shape[1]), dtype=np.float32)
            for i in range(f.shape[1]):
                o[:, i] = np.interp(tu, self.t, f[:, i])
            out[fld] = o
        return tu, out

    def event_windows(self, pre_s=10.0, post_s=200.0, fields=FIELDS):
        """for each catalog event, the records between onset - pre_s and
        arrest + post_s (seconds), as a list of dicts with t[s] relative
        to onset"""
        wins = []
        if self.cat is None:
            return wins
        ts = self.t * SEC_PER_YEAR
        for ev in self.cat:
            t_on, t_ar = ev[1] * SEC_PER_YEAR, ev[2] * SEC_PER_YEAR
            sel = (ts >= t_on - pre_s) & (ts <= t_ar + post_s)
            w = {"ievent": int(ev[0]), "t": ts[sel] - t_on, "Mw": ev[12]}
            for fld in fields:
                w[fld] = getattr(self, fld)[sel]
            wins.append(w)
        return wins


def _load_cat(fn):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")   # empty catalog files are legitimate
            c = np.loadtxt(fn, ndmin=2)
    except (OSError, ValueError):
        return None
    if c.size == 0:
        return None
    return c


def cat_columns(cat, group=False):
    """column indices for the two catalog layouts"""
    if group:
        return dict(ev=0, onset=1, arrest=2, ncell=5, area=6, mean_slip=7, max_slip=8,
                    mean_drop=9, max_drop=10, peak_v=11, M0=12, Mw=13)
    return dict(ev=0, onset=1, arrest=2, duration=3, ncell=4, area=5, mean_slip=6, max_slip=7,
                mean_drop=8, max_drop=9, peak_v=10, M0=11, Mw=12)


def export_training(run, fname, dt_yr=0.5, pre_s=10.0, post_s=300.0):
    """write a compressed npz with the uniform interseismic series, the
    event windows (ragged, stored as object arrays) and the catalog"""
    tu, u = run.resample(dt_yr)
    wins = run.event_windows(pre_s, post_s)
    np.savez_compressed(fname, t_uniform=tu, group=run.group, x=run.x, y=run.y, a=run.a, b=run.b,
                        params=json.dumps(run.par), catalog=run.cat if run.cat is not None else np.zeros((0, 13)),
                        n_events=len(wins),
                        event_t=np.array([w["t"] for w in wins], dtype=object),
                        **{"u_" + k: v for k, v in u.items()},
                        **{"ev_" + k: np.array([w[k] for w in wins], dtype=object)
                           for k in FIELDS})


if __name__ == "__main__":
    import sys
    r = Run(sys.argv[1] if len(sys.argv) > 1 else ".")
    print("%s: %i cells in %i group(s), %i records to %.1f yr, %i events" %
          (r.dir, len(r.x), r.ngroup, r.nrec, r.t[-1], 0 if r.cat is None else len(r.cat)))
