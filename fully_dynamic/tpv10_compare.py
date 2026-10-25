"""SCEC TPV10 comparison: reads the SCEC time-series files (HTML pages as downloaded from the
SCEC code-validation site, e.g. Kaneko et al. SPECFEM3D 100 m and Dunham WENORUP-2D 18.75 m)
and the on-fault station files written by `fd3d_test tpv10 <dx> <strike> <q2d> <margin> <prefix>`
(columns: time, up-dip slip, up-dip slip rate, up-dip shear stress in MPa, normal stress in MPa,
SCEC signs), and prints rupture time (slip rate > 1 mm/s), peak slip rate, and slip at a few
times for each station.
usage: tpv10_compare.py <scec files ...> -- <our prefixes ...>"""
import sys, re, html, numpy as np

def read_scec(path):
    s = open(path, errors="replace").read()
    t = html.unescape(re.sub(r"<[^>]+>", "\n", s))
    lines = [l for l in t.split("\n") if l.strip()]
    hdr = {}
    for l in lines[:80]:
        m = re.match(r"\s*#\s*(\w+)\s*=\s*(.*)", l)
        if m: hdr[m.group(1)] = m.group(2).strip()
    rows = [l.split() for l in lines if re.match(r"^\s*[-0-9.eE+]+(\s+[-0-9.eE+]+){3,}\s*$", l)]
    return hdr, np.array(rows, float)

def downdip_km(loc):
    m = re.search(r"dp(\d+)", loc)
    if m: return int(m.group(1))/10.0
    m = re.search(r"([\d.]+)\s*km down-dip", loc)
    return float(m.group(1)) if m else None

def metrics(t, sl, sr, times=(5.0, 8.0, 12.0, 15.0)):
    on = np.abs(sr) > 1e-3
    tr = t[np.argmax(on)] if on.any() else np.nan
    return tr, np.abs(sr).max(), [np.interp(x, t, sl) for x in times]

if __name__ == "__main__":
    args = sys.argv[1:]; k = args.index("--") if "--" in args else len(args)
    scec, ours = args[:k], args[k + 1:]
    for path in scec:
        hdr, a = read_scec(path)
        if a.shape[1] != 8: continue                     # on-fault files only
        dd = downdip_km(hdr.get("location", ""))
        tr, pk, sl = metrics(a[:, 0], np.abs(a[:, 4]), a[:, 5])
        print(f"{hdr.get('author','?')[:24]:24s} {hdr.get('code','?')[:16]:16s} {dd:5.1f} km: rupture {tr:6.3f} s, peak {pk:6.3f} m/s, slip at 5/8/12/15 s "
              + " ".join(f"{x:.3f}" for x in sl))
        for p in ours:
            try: b = np.loadtxt(f"{p}_st{int(round(dd*1000)):05d}.dat")
            except OSError: continue
            tr, pk, sl = metrics(b[:, 0], np.abs(b[:, 1]), b[:, 2])
            print(f"   {p:21s} {'':16s}{dd:5.1f} km: rupture {tr:6.3f} s, peak {pk:6.3f} m/s, slip at 5/8/12/15 s " + " ".join(f"{x:.3f}" for x in sl))
