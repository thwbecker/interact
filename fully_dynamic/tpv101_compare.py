"""SCEC TPV101 comparison: reads SCEC time-series files (HTML pages from the SCEC
code-validation site, e.g. Kaneko et al. SPECFEM3D 100 m: faultst000dp075, faultst-090dp075,
faultst120dp030) and the station files written by
`FD3D_STFILES=1 fd3d_test tpv101 <dx> - <tend>` (tpv101_<station>.dat: time, slip, slip rate,
shear stress in MPa), and prints rupture time (slip rate > 1 mm/s), peak slip rate and slip at
3, 6, 9 and 15 s. Stations at negative along-strike distance are compared with ours at the
mirrored position (the problem is symmetric).
usage: tpv101_compare.py <scec files ...> -- <our run directories ...>"""
import sys, re, html, os, numpy as np

def read_scec(path):
    s = open(path, errors="replace").read()
    t = html.unescape(re.sub(r"<[^>]+>", "\n", s))
    lines = [l for l in t.split("\n") if l.strip()]
    hdr = {}
    for l in lines[:80]:
        m = re.match(r"\s*#\s*(\w+)\s*=\s*(.*)", l)
        if m: hdr[m.group(1)] = m.group(2).strip()
    rows = [l.split() for l in lines if re.match(r"^\s*[-0-9.eE+]+(\s+[-0-9.eE+]+){3,}\s*$", l)]
    return hdr, (np.array(rows, float) if rows else None)

def metrics(t, sl, sr, times=(3.0, 6.0, 9.0, 15.0)):
    on = np.abs(sr) > 1e-3
    tr = t[np.argmax(on)] if on.any() else np.nan
    return tr, np.abs(sr).max(), [np.interp(x, t, sl) for x in times]

if __name__ == "__main__":
    args = sys.argv[1:]; k = args.index("--") if "--" in args else len(args)
    for path in args[:k]:
        hdr, a = read_scec(path)
        if a is None: print(f"{path}: no data (a login page?)"); continue
        loc = hdr.get("location", ""); m = re.search(r"st(-?\d+)dp(\d+)", loc)
        if m:
            st = f"st{abs(int(m.group(1))):03d}dp{m.group(2)}"
        else:   # "on fault -9.00 km along strike,  7.50 km down-dip" (Dunham's files)
            m = re.search(r"(-?[\d.]+)\s*km along strike,\s*(-?[\d.]+)\s*km down-dip", loc)
            xs, zs = round(float(m.group(1))*10), round(float(m.group(2))*10)
            st = f"st{abs(xs):03d}dp{zs:03d}"; loc = f"faultst{xs:03d}dp{zs:03d}"
        tr, pk, sl = metrics(a[:, 0], np.abs(a[:, 1]), a[:, 2])
        who = (hdr.get('author', '?').split(',')[0].split()[-1] + ' ' + hdr.get('code', '?'))[:22]
        print(f"{loc:20s} {who:22s}: rupture {tr:6.3f} s, peak {pk:6.3f} m/s, slip at 3/6/9/15 s " + " ".join(f"{x:.3f}" for x in sl))
        for d in args[k + 1:]:
            f = os.path.join(d, f"tpv101_{st}.dat")
            if not os.path.exists(f): continue
            b = np.loadtxt(f); tr, pk, sl = metrics(b[:, 0], np.abs(b[:, 1]), b[:, 2])
            print(f"   {d:40s}: rupture {tr:6.3f} s, peak {pk:6.3f} m/s, slip at 3/6/9/15 s " + " ".join(f"{x:.3f}" for x in sl))
