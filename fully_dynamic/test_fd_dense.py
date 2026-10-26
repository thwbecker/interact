"""Check the C dense (non-parallel) pair tables (src/testing/rsf_solve/test_fd_dense.c)
against the prototype's per-receiver tables inplane.W_table / antiplane.W_table_ap.
usage: test_fd_dense.py <angle_deg> <mode 0|1> <file written by test_fd_dense>"""
import sys, numpy as np
import inplane as ip, antiplane as ap
ang = np.radians(float(sys.argv[1])); mode = int(sys.argv[2]); fn = sys.argv[3]
alpha, beta, rho, dx = 6000.0, 3464.0, 2670.0, 100.0; dt = 0.5*dx/beta; ns, nr = 6, 5
ss = np.array([1.0, 0.0]); sr = np.array([np.cos(ang), np.sin(ang)])
Xs = np.array([[(j + 0.5)*dx, 0.0] for j in range(ns)])
Xr = np.array([[ns*dx + (i + 0.5)*dx*sr[0], (i + 0.5)*dx*sr[1]] for i in range(nr)])
d = np.loadtxt(fn)
nlag = 61
ref = {}
for j in range(ns):
    if mode == 0: Ws, Wn, Ks, Kn = ip.W_table(Xs[j], ss, 0.5*dx, Xr, sr, nlag, dt, alpha, beta, rho)
    else:
        out = ap.W_table_ap(Xs[j], ss, 0.5*dx, Xr, sr, nlag, dt, alpha, beta, rho)
        Ws, Wn, Ks, Kn = out[0], out[1], out[2], out[3]
    for i in range(nr):
        ref[(0, i, j)] = (Ks[i], Kn[i])
        for l in (1, 2, 5, 20, 60): ref[(l, i, j)] = (Ws[i, l - 1], Wn[i, l - 1])
err = 0.0; scale = 0.0
for row in d:
    l, i, j = int(row[0]), int(row[1]), int(row[2])
    r = ref[(l, i, j)]
    err = max(err, abs(row[3] - r[0]), abs(row[4] - r[1]) if mode == 0 else 0.0)
    scale = max(scale, abs(r[0]))
print(f"test_fd_dense.py: angle {sys.argv[1]} deg, mode {mode}: {len(d)} entries, max |C - prototype| {err:.3e} of {scale:.3e}: {'PASS' if err <= 1e-9*scale else 'FAIL'}")
