"""Checks of the DWN half-plane kernel (halfplane_dwn.py, fd_plan.md M10): dynamic, against the
virtual-surface prototype (densesys,
kinematic: unit slip rate on one fault element, zero on another; surface solved each step).
Usage: check_halfplane_dwn.py [lags] [fault top depth]; DWN_N sets the number of frequencies;
TDD_LIB the tdd library of fs_static2d.py."""
import sys, os, time
import numpy as np
os.environ.setdefault('TDD_LIB', '/tmp/tdd/libtdd.so')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import halfplane_dwn as hd
import densesys as ds
al, be, rho = 6000., 3464., 2670.
hf = 250.0; a = hf/2; dt = 0.5*hf/be; L = int(sys.argv[1]) if len(sys.argv) > 1 else 120
dip = np.radians(30.0)
sd = np.array([np.cos(dip), -np.sin(dip)])            # down-dip tangent, y up
top = float(sys.argv[2]) if len(sys.argv) > 2 else 750.0
Xs = np.array([0.0, -top]) + 0.5*hf*sd                 # source element (fault element 0)
Xr = np.array([0.0, -top]) + 2.5*hf*sd                 # receiver element (2 elements further down)
def build(surf, hs=125.0, Ls=10e3):
    X = [Xs, Xr]; s = [sd, sd]; aa = [a, a]; kind = [0, 0]
    if surf:
        n = int(2*Ls/hs)
        for j in range(n):
            X.append(np.array([-Ls + (j + 0.5)*hs, 0.0])); s.append(np.array([1.0, 0.0])); aa.append(hs/2); kind.append(1)
    return ds.DenseSys(np.array(X), np.array(s), np.array(aa), np.array(kind), dt, al, be, rho, L)
res = {}
for surf in (False, True):
    t0 = time.time()
    d = build(surf)
    m = ds.Marcher(d, nl=False, vf_fun=lambda t, f: np.array([1.0, 0.0]))
    tau = []
    for n in range(L):
        m.step(); tau.append(m.tau[d.idx[(1, 0)]])
    res[surf] = np.array(tau)
    print('densesys surface %s: %.0f s' % (surf, time.time() - t0), flush=True)
ts = dt*np.arange(1, L + 1)
# DWN: (x, z) frame, z = -y
def dwn(part):
    S = hd.ramp_stress(Xr[0], -Xr[1], Xs[0], -Xs[1], (sd[0], -sd[1]), a, ts, al, be, rho, part=part,
                       T=2.5*ts.max(), N=int(os.environ.get("DWN_N", "2048")), wI_T=18.0)
    # traction on the receiver: normal (y up frame) n = rot90(s) = (-s_y, s_x); stress tensor in the y-up frame
    sxx, sxz, szz = S[:, 0], S[:, 1], S[:, 2]
    sxy = -sxz; syy = szz
    n = np.array([-sd[1], sd[0]])
    tx = sxx*n[0] + sxy*n[1]; ty = sxy*n[0] + syy*n[1]
    return tx*sd[0] + ty*sd[1]
t0 = time.time()
Df = dwn("full"); Dt = dwn("total")
print('dwn %.0f s' % (time.time() - t0))
# calibrate the sign with the full plane, then compare the surface effect (difference) and the totals
sgn = np.sign(np.dot(Df, res[False]))
print('full plane: sign %+d, largest |DWN - densesys| %.3e of %.3e' % (sgn, np.abs(sgn*Df - res[False]).max(), np.abs(res[False]).max()))
eff_d = res[True] - res[False]; eff_w = sgn*(Dt - Df)
for n in range(0, L, max(1, L//12)):
    print('t %.3f  surface effect: densesys %+.4e  DWN %+.4e   total densesys %+.4e DWN %+.4e' % (ts[n], eff_d[n], eff_w[n], res[True][n], sgn*Dt[n]))
print('largest |difference of totals| %.3e of %.3e; of the surface effect %.3e of %.3e' % (np.abs(res[True] - sgn*Dt).max(), np.abs(res[True]).max(), np.abs(eff_d - eff_w).max(), np.abs(eff_d).max()))
