"""Stability of the virtual surface per wavenumber (fd_plan.md, step B4): for the
x-slip -> xz-traction tables written by `fd3d_test surfsym <cfl> <m> <file>`, the
traction-free recursion of one Fourier mode is iterated (second-order NL scheme and
the first-order scheme) and the growth factor per step reported for every
wavenumber in the file. usage: fs3d_surf_stability.py <file> [nsteps]
Modes are indexed (i/8, j/8) of the Nyquist wavenumber along strike and across; the
checkerboard is (8/8, 4/8) i.e. i = 8, j = 4."""
import sys, numpy as np

def simulate(k0, w, eta, nl, nsteps=2000):
    L = len(w); K1 = w[0] + k0
    V = np.zeros(nsteps); m = np.zeros(nsteps); V[0] = 1.0; m[0] = 0.5 if nl else 1.0; delta = m[0]
    for n in range(1, nsteps):
        lm = min(L, n + 1)
        hist = np.dot(w[1:lm], m[n - 1 - np.arange(lm - 1)]) if lm > 1 else 0.0   # sum_{l=2..lm} W(l) m(n+1-l)
        if nl:
            B = 0.5*(K1 + eta) - eta; C = 0.5*(K1 + eta)
            V[n] = -(k0*delta + hist + C*V[n-1])/B; m[n] = 0.5*(V[n-1] + V[n])
        else:
            V[n] = -(k0*delta + hist)/K1; m[n] = V[n]
        delta += m[n]
        if not np.isfinite(V[n]) or abs(V[n]) > 1e100: return np.inf
    a = np.abs(V[nsteps//2:]).max(); b = max(np.abs(V[:nsteps//2]).max(), 1e-300)
    return (a/b)**(2.0/nsteps)

def analyse(path, nsteps=2000, verbose=True):
    head = open(path).readline()
    eta = float(head.split("eta")[1].split(",")[0]); ratio = float(head.split("cS dt/hs")[1].split(",")[0])
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    res = []
    for r in rows:
        i, j = int(r[0]), int(r[1]); k0 = float(r[2]); w = np.array(r[3:], float)
        res.append((i, j, simulate(k0, w, eta, True, nsteps), simulate(k0, w, eta, False, nsteps)))
    res = np.array(res)
    wn = res[np.argmax(res[:, 2])]; w1 = res[np.argmax(res[:, 3])]
    cb = res[(res[:, 0] == 8) & (res[:, 1] == 4)][0]
    if verbose:
        print(f"cS dt/hs {ratio:.2f}: worst NL growth {wn[2]:.4f} per step at (i, j) = ({int(wn[0])}, {int(wn[1])}); "
              f"worst first order {w1[3]:.4f} at ({int(w1[0])}, {int(w1[1])}); checkerboard NL {cb[2]:.4f}, first order {cb[3]:.4f}")
    return res

if __name__ == "__main__":
    analyse(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2000)
