#!/usr/bin/env python3
"""
mn19_exact_chain.py: independent integrator for the Miyake & Noda (2019)
uniform-Maxwell antiplane patch model, cross-checking rsf_solve -ve_mode 1.

Kernel: K(t) = K0 exp(-t/tc) exactly, via ONE memory state per cell:
  tau_i = tau0_i + sum_j K0_ij S_j,   S_j' = (V_j - Vpl) - S_j/tc
(K0 is the static antiplane kernel of the same buried finite strip used by
gen_mn19.py: cells along depth around 30 km, free surface included, which
is negligible at 52R distance).  tc = 0 is not allowed; tc <= 0 on the
command line means ELASTIC (S ignored, tau = tau0 + K0 D).

Friction: aging law with M&N parameters (a = 0.02 uniform, b smoothed
boxcar 0.025 in the patch, 0.001 floor outside per their Fig. 2,
sigma = 100 MPa, f0 = 0.6 at V0 = Vpl = 1e-9 m/s), sinh-regularized like
rsf_solve, radiation damping eta = mu/(2 cs), cs = 3000 m/s.  IC: their
Eq. (23) perturbation, Vini = 0.999 Vpl with the 1 percent gradient.

usage: mn19_exact_chain.py tc_yr tend_yr [ds_m] [RcR] [lenR]
prints onsets (vmax > 0.1 m/s with 0.05 hysteresis) and intervals.
"""
import sys
import numpy as np

tc_yr = float(sys.argv[1])
tend  = float(sys.argv[2])
ds    = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
RcR   = float(sys.argv[4]) if len(sys.argv) > 4 else 0.4
lenR  = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0

yr  = 3.15576e7
G   = 30.0e9; cs = 3000.0; eta = G/(2*cs)
sig = 100.0e6; f0 = 0.6; Vpl = 1.0e-9; V0 = Vpl
R   = 500.0
a0  = 0.02; bmax = a0/0.8; b_out = 0.001
L   = RcR*R*np.pi*(bmax*sig/G)*((bmax-a0)/bmax)**2
zdep = 30.0e3; wtr = 0.1*R
tc  = tc_yr*yr
elastic = tc_yr <= 0.0

n = int(round(lenR*R/ds))
x = (np.arange(n)+0.5)*ds - lenR*R/2.0
z = zdep + x
ed = zdep + np.arange(n+1)*ds - lenR*R/2.0
ax = np.abs(x)
w = np.where(ax < R-wtr, 1.0,
    np.where(ax > R+wtr, 0.0, 0.5*(1.0+np.cos(np.pi*(ax-(R-wtr))/(2*wtr)))))
b = b_out + (bmax-b_out)*w
a = np.full(n, a0)
patch = w > 0.5

def F0(zz, zp): return 1.0/(zz-zp)-1.0/(zz+zp)
K0 = np.zeros((n, n))
for i in range(n):
    K0[i, :] = (G/(2*np.pi))*(F0(z[i], ed[1:])-F0(z[i], ed[:-1]))
if K0[0, 0] > 0:
    K0 = -K0

tau0 = sig*f0 + eta*Vpl            # steady state at Vpl (V0 = Vpl)
# IC per their Eq. (23): ln(Vpl th/L) = eps(x), theta above steady state
eps = np.log(1.0/0.999)*(1.0 + 0.01*x/(2*R))
lth0 = np.log(L/Vpl) + eps

def vsolve(tau, lth, Vg):
    Phi = f0 + b*(np.log(V0*np.exp(lth)/L))
    E = np.exp(Phi/a)
    V = Vg.copy()
    for _ in range(80):
        s = V/(2*V0)*E
        f = a*sig*np.arcsinh(s) + eta*V - tau
        df = a*sig*E/(2*V0)/np.sqrt(1.0+s*s) + eta
        Vn = V - f/df
        V = np.where(Vn > 0, Vn, V*0.5)
        if np.max(np.abs(f/(a*sig))) < 1e-12:
            break
    return V

def rhs(t, y, Vg):
    lth = y[:n]
    if elastic:
        D = y[n:2*n]
        tau = tau0 + K0@D
    else:
        S = y[n:2*n]
        tau = tau0 + K0@S
    V = vsolve(tau, lth, Vg)
    th = np.exp(lth)
    dlth = 1.0/th - V/L            # aging law, d ln(theta)/dt
    q = V - Vpl
    dy = np.empty_like(y)
    dy[:n] = dlth
    dy[n:2*n] = q if elastic else q - S/tc
    return dy, V

y0 = np.concatenate([lth0, np.zeros(n)])
rtol, atol_lth, atol_D = 1e-6, 1e-8, 1e-4
t, y = 0.0, y0.copy()
Vg = np.full(n, Vpl)
dt = 1.0
tend_s = tend*yr
onsets = []
slipping = False
k1, V1 = rhs(t, y, Vg)
nstep = naccept = 0
while t < tend_s and nstep < 8_000_000:
    nstep += 1
    k2, V2 = rhs(t+0.5*dt, y+0.5*dt*k1, V1)
    k3, V3 = rhs(t+0.75*dt, y+0.75*dt*k2, V2)
    ynew = y + dt*(2.0*k1 + 3.0*k2 + 4.0*k3)/9.0
    k4, V4 = rhs(t+dt, ynew, V3)
    yerr = dt*(-5.0*k1/72.0 + k2/12.0 + k3/9.0 - k4/8.0)
    sc = np.empty_like(y)
    sc[:n] = atol_lth + rtol*np.abs(ynew[:n])
    sc[n:] = atol_D + rtol*np.abs(ynew[n:])
    err = np.sqrt(np.mean((yerr/sc)**2))
    if err <= 1.0:
        t += dt; y = ynew; k1, V1 = k4, V4
        naccept += 1
        vmx = V1.max()
        if (not slipping) and vmx > 0.1:
            slipping = True; onsets.append(t/yr)
        elif slipping and vmx < 0.05:
            slipping = False
    dt = dt*min(4.0, max(0.2, 0.9*err**(-1.0/3.0)))
    dt = min(dt, tend_s-t+1.0) if t < tend_s else dt

onsets = np.array(onsets)
print(f"tc={tc_yr} ds={ds} RcR={RcR} lenR={lenR} n={n} L={L:.4e}: "
      f"{len(onsets)} onsets, {naccept}/{nstep} steps")
print(" onsets [yr]  :", " ".join(f"{x:.3f}" for x in onsets))
if len(onsets) > 2:
    print(" intervals[yr]:", " ".join(f"{d:.3f}" for d in np.diff(onsets)))
