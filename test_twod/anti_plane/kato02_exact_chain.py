#!/usr/bin/env python3
"""
kato02_exact_chain.py: independent integrator for the Kato (2002) model,
using his hereditary kernel (his Eqs 1-5) via Erlang memory chains, with
NO Prony fitting.  Cross-check for rsf_solve's -ve_mode 2 machinery, and
a laboratory for testing implementation variants of his Eq. (2)-(3).

  tau_i(t) = tau0_i + sum_j [ K0_ij D_j + sum_m Km_ij S_mj ] (+ damping in
  the friction solve),   D_j' = V_j - Vpl,
  S_mj = int_0^t A_m((t-s)/tr) (V_j - Vpl) ds  exactly via the chain
  S_1' = (D - S_1)/tr,  S_m' = (S_{m-1} - S_m)/tr   (A_m' = (A_{m-1}-A_m)/tr)

M = Kato's image truncation (10 in the paper); larger M -> exact layered
kernel; M = 0 elastic.  Friction: composite Kato-Tullis law with the same
sinh regularization as rsf_solve; V per cell by Newton.  Integrator:
adaptive Bogacki-Shampine RK3(2), like rsf_solve's -ts_rk_type 3bs.

Kernel modes (the diagnostic ladder for his Table 1):
  exact   correct Erlang weights A_m(x) = 1 - e^-x sum_{l=0}^{m-1} x^l/l!
          (Bonafede et al. 1984); this is what his Eq. (2) MEANS
  kato3   his Eq. (3) AS PRINTED: the sum starts at l = 1, so
          A_m(x) = P(m,x) + e^-x, A_m(0) = 1 (all images instantly on),
          A_1 == 1 always; implemented exactly (chains plus one shared
          exponential state)
  static  A_m == 1: permanently relaxed truncated kernel, no time
          dependence at all (the tr-independent limit of kato3)
  current the history-free shortcut tau = [K0 + sum_m A_m(t/tr) Km] D(t)
          (current kernel times current deficit, no convolution)
  single  correct Erlang weights but only image family A per order
          (F_m = 1/(z-2mh-z') - 1/(z+2mh+z'), depths 2mh+z' and mirror),
          the analog of the missing-family bug this project found in the
          repo kernels
  singleB the complementary family B only (depths 2mh-z' and mirror,
          F_m = -1/(z-2mh+z') + 1/(z+2mh-z')); for m = 1 these images
          hug the plate bottom at 2h-z', the near-field substrate
          reflection of the deep creep zone

usage: kato02_exact_chain.py M tr_yr tend_yr [n] [dc] [mode] [probe_file]
       probe_file: if given, writes "t_yr vmax v(15km) v(18km)" time series
"""
import sys
import numpy as np
from scipy.special import lambertw, gammainc
from scipy.optimize import brentq

M      = int(sys.argv[1])
tr_yr  = float(sys.argv[2])
tend   = float(sys.argv[3])
n      = int(sys.argv[4])   if len(sys.argv) > 4 else 80
dc     = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
mode   = sys.argv[6]        if len(sys.argv) > 6 else "exact"
probe  = sys.argv[7]        if len(sys.argv) > 7 else None
assert mode in ("exact", "kato3", "static", "current", "single", "singleB")

yr   = 3.15576e7
G    = 30.0e9; cs = 3270.0; eta = G/(2*cs)
Vpl  = 35.0e-3/yr; Vc = 1e-8; V0 = 1e-6; f0 = 0.6
h    = 20e3; dz = h/n; tr = tr_yr*yr
rho, rhow, grav = 2800.0, 1000.0, 9.8

zc = (np.arange(n)+0.5)*dz
ed = np.arange(n+1)*dz
# Kato Fig. 2, DIGITIZED (same as gen_kato02.py; b floored at 0.001 where
# his a-b = a means b = 0, and a-b kept consistent with the floor)
zk_amb = np.array([0.0, 2.5, 5.0, 10.0, 12.0, 20.0])*1e3
va_amb = np.array([0.004, 0.0, -0.004, -0.004, 0.0121, 0.0345])
amb = np.interp(zc, zk_amb, va_amb)
a   = np.where(zc <= 12e3, 0.0121, 0.0121 + (0.0345-0.0121)*(zc-12e3)/8e3)
b   = np.maximum(a - amb, 0.001)
amb = a - b
sig = (rho-rhow)*grav*zc

def F0(z, zp): return 1.0/(z-zp)-1.0/(z+zp)
def Fm(z, zp, m):
    if mode == "single":
        return 1.0/(z-2*m*h-zp)-1.0/(z+2*m*h+zp)
    if mode == "singleB":
        return -1.0/(z-2*m*h+zp)+1.0/(z+2*m*h-zp)
    return (1.0/(z-2*m*h-zp)-1.0/(z-2*m*h+zp)
           +1.0/(z+2*m*h-zp)-1.0/(z+2*m*h+zp))
K0 = np.zeros((n, n)); Km = np.zeros((max(M,1), n, n))
for i in range(n):
    zi = zc[i]
    K0[i, :] = (G/(2*np.pi))*(F0(zi, ed[1:])-F0(zi, ed[:-1]))
    for m in range(1, M+1):
        Km[m-1, i, :] = (G/(2*np.pi))*(Fm(zi, ed[1:], m)-Fm(zi, ed[:-1], m))
if K0[0, 0] > 0:          # orientation: self-stiffness must be negative
    K0 = -K0; Km = -Km
Ksum = Km[:M].sum(axis=0) if M else np.zeros((n, n))

Wg = np.real(lambertw(np.exp(-Vpl/Vc)))
mu_ss = f0 + amb*np.log(Vpl/V0) + b*Wg
tau_i0 = sig*mu_ss + eta*Vpl
lth0 = np.log(dc*np.exp(Wg)/Vpl)

def vsolve(tau, lth, Vg):
    Phi = f0 + b*(np.log(V0/dc) + lth)
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

# state layout: [lth (n), D (n), chains (M*n if chained), S_e (n, kato3 only)]
chained = mode in ("exact", "kato3", "single", "singleB")
nS = 2*n + (M*n if chained else 0) + (n if mode == "kato3" else 0)
mm = np.arange(1, M+1)

def rhs(t, y, Vg):
    lth = y[:n]; D = y[n:2*n]
    tau = tau_i0 + K0@D
    if M:
        if chained:
            S = y[2*n:2*n+M*n].reshape(M, n)
            tau = tau + np.einsum('mij,mj->i', Km, S)
            if mode == "kato3":
                tau = tau + Ksum@y[2*n+M*n:]
        elif mode == "static":
            tau = tau + Ksum@D
        elif mode == "current":
            w = gammainc(mm, max(t, 0.0)/tr)
            tau = tau + np.tensordot(w, Km, axes=(0, 0))@D
    V = vsolve(tau, lth, Vg)
    th = np.exp(lth)
    dlth = np.exp(-V/Vc)/th - (V/dc)*np.log(np.maximum(V*th/dc, 1e-300))
    q = V - Vpl
    dy = np.empty_like(y)
    dy[:n] = dlth; dy[n:2*n] = q
    if M and chained:
        S = y[2*n:2*n+M*n].reshape(M, n)
        dS = np.empty((M, n))
        dS[0] = (D - S[0])/tr
        for m in range(1, M):
            dS[m] = (S[m-1] - S[m])/tr
        dy[2*n:2*n+M*n] = dS.ravel()
        if mode == "kato3":
            dy[2*n+M*n:] = q - y[2*n+M*n:]/tr
    return dy, V

# IC: steady state at Vpl; nucleation kick 4-6 km at same tau (theta lowered)
y0 = np.zeros(nS)
y0[:n] = lth0
kick = (zc > 4e3) & (zc < 6e3)
for j in np.where(kick)[0]:
    g = lambda l: (a[j]*sig[j]*np.arcsinh(3e-2/(2*V0)*np.exp((f0+b[j]*(np.log(V0/dc)+l))/a[j]))
                   + eta*3e-2 - tau_i0[j])
    y0[j] = brentq(g, lth0-60, lth0+5)

# probe cells for postseismic slip-rate histories (Kato Fig. 4 depths)
i15 = int(np.argmin(np.abs(zc-15e3))); i18 = int(np.argmin(np.abs(zc-18e3)))
pf = open(probe, "w") if probe else None
if pf: pf.write("# t_yr vmax v(%.2fkm) v(%.2fkm)\n" % (zc[i15]/1e3, zc[i18]/1e3))
t_lastw = -1.0

# Bogacki-Shampine 3(2), PI controller, like rsf_solve 3bs
rtol, atol_lth, atol_D = 1e-6, 1e-8, 1e-4
t, y = 0.0, y0.copy()
Vg = np.full(n, Vpl); Vg[kick] = 3e-2
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
        if (not slipping) and vmx > 1e-2:
            slipping = True; onsets.append(t/yr)
        elif slipping and vmx < 0.5e-2:
            slipping = False
        if pf and (t/yr - t_lastw) > 2e-3:
            pf.write("%.6f %.4e %.4e %.4e\n" % (t/yr, vmx, V1[i15], V1[i18]))
            t_lastw = t/yr
    dt = dt*min(4.0, max(0.2, 0.9*err**(-1.0/3.0)))
    dt = min(dt, tend_s-t+1.0) if t < tend_s else dt

if pf: pf.close()
onsets = np.array(onsets)
print(f"M={M} tr={tr_yr} n={n} dc={dc} mode={mode}: {len(onsets)} onsets, "
      f"{naccept}/{nstep} steps")
print(" onsets [yr]  :", " ".join(f"{x:.3f}" for x in onsets))
if len(onsets) > 2:
    print(" intervals[yr]:", " ".join(f"{d:.3f}" for d in np.diff(onsets)))
