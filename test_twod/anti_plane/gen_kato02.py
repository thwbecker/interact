#!/usr/bin/env python3
"""
gen_kato02.py: generate the Kato (EPS 2002) replica setup for rsf_solve.

Model (Kato, Earth Planets Space 54, 1077-1083, 2002): infinitely long
vertical strike-slip fault (2-D antiplane) cutting an elastic layer of
thickness h over a Maxwell half-space, backslip loading at
Vpl = 35 mm/yr, composite Kato-Tullis friction (state_law 5,
Vc = 0.01 um/s = 1e-2 * v0 with -v0 1e-6), sigma_eff = (rho - rho_w) g z
with rho = 2.8, rho_w = 1.0 g/cm^3, G = 30 GPa, c_s = 3.27 km/s,
L = 5 cm.  His tr = 2 eta/G is TWICE our -ve_tmaxwell_yr (equal
rigidities).

The a(z) and a-b(z) profiles are DIGITIZED from his Fig. 2 (pixel
measurement of the PDF, tick-calibrated; anti-aliasing worth about
+-0.0002): a-b = +0.004 at the surface, zero at 2.5 km, -0.004 flat
from 5 to 10 km, rising through zero at 10.55 km to +0.0121 at 12 km;
a = 0.0121 constant above 12 km; below 12 km a-b = a (i.e. b = 0) with
a rising linearly to 0.0345 at 20 km.  Since the psi-based engine
divides by b, b is floored at B_FLOOR = 0.001 below 12 km (the state
term is then b ln(theta V/L) <= 1e-3, dynamically negligible against
a = 0.012-0.035 there).

usage: gen_kato02.py [ds_km] [h_km] [prefix]
       defaults      0.1      20     kato02
writes <prefix>_{geom.in,rsf.dat,ic.in,sigma.in}
group 0 = a-b < 0 cells (so rsf_catalog.g000.dat mean slip is Kato's us)
"""
import sys
import numpy as np

ds     = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1     # cell [km]
h_km   = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0    # layer/fault depth [km]
prefix = sys.argv[3]        if len(sys.argv) > 3 else "kato02"

# Kato 2002 material and friction constants
G    = 30.0e9                    # [Pa]
cs   = 3270.0                    # [m/s]
eta  = G/(2.0*cs)                # radiation damping [Pa s/m]
Vpl  = 35.0e-3/3.15576e7         # 35 mm/yr [m/s] = 1.109e-9
Vc   = 1.0e-8                    # composite-law cutoff, 0.01 um/s [m/s]
V0   = 1.0e-6                    # reference velocity V* (arbitrary in Kato)
f0   = 0.6                       # reference friction at V*
L    = 0.05                      # [m]
rho, rhow, grav = 2800.0, 1000.0, 9.8
B_FLOOR = 0.001

# Fig. 2, digitized (see header); piecewise linear in depth [km]
ZMAX   = 20.0
zk_amb = np.array([0.0,   2.5, 5.0,    10.0,   12.0,   ZMAX])
va_amb = np.array([0.004, 0.0, -0.004, -0.004, 0.0121, 0.0345])
def prof_a(z):
    return np.where(z <= 12.0, 0.0121, 0.0121 + (0.0345-0.0121)*(z-12.0)/8.0)

def W(x):                        # Lambert W on (0, 1], Newton
    w = 0.5*np.ones_like(x)
    for _ in range(50):
        w = w - (w*np.exp(w)-x)/(np.exp(w)*(1.0+w))
    return w

nz = int(round(h_km/ds))
half = ds*1e3/2.0
zn1, zn2, Vnuc = 4.0, 6.0, 3e-2   # nucleation seed in the VW zone
fg = open(f"{prefix}_geom.in","w")
fr = open(f"{prefix}_rsf.dat","w")
fi = open(f"{prefix}_ic.in","w")
fs = open(f"{prefix}_sigma.in","w")
for j in range(nz):
    zmid = (j+0.5)*ds                                  # depth [km]
    amb  = float(np.interp(zmid, zk_amb, va_amb))      # a-b
    a    = float(prof_a(np.array(zmid)))
    b    = max(a - amb, B_FLOOR)
    amb  = a - b                                       # keep a-b consistent
    sig  = (rho-rhow)*grav*zmid*1e3                    # [Pa]
    v    = Vnuc if (zn1 <= zmid <= zn2) else Vpl
    # composite-law steady state at Vpl:
    # mu_ss = f0 + (a-b) ln(Vpl/V0) + b W(exp(-Vpl/Vc))
    mu_ss = f0 + amb*np.log(Vpl/V0) + b*float(W(np.array([np.exp(-Vpl/Vc)]))[0])
    tau  = sig*mu_ss + eta*v
    grp  = 0 if amb < 0.0 else 1
    fg.write(f"0 {-zmid*1e3:.6e} 0 0 -90 {half:.6e} 0 {grp}\n")
    fr.write(f"{a:.6f} {b:.6f}\n")
    fi.write(f"{tau:.8e} {v:.6e}\n")
    fs.write(f"{sig:.8e}\n")
fg.close(); fr.close(); fi.close(); fs.close()

# resolution check: cell size vs 0.093 h* (Rice 1993),
# h* = (2/pi) G L / ((b-a) sigma), over the velocity-weakening zone
zz = (np.arange(nz)+0.5)*ds
amb = np.interp(zz, zk_amb, va_amb); sg = (rho-rhow)*grav*zz*1e3
vw = amb < -1e-6
hstar = np.full(nz, np.inf)
hstar[vw] = (2.0/np.pi)*G*L/(-amb[vw]*sg[vw])
print(f"gen_kato02: {nz} cells, ds {ds} km, h {h_km} km; "
      f"min 0.093 h* = {0.093*hstar.min():.1f} m vs cell {ds*1e3:.0f} m; "
      f"Vpl {Vpl:.4e} m/s, eta {eta:.4e}")
