#!/usr/bin/env python3
"""
Analytical check for the 2D half-plane displacement-discontinuity
implementation (Crouch and Starfield style elements) in interact.

Reference solution: a uniform-slip segment in the lower half-plane
(free surface on y = 0) is exactly a dislocation dipole, one edge
dislocation at each fault end. The half-plane edge dislocation is
built here from Muskhelishvili potentials; the full-space part is the
standard phi = g ln(z-z0), psi = conj(g) ln(z-z0) - g conj(z0)/(z-z0)
with g = G b / (i pi (1+kappa)), kappa = 3-4 nu, and the image
coefficients were obtained by solving the traction-free boundary
condition symbolically (sympy), not copied from a reference, so this
check is independent of the Crouch and Starfield transcription. The
implementation was validated against: a traction-free surface to
machine precision, the Burgers circuit recovering the slip vector,
the deep-dislocation limit matching the textbook full-space field,
and translation invariance. The displacement branch cut is placed on
the fault by evaluating the endpoint log difference as one ratio log.

Usage: analytic_compare.py <dip_deg> <half_length> [G] [nu] [slip]
reads geom-implied endpoints (fault from the origin, dipping
clockwise from horizontal as in run_test), disp.out and stress.out
from the current directory (uniform slip run), and reports relative
misfits. Sign convention: interact slip +1 in the strike slot of the
2D setup corresponds to slip -1 along the (zB - zA) tangent here.
G defaults to 1e4 (the compiled SHEAR_MODULUS of the standard 2D
build); adjust if built differently.
"""
import sys
import numpy as np

def segment_halfplane(zz, zA, zB, s, G, nu):
    kap = 3.0 - 4.0*nu
    t = (zB - zA)/abs(zB - zA)
    b = s*t
    g = G*b/(1j*np.pi*(1.0+kap))
    gr, gi = g.real, g.imag
    gc = np.conj(g)
    def coeffs(z0):
        x0, y0 = z0.real, z0.imag
        a1 = -gr - 1j*gi
        a2 = 2*gi*y0 + 1j*2*gr*y0
        b1 = -gr + 1j*gi
        b2 = (3*gi*y0 + gr*x0) + 1j*(gi*x0 + gr*y0)
        b3 = 2*y0*(gi*x0 + gr*y0) + 1j*2*y0*(gr*x0 - gi*y0)
        return a1, a2, b1, b2, b3
    zAc, zBc = np.conj(zA), np.conj(zB)
    aA1,aA2,bA1,bA2,bA3 = coeffs(zA)
    aB1,aB2,bB1,bB2,bB3 = coeffs(zB)
    d0B, d0A = zz - zB, zz - zA
    dcB, dcA = zz - zBc, zz - zAc
    L0 = np.log(d0B/d0A)	# branch cut on the fault segment
    Lc = np.log(dcB/dcA)	# cut on the mirror segment, outside
    phi  = g*L0 + aB1*Lc + aB2/dcB - aA2/dcA
    dphi = g/d0B - g/d0A + aB1/dcB - aA1/dcA - aB2/dcB**2 + aA2/dcA**2
    ddphi= -g/d0B**2 + g/d0A**2 - aB1/dcB**2 + aA1/dcA**2 \
	+ 2*aB2/dcB**3 - 2*aA2/dcA**3
    psi  = gc*L0 - g*zBc/d0B + g*zAc/d0A + bB1*Lc \
	+ bB2/dcB - bA2/dcA + bB3/dcB**2 - bA3/dcA**2
    dpsi = gc/d0B - gc/d0A + g*zBc/d0B**2 - g*zAc/d0A**2 \
	+ bB1/dcB - bA1/dcA - bB2/dcB**2 + bA2/dcA**2 \
	- 2*bB3/dcB**3 + 2*bA3/dcA**3
    u = (kap*phi - zz*np.conj(dphi) - np.conj(psi))/(2*G)
    s_sum = 4.0*np.real(dphi)
    s_dif = 2.0*(np.conj(zz)*ddphi + dpsi)
    sxx = 0.5*(s_sum - np.real(s_dif))
    syy = 0.5*(s_sum + np.real(s_dif))
    sxy = 0.5*np.imag(s_dif)
    return u, sxx, syy, sxy

def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        sys.exit(1)
    dip = float(sys.argv[1])
    hl = float(sys.argv[2])
    G = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e4
    nu = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
    slip = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    d = np.loadtxt('disp.out'); s = np.loadtxt('stress.out')
    zz = d[:,0] + 1j*d[:,1]
    zA = 0.0 + 0.0j
    zB = 2.0*hl*(np.cos(np.radians(dip)) - 1j*np.sin(np.radians(dip)))
    # exclude points too close to the fault plane and its two tips,
    # where the discrete elements and the ideal cut legitimately differ
    t = (zB-zA)/abs(zB-zA)
    proj = np.clip(((zz-zA)*np.conj(t)).real, 0.0, abs(zB-zA))
    dist = np.abs(zz - (zA + proj*t))
    seg = abs(zB - zA)
    mask = (dist > 0.03*seg) & (np.abs(zz-zA) > 0.06*seg) & (np.abs(zz-zB) > 0.06*seg)
    u, sxx, syy, sxy = segment_halfplane(zz, zA, zB, -slip, G, nu)
    ok = True
    tol = 1.0e-5
    for nm, a, b in (("sxx", s[:,3], sxx), ("sxy", s[:,4], sxy), ("syy", s[:,6], syy),
                     ("ux", d[:,3], u.real), ("uy", d[:,4], u.imag)):
        aa, bb = a[mask], b[mask]
        rel = np.linalg.norm(aa-bb)/np.linalg.norm(bb)
        print(f"analytic_compare: {nm:3s} relative L2 misfit {rel:.3e}")
        if not (rel < tol):
            ok = False
    print(f"analytic_compare: {'PASS' if ok else 'FAIL'} (tolerance {tol:g}; "
          f"{np.sum(mask)} of {len(zz)} grid points compared)")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
