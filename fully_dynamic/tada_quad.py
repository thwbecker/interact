"""
3-D full-space stress response to a constant slip rate on a quadrantal fault
(xi1 > 0, xi2 > 0, xi3 = 0) starting at t = 0, after Tada (2006, GJI 164,
653), Appendix A (eqs A12-A27), which restates Tada (2005):

    sigma_ij(x, t) = -mu/(2 cT) [D1 L_ij/1 + D2 L_ij/2 + D3 L_ij/3]

for slip u_k = D_k t H(t) H(xi1) H(xi2). L is dimensionless; at early time
inside the quadrant L_31/1 = 1 (radiation damping), at late time it grows
linearly with t (static stress of the growing slip).

Only shear slip in the plane is needed here (k = 1, 2); k = 2 follows from
k = 1 by the symmetry x1 <-> x2 (indices 1 <-> 2 swapped).
Vectorised with numpy; x1, x2 must not be exactly 0 (receivers at element
centres, quadrant corners at element edges).
"""
import numpy as np

def _H(x):
    return (x > 0).astype(float)

def _sgn(x):
    return np.sign(x)

def _sq(a):
    return np.sqrt(np.maximum(a, 0.0))

def g1(x1, x2, x3, t, c):
    chi1 = np.sqrt(x1**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
    A = (2*_H(x2)*_H(t - chi1/c) - _sgn(x2)*_H(t - R/c))
    return _sq(t**2 - chi1**2/c**2)**3 * A + (x2*(2*R**2 + chi1**2)/(2*R**3)*t**3 - 3*x2*chi1**2/(2*c**2*R)*t) * _H(t - R/c)

def g2(x1, x2, x3, t, c):
    chi1 = np.sqrt(x1**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
    A = (2*_H(x2)*_H(t - chi1/c) - _sgn(x2)*_H(t - R/c))
    return _sq(t**2 - chi1**2/c**2) * A + x2/R*t*_H(t - R/c)

def g3(x1, x2, x3, t, c):
    chi1 = np.sqrt(x1**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
    A = (2*_H(x2)*_H(t - chi1/c) - _sgn(x2)*_H(t - R/c))
    return np.arctan(_sq(t**2 - chi1**2/c**2)/(x1/c)) * A + np.arctan(x2/x1)*_H(t - R/c)

def g4(R, t, c):
    return (t**3 - R**2/c**2*t)*_H(t - R/c)

def g5(R, t, c):
    return (t - R/c)*_H(t - R/c)

def L_1(x1, x2, x3, t, cL, cT):
    """L_ij/1 for (ij) = 11, 12, 22, 31, 32, 33 (eqs A18-A23)"""
    p = cT/cL
    chi1 = np.sqrt(x1**2 + x3**2); chi2 = np.sqrt(x2**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
    dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT)
    dg4 = g4(R, t, cL) - g4(R, t, cT)
    g2L = g2(x1, x2, x3, t, cL); g2T = g2(x1, x2, x3, t, cT)
    g2sT = g2(x2, x1, x3, t, cT)
    pi = np.pi
    L11 = (-2*cT**3*x3*(3*x1**2 - x3**2)/(3*pi*chi1**6)*dg1 - cT**3*x1**2*x2*x3/(pi*R**5*chi1**2)*dg4
           - cT*x3/(pi*chi1**2)*((1 - 2*p**2 + 2*p**2*x1**2/chi1**2)*g2L + (1 - 2*x1**2/chi1**2)*g2T))
    L12 = cT**3*x1*x3/(pi*R**5)*dg4 - cT*x3/(2*pi*chi2**2)*g2sT
    L22 = cT**3*x2*x3/(pi*R**5)*dg4 - (1 - 2*p**2)*cT*x3/(pi*chi1**2)*g2L
    L31 = (_H(x1)*_H(x2)*_H(t - np.abs(x3)/cT) - (g3(x1, x2, x3, t, cT) + g3(x2, x1, x3, t, cT))/(2*pi)
           - 2*cT**3*x1*(3*x3**2 - x1**2)/(3*pi*chi1**6)*dg1 - cT**3*x1*x2*x3**2/(pi*R**5*chi1**2)*dg4
           - 2*cT*x1*x3**2/(pi*chi1**4)*(p**2*g2L - g2T) + cT*x2/(2*pi*chi2**2)*g2sT)
    L32 = cT**3/(pi*R**5)*(x3**2 - R**2/3)*dg4 + cT/(pi*R)*(2*p**2/3*g5(R, t, cL) - g5(R, t, cT)/6)
    L33 = (2*cT**3*x3*(3*x1**2 - x3**2)/(3*pi*chi1**6)*dg1 - cT**3*x2*x3**3/(pi*R**5*chi1**2)*dg4
           - cT*x3/(pi*chi1**2)*((1 - 2*p**2 + 2*p**2*x3**2/chi1**2)*g2L + (1 - 2*x3**2/chi1**2)*g2T))
    return dict(s11=L11, s12=L12, s22=L22, s31=L31, s32=L32, s33=L33)

def L_3(x1, x2, x3, t, cL, cT):
    """L_ij/3 (opening, eqs A24-A27; L_22/3 and L_32/3 by the x1 <-> x2 symmetry)"""
    p = cT/cL; pi = np.pi
    def l11(x1, x2):
        chi1 = np.sqrt(x1**2 + x3**2); chi2 = np.sqrt(x2**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
        dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT); dg4 = g4(R, t, cL) - g4(R, t, cT)
        return ((1 - 2*p**2)/p*_H(x1)*_H(x2)*_H(t - np.abs(x3)/cL)
                - (1 - 2*p**2)/(2*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL))
                - 2*cT**3*x1*(3*x3**2 - x1**2)/(3*pi*chi1**6)*dg1 - cT**3*x1*x2*x3**2/(pi*R**5*chi1**2)*dg4
                - 2*cT*x1*x3**2/(pi*chi1**4)*(p**2*g2(x1, x2, x3, t, cL) - g2(x1, x2, x3, t, cT))
                + (1 - 2*p**2)*cT*x2/(pi*chi2**2)*g2(x2, x1, x3, t, cL))
    chi1 = np.sqrt(x1**2 + x3**2); chi2 = np.sqrt(x2**2 + x3**2); R = np.sqrt(x1**2 + x2**2 + x3**2)
    dg4 = g4(R, t, cL) - g4(R, t, cT)
    L12 = cT**3/(pi*R**5)*(x3**2 - R**2/3)*dg4 + cT/(pi*R)*((1 - 4*p**2/3)*g5(R, t, cL) - 2.0/3.0*g5(R, t, cT))
    L33 = (1/p*_H(x1)*_H(x2)*_H(t - np.abs(x3)/cL) - 1/(2*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL))
           + 2*cT**3*x1*(3*x3**2 - x1**2)/(3*pi*chi1**6)*(g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT))
           + 2*cT**3*x2*(3*x3**2 - x2**2)/(3*pi*chi2**6)*(g1(x2, x1, x3, t, cL) - g1(x2, x1, x3, t, cT))
           + cT**3*x1*x2*x3**2/(pi*R**5)*(1/chi1**2 + 1/chi2**2)*dg4
           + cT*x1/(pi*chi1**2)*((2*(1 - p**2) + 2*p**2*x3**2/chi1**2)*g2(x1, x2, x3, t, cL) - 2*x3**2/chi1**2*g2(x1, x2, x3, t, cT))
           + cT*x2/(pi*chi2**2)*((2*(1 - p**2) + 2*p**2*x3**2/chi2**2)*g2(x2, x1, x3, t, cL) - 2*x3**2/chi2**2*g2(x2, x1, x3, t, cT)))
    L31 = L_1(x1, x2, x3, t, cL, cT)["s33"]
    L32 = L_1(x2, x1, x3, t, cL, cT)["s33"]
    return dict(s11=l11(x1, x2), s22=l11(x2, x1), s12=L12, s31=L31, s32=L32, s33=L33)

def L_2(x1, x2, x3, t, cL, cT):
    """L_ij/2 by the x1 <-> x2 symmetry: L_ij/2(x1, x2) = L_i'j'/1(x2, x1), 1' = 2, 2' = 1"""
    L = L_1(x2, x1, x3, t, cL, cT)
    return dict(s11=L["s22"], s12=L["s12"], s22=L["s11"], s31=L["s32"], s32=L["s31"], s33=L["s33"])

def rect_ramp_stress(x1, x2, x3, t, a1, a2, cL, cT, rho, comp=1):
    """stress at receiver (x1, x2, x3) from slip rate 1 in direction comp (1 or 2; 3 = opening)
    on the rectangle |xi1| < a1, |xi2| < a2 (xi3 = 0) since t = 0: four quadrants"""
    mu = rho*cT**2
    Lf = {1: L_1, 2: L_2, 3: L_3}[comp]
    out = None
    for (s1, c1) in ((+1, -a1), (-1, a1)):
        for (s2, c2) in ((+1, -a2), (-1, a2)):
            L = Lf(x1 - c1, x2 - c2, x3, t, cL, cT)
            sgn = s1*s2
            if out is None: out = {k: sgn*v for k, v in L.items()}
            else:
                for k in out: out[k] = out[k] + sgn*L[k]
    return {k: -mu/(2*cT)*v for k, v in out.items()}
