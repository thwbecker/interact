"""Checks of the 2-D in-plane opening kernels (inplane_open.py):
- early time at the segment centre: the plane P wave, sxx = -lambda/(2 alpha),
  syy = -rho alpha/2 per unit opening rate, sxy = 0;
- static limit (late-time slope) against Crouch and Starfield's constant
  displacement discontinuity, whose convention D = u(-) - u(+) is opposite to
  ours (the same overall sign relation holds for the existing shear kernel)."""
import numpy as np, sympy as sp
import inplane as ip, inplane_open as io
AL, BE, RHO = 6000.0, 3464.0, 2670.0
MU = RHO*BE**2; LAM = RHO*(AL**2 - 2*BE**2); NU = LAM/(2*(LAM + MU)); A = 50.0

def _cs():
    x, y, a, G, nu = sp.symbols("x y a G nu", real=True)
    f = -1/(4*sp.pi*(1 - nu))*(y*(sp.atan(y/(x - a)) - sp.atan(y/(x + a))) - (x - a)*sp.log(sp.sqrt((x - a)**2 + y**2))
                              + (x + a)*sp.log(sp.sqrt((x + a)**2 + y**2)))
    d = lambda *v: sp.diff(f, *v)
    opening = [2*G*(d(y, y) + y*d(y, y, y)), 2*G*(-y*d(x, y, y)), 2*G*(d(y, y) - y*d(y, y, y))]
    shear = [2*G*(2*d(x, y) + y*d(x, y, y)), 2*G*(d(y, y) + y*d(y, y, y)), 2*G*(-y*d(x, y, y))]
    lf = lambda e: sp.lambdify((x, y, a, G, nu), e)
    return [lf(e) for e in opening], [lf(e) for e in shear]

def test_early_time_plane_p_wave():
    s = io.segment_ramp_stress_open(0.0, 0.0, 0.001, A, AL, BE, RHO)
    assert abs(float(s[0]) + LAM/(2*AL)) < 1e-9*LAM/AL
    assert abs(float(s[1])) < 1e-6
    assert abs(float(s[2]) + RHO*AL/2) < 1e-9*RHO*AL

def test_static_limit_against_crouch_starfield():
    op, sh = _cs()
    for xr, yr in ((150., 0.), (80., 60.), (300., -200.), (10., 500.)):
        T1, T2 = 20.0, 40.0
        s1 = io.segment_ramp_stress_open(xr, yr, T1, A, AL, BE, RHO); s2 = io.segment_ramp_stress_open(xr, yr, T2, A, AL, BE, RHO)
        ref = max(abs(op[k](xr, yr, A, MU, NU)) for k in range(3))
        for k in range(3):
            slope = (float(s2[k]) - float(s1[k]))/(T2 - T1)
            assert abs(slope + op[k](xr, yr, A, MU, NU)) < 2e-5*ref
        s = ip.static_segment_stress(xr, yr, A, MU, NU)
        refs = max(abs(sh[k](xr, yr, A, MU, NU)) for k in range(3))
        for k in range(3):
            assert abs(float(s[k]) + sh[k](xr, yr, A, MU, NU)) < 1e-9*refs

if __name__ == "__main__":
    test_early_time_plane_p_wave(); test_static_limit_against_crouch_starfield(); print("test_inplane_open: PASS")
