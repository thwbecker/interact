"""Memory estimates (GB) for the 3-D virtual-surface tables in the design of
fd_plan.md ("Free surface in 3-D: design"). All couplings are Toeplitz along
strike (FFT, half spectrum, single precision complex = 8 bytes) and dense
across strike / down dip. L: lags, the S-wave crossing time of the farthest
pair (W vanishes exactly afterwards in 3-D) at dt = 0.5 h_f / cS."""
import numpy as np

def est(name, Lx, W, dip, hf, hs_fine, w_fine, ext, h_out, opening_only=False):
    nx = int(Lx/hf); nzf = int(W/hf)
    # surface rows across strike: a fine strip of width w_fine around the trace, coarse rows
    # out to +-ext beyond the fault's surface projection
    proj = W*np.cos(np.radians(dip))
    n_fine = int(w_fine/hs_fine); n_out = int(2*ext/h_out) + int(max(proj - w_fine, 0.0)/h_out)
    ny = n_fine + n_out
    Mx = int(2**np.ceil(np.log2(2*nx + 2*ext/hf)))          # along strike, including the surface margin
    Mh = Mx//2 + 1
    rmax = np.hypot(Lx + 2*ext, proj + 2*ext); L = int(rmax/(0.5*hf))
    cs = 1 if opening_only else 3                 # surface unknowns per element
    Mz = int(2**np.ceil(np.log2(2*nzf)))
    ff = Mh*Mz*2                                  # fault -> fault: 2-D Toeplitz (shear and normal traction)
    fs = Mh*nzf*ny*cs                             # fault slip -> surface tractions
    sf = Mh*ny*nzf*cs*2                           # surface sources -> fault shear and normal
    ss = Mh*ny*ny*(cs if opening_only else 5)     # coplanar: shear-shear 2x2 (3 independent) + opening-normal
    gb = (ff + fs + sf + ss)*L*8/1e9
    print(f"{name}: fault {nx} x {nzf} of {hf/1e3:g} km, surface {ny} rows ({n_fine} at {hs_fine:g} m), "
          f"FFT {Mx} along strike, {L} lags: tables about {gb:.1f} GB")
    return gb

if __name__ == "__main__":
    est("BP5-FD, mirror + opening-only surface layer", 100e3, 40e3, 90.0, 1e3, 250.0, 4e3, 40e3, 2e3, opening_only=True)
    est("BP5-FD, full virtual surface (no mirror)", 100e3, 40e3, 90.0, 1e3, 250.0, 4e3, 40e3, 2e3)
    est("dipping thrust 60 x 30 km, 30 deg, 1 km", 60e3, 30e3, 30.0, 1e3, 250.0, 6e3, 30e3, 2e3)
    est("dipping thrust 60 x 30 km, 30 deg, 500 m", 60e3, 30e3, 30.0, 500.0, 125.0, 6e3, 30e3, 1e3)
    est("small dipping thrust 20 x 10 km, 30 deg, 250 m", 20e3, 10e3, 30.0, 250.0, 62.5, 2e3, 10e3, 500.0)
