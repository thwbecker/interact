#!/usr/bin/env python3
"""
Convert a two-column depth[km] a-b table (zab_TzT.dat, a-b derived from a
subduction geotherm T(z) and laboratory a-b(T)) into the three-column
"depth_km a b" control-point format that rsf_input_gen_cascadia.py reads.

Only a-b controls frictional stability, but a and b are needed
separately: a sets the direct velocity effect and b enters the
nucleation length L_b = G D_c / (b sigma).  The table fixes only the
difference, so b is held at a constant lab value here and a = b + (a-b)
follows.  Holding b rather than a keeps L_b constant with depth, which
is the choice that makes the resolution census meaningful across the
whole fault.

The input table starts below the seafloor (at 3.1 km on the supplied
file) while the mesh has elements shallower than that.  Rather than let
the interpolation flat-extrapolate an arbitrary value up to y = 0, the
shallow end is continued explicitly to surface_amb over surface_z_km,
which keeps the updip region velocity strengthening as it should be.

Parameters are plain assignments below; no environment variables.
"""
import numpy as np

infile = "zab_TzT.dat"
outfile = "profile_cascadia_TzT.dat"
b_const = 0.015       # constant b.  Matches the Uphoff runs, so D_c,
                      # sigma and f0 carry over and L_b is unchanged
surface_z_km = 0.0    # continue the table up to this depth ...
surface_amb = None    # ... at this a-b.  None: hold the table's own
                      # shallowest value, which is already velocity
                      # strengthening on this file
decimate_km = 0.0     # thin the control points to roughly this spacing
                      # (0.0: keep every row).  The generator interpolates
                      # linearly either way, so this only shortens the file
max_z_km = 0.0        # drop control points below this depth (0.0: keep
                      # all).  The tail past the mesh bottom costs
                      # nothing and documents where a-b is heading

d = np.loadtxt(infile)
z, amb = d[:, 0], d[:, 1]
o = np.argsort(z)
z, amb = z[o], amb[o]

if z[0] > surface_z_km:
    z = np.hstack([surface_z_km, z])
    amb = np.hstack([amb[0] if surface_amb is None else surface_amb, amb])
if max_z_km > 0:
    k = z <= max_z_km
    z, amb = z[k], amb[k]
if decimate_km > 0:
    keep = [0]
    for i in range(1, len(z) - 1):
        if z[i] - z[keep[-1]] >= decimate_km:
            keep.append(i)
    keep.append(len(z) - 1)
    z, amb = z[keep], amb[keep]

a = b_const + amb
sgn = np.sign(amb)
cross = [z[i] + (z[i+1] - z[i])*(-amb[i])/(amb[i+1] - amb[i])
         for i in np.where(np.diff(sgn) != 0)[0]]

with open(outfile, "w") as f:
    f.write("# Depth-dependent rate and state parameters for the megathrust,\n")
    f.write("# converted from %s by ab2profile.py.  That table gives a-b(z)\n" % infile)
    f.write("# from a subduction geotherm and laboratory a-b(T); here b is held\n")
    f.write("# constant at %.4f and a = b + (a-b).\n#\n" % b_const)
    f.write("# a-b spans %+.5f to %+.5f over %.2f to %.2f km depth.\n"
            % (amb.min(), amb.max(), z[0], z[-1]))
    f.write("# Velocity weakening between %s km depth;\n"
            % " and ".join("%.2f" % c for c in cross))
    f.write("# most weakening a-b = %+.5f at %.2f km.\n" % (amb.min(), z[int(np.argmin(amb))]))
    if d[:, 0].min() > surface_z_km:
        f.write("# The table began at %.2f km; it is continued to %.2f km at a-b\n"
                "# = %+.5f so the shallowest elements are not flat-extrapolated.\n"
                % (d[:, 0].min(), surface_z_km, amb[0]))
    f.write("#\n# Linear interpolation between control points; end values extend beyond.\n")
    f.write("# depth_km        a         b\n")
    for zz, aa in zip(z, a):
        f.write("%10.4f %10.6f %10.6f\n" % (zz, aa, b_const))

print("wrote %s: %d control points, %.2f to %.2f km" % (outfile, len(z), z[0], z[-1]))
print("  a-b %+.5f to %+.5f, b held at %.4f so a spans %.6f to %.6f"
      % (amb.min(), amb.max(), b_const, a.min(), a.max()))
print("  velocity weakening between " + " and ".join("%.2f" % c for c in cross) + " km depth")
