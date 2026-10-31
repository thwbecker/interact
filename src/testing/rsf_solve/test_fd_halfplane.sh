#!/bin/bash
# rsf_solve -fd_mode 1 on a 2-D in-plane fault in a half plane (fd_plan.md M7 step 2c): a 30
# degree thrust of 4 km with elements of 200 m, its top 1 km deep; the virtual free surface is
# built from the geometry and its static limit (the fault operator with the surface in static
# equilibrium) must agree with interact's exact half-plane Is to 5e-3. Run from the repository
# root after make bin/rsf_solve; about a minute.
d=tmp_fd_halfplane
rm -rf $d
mkdir -p $d
cd $d || exit 1
python3 - << 'PYEOF'
import numpy as np
top, Lf, dx = 1000.0, 4e3, 200.0
dip = np.radians(30.0)
n = int(Lf/dx)
s = np.array([np.cos(dip), -np.sin(dip)])
c = np.array([0.0, -top])[None, :] + ((np.arange(n) + 0.5)*dx)[:, None]*s[None, :]
strike = 90.0 - np.degrees(np.arctan2(s[1], s[0]))
a = np.full(n, 0.010)
b = np.full(n, 0.015)
psi = 0.6 + b*np.log(1e-6/1e-9)
tau0 = a*50e6*np.arcsinh(1e-9/(2e-6)*np.exp(psi/a))
with open("geom.in", "w") as g, open("rsf.dat", "w") as r, open("ic.in", "w") as ic:
    for i in range(n):
        g.write(f"{c[i, 0]:.9e} {c[i, 1]:.9e} 0 {strike:.9f} 90 {dx/2:.6e} 0 0\n")
        r.write(f"{a[i]:.6f} {b[i]:.6f}\n")
        ic.write(f"{tau0[i]:.8e} 1e-9\n")
PYEOF
../bin/rsf_solve -geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -dc 0.016 -sigma_init 50e6 -v0 1e-6 -f0 0.6 \
  -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 -rtol 1e-6 -stop_time_yr 1e-6 -fd_mode 1 > run.log 2>&1
dev=`grep 'static gate: rel dev of Ks' run.log | sed 's/.*from Is \([0-9.e+-]*\).*/\1/'`
nsurf=`grep 'virtual free surface:' run.log | sed 's/.*m, \([0-9]*\) elements.*/\1/'`
if [ -n "$dev" ] && awk -v d=$dev 'BEGIN{exit !(d < 5e-3)}'; then
  echo "test_fd_halfplane: PASS (static gate $dev against the exact half plane, $nsurf surface elements)"
  exit 0
fi
echo "test_fd_halfplane: FAIL (static gate '$dev')"
exit 1
