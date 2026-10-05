#!/bin/bash
#
# run_sp_test.sh: exact test of the layer code's time dependence and
# relaxation-time convention against the Savage and Prescott (1978)
# image series for an infinitely long vertical strike-slip fault
# cutting the whole elastic plate over a Maxwell half-space of the
# same rigidity.  The series uses tau_a = 2 eta/mu; Johnson's tR is
# that same quantity, so the comparison below passes with tR = tau_a
# and fails by about 50 percent at early times with tR = tau_a/2.
#
# usage: run_sp_test.sh [nk] [H_km] [tR_yr]
#        defaults        1000  30     1.0
#
# The half-space is represented by the viscoelastic layer (H1..H2)
# plus the half-space below with tR2 = 1.002 tR1 (equal relaxation
# times are ill conditioned, see README).  nk = 100 (the original) is
# 7 to 12 percent low at t = 50 tR; nk = 1000 is 2 to 4 percent low.
#
nk=${1-1000}
H=${2-30}
tR=${3-1.0}
tR2=`awk -v t="$tR" 'BEGIN{printf "%.6f", t*1.002}'`
H2=`awk -v h="$H" 'BEGIN{printf "%g", 2*h}'`
L=`awk -v h="$H" 'BEGIN{printf "%g", 130*h}'`
nl=`awk -v l="$L" 'BEGIN{printf "%d", l/5}'`
printf "12 0\n32 0\n62 0\n100 0\n" > sp_sta.txt
for t in 2 10 50;do
    tt=`awk -v t="$t" -v r="$tR" 'BEGIN{printf "%g", t*r}'`
    ./pom layer -m $L $H $H 90 0 0 0 1 0 0 -H1 $H -H2 $H2 -nu 0.25 -tR1 $tR -tR2 $tR2 \
	-t $tt -nl $nl -nw 8 -nk $nk -rg 0 -sta sp_sta.txt -o sp_t$t.txt || exit 1
done
python3 - "$H" "$tR" <<'PY'
import sys, math, numpy as np
H = float(sys.argv[1]); tR = float(sys.argv[2]); D = H
def pcdf(nm1, x):            # e^{-x} sum_{m<=nm1} x^m/m!, stable
    s, term = 0.0, 1.0
    for m in range(nm1 + 1):
        s += term; term *= x / (m + 1)
    return math.exp(-x) * s
def sp(x, t, ta, nmax=600):
    u = np.zeros_like(x)
    for n in range(1, nmax):
        Fn = 1.0 - pcdf(n - 1, t / ta)
        u += Fn * (np.arctan((2 * n * H + D) / x) - np.arctan((2 * n * H - D) / x))
    return u / np.pi
ok = True
for t in (2, 10, 50):
    d = np.loadtxt(f'sp_t{t}.txt'); x = d[:, 0]; u = np.abs(d[:, 3])
    a = sp(x, t * tR, tR); b = sp(x, t * tR, tR / 2)
    ea = 100 * (u / a - 1); eb = 100 * (u / b - 1)
    print(f't = {t:2d} tR: pom/SP(tau_a = tR) - 1 [%]: ' + ' '.join(f'{v:+6.1f}' for v in ea)
          + f'   | with tau_a = tR/2: ' + ' '.join(f'{v:+6.1f}' for v in eb))
    ok &= np.abs(ea).max() < 6.0
print('Savage-Prescott test', 'PASSED (tR = 2 eta/mu confirmed)' if ok else 'FAILED')
PY
