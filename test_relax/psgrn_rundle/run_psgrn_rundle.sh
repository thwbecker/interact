#!/bin/bash
#
# run_psgrn_rundle.sh: PSGRN/PSCMP (Wang et al., 2020 version) runs for
# the Rundle (1982) thrust configuration at several fault lengths, with
# and without gravity, as the 3-D reference for the 2-D codes
# (ve_thrust_relax.py, inplane2d.py) and for the finite-fault-length
# question (Rundle's 2L = 200 km vs the 2-D limit).
#
# Model: elastic plate H = 30 km (vp 5.2223, vs 3.0151 km/s, rho 3300,
# i.e. mu = lam = 30 GPa) over a Maxwell half-space (vp 4.8676, vs
# 2.8103, rho 3800, eta = 9.468e17 Pa s, tau_M = eta/mu = 1.0 yr).
# Fault: 30 deg thrust, surface breaking, 30 km down-dip (bottom at
# 15 km), 1 m slip, strike east, dipping south; receivers at
# mid-strike from 90 km north (footwall) to 160 km south (hanging
# wall), 2.5 km spacing.  PSGRN relaxes the shear modulus at constant
# bulk modulus (same as ve_thrust_relax.py without --lam-const).
#
# usage: run_psgrn_rundle.sh [bindir] [lengths_km] [gravity_list] [twin_yr] [nt]
#        defaults          ../psgrn_pscmp/bin "200 600 2000" "0 1"  90    128
#
#   bindir   : where psgrn2020 and pscmp2020 are (tools/psgrn/install_psgrn
#              puts them in <install_dir>/bin)
#   lengths  : PSCMP fault lengths; 200 = Rundle's 2L = 6.7 H,
#              2000 = 2-D limit (67 H)
#   gravity  : PSGRN gravity factor, 0 = off, 1 = on
#   twin, nt : Green's function time window [yr] and number of time
#              samples (uniform; 128 over 90 yr gives dt = 0.7 yr, fine
#              for reading 10 and 90 yr; the time series start at t = 0,
#              which is the coseismic field)
#
# Output: grn_g<g>/ (Green's functions, large, reusable), out_g<g>_L<L>/
# with U_down.dat, U_north.dat, U_east.dat (rows = times, columns =
# receivers), then psgrn_rundle_compare.py for the figures and tables.
# Run from an empty scratch directory with a few GB free; PSGRN takes
# tens of minutes per gravity setting, each PSCMP a few minutes.
#
bindir=${1-../../tools/psgrn/psgrn_pscmp/bin}
lengths=${2-"200 600 2000"}
gravs=${3-"0 1"}
twin=${4-90}
nt=${5-128}
bindir=`cd "$bindir" && pwd` || { echo "$0: no bindir $1"; exit 1; }
[ -x $bindir/psgrn2020 ] && [ -x $bindir/pscmp2020 ] || { echo "$0: psgrn2020/pscmp2020 not in $bindir"; exit 1; }
tdays=`awk -v y="$twin" 'BEGIN{printf "%.1f", y*365.25}'`
KM=111.195

# receivers: 101 points along longitude 0 from lat +90 km to -160 km
lat1=`awk -v k=$KM 'BEGIN{printf "%.6f", 90.0/k}'`
lat2=`awk -v k=$KM 'BEGIN{printf "%.6f", -160.0/k}'`

for g in $gravs;do
    grn=grn_g$g
    if [ -d $grn ] && [ `ls $grn 2>/dev/null | wc -l` -ge 10 ];then
	echo "$0: $grn exists, skipping PSGRN"
    else
	mkdir -p $grn
	cat > psgrn_g$g.inp <<EOF
        0.0       1
  126   0.0    500.0   1.0
    8   1.0    15.0
 $nt   $tdays
 0.05
 $g.00
 './$grn/'
 'uz'  'ur'  'ut'
 'szz' 'srr' 'stt' 'szr' 'srt' 'stz'
 'tr'  'tt'  'rot' 'gd'  'gr'
  3
1       0.0       5.2223    3.0151    3300.0     0.0E+00    0.0E+00    1.000
2      30.0       5.2223    3.0151    3300.0     0.0E+00    0.0E+00    1.000
3      30.0       4.8676    2.8103    3800.0     0.0E+00    9.468E+17  1.000
EOF
	echo "$0: PSGRN gravity $g -> $grn (log psgrn_g$g.log)"
	echo psgrn_g$g.inp | $bindir/psgrn2020 > psgrn_g$g.log 2>&1 || { echo "$0: psgrn failed, see psgrn_g$g.log"; exit 1; }
    fi
    for L in $lengths;do
	out=out_g${g}_L$L
	[ -f $out/U_down.dat ] && { echo "$0: $out exists, skipping"; continue; }
	mkdir -p $out
	nst=`awk -v l=$L 'BEGIN{printf "%d", l/5}'`
	lon0=`awk -v l=$L -v k=$KM 'BEGIN{printf "%.6f", -(l/2)/k}'`
	{
	cat <<EOF
  1
 101
  ($lat1, 0.0000), ($lat2, 0.0000)
 0   0.0  0.0  -1.0
 0   0.700  0.500  90.000  30.000   90.000    0.1E+07   -0.5E+07   -1.0E+07
'./$out/'
  1                1               1
  'U_north.dat'    'U_east.dat'    'U_down.dat'
  0            0            0            0             0            0
  'S_nn.dat'   'S_ee.dat'   'S_dd.dat'   'S_ne.dat'    'S_ed.dat'   'S_dn.dat'
  0               0               0                0              0
  'Tilt_n.dat'    'Tilt_e.dat'    'Rotation.dat'   'geoid.dat'    'Gravity.dat'
  0
  1
 './$grn/'
 'uz'  'ur'  'ut'
 'szz' 'srr' 'stt' 'szr' 'srt' 'stz'
 'tr'  'tt'  'rot' 'gd'  'gr'
  1   0.0  0.0
 1    0.0 $lon0 0.00  $L.00  30.00  90.00  30.00  $nst  15  0.00
EOF
	# patches 5 km x 2 km, 1 m up-dip (thrust) slip: slip_ddip = -1
	awk -v n=$nst 'BEGIN{for(i=0;i<n;i++)for(j=0;j<15;j++)printf "   %8.2f %7.2f   0.0000  -1.0000   0.0000\n", 2.5+5*i, 1+2*j}'
	} > pscmp_g${g}_L$L.inp
	echo "$0: PSCMP gravity $g, L = $L km -> $out"
	echo pscmp_g${g}_L$L.inp | $bindir/pscmp2020 > pscmp_g${g}_L$L.log 2>&1 || { echo "$0: pscmp failed, see pscmp_g${g}_L$L.log"; exit 1; }
    done
done
echo "$0: done; now: python3 psgrn_rundle_compare.py . \"$lengths\" \"$gravs\""
