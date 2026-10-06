#!/bin/bash
#
# run_psgrn_convergence.sh: PSGRN/PSCMP numerical-parameter study for the
# 2-D-limit thrust case of run_psgrn_rundle.sh, to locate the 1.7 to 2.5
# percent by which PSCMP sits below the two 2-D codes at 10 and 90 tM.
#
# One Green's function set per variant, all without gravity unless the
# variant name ends in _g1 (the offset is the same with and without
# gravity, so the no-gravity case isolates the numerics):
#
#   base      as run_psgrn_rundle.sh: 126 distances 0-500 km, accuracy
#             0.05, 8 source depths 1-15 km
#   acc01     accuracy 0.01 (PSGRN's documented range is 0.1 to 0.01)
#   acc005    accuracy 0.005
#   nr252     252 distances
#   nr504     504 distances
#   nz29      29 source depths 0.5-14.5 km, i.e. exactly the patch centre
#             depths of the 5 x 2 km PSCMP patches (the base grid starts
#             at 1 km, so the shallowest patch row at 0.5 km is
#             extrapolated)
#   r1200     distances 0-1200 km (the L = 2000 km PSCMP fault has
#             patches up to 1000 km from the mid-strike receivers; the
#             base grid stops at 500 km)
#   all       252 distances 0-1200 km, accuracy 0.01, 29 depths
#   all_g1    the same with gravity on
#
# For each variant PSCMP is run with L = 2000 km and 5 x 2 km patches;
# for `all` additionally with 2.5 x 1 km patches (suffix _fine).  Each
# PSGRN run costs roughly (nr x nz / 1000) x base; `all` is about 11 x
# base.  Run the variants you want as arguments; default is the full
# list.  Then psgrn_convergence_compare.py.
#
# usage: run_psgrn_convergence.sh [bindir] [variants] [nt]
#        defaults   ../../tools/psgrn/psgrn_pscmp/bin  "base acc01 nr252 nz29 r1200 all acc005 nr504 all_g1"  128
#
bindir=${1-../../tools/psgrn/psgrn_pscmp/bin}
variants=${2-"base acc01 nr252 nz29 r1200 all acc005 nr504 all_g1"}
nt=${3-128}
bindir=`cd "$bindir" && pwd` || { echo "$0: no bindir $1"; exit 1; }
[ -x $bindir/psgrn2020 ] && [ -x $bindir/pscmp2020 ] || { echo "$0: psgrn2020/pscmp2020 not in $bindir"; exit 1; }
tdays=`awk 'BEGIN{printf "%.1f", 90*365.25}'`
KM=111.195
lat1=`awk -v k=$KM 'BEGIN{printf "%.6f", 90.0/k}'`
lat2=`awk -v k=$KM 'BEGIN{printf "%.6f", -160.0/k}'`
L=2000

pscmp_run ()	# pscmp_run <grn dir> <out dir> <patch length km> <patch width km>
{
    grn=$1; out=$2; ps=$3; pd=$4
    [ -f $out/U_down.dat ] && { echo "$0: $out exists, skipping"; return 0; }
    mkdir -p $out
    nst=`awk -v l=$L -v p=$ps 'BEGIN{printf "%d", l/p}'`
    ndi=`awk -v p=$pd 'BEGIN{printf "%d", 30/p}'`
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
 1    0.0 $lon0 0.00  $L.00  30.00  90.00  30.00  $nst  $ndi  0.00
EOF
    awk -v n=$nst -v m=$ndi -v ps=$ps -v pd=$pd 'BEGIN{for(i=0;i<n;i++)for(j=0;j<m;j++)printf "   %8.3f %7.3f   0.0000  -1.0000   0.0000\n", ps/2+ps*i, pd/2+pd*j}'
    } > pscmp_$out.inp
    echo "$0: PSCMP $out (patches $ps x $pd km)"
    echo pscmp_$out.inp | $bindir/pscmp2020 > pscmp_$out.log 2>&1 || { echo "$0: pscmp failed, see pscmp_$out.log"; exit 1; }
}

for v in $variants;do
    nr=126; r2=500.0; acc=0.05; nz=8; z1=1.0; z2=15.0; g=0
    case $v in
	base) ;;
	acc01) acc=0.01 ;;
	acc005) acc=0.005 ;;
	nr252) nr=252 ;;
	nr504) nr=504 ;;
	nz29) nz=29; z1=0.5; z2=14.5 ;;
	r1200) nr=302; r2=1200.0 ;;
	all) nr=605; r2=1200.0; acc=0.01; nz=29; z1=0.5; z2=14.5 ;;
	all_g1) nr=605; r2=1200.0; acc=0.01; nz=29; z1=0.5; z2=14.5; g=1 ;;
	*) echo "$0: unknown variant $v"; exit 1 ;;
    esac
    grn=grn_$v
    if [ -d $grn ] && [ `ls $grn 2>/dev/null | wc -l` -ge 10 ];then
	echo "$0: $grn exists, skipping PSGRN"
    else
	mkdir -p $grn
	cat > psgrn_$v.inp <<EOF
        0.0       1
  $nr   0.0    $r2   1.0
    $nz   $z1    $z2
 $nt   $tdays
 $acc
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
	echo "$0: PSGRN $v: nr $nr to $r2 km, accuracy $acc, nz $nz ($z1-$z2 km), gravity $g (log psgrn_$v.log)"
	echo psgrn_$v.inp | $bindir/psgrn2020 > psgrn_$v.log 2>&1 || { echo "$0: psgrn failed, see psgrn_$v.log"; exit 1; }
    fi
    pscmp_run $grn out_$v 5 2
    case $v in all|all_g1) pscmp_run $grn out_${v}_fine 2.5 1 ;; esac
done
echo "$0: done; now: python3 psgrn_convergence_compare.py . \"$variants\""
