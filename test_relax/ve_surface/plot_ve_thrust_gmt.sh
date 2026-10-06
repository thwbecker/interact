#!/bin/bash
#
# plot_ve_thrust_gmt.sh: GMT 6 figures from the <prefix>_g{0,1}_profiles.txt
# tables of ve_thrust_relax.py as produced by run_thrust_relax.sh.
#
# One PDF with four panels: rows u_x (positive towards the hanging wall)
# and u_z (positive up), columns gravity off and on; curves coseismic
# (black) and the Maxwell times (dark to light), relaxed limit red dashed
# where it exists (gravity on).  A second PDF shows the full x window.
# Displacements in units of slip, x in units of H (fault trace at 0,
# dips towards +x).
#
# usage: plot_ve_thrust_gmt.sh [prefix] [times] [H_km] [title] [near_H]
#        defaults             thrust "0.1 0.5 1 2 10 100" 40 "" 3
#
pre=${1-thrust}
times=${2-"0.1 0.5 1 2 10 100"}
H=${3-40}
title=${4-""}
near=${5-3}
command -v gmt > /dev/null || { echo "$0: gmt not found"; exit 1; }
set -- $times
nt=$(($# + 1))				# coseismic + Maxwell times
# colours, coseismic first
cols=(black 68/1/84 59/82/139 33/145/140 53/183/121 144/215/67 253/231/37 255/150/0)
labels=("coseismic")
for t in $times;do labels+=("$t t@-M@-");done

one_fig ()	# one_fig <xrange_H> <suffix>
{
    xr=$1; sfx=$2
    f0=${pre}_g0_profiles.txt; f1=${pre}_g1_profiles.txt
    [ -f $f0 ] && [ -f $f1 ] || { echo "$0: missing $f0 or $f1"; return 1; }
    xmax=`awk 'NR==2{print -$1}' $f0`
    [ "$xr" = "full" ] && xr=$xmax
    gmt begin ${pre}_$sfx pdf
    gmt subplot begin 2x2 -Fs13.5c/9c -M1.5c/1.2c -SRl -SCb -A -T"thrust $title, H = $H km, elastic layer over Maxwell half-space"
    for row in 0 1;do
	for g in 0 1;do
	    f=${pre}_g${g}_profiles.txt
	    if [ $row = 0 ];then
		c0=2; ylab="u@-x@- / slip (+ towards hanging wall)"
	    else
		c0=$((2 + nt)); ylab="u@-z@- / slip (+ up)"
	    fi
	    # y range from the data of the panel (all times)
	    # (gravity on: include the relaxed curve, column 3 + 2 nt + row)
	    yr=`awk -v c0=$c0 -v nt=$nt -v xr=$xr -v g=$g -v cr=$((3 + 2 * nt + row)) 'NR>1 && $1>=-xr && $1<=xr {for(i=c0+1;i<=c0+nt;i++){if(min==""||$i<min)min=$i;if(max==""||$i>max)max=$i}; if(g==1){if($cr<min)min=$cr;if($cr>max)max=$cr}} END{d=(max-min)*0.05; printf "%g/%g", min-d, max+d}' $f`
	    [ $g = 0 ] && gtitle="gravity off" || gtitle="gravity on (interface buoyancy)"
	    gmt subplot set $row,$g -A"$gtitle"
	    gmt basemap -R-$xr/$xr/$yr -JX? -Bxaf+l"x / H" -Byaf+l"$ylab" -BWSen
	    gmt plot -W0.25p,gray <<EOF
-$xr 0
$xr 0
EOF
	    gmt plot -W0.25p,gray <<EOF
0 -1e9
0 1e9
EOF
	    for i in `seq 0 $((nt - 1))`;do
		w=1p; [ $i = 0 ] && w=1.5p
		gmt plot $f -i0,$((c0 + i)) -W$w,${cols[$i]} -l"${labels[$i]}"
	    done
	    if [ $g = 1 ];then
		cr=$((2 + 2 * nt + row))
		gmt plot $f -i0,$cr -W1p,red,- -l"relaxed (@~m@~@-2@- @~\256@~ 0)"
	    fi
	    [ $row = 0 ] && [ $g = 1 ] && gmt legend -DjTR+o0.2c -F+p0.5p+gwhite --FONT_ANNOT_PRIMARY=7p
	done
    done
    gmt subplot end
    gmt end
    echo "$0: wrote ${pre}_$sfx.pdf"
}

one_fig $near near
one_fig full far
