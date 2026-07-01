#!/bin/bash
#
# turn the compact slip-rate field frames written by rsf_solve into GMT
# .grd (netCDF) files, one per fault group per frame, plus an optional
# PNG per frame.
#
# rsf_solve writes, when run with -field_step_interval N > 0:
#   rsf_geom.gGGG.dat                     static geometry per group, once
#                                         (header carries that group's -R/-I)
#   tmp_rsf/rsf_vel.gGGG.NNNNNN.bin        per-frame field per group: float32
#                                         (along_strike[m], down_dip[m], log10|v|[m/s])
#   rsf_vel.times                          index: frame step time[yr] time[s] log10(max|v|)
#
# frames are written every N accepted steps, so they densify through
# nucleation and rupture and thin out through the interseismic.
# each group/frame is a single xyz2grd call.  requires GMT 6 on PATH.
#
# parameters are hardcoded plain assignments below; edit them here.

times=rsf_vel.times              # frame index
outdir=tmp_rsf                   # where the .bin frames are and the .grd go
wdir=`pwd`
tdir=`basename $wdir`

make_png=2                       # 1: also render a PNG per group per frame, 0: just .grd 2: make png and delete grid

proj=-Jx.125/-.125                # down_dip increases downward (negative y height)

proj2=-JX12.5/3
y2min=-9.25;y2max=0.25
reg2=-R0/3000/$y2min/$y2max

# -R / -I are read per group from each rsf_geom.gGGG.dat header, so no
# region needs to be set here.  override below only if a header lacks a
# grid line (irregular geometry); leave empty to use the header value.
reg_override=
inc_override=

# --- nothing below here is a parameter -------------------------------------

if [ ! -f "$times" ]; then
    echo "$0: $times not found; run rsf_solve with -field_step_interval N > 0 first" 1>&2
    exit 1
fi
if [ "$make_png" -ne 0 ]; then
    gmt makecpt -Croma -T-11/0/0.5 -I > $outdir/rsf_vel.cpt
fi

# loop over fault groups
for geom in rsf_geom.g*.dat; do
    if [ ! -f "$geom" ]; then
	echo "$0: no rsf_geom.g*.dat geometry files found" 1>&2
	exit 1
    fi
    tag=`echo "$geom" | sed -n 's/^rsf_geom\.\(g[0-9]*\)\.dat$/\1/p'`
    reg=`awk '/^# gmt_region/{print $3; exit}' "$geom"`
    regp=`echo $reg | gawk -f reg2wesn.awk | gawk '{printf("-R%g/%g/%g/%g",$1/1e3,$2/1e3,$3/1e3,$4/1e3)}'`
    inc=`awk '/^# gmt_inc/{print $3; exit}' "$geom"`
    echo
    echo $reg $regp $inc
    echo
    if [ -n "$reg_override" ]; then reg=$reg_override; fi
    if [ -n "$inc_override" ]; then inc=$inc_override; fi
    if [ -z "$reg" ] || [ -z "$inc" ]; then
	echo "$0: $geom has no regular-grid -R/-I (irregular geometry); set reg_override/inc_override or grid with surface" 1>&2
	continue
    fi
    echo "$0: group $tag  $reg $inc"

    # loop frames for this group
    grep -v '^#' "$times" | while read frame step tyr tsec lmax rest; do
	tyrs=`echo $tyr | gawk '{printf("%06.1f",$1)}'`
	f=`printf "%06d" "$frame"`
	bin=$outdir/rsf_vel.$tag.$f.bin
	grd=$outdir/rsf_vel.$tag.$f.grd
	if [ ! -f "$bin" ]; then
	    continue
	fi
	gmt xyz2grd "$bin" -bi3f $reg $inc -G"$grd" -r -V
	grdedit $grd $regp -fc	# rescale to km
	grdinfo -C $grd
	if [ "$make_png" -ne 0 ]; then
	    ofile=$outdir/rsf_vel.$tag.$f.ps
	    
 	    grdimage "$grd" $regp -Ba5f1:"x [km]":/a5f1:"z [km]"::."$tag, step $f, t = $tyrs yr":WesN -Y5 $proj -C$outdir/rsf_vel.cpt  -P -K > $ofile
 	    psscale -E -C$outdir/rsf_vel.cpt -D6/-.25/4/.15h -B2/:"log@-10@-(|v| [m/s])": -O -K >> $ofile
	    psbasemap -Y-4 $proj2 $reg2 -Ba500f100:"time [yr]":/a2f.2:"log@-10@-(@~\341@~|v| [m/s]@~\361@~)":WeSn -O -K >> $ofile
	    cat <<EOF | psxy $reg2 $proj2 -W4,gray -O -K >> $ofile
$tyr $y2min
$tyr $y2max

EOF
	    gawk '{print($3,$5)}' rsf_monitor.dat  | psxy $reg2 $proj2 -W1,blue -O -K >> $ofile
	    psxy -T -O >> $ofile
	    modifybb $ofile 2

	    echo $0: written to $outdir/rsf_vel.$tag.$f.png
	    #eog  $outdir/rsf_vel.$tag.$f.png; exit
	fi
	if [ "$make_png" -eq 2 ]; then
	    rm $grd
	fi
    done
done

echo "$0: wrote $outdir/rsf_vel.gGGG.NNNNNN.grd"
echo "to animate one group, e.g.:  ffmpeg -framerate 12 -pattern_type glob -i '$outdir/rsf_vel.g000.*.png' rsf_vel.g000.$tdir.mp4"
rm  rsf_vel.g000.$tdir.mp4
ffmpeg -framerate 12 -pattern_type glob -i '$outdir/rsf_vel.g000.*.png' rsf_vel.g000.$tdir.mp4
cp rsf_vel.g000.$tdir.mp4 $HOME/Dropbox/tmp/
