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

make_png=1                       # 1: also render a PNG per group per frame, 0: just .grd
cpt_low=-12                      # color range, log10|v| [m/s]
cpt_high=1
cpt_step=0.5
proj=-JX15c/-7.5c                # down_dip increases downward (negative y height)

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
if [ "$make_png" -eq 1 ]; then
    gmt makecpt -Chot -T${cpt_low}/${cpt_high}/${cpt_step} -Z > $outdir/rsf_vel.cpt
fi

# loop over fault groups
for geom in rsf_geom.g*.dat; do
    if [ ! -f "$geom" ]; then
	echo "$0: no rsf_geom.g*.dat geometry files found" 1>&2
	exit 1
    fi
    tag=`echo "$geom" | sed -n 's/^rsf_geom\.\(g[0-9]*\)\.dat$/\1/p'`
    reg=`awk '/^# gmt_region/{print $3; exit}' "$geom"`
    inc=`awk '/^# gmt_inc/{print $3; exit}' "$geom"`
    if [ -n "$reg_override" ]; then reg=$reg_override; fi
    if [ -n "$inc_override" ]; then inc=$inc_override; fi
    if [ -z "$reg" ] || [ -z "$inc" ]; then
	echo "$0: $geom has no regular-grid -R/-I (irregular geometry); set reg_override/inc_override or grid with surface" 1>&2
	continue
    fi
    echo "$0: group $tag  $reg $inc"

    # loop frames for this group
    grep -v '^#' "$times" | while read frame step tyr tsec lmax rest; do
	f=`printf "%06d" "$frame"`
	bin=$outdir/rsf_vel.$tag.$f.bin
	grd=$outdir/rsf_vel.$tag.$f.grd
	if [ ! -f "$bin" ]; then
	    continue
	fi
	gmt xyz2grd "$bin" -bi3f $reg $inc -G"$grd"
	if [ "$make_png" -eq 1 ]; then
	    gmt begin $outdir/rsf_vel.$tag.$f png
	      gmt grdimage "$grd" -C$outdir/rsf_vel.cpt $proj -Baf -BWSen+t"$tag  step $step  t = $tyr yr"
	      gmt colorbar -C$outdir/rsf_vel.cpt -Bxaf+l"log@-10@-|v| [m/s]"
	    gmt end
	fi
    done
done

echo "$0: wrote $outdir/rsf_vel.gGGG.NNNNNN.grd"
echo "to animate one group, e.g.:  ffmpeg -framerate 12 -pattern_type glob -i '$outdir/rsf_vel.g000.*.png' rsf_vel.g000.mp4"
