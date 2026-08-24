#!/bin/bash
#
# reproduce run_postseismic_example.m with the C port and plot the
# cumulative postseismic displacement and velocity vectors with GMT
#
# all parameters are set here as plain assignments
for model in strike thrust normal;do
    if [ $model = strike ];then
	# fault: L W depth dip strike east north strike-slip dip-slip tensile
	fault="200 20 20 90 0 0 0 1  0 0"
    elif [ $model = thrust   ];then
	fault="200 20 20 60 0 0 0 0  1 0"
    elif [ $model = normal ];then
	fault="200 20 20 30 0 0 0 0 -1 0"
    fi
    echo $0: model $model fault: $fault
    # layered model
    h1=20
    h2=40
    nu=0.25
    tr1=25
    tr2=25
    # source quadrature (see make_movie.sh; 0 0 restores original defaults)
    nl=40
    nw=8
    # output grid (km)
    grid="-300 300 24 -300 300 24"
    # times (years)
    times="1 5 10 50 100"
    # vector scale: cm of arrow per m of displacement, per m/yr of velocity
    dscale=3
    vscale=200
    
    outdir=$model"_out"
    
    mkdir -p $outdir
    
    for t in $times; do
	./pom layer -nl $nl -nw $nw -m $fault -H1 $h1 -H2 $h2 -nu $nu -tR1 $tr1 -tR2 $tr2 \
	      -grid $grid -t $t -o $outdir/disp_$t.txt
	./pom layer -vel -nl $nl -nw $nw -m $fault -H1 $h1 -H2 $h2 -nu $nu -tR1 $tr1 -tR2 $tr2 \
	      -grid $grid -t $t -o $outdir/vel_$t.txt
    done
    
    region=-R-300/300/-300/300
    proj=-JX12c
    
    plot_panel () {
	# $1: file prefix (disp or vel), $2: scale, $3: title
	local pre=$1
	local scale=$2
	local title=$3
	gmt begin $outdir/${pre}_vectors png
	gmt basemap $region $proj -Bafg100 -B+t"$title"
	# fault trace: vertical strike-slip fault along north at x=0
	printf "0 -100\n0 100\n" | gmt plot -W2p,red
	local colors="black blue darkgreen orange purple"
	local c
	local t
	set -- $times
	for c in $colors; do
	    t=$1; shift
	    awk -v s=$scale 'NR>1 {
	    l = sqrt($3*$3+ $4*$4)*s;
	    if (l > 0.02) print $1, $2, atan2($4,$3)*180/3.14159265, l
	}' $outdir/${pre}_$t.txt     | gmt plot -Sv0.15c+e -W0.5p,$c -G$c
	done
	gmt end
    }
    
    plot_panel disp $dscale "cumulative postseismic displacement, t = $times yr"
    plot_panel vel $vscale "postseismic velocity, t = $times yr"
    
    echo "wrote $outdir/disp_vectors.png and $outdir/vel_vectors.png"
done
