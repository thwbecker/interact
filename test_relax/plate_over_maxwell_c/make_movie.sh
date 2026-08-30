#!/bin/bash
#
# make a movie of postseismic surface deformation with GMT and ffmpeg
#
# for each time step: compute the deformation field with the C port,
# plot vertical deformation as a color field and horizontal deformation
# as vectors, then assemble the frames into mp4 (and gif)
#
# all parameters are set here as plain assignments
for  model in strike thrust normal;do
    if [ $model = strike ];then
	# fault: L W depth dip strike east north strike-slip dip-slip tensile
	fault="200 20 20 90 0 0 0 1  0 0"
    elif [ $model = thrust   ];then
	fault="200 20 20 60 0 0 0 0  1 0"
    elif [ $model = normal ];then
	fault="200 20 20 30 0 0 0 0 -1 0"
    fi


    # layered model
    h1=20
    h2=40
    nu=0.25
    tr1=25
    tr2=25
    # what to show: disp or vel
    field=disp
    # grid for the color field (km); nx=ny=61 gives 10 km spacing
    xmin=-300
    xmax=300
    ymin=-300
    ymax=300
    #ngrid=61
    ngrid=122
    
    # vector decimation: keep every nth grid node in each direction
    #vdec=4
    vdec=8
    # source quadrature; 40 and 8 suppress the near-fault checkerboard from
    # the discrete point sources, at about 32 times the cost of the
    # original defaults (set both to 0 for the original behavior)
    nl=40
    nw=8
    # times: tstart, tend, dt (years)
    tstart=2
    tend=100
    dt=2
    # vector scale (cm of arrow per m); adjust for vel
    dscale=5
    # frame rate of the movie
    fps=8
    
    outdir=$model"_movie_out"
    
    mkdir -p $outdir
    rm -f $outdir/frame_*.png
# field files are kept and reused if present; delete the movie_out
# directory when changing fault or model parameters
    
    dx=$(echo "($xmax - $xmin) / ($ngrid - 1)" | bc -l)
    region=-R$xmin/$xmax/$ymin/$ymax
    proj=-JX8
    
    velflag=""
    label="displacement"
    unit="m"
    if [ $field = "vel" ]; then
	velflag="-vel"
	label="velocity"
	unit="m/yr"
    fi
    
    # pass 1: compute all fields
    tlist=$(seq $tstart $dt $tend)
    j=0
    for t in $tlist; do
	if [ ! -s $outdir/field_$t.txt ]; then
	    (./pom layer $velflag -m $fault -H1 $h1 -H2 $h2 -nu $nu -tR1 $tr1 -tR2 $tr2 -nl $nl -nw $nw -grid $xmin $xmax $ngrid $ymin $ymax $ngrid -t $t  -o $outdir/field_$t.txt) &
	    echo "computing t = $t"
	    ((j=j+1))
	    if [ $j -ge $NR_CPUS ];then
		wait
		j=0
	    fi
	fi
    done
    wait
    
    # pass 2: global symmetric color range for the vertical component
    zmax=$(awk 'FNR>1 {a=($5<0)?-$5:$5; if(a>m) m=a} END {printf "%.6f", m*1.05}' $outdir/field_*.txt)
    echo "vertical color range: +/- $zmax $unit"
    zinc=`echo $zmax | gawk '{print($1*2/21)}'`
    gmt makecpt -Cvik -T-$zmax/$zmax > $outdir/z.cpt

    # pass 3: plot frames
    i=0
    for t in $tlist; do
	frame=$(printf "%s/frame_%04d" $outdir $i)
	gmt begin $frame png
	awk 'NR>1 {print $1, $2, $5}' $outdir/field_$t.txt | \
	    gmt xyz2grd $region -I$dx -G$outdir/uz.grd
	gmt grdimage $outdir/uz.grd $region $proj -C$outdir/z.cpt
	gmt basemap $region $proj -Baf -BWSen+t"postseismic $label, t = $t yr"
	# fault trace (strike 0 fault centered on the origin)
	printf "0 -100\n0 100\n" | gmt plot -W2p,black
	# horizontal vectors, decimated
	awk -v s=$dscale -v d=$vdec -v n=$ngrid 'NR>1 {
	r = int((NR-2)/n); c = (NR-2)%n;	
	if (r%d==0 && c%d==0) {			
	    l = sqrt($3*$3+$4*$4)*s;		
	    if (l > 0.02) print $1, $2, atan2($4,$3)*180/3.14159265, l
	}
	    }' $outdir/field_$t.txt | gmt plot -Sv0.2c+e -W1p,black -Gblack
	gmt colorbar -C$outdir/z.cpt -DJBC+w8c/0.4c+h -Baf+l"vertical $label ($unit)"
	gmt end
	i=$((i+1))
	echo "frame $i / $(echo $tlist | wc -w)"
	rm $outdir/field_$t.txt 
    done
    
    # pass 4: assemble
    #ffmpeg -y -loglevel error -framerate $fps -i $outdir/frame_%04d.png \
#	   -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
#	   -c:v libx264 -pix_fmt yuv420p $outdir/${field}_movie.mp4
    ffmpeg -y -loglevel error -framerate $fps -i $outdir/frame_%04d.png \
	   -vf "scale=640:-1" $outdir/$model"_"$field"_movie.gif"
    rm $outdir/frame*.png 
    echo wrote $outdir/$model"_"$field"_movie.gif"
    cp  $outdir/$model"_"$field"_movie.gif" $HOME/Dropbox/tmp/
    
done
