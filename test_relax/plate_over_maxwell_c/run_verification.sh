#!/bin/bash
#
# verify the C port against the original MATLAB codes run under Octave
#
# requires: octave, the four original .m files, and roots.m (local
# wrapper restoring MATLAB root ordering) in this directory
#
# prints, for each test case, the maximum absolute difference of each
# component scaled by the maximum absolute value of that component

# test stations, matching gen_reference.m
stafile=verify_sta.txt
cat > $stafile << EOF
-150 -80
-50 40
0 60
50 10
150 -60
-100 120
100 -120
30 200
-30 -200
250 0
EOF

make pom || exit 1

echo "running original codes under Octave (this can take a few minutes)"
octave -q gen_reference.m 2>&1 | grep -v 'shadows a core'

compare () {
    # $1: reference file, $2: C output file, $3: label
    paste $1 <(grep -v '^#' $2) | awk -v lab="$3" '
	{
	    for (c = 0; c < 3; c++) {
		a = $(3+c); b = $(8+c);
		d = a-b; if (d<0) d=-d;
		if (d > dmax[c]) dmax[c] = d;
		aa = (a<0)?-a:a;
		if (aa > amax[c]) amax[c] = aa;
	    }
	}
	END {
	    printf "%-28s", lab;
	    for (c = 0; c < 3; c++) {
		s = (amax[c] > 0) ? amax[c] : 1.0;
		printf "  %.2e", dmax[c]/s;
	    }
	    printf "   (east north up)\n";
	}'
}

echo
echo "case                          scaled max difference"

./pom layer -m 200 20 20 90 0 0 0 1 0 0 -t 5 -tR1 25 -tR2 100 \
      -sta $stafile -o cmp1.txt
compare ref_L_disp_ss.txt cmp1.txt "layer disp strike-slip"

./pom layer -vel -m 200 20 20 90 0 0 0 1 0 0 -t 5 -tR1 25 -tR2 100 \
      -sta $stafile -o cmp2.txt
compare ref_L_vel_ss.txt cmp2.txt "layer vel strike-slip"

./pom layer -m 200 20 20 30 0 0 0 0 1 0 -t 5 -tR1 25 -tR2 100 \
      -sta $stafile -o cmp3.txt
compare ref_L_disp_ds30.txt cmp3.txt "layer disp dip-slip 30"

./pom layer -m 100 30 25 60 40 15 -10 0 0 1 -t 10 -tR1 25 -tR2 100 \
      -sta $stafile -o cmp4.txt
compare ref_L_disp_ten.txt cmp4.txt "layer disp tensile oblique"

./pom layer -m 200 20 20 70 20 5 5 1 0.5 0 -t 5 -tR1 25 -tR2 25 \
      -sta $stafile -o cmp5.txt
compare ref_L_disp_mix_eq.txt cmp5.txt "layer disp mixed, tR1=tR2"

i=1
for t in 1 5 10; do
    ./pom maxwell -m 200 20 20 90 0 0 0 1 0.3 0.1 -t $t -H 20 -mulam 1 \
	  -tR 25 -sta $stafile -o cmpM$i.txt
    compare ref_M_t$i.txt cmpM$i.txt "maxwell disp t=$t"
    i=$((i+1))
done

echo
echo "note: the tR1=tR2 case is intrinsically ill conditioned (nearly"
echo "repeated denominator roots); differences at the 1e-5 level there"
echo "reflect the formulation, not the translation (see README)"
