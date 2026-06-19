#!/bin/bash
#
# cluster_compression_test.sh
#
# Two related questions about H-matrix compression of interact operators,
# on a panel of makefault test geometries:
#
#   1. CLUSTERING: does isolating faults from each other before clustering
#      help? For each two-fault case we dump the dense operator and the
#      patch coordinates (compress_interaction_matrix -dump_matrix /
#      -dump_coords) and run cluster_compare.py, which builds a simple
#      H-matrix two ways at matched accuracy: a single geometric tree over
#      all patches (joint, what the real backends do) versus splitting by
#      fault first and clustering geometrically within each (split, the
#      block H-matrix F1F1/F2F2/F1F2/F2F1). It reports the storage of each.
#
#   2. COMPRESSIBILITY: how much do the actual backends (HTOOL, HACApK,
#      hmmvp) compress each geometry? We run them on the same operators and
#      tabulate their compression ratio. The backends cluster on centroids
#      only and cannot split by fault, so this is the joint number; the
#      clustering table above is what a fault-aware tree could add.
#
# The cases are chosen to bracket the effect: a single fault, two parallel
# same-orientation faults, two different-orientation faults that diverge in
# space, two perpendicular faults intersecting through one box (the genuine
# interleaving case), and two different-orientation faults far apart.
#
# All numbers are specific to these geometries, this accuracy band, and the
# admissibility/leaf settings; a different operator, size, tolerance, or a
# production library may shift them. See cluster_compression.md.
#
# IMPORTANT: compress_interaction_matrix must be built from a source that
# has the -dump_matrix / -dump_coords flags. If it is not, the dump writes
# nothing and this script stops with a clear message; rebuild the binary
# (make bin/compress_interaction_matrix) and rerun.
#
# Usage:
#   ./cluster_compression_test.sh [SCALE] [np]
#     SCALE : patch-count knob, strike segments per fault (default 24);
#             dip segments scale with it. Bump it (e.g. 48, 64) to push the
#             backends into a regime where they compress meaningfully.
#     np    : MPI ranks for the backend survey (default 1). The dumps always
#             use 1 rank so the full matrix is local.
#
# All other parameters are set as plain assignments below; edit them here.
# This script intentionally does NOT read parameters from environment
# variables.

set -u

# ---- parameters (edit here; no environment variables) ----------------------
SCALE=24                  # default strike segments per fault; arg 1 overrides
NP=1                      # default MPI ranks for the survey; arg 2 overrides
[ $# -ge 1 ] && SCALE=$1
[ $# -ge 2 ] && NP=$2

ETA=2.0                   # admissibility constant for cluster_compare.py
TOL=1e-5                  # low-rank relative tolerance
LEAF=8                    # leaf cluster size
BACKENDS="1 3 4"          # survey backends: 1=HTOOL 3=HACApK 4=hmmvp (0=dense)
DO_SURVEY=1               # 1 run the backend survey, 0 skip
KEEP=0                    # 1 keep the per-case dumps, 0 clean up

# tools, relative to the directory this script is run from (repo subdir)
MF=../bin/makefault
CIM=../bin/compress_interaction_matrix
CCCOMP=cluster_compare.py
MPIRUN=mpirun

for f in "$MF" "$CIM" "$CCCOMP"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

# strike/dip segment counts (dip ~ a third of strike, min 6)
NS=$SCALE
ND=$(( SCALE/3 )); [ $ND -lt 6 ] && ND=6
# half-length / half-width of each fault (long and thin gives far-field reach)
L=14; W=6
WORK=$(mktemp -d)
trap '[ "$KEEP" = 1 ] || rm -rf "$WORK"' EXIT

mkfault(){ # x y z strike dip grp  -> stdout patches
  $MF -n $NS -m $ND -l $L -w $W -x "$1" -y "$2" -z "$3" -strike "$4" -dip "$5" -grp "$6" 2>/dev/null
}

# ---- build the panel -------------------------------------------------------
# each case: name + the geometry file; two-fault cases get a clustering test
declare -a CASES

# single long fault (baseline far-field, one group)
{ mkfault 0 0 -8 0 90 0; }                         > "$WORK/single.in";   CASES+=("single")
# two parallel, SAME orientation, close (geometric tree separates trivially)
{ mkfault -1 0 -8 0 90 0; mkfault 1 0 -8 0 90 1; } > "$WORK/parallel.in"; CASES+=("parallel")
# two DIFFERENT orientation but spatially divergent (dip 90 vs dip 45)
{ mkfault 0 0 -8 0 90 0; mkfault 0 0 -8 0 45 1; }  > "$WORK/divergent.in";CASES+=("divergent")
# two PERPENDICULAR faults intersecting through one box (genuine interleaving)
{ mkfault 0 0 -8 0 90 0; mkfault 0 0 -8 0 0  1; }  > "$WORK/cross.in";    CASES+=("cross")
# two different-orientation faults FAR apart (both trees separate them)
{ mkfault 0 0 -8 0 90 0; mkfault 0 60 -8 0 45 1; } > "$WORK/separated.in";CASES+=("separated")

echo "panel scale: NS=$NS ND=$ND per fault; admissibility eta=$ETA tol=$TOL leaf=$LEAF"
echo

# ---- 1. clustering comparison (two-fault cases) ----------------------------
echo "=== 1. clustering: joint geometric tree vs split-by-fault (matched accuracy) ==="
printf "%-11s %6s %10s %10s %9s  %s\n" case N joint_x split_x gain verdict
for c in "${CASES[@]}"; do
  g="$WORK/$c.in"
  ngrp=$(awk '{print $8}' "$g" | sort -u | wc -l)
  [ "$ngrp" -lt 2 ] && continue   # clustering test needs two faults

  # dump the dense operator and the coordinates. keep the run log so that a
  # missing dump can be explained rather than crashing cluster_compare.py.
  dlog="$WORK/$c.dumplog"
  $MPIRUN -np 1 "$CIM" -geom_file "$g" -use_hmatrix 0 \
      -dump_matrix "$WORK/$c.bin" -dump_coords "$WORK/$c.coords" > "$dlog" 2>&1
  if [ ! -s "$WORK/$c.bin.info" ] || [ ! -s "$WORK/$c.bin" ]; then
    echo >&2
    echo "ERROR: the dump for case '$c' produced no matrix file." >&2
    echo "       Expected $WORK/$c.bin and .info, but they are missing or empty." >&2
    echo "       The likely cause is that" >&2
    echo "         $CIM" >&2
    echo "       was not built with the -dump_matrix / -dump_coords flags." >&2
    echo "       Rebuild it from the updated src/compress_interaction_matrix.c:" >&2
    echo "         make bin/compress_interaction_matrix" >&2
    echo "       then rerun. Last lines of the dump run:" >&2
    tail -n 6 "$dlog" | sed 's/^/         | /' >&2
    exit 1
  fi

  out=$(python3 "$CCCOMP" -m "$WORK/$c.bin" -c "$WORK/$c.coords" -tol "$TOL" -eta "$ETA" -leaf "$LEAF")
  N=$(echo "$out"   | awk '/^operator:/{print $2}')
  jx=$(echo "$out"  | awk '$1=="joint"{print $3}')
  sx=$(echo "$out"  | awk '$1=="split"{print $3}')
  if [ -z "$jx" ] || [ -z "$sx" ]; then
    Nshow="$N"; [ -z "$Nshow" ] && Nshow="?"
    printf "%-11s %6s %10s %10s %9s  %s\n" "$c" "$Nshow" "?" "?" "?" "cluster_compare failed"
    continue
  fi
  jv=${jx%x}; sv=${sx%x}
  gain=$(awk -v j="$jv" -v s="$sv" 'BEGIN{printf "%+.1f%%", 100*(s-j)/j}')
  verdict=$(awk -v j="$jv" -v s="$sv" 'BEGIN{d=100*(s-j)/j; if(d<-1)print "split worse"; else if(d>1)print "split helps"; else print "geometric ok"}')
  printf "%-11s %6s %10s %10s %9s  %s\n" "$c" "$N" "$jx" "$sx" "$gain" "$verdict"
done
echo
echo "  joint_x/split_x = dense/stored compression; gain = (split-joint)/joint."
echo "  positive gain means isolating the faults stores less (helps)."

# ---- 2. backend compressibility survey -------------------------------------
[ "$DO_SURVEY" = 1 ] || { echo; echo "(backend survey skipped: DO_SURVEY=0)"; exit 0; }
echo
echo "=== 2. backend compressibility on the joint operator (np=$NP) ==="
name_of(){ case $1 in 0) echo dense;; 1) echo htool;; 3) echo hacapk;; 4) echo hmmvp;; *) echo "hm$1";; esac; }
printf "%-11s %6s" case N; for b in $BACKENDS; do printf " %9s" "$(name_of "$b")"; done; printf "\n"
for c in "${CASES[@]}"; do
  g="$WORK/$c.in"
  N=$(wc -l < "$g")
  printf "%-11s %6s" "$c" "$N"
  for b in $BACKENDS; do
    log=$($MPIRUN -np "$NP" "$CIM" -geom_file "$g" -use_hmatrix "$b" 2>&1)
    case $b in
      1) cx=$(echo "$log" | awk '/compression ratio:/{print $3}' | tail -1);;
      4) cx=$(echo "$log" | grep -oE 'compression ratio [0-9.]+' | awk '{print $3}' | tail -1);;
      3) pct=$(echo "$log" | awk -F= '/Memory compression v.s. dense/{print $2}' | awk '{print $1}' | tail -1)
         cx=$(awk -v p="$pct" 'BEGIN{if(p>0)printf "%.2f",100/p; else print "na"}');;
      *) cx="na";;
    esac
    [ -z "$cx" ] && cx="na"
    printf " %9s" "$cx"
  done
  printf "\n"
done
echo
echo "  values are compression ratio (dense/H-matrix); >1 compresses, <=1 near-dense."
echo "  small N compresses little for every backend; raise SCALE to see the trend."
